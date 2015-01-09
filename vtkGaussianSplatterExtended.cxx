/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkGaussianSplatterExtended.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkGaussianSplatterExtended.h"

#include "vtkDoubleArray.h"
#include "vtkImageData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkPointData.h"

#include "vtkMultiThreader.h"
#include "vtkCriticalSection.h"

#include "vtkUnstructuredGrid.h"
#include "vtkStructuredGrid.h"
#include "vtkRectilinearGrid.h"

#include <math.h>

vtkStandardNewMacro(vtkGaussianSplatterExtended);

struct vtkGaussianSplatterExtendedThreadInfo
{
  vtkGaussianSplatterExtended *Modeller;
  vtkDataSet               **Input;
};

// Construct object with dimensions=(50,50,50); automatic computation of
// bounds; a splat radius of 0.1; an exponent factor of -5; and normal and
// scalar warping turned on.
vtkGaussianSplatterExtended::vtkGaussianSplatterExtended()
{
  this->SampleDimensions[0] = 50;
  this->SampleDimensions[1] = 50;
  this->SampleDimensions[2] = 50;

  this->Radius = 0.1;
  this->ExponentFactor = -5.0;

  this->ModelBounds[0] = 0.0;
  this->ModelBounds[1] = 0.0;
  this->ModelBounds[2] = 0.0;
  this->ModelBounds[3] = 0.0;
  this->ModelBounds[4] = 0.0;
  this->ModelBounds[5] = 0.0;

  this->NormalWarping = 1;
  this->Eccentricity = 2.5;

  this->ScalarWarping = 1;
  this->ScaleFactor = 1.0;

  this->Capping = 1;
  this->CapValue = 0.0;

  this->AccumulationMode = VTK_ACCUMULATION_MODE_MAX;
  this->NullValue = 0.0;


  this->Threader        = vtkMultiThreader::New();
  this->NumberOfThreads = this->Threader->GetNumberOfThreads();
  this->ProgressMutex = new vtkSimpleCriticalSection;
}

//----------------------------------------------------------------------------
vtkGaussianSplatterExtended::~vtkGaussianSplatterExtended()
{
  if (this->Threader)
    {
    this->Threader->Delete();
    }

  if (this->ProgressMutex)
    {
    delete this->ProgressMutex;
    }
}

//----------------------------------------------------------------------------
int vtkGaussianSplatterExtended::RequestInformation (
  vtkInformation * vtkNotUsed(request),
  vtkInformationVector ** vtkNotUsed( inputVector ),
  vtkInformationVector *outputVector)
{
  // get the info objects
  vtkInformation* outInfo = outputVector->GetInformationObject(0);

  // use model bounds if set
  this->Origin[0] = 0;
  this->Origin[1] = 0;
  this->Origin[2] = 0;
  if ( this->ModelBounds[0] < this->ModelBounds[1] &&
       this->ModelBounds[2] < this->ModelBounds[3] &&
       this->ModelBounds[4] < this->ModelBounds[5] )
    {
    this->Origin[0] = this->ModelBounds[0];
    this->Origin[1] = this->ModelBounds[2];
    this->Origin[2] = this->ModelBounds[4];
    }

  outInfo->Set(vtkDataObject::ORIGIN(), this->Origin, 3);

  int i;
  for (i=0; i<3; i++)
    {
    this->Spacing[i] = (this->ModelBounds[2*i+1] - this->ModelBounds[2*i])
      / (this->SampleDimensions[i] - 1);
    if ( this->Spacing[i] <= 0.0 )
      {
      this->Spacing[i] = 1.0;
      }
    }
  outInfo->Set(vtkDataObject::SPACING(),this->Spacing,3);

  outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(),
               0, this->SampleDimensions[0] - 1,
               0, this->SampleDimensions[1] - 1,
               0, this->SampleDimensions[2] - 1);
  vtkDataObject::SetPointDataActiveScalarInfo(outInfo, VTK_DOUBLE, 1);
  return 1;
}

//----------------------------------------------------------------------------
int vtkGaussianSplatterExtended::RequestData(
  vtkInformation* vtkNotUsed( request ),
  vtkInformationVector** inputVector,
  vtkInformationVector* outputVector)
{
  // get the data object
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  vtkImageData *output = vtkImageData::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  output->SetExtent(
    outInfo->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT()));
  output->AllocateScalars(outInfo);

  vtkIdType numPts;

  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkDataSet *input = vtkDataSet::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));
  int sliceSize=this->SampleDimensions[0]*this->SampleDimensions[1];

  // Check that the input type is supported
  int dataObjectType = input->GetDataObjectType();
  if (dataObjectType != VTK_STRUCTURED_GRID &&
      dataObjectType != VTK_UNSTRUCTURED_GRID &&
      dataObjectType != VTK_RECTILINEAR_GRID)
    {
    vtkErrorMacro(<< "vtkPartialVolumeModeller expects an input data set of type "
                  << "vtkStructuredGrid, vtkUnstructuredGrid, or vtkRectilinearGrid.");
    return 0;
    }


  vtkGaussianSplatterExtendedThreadInfo info;
  info.Modeller = this;

  numPts=input->GetNumberOfPoints();

  // Set the number of threads to use
  if ( this->NumberOfThreads > numPts)
    {
    this->Threader->SetNumberOfThreads(numPts);
    }
  else
    {
    this->Threader->SetNumberOfThreads(this->NumberOfThreads);
    }

  //[[[deep copy likely not needed]]]
  info.Input = new vtkDataSet*[this->Threader->GetNumberOfThreads()];
  for (int threadId = 0; threadId < this->Threader->GetNumberOfThreads(); threadId++)
    {
    switch (dataObjectType)
      {
      case VTK_STRUCTURED_GRID:
        info.Input[threadId] = vtkStructuredGrid::New();
        break;

      case VTK_UNSTRUCTURED_GRID:
        info.Input[threadId] = vtkUnstructuredGrid::New();
        break;

      case VTK_RECTILINEAR_GRID:
        info.Input[threadId] = vtkRectilinearGrid::New();
        break;
      }
    info.Input[threadId]->DeepCopy(input);
    }

  this->Threader->SetSingleMethod( vtkGaussianSplatterExtended::ThreadedExecute,
    (void *)&info);
  this->TotalProgress = 0.0;
  this->Threader->SingleMethodExecute();

  // Clean up.
  for (int threadId = 0; threadId < this->Threader->GetNumberOfThreads(); threadId++)
    {
    info.Input[threadId]->Delete();
    }
  delete[] info.Input;

  return 1;
}

//----------------------------------------------------------------------------
VTK_THREAD_RETURN_TYPE vtkGaussianSplatterExtended::ThreadedExecute( void *arg )
{
  int threadId = ((vtkMultiThreader::ThreadInfo *)(arg))->ThreadID;
  int threadCount = ((vtkMultiThreader::ThreadInfo *)(arg))->NumberOfThreads;
  vtkGaussianSplatterExtendedThreadInfo *userData = (vtkGaussianSplatterExtendedThreadInfo *)
    (((vtkMultiThreader::ThreadInfo *)(arg))->UserData);

  vtkDataSet *input = userData->Input[threadId];
  vtkImageData *output = userData->Modeller->GetOutput();

  vtkIdType numPts, numNewPts, ptId, idx, i;
  int j, k;
  int min[3], max[3];
  vtkPointData *pd;
  vtkDataArray *inNormals=NULL;
  double loc[3], dist2, cx[3];
  vtkDoubleArray *newScalars =
    vtkDoubleArray::SafeDownCast(output->GetPointData()->GetScalars());
  newScalars->SetName("SplatterValues");


  vtkDebugMacro(<< "Splatting data");

  //  Make sure points are available
  //
  if ( (numPts=input->GetNumberOfPoints()) < 1 )
    {
    vtkDebugMacro(<<"No points to splat!");
    return VTK_THREAD_RETURN_VALUE;
    }

  //  Compute the radius of influence of the points.  If an
  //  automatically generated bounding box has been generated, increase
  //  its size slightly to acoomodate the radius of influence.
  //
  this->Eccentricity2 = this->Eccentricity * this->Eccentricity;

  numNewPts = this->SampleDimensions[0] * this->SampleDimensions[1] *
              this->SampleDimensions[2];
  for (i=0; i<numNewPts; i++)
    {
    newScalars->SetTuple(i,&this->NullValue);
    }
  this->Visited = new char[numNewPts];
  for (i=0; i < numNewPts; i++)
    {
    this->Visited[i] = 0;
    }

  output->SetDimensions(this->GetSampleDimensions());
  this->ComputeModelBounds(input,output, outInfo);

  //decide which array to splat, if any
  pd = input->GetPointData();
  int association = vtkDataObject::FIELD_ASSOCIATION_POINTS;
  vtkDataArray *inScalars=this->GetInputArrayToProcess
    (0, inputVector, association);
  if (!inScalars)
    {
    inScalars = pd->GetScalars();
    }

  //  Set up function pointers to sample functions
  //
  if ( this->NormalWarping && (inNormals=pd->GetNormals()) != NULL )
    {
    this->Sample = &vtkGaussianSplatterExtended::EccentricGaussian;
    }
  else
    {
    this->Sample = &vtkGaussianSplatterExtended::Gaussian;
    }

  if ( this->ScalarWarping && inScalars != NULL )
    {
    this->SampleFactor = &vtkGaussianSplatterExtended::ScalarSampling;
    }
  else
    {
    this->SampleFactor = &vtkGaussianSplatterExtended::PositionSampling;
    this->S = 0.0; //position sampling does not require S to be defined
                   //but this makes purify happy.
    }

  // Traverse all points - splatting each into the volume.
  // For each point, determine which voxel it is in.  Then determine
  // the subvolume that the splat is contained in, and process that.
  //
  int abortExecute=0;
  vtkIdType progressInterval = numPts/100 + 1;
  for (ptId=0; ptId < numPts && !abortExecute; ptId++)
    {
    if ( ! (ptId % progressInterval) )
      {
      vtkDebugMacro(<<"Inserting point #" << ptId);
      this->UpdateProgress (static_cast<double>(ptId)/numPts);
      abortExecute = this->GetAbortExecute();
      }

    userData->Modeller->UpdateThreadProgress(voxelProgressWeight*count);
    if (threadId == 0)
	{
	userData->Modeller->UpdateProgress( userData->Modeller->TotalProgress );
	}

    this->P = input->GetPoint(ptId);
    if ( inNormals != NULL )
      {
      this->N = inNormals->GetTuple(ptId);
      }
    if ( inScalars != NULL )
      {
      this->S = inScalars->GetComponent(ptId,0);
      }

    // Determine the voxel that the point is in
    for (i=0; i<3; i++)
      {
      loc[i] = (this->P[i] - this->Origin[i]) / this->Spacing[i];
      }

    if (this->Radius == 0)
	{
        char splatit= 1;
	// Determine splat footprint
	for (i=0; i<3; i++)
	    {
	    min[i] = static_cast<int>(round(loc[i]));
	    if ( min[i] < 0 || min[i] >= this->SampleDimensions[i] )
		{
		splatit= 0;
		}
	    } 
	if (splatit){
	    idx = min[0] + min[1]*this->SampleDimensions[0] + min[2]*sliceSize;
	    this->SetScalar(idx,0, newScalars);
	    }                
	}
    else {
        // Determine splat footprint
        for (i=0; i<3; i++)
            {
            min[i] = static_cast<int>(round(loc[i] - this->SplatDistance[i]));
            max[i] = static_cast<int>(round(loc[i] + this->SplatDistance[i]));
      if ( min[i] < 0 )
        {
        min[i] = 0;
        }
      if ( max[i] >= this->SampleDimensions[i] )
        {
        max[i] = this->SampleDimensions[i] - 1;
        }
      }

    // Loop over all sample points in volume within footprint and
    // evaluate the splat
    for (k=min[2]; k<=max[2]; k++)
      {
      cx[2] = this->Origin[2] + this->Spacing[2]*k;
      for (j=min[1]; j<=max[1]; j++)
        {
        cx[1] = this->Origin[1] + this->Spacing[1]*j;
        for (i=min[0]; i<=max[0]; i++)
          {
          cx[0] = this->Origin[0] + this->Spacing[0]*i;
          if ( (dist2=(this->*Sample)(cx)) <= this->Radius2 )
            {
            idx = i + j*this->SampleDimensions[0] + k*sliceSize;
            this->SetScalar(idx,dist2, newScalars);
            }//if within splat radius
          }
        }
      }//within splat footprint
     }//Radius != 0
    }//for all input points

  // If capping is turned on, set the distances of the outside of the volume
  // to the CapValue.
  //
  if ( this->Capping )
    {
    this->Cap(newScalars);
    }

  vtkDebugMacro(<< "Splatted " << input->GetNumberOfPoints() << " points");

  // Update self and release memeory
  //
  delete [] this->Visited;

  return VTK_THREAD_RETURN_VALUE;
}

//----------------------------------------------------------------------------
// Compute the size of the sample bounding box automatically from the
// input data.
void vtkGaussianSplatterExtended::ComputeModelBounds(vtkDataSet *input,
                                             vtkImageData *output,
                                             vtkInformation *outInfo)
{
  double *bounds, maxDist;
  int i, adjustBounds=0;

  // compute model bounds if not set previously
  if ( this->ModelBounds[0] >= this->ModelBounds[1] ||
       this->ModelBounds[2] >= this->ModelBounds[3] ||
       this->ModelBounds[4] >= this->ModelBounds[5] )
    {
    adjustBounds = 1;
    bounds = input->GetBounds();
    }
  else
    {
    bounds = this->ModelBounds;
    }

  for (maxDist=0.0, i=0; i<3; i++)
    {
    if ( (bounds[2*i+1] - bounds[2*i]) > maxDist )
      {
      maxDist = bounds[2*i+1] - bounds[2*i];
      }
    }
  maxDist *= this->Radius;
  this->Radius2 = maxDist * maxDist;

  // adjust bounds so model fits strictly inside (only if not set previously)
  if ( adjustBounds )
    {
    for (i=0; i<3; i++)
      {
      this->ModelBounds[2*i] = bounds[2*i] - maxDist;
      this->ModelBounds[2*i+1] = bounds[2*i+1] + maxDist;
      }
    }

  // Set volume origin and data spacing
  outInfo->Set(vtkDataObject::ORIGIN(),
               this->ModelBounds[0],this->ModelBounds[2],
               this->ModelBounds[4]);
  memcpy(this->Origin,outInfo->Get(vtkDataObject::ORIGIN()), sizeof(double)*3);
  output->SetOrigin(this->Origin);

  for (i=0; i<3; i++)
    {
    this->Spacing[i] = (this->ModelBounds[2*i+1] - this->ModelBounds[2*i])
      / (this->SampleDimensions[i] - 1);
    if ( this->Spacing[i] <= 0.0 )
      {
      this->Spacing[i] = 1.0;
      }
    }
  outInfo->Set(vtkDataObject::SPACING(),this->Spacing,3);
  output->SetSpacing(this->Spacing);

  // Determine the splat propagation distance...used later
  for (i=0; i<3; i++)
    {
    this->SplatDistance[i] = maxDist / this->Spacing[i];
    }
}

// Set the dimensions of the sampling structured point set.
void vtkGaussianSplatterExtended::SetSampleDimensions(int i, int j, int k)
{
  int dim[3];

  dim[0] = i;
  dim[1] = j;
  dim[2] = k;

  this->SetSampleDimensions(dim);
}

//----------------------------------------------------------------------------
void vtkGaussianSplatterExtended::SetSampleDimensions(int dim[3])
{
  int dataDim, i;

  vtkDebugMacro(<< " setting SampleDimensions to (" << dim[0] << ","
                << dim[1] << "," << dim[2] << ")");

  if (dim[0] != this->SampleDimensions[0] ||
      dim[1] != this->SampleDimensions[1] ||
      dim[2] != this->SampleDimensions[2] )
    {
    if ( dim[0]<1 || dim[1]<1 || dim[2]<1 )
      {
      vtkErrorMacro (<< "Bad Sample Dimensions, retaining previous values");
      return;
      }

    for (dataDim=0, i=0; i<3 ; i++)
      {
      if (dim[i] > 1)
        {
        dataDim++;
        }
      }

    if ( dataDim  < 3 )
      {
      vtkErrorMacro(<<"Sample dimensions must define a volume!");
      return;
      }

    for ( i=0; i<3; i++)
      {
      this->SampleDimensions[i] = dim[i];
      }

    this->Modified();
    }
}

//----------------------------------------------------------------------------
void vtkGaussianSplatterExtended::Cap(vtkDoubleArray *s)
{
  int i,j,k;
  vtkIdType idx;
  int d01=this->SampleDimensions[0]*this->SampleDimensions[1];

  // i-j planes
  //k = 0;
  for (j=0; j<this->SampleDimensions[1]; j++)
    {
    for (i=0; i<this->SampleDimensions[0]; i++)
      {
      s->SetTuple(i+j*this->SampleDimensions[0], &this->CapValue);
      }
    }
  k = this->SampleDimensions[2] - 1;
  idx = k*d01;
  for (j=0; j<this->SampleDimensions[1]; j++)
    {
    for (i=0; i<this->SampleDimensions[0]; i++)
      {
      s->SetTuple(idx+i+j*this->SampleDimensions[0], &this->CapValue);
      }
    }
  // j-k planes
  //i = 0;
  for (k=0; k<this->SampleDimensions[2]; k++)
    {
    for (j=0; j<this->SampleDimensions[1]; j++)
      {
      s->SetTuple(j*this->SampleDimensions[0]+k*d01, &this->CapValue);
      }
    }
  i = this->SampleDimensions[0] - 1;
  for (k=0; k<this->SampleDimensions[2]; k++)
    {
    for (j=0; j<this->SampleDimensions[1]; j++)
      {
      s->SetTuple(i+j*this->SampleDimensions[0]+k*d01, &this->CapValue);
      }
    }
  // i-k planes
  //j = 0;
  for (k=0; k<this->SampleDimensions[2]; k++)
    {
    for (i=0; i<this->SampleDimensions[0]; i++)
      {
      s->SetTuple(i+k*d01, &this->CapValue);
      }
    }
  j = this->SampleDimensions[1] - 1;
  idx = j*this->SampleDimensions[0];
  for (k=0; k<this->SampleDimensions[2]; k++)
    {
    for (i=0; i<this->SampleDimensions[0]; i++)
      {
      s->SetTuple(idx+i+k*d01, &this->CapValue);
      }
    }
}

//----------------------------------------------------------------------------
//
//  Gaussian sampling
//
double vtkGaussianSplatterExtended::Gaussian (double cx[3])
{
  return ((cx[0]-P[0])*(cx[0]-P[0]) + (cx[1]-P[1])*(cx[1]-P[1]) +
          (cx[2]-P[2])*(cx[2]-P[2]) );
}

//----------------------------------------------------------------------------
//
//  Ellipsoidal Gaussian sampling
//
double vtkGaussianSplatterExtended::EccentricGaussian (double cx[3])
{
  double   v[3], r2, z2, rxy2, mag;

  v[0] = cx[0] - this->P[0];
  v[1] = cx[1] - this->P[1];
  v[2] = cx[2] - this->P[2];

  r2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];

  if ( (mag=this->N[0]*this->N[0]+
            this->N[1]*this->N[1]+
            this->N[2]*this->N[2]) != 1.0  )
    {
    if ( mag == 0.0 )
      {
      mag = 1.0;
      }
    else
      {
      mag = sqrt(mag);
      }
    }

  z2 = (v[0]*this->N[0] + v[1]*this->N[1] + v[2]*this->N[2])/mag;
  z2 = z2*z2;

  rxy2 = r2 - z2;

  return (rxy2/this->Eccentricity2 + z2);
}

//----------------------------------------------------------------------------
void vtkGaussianSplatterExtended::SetScalar(int idx, double dist2,
                                    vtkDoubleArray *newScalars)
{
  double v;
  if (this->Radius2 == 0)
      v = (this->*SampleFactor)(this->S);
  else
      v = (this->*SampleFactor)(this->S) * exp(
    static_cast<double>
    (this->ExponentFactor*(dist2)/(this->Radius2)));

  if ( ! this->Visited[idx] )
    {
    this->Visited[idx] = 1;
    newScalars->SetTuple(idx,&v);
    }
  else
    {
    double s = newScalars->GetValue(idx);
    switch (this->AccumulationMode)
      {
      case VTK_ACCUMULATION_MODE_MIN:
        newScalars->SetTuple(idx,(s < v ? &s : &v));
        break;
      case VTK_ACCUMULATION_MODE_MAX:
        newScalars->SetTuple(idx,(s > v ? &s : &v));
        break;
      case VTK_ACCUMULATION_MODE_SUM:
        s += v;
        newScalars->SetTuple(idx,&s);
        break;
      }
    }//not first visit
}

//----------------------------------------------------------------------------
const char *vtkGaussianSplatterExtended::GetAccumulationModeAsString()
{
  if ( this->AccumulationMode == VTK_ACCUMULATION_MODE_MIN )
    {
    return "Minimum";
    }
  else if ( this->AccumulationMode == VTK_ACCUMULATION_MODE_MAX )
    {
    return "Maximum";
    }
  else //if ( this->AccumulationMode == VTK_ACCUMULATION_MODE_SUM )
    {
    return "Sum";
    }
}

//----------------------------------------------------------------------------
void vtkGaussianSplatterExtended::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "Sample Dimensions: ("
               << this->SampleDimensions[0] << ", "
               << this->SampleDimensions[1] << ", "
               << this->SampleDimensions[2] << ")\n";

  os << indent << "Radius: " << this->Radius << "\n";
  os << indent << "Exponent Factor: " << this->ExponentFactor << "\n";

  os << indent << "ModelBounds: \n";
  os << indent << "  Xmin,Xmax: (" << this->ModelBounds[0]
     << ", " << this->ModelBounds[1] << ")\n";
  os << indent << "  Ymin,Ymax: (" << this->ModelBounds[2]
     << ", " << this->ModelBounds[3] << ")\n";
  os << indent << "  Zmin,Zmax: (" << this->ModelBounds[4]
     << ", " << this->ModelBounds[5] << ")\n";

  os << indent << "Normal Warping: "
     << (this->NormalWarping ? "On\n" : "Off\n");
  os << indent << "Eccentricity: " << this->Eccentricity << "\n";

  os << indent << "Scalar Warping: "
     << (this->ScalarWarping ? "On\n" : "Off\n");
  os << indent << "Scale Factor: " << this->ScaleFactor << "\n";

  os << indent << "Capping: " << (this->Capping ? "On\n" : "Off\n");
  os << indent << "Cap Value: " << this->CapValue << "\n";

  os << indent << "Accumulation Mode: "
     << this->GetAccumulationModeAsString() << "\n";

  os << indent << "Null Value: " << this->NullValue << "\n";
}

//----------------------------------------------------------------------------
int vtkGaussianSplatterExtended::FillInputPortInformation(
  int vtkNotUsed(port), vtkInformation* info)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
  return 1;
}
