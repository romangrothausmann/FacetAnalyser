
#include "FacetAnalyser.h"

#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"

#include <vtkSmartPointer.h>

#include <vtkMath.h>
#include <vtkPolyDataNormals.h>
#include <vtkMeshQuality.h>
#include <vtkGaussianSplatter.h>
#include <vtkImageCast.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkPolyData.h>
#include <vtkIdTypeArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkPlanes.h>
#include <vtkHull.h>

#include <itkVTKImageToImageFilter.h>
#include <itkShiftScaleImageFilter.h>
#include <itkHMinimaImageFilter.h>
#include <itkRegionalMinimaImageFilter.h>
#include <itkConnectedComponentImageFilter.h>
#include <itkMorphologicalWatershedFromMarkersImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkAddImageFilter.h>
#include <itkGradientMagnitudeImageFilter.h>
#include <itkChangeLabelImageFilter.h>
#include <itkMaskImageFilter.h>

#include <itkStatisticsLabelObject.h>
#include <itkLabelMap.h>
#include <itkLabelImageToStatisticsLabelMapFilter.h>

#define SMB 1.1 //splatter boundaries: make it bigger than 1 to catch the centres!
#define msigma 2 //render gauss upto msigma*sigma; 1~68%; 2~95%; 3~99%


vtkStandardNewMacro(FacetAnalyser);

//----------------------------------------------------------------------------
// Description:
FacetAnalyser::FacetAnalyser(){
    this->SetNumberOfOutputPorts(2);

    this->SampleSize= 101;
    this->AngleUncertainty= 10;
    this->MinTrianglesPerFacet= 10;
    }

//----------------------------------------------------------------------------

int FacetAnalyser::RequestData(
    vtkInformation *vtkNotUsed(request),
    vtkInformationVector **inputVector,
    vtkInformationVector *outputVector)
    {
    // get the info objects
    vtkInformation *inInfo0 = inputVector[0]->GetInformationObject(0);
    //vtkInformation *inInfo1 = inputVector[0]->GetInformationObject(1);
    vtkInformation *outInfo0 = outputVector->GetInformationObject(0);
    vtkInformation *outInfo1 = outputVector->GetInformationObject(1);

    // get the input and ouptut
    vtkPolyData *input = vtkPolyData::SafeDownCast(
        inInfo0->Get(vtkDataObject::DATA_OBJECT()));
    vtkPolyData *output0 = vtkPolyData::SafeDownCast(
        outInfo0->Get(vtkDataObject::DATA_OBJECT()));
    vtkPolyData *output1 = vtkPolyData::SafeDownCast(
        outInfo1->Get(vtkDataObject::DATA_OBJECT()));


    double da= this->AngleUncertainty / 180.0 * vtkMath::Pi(); 
    double f= 1/2./sin(da)/sin(da); //sin(da) corresponds to sigma
    double R= msigma / double(SMB) * sqrt(1/2./f);

    vtkSmartPointer<vtkPolyDataNormals> PDnormals0= vtkSmartPointer<vtkPolyDataNormals>::New();
    PDnormals0->SetInputData(input);
    PDnormals0->ComputePointNormalsOff(); 
    PDnormals0->ComputeCellNormalsOn();
    PDnormals0->Update();

    vtkSmartPointer<vtkMeshQuality> cellArea= vtkSmartPointer<vtkMeshQuality>::New();
    cellArea->SetInputConnection(PDnormals0->GetOutputPort());
    cellArea->SetTriangleQualityMeasureToArea();
    cellArea->SaveCellQualityOn();
    cellArea->Update();

    vtkDataArray *normals= PDnormals0->GetOutput()->GetCellData()->GetNormals();
    vtkDoubleArray *areas= vtkDoubleArray::SafeDownCast(cellArea->GetOutput()->GetCellData()->GetArray("Quality"));

    vtkIdType NumPolyDataPoints= normals->GetNumberOfTuples();

    ////vtkMeshQuality has no sum (only min, max, avg) so calc it here:
    double totalPolyDataArea= 0;
    for(vtkIdType k= 0; k < NumPolyDataPoints; k++)
        totalPolyDataArea+= areas->GetValue(k);

    ////regard normals coords as point coords
    vtkSmartPointer<vtkPoints> Points = vtkSmartPointer<vtkPoints>::New();
    for(vtkIdType i= 0; i < NumPolyDataPoints; i++)
        Points->InsertNextPoint(normals->GetTuple3(i));

    vtkSmartPointer<vtkPolyData> polydata0 = vtkSmartPointer<vtkPolyData>::New();
    polydata0->SetPoints(Points);
    polydata0->GetPointData()->SetScalars(areas);

    vtkSmartPointer<vtkGaussianSplatter> Splatter = vtkSmartPointer<vtkGaussianSplatter>::New();
    Splatter->SetInputData(polydata0);
    Splatter->SetSampleDimensions(this->SampleSize,this->SampleSize,this->SampleSize); //set the resolution of the final! volume
    Splatter->SetModelBounds(-SMB,SMB, -SMB,SMB, -SMB,SMB);
    Splatter->SetExponentFactor(-f); //GaussianSplat decay value
    Splatter->SetRadius(R); //GaussianSplat truncated outside Radius
    Splatter->SetAccumulationModeToSum();
    Splatter->ScalarWarpingOn(); //use individual weights
    Splatter->SetScaleFactor(1/totalPolyDataArea);
    Splatter->CappingOff(); //don't pad with 0
    Splatter->Update();
 
    vtkSmartPointer<vtkImageCast> cast = vtkSmartPointer<vtkImageCast>::New();
    cast->SetInputConnection(Splatter->GetOutputPort());
    cast->SetOutputScalarTypeToDouble();
    cast->Update(); //seems to be essential for vtk-6.1.0 + itk-4.5.1


/////////////////////////////////////////////////
///////////////// going to ITK //////////////////
/////////////////////////////////////////////////


    const int dim = 3;

    typedef double PixelType; //vtk 6.1.0 + itk 4.5.1 seem to use double instead of float
    typedef itk::Image< PixelType, dim >    ImageType;
    typedef ImageType GreyImageType;

    typedef itk::VTKImageToImageFilter<ImageType> ConnectorType;

    ConnectorType::Pointer vtkitkf = ConnectorType::New();
    vtkitkf->SetInput(cast->GetOutput()); //NOT GetOutputPort()!!!
    vtkitkf->Update();

/////////////////////// creating the sphere mask
//also outer radius to ensure proper probability interpretation!!!

    typedef  unsigned char MPixelType;

    const unsigned int Dimension = dim;

    typedef itk::Image<MPixelType,  Dimension>   MImageType;


    typedef itk::ShiftScaleImageFilter<GreyImageType, GreyImageType> SSType;
    SSType::Pointer ss = SSType::New();
    ss->SetScale(-1); //invert by mul. with -1
    ss->SetInput(vtkitkf->GetOutput());
    ss->Update();

    typedef unsigned short LabelType;
    typedef itk::Image<LabelType,  dim>   LImageType;

    //use ws from markers because the label image is needed twice
    bool ws1_conn= true;  //fully connected borders!
    bool ws2_conn= false; //finer borders??? and no lost labels! (see mws pdf)
    bool ws3_conn= ws2_conn; 

    typedef itk::HMinimaImageFilter<GreyImageType, GreyImageType> HMType; //seems for hmin in-type==out-type!!!
    HMType::Pointer hm= HMType::New();
    hm->SetHeight(this->MinTrianglesPerFacet / NumPolyDataPoints);
    hm->SetFullyConnected(ws1_conn);
    hm->SetInput(ss->GetOutput());

    typedef itk::RegionalMinimaImageFilter<GreyImageType, MImageType> RegMinType;
    RegMinType::Pointer rm = RegMinType::New();
    rm->SetFullyConnected(ws1_conn);
    rm->SetInput(hm->GetOutput());

    // connected component labelling
    typedef itk::ConnectedComponentImageFilter<MImageType, LImageType> CCType;
    CCType::Pointer labeller = CCType::New();
    labeller->SetFullyConnected(ws1_conn);
    labeller->SetInput(rm->GetOutput());

    typedef itk::MorphologicalWatershedFromMarkersImageFilter<GreyImageType, LImageType> MWatershedType;
    MWatershedType::Pointer ws1 = MWatershedType::New();
    ws1->SetMarkWatershedLine(true); //use borders as marker in sd. ws
    ws1->SetFullyConnected(ws1_conn); //true reduces amount of watersheds
    //ws1->SetLevel(1 / double(ndp) * atof(argv[4])); //if 0: hminima skipted?
    ws1->SetInput(ss->GetOutput());
    ws1->SetMarkerImage(labeller->GetOutput());
    ws1->Update(); //whith out this update label one is lost for facet_holger particle0195!!! Why???

    // extract the watershed lines and combine with the orginal markers
    typedef itk::BinaryThresholdImageFilter<LImageType, LImageType> ThreshType;
    ThreshType::Pointer th = ThreshType::New();
    th->SetUpperThreshold(0);
    th->SetOutsideValue(0);
    // set the inside value to the number of markers + 1
    th->SetInsideValue(labeller->GetObjectCount() + 1);
    th->SetInput(ws1->GetOutput());

    // Add the marker image to the watershed line image
    typedef itk::AddImageFilter<LImageType, LImageType, LImageType> AddType;
    AddType::Pointer adder1= AddType::New();
    adder1->SetInput1(th->GetOutput());
    adder1->SetInput2(labeller->GetOutput());

    // compute a gradient
    typedef itk::GradientMagnitudeImageFilter<GreyImageType, GreyImageType> GMType;
    GMType::Pointer gm1= GMType::New();
    gm1->SetInput(ss->GetOutput());

    // Now apply sd. watershed
    MWatershedType::Pointer ws2 = MWatershedType::New();
    ws2->SetMarkWatershedLine(false); //no use for a border in sd. stage
    ws2->SetFullyConnected(ws2_conn); 
    //ws->SetLevel(1 / double(ndp) * atof(argv[4])); //if 0: hminima skipted?
    ws2->SetInput(gm1->GetOutput());
    ws2->SetMarkerImage(adder1->GetOutput());

    // delete the background label
    typedef itk::ChangeLabelImageFilter<LImageType, LImageType> ChangeLabType;
    ChangeLabType::Pointer ch1= ChangeLabType::New();
    ch1->SetInput(ws2->GetOutput());
    ch1->SetChange(labeller->GetObjectCount() + 1, 0);

    // combine the markers again
    //this result is not the same as ws2 because the bg-label has grown!
    AddType::Pointer adder2 = AddType::New();
    adder2->SetInput1(th->GetOutput());
    adder2->SetInput2(ch1->GetOutput());

    GMType::Pointer gm2 = GMType::New();
    gm2->SetInput(gm1->GetOutput());

    MWatershedType::Pointer ws3 = MWatershedType::New();
    ws3->SetFullyConnected(ws3_conn);
    ws3->SetInput(gm2->GetOutput());
    ws3->SetMarkerImage(adder2->GetOutput());
    ws3->SetMarkWatershedLine(false);			

    // delete the background label again
    ChangeLabType::Pointer ch2= ChangeLabType::New();
    ch2->SetInput(ws3->GetOutput());
    ch2->SetChange(labeller->GetObjectCount() + 1, 0);

////////////////////////Now label and grow the facet reagions... done.


    ////spalter only single points with weights    
    vtkGaussianSplatter *Splatter2 = vtkGaussianSplatter::New();
    Splatter2->SetInputData(polydata0);
    Splatter2->SetSampleDimensions(this->SampleSize,this->SampleSize,this->SampleSize); //set the resolution of the final! volume
    Splatter2->SetModelBounds(-SMB,SMB, -SMB,SMB, -SMB,SMB);
    Splatter2->SetExponentFactor(-1); //GaussianSplat decay value
    Splatter2->SetRadius(0); //only splat single points
    Splatter2->SetAccumulationModeToSum();
    Splatter2->ScalarWarpingOn(); //use individual weights
    Splatter2->SetScaleFactor(1/totalPolyDataArea);
    Splatter2->CappingOff(); //don't pad with 0
    Splatter2->Update();

    vtkImageCast* cast2 = vtkImageCast::New();
    cast2->SetInputConnection(Splatter2->GetOutputPort());
    //cast->SetOutputScalarTypeToUnsignedChar();
    //cast2->SetOutputScalarTypeToFloat();
    cast2->SetOutputScalarTypeToDouble();
    cast2->Update(); //seems to be essential for vtk-6.1.0 + itk-4.5.1

    ConnectorType::Pointer vtkitkf2 = ConnectorType::New();
    vtkitkf2->SetInput(cast2->GetOutput()); //NOT GetOutputPort()!!!
    vtkitkf2->Update();

    typedef itk::MaskImageFilter<LImageType, GreyImageType, LImageType> MaskType2;
    MaskType2::Pointer mask2 = MaskType2::New();
    mask2->SetInput1(ch2->GetOutput());
    mask2->SetInput2(vtkitkf2->GetOutput()); //mask
    mask2->Update();


    /////////////create first output/////////////

    ///////////probe each splat point value
    vtkDoubleArray* splatValues= vtkDoubleArray::SafeDownCast(Splatter->GetOutput()->GetPointData()->GetScalars());

    vtkSmartPointer<vtkIdTypeArray> fId= vtkSmartPointer<vtkIdTypeArray>::New();
    fId->SetName("FacetIds");
    fId->SetNumberOfComponents(1);

    vtkSmartPointer<vtkDoubleArray> fPb= vtkSmartPointer<vtkDoubleArray>::New();
    fPb->SetName("FacetProbabilities");
    fPb->SetNumberOfComponents(1);

    for(vtkIdType k= 0; k < NumPolyDataPoints; k++){
        double pp[3];
        Points->GetPoint(k, pp); 
        vtkIdType pi[3];
        vtkIdType idx= Splatter->ProbePoint(pp, pi);

        if (idx < 0){
            vtkErrorMacro(<< "Prope Point: " << pp[0] << ";" << pp[1] << ";" << pp[2] << " not insied sample data");
            return VTK_ERROR;
            }

        ////get probability value
        double tv= splatValues->GetValue(idx);

        ////get label value
        LImageType::IndexType lidx;
        lidx[0]= pi[0];
        lidx[1]= pi[1];
        lidx[2]= pi[2];
        LabelType tl= mask2->GetOutput()->GetPixel(lidx);
        //LabelType tl= mask2->GetOutput()->GetPixel(dynamic_cast<LImageType::IndexType>(idx));//does not work

        fId->InsertNextValue(tl);
        fPb->InsertNextValue(tv);
        }

    // Copy original points and point data
    output0->CopyStructure(input);
    output0->GetPointData()->PassData(input->GetPointData());
    output0->GetCellData()->PassData(input->GetCellData());
    output0->GetCellData()->AddArray(fId);
    output0->GetCellData()->AddArray(fPb);



    /////////////create data for second output/////////////

    vtkSmartPointer<vtkDoubleArray> facetNormals= vtkSmartPointer<vtkDoubleArray>::New();
    facetNormals->SetNumberOfComponents(3);
    facetNormals->SetName ("facetNormals");

    vtkSmartPointer<vtkDoubleArray> relFacetSizes= vtkSmartPointer<vtkDoubleArray>::New();
    relFacetSizes->SetNumberOfComponents(1);
    relFacetSizes->SetName ("relFacetSize");

    vtkSmartPointer<vtkDoubleArray> absFacetSizes= vtkSmartPointer<vtkDoubleArray>::New();
    absFacetSizes->SetNumberOfComponents(1);
    absFacetSizes->SetName ("absFacetSize");

    vtkSmartPointer<vtkPoints> facetNormalPoints = vtkSmartPointer<vtkPoints>::New();
    typedef itk::StatisticsLabelObject< LabelType, dim > LabelObjectType;
    typedef itk::LabelMap< LabelObjectType > LabelMapType;

    typedef itk::LabelImageToStatisticsLabelMapFilter<LImageType, GreyImageType, LabelMapType> ConverterType;
    ConverterType::Pointer converter = ConverterType::New();
    converter->SetInput(mask2->GetOutput());
    converter->SetFeatureImage(vtkitkf2->GetOutput()); //this should be the single splat grey image to be exact!
    //converter->SetFullyConnected(true); //true: 26-connectivity; false: 6-connectivity
    converter->SetComputePerimeter(false);
    converter->Update();

    LabelMapType::Pointer labelMap = converter->GetOutput();
    const LabelObjectType * labelObject;
    for( unsigned int label=1; label<=labelMap->GetNumberOfLabelObjects(); label++ )
        {
        try{
            labelObject = labelMap->GetLabelObject( label );
            }
        catch( itk::ExceptionObject exp ){
            if (strstr(exp.GetDescription(), "No label object with label")){
                std::cerr << "Missing label: " << label << std::endl;
                continue;
                }
            }
  
        double c[3], fw;
        //// Since GetCenterOfGravity(), GetCentroid() are in 3D the centres can lie off the unit sphere!
        //// The inverse of the length of the centre vector is therefor a measure of how concentrated/dispersed the label is and therefor how destinct a facet is
        c[0]= labelObject->GetCenterOfGravity()[0];
        c[1]= labelObject->GetCenterOfGravity()[1];
        c[2]= labelObject->GetCenterOfGravity()[2];
        fw= labelObject->GetSum();
        facetNormals->InsertNextTuple(c);
        relFacetSizes->InsertNextValue(fw);
        absFacetSizes->InsertNextValue(fw * totalPolyDataArea);
        facetNormalPoints->InsertNextPoint(0,0,0);
        }


    /////////////create second output/////////////

    vtkSmartPointer<vtkPlanes> planes= vtkSmartPointer<vtkPlanes>::New();
    planes->SetPoints(facetNormalPoints);
    planes->SetNormals(facetNormals);

    //vtkSmartPointer<vtkPolyData> polydata1 = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkHull> hull= vtkSmartPointer<vtkHull>::New();
    hull->SetPlanes(planes);
    hull->SetInputData(input);
    //hull->GenerateHull(polydata1, -20, 20, -20, 20, -20, 20)
    hull->Update();

    output1->ShallowCopy(hull->GetOutput());
    output1->GetCellData()->SetNormals(facetNormals);
    output1->GetCellData()->AddArray(relFacetSizes);
    output1->GetCellData()->AddArray(absFacetSizes);

    static const std::string fn= "test";


    //unsigned int nfp= facetNormals->GetNumberOfPoints();
    vtkIdType nfp= facetNormals->GetNumberOfTuples();

    //fstream of;
    //of.open((fn + ".adat").data(), ios::out);
    //of << "#p0_x\tp0_y\tp0_z\tp1_x\tp1_y\tp1_z\ta" << endl;
    FILE *of;
    of= fopen ((fn + ".adat").data(),"w");
    fprintf(of, "#i1\ti2\tp0_x\tp0_y\tp0_z\tp1_x\tp1_y\tp1_z\tangle\ta_weight\n");

    if(nfp > 1){

        double angle, aw;
        double p0[3], p1[3];
        //unsigned int u,v;
        vtkIdType u,v;

        for(u= 0; u < nfp - 1 ; u++){
            for(v= u + 1; v < nfp; v++){
                facetNormals->GetTuple(u, p0); 
                facetNormals->GetTuple(v, p1);
                vtkMath::Normalize(p0);
                vtkMath::Normalize(p1);
                angle= acos(vtkMath::Dot(p0, p1)) * 180.0 / vtkMath::Pi();
                aw= 2 / (1/relFacetSizes->GetTuple1(u) + 1/relFacetSizes->GetTuple1(v)); //harmonic mean because its the "smallest" of the Pythagorean means: http://en.wikipedia.org/wiki/Pythagorean_means
                //cout << angle << "; ";
                fprintf(of, "%lld\t%lld\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", u+1 , v+1, p0[0], p0[1], p0[2], p1[0], p1[1], p1[2], angle, aw);
                }
            }

        //cout << endl;
        }
    //of.close();
    fclose(of);




    return 1;
    }

//----------------------------------------------------------------------------
void FacetAnalyser::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "SampleSize: " << this->SampleSize << endl;
  os << indent << "AngleUncertainty: " << this->AngleUncertainty << endl;
  os << indent << "MinTrianglesPerFacet: " << this->MinTrianglesPerFacet << endl;
  os << indent << endl;
}
