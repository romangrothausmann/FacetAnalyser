
#include "FacetAnalyser.h"

#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"

#include <vtkSmartPointer.h>

#include <vtkCallbackCommand.h>
#include <vtkCommand.h>

#include <vtkMath.h>
#include <vtkTriangleFilter.h> 
#include <vtkPolyDataNormals.h>
#include <vtkMeshQuality.h>
#include "vtkGaussianSplatterExtended.h"
#include <vtkImageCast.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkPolyData.h>
#include <vtkIdTypeArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellCenters.h>
#include <vtkPlanes.h>
#include <vtkHull.h>
#include <vtkCleanPolyData.h>
#include <vtkCellArray.h>

#include <itkCommand.h>

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
    this->SetNumberOfOutputPorts(3);

    this->SampleSize= 101;
    this->AngleUncertainty= 10;
    this->SplatRadius= 0;
    this->MinRelFacetSize= 0.001;
    this->NumberOfExtraWS= 2;
    this->OuterHull= 0;
    this->AreaWeight= 0;
    }

//----------------------------------------------------------------------------
void FilterEventHandlerVTK(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData){

    vtkAlgorithm *filter= static_cast<vtkAlgorithm*>(caller);

    switch(eventId){
	case vtkCommand::ProgressEvent:
	    fprintf(stderr, "\r%s progress: %5.1f%%", filter->GetClassName(), 100.0 * filter->GetProgress());//stderr is flushed directly
	    break;
	case vtkCommand::EndEvent:
	    std::cerr << std::endl << std::flush;   
	    break;
	}
    }

void FilterEventHandlerITK(itk::Object *caller, const itk::EventObject &event, void*){

    const itk::ProcessObject* filter = static_cast<const itk::ProcessObject*>(caller);

    if(itk::ProgressEvent().CheckEvent(&event))
	fprintf(stderr, "\r%s progress: %5.1f%%", filter->GetNameOfClass(), 100.0 * filter->GetProgress());//stderr is flushed directly
    else if(itk::EndEvent().CheckEvent(&event))
	std::cerr << std::endl << std::flush;   
    }

//----------------------------------------------------------------------------
int FacetAnalyser::RequestData(
    vtkInformation *vtkNotUsed(request),
    vtkInformationVector **inputVector,
    vtkInformationVector *outputVector)
    {
    // get the info objects
    vtkInformation *inInfo0 = inputVector[0]->GetInformationObject(0);
    vtkInformation *outInfo0 = outputVector->GetInformationObject(0);
    vtkInformation *outInfo1 = outputVector->GetInformationObject(1);
    vtkInformation *outInfo2 = outputVector->GetInformationObject(2);

    // get the input
    vtkSmartPointer<vtkPolyData> tinput= vtkSmartPointer<vtkPolyData>::New();

	{//scoped to ensure only tinput is used later on
	vtkPolyData *input = vtkPolyData::SafeDownCast(
	    inInfo0->Get(vtkDataObject::DATA_OBJECT()));

	////triangulation needed for vtkMeshQuality, which cannot compute the area for general polygons!
	vtkSmartPointer<vtkTriangleFilter> triangulate= vtkSmartPointer<vtkTriangleFilter>::New();
	triangulate->SetInputData(input);
	triangulate->PassLinesOff();
	triangulate->PassVertsOff();
	triangulate->Update();

	tinput->ShallowCopy(triangulate->GetOutput());
	triangulate->ReleaseDataFlagOn();
	}

    // get the ouptut
    vtkPolyData *output0 = vtkPolyData::SafeDownCast(
        outInfo0->Get(vtkDataObject::DATA_OBJECT()));
    vtkPolyData *output1 = vtkPolyData::SafeDownCast(
        outInfo1->Get(vtkDataObject::DATA_OBJECT()));
    vtkPolyData *output2 = vtkPolyData::SafeDownCast(
        outInfo2->Get(vtkDataObject::DATA_OBJECT()));

    vtkSmartPointer<vtkCallbackCommand> eventCallbackVTK = vtkSmartPointer<vtkCallbackCommand>::New();
    eventCallbackVTK->SetCallback(FilterEventHandlerVTK);

    this->UpdateProgress(0.0);

    double da= this->AngleUncertainty / 180.0 * vtkMath::Pi(); 
    double f= 1/2./sin(da)/sin(da); //sin(da) corresponds to sigma

    double R;
    if(this->SplatRadius)
	R= this->SplatRadius;
    else
	R= msigma / double(SMB) * sqrt(1/2./f);

    vtkSmartPointer<vtkPolyDataNormals> PDnormals0= vtkSmartPointer<vtkPolyDataNormals>::New();
    PDnormals0->SetInputData(tinput);
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

    this->UpdateProgress(0.1);

    vtkSmartPointer<vtkGaussianSplatterExtended> Splatter = vtkSmartPointer<vtkGaussianSplatterExtended>::New();
    Splatter->SetInputData(polydata0);
    Splatter->SetSampleDimensions(this->SampleSize,this->SampleSize,this->SampleSize); //set the resolution of the final! volume
    Splatter->SetModelBounds(-SMB,SMB, -SMB,SMB, -SMB,SMB);
    Splatter->SetExponentFactor(-f); //GaussianSplat decay value
    Splatter->SetRadius(R); //GaussianSplat truncated outside Radius
    Splatter->SetAccumulationModeToSum();
    Splatter->ScalarWarpingOn(); //use individual weights
    Splatter->SetScaleFactor(1/totalPolyDataArea);
    Splatter->CappingOff(); //don't pad with 0
    Splatter->AddObserver(vtkCommand::ProgressEvent, eventCallbackVTK);
    Splatter->AddObserver(vtkCommand::EndEvent, eventCallbackVTK);
    Splatter->Update();
 
    this->UpdateProgress(0.2);

    vtkSmartPointer<vtkImageCast> cast = vtkSmartPointer<vtkImageCast>::New();
    cast->SetInputConnection(Splatter->GetOutputPort());
    cast->SetOutputScalarTypeToDouble();
    cast->Update(); //seems to be essential for vtk-6.1.0 + itk-4.5.1


/////////////////////////////////////////////////
///////////////// going to ITK //////////////////
/////////////////////////////////////////////////


    const unsigned int Dimension = 3;

    bool ws0_conn= true;//[[[doc difference]]]//true reduces amount of watersheds
    bool ws_conn= false;//[[[doc difference]]]

    typedef double         GreyType;
    typedef unsigned short LabelType;
    typedef unsigned char  MaskType;

    typedef itk::Image<GreyType,  Dimension>   GreyImageType;
    typedef itk::Image<LabelType, Dimension>   LabelImageType;
    typedef itk::Image<MaskType,  Dimension>   MaskImageType;

    itk::CStyleCommand::Pointer eventCallbackITK;
    eventCallbackITK = itk::CStyleCommand::New();
    eventCallbackITK->SetCallback(FilterEventHandlerITK);


    typedef itk::VTKImageToImageFilter<GreyImageType> ConnectorType;
    ConnectorType::Pointer vtkitkf = ConnectorType::New();
    vtkitkf->SetInput(cast->GetOutput()); //NOT GetOutputPort()!!!
    vtkitkf->Update();

    this->UpdateProgress(0.3);

    typedef itk::ShiftScaleImageFilter<GreyImageType, GreyImageType> SSType;
    SSType::Pointer ss = SSType::New();
    ss->SetScale(-1); //invert by mul. with -1
    ss->SetInput(vtkitkf->GetOutput());
    ss->Update();

    typedef itk::HMinimaImageFilter<GreyImageType, GreyImageType> HMType; //seems for hmin in-type==out-type!!!
    HMType::Pointer hm= HMType::New();
    hm->SetHeight(this->MinRelFacetSize);
    hm->SetFullyConnected(ws0_conn);
    hm->SetInput(ss->GetOutput());
    hm->Update();

    typedef itk::RegionalMinimaImageFilter<GreyImageType, MaskImageType> RegMinType;
    RegMinType::Pointer rm = RegMinType::New();
    rm->SetFullyConnected(ws0_conn);
    rm->SetInput(hm->GetOutput());
    rm->Update();

    // connected component labelling
    typedef itk::ConnectedComponentImageFilter<MaskImageType, LabelImageType> CCType;
    CCType::Pointer labeller = CCType::New();
    labeller->SetFullyConnected(ws0_conn);
    labeller->SetInput(rm->GetOutput());
    labeller->Update();

    this->UpdateProgress(0.4);

    typedef itk::MorphologicalWatershedFromMarkersImageFilter<GreyImageType, LabelImageType> MWatershedType;
    MWatershedType::Pointer ws = MWatershedType::New();
    ws->SetMarkWatershedLine(this->NumberOfExtraWS); //use borders if higher order WS are wanted
    ws->SetFullyConnected(ws0_conn);
    ws->SetInput(ss->GetOutput());
    ws->SetMarkerImage(labeller->GetOutput());
    ws->AddObserver(itk::ProgressEvent(), eventCallbackITK);
    ws->AddObserver(itk::EndEvent(), eventCallbackITK);
    ws->Update(); 

    // extract the watershed lines and combine with the orginal markers
    typedef itk::BinaryThresholdImageFilter<LabelImageType, LabelImageType> ThreshType;
    ThreshType::Pointer th = ThreshType::New();
    th->SetUpperThreshold(0);
    th->SetOutsideValue(0);
    // set the inside value to the number of markers + 1
    th->SetInsideValue(labeller->GetObjectCount() + 1);
    th->SetInput(ws->GetOutput());
    th->Update();

    // to combine the markers again
    typedef itk::AddImageFilter<LabelImageType, LabelImageType, LabelImageType> AddType;
    AddType::Pointer adder = AddType::New();

    // to create gradient magnitude image
    typedef itk::GradientMagnitudeImageFilter<GreyImageType, GreyImageType> GMType;
    GMType::Pointer gm = GMType::New();

    // Add the marker image to the watershed line image
    adder->SetInput1(th->GetOutput());
    adder->SetInput2(labeller->GetOutput());
    adder->Update();

    // compute a gradient
    gm->SetInput(ss->GetOutput());
    gm->Update();

    this->UpdateProgress(0.5);

    LabelImageType::Pointer markerImg= adder->GetOutput();
    markerImg->DisconnectPipeline();
    GreyImageType::Pointer gradientImg= gm->GetOutput();
    gradientImg->DisconnectPipeline();
    LabelImageType::Pointer finalLabelImg= ws->GetOutput();
    finalLabelImg->DisconnectPipeline();

    ws->SetMarkWatershedLine(false); //no use for a border in higher stages
    ws->SetFullyConnected(ws_conn); 

    // to delete the background label
    typedef itk::ChangeLabelImageFilter<LabelImageType, LabelImageType> ChangeLabType;
    ChangeLabType::Pointer ch= ChangeLabType::New();
    ch->SetChange(labeller->GetObjectCount() + 1, 0);


    for(char i= 0; i < this->NumberOfExtraWS; i++){

	// Now apply higher order watershed
	ws->SetInput(gradientImg);
	ws->SetMarkerImage(markerImg);
	ws->Update();

	// delete the background label
	ch->SetInput(ws->GetOutput());
	ch->Update();

	// combine the markers again
	adder->SetInput1(th->GetOutput());
	adder->SetInput2(ch->GetOutput());
	adder->Update();

	gm->SetInput(gradientImg);
	gm->Update();

	markerImg= adder->GetOutput();
	markerImg->DisconnectPipeline();
	gradientImg= gm->GetOutput();
	gradientImg->DisconnectPipeline();
	finalLabelImg= ch->GetOutput();
	finalLabelImg->DisconnectPipeline();
	}

    this->UpdateProgress(0.7);

////////////////////////Now label and grow the facet reagions... done.


    ////spalter only single points with weights    
    vtkGaussianSplatterExtended *Splatter2 = vtkGaussianSplatterExtended::New();
    Splatter2->SetInputData(polydata0);
    Splatter2->SetSampleDimensions(this->SampleSize,this->SampleSize,this->SampleSize); //set the resolution of the final! volume
    Splatter2->SetModelBounds(-SMB,SMB, -SMB,SMB, -SMB,SMB);
    //Splatter2->SetRadius(1.0/(this->SampleSize+1)); //only splat single points (nearly, as a value of 0 does not work correctly), works with standard vtkGaussianSplatter
    Splatter2->SetRadius(0); //only splat single voxels for each point, needs vtkGaussianSplatterExtended
    Splatter2->SetExponentFactor(0); //no decay
    Splatter2->SetAccumulationModeToSum();
    Splatter2->ScalarWarpingOn(); //use individual weights
    Splatter2->SetScaleFactor(1/totalPolyDataArea);
    Splatter2->CappingOff(); //don't pad with 0
    Splatter2->AddObserver(vtkCommand::ProgressEvent, eventCallbackVTK);
    Splatter2->AddObserver(vtkCommand::EndEvent, eventCallbackVTK);
    Splatter2->Update();

    this->UpdateProgress(0.80);

    vtkImageCast* cast2 = vtkImageCast::New();
    cast2->SetInputConnection(Splatter2->GetOutputPort());
    //cast->SetOutputScalarTypeToUnsignedChar();
    //cast2->SetOutputScalarTypeToFloat();
    cast2->SetOutputScalarTypeToDouble();
    cast2->Update(); //seems to be essential for vtk-6.1.0 + itk-4.5.1

    ConnectorType::Pointer vtkitkf2 = ConnectorType::New();
    vtkitkf2->SetInput(cast2->GetOutput()); //NOT GetOutputPort()!!!
    vtkitkf2->Update();

    typedef itk::MaskImageFilter<LabelImageType, GreyImageType, LabelImageType> MaskFilterType;
    MaskFilterType::Pointer mask = MaskFilterType::New();
    mask->SetInput1(finalLabelImg);
    mask->SetInput2(vtkitkf->GetOutput()); //mask
    mask->Update();

    typedef itk::StatisticsLabelObject<LabelType, Dimension> LabelObjectType;
    typedef itk::LabelMap<LabelObjectType> LabelMapType;
    typedef itk::LabelImageToStatisticsLabelMapFilter<LabelImageType, GreyImageType, LabelMapType> ConverterType;
    ConverterType::Pointer converter = ConverterType::New();
    converter->SetInput(mask->GetOutput());
    converter->SetFeatureImage(vtkitkf2->GetOutput()); //this should be the single splat grey image to be exact!
    //converter->SetFullyConnected(true); //true: 26-connectivity; false: 6-connectivity
    converter->SetComputePerimeter(false);
    converter->Update();

    this->UpdateProgress(0.85);

    LabelMapType::Pointer labelMap = converter->GetOutput();
    vtkIdType NumFacets= labelMap->GetNumberOfLabelObjects();//needed later on

    /////////////ITK work done/////////////


    /////////////create first output/////////////

    ///////////probe each splat point value
    vtkDoubleArray* splatValues= vtkDoubleArray::SafeDownCast(Splatter->GetOutput()->GetPointData()->GetScalars());

    vtkSmartPointer<vtkIdTypeArray> fId= vtkSmartPointer<vtkIdTypeArray>::New();
    fId->SetName("FacetIds");
    fId->SetNumberOfComponents(1);

    vtkSmartPointer<vtkDoubleArray> fPb= vtkSmartPointer<vtkDoubleArray>::New();
    fPb->SetName("FacetProbabilities");
    fPb->SetNumberOfComponents(1);

    vtkSmartPointer<vtkPoints> facetCenterPoints= vtkSmartPointer<vtkPoints>::New();
    facetCenterPoints->SetNumberOfPoints(NumFacets);//label 0 is for unfacetted regions, not counted
    for(vtkIdType i= 0; i < NumFacets; i++) facetCenterPoints->SetPoint(i, 0, 0, 0);

    vtkSmartPointer<vtkDoubleArray> facetCentersCounter= vtkSmartPointer<vtkDoubleArray>::New();
    facetCentersCounter->SetNumberOfTuples(NumFacets);//label 0 is for unfacetted regions, not counted
    facetCentersCounter->FillComponent(0, 0);

    vtkSmartPointer<vtkCellCenters> cellCenters= vtkSmartPointer<vtkCellCenters>::New();
    cellCenters->SetInputData(tinput);
    cellCenters->Update(); //The cell attributes will be associated with the points on output.

    for(vtkIdType k= 0; k < NumPolyDataPoints; k++){
        double pp[3];
        Points->GetPoint(k, pp); 
        vtkIdType pi[3];
        vtkIdType idx= ProbePoint(Splatter->GetOutput()->GetOrigin(), Splatter->GetOutput()->GetSpacing(), Splatter->GetSampleDimensions(), pp, pi);

        if (idx < 0){
            vtkErrorMacro(<< "Probe Point: " << pp[0] << ";" << pp[1] << ";" << pp[2] << " not insied sample data");
            return VTK_ERROR;
            }

        ////get probability value
        double tv= splatValues->GetValue(idx);

        ////get label value
        LabelImageType::IndexType lidx;
        lidx[0]= pi[0];
        lidx[1]= pi[1];
        lidx[2]= pi[2];
        LabelType tl= mask->GetOutput()->GetPixel(lidx);

        fId->InsertNextValue(tl);
        fPb->InsertNextValue(tv);

	
	vtkIdType fl= tl - 1;
	if(fl >= 0){
	    double cp[3], fp[3];
	    double w= 1;
	    if(this->AreaWeight)
		w= areas->GetValue(k);
	    cellCenters->GetOutput()->GetPoint(k, cp);
	    facetCenterPoints->GetPoint(fl, fp);
	    facetCenterPoints->SetPoint(fl, fp[0]+cp[0]*w, fp[1]+cp[1]*w, fp[2]+cp[2]*w);
	    facetCentersCounter->SetValue(fl, facetCentersCounter->GetValue(fl)+w);//facetCentersCounter[k]+= w;
	    }
        }

    // Copy original points and point data
    output0->CopyStructure(tinput);//needs to be the polydata the normals were computed for!
    output0->GetPointData()->PassData(tinput->GetPointData());
    output0->GetCellData()->PassData(tinput->GetCellData());
    output0->GetCellData()->AddArray(fId);
    output0->GetCellData()->AddArray(fPb);



    /////////////create data for second output/////////////

    vtkSmartPointer<vtkDoubleArray> facetNormals= vtkSmartPointer<vtkDoubleArray>::New();
    facetNormals->SetNumberOfComponents(3);
    facetNormals->SetName ("facetNormals");

    vtkSmartPointer<vtkIdTypeArray> hullFacetIds= vtkSmartPointer<vtkIdTypeArray>::New();
    hullFacetIds->SetNumberOfComponents(1);
    hullFacetIds->SetName ("FacetIds");

    vtkSmartPointer<vtkDoubleArray> relFacetSizes= vtkSmartPointer<vtkDoubleArray>::New();
    relFacetSizes->SetNumberOfComponents(1);
    relFacetSizes->SetName ("relFacetSize");

    vtkSmartPointer<vtkDoubleArray> absFacetSizes= vtkSmartPointer<vtkDoubleArray>::New();
    absFacetSizes->SetNumberOfComponents(1);
    absFacetSizes->SetName ("absFacetSize");

    for(unsigned int label= 1; label <= NumFacets; label++){//skipping bg label 0, ie the "unfacetted" regions
        const LabelObjectType* labelObject;
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
        hullFacetIds->InsertNextValue(label);
        relFacetSizes->InsertNextValue(fw);
        absFacetSizes->InsertNextValue(fw * totalPolyDataArea);

	double fp[3];
	facetCenterPoints->GetPoint(label-1, fp);
	facetCenterPoints->SetPoint(label-1, fp[0]/facetCentersCounter->GetValue(label-1), fp[1]/facetCentersCounter->GetValue(label-1), fp[2]/facetCentersCounter->GetValue(label-1));
        }


    this->UpdateProgress(0.90);

    /////////////create second output/////////////

    vtkSmartPointer<vtkPlanes> planes= vtkSmartPointer<vtkPlanes>::New();
    planes->SetPoints(facetCenterPoints);
    planes->SetNormals(facetNormals);

    vtkSmartPointer<vtkCleanPolyData> cleanFilter= vtkSmartPointer<vtkCleanPolyData>::New();
    vtkSmartPointer<vtkHull> hull= vtkSmartPointer<vtkHull>::New();

    hull->SetPlanes(planes);
    if(this->OuterHull){
	hull->SetInputData(tinput);
	hull->Update();
	cleanFilter->SetInputConnection(hull->GetOutputPort());
	}
    else {
	vtkSmartPointer<vtkPolyData> polydata1 = vtkSmartPointer<vtkPolyData>::New();
	hull->GenerateHull(polydata1, tinput->GetBounds());//replaced SetInputData and Update
	cleanFilter->SetInputData(polydata1);
	}

    cleanFilter->PointMergingOn();//this is why it's done
    cleanFilter->ConvertPolysToLinesOn();
    cleanFilter->SetTolerance(0.000001);//small tolerance needed for most hulls
    cleanFilter->ToleranceIsAbsoluteOff();//relative tolerance 
    cleanFilter->Update();

    vtkDataArray* facetCenters = facetCenterPoints->GetData();
    facetCenters->SetName("FacetCenters");

    output1->ShallowCopy(cleanFilter->GetOutput());
    output1->GetCellData()->SetNormals(facetNormals);
    output1->GetCellData()->AddArray(facetCenters);
    output1->GetCellData()->AddArray(relFacetSizes);
    output1->GetCellData()->AddArray(hullFacetIds);
    output1->GetCellData()->AddArray(absFacetSizes);

    ////some of the planes set as input for vtkHull can get lost
    ////so set face-analyses also as FieldData to output0
    output0->GetFieldData()->AddArray(facetNormals);
    output0->GetFieldData()->AddArray(facetCenters);
    output0->GetFieldData()->AddArray(hullFacetIds);
    output0->GetFieldData()->AddArray(relFacetSizes);
    output0->GetFieldData()->AddArray(absFacetSizes);

    this->UpdateProgress(0.95);

    /////////////create second output field data/////////////
    
    vtkSmartPointer<vtkIdTypeArray> cellPairingIds= vtkSmartPointer<vtkIdTypeArray>::New();
    cellPairingIds->SetNumberOfComponents(1);
    cellPairingIds->SetName ("cellPairingIds");

    vtkSmartPointer<vtkDoubleArray> interplanarAngles= vtkSmartPointer<vtkDoubleArray>::New();
    interplanarAngles->SetNumberOfComponents(1);
    interplanarAngles->SetName ("interplanarAngles");

    vtkSmartPointer<vtkDoubleArray> angleWeights= vtkSmartPointer<vtkDoubleArray>::New();
    angleWeights->SetNumberOfComponents(1);
    angleWeights->SetName ("angleWeights");


    /////////////additionaly, create third output with cell data/////////////

    vtkSmartPointer<vtkCellArray> lines= vtkSmartPointer<vtkCellArray>::New();

    vtkSmartPointer<vtkIdTypeArray> lcellPairingIds= vtkSmartPointer<vtkIdTypeArray>::New();
    lcellPairingIds->SetNumberOfComponents(1);
    lcellPairingIds->SetName ("cellPairingIds");

    vtkSmartPointer<vtkDoubleArray> linterplanarAngles= vtkSmartPointer<vtkDoubleArray>::New();
    linterplanarAngles->SetNumberOfComponents(1);
    linterplanarAngles->SetName ("interplanarAngles");

    vtkSmartPointer<vtkDoubleArray> langleWeights= vtkSmartPointer<vtkDoubleArray>::New();
    langleWeights->SetNumberOfComponents(1);
    langleWeights->SetName ("angleWeights");

    if(NumFacets > 1){
        output1->BuildCells();
        for(vtkIdType u= 0; u < NumFacets - 1 ; u++){
            for(vtkIdType v= u + 1; v < NumFacets; v++){
                double p0[3], p1[3];
                facetNormals->GetTuple(u, p0); 
                facetNormals->GetTuple(v, p1);
                vtkMath::Normalize(p0);
                vtkMath::Normalize(p1);
                double angle= acos(vtkMath::Dot(p0, p1)) * 180.0 / vtkMath::Pi();
                double aw= 2 / (1/relFacetSizes->GetTuple1(u) + 1/relFacetSizes->GetTuple1(v)); 
 
                cellPairingIds->InsertNextValue(CantorPairing(u+1,v+1));//for consistency: label 0 is for unfacetted regions in output0
                interplanarAngles->InsertNextValue(angle);
                angleWeights->InsertNextValue(aw);

                vtkIdType *pts0, *pts1;
                vtkIdType npts0, npts1, nosc;
                output1->GetCellPoints(u,npts0,pts0);
                output1->GetCellPoints(v,npts1,pts1);
                vtkSmartPointer<vtkIdList> idlst= vtkSmartPointer<vtkIdList>::New();
                if((nosc= findSharedPoints(pts0, pts1, npts0, npts1, idlst)) == 2){
                    //vtkSmartPointer<vtkLine> line= vtkSmartPointer<vtkLine>::New();
                    lines->InsertNextCell(2);
                    lines->InsertCellPoint(idlst->GetId(0));
                    lines->InsertCellPoint(idlst->GetId(1));

                    lcellPairingIds->InsertNextValue(CantorPairing(u+1,v+1));//for consistency: label 0 is for unfacetted regions in output0
                    linterplanarAngles->InsertNextValue(angle);
                    langleWeights->InsertNextValue(aw);
                    }
                else if(nosc>2){
                    vtkErrorMacro(<< "Adjacent cells share more than one edge. This is not handled! " << u << "; " << v << "; " << "; " << npts0 << "; " << npts1 << "; " << nosc);
                    //return VTK_ERROR; //not needed, mesh will just lack some lines
                    }
                }
            }
        }

    output1->GetFieldData()->AddArray(cellPairingIds);
    output1->GetFieldData()->AddArray(interplanarAngles);
    output1->GetFieldData()->AddArray(angleWeights);

    ////as output0 already containds fdat as FieldData also add adat
    ////so output0 contains complete analyses
    output0->GetFieldData()->AddArray(cellPairingIds);
    output0->GetFieldData()->AddArray(interplanarAngles);
    output0->GetFieldData()->AddArray(angleWeights);

    output2->SetPoints(output1->GetPoints());
    output2->SetLines(lines);
 
    output2->GetCellData()->AddArray(lcellPairingIds);
    output2->GetCellData()->AddArray(linterplanarAngles);
    output2->GetCellData()->AddArray(langleWeights);

    this->UpdateProgress(1);

    return 1;
    }

//----------------------------------------------------------------------------
void FacetAnalyser::PrintSelf(ostream& os, vtkIndent indent){
    this->Superclass::PrintSelf(os,indent);

    os << indent << "SampleSize: " << this->SampleSize << endl;
    os << indent << "AngleUncertainty: " << this->AngleUncertainty << endl;
    os << indent << "MinRelFacetSize: " << this->MinRelFacetSize << endl;
    os << indent << endl;
    }


////Cantor pairing function
////http://de.wikipedia.org/wiki/Cantorsche_Paarungsfunktion#Implementierung_der_Berechnungen_in_Java

vtkIdType FacetAnalyser::CantorPairing(vtkIdType x, vtkIdType y){
    return (x+y)*(x+y+1)/2 + y;
    }

// vtkIdType FacetAnalyser::computeX(vtkIdType z){
//     vtkIdType j  = (vtkIdType) floor(sqrt(0.25 + 2*z) - 0.5);
//     return j - (z - j*(j+1)/2);
//     }

// vtkIdType FacetAnalyser::computeY(vtkIdType z){
//     vtkIdType j  = (vtkIdType) floor(sqrt(0.25 + 2*z) - 0.5);
//     return z - j*(j+1)/2;
//     }


vtkIdType FacetAnalyser::findSharedPoints(vtkIdType* pts0, vtkIdType* pts1, vtkIdType npts0, vtkIdType npts1, vtkIdList* ptIds){
    
    for (vtkIdType i= 0; i<npts0; i++)
        for (vtkIdType j= 0; j<npts1; j++)
            if(pts0[i] == pts1[j])
                ptIds->InsertNextId(pts0[i]);
    return(ptIds->GetNumberOfIds());
    }

vtkIdType FacetAnalyser::ProbePoint(const double Origin[3], const double Spacing[3], const int SampleDimensions[3], const double *pp, vtkIdType *pi){
    int i;
    double x[3];

    // Determine the voxel that the point is in
    for (i=0; i<3; i++)
        {
        x[i] = round((pp[i] - Origin[i]) / Spacing[i]); //round(), int(), ceil() or floor()???
        //cout << x[i] << endl;
        //check if point is inside output data
        if ( x[i] < 0 || x[i] >= SampleDimensions[i] )
            {
            return(-1);
            }
        pi[i]= x[i];
        }

    return(x[0] + x[1]*SampleDimensions[0] + x[2]*SampleDimensions[0]*SampleDimensions[1]);
    }
