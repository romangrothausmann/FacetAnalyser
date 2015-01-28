////test program for the vtk-itk-filter: FacetAnalyser

#include <vtkSmartPointer.h>

#include <vtkDataSetReader.h>
#include "FacetAnalyser.h"
#include <vtkFieldData.h>


int main(int argc, char* argv[]){

    if( argc != 2 )
	{
	std::cerr << "Usage: " << argv[0];
	std::cerr << " inputMesh";
	std::cerr << std::endl;  
	return EXIT_FAILURE;
	}

    if(!(strcasestr(argv[1],".vtk"))) {
	std::cerr << "The input should end with .vtk" << std::endl; 
	return EXIT_FAILURE;
	}

    unsigned int SampleSize= 51;
    double AngleUncertainty= 10;
    double SplatRadius= 0.2;
    double MinRelFacetSize= 0.001;
    unsigned int NumberOfExtraWS= 2;
    bool OuterHull= false;
    bool AreaWeight= true;

    //vtkObject::SetGlobalWarningDisplay(1);

    vtkSmartPointer<vtkDataSetReader> reader0 = vtkSmartPointer<vtkDataSetReader>::New();
 
    reader0->SetFileName("8faced-rhombic-dodecahedron_twinned_simplified.vtk");
    reader0->Update();

    vtkSmartPointer<FacetAnalyser> filter = vtkSmartPointer<FacetAnalyser>::New(); 

    filter->SetInputConnection(reader0->GetOutputPort());
    filter->SetSampleSize(SampleSize);
    filter->SetAngleUncertainty(AngleUncertainty);
    filter->SetSplatRadius(SplatRadius);
    filter->SetMinRelFacetSize(MinRelFacetSize);
    filter->SetNumberOfExtraWS(NumberOfExtraWS);
    filter->SetOuterHull(OuterHull);
    filter->SetAreaWeight(AreaWeight);

    filter->Update();


    ///////////// testing ///////////
    if(!filter->GetOutput()->GetNumberOfPoints())
        return EXIT_FAILURE;

    int N= filter->GetOutput(0)->GetFieldData()->GetArray("interplanarAngles")->GetNumberOfTuples();

    for(int i = 0; i < N; i++){
	// Test the existence of interplanarAngles
	std::cout << filter->GetOutput(0)->GetFieldData()->GetArray("interplanarAngles")->GetTuple1(i)
		  << " "
		  << std::endl;
	}

    return EXIT_SUCCESS;
 
    }


