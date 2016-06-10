////example program for the vtk-itk-filter: FacetAnalyser

#include <vtkSmartPointer.h>

#include <vtkCallbackCommand.h>
#include <vtkCommand.h>

#include <vtkXMLPolyDataReader.h>
#include "FacetAnalyser.h"
#include <vtkXMLPolyDataWriter.h>



#define P_VERBOSE 1

#define VTK_CREATE(type, name) vtkSmartPointer<type> name = vtkSmartPointer<type>::New()


void FilterEventHandler(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData){

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

int main(int argc, char* argv[]){


    if( argc < 8 )
	{
	std::cerr << "Usage: " << argv[0];
	std::cerr << " inputMesh";
	std::cerr << " SampleSize AngleUncertainty SplatRadius MinRelFacetSize NumberOfExtraWS";
	std::cerr << " outputMesh0 [outputMesh1] [outputMesh2] ";
	std::cerr << std::endl;  
	return EXIT_FAILURE;
	}

    if(!(strcasestr(argv[1],".vtp"))) {
	std::cerr << "The input should end with .vtp" << std::endl; 
	return -1;
	}

    if(!(strcasestr(argv[7],".vtp"))) {
	std::cerr << "The output should end with .vtp" << std::endl; 
	return -1;
	}

    if(argc >= 9)
	if(!(strcasestr(argv[8],".vtp"))) {
	    std::cerr << "The output should end with .vtp" << std::endl; 
	    return -1;
	    }

    if(argc >= 10)
	if(!(strcasestr(argv[9],".vtp"))) {
	    std::cerr << "The output should end with .vtp" << std::endl; 
	    return -1;
	    }


    vtkSmartPointer<vtkCallbackCommand> eventCallbackVTK = vtkSmartPointer<vtkCallbackCommand>::New();
    eventCallbackVTK->SetCallback(FilterEventHandler);


    unsigned int SampleSize= atoi(argv[2]);
    double AngleUncertainty= atof(argv[3]);
    double SplatRadius= atof(argv[4]);
    double MinRelFacetSize= atof(argv[5]);
    unsigned char NumberOfExtraWS= atoi(argv[6]);

    //vtkObject::SetGlobalWarningDisplay(1);

    vtkSmartPointer<vtkXMLPolyDataReader> reader0 = vtkSmartPointer<vtkXMLPolyDataReader>::New();
 
    reader0->SetFileName(argv[1]);
    reader0->Update();

    vtkSmartPointer<FacetAnalyser> filter = vtkSmartPointer<FacetAnalyser>::New(); 

    filter->SetInputConnection(reader0->GetOutputPort());
    filter->SetSampleSize(SampleSize);
    filter->SetAngleUncertainty(AngleUncertainty);
    filter->SetSplatRadius(SplatRadius);
    filter->SetMinRelFacetSize(MinRelFacetSize);
    filter->SetNumberOfExtraWS(NumberOfExtraWS);
    if(P_VERBOSE) filter->AddObserver(vtkCommand::ProgressEvent, eventCallbackVTK);
    if(P_VERBOSE) filter->AddObserver(vtkCommand::EndEvent, eventCallbackVTK);
    filter->Update();


    ///////////// testing ///////////
    if(!filter->GetOutput()->GetNumberOfPoints())
        return EXIT_FAILURE;


    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(argv[7]);
    writer->SetInputConnection(filter->GetOutputPort(0));
    writer->Update();

    if(argc >= 9){
	writer->SetFileName(argv[8]);
	writer->SetInputConnection(filter->GetOutputPort(1));
	writer->Update();
	}
    
    if(argc >= 10){
	writer->SetFileName(argv[9]);
	writer->SetInputConnection(filter->GetOutputPort(2));
	writer->Update();
	}
    
    return EXIT_SUCCESS;
 
    }


