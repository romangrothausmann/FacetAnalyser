// .NAME FacetAnalyser - sample filter for paraview
// .SECTION Description
// This filter demonstrates importing filters to paraview


#ifndef __FacetAnalyser_h
#define __FacetAnalyser_h

#include "vtkPolyDataAlgorithm.h"
#include <vtkImplicitFunction.h>

class VTK_EXPORT FacetAnalyser : public vtkPolyDataAlgorithm
{
public:
  static FacetAnalyser *New();
  vtkTypeMacro(FacetAnalyser,vtkPolyDataAlgorithm);
  //vtkCxxRevisionMacro(FacetAnalyser, "$Revision$");//http://www.paraview.org/Wiki/ParaView/Plugin_HowTo#Undefined_symbol_ZTV12vtkYourFilter

  void PrintSelf(ostream& os, vtkIndent indent);

  vtkSetMacro(SampleSize, unsigned int);
  vtkGetMacro(SampleSize, unsigned int);
  vtkSetMacro(AngleUncertainty, double);
  vtkGetMacro(AngleUncertainty, double);
  /* vtkSetMacro(MinTrianglesPerFacet, unsigned int); */
  /* vtkGetMacro(MinTrianglesPerFacet, unsigned int); */
  vtkSetMacro(MinTrianglesPerFacet, double);
  vtkGetMacro(MinTrianglesPerFacet, double);

protected:
  FacetAnalyser();
  ~FacetAnalyser() {};

  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

  unsigned int SampleSize, MinTrianglesPerFacet;
  double AngleUncertainty;

private:
  FacetAnalyser(const FacetAnalyser&);  // Not implemented.
  void operator=(const FacetAnalyser&);  // Not implemented.

};

#endif


