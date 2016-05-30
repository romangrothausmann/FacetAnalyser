// .NAME FacetAnalyser - sample filter for paraview
// .SECTION Description
// This filter demonstrates importing filters to paraview


#ifndef __FacetAnalyser_h
#define __FacetAnalyser_h

#include <vtkPolyDataAlgorithm.h>

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
  vtkSetMacro(SplatRadius, double);
  vtkGetMacro(SplatRadius, double);
  vtkSetMacro(MinRelFacetSize, double);
  vtkGetMacro(MinRelFacetSize, double);
  vtkSetMacro(NumberOfExtraWS, unsigned char);
  vtkGetMacro(NumberOfExtraWS, unsigned char);
  vtkSetMacro(OuterHull, bool);
  vtkGetMacro(OuterHull, bool);
  vtkSetMacro(AreaWeight, bool);
  vtkGetMacro(AreaWeight, bool);

protected:
  FacetAnalyser();
  ~FacetAnalyser() {};

  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

  unsigned int SampleSize;
  double AngleUncertainty, SplatRadius, MinRelFacetSize;
  bool OuterHull, AreaWeight;
  unsigned char NumberOfExtraWS;

private:
  FacetAnalyser(const FacetAnalyser&);  // Not implemented.
  void operator=(const FacetAnalyser&);  // Not implemented.

  vtkIdType CantorPairing(vtkIdType x, vtkIdType y);
  // vtkIdType computeX(vtkIdType z);
  // vtkIdType computeY(vtkIdType z);

  vtkIdType findSharedPoints(vtkIdType* pts0, vtkIdType* pts1, vtkIdType npts0, vtkIdType npts1, vtkIdList* ptIds);
  vtkIdType ProbePoint(const double Origin[3], const double Spacing[3], const int SampleDimensions[3], const double *pp, vtkIdType *pi);

};

#endif


