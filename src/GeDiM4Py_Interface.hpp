#ifndef __GeDiM4Py_Interface_H
#define __GeDiM4Py_Interface_H

#include <Python.h>
#include <iostream>

#include "GeDiM4Py_Logic.hpp"

// ***************************************************************************
extern "C"
void GedimForPy_Initialize(PyObject* config);
extern "C"
void GedimForPy_CreateDomainSquare(PyObject* square);
extern "C"
PyObject* GedimForPy_Discretize(PyObject* discreteSpace,
                                double** dofsCoordinate,
                                double** strongsCoordinate);
extern "C"
void GedimForPy_AssembleStiffnessMatrix(GedimForPy::GeDiM4Py_Logic::K k,
                                        int* numTriplets,
                                        double** stiffnessTriplets);
extern "C"
void GedimForPy_AssembleForcingTerm(GedimForPy::GeDiM4Py_Logic::F f,
                                    int* size,
                                    double** forcingTerm);
extern "C"
void GedimForPy_CholeskySolver(const int nRows,
                               const int numNonZeros,
                               const double* pointerTriplets,
                               const double* pointerF,
                               double** solution);
// ***************************************************************************

namespace GedimForPy
{
  class GeDiM4Py_Interface final
  {
    public:
      static GedimForPy::InterfaceConfiguration InterfaceConfig;
      static GedimForPy::InterfaceData InterfaceData;
      static GedimForPy::Domain2D Domain;
      static GedimForPy::Domain2DMesh Mesh;
      static GedimForPy::DiscreteSpace Space;
      static GedimForPy::DiscreteProblemData ProblemData;

    private:
      template<typename T>
      static std::vector<T> ConvertToArray(PyObject* list);

    public:
      static InterfaceConfiguration ConvertConfiguration(PyObject* config);
      static Domain2D ConvertDomainSquare(InterfaceDataDAO& gedimData,
                                          PyObject* square);
      static DiscreteSpace ConvertDiscreteSpace(PyObject* discreteSpace);
      static PyObject* ConvertProblemData(GedimForPy::DiscreteProblemData& problemData,
                                          double*& dofsCoordinate,
                                          double*& strongsCoordinate);

      static void ConvertTriplets(const std::list<Eigen::Triplet<double>>& triplets,
                                  int& numTriplets,
                                  double*& convertedTriplets);
      static void ConvertMatrix(const Eigen::MatrixXd& matrix,
                                double*& convertedMatrix);
      static void ConvertArray(const Eigen::VectorXd& array,
                               int& size,
                               double*& convertedArray);
  };
}

#endif
