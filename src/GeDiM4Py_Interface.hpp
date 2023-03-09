#ifndef __GeDiM4Py_Interface_H
#define __GeDiM4Py_Interface_H

#include <Python.h>
#include <iostream>

#include "GeDiM4Py_Logic.hpp"

// ***************************************************************************
typedef const double* (*K_Py)(const int numPoints, const double* points);

extern "C"
void GedimForPy_Initialize(PyObject* config);
extern "C"
void GedimForPy_CreateDomainSquare(PyObject* square);
extern "C"
PyObject* GedimForPy_Discretize(PyObject* discreteSpace);
extern "C"
void GedimForPy_AssembleStiffnessMatrix(GedimForPy::GeDiM4Py_Logic::K k,
                                        double** stiffnessTriplets);
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
      static std::vector<T> ConvertArray(PyObject* list);

    public:
      static InterfaceConfiguration ConvertConfiguration(PyObject* config);
      static Domain2D ConvertDomainSquare(InterfaceDataDAO& gedimData,
                                          PyObject* square);
      static DiscreteSpace ConvertDiscreteSpace(PyObject* discreteSpace);
      static PyObject* ConvertProblemData(GedimForPy::DiscreteProblemData& problemData);

      static void ConvertTriplets(const std::list<Eigen::Triplet<double>>& triplets,
                                  double** convertedTriplets);
  };
}

#endif
