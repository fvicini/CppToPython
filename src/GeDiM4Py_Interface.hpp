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

    private:
      template<typename T>
      static std::vector<T> ConvertArray(PyObject* list);

    public:
      static InterfaceConfiguration ConvertConfiguration(PyObject* config);
      static Domain2D ConvertDomainSquare(InterfaceDataDAO& gedimData,
                                          PyObject* square);
  };
}

#endif
