#ifndef __GeDiM4Py_Interface_H
#define __GeDiM4Py_Interface_H

#include <Python.h>
#include <iostream>

#include "GeDiM4Py_Logic.hpp"

// ***************************************************************************
extern "C"
void GedimForPy_Initialize(PyObject* config);
// ***************************************************************************

namespace GedimForPy
{
  class GeDiM4Py_Interface final
  {
    public:
      static GedimForPy::GeDiM4Py_Logic_Configuration InterfaceConfig;
      static GedimForPy::InterfaceData InterfaceData;

    public:
      static GeDiM4Py_Logic_Configuration ConvertConfiguration(PyObject* config);
  };
}

#endif
