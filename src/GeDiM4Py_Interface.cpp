#include "GeDiM4Py_Interface.hpp"

// ***************************************************************************
void GedimForPy_Initialize(PyObject* config)
{
  if (!PyDict_Check(config))
  {
    std::cout<< "The input is not correct"<< std::endl;
    return;
  }

  PyObject* geometricTolerance = PyDict_GetItemString(config, "GeometricTolerance");

  std::cout<< "GeometricTolerance "<< PyFloat_AsDouble(geometricTolerance)<< std::endl;
}
// ***************************************************************************
