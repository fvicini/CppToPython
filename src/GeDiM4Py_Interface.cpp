#include "GeDiM4Py_Interface.hpp"

// ***************************************************************************
void GedimForPy_Initialize(PyObject* config)
{
  if (!PyDict_Check(config))
  {
    std::cout<< "The input is not correct"<< std::endl;
    return;
  }

  GedimForPy::GeDiM4Py_Logic_Configuration& configuration = GedimForPy::GeDiM4Py_Interface::InterfaceConfig;
  GedimForPy::InterfaceData& data = GedimForPy::GeDiM4Py_Interface::InterfaceData;

  configuration = GedimForPy::GeDiM4Py_Interface::ConvertConfiguration(config);

  GedimForPy::GeDiM4Py_Logic interface;
  interface.Initialize(configuration,
                       data);
}
// ***************************************************************************
namespace GedimForPy
{
  GeDiM4Py_Logic_Configuration GeDiM4Py_Interface::InterfaceConfig;
  InterfaceData GeDiM4Py_Interface::InterfaceData;
  // ***************************************************************************
  GeDiM4Py_Logic_Configuration GeDiM4Py_Interface::ConvertConfiguration(PyObject* config)
  {
    GeDiM4Py_Logic_Configuration configuration;

    configuration.GeometricTolerance = PyFloat_AsDouble(PyDict_GetItemString(config, "GeometricTolerance"));

    return configuration;
  }
  // ***************************************************************************
}
// ***************************************************************************
