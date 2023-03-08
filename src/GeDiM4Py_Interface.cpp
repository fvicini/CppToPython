#include "GeDiM4Py_Interface.hpp"

// ***************************************************************************
void GedimForPy_Initialize(PyObject* config)
{
  if (!PyDict_Check(config))
    throw std::runtime_error("The input is not correct");

  GedimForPy::InterfaceConfiguration& configuration = GedimForPy::GeDiM4Py_Interface::InterfaceConfig;
  GedimForPy::InterfaceData& data = GedimForPy::GeDiM4Py_Interface::InterfaceData;

  configuration = GedimForPy::GeDiM4Py_Interface::ConvertConfiguration(config);

  GedimForPy::GeDiM4Py_Logic::Initialize(configuration,
                                         data);
}
// ***************************************************************************
void GedimForPy_CreateDomainSquare(PyObject* square)
{
  if (!PyDict_Check(square))
    throw std::runtime_error("The input is not correct");

  GedimForPy::InterfaceConfiguration& configuration = GedimForPy::GeDiM4Py_Interface::InterfaceConfig;
  GedimForPy::InterfaceData& data = GedimForPy::GeDiM4Py_Interface::InterfaceData;
  GedimForPy::Domain2D& domain2D = GedimForPy::GeDiM4Py_Interface::Domain;
  GedimForPy::Domain2DMesh& mesh = GedimForPy::GeDiM4Py_Interface::Mesh;

  GedimForPy::InterfaceDataDAO gedimData = GedimForPy::InterfaceDataDAO(data);

  domain2D = GedimForPy::GeDiM4Py_Interface::ConvertDomainSquare(gedimData,
                                                                 square);
  mesh = GedimForPy::GeDiM4Py_Logic::CreateDomainMesh2D(domain2D,
                                                        gedimData);
}
// ***************************************************************************
namespace GedimForPy
{
  // ***************************************************************************
  template<typename T>
  std::vector<T> GeDiM4Py_Interface::ConvertArray(PyObject* list)
  {
    std::vector<T> result;

    PyList_Check(list);

    const unsigned int size = PyList_Size(list);
    result.resize(size);

    for (unsigned int i = 0; i < size; i++)
    {
      PyObject* value = PyList_GET_ITEM(list, i);

      if (PyLong_Check(value))
        result[i] = PyLong_AsLong(value);
      else if (PyFloat_Check(value))
        result[i] = PyFloat_AsDouble(value);
      else
        throw std::runtime_error("Unknown value type");
    }

    return result;
  }
  // ***************************************************************************
  InterfaceConfiguration GeDiM4Py_Interface::InterfaceConfig;
  InterfaceData GeDiM4Py_Interface::InterfaceData;
  Domain2D GeDiM4Py_Interface::Domain;
  Domain2DMesh GeDiM4Py_Interface::Mesh;
  // ***************************************************************************
  InterfaceConfiguration GeDiM4Py_Interface::ConvertConfiguration(PyObject* config)
  {
    InterfaceConfiguration configuration;

    configuration.GeometricTolerance = PyFloat_AsDouble(PyDict_GetItemString(config, "GeometricTolerance"));

    return configuration;
  }
  // ***************************************************************************
  Domain2D GeDiM4Py_Interface::ConvertDomainSquare(GedimForPy::InterfaceDataDAO& gedimData,
                                                   PyObject* square)
  {
    Domain2D domain;

    const double squareEdge = PyFloat_AsDouble(PyDict_GetItemString(square, "SquareEdge"));

    domain.Vertices = gedimData.GeometryUtilities().CreateSquare(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                                 squareEdge);
    domain.VerticesBoundaryCondition = ConvertArray<unsigned int>(PyDict_GetItemString(square, "VerticesBoundaryCondition"));
    domain.EdgesBoundaryCondition = ConvertArray<unsigned int>(PyDict_GetItemString(square, "EdgesBoundaryCondition"));
    domain.DiscretizationType = static_cast<Domain2D::DiscretizationTypes>(PyLong_AsLong(PyDict_GetItemString(square, "DiscretizationType")));
    domain.MeshCellsMaximumArea = PyFloat_AsDouble(PyDict_GetItemString(square, "MeshCellsMaximumArea"));

    return domain;
  }
  // ***************************************************************************
}
// ***************************************************************************