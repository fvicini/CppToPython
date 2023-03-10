#include "GeDiM4Py_Interface.hpp"
#include "MeshMatricesDAO.hpp"

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
PyObject* GedimForPy_Discretize(PyObject* discreteSpace)
{
  if (!PyDict_Check(discreteSpace))
    throw std::runtime_error("The input is not correct");

  GedimForPy::InterfaceConfiguration& configuration = GedimForPy::GeDiM4Py_Interface::InterfaceConfig;
  GedimForPy::InterfaceData& data = GedimForPy::GeDiM4Py_Interface::InterfaceData;
  GedimForPy::Domain2D& domain2D = GedimForPy::GeDiM4Py_Interface::Domain;
  GedimForPy::Domain2DMesh& mesh = GedimForPy::GeDiM4Py_Interface::Mesh;
  GedimForPy::DiscreteSpace& space = GedimForPy::GeDiM4Py_Interface::Space;
  GedimForPy::DiscreteProblemData& problemData = GedimForPy::GeDiM4Py_Interface::ProblemData;

  Gedim::MeshMatricesDAO meshDAO(mesh.Mesh);

  space = GedimForPy::GeDiM4Py_Interface::ConvertDiscreteSpace(discreteSpace);
  problemData = GedimForPy::GeDiM4Py_Logic::Discretize(meshDAO,
                                                       space);

  return GedimForPy::GeDiM4Py_Interface::ConvertProblemData(problemData);
}
// ***************************************************************************
void GedimForPy_AssembleStiffnessMatrix(GedimForPy::GeDiM4Py_Logic::K k,
                                        int* numTriplets,
                                        double** stiffnessTriplets)
{
  GedimForPy::InterfaceConfiguration& configuration = GedimForPy::GeDiM4Py_Interface::InterfaceConfig;
  GedimForPy::InterfaceData& data = GedimForPy::GeDiM4Py_Interface::InterfaceData;
  GedimForPy::Domain2D& domain2D = GedimForPy::GeDiM4Py_Interface::Domain;
  GedimForPy::Domain2DMesh& mesh = GedimForPy::GeDiM4Py_Interface::Mesh;
  GedimForPy::DiscreteSpace& space = GedimForPy::GeDiM4Py_Interface::Space;
  GedimForPy::DiscreteProblemData& problemData = GedimForPy::GeDiM4Py_Interface::ProblemData;

  Gedim::MeshMatricesDAO meshDAO(mesh.Mesh);

  GedimForPy::GeDiM4Py_Interface::ConvertTriplets(GedimForPy::GeDiM4Py_Logic::AssembleStiffnessMatrix(k,
                                                                                                      meshDAO,
                                                                                                      mesh.Cell2DsMap,
                                                                                                      problemData),
                                                  *numTriplets,
                                                  stiffnessTriplets);
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
  GedimForPy::InterfaceConfiguration GeDiM4Py_Interface::InterfaceConfig;
  GedimForPy::InterfaceData GeDiM4Py_Interface::InterfaceData;
  GedimForPy::Domain2D GeDiM4Py_Interface::Domain;
  GedimForPy::Domain2DMesh GeDiM4Py_Interface::Mesh;
  GedimForPy::DiscreteSpace GeDiM4Py_Interface::Space;
  GedimForPy::DiscreteProblemData GeDiM4Py_Interface::ProblemData;
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
  DiscreteSpace GeDiM4Py_Interface::ConvertDiscreteSpace(PyObject* discreteSpace)
  {
    DiscreteSpace space;

    space.Order = PyLong_AsLong(PyDict_GetItemString(discreteSpace, "Order"));
    space.Type = static_cast<DiscreteSpace::Types>(PyLong_AsLong(PyDict_GetItemString(discreteSpace, "Type")));

    std::vector<unsigned int> boundaryConditionsType = ConvertArray<unsigned int>(PyDict_GetItemString(discreteSpace, "BoundaryConditionsType"));
    space.BoundaryConditionsType.resize(boundaryConditionsType.size());
    for (unsigned int b = 0; b < boundaryConditionsType.size(); b++)
      space.BoundaryConditionsType[b] = static_cast<DiscreteSpace::BoundaryConditionTypes>(boundaryConditionsType[b]);

    return space;
  }
  // ***************************************************************************
  PyObject* GeDiM4Py_Interface::ConvertProblemData(DiscreteProblemData& problemData)
  {
    PyObject* problem = PyDict_New();

    PyDict_SetItemString(problem, "NumberDOFs", Py_BuildValue("i", problemData.NumberDOFs));
    PyDict_SetItemString(problem, "NumberStrongs", Py_BuildValue("i", problemData.NumberStrongs));

    return problem;
  }
  // ***************************************************************************
  void GeDiM4Py_Interface::ConvertTriplets(const std::list<Eigen::Triplet<double>>& triplets,
                                           int& numTriplets,
                                           double** convertedTriplets)
  {
    numTriplets = triplets.size();
    *convertedTriplets = new double[3 * numTriplets];

    Eigen::Map<Eigen::MatrixXd> tripl(*convertedTriplets, 3, numTriplets);

    unsigned int t = 0;
    for (const Eigen::Triplet<double>& triplet : triplets)
      tripl.col(t++)<< triplet.row(), triplet.col(), triplet.value();
  }
  // ***************************************************************************
}
// ***************************************************************************
