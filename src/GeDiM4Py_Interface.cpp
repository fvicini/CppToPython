#include "GeDiM4Py_Interface.hpp"
#include "FileTextReader.hpp"
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
PyObject* GedimForPy_CreateDomainSquare(PyObject* square,
                                        double** coordinates)
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

  Gedim::MeshMatricesDAO meshDAO(mesh.Mesh);

  return GedimForPy::GeDiM4Py_Interface::ConvertMesh(meshDAO,
                                                     *coordinates);
}
// ***************************************************************************
PyObject* GedimForPy_CreateDomainRectangle(PyObject* rectangle,
                                           double** coordinates)
{
  if (!PyDict_Check(rectangle))
    throw std::runtime_error("The input is not correct");

  GedimForPy::InterfaceConfiguration& configuration = GedimForPy::GeDiM4Py_Interface::InterfaceConfig;
  GedimForPy::InterfaceData& data = GedimForPy::GeDiM4Py_Interface::InterfaceData;
  GedimForPy::Domain2D& domain2D = GedimForPy::GeDiM4Py_Interface::Domain;
  GedimForPy::Domain2DMesh& mesh = GedimForPy::GeDiM4Py_Interface::Mesh;

  GedimForPy::InterfaceDataDAO gedimData = GedimForPy::InterfaceDataDAO(data);

  domain2D = GedimForPy::GeDiM4Py_Interface::ConvertDomainRectangle(gedimData,
                                                                    rectangle);
  mesh = GedimForPy::GeDiM4Py_Logic::CreateDomainMesh2D(domain2D,
                                                        gedimData);

  Gedim::MeshMatricesDAO meshDAO(mesh.Mesh);

  return GedimForPy::GeDiM4Py_Interface::ConvertMesh(meshDAO,
                                                     *coordinates);
}
// ***************************************************************************
PyObject* GedimForPy_ImportDomainMesh2D(double** coordinates)
{
  GedimForPy::InterfaceConfiguration& configuration = GedimForPy::GeDiM4Py_Interface::InterfaceConfig;
  GedimForPy::InterfaceData& data = GedimForPy::GeDiM4Py_Interface::InterfaceData;
  GedimForPy::ImportMesh2D& domain2D = GedimForPy::GeDiM4Py_Interface::DomainImported;
  GedimForPy::Domain2DMesh& mesh = GedimForPy::GeDiM4Py_Interface::Mesh;

  GedimForPy::InterfaceDataDAO gedimData = GedimForPy::InterfaceDataDAO(data);

  domain2D = GedimForPy::GeDiM4Py_Interface::ConvertDomainImport(gedimData);
  mesh = GedimForPy::GeDiM4Py_Logic::ImportDomainMesh2D(domain2D,
                                                        gedimData);

  Gedim::MeshMatricesDAO meshDAO(mesh.Mesh);

  return GedimForPy::GeDiM4Py_Interface::ConvertMesh(meshDAO,
                                                     *coordinates);
}
// ***************************************************************************
PyObject* GedimForPy_Discretize(PyObject* discreteSpace,
                                double** dofsCoordinate,
                                double** strongsCoordinate)
{
  if (!PyDict_Check(discreteSpace))
    throw std::runtime_error("The input is not correct");

  GedimForPy::InterfaceConfiguration& configuration = GedimForPy::GeDiM4Py_Interface::InterfaceConfig;
  GedimForPy::InterfaceData& data = GedimForPy::GeDiM4Py_Interface::InterfaceData;
  GedimForPy::Domain2D& domain2D = GedimForPy::GeDiM4Py_Interface::Domain;
  GedimForPy::Domain2DMesh& mesh = GedimForPy::GeDiM4Py_Interface::Mesh;
  std::vector<GedimForPy::DiscreteSpace>& spaces = GedimForPy::GeDiM4Py_Interface::Spaces;
  std::vector<GedimForPy::DiscreteProblemData>& problemsData = GedimForPy::GeDiM4Py_Interface::ProblemsData;

  Gedim::MeshMatricesDAO meshDAO(mesh.Mesh);

  spaces.push_back(GedimForPy::GeDiM4Py_Interface::ConvertDiscreteSpace(discreteSpace));
  problemsData.push_back(GedimForPy::GeDiM4Py_Logic::Discretize(meshDAO,
                                                                mesh.MeshGeometricData,
                                                                spaces.back()));

  return GedimForPy::GeDiM4Py_Interface::ConvertProblemData(problemsData.size() - 1,
                                                            problemsData.back(),
                                                            *dofsCoordinate,
                                                            *strongsCoordinate);
}
// ***************************************************************************
void GedimForPy_AssembleStiffnessMatrix(const int trialSpaceIndex,
                                        const int testSpaceIndex,
                                        GedimForPy::GeDiM4Py_Logic::A a,
                                        int* numStiffnessTriplets,
                                        double** stiffnessTriplets,
                                        int* numStiffnessStrongTriplets,
                                        double** stiffnessStrongTriplets)
{
  GedimForPy::InterfaceConfiguration& configuration = GedimForPy::GeDiM4Py_Interface::InterfaceConfig;
  GedimForPy::InterfaceData& data = GedimForPy::GeDiM4Py_Interface::InterfaceData;
  GedimForPy::Domain2D& domain2D = GedimForPy::GeDiM4Py_Interface::Domain;
  GedimForPy::Domain2DMesh& mesh = GedimForPy::GeDiM4Py_Interface::Mesh;
  std::vector<GedimForPy::DiscreteSpace>& spaces = GedimForPy::GeDiM4Py_Interface::Spaces;
  std::vector<GedimForPy::DiscreteProblemData>& problemsData = GedimForPy::GeDiM4Py_Interface::ProblemsData;

  Gedim::MeshMatricesDAO meshDAO(mesh.Mesh);

  std::list<Eigen::Triplet<double>> stiffness, stiffnessStrong;
  GedimForPy::GeDiM4Py_Logic::AssembleStiffnessMatrix(a,
                                                      meshDAO,
                                                      mesh.Cell2DsMap,
                                                      problemsData[trialSpaceIndex],
                                                      stiffness,
                                                      stiffnessStrong);

  GedimForPy::GeDiM4Py_Interface::ConvertTriplets(stiffness,
                                                  *numStiffnessTriplets,
                                                  *stiffnessTriplets);

  GedimForPy::GeDiM4Py_Interface::ConvertTriplets(stiffnessStrong,
                                                  *numStiffnessStrongTriplets,
                                                  *stiffnessStrongTriplets);
}
// ***************************************************************************
void GedimForPy_AssembleAnisotropicStiffnessMatrix(const int trialSpaceIndex,
                                                   const int testSpaceIndex,
                                                   GedimForPy::GeDiM4Py_Logic::A a,
                                                   int* numStiffnessTriplets,
                                                   double** stiffnessTriplets,
                                                   int* numStiffnessStrongTriplets,
                                                   double** stiffnessStrongTriplets)
{
  GedimForPy::InterfaceConfiguration& configuration = GedimForPy::GeDiM4Py_Interface::InterfaceConfig;
  GedimForPy::InterfaceData& data = GedimForPy::GeDiM4Py_Interface::InterfaceData;
  GedimForPy::Domain2D& domain2D = GedimForPy::GeDiM4Py_Interface::Domain;
  GedimForPy::Domain2DMesh& mesh = GedimForPy::GeDiM4Py_Interface::Mesh;
  std::vector<GedimForPy::DiscreteSpace>& spaces = GedimForPy::GeDiM4Py_Interface::Spaces;
  std::vector<GedimForPy::DiscreteProblemData>& problemsData = GedimForPy::GeDiM4Py_Interface::ProblemsData;

  Gedim::MeshMatricesDAO meshDAO(mesh.Mesh);

  std::list<Eigen::Triplet<double>> stiffness, stiffnessStrong;
  GedimForPy::GeDiM4Py_Logic::AssembleAnisotropicStiffnessMatrix(a,
                                                                 meshDAO,
                                                                 mesh.Cell2DsMap,
                                                                 problemsData[trialSpaceIndex],
                                                                 stiffness,
                                                                 stiffnessStrong);

  GedimForPy::GeDiM4Py_Interface::ConvertTriplets(stiffness,
                                                  *numStiffnessTriplets,
                                                  *stiffnessTriplets);

  GedimForPy::GeDiM4Py_Interface::ConvertTriplets(stiffnessStrong,
                                                  *numStiffnessStrongTriplets,
                                                  *stiffnessStrongTriplets);
}
// ***************************************************************************
void GedimForPy_AssembleAdvectionMatrix(const int trialSpaceIndex,
                                        const int testSpaceIndex,
                                        GedimForPy::GeDiM4Py_Logic::B b,
                                        int* numAdvectionTriplets,
                                        double** advectionTriplets,
                                        int* numAdvectionStrongTriplets,
                                        double** advectionStrongTriplets)
{
  GedimForPy::InterfaceConfiguration& configuration = GedimForPy::GeDiM4Py_Interface::InterfaceConfig;
  GedimForPy::InterfaceData& data = GedimForPy::GeDiM4Py_Interface::InterfaceData;
  GedimForPy::Domain2D& domain2D = GedimForPy::GeDiM4Py_Interface::Domain;
  GedimForPy::Domain2DMesh& mesh = GedimForPy::GeDiM4Py_Interface::Mesh;
  std::vector<GedimForPy::DiscreteSpace>& spaces = GedimForPy::GeDiM4Py_Interface::Spaces;
  std::vector<GedimForPy::DiscreteProblemData>& problemsData = GedimForPy::GeDiM4Py_Interface::ProblemsData;

  Gedim::MeshMatricesDAO meshDAO(mesh.Mesh);

  std::list<Eigen::Triplet<double>> advection, advectionStrong;
  GedimForPy::GeDiM4Py_Logic::AssembleAdvectionMatrix(b,
                                                      meshDAO,
                                                      mesh.Cell2DsMap,
                                                      problemsData[trialSpaceIndex],
                                                      problemsData[testSpaceIndex],
                                                      advection,
                                                      advectionStrong);

  GedimForPy::GeDiM4Py_Interface::ConvertTriplets(advection,
                                                  *numAdvectionTriplets,
                                                  *advectionTriplets);

  GedimForPy::GeDiM4Py_Interface::ConvertTriplets(advectionStrong,
                                                  *numAdvectionStrongTriplets,
                                                  *advectionStrongTriplets);
}
// ***************************************************************************
void GedimForPy_AssembleReactionMatrix(const int trialSpaceIndex,
                                       const int testSpaceIndex,
                                       GedimForPy::GeDiM4Py_Logic::C c,
                                       int* numReactionTriplets,
                                       double** reactionTriplets,
                                       int* numReactionStrongTriplets,
                                       double** reactionStrongTriplets)
{
  GedimForPy::InterfaceConfiguration& configuration = GedimForPy::GeDiM4Py_Interface::InterfaceConfig;
  GedimForPy::InterfaceData& data = GedimForPy::GeDiM4Py_Interface::InterfaceData;
  GedimForPy::Domain2D& domain2D = GedimForPy::GeDiM4Py_Interface::Domain;
  GedimForPy::Domain2DMesh& mesh = GedimForPy::GeDiM4Py_Interface::Mesh;
  std::vector<GedimForPy::DiscreteSpace>& spaces = GedimForPy::GeDiM4Py_Interface::Spaces;
  std::vector<GedimForPy::DiscreteProblemData>& problemsData = GedimForPy::GeDiM4Py_Interface::ProblemsData;

  Gedim::MeshMatricesDAO meshDAO(mesh.Mesh);

  std::list<Eigen::Triplet<double>> reaction, reactionStrong;
  GedimForPy::GeDiM4Py_Logic::AssembleReactionMatrix(c,
                                                     meshDAO,
                                                     mesh.Cell2DsMap,
                                                     problemsData[trialSpaceIndex],
                                                     reaction,
                                                     reactionStrong);

  GedimForPy::GeDiM4Py_Interface::ConvertTriplets(reaction,
                                                  *numReactionTriplets,
                                                  *reactionTriplets);

  GedimForPy::GeDiM4Py_Interface::ConvertTriplets(reactionStrong,
                                                  *numReactionStrongTriplets,
                                                  *reactionStrongTriplets);
}
// ***************************************************************************
void GedimForPy_AssembleForcingTerm(const int testSpaceIndex,
                                    GedimForPy::GeDiM4Py_Logic::F f,
                                    int* size,
                                    double** forcingTerm)
{
  GedimForPy::InterfaceConfiguration& configuration = GedimForPy::GeDiM4Py_Interface::InterfaceConfig;
  GedimForPy::InterfaceData& data = GedimForPy::GeDiM4Py_Interface::InterfaceData;
  GedimForPy::Domain2D& domain2D = GedimForPy::GeDiM4Py_Interface::Domain;
  GedimForPy::Domain2DMesh& mesh = GedimForPy::GeDiM4Py_Interface::Mesh;
  std::vector<GedimForPy::DiscreteSpace>& spaces = GedimForPy::GeDiM4Py_Interface::Spaces;
  std::vector<GedimForPy::DiscreteProblemData>& problemsData = GedimForPy::GeDiM4Py_Interface::ProblemsData;

  Gedim::MeshMatricesDAO meshDAO(mesh.Mesh);

  GedimForPy::GeDiM4Py_Interface::ConvertArray(GedimForPy::GeDiM4Py_Logic::AssembleForcingTerm(f,
                                                                                               meshDAO,
                                                                                               mesh.Cell2DsMap,
                                                                                               problemsData[testSpaceIndex]),
                                               *size,
                                               *forcingTerm);
}
// ***************************************************************************
void GedimForPy_AssembleStrongSolution(const int trialSpaceIndex,
                                       GedimForPy::GeDiM4Py_Logic::Strong g,
                                       int marker,
                                       int* size,
                                       double** solutionStrong)
{
  GedimForPy::InterfaceConfiguration& configuration = GedimForPy::GeDiM4Py_Interface::InterfaceConfig;
  GedimForPy::InterfaceData& data = GedimForPy::GeDiM4Py_Interface::InterfaceData;
  GedimForPy::Domain2D& domain2D = GedimForPy::GeDiM4Py_Interface::Domain;
  GedimForPy::Domain2DMesh& mesh = GedimForPy::GeDiM4Py_Interface::Mesh;
  std::vector<GedimForPy::DiscreteSpace>& spaces = GedimForPy::GeDiM4Py_Interface::Spaces;
  std::vector<GedimForPy::DiscreteProblemData>& problemsData = GedimForPy::GeDiM4Py_Interface::ProblemsData;

  Gedim::MeshMatricesDAO meshDAO(mesh.Mesh);

  GedimForPy::GeDiM4Py_Interface::ConvertArray(GedimForPy::GeDiM4Py_Logic::AssembleStrongSolution(g,
                                                                                                  marker,
                                                                                                  meshDAO,
                                                                                                  mesh.Cell2DsMap,
                                                                                                  problemsData[trialSpaceIndex]),
                                               *size,
                                               *solutionStrong);
}
// ***************************************************************************
void GedimForPy_AssembleWeakTerm(const int testSpaceIndex,
                                 GedimForPy::GeDiM4Py_Logic::Weak g,
                                 int marker,
                                 int* size,
                                 double** weakTerm)
{
  GedimForPy::InterfaceConfiguration& configuration = GedimForPy::GeDiM4Py_Interface::InterfaceConfig;
  GedimForPy::InterfaceData& data = GedimForPy::GeDiM4Py_Interface::InterfaceData;
  GedimForPy::Domain2D& domain2D = GedimForPy::GeDiM4Py_Interface::Domain;
  GedimForPy::Domain2DMesh& mesh = GedimForPy::GeDiM4Py_Interface::Mesh;
  std::vector<GedimForPy::DiscreteSpace>& spaces = GedimForPy::GeDiM4Py_Interface::Spaces;
  std::vector<GedimForPy::DiscreteProblemData>& problemsData = GedimForPy::GeDiM4Py_Interface::ProblemsData;

  Gedim::MeshMatricesDAO meshDAO(mesh.Mesh);

  GedimForPy::GeDiM4Py_Interface::ConvertArray(GedimForPy::GeDiM4Py_Logic::AssembleWeakTerm(g,
                                                                                            marker,
                                                                                            meshDAO,
                                                                                            mesh.MeshGeometricData.Cell2DsVertices,
                                                                                            mesh.MeshGeometricData.Cell2DsEdgeLengths,
                                                                                            mesh.MeshGeometricData.Cell2DsEdgeTangents,
                                                                                            mesh.Cell2DsMap,
                                                                                            problemsData[testSpaceIndex]),
                                               *size,
                                               *weakTerm);
}
// ***************************************************************************
void GedimForPy_CholeskySolver(const int nRows,
                               const int numNonZeros,
                               const double* pointerTriplets,
                               const double* pointerF,
                               double** solution)
{
  const Eigen::SparseMatrix<double> A = GedimForPy::GeDiM4Py_Interface::ConvertSparseMatrix(nRows,
                                                                                            nRows,
                                                                                            numNonZeros,
                                                                                            pointerTriplets);

  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>, Eigen::Lower> linearSolver;
  linearSolver.compute(A);

  const Eigen::VectorXd sol = linearSolver.solve(Eigen::Map<const Eigen::VectorXd>(pointerF, nRows));

  int size = 0;
  GedimForPy::GeDiM4Py_Interface::ConvertArray(sol,
                                               size,
                                               *solution);
}
// ***************************************************************************
void GedimForPy_LUSolver(const int nRows,
                         const int numNonZeros,
                         const double* pointerTriplets,
                         const double* pointerF,
                         double** solution)
{
  const Eigen::SparseMatrix<double> A = GedimForPy::GeDiM4Py_Interface::ConvertSparseMatrix(nRows,
                                                                                            nRows,
                                                                                            numNonZeros,
                                                                                            pointerTriplets);

  Eigen::SparseLU<Eigen::SparseMatrix<double>> linearSolver;
  linearSolver.compute(A);

  const Eigen::VectorXd sol = linearSolver.solve(Eigen::Map<const Eigen::VectorXd>(pointerF, nRows));

  int size = 0;
  GedimForPy::GeDiM4Py_Interface::ConvertArray(sol,
                                               size,
                                               *solution);
}
// ***************************************************************************
double GedimForPy_ComputeErrorL2(const int trialSpaceIndex,
                                 GedimForPy::GeDiM4Py_Logic::Exact u,
                                 const double* pointerNumericSolution,
                                 const double* pointerStrongSolution)
{
  GedimForPy::InterfaceConfiguration& configuration = GedimForPy::GeDiM4Py_Interface::InterfaceConfig;
  GedimForPy::InterfaceData& data = GedimForPy::GeDiM4Py_Interface::InterfaceData;
  GedimForPy::Domain2D& domain2D = GedimForPy::GeDiM4Py_Interface::Domain;
  GedimForPy::Domain2DMesh& mesh = GedimForPy::GeDiM4Py_Interface::Mesh;
  std::vector<GedimForPy::DiscreteSpace>& spaces = GedimForPy::GeDiM4Py_Interface::Spaces;
  std::vector<GedimForPy::DiscreteProblemData>& problemsData = GedimForPy::GeDiM4Py_Interface::ProblemsData;

  Gedim::MeshMatricesDAO meshDAO(mesh.Mesh);

  const Eigen::VectorXd cell2DsErrorL2 = GedimForPy::GeDiM4Py_Logic::ComputeErrorL2(u,
                                                                                    Eigen::Map<const Eigen::VectorXd>(pointerNumericSolution,
                                                                                                                      problemsData[trialSpaceIndex].NumberDOFs),
                                                                                    Eigen::Map<const Eigen::VectorXd>(pointerStrongSolution,
                                                                                                                      problemsData[trialSpaceIndex].NumberStrongs),
                                                                                    meshDAO,
                                                                                    mesh.Cell2DsMap,
                                                                                    problemsData[trialSpaceIndex]);
  return sqrt(cell2DsErrorL2.sum());
}
// ***************************************************************************
double GedimForPy_ComputeErrorL2_LastSpace(GedimForPy::GeDiM4Py_Logic::Exact u,
                                           const double* pointerNumericSolution,
                                           const double* pointerStrongSolution)
{
  return GedimForPy_ComputeErrorL2(GedimForPy::GeDiM4Py_Interface::ProblemsData.size() - 1,
                                   u,
                                   pointerNumericSolution,
                                   pointerStrongSolution);
}
// ***************************************************************************
double GedimForPy_ComputeErrorH1(const int trialSpaceIndex,
                                 GedimForPy::GeDiM4Py_Logic::ExactDerivative uDer,
                                 const double* pointerNumericSolution,
                                 const double* pointerStrongSolution)
{
  GedimForPy::InterfaceConfiguration& configuration = GedimForPy::GeDiM4Py_Interface::InterfaceConfig;
  GedimForPy::InterfaceData& data = GedimForPy::GeDiM4Py_Interface::InterfaceData;
  GedimForPy::Domain2D& domain2D = GedimForPy::GeDiM4Py_Interface::Domain;
  GedimForPy::Domain2DMesh& mesh = GedimForPy::GeDiM4Py_Interface::Mesh;
  std::vector<GedimForPy::DiscreteSpace>& spaces = GedimForPy::GeDiM4Py_Interface::Spaces;
  std::vector<GedimForPy::DiscreteProblemData>& problemsData = GedimForPy::GeDiM4Py_Interface::ProblemsData;

  Gedim::MeshMatricesDAO meshDAO(mesh.Mesh);

  const Eigen::VectorXd cell2DsErrorH1 = GedimForPy::GeDiM4Py_Logic::ComputeErrorH1(uDer,
                                                                                    Eigen::Map<const Eigen::VectorXd>(pointerNumericSolution,
                                                                                                                      problemsData[trialSpaceIndex].NumberDOFs),
                                                                                    Eigen::Map<const Eigen::VectorXd>(pointerStrongSolution,
                                                                                                                      problemsData[trialSpaceIndex].NumberStrongs),
                                                                                    meshDAO,
                                                                                    mesh.Cell2DsMap,
                                                                                    problemsData[trialSpaceIndex]);
  return sqrt(cell2DsErrorH1.sum());
}
// ***************************************************************************
double GedimForPy_ComputeErrorH1_LastSpace(GedimForPy::GeDiM4Py_Logic::ExactDerivative uDer,
                                           const double* pointerNumericSolution,
                                           const double* pointerStrongSolution)
{
  return GedimForPy_ComputeErrorH1(GedimForPy::GeDiM4Py_Interface::ProblemsData.size() - 1,
                                   uDer,
                                   pointerNumericSolution,
                                   pointerStrongSolution);
}
// ***************************************************************************
void GedimForPy_ExportSolution(const int trialSpaceIndex,
                               GedimForPy::GeDiM4Py_Logic::Exact u,
                               const double* pointerNumericSolution,
                               const double* pointerStrongSolution)
{
  GedimForPy::InterfaceConfiguration& configuration = GedimForPy::GeDiM4Py_Interface::InterfaceConfig;
  GedimForPy::InterfaceData& data = GedimForPy::GeDiM4Py_Interface::InterfaceData;
  GedimForPy::Domain2D& domain2D = GedimForPy::GeDiM4Py_Interface::Domain;
  GedimForPy::Domain2DMesh& mesh = GedimForPy::GeDiM4Py_Interface::Mesh;
  std::vector<GedimForPy::DiscreteSpace>& spaces = GedimForPy::GeDiM4Py_Interface::Spaces;
  std::vector<GedimForPy::DiscreteProblemData>& problemsData = GedimForPy::GeDiM4Py_Interface::ProblemsData;

  Gedim::MeshMatricesDAO meshDAO(mesh.Mesh);

  GedimForPy::GeDiM4Py_Logic::ExportSolution(u,
                                             Eigen::Map<const Eigen::VectorXd>(pointerNumericSolution,
                                                                               problemsData[trialSpaceIndex].NumberDOFs),
                                             Eigen::Map<const Eigen::VectorXd>(pointerStrongSolution,
                                                                               problemsData[trialSpaceIndex].NumberStrongs),
                                             meshDAO,
                                             mesh.Cell2DsMap,
                                             problemsData[trialSpaceIndex]);
}
// ***************************************************************************
void GedimForPy_ExportSolution_LastSpace(GedimForPy::GeDiM4Py_Logic::Exact u,
                                         const double* pointerNumericSolution,
                                         const double* pointerStrongSolution)
{
  GedimForPy_ExportSolution(GedimForPy::GeDiM4Py_Interface::ProblemsData.size() - 1,
                            u,
                            pointerNumericSolution,
                            pointerStrongSolution);
}
// ***************************************************************************
namespace GedimForPy
{
  // ***************************************************************************
  template<typename T>
  std::vector<T> GeDiM4Py_Interface::ConvertToArray(PyObject* list)
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
  GedimForPy::ImportMesh2D GeDiM4Py_Interface::DomainImported;
  GedimForPy::Domain2DMesh GeDiM4Py_Interface::Mesh;
  std::vector<GedimForPy::DiscreteSpace> GeDiM4Py_Interface::Spaces;
  std::vector<GedimForPy::DiscreteProblemData> GeDiM4Py_Interface::ProblemsData;
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
    domain.VerticesBoundaryCondition = ConvertToArray<unsigned int>(PyDict_GetItemString(square, "VerticesBoundaryCondition"));
    domain.EdgesBoundaryCondition = ConvertToArray<unsigned int>(PyDict_GetItemString(square, "EdgesBoundaryCondition"));
    domain.DiscretizationType = static_cast<Domain2D::DiscretizationTypes>(PyLong_AsLong(PyDict_GetItemString(square, "DiscretizationType")));
    domain.MeshCellsMaximumArea = PyFloat_AsDouble(PyDict_GetItemString(square, "MeshCellsMaximumArea"));

    return domain;
  }
  // ***************************************************************************
  Domain2D GeDiM4Py_Interface::ConvertDomainRectangle(GedimForPy::InterfaceDataDAO& gedimData,
                                                      PyObject* rectangle)
  {
    Domain2D domain;

    const double rectangleBase = PyFloat_AsDouble(PyDict_GetItemString(rectangle, "RectangleBase"));
    const double rectangleHeight = PyFloat_AsDouble(PyDict_GetItemString(rectangle, "RectangleHeight"));

    domain.Vertices = gedimData.GeometryUtilities().CreateRectangle(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                                    rectangleBase,
                                                                    rectangleHeight);
    domain.VerticesBoundaryCondition = ConvertToArray<unsigned int>(PyDict_GetItemString(rectangle, "VerticesBoundaryCondition"));
    domain.EdgesBoundaryCondition = ConvertToArray<unsigned int>(PyDict_GetItemString(rectangle, "EdgesBoundaryCondition"));
    domain.DiscretizationType = static_cast<Domain2D::DiscretizationTypes>(PyLong_AsLong(PyDict_GetItemString(rectangle, "DiscretizationType")));
    domain.MeshCellsMaximumArea = PyFloat_AsDouble(PyDict_GetItemString(rectangle, "MeshCellsMaximumArea"));

    return domain;
  }
  // ***************************************************************************
  ImportMesh2D GeDiM4Py_Interface::ConvertDomainImport(InterfaceDataDAO& gedimData)
  {
    ImportMesh2D domain;

    std::vector<std::string> inputFolders;
    Gedim::FileReader csvFileReader("./ImportMesh.csv");

    if (!csvFileReader.Open())
      throw runtime_error("File ImportMesh.csv not found in folder ./");

    csvFileReader.GetAllLines(inputFolders);
    csvFileReader.Close();

    unsigned int numInputFolders = inputFolders.size() - 1;
    if (numInputFolders != 1)
      throw runtime_error("File ImportMesh.csv wrong");

    domain.InputFolder = inputFolders[1];
    domain.Separator = ';';

    return domain;
  }
  // ***************************************************************************
  PyObject* GeDiM4Py_Interface::ConvertMesh(const Gedim::IMeshDAO& mesh,
                                            double*& coordinates)
  {
    PyObject* meshInfo = PyDict_New();

    PyDict_SetItemString(meshInfo, "NumberCell0Ds", Py_BuildValue("i", mesh.Cell0DTotalNumber()));
    PyDict_SetItemString(meshInfo, "NumberCell1Ds", Py_BuildValue("i", mesh.Cell1DTotalNumber()));
    PyDict_SetItemString(meshInfo, "NumberCell2Ds", Py_BuildValue("i", mesh.Cell2DTotalNumber()));

    ConvertMatrix(mesh.Cell0DsCoordinates(),
                  coordinates);

    return meshInfo;
  }
  // ***************************************************************************
  DiscreteSpace GeDiM4Py_Interface::ConvertDiscreteSpace(PyObject* discreteSpace)
  {
    DiscreteSpace space;

    space.Order = PyLong_AsLong(PyDict_GetItemString(discreteSpace, "Order"));
    space.Type = static_cast<DiscreteSpace::Types>(PyLong_AsLong(PyDict_GetItemString(discreteSpace, "Type")));

    std::vector<unsigned int> boundaryConditionsType = ConvertToArray<unsigned int>(PyDict_GetItemString(discreteSpace, "BoundaryConditionsType"));
    space.BoundaryConditionsType.resize(boundaryConditionsType.size());
    for (unsigned int b = 0; b < boundaryConditionsType.size(); b++)
      space.BoundaryConditionsType[b] = static_cast<DiscreteSpace::BoundaryConditionTypes>(boundaryConditionsType[b]);

    return space;
  }
  // ***************************************************************************
  PyObject* GeDiM4Py_Interface::ConvertProblemData(const unsigned int& spaceIndex,
                                                   DiscreteProblemData& problemData,
                                                   double*& dofsCoordinate,
                                                   double*& strongsCoordinate)
  {
    PyObject* problem = PyDict_New();

    PyDict_SetItemString(problem, "SpaceIndex", Py_BuildValue("i", spaceIndex));
    PyDict_SetItemString(problem, "NumberDOFs", Py_BuildValue("i", problemData.NumberDOFs));
    PyDict_SetItemString(problem, "NumberStrongs", Py_BuildValue("i", problemData.NumberStrongs));
    PyDict_SetItemString(problem, "H", Py_BuildValue("d", problemData.H));

    ConvertMatrix(problemData.DOFsCoordinate, dofsCoordinate);
    ConvertMatrix(problemData.StrongsCoordinate, strongsCoordinate);

    return problem;
  }
  // ***************************************************************************
  void GeDiM4Py_Interface::ConvertTriplets(const std::list<Eigen::Triplet<double>>& triplets,
                                           int& numTriplets,
                                           double*& convertedTriplets)
  {
    numTriplets = triplets.size();
    convertedTriplets = new double[3 * numTriplets];

    Eigen::Map<Eigen::MatrixXd> tripl(convertedTriplets, 3, numTriplets);

    unsigned int t = 0;
    for (const Eigen::Triplet<double>& triplet : triplets)
      tripl.col(t++)<< triplet.row(), triplet.col(), triplet.value();
  }
  // ***************************************************************************
  Eigen::SparseMatrix<double> GeDiM4Py_Interface::ConvertSparseMatrix(const int& nRows,
                                                                      const int& nCols,
                                                                      const int& numNonZeros,
                                                                      const double* pointerTriplets)
  {
    Eigen::Map<const Eigen::MatrixXd> triplets(pointerTriplets, 3, numNonZeros);

    std::vector<Eigen::Triplet<double>> matTriplets(numNonZeros);
    for (unsigned int t = 0; t < numNonZeros; t++)
      matTriplets[t] = Eigen::Triplet<double>(triplets(0, t), triplets(1, t), triplets(2, t));

    Eigen::SparseMatrix<double> mat(nRows,
                                    nCols);
    mat.setFromTriplets(matTriplets.begin(),
                        matTriplets.end());
    mat.makeCompressed();

    return mat;
  }
  // ***************************************************************************
  void GeDiM4Py_Interface::ConvertMatrix(const Eigen::MatrixXd& matrix,
                                         double*& convertedMatrix)
  {
    unsigned int size = matrix.rows() * matrix.cols();
    convertedMatrix = new double[size];

    Eigen::Map<Eigen::MatrixXd> mat(convertedMatrix, matrix.rows(), matrix.cols());
    mat = matrix;
  }
  // ***************************************************************************
  void GeDiM4Py_Interface::ConvertArray(const Eigen::VectorXd& array,
                                        int& size,
                                        double*& convertedArray)
  {
    size = array.size();
    convertedArray = new double[size];

    Eigen::Map<Eigen::VectorXd> arr(convertedArray, size);
    arr = array;
  }
  // ***************************************************************************
}
// ***************************************************************************
