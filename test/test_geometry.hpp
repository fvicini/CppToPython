#ifndef __TEST_GEOMETRY_H
#define __TEST_GEOMETRY_H

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>

#include "FileTextReader.hpp"
#include "GeDiM4Py_Logic.hpp"
#include "MeshMatricesDAO.hpp"
#include "VTKUtilities.hpp"
#include "test_Poisson.hpp"
#include "test_heat_conductivity.hpp"
#include "test_Stokes.hpp"
#include "test_Burger.hpp"

#define ACTIVE_CHECK 0

namespace UnitTesting
{
  // ***************************************************************************
  TEST(TestGeometry, Test_SquareLaplace)
  {
    const std::string exportFolder = "./Export/TestGeometry/Test_SquareLaplace";
    Gedim::Output::CreateFolder(exportFolder);

    GedimForPy::InterfaceConfiguration interfaceConfig;
    interfaceConfig.GeometricTolerance = 1.0e-8;

    GedimForPy::InterfaceData data;
    GedimForPy::InterfaceDataDAO gedimData(data);

    GedimForPy::GeDiM4Py_Logic interface;

    ASSERT_NO_THROW(interface.Initialize(interfaceConfig,
                                         data));

    const std::vector<double> meshSize = { 0.1 };
    const unsigned int order = 2;

    for (unsigned int m = 0; m < meshSize.size(); m++)
    {
      GedimForPy::Domain2D domain;
      domain.Vertices = gedimData.GeometryUtilities().CreateSquare(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                                   1.0);
      domain.VerticesBoundaryCondition = { 1, 1, 1, 1 };
      domain.EdgesBoundaryCondition = { 1, 2, 1, 3 };
      domain.DiscretizationType = GedimForPy::Domain2D::DiscretizationTypes::Triangular;
      domain.MeshCellsMaximumArea = meshSize[m];

      GedimForPy::Domain2DMesh mesh = GedimForPy::GeDiM4Py_Logic::CreateDomainMesh2D(domain,
                                                                                     gedimData);

#if ACTIVE_CHECK == 1
      ASSERT_EQ(13, mesh.Mesh.NumberCell0D);
      ASSERT_EQ(28, mesh.Mesh.NumberCell1D);
      ASSERT_EQ(16, mesh.Mesh.NumberCell2D);
#endif

      Gedim::MeshMatricesDAO meshDAO(mesh.Mesh);

      // export
      {
        {
          Gedim::VTKUtilities exporter;
          exporter.AddPolygon(domain.Vertices);
          exporter.Export(exportFolder +
                          "/Domain.vtu");
        }

        {
          Gedim::VTKUtilities exporter;
          gedimData.MeshUtilities().ExportMeshToVTU(meshDAO,
                                                    exportFolder,
                                                    "Mesh");
        }
      }

      GedimForPy::DiscreteSpace discreteSpace;
      discreteSpace.Type = GedimForPy::DiscreteSpace::Types::FEM;
      discreteSpace.Order = order;
      discreteSpace.BoundaryConditionsType = { GedimForPy::DiscreteSpace::BoundaryConditionTypes::None,
                                               GedimForPy::DiscreteSpace::BoundaryConditionTypes::Strong,
                                               GedimForPy::DiscreteSpace::BoundaryConditionTypes::Weak,
                                               GedimForPy::DiscreteSpace::BoundaryConditionTypes::Weak };

      GedimForPy::DiscreteProblemData problemData = GedimForPy::GeDiM4Py_Logic::Discretize(meshDAO,
                                                                                           mesh.MeshGeometricData,
                                                                                           discreteSpace);

#if ACTIVE_CHECK == 1
      ASSERT_EQ(31, problemData.NumberDOFs);
      ASSERT_EQ(10, problemData.NumberStrongs);
      ASSERT_EQ(13, problemData.Cell0Ds_DOF.size());
      ASSERT_EQ(28, problemData.Cell1Ds_DOF.size());
#endif

      // export
      {
        {
          std::vector<double> cell0Ds_DOFType(meshDAO.Cell0DTotalNumber(), 0.0);
          std::vector<double> cell0Ds_DOFGlobalIndex(meshDAO.Cell0DTotalNumber(), 0.0);
          std::vector<double> cell1Ds_DOFType(meshDAO.Cell1DTotalNumber(), 0.0);
          std::vector<double> cell1Ds_DOFGlobalIndex(meshDAO.Cell1DTotalNumber(), 0.0);

          for (unsigned int p = 0; p < problemData.Cell0Ds_DOF.size(); p++)
          {
            const GedimForPy::DiscreteProblemData::DOF& dof = problemData.Cell0Ds_DOF[p];
            cell0Ds_DOFType[p] = (unsigned int)dof.Type;
            cell0Ds_DOFGlobalIndex[p] = dof.Global_Index;
          }

          for (unsigned int p = 0; p < problemData.Cell1Ds_DOF.size(); p++)
          {
            const GedimForPy::DiscreteProblemData::DOF& dof = problemData.Cell1Ds_DOF[p];
            cell1Ds_DOFType[p] = (unsigned int)dof.Type;
            cell1Ds_DOFGlobalIndex[p] = dof.Global_Index;
          }

          Gedim::VTKUtilities exporter;
          exporter.AddSegments(meshDAO.Cell0DsCoordinates(),
                               meshDAO.Cell1DsExtremes(),
                               {
                                 {
                                   "cell0Ds_DOFType",
                                   Gedim::VTPProperty::Formats::Points,
                                   static_cast<unsigned int>(cell0Ds_DOFType.size()),
                                   cell0Ds_DOFType.data()
                                 },
                                 {
                                   "cell0Ds_DOFGlobalIndex",
                                   Gedim::VTPProperty::Formats::Points,
                                   static_cast<unsigned int>(cell0Ds_DOFGlobalIndex.size()),
                                   cell0Ds_DOFGlobalIndex.data()
                                 },
                                 {
                                   "cell1Ds_DOFType",
                                   Gedim::VTPProperty::Formats::Cells,
                                   static_cast<unsigned int>(cell1Ds_DOFType.size()),
                                   cell1Ds_DOFType.data()
                                 },
                                 {
                                   "cell1Ds_DOFGlobalIndex",
                                   Gedim::VTPProperty::Formats::Cells,
                                   static_cast<unsigned int>(cell1Ds_DOFGlobalIndex.size()),
                                   cell1Ds_DOFGlobalIndex.data()
                                 }
                               });
          exporter.Export(exportFolder +
                          "/DOFs.vtu");
        }
      }

      std::list<Eigen::Triplet<double>> stiffnessTriplets, stiffnessStrongTriplets;
      GedimForPy::GeDiM4Py_Logic::AssembleStiffnessMatrix(Poisson::DiffusionTerm,
                                                          meshDAO,
                                                          mesh.Cell2DsMap,
                                                          problemData,
                                                          stiffnessTriplets,
                                                          stiffnessStrongTriplets);
      std::list<Eigen::Triplet<double>> advectionTriplets, advectionStrongTriplets;
      GedimForPy::GeDiM4Py_Logic::AssembleAdvectionMatrix(Poisson::AdvectionTerm,
                                                          meshDAO,
                                                          mesh.Cell2DsMap,
                                                          problemData,
                                                          problemData,
                                                          advectionTriplets,
                                                          advectionStrongTriplets);
      std::list<Eigen::Triplet<double>> reactionTriplets, reactionStrongTriplets;
      GedimForPy::GeDiM4Py_Logic::AssembleReactionMatrix(Poisson::ReactionTerm,
                                                         meshDAO,
                                                         mesh.Cell2DsMap,
                                                         problemData,
                                                         reactionTriplets,
                                                         reactionStrongTriplets);

      const Eigen::VectorXd forcingTerm = GedimForPy::GeDiM4Py_Logic::AssembleForcingTerm(Poisson::ForcingTerm,
                                                                                          meshDAO,
                                                                                          mesh.Cell2DsMap,
                                                                                          problemData);

      const Eigen::VectorXd weakTerm_Right = GedimForPy::GeDiM4Py_Logic::AssembleWeakTerm(Poisson::WeakTerm_Right,
                                                                                          2,
                                                                                          meshDAO,
                                                                                          mesh.MeshGeometricData.Cell2DsVertices,
                                                                                          mesh.MeshGeometricData.Cell2DsEdgeLengths,
                                                                                          mesh.MeshGeometricData.Cell2DsEdgeTangents,
                                                                                          mesh.Cell2DsMap,
                                                                                          problemData);
      const Eigen::VectorXd weakTerm_Left = GedimForPy::GeDiM4Py_Logic::AssembleWeakTerm(Poisson::WeakTerm_Left,
                                                                                         3,
                                                                                         meshDAO,
                                                                                         mesh.MeshGeometricData.Cell2DsVertices,
                                                                                         mesh.MeshGeometricData.Cell2DsEdgeLengths,
                                                                                         mesh.MeshGeometricData.Cell2DsEdgeTangents,
                                                                                         mesh.Cell2DsMap,
                                                                                         problemData);

      const Eigen::VectorXd solutionStrong = GedimForPy::GeDiM4Py_Logic::AssembleStrongSolution(Poisson::ExactSolution,
                                                                                                1,
                                                                                                meshDAO,
                                                                                                mesh.Cell2DsMap,
                                                                                                problemData);

      Eigen::SparseMatrix<double> stiffness(problemData.NumberDOFs,
                                            problemData.NumberDOFs);
      stiffness.setFromTriplets(stiffnessTriplets.begin(),
                                stiffnessTriplets.end());
      stiffness.makeCompressed();
      stiffnessTriplets.clear();

      Eigen::SparseMatrix<double> advection(problemData.NumberDOFs,
                                            problemData.NumberDOFs);
      advection.setFromTriplets(advectionTriplets.begin(),
                                advectionTriplets.end());
      advection.makeCompressed();
      advectionTriplets.clear();

      Eigen::SparseMatrix<double> reaction(problemData.NumberDOFs,
                                           problemData.NumberDOFs);
      reaction.setFromTriplets(reactionTriplets.begin(),
                               reactionTriplets.end());
      reaction.makeCompressed();
      reactionTriplets.clear();

      Eigen::SparseMatrix<double> stiffnessStrong(problemData.NumberDOFs,
                                                  problemData.NumberStrongs);
      stiffnessStrong.setFromTriplets(stiffnessStrongTriplets.begin(),
                                      stiffnessStrongTriplets.end());
      stiffnessStrong.makeCompressed();
      stiffnessStrongTriplets.clear();

      Eigen::SparseMatrix<double> advectionStrong(problemData.NumberDOFs,
                                                  problemData.NumberStrongs);
      advectionStrong.setFromTriplets(advectionStrongTriplets.begin(),
                                      advectionStrongTriplets.end());
      advectionStrong.makeCompressed();
      advectionStrongTriplets.clear();

      Eigen::SparseMatrix<double> reactionStrong(problemData.NumberDOFs,
                                                 problemData.NumberStrongs);
      reactionStrong.setFromTriplets(reactionStrongTriplets.begin(),
                                     reactionStrongTriplets.end());
      reactionStrong.makeCompressed();
      reactionStrongTriplets.clear();

      Eigen::SparseLU<Eigen::SparseMatrix<double>> linearSolver;
      linearSolver.compute(stiffness + advection + reaction);

      const Eigen::VectorXd solution = linearSolver.solve(forcingTerm -
                                                          (stiffnessStrong +
                                                           advectionStrong +
                                                           reactionStrong) * solutionStrong +
                                                          weakTerm_Left +
                                                          weakTerm_Right);

      const Eigen::VectorXd cell2DsErrorL2 = GedimForPy::GeDiM4Py_Logic::ComputeErrorL2(Poisson::ExactSolution,
                                                                                        solution,
                                                                                        solutionStrong,
                                                                                        meshDAO,
                                                                                        mesh.Cell2DsMap,
                                                                                        problemData);
      const Eigen::VectorXd cell2DsErrorH1 = GedimForPy::GeDiM4Py_Logic::ComputeErrorH1(Poisson::ExactDerivativeSolution,
                                                                                        solution,
                                                                                        solutionStrong,
                                                                                        meshDAO,
                                                                                        mesh.Cell2DsMap,
                                                                                        problemData);


#if ACTIVE_CHECK == 0
      std::cerr.precision(16);
      std::cerr<< std::scientific<< "dofs"<< ","<< "h"<< ","<< "errorL2"<< ","<< "errorH1"<< std::endl;
      std::cerr<< std::scientific<< problemData.NumberDOFs<< ","<< problemData.H<< ","<< sqrt(cell2DsErrorL2.sum())<< ","<< sqrt(cell2DsErrorH1.sum())<< std::endl;
#endif

      // export
      {
        GedimForPy::GeDiM4Py_Logic::ExportSolution(Poisson::ExactSolution,
                                                   solution,
                                                   solutionStrong,
                                                   meshDAO,
                                                   problemData,
                                                   {
                                                     exportFolder,
                                                     "Solution"
                                                   });

        {
          std::vector<double> cell0Ds_numeric_solution(meshDAO.Cell0DTotalNumber(),
                                                       0.0);
          const Eigen::MatrixXd coordinates = meshDAO.Cell0DsCoordinates();

          const double* cell0Ds_exact_solution = Poisson::ExactSolution(coordinates.cols(),
                                                                        coordinates.data());

          for (unsigned int p = 0; p < meshDAO.Cell0DTotalNumber(); p++)
          {
            const GedimForPy::DiscreteProblemData::DOF& dof = problemData.Cell0Ds_DOF[p];

            switch (dof.Type)
            {
              case GedimForPy::DiscreteProblemData::DOF::Types::DOF:
                cell0Ds_numeric_solution[p] = solution[dof.Global_Index];
                break;
              case GedimForPy::DiscreteProblemData::DOF::Types::Strong:
                cell0Ds_numeric_solution[p] = solutionStrong[dof.Global_Index];
                break;
              default:
                throw std::runtime_error("DOF Type " +
                                         std::to_string((unsigned int)dof.Type) +
                                         " not supported");
            }
          }

          Gedim::VTKUtilities exporter;
          exporter.AddPolygons(meshDAO.Cell0DsCoordinates(),
                               meshDAO.Cell2DsVertices(),
                               {
                                 {
                                   "cell0Ds_exact_solution",
                                   Gedim::VTPProperty::Formats::Points,
                                   static_cast<unsigned int>(coordinates.cols()),
                                   cell0Ds_exact_solution
                                 },
                                 {
                                   "cell0Ds_numeric_solution",
                                   Gedim::VTPProperty::Formats::Points,
                                   static_cast<unsigned int>(cell0Ds_numeric_solution.size()),
                                   cell0Ds_numeric_solution.data()
                                 },
                                 {
                                   "cell2Ds_errorL2",
                                   Gedim::VTPProperty::Formats::Cells,
                                   static_cast<unsigned int>(cell2DsErrorL2.size()),
                                   cell2DsErrorL2.data()
                                 },
                                 {
                                   "cell2Ds_errorH1",
                                   Gedim::VTPProperty::Formats::Cells,
                                   static_cast<unsigned int>(cell2DsErrorH1.size()),
                                   cell2DsErrorH1.data()
                                 }
                               });
          exporter.Export(exportFolder +
                          "/Solution.vtu");

          delete[] cell0Ds_exact_solution;
        }
      }
    }

    gedimData.Destroy();
  }
  // ***************************************************************************
  TEST(TestGeometry, Test_SquareLaplace_ImportedMesh)
  {
    const std::string exportFolder = "./Export/TestGeometry/Test_SquareLaplace_ImportedMesh";
    Gedim::Output::CreateFolder(exportFolder);

    GedimForPy::InterfaceConfiguration interfaceConfig;
    interfaceConfig.GeometricTolerance = 1.0e-8;

    GedimForPy::InterfaceData data;
    GedimForPy::InterfaceDataDAO gedimData(data);

    GedimForPy::GeDiM4Py_Logic interface;

    ASSERT_NO_THROW(interface.Initialize(interfaceConfig,
                                         data));

    const std::vector<std::string> meshImport = {
      "/home/geoscore/Desktop/GEO++/Courses/CppToPython/Meshes/Mesh1"
    };
    const unsigned int order = 2;

    for (unsigned int m = 0; m < meshImport.size(); m++)
    {
      GedimForPy::ImportMesh2D domain;
      domain.InputFolder = meshImport.at(m);
      domain.Separator = ';';

      GedimForPy::Domain2DMesh mesh = GedimForPy::GeDiM4Py_Logic::ImportDomainMesh2D(domain,
                                                                                     gedimData);

      Gedim::MeshMatricesDAO meshDAO(mesh.Mesh);

      // export
      {
        {
          Gedim::VTKUtilities exporter;
          gedimData.MeshUtilities().ExportMeshToVTU(meshDAO,
                                                    exportFolder,
                                                    "Mesh");
        }
      }

#if ACTIVE_CHECK == 1
      ASSERT_EQ(197, mesh.Mesh.NumberCell0D);
      ASSERT_EQ(537, mesh.Mesh.NumberCell1D);
      ASSERT_EQ(341, mesh.Mesh.NumberCell2D);
#endif

      GedimForPy::DiscreteSpace discreteSpace;
      discreteSpace.Type = GedimForPy::DiscreteSpace::Types::FEM;
      discreteSpace.Order = order;
      discreteSpace.BoundaryConditionsType = { GedimForPy::DiscreteSpace::BoundaryConditionTypes::None,
                                               GedimForPy::DiscreteSpace::BoundaryConditionTypes::Weak,
                                               GedimForPy::DiscreteSpace::BoundaryConditionTypes::Weak,
                                               GedimForPy::DiscreteSpace::BoundaryConditionTypes::Strong };

      GedimForPy::DiscreteProblemData problemData = GedimForPy::GeDiM4Py_Logic::Discretize(meshDAO,
                                                                                           mesh.MeshGeometricData,
                                                                                           discreteSpace);

#if ACTIVE_CHECK == 1
      ASSERT_EQ(711, problemData.NumberDOFs);
      ASSERT_EQ(23, problemData.NumberStrongs);
      ASSERT_EQ(197, problemData.Cell0Ds_DOF.size());
      ASSERT_EQ(537, problemData.Cell1Ds_DOF.size());
#endif

      // export
      {
        {
          std::vector<double> cell0Ds_DOFType(meshDAO.Cell0DTotalNumber(), 0.0);
          std::vector<double> cell0Ds_DOFGlobalIndex(meshDAO.Cell0DTotalNumber(), 0.0);
          std::vector<double> cell1Ds_DOFType(meshDAO.Cell1DTotalNumber(), 0.0);
          std::vector<double> cell1Ds_DOFGlobalIndex(meshDAO.Cell1DTotalNumber(), 0.0);

          for (unsigned int p = 0; p < problemData.Cell0Ds_DOF.size(); p++)
          {
            const GedimForPy::DiscreteProblemData::DOF& dof = problemData.Cell0Ds_DOF[p];
            cell0Ds_DOFType[p] = (unsigned int)dof.Type;
            cell0Ds_DOFGlobalIndex[p] = dof.Global_Index;
          }

          for (unsigned int p = 0; p < problemData.Cell1Ds_DOF.size(); p++)
          {
            const GedimForPy::DiscreteProblemData::DOF& dof = problemData.Cell1Ds_DOF[p];
            cell1Ds_DOFType[p] = (unsigned int)dof.Type;
            cell1Ds_DOFGlobalIndex[p] = dof.Global_Index;
          }

          Gedim::VTKUtilities exporter;
          exporter.AddSegments(meshDAO.Cell0DsCoordinates(),
                               meshDAO.Cell1DsExtremes(),
                               {
                                 {
                                   "cell0Ds_DOFType",
                                   Gedim::VTPProperty::Formats::Points,
                                   static_cast<unsigned int>(cell0Ds_DOFType.size()),
                                   cell0Ds_DOFType.data()
                                 },
                                 {
                                   "cell0Ds_DOFGlobalIndex",
                                   Gedim::VTPProperty::Formats::Points,
                                   static_cast<unsigned int>(cell0Ds_DOFGlobalIndex.size()),
                                   cell0Ds_DOFGlobalIndex.data()
                                 },
                                 {
                                   "cell1Ds_DOFType",
                                   Gedim::VTPProperty::Formats::Cells,
                                   static_cast<unsigned int>(cell1Ds_DOFType.size()),
                                   cell1Ds_DOFType.data()
                                 },
                                 {
                                   "cell1Ds_DOFGlobalIndex",
                                   Gedim::VTPProperty::Formats::Cells,
                                   static_cast<unsigned int>(cell1Ds_DOFGlobalIndex.size()),
                                   cell1Ds_DOFGlobalIndex.data()
                                 }
                               });
          exporter.Export(exportFolder +
                          "/DOFs.vtu");
        }
      }

      std::list<Eigen::Triplet<double>> stiffnessTriplets, stiffnessStrongTriplets;
      GedimForPy::GeDiM4Py_Logic::AssembleStiffnessMatrix(HeatConductivity::DiffusionTerm,
                                                          meshDAO,
                                                          mesh.Cell2DsMap,
                                                          problemData,
                                                          stiffnessTriplets,
                                                          stiffnessStrongTriplets);

      const Eigen::VectorXd weakTerm_Down = GedimForPy::GeDiM4Py_Logic::AssembleWeakTerm(HeatConductivity::WeakTerm_Down,
                                                                                         1,
                                                                                         meshDAO,
                                                                                         mesh.MeshGeometricData.Cell2DsVertices,
                                                                                         mesh.MeshGeometricData.Cell2DsEdgeLengths,
                                                                                         mesh.MeshGeometricData.Cell2DsEdgeTangents,
                                                                                         mesh.Cell2DsMap,
                                                                                         problemData);

      Eigen::SparseMatrix<double> stiffness(problemData.NumberDOFs,
                                            problemData.NumberDOFs);
      stiffness.setFromTriplets(stiffnessTriplets.begin(),
                                stiffnessTriplets.end());
      stiffness.makeCompressed();
      stiffnessTriplets.clear();

      Eigen::SparseLU<Eigen::SparseMatrix<double>> linearSolver;
      linearSolver.compute(stiffness);

      const Eigen::VectorXd solution = linearSolver.solve(weakTerm_Down);

      // export
      {
        GedimForPy::GeDiM4Py_Logic::ExportSolution(Poisson::ExactSolution,
                                                   solution,
                                                   Eigen::VectorXd::Zero(problemData.NumberStrongs),
                                                   meshDAO,
                                                   problemData,
                                                   {
                                                     exportFolder,
                                                     "Solution"
                                                   });

        {
          std::vector<double> cell0Ds_numeric_solution(meshDAO.Cell0DTotalNumber(),
                                                       0.0);
          const Eigen::MatrixXd coordinates = meshDAO.Cell0DsCoordinates();

          const double* cell0Ds_diffusion = HeatConductivity::DiffusionTerm(coordinates.cols(),
                                                                            coordinates.data());

          for (unsigned int p = 0; p < meshDAO.Cell0DTotalNumber(); p++)
          {
            const GedimForPy::DiscreteProblemData::DOF& dof = problemData.Cell0Ds_DOF[p];

            switch (dof.Type)
            {
              case GedimForPy::DiscreteProblemData::DOF::Types::DOF:
                cell0Ds_numeric_solution[p] = solution[dof.Global_Index];
                break;
              case GedimForPy::DiscreteProblemData::DOF::Types::Strong:
                cell0Ds_numeric_solution[p] = 0.0;
                break;
              default:
                throw std::runtime_error("DOF Type " +
                                         std::to_string((unsigned int)dof.Type) +
                                         " not supported");
            }
          }

          Gedim::VTKUtilities exporter;
          exporter.AddPolygons(meshDAO.Cell0DsCoordinates(),
                               meshDAO.Cell2DsVertices(),
                               {
                                 {
                                   "cell0Ds_diffusion",
                                   Gedim::VTPProperty::Formats::Points,
                                   static_cast<unsigned int>(coordinates.cols()),
                                   cell0Ds_diffusion
                                 },
                                 {
                                   "cell0Ds_numeric_solution",
                                   Gedim::VTPProperty::Formats::Points,
                                   static_cast<unsigned int>(cell0Ds_numeric_solution.size()),
                                   cell0Ds_numeric_solution.data()
                                 }
                               });
          exporter.Export(exportFolder +
                          "/Solution.vtu");

          delete[] cell0Ds_diffusion;
        }
      }
    }

    gedimData.Destroy();
  }
  // ***************************************************************************
  TEST(TestGeometry, Test_SquareStokes)
  {
    const std::string exportFolder = "./Export/TestGeometry/Test_SquareStokes";
    Gedim::Output::CreateFolder(exportFolder);

    GedimForPy::InterfaceConfiguration interfaceConfig;
    interfaceConfig.GeometricTolerance = 1.0e-8;

    GedimForPy::InterfaceData data;
    GedimForPy::InterfaceDataDAO gedimData(data);

    GedimForPy::GeDiM4Py_Logic interface;

    ASSERT_NO_THROW(interface.Initialize(interfaceConfig,
                                         data));

    const std::vector<double> meshSize = { 0.1 };

    for (unsigned int m = 0; m < meshSize.size(); m++)
    {
      GedimForPy::Domain2DMesh mesh;

      // create with triangle
      GedimForPy::Domain2D domain;
      domain.Vertices = gedimData.GeometryUtilities().CreateSquare(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                                   1.0);
      domain.VerticesBoundaryCondition = { 1, 1, 1, 1 };
      domain.EdgesBoundaryCondition = { 2, 3, 4, 5 };
      domain.DiscretizationType = GedimForPy::Domain2D::DiscretizationTypes::Triangular;
      domain.MeshCellsMaximumArea = meshSize[m];

      mesh = GedimForPy::GeDiM4Py_Logic::CreateDomainMesh2D(domain,
                                                            gedimData);

      // export
      {
        {
          Gedim::VTKUtilities exporter;
          exporter.AddPolygon(domain.Vertices);
          exporter.Export(exportFolder +
                          "/Domain.vtu");
        }
      }


      Gedim::MeshMatricesDAO meshDAO(mesh.Mesh);

      // export
      {
        {
          Gedim::VTKUtilities exporter;
          gedimData.MeshUtilities().ExportMeshToVTU(meshDAO,
                                                    exportFolder,
                                                    "Mesh");
        }
      }

      GedimForPy::DiscreteSpace speed_DiscreteSpace;
      speed_DiscreteSpace.Type = GedimForPy::DiscreteSpace::Types::FEM;
      speed_DiscreteSpace.Order = 2;
      speed_DiscreteSpace.BoundaryConditionsType = { GedimForPy::DiscreteSpace::BoundaryConditionTypes::None,
                                                     GedimForPy::DiscreteSpace::BoundaryConditionTypes::Strong,
                                                     GedimForPy::DiscreteSpace::BoundaryConditionTypes::Strong,
                                                     GedimForPy::DiscreteSpace::BoundaryConditionTypes::Strong,
                                                     GedimForPy::DiscreteSpace::BoundaryConditionTypes::Strong,
                                                     GedimForPy::DiscreteSpace::BoundaryConditionTypes::Strong };
      GedimForPy::DiscreteSpace pressure_DiscreteSpace;
      pressure_DiscreteSpace.Type = GedimForPy::DiscreteSpace::Types::FEM;
      pressure_DiscreteSpace.Order = 1;
      pressure_DiscreteSpace.BoundaryConditionsType = { GedimForPy::DiscreteSpace::BoundaryConditionTypes::None,
                                                        GedimForPy::DiscreteSpace::BoundaryConditionTypes::Strong,
                                                        GedimForPy::DiscreteSpace::BoundaryConditionTypes::None,
                                                        GedimForPy::DiscreteSpace::BoundaryConditionTypes::None,
                                                        GedimForPy::DiscreteSpace::BoundaryConditionTypes::None,
                                                        GedimForPy::DiscreteSpace::BoundaryConditionTypes::None };

      GedimForPy::DiscreteProblemData speed_problemData = GedimForPy::GeDiM4Py_Logic::Discretize(meshDAO,
                                                                                                 mesh.MeshGeometricData,
                                                                                                 speed_DiscreteSpace);
      GedimForPy::DiscreteProblemData pressure_problemData = GedimForPy::GeDiM4Py_Logic::Discretize(meshDAO,
                                                                                                    mesh.MeshGeometricData,
                                                                                                    pressure_DiscreteSpace);

      // export
      {
        {
          std::vector<double> cell0Ds_DOFType(meshDAO.Cell0DTotalNumber(), 0.0);
          std::vector<double> cell0Ds_DOFGlobalIndex(meshDAO.Cell0DTotalNumber(), 0.0);
          std::vector<double> cell1Ds_DOFType(meshDAO.Cell1DTotalNumber(), 0.0);
          std::vector<double> cell1Ds_DOFGlobalIndex(meshDAO.Cell1DTotalNumber(), 0.0);

          for (unsigned int p = 0; p < speed_problemData.Cell0Ds_DOF.size(); p++)
          {
            const GedimForPy::DiscreteProblemData::DOF& dof = speed_problemData.Cell0Ds_DOF[p];
            cell0Ds_DOFType[p] = (unsigned int)dof.Type;
            cell0Ds_DOFGlobalIndex[p] = dof.Global_Index;
          }

          for (unsigned int p = 0; p < speed_problemData.Cell1Ds_DOF.size(); p++)
          {
            const GedimForPy::DiscreteProblemData::DOF& dof = speed_problemData.Cell1Ds_DOF[p];
            cell1Ds_DOFType[p] = (unsigned int)dof.Type;
            cell1Ds_DOFGlobalIndex[p] = dof.Global_Index;
          }

          Gedim::VTKUtilities exporter;
          exporter.AddSegments(meshDAO.Cell0DsCoordinates(),
                               meshDAO.Cell1DsExtremes(),
                               {
                                 {
                                   "cell0Ds_DOFType",
                                   Gedim::VTPProperty::Formats::Points,
                                   static_cast<unsigned int>(cell0Ds_DOFType.size()),
                                   cell0Ds_DOFType.data()
                                 },
                                 {
                                   "cell0Ds_DOFGlobalIndex",
                                   Gedim::VTPProperty::Formats::Points,
                                   static_cast<unsigned int>(cell0Ds_DOFGlobalIndex.size()),
                                   cell0Ds_DOFGlobalIndex.data()
                                 },
                                 {
                                   "cell1Ds_DOFType",
                                   Gedim::VTPProperty::Formats::Cells,
                                   static_cast<unsigned int>(cell1Ds_DOFType.size()),
                                   cell1Ds_DOFType.data()
                                 },
                                 {
                                   "cell1Ds_DOFGlobalIndex",
                                   Gedim::VTPProperty::Formats::Cells,
                                   static_cast<unsigned int>(cell1Ds_DOFGlobalIndex.size()),
                                   cell1Ds_DOFGlobalIndex.data()
                                 }
                               });
          exporter.Export(exportFolder +
                          "/speed_DOFs.vtu");
        }
      }

      std::list<Eigen::Triplet<double>> stiffness_dx_Triplets, stiffnessStrong_dx_Triplets;
      GedimForPy::GeDiM4Py_Logic::AssembleStiffnessMatrix(Stokes::ViscosityTerm,
                                                          meshDAO,
                                                          mesh.Cell2DsMap,
                                                          speed_problemData,
                                                          stiffness_dx_Triplets,
                                                          stiffnessStrong_dx_Triplets);
      std::list<Eigen::Triplet<double>> stiffness_dy_Triplets, stiffnessStrong_dy_Triplets;
      GedimForPy::GeDiM4Py_Logic::AssembleStiffnessMatrix(Stokes::ViscosityTerm,
                                                          meshDAO,
                                                          mesh.Cell2DsMap,
                                                          speed_problemData,
                                                          stiffness_dy_Triplets,
                                                          stiffnessStrong_dy_Triplets);
      std::list<Eigen::Triplet<double>> advection_dx_Triplets, advectionStrong_dx_Triplets;
      GedimForPy::GeDiM4Py_Logic::AssembleAdvectionMatrix(Stokes::AdvectionTerm_1,
                                                          meshDAO,
                                                          mesh.Cell2DsMap,
                                                          speed_problemData,
                                                          pressure_problemData,
                                                          advection_dx_Triplets,
                                                          advectionStrong_dx_Triplets);
      std::list<Eigen::Triplet<double>> advection_dy_Triplets, advectionStrong_dy_Triplets;
      GedimForPy::GeDiM4Py_Logic::AssembleAdvectionMatrix(Stokes::AdvectionTerm_2,
                                                          meshDAO,
                                                          mesh.Cell2DsMap,
                                                          speed_problemData,
                                                          pressure_problemData,
                                                          advection_dy_Triplets,
                                                          advectionStrong_dy_Triplets);

      const Eigen::VectorXd forcingTerm_1 = GedimForPy::GeDiM4Py_Logic::AssembleForcingTerm(Stokes::ForcingTerm_1,
                                                                                            meshDAO,
                                                                                            mesh.Cell2DsMap,
                                                                                            speed_problemData);
      const Eigen::VectorXd forcingTerm_2 = GedimForPy::GeDiM4Py_Logic::AssembleForcingTerm(Stokes::ForcingTerm_2,
                                                                                            meshDAO,
                                                                                            mesh.Cell2DsMap,
                                                                                            speed_problemData);

      const Eigen::VectorXd pressure_solutionStrong = GedimForPy::GeDiM4Py_Logic::AssembleStrongSolution(Stokes::ExactPressureSolution,
                                                                                                         1,
                                                                                                         meshDAO,
                                                                                                         mesh.Cell2DsMap,
                                                                                                         pressure_problemData);

      std::list<Eigen::Triplet<double>> saddlePoint_Triplets;
      for (const Eigen::Triplet<double>& triplet : stiffness_dx_Triplets)
      {
        saddlePoint_Triplets.push_back(Eigen::Triplet<double>(triplet.row(),
                                                              triplet.col(),
                                                              triplet.value()));
      }
      stiffness_dx_Triplets.clear();
      for (const Eigen::Triplet<double>& triplet : stiffness_dy_Triplets)
      {
        saddlePoint_Triplets.push_back(Eigen::Triplet<double>(speed_problemData.NumberDOFs + triplet.row(),
                                                              speed_problemData.NumberDOFs + triplet.col(),
                                                              triplet.value()));
      }
      stiffness_dy_Triplets.clear();
      for (const Eigen::Triplet<double>& triplet : advection_dx_Triplets)
      {
        saddlePoint_Triplets.push_back(Eigen::Triplet<double>(triplet.col(),
                                                              2 * speed_problemData.NumberDOFs + triplet.row(),
                                                              - triplet.value()));
        saddlePoint_Triplets.push_back(Eigen::Triplet<double>(2 * speed_problemData.NumberDOFs + triplet.row(),
                                                              triplet.col(),
                                                              - triplet.value()));
      }
      advection_dx_Triplets.clear();
      for (const Eigen::Triplet<double>& triplet : advection_dy_Triplets)
      {
        saddlePoint_Triplets.push_back(Eigen::Triplet<double>(speed_problemData.NumberDOFs + triplet.col(),
                                                              2 * speed_problemData.NumberDOFs + triplet.row(),
                                                              - triplet.value()));
        saddlePoint_Triplets.push_back(Eigen::Triplet<double>(2 * speed_problemData.NumberDOFs + triplet.row(),
                                                              speed_problemData.NumberDOFs + triplet.col(),
                                                              - triplet.value()));
      }
      advection_dy_Triplets.clear();

      Eigen::SparseMatrix<double> saddlePoint_matrix(2 * speed_problemData.NumberDOFs +
                                                     pressure_problemData.NumberDOFs,
                                                     2 * speed_problemData.NumberDOFs +
                                                     pressure_problemData.NumberDOFs);
      saddlePoint_matrix.setFromTriplets(saddlePoint_Triplets.begin(),
                                         saddlePoint_Triplets.end());
      saddlePoint_matrix.makeCompressed();
      saddlePoint_Triplets.clear();

      Eigen::VectorXd saddlePoint_forcingTerm = Eigen::VectorXd::Zero(2 * speed_problemData.NumberDOFs +
                                                                      pressure_problemData.NumberDOFs);
      saddlePoint_forcingTerm.segment(0, speed_problemData.NumberDOFs) = forcingTerm_1;
      saddlePoint_forcingTerm.segment(speed_problemData.NumberDOFs, speed_problemData.NumberDOFs) = forcingTerm_2;

      Eigen::SparseLU<Eigen::SparseMatrix<double>> linearSolver;
      linearSolver.compute(saddlePoint_matrix);

      const Eigen::VectorXd solution = linearSolver.solve(saddlePoint_forcingTerm);

      const Eigen::VectorXd pressure_cell2DsErrorL2 = GedimForPy::GeDiM4Py_Logic::ComputeErrorL2(Stokes::ExactPressureSolution,
                                                                                                 solution.segment(2 * speed_problemData.NumberDOFs, pressure_problemData.NumberDOFs),
                                                                                                 pressure_solutionStrong,
                                                                                                 meshDAO,
                                                                                                 mesh.Cell2DsMap,
                                                                                                 pressure_problemData);


#if ACTIVE_CHECK == 0
      std::cerr.precision(16);
      std::cerr<< std::scientific<< "dofs"<< ","<< "h"<< ","<< "errorL2"<< std::endl;
      std::cerr<< std::scientific<< pressure_problemData.NumberDOFs<< ","<< pressure_problemData.H<< ","<< sqrt(pressure_cell2DsErrorL2.sum())<< std::endl;
#endif

      // export
      {
        {
          std::vector<double> pressure_cell0Ds_numeric_solution(meshDAO.Cell0DTotalNumber(),
                                                                0.0);
          std::vector<double> speed_1_cell0Ds_numeric_solution(meshDAO.Cell0DTotalNumber(),
                                                               0.0);
          std::vector<double> speed_2_cell0Ds_numeric_solution(meshDAO.Cell0DTotalNumber(),
                                                               0.0);

          const Eigen::MatrixXd coordinates = meshDAO.Cell0DsCoordinates();

          const double* pressure_cell0Ds_exact_solution = Stokes::ExactPressureSolution(coordinates.cols(),
                                                                                        coordinates.data());
          const double* speed_1_cell0Ds_exact_solution = Stokes::ExactSpeedSolution_1(coordinates.cols(),
                                                                                      coordinates.data());
          const double* speed_2_cell0Ds_exact_solution = Stokes::ExactSpeedSolution_2(coordinates.cols(),
                                                                                      coordinates.data());

          for (unsigned int p = 0; p < meshDAO.Cell0DTotalNumber(); p++)
          {
            const GedimForPy::DiscreteProblemData::DOF& speed_dof = speed_problemData.Cell0Ds_DOF[p];
            const GedimForPy::DiscreteProblemData::DOF& pressure_dof = pressure_problemData.Cell0Ds_DOF[p];

            switch (speed_dof.Type)
            {
              case GedimForPy::DiscreteProblemData::DOF::Types::DOF:
                speed_1_cell0Ds_numeric_solution[p] = solution[speed_dof.Global_Index];
                speed_2_cell0Ds_numeric_solution[p] = solution[speed_problemData.NumberDOFs + speed_dof.Global_Index];
                break;
              case GedimForPy::DiscreteProblemData::DOF::Types::Strong:
                speed_1_cell0Ds_numeric_solution[p] = 0.0;
                speed_2_cell0Ds_numeric_solution[p] = 0.0;
                break;
              default:
                throw std::runtime_error("DOF Type " +
                                         std::to_string((unsigned int)speed_dof.Type) +
                                         " not supported");
            }

            switch (pressure_dof.Type)
            {
              case GedimForPy::DiscreteProblemData::DOF::Types::DOF:
                pressure_cell0Ds_numeric_solution[p] = solution[2 * speed_problemData.NumberDOFs + pressure_dof.Global_Index];
                break;
              case GedimForPy::DiscreteProblemData::DOF::Types::Strong:
                pressure_cell0Ds_numeric_solution[p] = pressure_solutionStrong[pressure_dof.Global_Index];
                break;
              default:
                throw std::runtime_error("DOF Type " +
                                         std::to_string((unsigned int)pressure_dof.Type) +
                                         " not supported");
            }
          }

          GedimForPy::GeDiM4Py_Logic::ExportSolution(Stokes::ExactPressureSolution,
                                                     solution.segment(2 * speed_problemData.NumberDOFs,
                                                                      pressure_problemData.NumberDOFs),
                                                     pressure_solutionStrong,
                                                     meshDAO,
                                                     pressure_problemData,
                                                     {
                                                       exportFolder,
                                                       "Pressure"
                                                     });
          GedimForPy::GeDiM4Py_Logic::ExportSolution(Stokes::ExactSpeedSolution_1,
                                                     solution.segment(0,
                                                                      speed_problemData.NumberDOFs),
                                                     Eigen::VectorXd::Zero(speed_problemData.NumberStrongs),
                                                     meshDAO,
                                                     speed_problemData,
                                                     {
                                                       exportFolder,
                                                       "Speed_1"
                                                     });
          GedimForPy::GeDiM4Py_Logic::ExportSolution(Stokes::ExactSpeedSolution_2,
                                                     solution.segment(speed_problemData.NumberDOFs,
                                                                      speed_problemData.NumberDOFs),
                                                     Eigen::VectorXd::Zero(speed_problemData.NumberStrongs),
                                                     meshDAO,
                                                     speed_problemData,
                                                     {
                                                       exportFolder,
                                                       "Speed_2"
                                                     });

          Gedim::VTKUtilities exporter;
          exporter.AddPolygons(meshDAO.Cell0DsCoordinates(),
                               meshDAO.Cell2DsVertices(),
                               {
                                 {
                                   "pressure_cell0Ds_exact_solution",
                                   Gedim::VTPProperty::Formats::Points,
                                   static_cast<unsigned int>(coordinates.cols()),
                                   pressure_cell0Ds_exact_solution
                                 },
                                 {
                                   "pressure_cell0Ds_numeric_solution",
                                   Gedim::VTPProperty::Formats::Points,
                                   static_cast<unsigned int>(pressure_cell0Ds_numeric_solution.size()),
                                   pressure_cell0Ds_numeric_solution.data()
                                 },
                                 {
                                   "pressure_cell2Ds_errorL2",
                                   Gedim::VTPProperty::Formats::Cells,
                                   static_cast<unsigned int>(pressure_cell2DsErrorL2.size()),
                                   pressure_cell2DsErrorL2.data()
                                 },
                                 {
                                   "speed_1_cell0Ds_exact_solution",
                                   Gedim::VTPProperty::Formats::Points,
                                   static_cast<unsigned int>(coordinates.cols()),
                                   speed_1_cell0Ds_exact_solution
                                 },
                                 {
                                   "speed_1_cell0Ds_numeric_solution",
                                   Gedim::VTPProperty::Formats::Points,
                                   static_cast<unsigned int>(speed_1_cell0Ds_numeric_solution.size()),
                                   speed_1_cell0Ds_numeric_solution.data()
                                 },
                                 {
                                   "speed_2_cell0Ds_exact_solution",
                                   Gedim::VTPProperty::Formats::Points,
                                   static_cast<unsigned int>(coordinates.cols()),
                                   speed_2_cell0Ds_exact_solution
                                 },
                                 {
                                   "speed_2_cell0Ds_numeric_solution",
                                   Gedim::VTPProperty::Formats::Points,
                                   static_cast<unsigned int>(speed_2_cell0Ds_numeric_solution.size()),
                                   speed_2_cell0Ds_numeric_solution.data()
                                 }
                               });
          exporter.Export(exportFolder +
                          "/solution.vtu");

          delete[] pressure_cell0Ds_exact_solution;
        }
      }
    }

    gedimData.Destroy();
  }
  // ***************************************************************************
  TEST(TestGeometry, Test_SquareBurger)
  {
    const std::string exportFolder = "./Export/TestGeometry/Test_SquareBurger";
    Gedim::Output::CreateFolder(exportFolder);
    const std::string exportSolutionFolder = exportFolder + "/Solution";
    Gedim::Output::CreateFolder(exportSolutionFolder);

    GedimForPy::InterfaceConfiguration interfaceConfig;
    interfaceConfig.GeometricTolerance = 1.0e-8;

    GedimForPy::InterfaceData data;
    GedimForPy::InterfaceDataDAO gedimData(data);

    GedimForPy::GeDiM4Py_Logic interface;

    ASSERT_NO_THROW(interface.Initialize(interfaceConfig,
                                         data));

    const std::vector<double> meshSize = { 0.01 };
    const unsigned int order = 1;

    for (unsigned int m = 0; m < meshSize.size(); m++)
    {
      GedimForPy::Domain2D domain;
      domain.Vertices = gedimData.GeometryUtilities().CreateSquare(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                                   1.0);
      domain.VerticesBoundaryCondition = { 1, 1, 1, 1 };
      domain.EdgesBoundaryCondition = { 1, 1, 1, 1 };
      domain.DiscretizationType = GedimForPy::Domain2D::DiscretizationTypes::Triangular;
      domain.MeshCellsMaximumArea = meshSize[m];

      GedimForPy::Domain2DMesh mesh = GedimForPy::GeDiM4Py_Logic::CreateDomainMesh2D(domain,
                                                                                     gedimData);
      Gedim::MeshMatricesDAO meshDAO(mesh.Mesh);

      // export
      {
        {
          Gedim::VTKUtilities exporter;
          exporter.AddPolygon(domain.Vertices);
          exporter.Export(exportFolder +
                          "/Domain.vtu");
        }

        {
          Gedim::VTKUtilities exporter;
          gedimData.MeshUtilities().ExportMeshToVTU(meshDAO,
                                                    exportFolder,
                                                    "Mesh");
        }
      }

      GedimForPy::DiscreteSpace discreteSpace;
      discreteSpace.Type = GedimForPy::DiscreteSpace::Types::FEM;
      discreteSpace.Order = order;
      discreteSpace.BoundaryConditionsType = {
        GedimForPy::DiscreteSpace::BoundaryConditionTypes::None,
        GedimForPy::DiscreteSpace::BoundaryConditionTypes::Strong
      };

      GedimForPy::DiscreteProblemData problemData = GedimForPy::GeDiM4Py_Logic::Discretize(meshDAO,
                                                                                           mesh.MeshGeometricData,
                                                                                           discreteSpace);
      // export
      {
        {
          std::vector<double> cell0Ds_DOFType(meshDAO.Cell0DTotalNumber(), 0.0);
          std::vector<double> cell0Ds_DOFGlobalIndex(meshDAO.Cell0DTotalNumber(), 0.0);
          std::vector<double> cell1Ds_DOFType(meshDAO.Cell1DTotalNumber(), 0.0);
          std::vector<double> cell1Ds_DOFGlobalIndex(meshDAO.Cell1DTotalNumber(), 0.0);

          for (unsigned int p = 0; p < problemData.Cell0Ds_DOF.size(); p++)
          {
            const GedimForPy::DiscreteProblemData::DOF& dof = problemData.Cell0Ds_DOF[p];
            cell0Ds_DOFType[p] = (unsigned int)dof.Type;
            cell0Ds_DOFGlobalIndex[p] = dof.Global_Index;
          }

          for (unsigned int p = 0; p < problemData.Cell1Ds_DOF.size(); p++)
          {
            const GedimForPy::DiscreteProblemData::DOF& dof = problemData.Cell1Ds_DOF[p];
            cell1Ds_DOFType[p] = (unsigned int)dof.Type;
            cell1Ds_DOFGlobalIndex[p] = dof.Global_Index;
          }

          Gedim::VTKUtilities exporter;
          exporter.AddSegments(meshDAO.Cell0DsCoordinates(),
                               meshDAO.Cell1DsExtremes(),
                               {
                                 {
                                   "cell0Ds_DOFType",
                                   Gedim::VTPProperty::Formats::Points,
                                   static_cast<unsigned int>(cell0Ds_DOFType.size()),
                                   cell0Ds_DOFType.data()
                                 },
                                 {
                                   "cell0Ds_DOFGlobalIndex",
                                   Gedim::VTPProperty::Formats::Points,
                                   static_cast<unsigned int>(cell0Ds_DOFGlobalIndex.size()),
                                   cell0Ds_DOFGlobalIndex.data()
                                 },
                                 {
                                   "cell1Ds_DOFType",
                                   Gedim::VTPProperty::Formats::Cells,
                                   static_cast<unsigned int>(cell1Ds_DOFType.size()),
                                   cell1Ds_DOFType.data()
                                 },
                                 {
                                   "cell1Ds_DOFGlobalIndex",
                                   Gedim::VTPProperty::Formats::Cells,
                                   static_cast<unsigned int>(cell1Ds_DOFGlobalIndex.size()),
                                   cell1Ds_DOFGlobalIndex.data()
                                 }
                               });
          exporter.Export(exportFolder +
                          "/DOFs.vtu");
        }
      }

      Eigen::VectorXd u_k = Eigen::VectorXd::Zero(problemData.NumberDOFs);

      double residual_norm = 1.0, solution_norm = 1.0;
      const double newton_tol = 1.0e-6;
      const unsigned int max_iterations = 20;
      int num_iteration = 1;

      const Eigen::VectorXd u_strong = Eigen::VectorXd::Zero(problemData.NumberStrongs);

      while (residual_norm > newton_tol * solution_norm &&
             num_iteration < max_iterations)
      {
        std::list<Eigen::Triplet<double>> J_stiffnessTriplets, J_stiffnessStrongTriplets;
        GedimForPy::GeDiM4Py_Logic::AssembleStiffnessMatrix(Burger::DiffusionTerm,
                                                            meshDAO,
                                                            mesh.Cell2DsMap,
                                                            problemData,
                                                            J_stiffnessTriplets,
                                                            J_stiffnessStrongTriplets);
        std::list<Eigen::Triplet<double>> J_reactionTriplets, J_reactionStrongTriplets;
        GedimForPy::GeDiM4Py_Logic::AssembleNonLinearReactionMatrix(Burger::ReactionTerm,
                                                                    Burger::NonLinear_Reaction,
                                                                    meshDAO,
                                                                    mesh.Cell2DsMap,
                                                                    problemData,
                                                                    u_k,
                                                                    u_strong,
                                                                    J_reactionTriplets,
                                                                    J_reactionStrongTriplets);
        std::list<Eigen::Triplet<double>> J_advectionTriplets, J_advectionStrongTriplets;
        GedimForPy::GeDiM4Py_Logic::AssembleNonLinearAdvectionMatrix(Burger::AdvectionTerm,
                                                                     Burger::NonLinear_Advection,
                                                                     meshDAO,
                                                                     mesh.Cell2DsMap,
                                                                     problemData,
                                                                     problemData,
                                                                     u_k,
                                                                     u_strong,
                                                                     J_advectionTriplets,
                                                                     J_advectionStrongTriplets);

        const Eigen::VectorXd J_forcingTerm_g = GedimForPy::GeDiM4Py_Logic::AssembleForcingTerm(Burger::ForcingTerm,
                                                                                                meshDAO,
                                                                                                mesh.Cell2DsMap,
                                                                                                problemData);
        const Eigen::VectorXd J_forcingTerm_der_v = GedimForPy::GeDiM4Py_Logic::AssembleNonLinearDerivativeForcingTerm(Burger::OnesDerivative,
                                                                                                                       Burger::NonLinear_f_der_v,
                                                                                                                       meshDAO,
                                                                                                                       mesh.Cell2DsMap,
                                                                                                                       problemData,
                                                                                                                       u_k,
                                                                                                                       u_strong);
        const Eigen::VectorXd J_forcingTerm_v = GedimForPy::GeDiM4Py_Logic::AssembleNonLinearForcingTerm(Burger::Ones,
                                                                                                         Burger::NonLinear_f_v,
                                                                                                         meshDAO,
                                                                                                         mesh.Cell2DsMap,
                                                                                                         problemData,
                                                                                                         u_k,
                                                                                                         u_strong);


        Eigen::SparseMatrix<double> J_stiffness(problemData.NumberDOFs,
                                                problemData.NumberDOFs);
        J_stiffness.setFromTriplets(J_stiffnessTriplets.begin(),
                                    J_stiffnessTriplets.end());
        J_stiffness.makeCompressed();
        J_stiffnessTriplets.clear();

        Eigen::SparseMatrix<double> J_reaction(problemData.NumberDOFs,
                                               problemData.NumberDOFs);
        J_reaction.setFromTriplets(J_reactionTriplets.begin(),
                                   J_reactionTriplets.end());
        J_reaction.makeCompressed();
        J_reactionTriplets.clear();

        Eigen::SparseMatrix<double> J_advection(problemData.NumberDOFs,
                                                problemData.NumberDOFs);
        J_advection.setFromTriplets(J_advectionTriplets.begin(),
                                    J_advectionTriplets.end());
        J_advection.makeCompressed();
        J_advectionTriplets.clear();


        Eigen::SparseLU<Eigen::SparseMatrix<double>> linearSolver;
        linearSolver.compute(J_stiffness +
                             J_reaction +
                             J_advection);

        const Eigen::VectorXd du = linearSolver.solve(J_forcingTerm_g -
                                                      J_forcingTerm_der_v -
                                                      J_forcingTerm_v);
        u_k = u_k + du;

        const Eigen::VectorXd cell2DsErrorL2 = GedimForPy::GeDiM4Py_Logic::ComputeErrorL2(Burger::ExactSolution,
                                                                                          u_k,
                                                                                          u_strong,
                                                                                          meshDAO,
                                                                                          mesh.Cell2DsMap,
                                                                                          problemData);
        const Eigen::VectorXd cell2DsErrorH1 = GedimForPy::GeDiM4Py_Logic::ComputeErrorH1(Burger::ExactDerivativeSolution,
                                                                                          u_k,
                                                                                          u_strong,
                                                                                          meshDAO,
                                                                                          mesh.Cell2DsMap,
                                                                                          problemData);
        const Eigen::VectorXd cell2DsduNormL2 = GedimForPy::GeDiM4Py_Logic::ComputeErrorL2(Burger::ZeroSolution,
                                                                                           du,
                                                                                           Eigen::VectorXd::Zero(problemData.NumberStrongs),
                                                                                           meshDAO,
                                                                                           mesh.Cell2DsMap,
                                                                                           problemData);
        const Eigen::VectorXd cell2DsNormL2 = GedimForPy::GeDiM4Py_Logic::ComputeErrorL2(Burger::ZeroSolution,
                                                                                         u_k,
                                                                                         u_strong,
                                                                                         meshDAO,
                                                                                         mesh.Cell2DsMap,
                                                                                         problemData);
        const Eigen::VectorXd cell2DsNormH1 = GedimForPy::GeDiM4Py_Logic::ComputeErrorH1(Burger::ZeroDerivativeSolution,
                                                                                         u_k,
                                                                                         u_strong,
                                                                                         meshDAO,
                                                                                         mesh.Cell2DsMap,
                                                                                         problemData);
        solution_norm = std::sqrt(cell2DsNormL2.sum());
        residual_norm = std::sqrt(cell2DsduNormL2.sum());

        // export
        {
          GedimForPy::GeDiM4Py_Logic::ExportSolution(Burger::ExactSolution,
                                                     u_k,
                                                     u_strong,
                                                     meshDAO,
                                                     problemData,
                                                     {
                                                       exportSolutionFolder,
                                                       "Solution"
                                                     });

          {
            std::vector<double> cell0Ds_numeric_solution(meshDAO.Cell0DTotalNumber(),
                                                         0.0);
            const Eigen::MatrixXd coordinates = meshDAO.Cell0DsCoordinates();

            const double* cell0Ds_exact_solution = Burger::ExactSolution(coordinates.cols(),
                                                                         coordinates.data());

            for (unsigned int p = 0; p < meshDAO.Cell0DTotalNumber(); p++)
            {
              const GedimForPy::DiscreteProblemData::DOF& dof = problemData.Cell0Ds_DOF[p];

              switch (dof.Type)
              {
                case GedimForPy::DiscreteProblemData::DOF::Types::DOF:
                  cell0Ds_numeric_solution[p] = u_k[dof.Global_Index];
                  break;
                case GedimForPy::DiscreteProblemData::DOF::Types::Strong:
                  cell0Ds_numeric_solution[p] = u_strong[dof.Global_Index];
                  break;
                default:
                  throw std::runtime_error("DOF Type " +
                                           std::to_string((unsigned int)dof.Type) +
                                           " not supported");
              }
            }

            Gedim::VTKUtilities exporter;
            exporter.AddPolygons(meshDAO.Cell0DsCoordinates(),
                                 meshDAO.Cell2DsVertices(),
                                 {
                                   {
                                     "cell0Ds_exact_solution",
                                     Gedim::VTPProperty::Formats::Points,
                                     static_cast<unsigned int>(coordinates.cols()),
                                     cell0Ds_exact_solution
                                   },
                                   {
                                     "cell0Ds_numeric_solution",
                                     Gedim::VTPProperty::Formats::Points,
                                     static_cast<unsigned int>(cell0Ds_numeric_solution.size()),
                                     cell0Ds_numeric_solution.data()
                                   },
                                   {
                                     "cell2Ds_errorL2",
                                     Gedim::VTPProperty::Formats::Cells,
                                     static_cast<unsigned int>(cell2DsErrorL2.size()),
                                     cell2DsErrorL2.data()
                                   },
                                   {
                                     "cell2Ds_errorH1",
                                     Gedim::VTPProperty::Formats::Cells,
                                     static_cast<unsigned int>(cell2DsErrorH1.size()),
                                     cell2DsErrorH1.data()
                                   }
                                 });
            exporter.Export(exportSolutionFolder +
                            "/Solution_" + std::to_string(num_iteration) + ".vtu");

            delete[] cell0Ds_exact_solution;
          }
        }


#if ACTIVE_CHECK == 0
        std::cerr.precision(3);
        std::cerr<< std::scientific<< "dofs"<< ","
                 << "h"<< ","<< "errorL2"<< ","
                 << "errorH1"<< "normL2"<< ","
                 << "normH1"<< std::endl;
        std::cerr<< std::scientific<< problemData.NumberDOFs<< ","
                 << problemData.H<< ","<< sqrt(cell2DsErrorL2.sum())<< ","
                 << sqrt(cell2DsErrorH1.sum())<< ","
                 << sqrt(cell2DsNormL2.sum())<< ","
                 << sqrt(cell2DsNormH1.sum())
                 << std::endl;

        std::cout.precision(3);
        std::cout<< std::scientific<<
                    " Newton it "<< num_iteration<< " / "<< max_iterations<<
                    " residual "<< residual_norm<< " / "<< newton_tol * solution_norm<< std::endl;
#endif

        num_iteration++;
      }
    }

    gedimData.Destroy();
  }
}

#endif // __TEST_GEOMETRY_H
