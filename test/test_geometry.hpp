#ifndef __TEST_GEOMETRY_H
#define __TEST_GEOMETRY_H

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>

#include "GeDiM4Py_Logic.hpp"
#include "MeshMatricesDAO.hpp"
#include "VTKUtilities.hpp"
#include "test_Poisson.hpp"
#include "test_heat_conductivity.hpp"

#define ACTIVE_CHECK 1

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
}

#endif // __TEST_GEOMETRY_H
