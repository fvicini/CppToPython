#ifndef __TEST_GEOMETRY_H
#define __TEST_GEOMETRY_H

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>

#include "GeDiM4Py_Logic.hpp"
#include "MeshMatricesDAO.hpp"
#include "VTKUtilities.hpp"
#include "test_Poisson.hpp"

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

    GedimForPy::Domain2D domain;
    domain.Vertices = gedimData.GeometryUtilities().CreateSquare(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                                 1.0);
    domain.VerticesBoundaryCondition = std::vector<unsigned int> { 1, 1, 1, 1 };
    domain.EdgesBoundaryCondition = std::vector<unsigned int> { 1, 1, 1, 1 };
    domain.DiscretizationType = GedimForPy::Domain2D::DiscretizationTypes::Triangular;
    domain.MeshCellsMaximumArea = 0.1;

    GedimForPy::Domain2DMesh mesh = GedimForPy::GeDiM4Py_Logic::CreateDomainMesh2D(domain,
                                                                                   gedimData);

    ASSERT_EQ(13, mesh.Mesh.NumberCell0D);
    ASSERT_EQ(28, mesh.Mesh.NumberCell1D);
    ASSERT_EQ(16, mesh.Mesh.NumberCell2D);

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
    discreteSpace.Order = 2;
    discreteSpace.BoundaryConditionsType = { GedimForPy::DiscreteSpace::BoundaryConditionTypes::None,
                                             GedimForPy::DiscreteSpace::BoundaryConditionTypes::Strong };

    GedimForPy::DiscreteProblemData problemData = GedimForPy::GeDiM4Py_Logic::Discretize(meshDAO,
                                                                                         discreteSpace);

    ASSERT_EQ(25, problemData.NumberDOFs);
    ASSERT_EQ(16, problemData.NumberStrongs);
    ASSERT_EQ(13, problemData.Cell0Ds_DOF.size());
    ASSERT_EQ(28, problemData.Cell1Ds_DOF.size());

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

    std::list<Eigen::Triplet<double>> stiffnessTriplets = GedimForPy::GeDiM4Py_Logic::AssembleStiffnessMatrix(Poisson::DiffusionTerm,
                                                                                                              meshDAO,
                                                                                                              mesh.Cell2DsMap,
                                                                                                              problemData);
    const Eigen::VectorXd forcingTerm = GedimForPy::GeDiM4Py_Logic::AssembleForcingTerm(Poisson::ForcingTerm,
                                                                                        meshDAO,
                                                                                        mesh.Cell2DsMap,
                                                                                        problemData);

    Eigen::SparseMatrix<double> stiffness(problemData.NumberDOFs,
                                          problemData.NumberDOFs);
    stiffness.setFromTriplets(stiffnessTriplets.begin(),
                              stiffnessTriplets.end());
    stiffness.makeCompressed();
    stiffnessTriplets.clear();
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>, Eigen::Lower> linearSolver;
    linearSolver.compute(stiffness);

    Eigen::VectorXd solution = linearSolver.solve(forcingTerm);

    std::cerr.precision(16);
    std::cerr<< std::scientific<< solution.transpose()<< std::endl;

    // export
    {
      {
        std::vector<double> cell0Ds_numeric_solution(meshDAO.Cell0DTotalNumber(),
                                                     0.0);
        Eigen::VectorXd cell0Ds_exact_solution = Poisson::ExactSolution(meshDAO.Cell0DsCoordinates());

        for (unsigned int p = 0; p < meshDAO.Cell0DTotalNumber(); p++)
        {
          const GedimForPy::DiscreteProblemData::DOF& dof = problemData.Cell0Ds_DOF[p];

          if (dof.Type != GedimForPy::DiscreteProblemData::DOF::Types::DOF)
            continue;

          cell0Ds_numeric_solution[p] = solution[dof.Global_Index];
        }

        Gedim::VTKUtilities exporter;
        exporter.AddPolygons(meshDAO.Cell0DsCoordinates(),
                             meshDAO.Cell2DsVertices(),
                             {
                               {
                                 "cell0Ds_exact_solution",
                                 Gedim::VTPProperty::Formats::Points,
                                 static_cast<unsigned int>(cell0Ds_exact_solution.size()),
                                 cell0Ds_exact_solution.data()
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
      }
    }

  }
  // ***************************************************************************
}

#endif // __TEST_GEOMETRY_H
