#ifndef __TEST_HEAT_CONDUCTIVITY_H
#define __TEST_HEAT_CONDUCTIVITY_H

#include "test_utilities.hpp"

namespace UnitTesting
{
  class HeatConductivity final
  {
    public:
      static constexpr double R() { return 0.5; }
      static constexpr double K() { return 6.68; }
      static constexpr double G() { return 0.94; }
      // ***************************************************************************
      static double* DiffusionTerm(const int numPoints, const double* points)
      {
        double* values = new double[numPoints];

        Eigen::Map<const Eigen::MatrixXd> matPoints(points, 3, numPoints);
        Eigen::Map<Eigen::VectorXd> matValues(values, numPoints);
        matValues.setOnes();

        Eigen::Map<Eigen::VectorXd> vecValues(values, numPoints);

        for (unsigned int p = 0; p < numPoints; p++)
        {
          if (matPoints(0, p) * matPoints(0, p) + matPoints(1, p) * matPoints(1, p) <= (R() * R() + 1.0e-16))
            matValues[p] = K();
        }

        return values;
      }
      // ***************************************************************************
      static double* WeakTerm_Down(const int numPoints, const double* points)
      {
        double* values = new double[numPoints];

        Eigen::Map<Eigen::VectorXd> vecValues(values, numPoints);
        vecValues.setConstant(G());

        return values;
      }
      // ***************************************************************************
      static double* ExactSolution(const int numPoints, const double* points)
      {
        double* values = new double[numPoints];

        Eigen::Map<const Eigen::MatrixXd> matPoints(points, 3, numPoints);
        Eigen::Map<Eigen::VectorXd> vecValues(values, numPoints);

        vecValues = 16.0 * (matPoints.row(1).array() * (1.0 - matPoints.row(1).array()) *
                            matPoints.row(0).array() * (1.0 - matPoints.row(0).array())) + 1.1;
        return values;
      }
      // ***************************************************************************
  };
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
        GedimForPy::GeDiM4Py_Logic::ExportSolution(HeatConductivity::ExactSolution,
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
}

#endif // __TEST_HEAT_CONDUCTIVITY_H
