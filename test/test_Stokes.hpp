#ifndef __TEST_STOKES_H
#define __TEST_STOKES_H

#include "test_utilities.hpp"

namespace UnitTesting
{
  class Stokes final
  {
    public:
      static constexpr double v() { return 1.0; }

      // ***************************************************************************
      static double* ViscosityTerm(const int numPoints, const double* points)
      {
        double* values = new double[numPoints];

        Eigen::Map<Eigen::VectorXd> matValues(values, numPoints);
        matValues.setConstant(v());

        return values;
      }
      // ***************************************************************************
      static double* AdvectionTerm_1(const int numPoints, const double* points)
      {
        double* values = new double[2 * numPoints];

        Eigen::Map<Eigen::MatrixXd> matValues(values, 2, numPoints);
        matValues.row(0).setOnes();
        matValues.row(1).setZero();

        return values;
      }
      // ***************************************************************************
      static double* AdvectionTerm_2(const int numPoints, const double* points)
      {
        double* values = new double[2 * numPoints];

        Eigen::Map<Eigen::MatrixXd> matValues(values, 2, numPoints);
        matValues.row(0).setZero();
        matValues.row(1).setOnes();

        return values;
      }
      // ***************************************************************************
      static double* ForcingTerm_1(const int numPoints, const double* points)
      {
        double* values = new double[numPoints];

        Eigen::Map<const Eigen::MatrixXd> matPoints(points, 3, numPoints);
        Eigen::Map<Eigen::VectorXd> vecValues(values, numPoints);
        vecValues = - (+ 8.0 * M_PI * M_PI * cos(4.0 * M_PI * matPoints.row(0).array()) -
                     4.0 * M_PI * M_PI) *
                    sin(2.0 * M_PI * matPoints.row(1).array()) *
                    cos(2.0 * M_PI * matPoints.row(1).array()) +
                    (+ 2.0 * M_PI *
                     cos(2.0 * M_PI * matPoints.row(0).array()) *
                     cos(2.0 * M_PI * matPoints.row(1).array()));

        return values;
      }
      // ***************************************************************************
      static double* ForcingTerm_2(const int numPoints, const double* points)
      {
        double* values = new double[numPoints];

        Eigen::Map<const Eigen::MatrixXd> matPoints(points, 3, numPoints);
        Eigen::Map<Eigen::VectorXd> vecValues(values, numPoints);
        vecValues = - (- 8.0 * M_PI * M_PI * cos(4.0 * M_PI * matPoints.row(1).array()) +
                     4.0 * M_PI * M_PI) *
                    sin(2.0 * M_PI * matPoints.row(0).array()) *
                    cos(2.0 * M_PI * matPoints.row(0).array()) +
                    (- 2.0 * M_PI *
                     sin(2.0 * M_PI * matPoints.row(0).array()) *
                     sin(2.0 * M_PI * matPoints.row(1).array()));

        return values;
      }
      // ***************************************************************************
      static double* ExactPressureSolution(const int numPoints,
                                           const double* points)
      {
        double* values = new double[numPoints];

        Eigen::Map<const Eigen::MatrixXd> matPoints(points, 3, numPoints);
        Eigen::Map<Eigen::VectorXd> vecValues(values, numPoints);

        vecValues = sin(2.0 * M_PI * matPoints.row(0).array()) *
                    cos(2.0 * M_PI * matPoints.row(1).array());
        return values;
      }
      // ***************************************************************************
      static double* ExactSpeedSolution_1(const int numPoints,
                                          const double* points)
      {
        double* values = new double[numPoints];

        Eigen::Map<const Eigen::MatrixXd> matPoints(points, 3, numPoints);
        Eigen::Map<Eigen::VectorXd> vecValues(values, numPoints);

        vecValues = +0.5 *
                    sin(2.0 * M_PI * matPoints.row(0).array()) *
                    sin(2.0 * M_PI * matPoints.row(0).array()) *
                    sin(2.0 * M_PI * matPoints.row(1).array()) *
                    cos(2.0 * M_PI * matPoints.row(1).array());

        return values;
      }
      // ***************************************************************************
      static double* ExactSpeedSolution_2(const int numPoints,
                                          const double* points)
      {
        double* values = new double[numPoints];

        Eigen::Map<const Eigen::MatrixXd> matPoints(points, 3, numPoints);
        Eigen::Map<Eigen::VectorXd> vecValues(values, numPoints);

        vecValues = -0.5 *
                    sin(2.0 * M_PI * matPoints.row(1).array()) *
                    sin(2.0 * M_PI * matPoints.row(1).array()) *
                    sin(2.0 * M_PI * matPoints.row(0).array()) *
                    cos(2.0 * M_PI * matPoints.row(0).array());
        return values;
      }
      // ***************************************************************************
      static double* ExactPressureDerivativeSolution(const int direction,
                                                     const int numPoints,
                                                     const double* points)
      {
        double* values = new double[numPoints];

        Eigen::Map<const Eigen::MatrixXd> matPoints(points, 3, numPoints);
        Eigen::Map<Eigen::VectorXd> vecValues(values, numPoints);

        if(direction == 0)
          vecValues = +2.0 * M_PI *
                      cos(2.0 * M_PI * matPoints.row(0).array()) *
                      cos(2.0 * M_PI * matPoints.row(1).array());
        else if (direction == 1)
          vecValues = -2.0 * M_PI *
                      sin(2.0 * M_PI * matPoints.row(0).array()) *
                      sin(2.0 * M_PI * matPoints.row(1).array());
        else if (direction == 2)
          vecValues.setZero();
        else
          throw std::runtime_error("Error on direction");

        return values;
      }
      // ***************************************************************************
  };
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
}

#endif // __TEST_STOKES_H
