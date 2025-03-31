#ifndef __TEST_NAVIER_STOKES_H
#define __TEST_NAVIER_STOKES_H

#include "test_utilities.hpp"

namespace UnitTesting
{
  // ***************************************************************************
  class NavierStokes final
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
      static int A_transform_triplet_row(const int row, const int col)
      { return row; }
      static int A_transform_triplet_col(const int row, const int col)
      { return col; }
      static double A_transform_triplet_value(const double value)
      { return value; }
      // ***************************************************************************
      static int BT_transform_triplet_row(const int row, const int col)
      { return col; }
      static int BT_transform_triplet_col(const int row, const int col)
      { return row; }
      static double BT_transform_triplet_value(const double value)
      { return -1.0 * value; }
      // ***************************************************************************
      static int B_transform_triplet_row(const int row, const int col)
      { return row; }
      static int B_transform_triplet_col(const int row, const int col)
      { return col; }
      static double B_transform_triplet_value(const double value)
      { return -1.0 * value; }
      // ***************************************************************************
      static double* ZeroSolution(const int numPoints, const double* points)
      {
        double* values = new double[numPoints];

        Eigen::Map<Eigen::VectorXd> vecValues(values, numPoints);
        vecValues.setZero();

        return values;
      }
      // ***************************************************************************
      static double* ZeroDerivativeSolution(const int direction,
                                            const int numPoints,
                                            const double* points)
      {
        double* values = new double[numPoints];

        Eigen::Map<Eigen::VectorXd> vecValues(values, numPoints);
        vecValues.setZero();

        return values;
      }
      // ***************************************************************************
      static double* Ones(const int numPoints,
                          const double* points)
      {
        double* values = new double[numPoints];

        Eigen::Map<Eigen::VectorXd> matValues(values, numPoints);
        matValues.setConstant(1.0);

        return values;
      }
      // ***************************************************************************
      static double* OnesDerivative(const int numPoints,
                                    const double* points)
      {
        double* values = new double[2 * numPoints];

        Eigen::Map<Eigen::MatrixXd> matValues(values,
                                              2,
                                              numPoints);
        matValues.setOnes();

        return values;
      }
      // ***************************************************************************
      static double* NonLinear_double_dot_product(const int numPoints,
                                                  const double* points,
                                                  const double* u,
                                                  const double* u_x,
                                                  const double* u_y)
      {
        double* values = new double[2 * numPoints];

        Eigen::Map<Eigen::MatrixXd> matValues(values,
                                              2,
                                              numPoints);
        matValues.row(0)<< Eigen::Map<const Eigen::VectorXd>(u_x,
                                                             numPoints).transpose();
        matValues.row(1)<< Eigen::Map<const Eigen::VectorXd>(u_y,
                                                             numPoints).transpose();

        return values;
      }
      // ***************************************************************************
  };
  // ***************************************************************************
  class NavierStokes_T1 final
  {
    public:
      // ***************************************************************************
      static double* ForcingTerm_1(const int numPoints, const double* points)
      {
        double* values = new double[numPoints];

        Eigen::Map<Eigen::VectorXd> vecValues(values, numPoints);
        vecValues.setConstant(1.0);

        return values;
      }
      // ***************************************************************************
      static double* ForcingTerm_2(const int numPoints, const double* points)
      {
        double* values = new double[numPoints];

        Eigen::Map<Eigen::VectorXd> vecValues(values, numPoints);
        vecValues.setConstant(1.0);

        return values;
      }
      // ***************************************************************************
      static double* ExactPressureSolution(const int numPoints,
                                           const double* points)
      {
        double* values = new double[numPoints];

        Eigen::Map<const Eigen::MatrixXd> matPoints(points, 3, numPoints);
        Eigen::Map<Eigen::VectorXd> vecValues(values, numPoints);

        vecValues = matPoints.row(0).array() +
                    matPoints.row(1).array() -
                    1.0;
        return values;
      }
      // ***************************************************************************
      static double* ExactSpeedSolution_1(const int numPoints,
                                          const double* points)
      {
        double* values = new double[numPoints];

        Eigen::Map<Eigen::VectorXd> vecValues(values, numPoints);
        vecValues.setZero();

        return values;
      }
      // ***************************************************************************
      static double* ExactSpeedSolution_2(const int numPoints,
                                          const double* points)
      {
        double* values = new double[numPoints];

        Eigen::Map<Eigen::VectorXd> vecValues(values, numPoints);
        vecValues.setZero();

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
          vecValues.setOnes();
        else if (direction == 1)
          vecValues.setOnes();
        else if (direction == 2)
          vecValues.setZero();
        else
          throw std::runtime_error("Error on direction");

        return values;
      }
      // ***************************************************************************
  };
  // ***************************************************************************
  TEST(TestGeometry, Test_NavierStokes_T1)
  {
    const std::string exportFolder = "./Export/TestGeometry/Test_NavierStokes_T1";
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

      Eigen::VectorXd sol_k = Eigen::VectorXd::Zero(2 * speed_problemData.NumberDOFs +
                                                    pressure_problemData.NumberDOFs);

      double residual_norm = 1.0, solution_norm = 1.0;
      const double newton_tol = 1.0e-6;
      const unsigned int max_iterations = 20;
      int num_iteration = 1;

      Eigen::VectorXd sol_strong = Eigen::VectorXd::Zero(2 * speed_problemData.NumberStrongs +
                                                         pressure_problemData.NumberStrongs);

      if (speed_problemData.NumberStrongs > 0)
      {
        sol_strong.segment(0,
                           speed_problemData.NumberStrongs) =
            GedimForPy::GeDiM4Py_Logic::AssembleStrongSolution(NavierStokes_T1::ExactSpeedSolution_1,
                                                               1,
                                                               meshDAO,
                                                               mesh.Cell2DsMap,
                                                               speed_problemData);
        sol_strong.segment(speed_problemData.NumberStrongs,
                           speed_problemData.NumberStrongs) =
            GedimForPy::GeDiM4Py_Logic::AssembleStrongSolution(NavierStokes_T1::ExactSpeedSolution_2,
                                                               1,
                                                               meshDAO,
                                                               mesh.Cell2DsMap,
                                                               speed_problemData);
      }

      if (pressure_problemData.NumberStrongs > 0)
      {
        sol_strong.segment(2 * speed_problemData.NumberStrongs,
                           pressure_problemData.NumberStrongs) =
            GedimForPy::GeDiM4Py_Logic::AssembleStrongSolution(NavierStokes_T1::ExactPressureSolution,
                                                               1,
                                                               meshDAO,
                                                               mesh.Cell2DsMap,
                                                               pressure_problemData);
      }



      const Eigen::VectorXd u_x_strong = sol_strong.segment(0,
                                                            speed_problemData.NumberStrongs);
      const Eigen::VectorXd u_y_strong = sol_strong.segment(speed_problemData.NumberStrongs,
                                                            speed_problemData.NumberStrongs);
      const Eigen::VectorXd p_strong = sol_strong.segment(2 * speed_problemData.NumberStrongs,
                                                          pressure_problemData.NumberStrongs);

      Eigen::VectorXd speed_x_cell2DsErrorL2 = GedimForPy::GeDiM4Py_Logic::ComputeErrorL2(NavierStokes_T1::ExactPressureSolution,
                                                                                          sol_k.segment(0, speed_problemData.NumberDOFs),
                                                                                          sol_strong.segment(0, speed_problemData.NumberStrongs),
                                                                                          meshDAO,
                                                                                          mesh.Cell2DsMap,
                                                                                          speed_problemData);
      Eigen::VectorXd speed_y_cell2DsErrorL2 = GedimForPy::GeDiM4Py_Logic::ComputeErrorL2(NavierStokes_T1::ExactPressureSolution,
                                                                                          sol_k.segment(speed_problemData.NumberDOFs, speed_problemData.NumberDOFs),
                                                                                          sol_strong.segment(speed_problemData.NumberStrongs, speed_problemData.NumberStrongs),
                                                                                          meshDAO,
                                                                                          mesh.Cell2DsMap,
                                                                                          speed_problemData);
      Eigen::VectorXd pressure_cell2DsErrorL2 = GedimForPy::GeDiM4Py_Logic::ComputeErrorL2(NavierStokes_T1::ExactPressureSolution,
                                                                                           sol_k.segment(2 * speed_problemData.NumberDOFs, pressure_problemData.NumberDOFs),
                                                                                           p_strong,
                                                                                           meshDAO,
                                                                                           mesh.Cell2DsMap,
                                                                                           pressure_problemData);

      while (residual_norm > newton_tol * solution_norm &&
             num_iteration < max_iterations)
      {
        const Eigen::VectorXd u_x_k = sol_k.segment(0,
                                                    speed_problemData.NumberDOFs);
        const Eigen::VectorXd u_y_k = sol_k.segment(speed_problemData.NumberDOFs,
                                                    speed_problemData.NumberDOFs);
        const Eigen::VectorXd p_k = sol_k.segment(2 * speed_problemData.NumberDOFs,
                                                  pressure_problemData.NumberDOFs);

        std::list<Eigen::Triplet<double>> J_saddle_Point_triplets;
        Eigen::VectorXd J_saddlePoint_f = Eigen::VectorXd::Zero(2 * speed_problemData.NumberDOFs +
                                                                pressure_problemData.NumberDOFs);

        {
          std::list<Eigen::Triplet<double>> J_stiffness_dx_Triplets, J_stiffnessStrong_dx_Triplets;
          GedimForPy::GeDiM4Py_Logic::AssembleStiffnessMatrix(NavierStokes::ViscosityTerm,
                                                              meshDAO,
                                                              mesh.Cell2DsMap,
                                                              speed_problemData,
                                                              J_stiffness_dx_Triplets,
                                                              J_stiffnessStrong_dx_Triplets);

          GedimForPy::GeDiM4Py_Logic::ShiftTriplets(J_stiffness_dx_Triplets,
                                                    0,
                                                    0,
                                                    NavierStokes::A_transform_triplet_row,
                                                    NavierStokes::A_transform_triplet_col,
                                                    NavierStokes::A_transform_triplet_value,
                                                    J_saddle_Point_triplets);
        }

        {
          std::list<Eigen::Triplet<double>> J_stiffness_dy_Triplets, J_stiffnessStrong_dy_Triplets;
          GedimForPy::GeDiM4Py_Logic::AssembleStiffnessMatrix(NavierStokes::ViscosityTerm,
                                                              meshDAO,
                                                              mesh.Cell2DsMap,
                                                              speed_problemData,
                                                              J_stiffness_dy_Triplets,
                                                              J_stiffnessStrong_dy_Triplets);

          GedimForPy::GeDiM4Py_Logic::ShiftTriplets(J_stiffness_dy_Triplets,
                                                    speed_problemData.NumberDOFs,
                                                    speed_problemData.NumberDOFs,
                                                    NavierStokes::A_transform_triplet_row,
                                                    NavierStokes::A_transform_triplet_col,
                                                    NavierStokes::A_transform_triplet_value,
                                                    J_saddle_Point_triplets);
        }

        {
          std::list<Eigen::Triplet<double>> J_advection_dx_Triplets, J_advectionStrong_dx_Triplets;
          GedimForPy::GeDiM4Py_Logic::AssembleAdvectionMatrix(NavierStokes::AdvectionTerm_1,
                                                              meshDAO,
                                                              mesh.Cell2DsMap,
                                                              speed_problemData,
                                                              pressure_problemData,
                                                              J_advection_dx_Triplets,
                                                              J_advectionStrong_dx_Triplets);

          GedimForPy::GeDiM4Py_Logic::ShiftTriplets(J_advection_dx_Triplets,
                                                    0,
                                                    2.0 * speed_problemData.NumberDOFs,
                                                    NavierStokes::BT_transform_triplet_row,
                                                    NavierStokes::BT_transform_triplet_col,
                                                    NavierStokes::BT_transform_triplet_value,
                                                    J_saddle_Point_triplets);
          GedimForPy::GeDiM4Py_Logic::ShiftTriplets(J_advection_dx_Triplets,
                                                    2 * speed_problemData.NumberDOFs,
                                                    0,
                                                    NavierStokes::B_transform_triplet_row,
                                                    NavierStokes::B_transform_triplet_col,
                                                    NavierStokes::B_transform_triplet_value,
                                                    J_saddle_Point_triplets);

        }

        {
          std::list<Eigen::Triplet<double>> J_advection_dy_Triplets, J_advectionStrong_dy_Triplets;
          GedimForPy::GeDiM4Py_Logic::AssembleAdvectionMatrix(NavierStokes::AdvectionTerm_2,
                                                              meshDAO,
                                                              mesh.Cell2DsMap,
                                                              speed_problemData,
                                                              pressure_problemData,
                                                              J_advection_dy_Triplets,
                                                              J_advectionStrong_dy_Triplets);
          GedimForPy::GeDiM4Py_Logic::ShiftTriplets(J_advection_dy_Triplets,
                                                    speed_problemData.NumberDOFs,
                                                    2.0 * speed_problemData.NumberDOFs,
                                                    NavierStokes::BT_transform_triplet_row,
                                                    NavierStokes::BT_transform_triplet_col,
                                                    NavierStokes::BT_transform_triplet_value,
                                                    J_saddle_Point_triplets);
          GedimForPy::GeDiM4Py_Logic::ShiftTriplets(J_advection_dy_Triplets,
                                                    2 * speed_problemData.NumberDOFs,
                                                    speed_problemData.NumberDOFs,
                                                    NavierStokes::B_transform_triplet_row,
                                                    NavierStokes::B_transform_triplet_col,
                                                    NavierStokes::B_transform_triplet_value,
                                                    J_saddle_Point_triplets);

        }

        {
          const Eigen::VectorXd J_forcingTerm_f_1 = GedimForPy::GeDiM4Py_Logic::AssembleForcingTerm(NavierStokes_T1::ForcingTerm_1,
                                                                                                    meshDAO,
                                                                                                    mesh.Cell2DsMap,
                                                                                                    speed_problemData);
          const Eigen::VectorXd J_forcingTerm_f_2 = GedimForPy::GeDiM4Py_Logic::AssembleForcingTerm(NavierStokes_T1::ForcingTerm_2,
                                                                                                    meshDAO,
                                                                                                    mesh.Cell2DsMap,
                                                                                                    speed_problemData);

          J_saddlePoint_f.segment(0, speed_problemData.NumberDOFs) += J_forcingTerm_f_1;
          J_saddlePoint_f.segment(speed_problemData.NumberDOFs, speed_problemData.NumberDOFs) += J_forcingTerm_f_2;
        }


        {
          const Eigen::VectorXd J_forcingTerm_double_dot_u_x = GedimForPy::GeDiM4Py_Logic::AssembleNonLinearDerivativeForcingTerm(NavierStokes::OnesDerivative,
                                                                                                                                  NavierStokes::NonLinear_double_dot_product,
                                                                                                                                  meshDAO,
                                                                                                                                  mesh.Cell2DsMap,
                                                                                                                                  speed_problemData,
                                                                                                                                  u_x_k,
                                                                                                                                  u_x_strong);
          const Eigen::VectorXd J_forcingTerm_double_dot_u_y = GedimForPy::GeDiM4Py_Logic::AssembleNonLinearDerivativeForcingTerm(NavierStokes::OnesDerivative,
                                                                                                                                  NavierStokes::NonLinear_double_dot_product,
                                                                                                                                  meshDAO,
                                                                                                                                  mesh.Cell2DsMap,
                                                                                                                                  speed_problemData,
                                                                                                                                  u_y_k,
                                                                                                                                  u_y_strong);
          J_saddlePoint_f.segment(0, speed_problemData.NumberDOFs) += J_forcingTerm_double_dot_u_x;
          J_saddlePoint_f.segment(speed_problemData.NumberDOFs, speed_problemData.NumberDOFs) += J_forcingTerm_double_dot_u_y;
        }

        Eigen::SparseMatrix<double> J_saddle_point(2 * speed_problemData.NumberDOFs +
                                                   pressure_problemData.NumberDOFs,
                                                   2 * speed_problemData.NumberDOFs +
                                                   pressure_problemData.NumberDOFs);

        J_saddle_point.setFromTriplets(J_saddle_Point_triplets.begin(),
                                       J_saddle_Point_triplets.end());
        J_saddle_point.makeCompressed();
        J_saddle_Point_triplets.clear();



        Eigen::SparseLU<Eigen::SparseMatrix<double>> linearSolver;
        linearSolver.compute(J_saddle_point);

        const Eigen::VectorXd d_sol = linearSolver.solve(J_saddlePoint_f);
        sol_k = sol_k + d_sol;

        speed_x_cell2DsErrorL2 = GedimForPy::GeDiM4Py_Logic::ComputeErrorL2(NavierStokes_T1::ExactPressureSolution,
                                                                            u_x_k,
                                                                            u_x_strong,
                                                                            meshDAO,
                                                                            mesh.Cell2DsMap,
                                                                            speed_problemData);
        speed_y_cell2DsErrorL2 = GedimForPy::GeDiM4Py_Logic::ComputeErrorL2(NavierStokes_T1::ExactPressureSolution,
                                                                            u_y_k,
                                                                            u_y_strong,
                                                                            meshDAO,
                                                                            mesh.Cell2DsMap,
                                                                            speed_problemData);
        pressure_cell2DsErrorL2 = GedimForPy::GeDiM4Py_Logic::ComputeErrorL2(NavierStokes_T1::ExactPressureSolution,
                                                                             p_k,
                                                                             p_strong,
                                                                             meshDAO,
                                                                             mesh.Cell2DsMap,
                                                                             pressure_problemData);

        const Eigen::VectorXd speed_x_cell2DsNormL2 = GedimForPy::GeDiM4Py_Logic::ComputeErrorL2(NavierStokes::ZeroSolution,
                                                                                                 u_x_k,
                                                                                                 u_x_strong,
                                                                                                 meshDAO,
                                                                                                 mesh.Cell2DsMap,
                                                                                                 speed_problemData);
        const Eigen::VectorXd speed_y_cell2DsNormL2 = GedimForPy::GeDiM4Py_Logic::ComputeErrorL2(NavierStokes::ZeroSolution,
                                                                                                 u_y_k,
                                                                                                 u_y_strong,
                                                                                                 meshDAO,
                                                                                                 mesh.Cell2DsMap,
                                                                                                 speed_problemData);
        const Eigen::VectorXd pressure_cell2DsNormL2 = GedimForPy::GeDiM4Py_Logic::ComputeErrorL2(NavierStokes::ZeroSolution,
                                                                                                  p_k,
                                                                                                  p_strong,
                                                                                                  meshDAO,
                                                                                                  mesh.Cell2DsMap,
                                                                                                  pressure_problemData);

        const Eigen::VectorXd speed_x_cell2Ds_d_sol_NormL2 = GedimForPy::GeDiM4Py_Logic::ComputeErrorL2(NavierStokes::ZeroSolution,
                                                                                                        d_sol.segment(0, speed_problemData.NumberDOFs),
                                                                                                        u_x_strong,
                                                                                                        meshDAO,
                                                                                                        mesh.Cell2DsMap,
                                                                                                        speed_problemData);
        const Eigen::VectorXd speed_y_cell2Ds_d_sol_NormL2 = GedimForPy::GeDiM4Py_Logic::ComputeErrorL2(NavierStokes::ZeroSolution,
                                                                                                        d_sol.segment(speed_problemData.NumberDOFs, speed_problemData.NumberDOFs),
                                                                                                        u_y_strong,
                                                                                                        meshDAO,
                                                                                                        mesh.Cell2DsMap,
                                                                                                        speed_problemData);
        const Eigen::VectorXd pressure_cell2Ds_d_sol_NormL2 = GedimForPy::GeDiM4Py_Logic::ComputeErrorL2(NavierStokes::ZeroSolution,
                                                                                                         d_sol.segment(2 * speed_problemData.NumberDOFs, pressure_problemData.NumberDOFs),
                                                                                                         p_strong,
                                                                                                         meshDAO,
                                                                                                         mesh.Cell2DsMap,
                                                                                                         pressure_problemData);


        solution_norm = std::sqrt(speed_x_cell2DsNormL2.sum() + speed_y_cell2DsNormL2.sum() + pressure_cell2DsNormL2.sum());
        residual_norm = std::sqrt(speed_x_cell2Ds_d_sol_NormL2.sum() + speed_y_cell2Ds_d_sol_NormL2.sum() + pressure_cell2Ds_d_sol_NormL2.sum());


#if ACTIVE_CHECK == 0
        std::cerr.precision(3);
        std::cerr<< std::scientific
                 << "dofs"<< ","
                 << "h"<< ","
                 << "speed_x_errorL2"<< ","
                 << "speed_y_errorL2"<< ","
                 << "pressure_errorL2"<< ","
                 << "speed_x_normL2"<< ","
                 << "speed_y_normL2"<< ","
                 << "pressure_normL2"
                 << std::endl;
        std::cerr<< std::scientific
                 << 2.0 * speed_problemData.NumberDOFs + pressure_problemData.NumberDOFs<< ","
                 << speed_problemData.H<< ","
                 << sqrt(speed_x_cell2DsErrorL2.sum())<< ","
                 << sqrt(speed_y_cell2DsErrorL2.sum())<< ","
                 << sqrt(pressure_cell2DsErrorL2.sum())<< ","
                 << sqrt(speed_x_cell2DsNormL2.sum())<< ","
                 << sqrt(speed_y_cell2DsNormL2.sum())<< ","
                 << sqrt(pressure_cell2DsNormL2.sum())
                 << std::endl;

        std::cout.precision(3);
        std::cout<< std::scientific<<
                    " Newton it "<< num_iteration<< " / "<< max_iterations<<
                    " residual "<< residual_norm<< " / "<< newton_tol * solution_norm<< std::endl;
#endif

        num_iteration++;
      }


      // export
      {
        {
          std::vector<double> pressure_cell0Ds_numeric_solution(meshDAO.Cell0DTotalNumber(),
                                                                0.0);
          std::vector<double> speed_x_cell0Ds_numeric_solution(meshDAO.Cell0DTotalNumber(),
                                                               0.0);
          std::vector<double> speed_y_cell0Ds_numeric_solution(meshDAO.Cell0DTotalNumber(),
                                                               0.0);

          const Eigen::MatrixXd coordinates = meshDAO.Cell0DsCoordinates();

          const double* pressure_cell0Ds_exact_solution = NavierStokes_T1::ExactPressureSolution(coordinates.cols(),
                                                                                                 coordinates.data());
          const double* speed_x_cell0Ds_exact_solution = NavierStokes_T1::ExactSpeedSolution_1(coordinates.cols(),
                                                                                               coordinates.data());
          const double* speed_y_cell0Ds_exact_solution = NavierStokes_T1::ExactSpeedSolution_2(coordinates.cols(),
                                                                                               coordinates.data());

          for (unsigned int p = 0; p < meshDAO.Cell0DTotalNumber(); p++)
          {
            const GedimForPy::DiscreteProblemData::DOF& speed_dof = speed_problemData.Cell0Ds_DOF[p];
            const GedimForPy::DiscreteProblemData::DOF& pressure_dof = pressure_problemData.Cell0Ds_DOF[p];

            switch (speed_dof.Type)
            {
              case GedimForPy::DiscreteProblemData::DOF::Types::DOF:
                speed_x_cell0Ds_numeric_solution[p] = sol_k[speed_dof.Global_Index];
                speed_y_cell0Ds_numeric_solution[p] = sol_k[speed_problemData.NumberDOFs + speed_dof.Global_Index];
                break;
              case GedimForPy::DiscreteProblemData::DOF::Types::Strong:
                speed_x_cell0Ds_numeric_solution[p] = sol_strong[speed_dof.Global_Index];
                speed_y_cell0Ds_numeric_solution[p] = sol_strong[speed_problemData.NumberStrongs + speed_dof.Global_Index];
                break;
              default:
                throw std::runtime_error("DOF Type " +
                                         std::to_string((unsigned int)speed_dof.Type) +
                                         " not supported");
            }

            switch (pressure_dof.Type)
            {
              case GedimForPy::DiscreteProblemData::DOF::Types::DOF:
                pressure_cell0Ds_numeric_solution[p] = sol_k[2 * speed_problemData.NumberDOFs + pressure_dof.Global_Index];
                break;
              case GedimForPy::DiscreteProblemData::DOF::Types::Strong:
                pressure_cell0Ds_numeric_solution[p] = sol_strong[2 * speed_problemData.NumberStrongs + pressure_dof.Global_Index];
                break;
              default:
                throw std::runtime_error("DOF Type " +
                                         std::to_string((unsigned int)pressure_dof.Type) +
                                         " not supported");
            }
          }

          GedimForPy::GeDiM4Py_Logic::ExportSolution(NavierStokes_T1::ExactPressureSolution,
                                                     sol_k.segment(2 * speed_problemData.NumberDOFs,
                                                                   pressure_problemData.NumberDOFs),
                                                     p_strong,
                                                     meshDAO,
                                                     pressure_problemData,
                                                     {
                                                       exportFolder,
                                                       "Pressure"
                                                     });
          GedimForPy::GeDiM4Py_Logic::ExportSolution(NavierStokes_T1::ExactSpeedSolution_1,
                                                     sol_k.segment(0,
                                                                   speed_problemData.NumberDOFs),
                                                     u_x_strong,
                                                     meshDAO,
                                                     speed_problemData,
                                                     {
                                                       exportFolder,
                                                       "Speed_1"
                                                     });
          GedimForPy::GeDiM4Py_Logic::ExportSolution(NavierStokes_T1::ExactSpeedSolution_2,
                                                     sol_k.segment(speed_problemData.NumberDOFs,
                                                                   speed_problemData.NumberDOFs),
                                                     u_y_strong,
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
                                   "speed_x_cell2Ds_errorL2",
                                   Gedim::VTPProperty::Formats::Cells,
                                   static_cast<unsigned int>(speed_x_cell2DsErrorL2.size()),
                                   speed_x_cell2DsErrorL2.data()
                                 },
                                 {
                                   "speed_y_cell2Ds_errorL2",
                                   Gedim::VTPProperty::Formats::Cells,
                                   static_cast<unsigned int>(speed_y_cell2DsErrorL2.size()),
                                   speed_y_cell2DsErrorL2.data()
                                 },
                                 {
                                   "pressure_cell2Ds_errorL2",
                                   Gedim::VTPProperty::Formats::Cells,
                                   static_cast<unsigned int>(pressure_cell2DsErrorL2.size()),
                                   pressure_cell2DsErrorL2.data()
                                 },
                                 {
                                   "speed_x_cell0Ds_exact_solution",
                                   Gedim::VTPProperty::Formats::Points,
                                   static_cast<unsigned int>(coordinates.cols()),
                                   speed_x_cell0Ds_exact_solution
                                 },
                                 {
                                   "speed_x_cell0Ds_numeric_solution",
                                   Gedim::VTPProperty::Formats::Points,
                                   static_cast<unsigned int>(speed_x_cell0Ds_numeric_solution.size()),
                                   speed_x_cell0Ds_numeric_solution.data()
                                 },
                                 {
                                   "speed_y_cell0Ds_exact_solution",
                                   Gedim::VTPProperty::Formats::Points,
                                   static_cast<unsigned int>(coordinates.cols()),
                                   speed_y_cell0Ds_exact_solution
                                 },
                                 {
                                   "speed_y_cell0Ds_numeric_solution",
                                   Gedim::VTPProperty::Formats::Points,
                                   static_cast<unsigned int>(speed_y_cell0Ds_numeric_solution.size()),
                                   speed_y_cell0Ds_numeric_solution.data()
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
