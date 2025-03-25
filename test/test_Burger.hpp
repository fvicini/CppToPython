#ifndef __TEST_BURGER_H
#define __TEST_BURGER_H

#include "test_utilities.hpp"

namespace UnitTesting
{
  class Burger final
  {
    public:
      // ***************************************************************************
      static double* DiffusionTerm(const int numPoints,
                                   const double* points)
      {
        double* values = new double[numPoints];

        Eigen::Map<Eigen::VectorXd> matValues(values, numPoints);
        matValues.setConstant(1.0);

        return values;
      }
      // ***************************************************************************
      static double* AdvectionTerm(const int numPoints,
                                   const double* points)
      {
        double* values = new double[2 * numPoints];

        Eigen::Map<Eigen::MatrixXd> matValues(values, 2, numPoints);
        matValues.row(0).setConstant(1.0);
        matValues.row(1).setConstant(1.0);

        return values;
      }
      // ***************************************************************************
      static double* ReactionTerm(const int numPoints,
                                  const double* points)
      {
        double* values = new double[numPoints];

        Eigen::Map<Eigen::VectorXd> vecValues(values, numPoints);
        vecValues.setConstant(1.0);

        return values;
      }
      // ***************************************************************************
      static double* NonLinear_f_der_v(const int numPoints,
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
      static double* NonLinear_f_v(const int numPoints,
                                   const double* points,
                                   const double* u,
                                   const double* u_x,
                                   const double* u_y)
      {
        double* values = new double[numPoints];

        Eigen::Map<Eigen::VectorXd> matValues(values,
                                              numPoints);
        matValues = Eigen::Map<const Eigen::VectorXd>(u,
                                                      numPoints).array() *
                    (Eigen::Map<const Eigen::VectorXd>(u_x,
                                                       numPoints) +
                     Eigen::Map<const Eigen::VectorXd>(u_y,
                                                       numPoints)).array();

        return values;
      }
      // ***************************************************************************
      static double* NonLinear_Reaction(const int numPoints,
                                        const double* points,
                                        const double* u,
                                        const double* u_x,
                                        const double* u_y)
      {
        double* values = new double[numPoints];

        Eigen::Map<Eigen::VectorXd> matValues(values,
                                              numPoints);
        matValues = Eigen::Map<const Eigen::VectorXd>(u_x,
                                                      numPoints) +
                    Eigen::Map<const Eigen::VectorXd>(u_y,
                                                      numPoints);

        return values;
      }
      // ***************************************************************************
      static double* NonLinear_Advection(const int numPoints,
                                         const double* points,
                                         const double* u,
                                         const double* u_x,
                                         const double* u_y)
      {
        double* values = new double[numPoints];

        Eigen::Map<Eigen::VectorXd> matValues(values,
                                              numPoints);
        matValues = Eigen::Map<const Eigen::VectorXd>(u,
                                                      numPoints);

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
      static double* ForcingTerm(const int numPoints, const double* points)
      {
        double* values = new double[numPoints];

        Eigen::Map<const Eigen::MatrixXd> matPoints(points, 3, numPoints);
        Eigen::Map<Eigen::VectorXd> vecValues(values, numPoints);
        vecValues = 32.0 * (matPoints.row(1).array() * (1.0 - matPoints.row(1).array()) +
                            matPoints.row(0).array() * (1.0 - matPoints.row(0).array())).array() +
                    16.0 * ((1.0 - 2.0 * matPoints.row(0).array()) *
                            matPoints.row(1).array() *
                            (1.0 - matPoints.row(1).array()) +
                            (1.0 - 2.0 * matPoints.row(1).array()) *
                            matPoints.row(0).array() *
                            (1.0 - matPoints.row(0).array())).array() *
                    16.0 * (matPoints.row(1).array() *
                            (1.0 - matPoints.row(1).array()) *
                            matPoints.row(0).array() *
                            (1.0 - matPoints.row(0).array()));
        return values;
      }
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
      static double* ExactSolution(const int numPoints, const double* points)
      {
        double* values = new double[numPoints];

        Eigen::Map<const Eigen::MatrixXd> matPoints(points, 3, numPoints);
        Eigen::Map<Eigen::VectorXd> vecValues(values, numPoints);

        vecValues = 16.0 * (matPoints.row(1).array() * (1.0 - matPoints.row(1).array()) *
                            matPoints.row(0).array() * (1.0 - matPoints.row(0).array()));
        return values;
      }
      // ***************************************************************************
      static double* ExactDerivativeSolution(const int direction,
                                             const int numPoints,
                                             const double* points)
      {
        double* values = new double[numPoints];

        Eigen::Map<const Eigen::MatrixXd> matPoints(points, 3, numPoints);
        Eigen::Map<Eigen::VectorXd> vecValues(values, numPoints);

        if(direction == 0)
          vecValues = 16.0 * (1.0 - 2.0 * matPoints.row(0).array()) *
                      matPoints.row(1).array() *
                      (1.0 - matPoints.row(1).array());
        else if (direction == 1)
          vecValues = 16.0 * (1.0 - 2.0 * matPoints.row(1).array()) *
                      matPoints.row(0).array() *
                      (1.0 - matPoints.row(0).array());
        else if (direction == 2)
          vecValues.setZero();
        else
          throw std::runtime_error("Error on direction");

        return values;
      }
      // ***************************************************************************
  };
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
      const unsigned int max_iterations = 5;
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
                 << "errorH1"<< ","<< "residual"<< std::endl;
        std::cerr<< std::scientific<< problemData.NumberDOFs<< ","
                 << problemData.H<< ","<< sqrt(cell2DsErrorL2.sum()) / sqrt(cell2DsNormL2.sum())<< ","
                 << sqrt(cell2DsErrorH1.sum()) / sqrt(cell2DsNormH1.sum())<< ","
                 << residual_norm / solution_norm
                 << std::endl;

        std::cout.precision(3);
        std::cout<< std::scientific<<
                    " Newton it "<< num_iteration<< " / "<< max_iterations<<
                    " residual "<< residual_norm<< " / "<< newton_tol * solution_norm<< std::endl;
#endif

        num_iteration++;
      }

#if ACTIVE_CHECK == 1
      ASSERT_EQ(5, num_iteration);
#endif
    }

    gedimData.Destroy();
  }
}

#endif // __TEST_BURGER_H
