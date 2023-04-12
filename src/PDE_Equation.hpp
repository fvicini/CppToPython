#ifndef __PDE_Equation_H
#define __PDE_Equation_H

#include "Eigen/Eigen"

namespace GedimForPy
{
  /// \brief 2D Primal Conforming Constant Lagrange Element Degree variable
  class PDE_Equation final
  {
    public:
      static Eigen::MatrixXd ComputeStiffnessMatrix(const unsigned int& numCellDofs,
                                                    const std::vector<Eigen::MatrixXd>& basisFunctionDerivativeValues,
                                                    const Eigen::MatrixXd& diffusionTermValues,
                                                    const Eigen::VectorXd& quadratureWeights);
      static Eigen::MatrixXd ComputeCellAdvectionMatrix(const unsigned int& num_Cell_trialDofs,
                                                        const unsigned int& num_Cell_testDofs,
                                                        const Eigen::MatrixXd& advectionTermValues,
                                                        const Eigen::MatrixXd& test_basisFunctionValues,
                                                        const std::vector<Eigen::MatrixXd>& trial_basisFunctionDerivativeValues,
                                                        const Eigen::VectorXd& quadratureWeights);
      static Eigen::MatrixXd ComputeCellReactionMatrix(const Eigen::VectorXd& reactionTermValues,
                                                       const Eigen::MatrixXd& basisFunctionValues,
                                                       const Eigen::VectorXd& quadratureWeights);
      static Eigen::VectorXd ComputeCellForcingTerm(const Eigen::VectorXd& forcingTermValues,
                                                    const Eigen::MatrixXd& basisFunctionValues,
                                                    const Eigen::VectorXd& quadratureWeights);
      static Eigen::VectorXd ComputeCellWeakTerm(const Eigen::VectorXd& weakTermValues,
                                                 const Eigen::MatrixXd& basisFunctionValues,
                                                 const Eigen::VectorXd& quadratureWeights);

      static Eigen::MatrixXd ComputeStokes_StiffnessMatrix(const unsigned int& speed_singleComponent_numDofs,
                                                           const std::vector<Eigen::MatrixXd>& speed_basisFunctionDerivativeValues,
                                                           const Eigen::VectorXd& viscosityTermValues,
                                                           const Eigen::VectorXd& quadratureWeights);
      static Eigen::MatrixXd ComputeStokes_AdvectionMatrix(const unsigned int& speed_singleComponent_numDofs,
                                                           const unsigned int& pressure_numDofs,
                                                           const Eigen::MatrixXd& pressure_basisFunctionValues,
                                                           const std::vector<Eigen::MatrixXd>& speed_basisFunctionDerivativeValues,
                                                           const Eigen::VectorXd& quadratureWeights);
      static Eigen::VectorXd ComputeStokes_ForcingTerm(const unsigned int& speed_singleComponent_numDofs,
                                                       const std::vector<Eigen::VectorXd>& forcingTermValues,
                                                       const Eigen::MatrixXd& speed_basisFunctionValues,
                                                       const Eigen::VectorXd& quadratureWeights);
  };
}

#endif // __PDE_Equation_H
