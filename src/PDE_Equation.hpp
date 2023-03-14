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
                                                    const Eigen::VectorXd& diffusionTermValues,
                                                    const Eigen::VectorXd& quadratureWeights);
      static Eigen::MatrixXd ComputeCellAdvectionMatrix(const unsigned int& numCellDofs,
                                                        const Eigen::MatrixXd& advectionTermValues,
                                                        const Eigen::MatrixXd& basisFunctionValues,
                                                        const std::vector<Eigen::MatrixXd>& basisFunctionDerivativeValues,
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
  };
}

#endif // __PDE_Equation_H
