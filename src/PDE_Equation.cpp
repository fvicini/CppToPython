#include "PDE_Equation.hpp"

using namespace std;
using namespace Eigen;

namespace GedimForPy
{
  // ***************************************************************************
  Eigen::MatrixXd PDE_Equation::ComputeStiffnessMatrix(const unsigned int& numCellDofs,
                                                       const std::vector<Eigen::MatrixXd>& basisFunctionDerivativeValues,
                                                       const VectorXd& diffusionTermValues,
                                                       const Eigen::VectorXd& quadratureWeights)
  {
    MatrixXd cellMatrixA;
    cellMatrixA.setZero(numCellDofs,
                        numCellDofs);

    for(unsigned int i = 0; i < basisFunctionDerivativeValues.size(); i++)
      cellMatrixA += basisFunctionDerivativeValues[i].transpose() *
                     quadratureWeights.cwiseProduct(diffusionTermValues).asDiagonal() *
                     basisFunctionDerivativeValues[i];

    return cellMatrixA;
  }
  // ***************************************************************************
  MatrixXd PDE_Equation::ComputeCellAdvectionMatrix(const unsigned int& numCellDofs,
                                                    const Eigen::MatrixXd& advectionTermValues,
                                                    const Eigen::MatrixXd& basisFunctionValues,
                                                    const std::vector<Eigen::MatrixXd>& basisFunctionDerivativeValues,
                                                    const Eigen::VectorXd& quadratureWeights)
  {
    Eigen::MatrixXd cellMatrixA;
    cellMatrixA.setZero(numCellDofs,
                        numCellDofs);

    for(unsigned int i = 0; i < basisFunctionDerivativeValues.size(); i++)
      cellMatrixA += basisFunctionValues.transpose() *
                     quadratureWeights.cwiseProduct(advectionTermValues.row(i).transpose()).asDiagonal() *
                     basisFunctionDerivativeValues[i];

    return cellMatrixA;
  }
  // ***************************************************************************
  MatrixXd PDE_Equation::ComputeCellReactionMatrix(const Eigen::VectorXd& reactionTermValues,
                                                   const Eigen::MatrixXd& basisFunctionValues,
                                                   const Eigen::VectorXd& quadratureWeights)
  {
    return basisFunctionValues.transpose() *
        quadratureWeights.cwiseProduct(reactionTermValues).asDiagonal() *
        basisFunctionValues;
  }
  // ***************************************************************************
  Eigen::VectorXd PDE_Equation::ComputeCellForcingTerm(const Eigen::VectorXd& forcingTermValues,
                                                       const Eigen::MatrixXd& basisFunctionValues,
                                                       const Eigen::VectorXd& quadratureWeights)
  {
    return
        basisFunctionValues.transpose() *
        quadratureWeights.asDiagonal() *
        forcingTermValues;
  }
  // ***************************************************************************
  VectorXd PDE_Equation::ComputeCellWeakTerm(const Eigen::VectorXd& weakTermValues,
                                             const Eigen::MatrixXd& basisFunctionValues,
                                             const Eigen::VectorXd& quadratureWeights)
  {
    return
        basisFunctionValues.transpose() *
        quadratureWeights.asDiagonal() *
        weakTermValues;
  }
  // ***************************************************************************
}
