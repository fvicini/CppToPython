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
}
