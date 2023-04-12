#include "PDE_Equation.hpp"

using namespace std;
using namespace Eigen;

namespace GedimForPy
{
  // ***************************************************************************
  Eigen::MatrixXd PDE_Equation::ComputeStiffnessMatrix(const unsigned int& numCellDofs,
                                                       const std::vector<Eigen::MatrixXd>& basisFunctionDerivativeValues,
                                                       const Eigen::MatrixXd& diffusionTermValues,
                                                       const Eigen::VectorXd& quadratureWeights)
  {
    MatrixXd cellMatrixA;
    cellMatrixA.setZero(numCellDofs,
                        numCellDofs);

    for(unsigned int i = 0; i < basisFunctionDerivativeValues.size(); i++)
    {
      for(unsigned int j = 0; j < basisFunctionDerivativeValues.size(); j++)
      {
        cellMatrixA += basisFunctionDerivativeValues[i].transpose() *
                       quadratureWeights.cwiseProduct(diffusionTermValues.row(i + j).transpose()).asDiagonal() *
                       basisFunctionDerivativeValues[j];
      }
    }

    return cellMatrixA;
  }
  // ***************************************************************************
  MatrixXd PDE_Equation::ComputeCellAdvectionMatrix(const unsigned int& num_Cell_trialDofs,
                                                    const unsigned int& num_Cell_testDofs,
                                                    const Eigen::MatrixXd& advectionTermValues,
                                                    const Eigen::MatrixXd& test_basisFunctionValues,
                                                    const std::vector<Eigen::MatrixXd>& trial_basisFunctionDerivativeValues,
                                                    const Eigen::VectorXd& quadratureWeights)
  {
    Eigen::MatrixXd cellMatrixA;
    cellMatrixA.setZero(num_Cell_testDofs,
                        num_Cell_trialDofs);

    for(unsigned int i = 0; i < trial_basisFunctionDerivativeValues.size(); i++)
      cellMatrixA += test_basisFunctionValues.transpose() *
                     quadratureWeights.cwiseProduct(advectionTermValues.row(i).transpose()).asDiagonal() *
                     trial_basisFunctionDerivativeValues[i];

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
  MatrixXd PDE_Equation::ComputeStokes_StiffnessMatrix(const unsigned int& speed_singleComponent_numDofs,
                                                       const std::vector<Eigen::MatrixXd>& speed_basisFunctionDerivativeValues,
                                                       const Eigen::VectorXd& viscosityTermValues,
                                                       const Eigen::VectorXd& quadratureWeights)
  {
    MatrixXd cellMatrixA = MatrixXd::Zero(2 * speed_singleComponent_numDofs,
                                          2 * speed_singleComponent_numDofs);

    for (unsigned int i = 0; i < speed_singleComponent_numDofs; i++)
    {
      for (unsigned int j = 0; j < speed_singleComponent_numDofs; j++)
      {
        cellMatrixA(i, j) +=
            speed_basisFunctionDerivativeValues[0].col(i).transpose() *
            quadratureWeights.cwiseProduct(viscosityTermValues).asDiagonal() *
            speed_basisFunctionDerivativeValues[0].col(j);

        cellMatrixA(speed_singleComponent_numDofs + i,
                    speed_singleComponent_numDofs + j) +=
            speed_basisFunctionDerivativeValues[1].col(i).transpose() *
            quadratureWeights.cwiseProduct(viscosityTermValues).asDiagonal() *
            speed_basisFunctionDerivativeValues[1].col(j);
      }
    }

    return cellMatrixA;
  }
  // ***************************************************************************
  MatrixXd PDE_Equation::ComputeStokes_AdvectionMatrix(const unsigned int& speed_singleComponent_numDofs,
                                                       const unsigned int& pressure_numDofs,
                                                       const Eigen::MatrixXd& pressure_basisFunctionValues,
                                                       const std::vector<Eigen::MatrixXd>& speed_basisFunctionDerivativeValues,
                                                       const Eigen::VectorXd& quadratureWeights)
  {
    Eigen::MatrixXd cellMatrixB = MatrixXd::Zero(pressure_numDofs,
                                                 2 * speed_singleComponent_numDofs);

    for (unsigned int i = 0; i < pressure_numDofs; i++)
    {
      for (unsigned int j = 0; j < speed_singleComponent_numDofs; j++)
      {
        cellMatrixB(i, j) +=
            pressure_basisFunctionValues.col(i).transpose() *
            quadratureWeights.asDiagonal() *
            speed_basisFunctionDerivativeValues[0].col(j);

        cellMatrixB(i,
                    speed_singleComponent_numDofs + j) +=
            pressure_basisFunctionValues.col(i).transpose() *
            quadratureWeights.asDiagonal() *
            speed_basisFunctionDerivativeValues[1].col(j);
      }
    }

    return cellMatrixB;
  }
  // ***************************************************************************
  VectorXd PDE_Equation::ComputeStokes_ForcingTerm(const unsigned int& speed_singleComponent_numDofs,
                                                   const std::vector<Eigen::VectorXd>& forcingTermValues,
                                                   const Eigen::MatrixXd& speed_basisFunctionValues,
                                                   const Eigen::VectorXd& quadratureWeights)
  {
    Eigen::VectorXd cellForcingTerm = VectorXd::Zero(2 * speed_singleComponent_numDofs);

    for (unsigned int i = 0; i < speed_singleComponent_numDofs; i++)
    {
      cellForcingTerm(i) +=
          speed_basisFunctionValues.col(i).transpose() *
          quadratureWeights.asDiagonal() *
          forcingTermValues[0];

      cellForcingTerm(speed_singleComponent_numDofs + i) +=
          speed_basisFunctionValues.col(i).transpose() *
          quadratureWeights.asDiagonal() *
          forcingTermValues[1];
    }

    return cellForcingTerm;
  }
  // ***************************************************************************
}
