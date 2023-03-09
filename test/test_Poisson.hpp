#ifndef __TEST_POISSON_H
#define __TEST_POISSON_H

#include "Eigen/Eigen"

namespace UnitTesting
{
  class Poisson final
  {
    public:
      // ***************************************************************************
      static Eigen::VectorXd DiffusionTerm(const Eigen::MatrixXd& points)
      {
        return Eigen::VectorXd::Constant(points.cols(), 1.0);
      }
      // ***************************************************************************
      static Eigen::VectorXd ForcingTerm(const Eigen::MatrixXd& points)
      {
        return 32.0 * (points.row(1).array() * (1.0 - points.row(1).array()) +
                       points.row(0).array() * (1.0 - points.row(0).array()));
      }
      // ***************************************************************************
      static Eigen::VectorXd ExactSolution(const Eigen::MatrixXd& points)
      {
        return 16.0 * (points.row(1).array() * (1.0 - points.row(1).array()) *
                       points.row(0).array() * (1.0 - points.row(0).array()));
      }
      // ***************************************************************************
      static Eigen::VectorXd ExactDerivativeSolution(const unsigned int& direction,
                                                     const Eigen::MatrixXd& points)
      {
        if(direction == 0)
          return 16.0 * (1.0 - 2.0 * points.row(0).array()) * points.row(1).array() * (1.0 - points.row(1).array());
        else if (direction == 1)
          return 16.0 * (1.0 - 2.0 * points.row(1).array()) * points.row(0).array() * (1.0 - points.row(0).array());
        else if (direction == 2)
          return Eigen::VectorXd::Zero(points.cols());
        else
          throw std::runtime_error("Error on direction");
      }
      // ***************************************************************************
  };
}

#endif // __TEST_POISSON_H
