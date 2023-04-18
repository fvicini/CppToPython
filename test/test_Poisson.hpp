#ifndef __TEST_POISSON_H
#define __TEST_POISSON_H

#include "Eigen/Eigen"

namespace UnitTesting
{
  class Poisson final
  {
    public:
      static constexpr double A() { return 10.0; }
      static constexpr double B() { return 0.1; }
      static constexpr double C() { return 2.0; }

      // ***************************************************************************
      static double* DiffusionTerm(const int numPoints, const double* points)
      {
        double* values = new double[numPoints];

        Eigen::Map<Eigen::VectorXd> matValues(values, numPoints);
        matValues.setConstant(A());

        return values;
      }
      // ***************************************************************************
      static double* AdvectionTerm(const int numPoints, const double* points)
      {
        double* values = new double[2 * numPoints];

        Eigen::Map<Eigen::MatrixXd> matValues(values, 2, numPoints);
        matValues.row(0).setConstant(B());
        matValues.row(1).setConstant(B());

        return values;
      }
      // ***************************************************************************
      static double* ReactionTerm(const int numPoints, const double* points)
      {
        double* values = new double[numPoints];

        Eigen::Map<Eigen::VectorXd> vecValues(values, numPoints);
        vecValues.setConstant(C());

        return values;
      }
      // ***************************************************************************
      static double* ForcingTerm(const int numPoints, const double* points)
      {
        double* values = new double[numPoints];

        Eigen::Map<const Eigen::MatrixXd> matPoints(points, 3, numPoints);
        Eigen::Map<Eigen::VectorXd> vecValues(values, numPoints);
        vecValues = A() * 32.0 * (matPoints.row(1).array() * (1.0 - matPoints.row(1).array()) +
                                  matPoints.row(0).array() * (1.0 - matPoints.row(0).array())) +
                    B() * 16.0 * (1.0 - 2.0 * matPoints.row(0).array()) *
                    matPoints.row(1).array() * (1.0 - matPoints.row(1).array()) +
                    B() * 16.0 * (1.0 - 2.0 * matPoints.row(1).array()) *
                    matPoints.row(0).array() * (1.0 - matPoints.row(0).array()) +
                    C() * 16.0 * (matPoints.row(1).array() *
                                  (1.0 - matPoints.row(1).array()) *
                                  matPoints.row(0).array() *
                                  (1.0 - matPoints.row(0).array())) +
                    C() * 1.1;
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
      static double* WeakTerm_Right(const int numPoints, const double* points)
      {
        double* values = new double[numPoints];

        Eigen::Map<const Eigen::MatrixXd> matPoints(points, 3, numPoints);
        Eigen::Map<Eigen::VectorXd> vecValues(values, numPoints);

        vecValues = A() * 16.0 * (1.0 - 2.0 * matPoints.row(0).array()) *
                    matPoints.row(1).array() *
                    (1.0 - matPoints.row(1).array());

        return values;
      }
      // ***************************************************************************
      static double* WeakTerm_Left(const int numPoints, const double* points)
      {
        double* values = new double[numPoints];

        Eigen::Map<const Eigen::MatrixXd> matPoints(points, 3, numPoints);
        Eigen::Map<Eigen::VectorXd> vecValues(values, numPoints);

        vecValues = - A() * 16.0 * (1.0 - 2.0 * matPoints.row(0).array()) *
                    matPoints.row(1).array() *
                    (1.0 - matPoints.row(1).array());

        return values;
      }
      // ***************************************************************************
  };
}

#endif // __TEST_POISSON_H
