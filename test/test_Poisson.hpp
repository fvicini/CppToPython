#ifndef __TEST_POISSON_H
#define __TEST_POISSON_H

#include "Eigen/Eigen"

namespace UnitTesting
{
  class Poisson final
  {
    public:
      // ***************************************************************************
      static double* DiffusionTerm(const int numPoints, const double* points)
      {
        double* values = new double[numPoints];

        Eigen::Map<Eigen::VectorXd> vecValues(values, numPoints);
        vecValues.setConstant(10.0);

        return values;
      }
      // ***************************************************************************
      static double* AdvectionTerm(const int numPoints, const double* points)
      {
        double* values = new double[2 * numPoints];

        Eigen::Map<Eigen::MatrixXd> matValues(values, 2, numPoints);
        matValues.row(0).setConstant(0.1);
        matValues.row(1).setConstant(0.1);

        return values;
      }
      // ***************************************************************************
      static double* ReactionTerm(const int numPoints, const double* points)
      {
        double* values = new double[numPoints];

        Eigen::Map<Eigen::VectorXd> vecValues(values, numPoints);
        vecValues.setConstant(2.0);

        return values;
      }
      // ***************************************************************************
      static double* ForcingTerm(const int numPoints, const double* points)
      {
        double* values = new double[numPoints];

        Eigen::Map<const Eigen::MatrixXd> matPoints(points, 3, numPoints);
        Eigen::Map<Eigen::VectorXd> vecValues(values, numPoints);
        vecValues = 10.0 * 32.0 * (matPoints.row(1).array() * (1.0 - matPoints.row(1).array()) +
                                   matPoints.row(0).array() * (1.0 - matPoints.row(0).array())) +
                    0.1 * 16.0 * (1.0 - 2.0 * matPoints.row(0).array()) *
                    matPoints.row(1).array() *
                    (1.0 - matPoints.row(1).array()) +
                    0.1 * 6.0 * (1.0 - 2.0 * matPoints.row(1).array()) *
                    matPoints.row(0).array() *
                    (1.0 - matPoints.row(0).array()) +
                    2.0 * 16.0 * (matPoints.row(1).array() *
                                  (1.0 - matPoints.row(1).array()) *
                                  matPoints.row(0).array() *
                                  (1.0 - matPoints.row(0).array())) +
                    2.0 * 1.1;
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
      static double* ExactDerivativeSolution(const unsigned int& direction,
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

        vecValues = 10.0 * 16.0 * (1.0 - 2.0 * matPoints.row(0).array()) *
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

        vecValues = - 10.0 * 16.0 * (1.0 - 2.0 * matPoints.row(0).array()) *
                    matPoints.row(1).array() *
                    (1.0 - matPoints.row(1).array());

        return values;
      }
      // ***************************************************************************
  };
}

#endif // __TEST_POISSON_H
