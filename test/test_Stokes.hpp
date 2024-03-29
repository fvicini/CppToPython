#ifndef __TEST_STOKES_H
#define __TEST_STOKES_H

#include "Eigen/Eigen"
#include <iostream>

namespace UnitTesting
{
  class Stokes final
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
      static double* ForcingTerm_1(const int numPoints, const double* points)
      {
        double* values = new double[numPoints];

        Eigen::Map<const Eigen::MatrixXd> matPoints(points, 3, numPoints);
        Eigen::Map<Eigen::VectorXd> vecValues(values, numPoints);
        vecValues = - (+ 8.0 * M_PI * M_PI * cos(4.0 * M_PI * matPoints.row(0).array()) -
                     4.0 * M_PI * M_PI) *
                    sin(2.0 * M_PI * matPoints.row(1).array()) *
                    cos(2.0 * M_PI * matPoints.row(1).array()) +
                    (+ 2.0 * M_PI *
                     cos(2.0 * M_PI * matPoints.row(0).array()) *
                     cos(2.0 * M_PI * matPoints.row(1).array()));

        return values;
      }
      // ***************************************************************************
      static double* ForcingTerm_2(const int numPoints, const double* points)
      {
        double* values = new double[numPoints];

        Eigen::Map<const Eigen::MatrixXd> matPoints(points, 3, numPoints);
        Eigen::Map<Eigen::VectorXd> vecValues(values, numPoints);
        vecValues = - (- 8.0 * M_PI * M_PI * cos(4.0 * M_PI * matPoints.row(1).array()) +
                     4.0 * M_PI * M_PI) *
                    sin(2.0 * M_PI * matPoints.row(0).array()) *
                    cos(2.0 * M_PI * matPoints.row(0).array()) +
                    (- 2.0 * M_PI *
                     sin(2.0 * M_PI * matPoints.row(0).array()) *
                     sin(2.0 * M_PI * matPoints.row(1).array()));

        return values;
      }
      // ***************************************************************************
      static double* ExactPressureSolution(const int numPoints,
                                           const double* points)
      {
        double* values = new double[numPoints];

        Eigen::Map<const Eigen::MatrixXd> matPoints(points, 3, numPoints);
        Eigen::Map<Eigen::VectorXd> vecValues(values, numPoints);

        vecValues = sin(2.0 * M_PI * matPoints.row(0).array()) *
                    cos(2.0 * M_PI * matPoints.row(1).array());
        return values;
      }
      // ***************************************************************************
      static double* ExactSpeedSolution_1(const int numPoints,
                                          const double* points)
      {
        double* values = new double[numPoints];

        Eigen::Map<const Eigen::MatrixXd> matPoints(points, 3, numPoints);
        Eigen::Map<Eigen::VectorXd> vecValues(values, numPoints);

        vecValues = +0.5 *
                    sin(2.0 * M_PI * matPoints.row(0).array()) *
                    sin(2.0 * M_PI * matPoints.row(0).array()) *
                    sin(2.0 * M_PI * matPoints.row(1).array()) *
                    cos(2.0 * M_PI * matPoints.row(1).array());

        return values;
      }
      // ***************************************************************************
      static double* ExactSpeedSolution_2(const int numPoints,
                                          const double* points)
      {
        double* values = new double[numPoints];

        Eigen::Map<const Eigen::MatrixXd> matPoints(points, 3, numPoints);
        Eigen::Map<Eigen::VectorXd> vecValues(values, numPoints);

        vecValues = -0.5 *
                    sin(2.0 * M_PI * matPoints.row(1).array()) *
                    sin(2.0 * M_PI * matPoints.row(1).array()) *
                    sin(2.0 * M_PI * matPoints.row(0).array()) *
                    cos(2.0 * M_PI * matPoints.row(0).array());
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
          vecValues = +2.0 * M_PI *
                      cos(2.0 * M_PI * matPoints.row(0).array()) *
                      cos(2.0 * M_PI * matPoints.row(1).array());
        else if (direction == 1)
          vecValues = -2.0 * M_PI *
                      sin(2.0 * M_PI * matPoints.row(0).array()) *
                      sin(2.0 * M_PI * matPoints.row(1).array());
        else if (direction == 2)
          vecValues.setZero();
        else
          throw std::runtime_error("Error on direction");

        return values;
      }
      // ***************************************************************************
  };
}

#endif // __TEST_STOKES_H
