#ifndef __TEST_BURGER_H
#define __TEST_BURGER_H

#include "Eigen/Eigen"

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
        matValues = Eigen::Map<const Eigen::VectorXd>(u_x,
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
                            matPoints.row(0).array() * (1.0 - matPoints.row(0).array())) +
                    16.0 * ((1.0 - 2.0 * matPoints.row(0).array()) *
                            matPoints.row(1).array() *
                            (1.0 - matPoints.row(1).array()) +
                            (1.0 - 2.0 * matPoints.row(1).array()) *
                            matPoints.row(0).array() *
                            (1.0 - matPoints.row(0).array())) *
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
}

#endif // __TEST_BURGER_H
