#ifndef __TEST_NONLINEAR_POISSON_H
#define __TEST_NONLINEAR_POISSON_H

#include "Eigen/Eigen"

namespace UnitTesting
{
  class NonLinearPoisson final
  {
    public:
      static constexpr double MU_1() { return 3.34; }
      static constexpr double MU_2() { return 4.45; }
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
        matValues = MU_1() / MU_2() *
                    ((MU_2() *
                      Eigen::Map<const Eigen::VectorXd>(u,
                                                        numPoints).array()).exp() - 1.0);

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
        matValues = MU_1() *
                    (MU_2() *
                     Eigen::Map<const Eigen::VectorXd>(u,
                                                       numPoints).array()).exp();

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
        vecValues = 100.0 *
                    (2.0 * M_PI * matPoints.row(0).array()).sin() *
                    (2.0 * M_PI * matPoints.row(1).array()).cos();
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
  };
}

#endif // __TEST_NONLINEAR_POISSON_H
