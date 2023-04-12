#ifndef __TEST_HEAT_CONDUCTIVITY_H
#define __TEST_HEAT_CONDUCTIVITY_H

#include "Eigen/Eigen"

namespace UnitTesting
{
  class HeatConductivity final
  {
    public:
      static constexpr double R() { return 0.5; }
      static constexpr double K() { return 6.68; }
      static constexpr double G() { return 0.94; }
      // ***************************************************************************
      static double* DiffusionTerm(const int numPoints, const double* points)
      {
        double* values = new double[3 * numPoints];

        Eigen::Map<const Eigen::MatrixXd> matPoints(points, 3, numPoints);
        Eigen::Map<Eigen::MatrixXd> matValues(values, 3, numPoints);
        matValues.row(0).setOnes();
        matValues.row(1).setZero();
        matValues.row(2).setOnes();

        Eigen::Map<Eigen::VectorXd> vecValues(values, numPoints);

        for (unsigned int p = 0; p < numPoints; p++)
        {
          if (matPoints(0, p) * matPoints(0, p) + matPoints(1, p) * matPoints(1, p) <= (R() * R() + 1.0e-16))
          {
            matValues(0, p) = K();
            matValues(2, p) = K();
          }
        }

        return values;
      }
      // ***************************************************************************
      static double* WeakTerm_Down(const int numPoints, const double* points)
      {
        double* values = new double[numPoints];

        Eigen::Map<Eigen::VectorXd> vecValues(values, numPoints);
        vecValues.setConstant(G());

        return values;
      }
      // ***************************************************************************
  };
}

#endif // __TEST_HEAT_CONDUCTIVITY_H
