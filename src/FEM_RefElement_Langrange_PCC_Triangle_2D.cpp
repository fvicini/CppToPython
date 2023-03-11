#include "FEM_RefElement_Langrange_PCC_Triangle_2D.hpp"
#include "Quadrature_Gauss2D_Triangle.hpp"
#include "Quadrature_Gauss1D.hpp"

using namespace std;
using namespace Eigen;

namespace GedimForPy
{
  // ***************************************************************************
  FEM_RefElement_Langrange_PCC_Triangle_2D::LocalSpace FEM_RefElement_Langrange_PCC_Triangle_2D::Compute(const unsigned int& order) const
  {
    if (order < 1 && order > 2)
      throw std::runtime_error("DiscreteSpace Order " +
                               std::to_string(order) +
                               " not supported");

    LocalSpace localSpace;

    localSpace.Order = order;

    if (order == 0)
    {
      localSpace.NumberDofs0D = 0;
      localSpace.NumberDofs1D = 0;
      localSpace.NumberDofs2D = 1;
      localSpace.NumberBasisFunctions = 1;
      localSpace.ReferenceElement.DofPositions.setZero(3, localSpace.NumberBasisFunctions);
      localSpace.ReferenceElement.DofTypes.setZero(3, localSpace.NumberBasisFunctions);

      localSpace.ReferenceElement.DofPositions.col(0)<< 1.0 / 3.0, 1.0 / 3.0, 0.0;
      return localSpace;
    }

    localSpace.NumberBasisFunctions = (order + 1) * (order + 2) / 2;

    vector<unsigned int> nodeDofs = { 0, order, localSpace.NumberBasisFunctions - 1 };
    list<unsigned int> cellDofs;
    vector<list<unsigned int>> edgeDofs(3);

    Eigen::MatrixXd localDofPositions = Eigen::MatrixXd::Zero(3, localSpace.NumberBasisFunctions);
    Eigen::MatrixXi localDofTypes = Eigen::MatrixXi::Zero(3, localSpace.NumberBasisFunctions);

    const double h = 1.0 / order;
    unsigned int dof = 0;
    for (unsigned int i = 0; i < order + 1; i++)
    {
      for (unsigned int j = 0; j < order + 1 - i; j++)
      {
        if (i == 0)
        {
          if (j > 0 && j < order - i)
            edgeDofs[0].push_back(dof);
        }
        else if (i < order && (j == 0 || j == order - i))
        {
          if (j == 0)
            edgeDofs[2].push_front(dof);
          else
            edgeDofs[1].push_back(dof);
        }
        else if (i < order)
          cellDofs.push_back(dof);

        localDofPositions.col(dof)<< (double)j * h, (double)i * h, 0.0;
        localDofTypes.col(dof)<< order - i - j, j, i;

        dof++;
      }
    }

    // Reordering Dofs using convention [point, edge, cell]
    localSpace.NumberDofs0D = nodeDofs.size() / 3;
    localSpace.NumberDofs1D = edgeDofs[0].size();
    localSpace.NumberDofs2D = cellDofs.size();
    localSpace.ReferenceElement.DofPositions.setZero(3, localSpace.NumberBasisFunctions);
    localSpace.ReferenceElement.DofTypes.setZero(3, localSpace.NumberBasisFunctions);

    dof = 0;
    for (const unsigned int dofIndex : nodeDofs)
    {
      localSpace.ReferenceElement.DofPositions.col(dof)<< localDofPositions.col(dofIndex);
      localSpace.ReferenceElement.DofTypes.col(dof)<< localDofTypes.col(dofIndex);
      dof++;
    }
    for (unsigned int e = 0; e < 3; e++)
    {
      for (const unsigned int dofIndex : edgeDofs.at(e))
      {
        localSpace.ReferenceElement.DofPositions.col(dof)<< localDofPositions.col(dofIndex);
        localSpace.ReferenceElement.DofTypes.col(dof)<< localDofTypes.col(dofIndex);
        dof++;
      }
    }
    for (const unsigned int dofIndex : cellDofs)
    {
      localSpace.ReferenceElement.DofPositions.col(dof)<< localDofPositions.col(dofIndex);
      localSpace.ReferenceElement.DofTypes.col(dof)<< localDofTypes.col(dofIndex);
      dof++;
    }

    localSpace.Dof0DsIndex.resize(4, 0);
    for (unsigned int v = 0; v < 3; v++)
      localSpace.Dof0DsIndex[v + 1] = localSpace.Dof0DsIndex[v] + localSpace.NumberDofs0D;

    localSpace.Dof1DsIndex.resize(4, localSpace.Dof0DsIndex[3]);
    for (unsigned int e = 0; e < 3; e++)
      localSpace.Dof1DsIndex[e + 1] = localSpace.Dof1DsIndex[e] + localSpace.NumberDofs1D;

    Gedim::Quadrature_Gauss2D_Triangle::FillPointsAndWeights(2 * order,
                                                             localSpace.ReferenceElement.InternalQuadrature.Points,
                                                             localSpace.ReferenceElement.InternalQuadrature.Weights);
    Gedim::Quadrature_Gauss1D::FillPointsAndWeights(2 * order,
                                                    localSpace.ReferenceElement.BorderQuadrature.Points,
                                                    localSpace.ReferenceElement.BorderQuadrature.Weights);

    return localSpace;
  }
  // ***************************************************************************
  MatrixXd FEM_RefElement_Langrange_PCC_Triangle_2D::Reference_BasisFunctions(const LocalSpace& localSpace,
                                                                              const Eigen::MatrixXd& points) const
  {
    switch (localSpace.Order)
    {
      case 0:
        return Eigen::VectorXd::Constant(points.cols(), 1.0);
      default:
      {
        const double h = 1.0 / localSpace.Order;
        const Eigen::ArrayXd x = points.row(0).transpose().array();
        const Eigen::ArrayXd y = points.row(1).transpose().array();
        Eigen::MatrixXd values = Eigen::MatrixXd::Ones(points.cols(),
                                                       localSpace.NumberBasisFunctions);

        for (unsigned int d = 0; d < localSpace.NumberBasisFunctions; d++)
        {
          const Eigen::Vector3i& dofType = localSpace.ReferenceElement.DofTypes.col(d);
          const Eigen::Vector3d& dofPosition = localSpace.ReferenceElement.DofPositions.col(d);

          // terms of equation 1 - x - y - t * h
          for (unsigned int t = 0; t < dofType[0]; t++)
          {
            values.col(d).array() *= (1.0 - x - y - t * h);
            values.col(d) /= (1.0 - dofPosition.x() - dofPosition.y() - t * h);
          }

          // terms of equation x - t * h
          for (unsigned int t = 0; t < dofType[1]; t++)
          {
            values.col(d).array() *= (x - t * h);
            values.col(d) /= (dofPosition.x() - t * h);
          }

          // terms of equation y - t * h
          for (unsigned int t = 0; t < dofType[2]; t++)
          {
            values.col(d).array() *= (y  - t * h);
            values.col(d) /= (dofPosition.y() - t * h);
          }
        }
        return values;
      }
    }
  }
  // ***************************************************************************
  std::vector<MatrixXd> FEM_RefElement_Langrange_PCC_Triangle_2D::Reference_BasisFunctionDerivatives(const LocalSpace& localSpace, const Eigen::MatrixXd& points) const
  {
    switch (localSpace.Order)
    {
      case 0:
        return vector<MatrixXd>(2, MatrixXd::Zero(points.cols(), localSpace.NumberBasisFunctions));
      default:
      {
        const double h = 1.0 / localSpace.Order;
        const Eigen::ArrayXd x = points.row(0).transpose().array();
        const Eigen::ArrayXd y = points.row(1).transpose().array();
        vector<MatrixXd> gradValues(2, MatrixXd::Zero(points.cols(), localSpace.NumberBasisFunctions));

        for (unsigned int d = 0; d < localSpace.NumberBasisFunctions; d++)
        {
          const Eigen::Vector3i& dofType = localSpace.ReferenceElement.DofTypes.col(d);
          const Eigen::Vector3d& dofPosition = localSpace.ReferenceElement.DofPositions.col(d);

          const unsigned int numProds = dofType[0] +
                                        dofType[1] +
                                        dofType[2];

          vector<Eigen::ArrayXd> prod_terms = vector<Eigen::ArrayXd>(numProds);
          vector<Eigen::Array2d> grad_terms = vector<Eigen::Array2d>(numProds);
          double denominator = 1.0;

          unsigned int dt = 0;
          // terms of equation 1 - x - y - t * h
          for (unsigned int t = 0; t < dofType[0]; t++)
          {
            prod_terms[dt] = (1.0 - x - y - t * h);
            grad_terms[dt]<< -1.0, -1.0;
            denominator *= (1.0 - dofPosition.x() - dofPosition.y() - t * h);
            dt++;
          }

          // terms of equation x - t * h
          for (unsigned int t = 0; t < dofType[1]; t++)
          {
            prod_terms[dt] = (x - t * h);
            grad_terms[dt]<< 1.0, 0.0;
            denominator *= (dofPosition.x() - t * h);
            dt++;
          }

          // terms of equation y - t * h
          for (unsigned int t = 0; t < dofType[2]; t++)
          {
            prod_terms[dt] = (y  - t * h);
            grad_terms[dt]<< 0.0, 1.0;
            denominator *= (dofPosition.y() - t * h);
            dt++;
          }

          for (unsigned int i = 0; i < numProds; i++)
          {
            Eigen::ArrayXd inner_prod = Eigen::ArrayXd::Ones(points.cols());
            for (unsigned int j = 0; j < numProds; j++)
            {
              if (i != j)
                inner_prod *= prod_terms[j];
            }

            gradValues[0].col(d).array() += inner_prod * grad_terms[i][0];
            gradValues[1].col(d).array() += inner_prod * grad_terms[i][1];
          }

          gradValues[0].col(d) /= denominator;
          gradValues[1].col(d) /= denominator;
        }

        return gradValues;
      }
    }
  }
  // ***************************************************************************
  MatrixXd FEM_RefElement_Langrange_PCC_Triangle_2D::BasisFunctionsOnPoints(const LocalSpace& localSpace,
                                                                            const Gedim::MapTriangle::MapTriangleData& mapData,
                                                                            const Eigen::MatrixXd& points) const
  {
    Gedim::MapTriangle mapTriangle;
    return Reference_BasisFunctions(localSpace,
                                    mapTriangle.FInv(mapData,
                                                     points));

  }
  // ***************************************************************************
  std::vector<MatrixXd> FEM_RefElement_Langrange_PCC_Triangle_2D::BasisFunctionDerivatives(const LocalSpace& localSpace,
                                                                                           const Gedim::MapTriangle::MapTriangleData& mapData,
                                                                                           const std::vector<Eigen::MatrixXd>& reference_values) const
  {
    std::vector<Eigen::MatrixXd> basisFunctionsDerivativeValues(2,
                                                                Eigen::MatrixXd::Zero(reference_values[0].rows(),
                                                                localSpace.NumberBasisFunctions));

    for(unsigned int i = 0; i < 2; i++)
    {
      basisFunctionsDerivativeValues[i] = mapData.BInv(i,i) * reference_values[i];
      for(unsigned int j = 0; j < i; j++)
      {
        basisFunctionsDerivativeValues[i] += mapData.BInv(j,i) * reference_values[j];
        basisFunctionsDerivativeValues[j] += mapData.BInv(i,j) * reference_values[i];
      }
    }

    return basisFunctionsDerivativeValues;
  }
  // ***************************************************************************
  std::vector<MatrixXd> FEM_RefElement_Langrange_PCC_Triangle_2D::BasisFunctionDerivativesOnPoints(const LocalSpace& localSpace,
                                                                                                   const Gedim::MapTriangle::MapTriangleData& mapData,
                                                                                                   const Eigen::MatrixXd& points) const
  {
    Gedim::MapTriangle mapTriangle;
    return BasisFunctionDerivatives(localSpace,
                                    mapData,
                                    Reference_BasisFunctionDerivatives(localSpace,
                                                                       mapTriangle.FInv(mapData,
                                                                                        points)));
  }
  // ***************************************************************************
}
