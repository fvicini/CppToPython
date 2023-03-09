#ifndef __FEM_RefElement_Langrange_PCC_Triangle_2D_H
#define __FEM_RefElement_Langrange_PCC_Triangle_2D_H

#include "Eigen/Eigen"
#include "MapTriangle.hpp"

namespace GedimForPy
{
  /// \brief 1D Primal Conforming Constant Lagrange Element Degree variable
  class FEM_RefElement_Langrange_PCC_Triangle_2D final
  {
    public:
      struct LocalSpace
      {
          unsigned int Order; ///< Order of the method
          unsigned int NumberBasisFunctions; ///< Number of total basis functions
          unsigned int NumberDofs0D; ///< Number of Dofs 0D
          unsigned int NumberDofs1D; ///< Number of Dofs 1D
          unsigned int NumberDofs2D; ///< Number of Dofs 2D
          std::vector<unsigned int> Dof0DsIndex; ///< local DOF index for each element 0D, size num0D + 1
          std::vector<unsigned int> Dof1DsIndex; ///< local DOF index for each element 1D, size num1D + 1
          std::vector<unsigned int> Dof2DsIndex; ///< local DOF index for each element 2D, size num2D + 1
          Eigen::MatrixXd ReferenceElementDofPositions; ///< reference element dof points
          Eigen::MatrixXi ReferenceElementDofTypes; ///< dof type [num oblique edges, num vertical edges, num horizontal edges]
      };

    public:
      FEM_RefElement_Langrange_PCC_Triangle_2D() {}
      ~FEM_RefElement_Langrange_PCC_Triangle_2D() {}

      LocalSpace Compute(const unsigned int& order) const;

      Eigen::MatrixXd Reference_BasisFunctions(const LocalSpace& localSpace,
                                               const Eigen::MatrixXd& points) const;

      std::vector<Eigen::MatrixXd> Reference_BasisFunctionDerivatives(const LocalSpace& localSpace,
                                                                      const Eigen::MatrixXd& points) const;

      inline Eigen::MatrixXd Map_BasisFunctions(const LocalSpace& localSpace,
                                                const Gedim::MapTriangle::MapTriangleData& mapData,
                                                const Eigen::MatrixXd& reference_values) const
      { return reference_values; }

      std::vector<Eigen::MatrixXd> Map_BasisFunctionDerivatives(const LocalSpace& localSpace,
                                                                const Gedim::MapTriangle::MapTriangleData& mapData,
                                                                const std::vector<Eigen::MatrixXd>& reference_values) const;

  };
}

#endif // __FEM_RefElement_Langrange_PCC_Triangle_2D_H
