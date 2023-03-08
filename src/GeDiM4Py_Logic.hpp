#ifndef __GeDiM4Py_Logic_H
#define __GeDiM4Py_Logic_H

#include "Eigen/Eigen"
#include <vector>

#include "MeshMatrices.hpp"
#include "MeshUtilities.hpp"

namespace GedimForPy
{
  struct InterfaceData final
  {
      Gedim::GeometryUtilitiesConfig* p_geometryUtilitiesConfig;
      Gedim::GeometryUtilities* p_geometryUtilities;
      Gedim::MeshUtilities* p_meshUtilities;
  };

  class InterfaceDataDAO final
  {
    private:
      InterfaceData& data;

      void Construct();
      void Destroy();

    public:
      InterfaceDataDAO(InterfaceData& data) : data(data) { Construct(); }
      ~InterfaceDataDAO() { Destroy(); }

      Gedim::GeometryUtilitiesConfig& GeometryUtilitiesConfig() { return *data.p_geometryUtilitiesConfig; }
      Gedim::GeometryUtilities& GeometryUtilities() { return *data.p_geometryUtilities; }
      Gedim::MeshUtilities& MeshUtilities() { return *data.p_meshUtilities; }
      const Gedim::GeometryUtilitiesConfig& GeometryUtilitiesConfig() const { return *data.p_geometryUtilitiesConfig; }
      const Gedim::GeometryUtilities& GeometryUtilities() const { return *data.p_geometryUtilities; }
      const Gedim::MeshUtilities& MeshUtilities() const { return *data.p_meshUtilities; }
  };

  struct GeDiM4Py_Logic_Configuration final
  {
      double GeometricTolerance = 1.0e-8;
  };

  struct Domain2D final
  {
      enum struct DiscretizationTypes
      {
        Unknown = 0,
        Triangular = 1
      };

      Eigen::MatrixXd Vertices;
      std::vector<unsigned int> VerticesBoundaryCondition = {};
      std::vector<unsigned int> EdgesBoundaryCondition = {};
      DiscretizationTypes DiscretizationType = DiscretizationTypes::Unknown;
      double MeshCellsMaximumArea = 0.0;
  };

  struct Domain2DMesh final
  {
      Gedim::MeshMatrices Mesh;
  };

  class GeDiM4Py_Logic final
  {
    private:
      InterfaceData& data;

    public:
      GeDiM4Py_Logic(InterfaceData& data);
      ~GeDiM4Py_Logic();

      void Initialize(const GeDiM4Py_Logic_Configuration& config);

      Domain2DMesh CreateDomainMesh2D(const Domain2D& domain);
  };

}

#endif
