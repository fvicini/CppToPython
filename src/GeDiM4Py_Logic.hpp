#ifndef __GeDiM4Py_Logic_H
#define __GeDiM4Py_Logic_H

#include "Eigen/Eigen"
#include <vector>

#include "MeshMatrices.hpp"
#include "MeshUtilities.hpp"
#include "FEM_RefElement_Langrange_PCC_Triangle_2D.hpp"

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

  struct InterfaceConfiguration final
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

  struct DiscreteSpace final
  {
      enum struct Types
      {
        Unknown = 0,
        FEM = 1
      };

      enum struct BoundaryConditionTypes
      {
        Unknown = 0,
        None = 1,
        Strong = 2,
        Weak = 3
      };

      unsigned int Order = 0;
      Types Type = Types::Unknown;
      std::vector<BoundaryConditionTypes> BoundaryConditionsType = {};
  };

  struct DiscreteProblemData final
  {
      struct DOF
      {
          enum struct Types
          {
            Unknwon = 0,
            Strong = 1,
            DOF = 2
          };

          Types Type = Types::Unknwon;
          unsigned int Global_Index = 0;
      };

      unsigned int NumberDOFs = 0;
      unsigned int NumberStrongs = 0;
      std::vector<DOF> Cell0Ds_DOF = {};
      std::vector<DOF> Cell1Ds_DOF = {};
      std::vector<std::vector<DOF*>> Cell2Ds_DOF = {};
      FEM_RefElement_Langrange_PCC_Triangle_2D::LocalSpace LocalSpace;
  };

  class GeDiM4Py_Logic final
  {
    public:
      typedef const Eigen::VectorXd (*K)(const Eigen::MatrixXd& points);

    private:

    public:
      GeDiM4Py_Logic();
      ~GeDiM4Py_Logic();

      static void Initialize(const InterfaceConfiguration& config,
                             InterfaceData& data);

      static Domain2DMesh CreateDomainMesh2D(const Domain2D& domain,
                                             InterfaceDataDAO& gedimData);

      static DiscreteProblemData Discretize(const Gedim::IMeshDAO& mesh,
                                            const DiscreteSpace& space);

      static std::list<Eigen::Triplet<double>> AssembleStiffnessMatrix(K k,
                                                                       const Gedim::IMeshDAO& mesh,
                                                                       const std::vector<Gedim::MapTriangle::MapTriangleData>& cell2DsMap,
                                                                       const DiscreteProblemData& problemData);
  };

}

#endif
