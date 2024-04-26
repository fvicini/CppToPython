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
      Gedim::GeometryUtilitiesConfig* p_geometryUtilitiesConfig = nullptr;
      Gedim::GeometryUtilities* p_geometryUtilities = nullptr;
      Gedim::MeshUtilities* p_meshUtilities = nullptr;
  };

  class InterfaceDataDAO final
  {
    private:
      InterfaceData& data;

    public:
      InterfaceDataDAO(InterfaceData& data) : data(data) { }
      ~InterfaceDataDAO() { }

      void Destroy();
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

  struct ImportMesh2D final
  {
      std::string InputFolder = "";
      char Separator = ';';
  };

  struct ExportData final
  {
      std::string ExportFolder;
      std::string FileName;
  };

  struct Domain2DMesh final
  {
      Gedim::MeshMatrices Mesh;
      std::vector<Gedim::MapTriangle::MapTriangleData> Cell2DsMap;
      Gedim::MeshUtilities::MeshGeometricData2D MeshGeometricData;
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

          struct BoundaryInfo
          {
              enum struct Types
              {
                Unknwon = 0,
                Strong = 1,
                Weak = 2,
                None = 3
              };

              Types Type = Types::Unknwon;
              unsigned int Marker = 0;
          };

          Types Type = Types::Unknwon;
          unsigned int Global_Index = 0;
          BoundaryInfo Boundary;
      };

      double H = 0.0;
      unsigned int NumberDOFs = 0;
      unsigned int NumberStrongs = 0;
      Eigen::MatrixXd DOFsCoordinate;
      Eigen::MatrixXd StrongsCoordinate;
      std::vector<DOF> Cell0Ds_DOF = {};
      std::vector<DOF> Cell1Ds_DOF = {};
      std::vector<std::vector<DOF*>> Cell2Ds_DOF = {};
      FEM_RefElement_Langrange_PCC_Triangle_2D::LocalSpace LocalSpace;
  };

  class GeDiM4Py_Logic final
  {
    public:
      typedef double* (*A)(const int numPoints, const double* points);
      typedef double* (*B)(const int numPoints, const double* points);
      typedef double* (*C)(const int numPoints, const double* points);
      typedef double* (*F)(const int numPoints, const double* points);
      typedef double* (*NNL)(const int numPoints, const double* points, const double* u, const double* u_x, const double* u_y);
      typedef double* (*Strong)(const int numPoints, const double* points);
      typedef double* (*Weak)(const int numPoints, const double* points);
      typedef double* (*Exact)(const int numPoints, const double* points);
      typedef double* (*ExactDerivative)(const int direction, const int numPoints, const double* points);

    private:

    public:
      GeDiM4Py_Logic();
      ~GeDiM4Py_Logic();

      static void Initialize(const InterfaceConfiguration& config,
                             InterfaceData& data);

      static Domain2DMesh CreateDomainMesh2D(const Domain2D& domain,
                                             InterfaceDataDAO& gedimData);

      static Domain2DMesh ImportDomainMesh2D(const ImportMesh2D& domain,
                                             InterfaceDataDAO& gedimData,
                                             const bool& checkMesh = false);

      static DiscreteProblemData Discretize(const Gedim::IMeshDAO& mesh,
                                            const Gedim::MeshUtilities::MeshGeometricData2D& meshGeometricData,
                                            const DiscreteSpace& space);

      static void AssembleStiffnessMatrix(A a,
                                          const Gedim::IMeshDAO& mesh,
                                          const std::vector<Gedim::MapTriangle::MapTriangleData>& cell2DsMap,
                                          const DiscreteProblemData& problemData,
                                          std::list<Eigen::Triplet<double> >& stiffnessTriplets,
                                          std::list<Eigen::Triplet<double> >& stiffnessStrongTriplets);
      static void AssembleNonLinearStiffnessMatrix(A a,
                                                   NNL non_linear_f,
                                                   const Gedim::IMeshDAO& mesh,
                                                   const std::vector<Gedim::MapTriangle::MapTriangleData>& cell2DsMap,
                                                   const DiscreteProblemData& problemData,
                                                   const Eigen::VectorXd& numeric_k,
                                                   const Eigen::VectorXd& strong_k,
                                                   std::list<Eigen::Triplet<double> >& stiffnessTriplets,
                                                   std::list<Eigen::Triplet<double> >& stiffnessStrongTriplets);
      static void AssembleAnisotropicStiffnessMatrix(A a,
                                                     const Gedim::IMeshDAO& mesh,
                                                     const std::vector<Gedim::MapTriangle::MapTriangleData>& cell2DsMap,
                                                     const DiscreteProblemData& problemData,
                                                     std::list<Eigen::Triplet<double> >& stiffnessTriplets,
                                                     std::list<Eigen::Triplet<double> >& stiffnessStrongTriplets);
      static void AssembleAdvectionMatrix(B b,
                                          const Gedim::IMeshDAO& mesh,
                                          const std::vector<Gedim::MapTriangle::MapTriangleData>& cell2DsMap,
                                          const DiscreteProblemData& trial_Functions,
                                          const DiscreteProblemData& test_Functions,
                                          std::list<Eigen::Triplet<double>>& advectionTriplets,
                                          std::list<Eigen::Triplet<double>>& advectionStrongTriplets);
      static void AssembleNonLinearAdvectionMatrix(B b,
                                                   NNL non_linear_f,
                                                   const Gedim::IMeshDAO& mesh,
                                                   const std::vector<Gedim::MapTriangle::MapTriangleData>& cell2DsMap,
                                                   const DiscreteProblemData& trial_Functions,
                                                   const DiscreteProblemData& test_Functions,
                                                   const Eigen::VectorXd& numeric_k,
                                                   const Eigen::VectorXd& strong_k,
                                                   std::list<Eigen::Triplet<double>>& advectionTriplets,
                                                   std::list<Eigen::Triplet<double>>& advectionStrongTriplets);
      static void AssembleReactionMatrix(C c,
                                         const Gedim::IMeshDAO& mesh,
                                         const std::vector<Gedim::MapTriangle::MapTriangleData>& cell2DsMap,
                                         const DiscreteProblemData& problemData,
                                         std::list<Eigen::Triplet<double>>& reactionTriplets,
                                         std::list<Eigen::Triplet<double>>& reactionStrongTriplets);
      static void AssembleNonLinearReactionMatrix(C c,
                                                  NNL non_linear_f,
                                                  const Gedim::IMeshDAO& mesh,
                                                  const std::vector<Gedim::MapTriangle::MapTriangleData>& cell2DsMap,
                                                  const DiscreteProblemData& problemData,
                                                  const Eigen::VectorXd& numeric_k,
                                                  const Eigen::VectorXd& strong_k,
                                                  std::list<Eigen::Triplet<double>>& reactionTriplets,
                                                  std::list<Eigen::Triplet<double>>& reactionStrongTriplets);


      static Eigen::VectorXd AssembleForcingTerm(F f,
                                                 const Gedim::IMeshDAO& mesh,
                                                 const std::vector<Gedim::MapTriangle::MapTriangleData>& cell2DsMap,
                                                 const DiscreteProblemData& problemData);
      static Eigen::VectorXd AssembleNonLinearForcingTerm(F f,
                                                          NNL non_linear_f,
                                                          const Gedim::IMeshDAO& mesh,
                                                          const std::vector<Gedim::MapTriangle::MapTriangleData>& cell2DsMap,
                                                          const DiscreteProblemData& problemData,
                                                          const Eigen::VectorXd& numeric_k,
                                                          const Eigen::VectorXd& strong_k);
      static Eigen::VectorXd AssembleNonLinearDerivativeForcingTerm(F f,
                                                                    NNL non_linear_f,
                                                                    const Gedim::IMeshDAO& mesh,
                                                                    const std::vector<Gedim::MapTriangle::MapTriangleData>& cell2DsMap,
                                                                    const DiscreteProblemData& problemData,
                                                                    const Eigen::VectorXd& numeric_k,
                                                                    const Eigen::VectorXd& strong_k);

      static Eigen::VectorXd AssembleStrongSolution(Strong g,
                                                    const unsigned int& marker,
                                                    const Gedim::IMeshDAO& mesh,
                                                    const std::vector<Gedim::MapTriangle::MapTriangleData>& cell2DsMap,
                                                    const DiscreteProblemData& problemData);

      static Eigen::VectorXd AssembleWeakTerm(Weak g,
                                              const unsigned int& marker,
                                              const Gedim::IMeshDAO& mesh,
                                              const std::vector<Eigen::MatrixXd>& cell2DsVertices,
                                              const std::vector<Eigen::VectorXd>& cell2DsEdgeLengths,
                                              const std::vector<Eigen::MatrixXd>& cell2DsEdgeTangents,
                                              const std::vector<Gedim::MapTriangle::MapTriangleData>& cell2DsMap,
                                              const DiscreteProblemData& problemData);

      static Eigen::VectorXd ComputeErrorL2(Exact u,
                                            const Eigen::VectorXd& numeric,
                                            const Eigen::VectorXd& strong,
                                            const Gedim::IMeshDAO& mesh,
                                            const std::vector<Gedim::MapTriangle::MapTriangleData>& cell2DsMap,
                                            const DiscreteProblemData& problemData);
      static Eigen::VectorXd ComputeErrorH1(ExactDerivative uDer,
                                            const Eigen::VectorXd& numeric,
                                            const Eigen::VectorXd& strong,
                                            const Gedim::IMeshDAO& mesh,
                                            const std::vector<Gedim::MapTriangle::MapTriangleData>& cell2DsMap,
                                            const DiscreteProblemData& problemData);

      static void ExportSolution(Exact u,
                                 const Eigen::VectorXd& numeric,
                                 const Eigen::VectorXd& strong,
                                 const Gedim::IMeshDAO& mesh,
                                 const DiscreteProblemData& problemData,
                                 const ExportData& exportData);
  };

}

#endif
