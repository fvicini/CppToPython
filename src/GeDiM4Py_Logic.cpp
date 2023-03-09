#include "GeDiM4Py_Logic.hpp"
#include "MeshMatricesDAO.hpp"
#include "MeshUtilities.hpp"

namespace GedimForPy
{
  // ***************************************************************************
  void InterfaceDataDAO::Construct()
  {
    data.p_geometryUtilitiesConfig = nullptr;
    data.p_geometryUtilities = nullptr;
    data.p_meshUtilities = nullptr;
  }
  // ***************************************************************************
  void InterfaceDataDAO::Destroy()
  {
    delete data.p_meshUtilities;
    delete data.p_geometryUtilitiesConfig;
    delete data.p_geometryUtilities;
  }
  // ***************************************************************************
  GeDiM4Py_Logic::GeDiM4Py_Logic()
  {
  }
  GeDiM4Py_Logic::~GeDiM4Py_Logic()
  {
  }
  // ***************************************************************************
  void GeDiM4Py_Logic::Initialize(const InterfaceConfiguration& config,
                                  InterfaceData& data)
  {
    data.p_geometryUtilitiesConfig = new Gedim::GeometryUtilitiesConfig();
    data.p_geometryUtilitiesConfig->Tolerance = config.GeometricTolerance;
    data.p_geometryUtilities = new Gedim::GeometryUtilities(*data.p_geometryUtilitiesConfig);
    data.p_meshUtilities = new Gedim::MeshUtilities();
  }
  // ***************************************************************************
  Domain2DMesh GeDiM4Py_Logic::CreateDomainMesh2D(const Domain2D& domain,
                                                  InterfaceDataDAO& gedimData)
  {
    Domain2DMesh mesh;
    Gedim::MeshMatricesDAO domainMesh(mesh.Mesh);

    switch (domain.DiscretizationType)
    {
      case Domain2D::DiscretizationTypes::Triangular:
      {
        gedimData.MeshUtilities().CreateTriangularMesh(domain.Vertices,
                                                       domain.MeshCellsMaximumArea,
                                                       domainMesh);
        gedimData.MeshUtilities().ChangePolygonMeshMarkers(domain.Vertices,
                                                           domain.VerticesBoundaryCondition,
                                                           domain.EdgesBoundaryCondition,
                                                           domainMesh);
      }
        break;
      default:
        throw std::runtime_error("MeshGenerator " +
                                 std::to_string((unsigned int)domain.DiscretizationType) +
                                 " not supported");
    }

    return mesh;
  }
  // ***************************************************************************
  DiscreteProblemData GeDiM4Py_Logic::Discretize(const Gedim::IMeshDAO& mesh,
                                                 const DiscreteSpace& space)
  {
    DiscreteProblemData problemData;

    if (space.Type != DiscreteSpace::Types::FEM)
      throw std::runtime_error("DiscreteSpace Type " +
                               std::to_string((unsigned int)space.Type) +
                               " not supported");

    if (space.Order < 1 && space.Order > 2)
      throw std::runtime_error("DiscreteSpace Order " +
                               std::to_string(space.Order) +
                               " not supported");

    problemData.NumberDOFs = 0;
    problemData.NumberStrongs = 0;
    problemData.Cell0Ds_DOF.resize(mesh.Cell0DTotalNumber());
    problemData.Cell1Ds_DOF.resize(mesh.Cell1DTotalNumber());

    for (unsigned int p = 0; p < mesh.Cell0DTotalNumber(); p++)
    {
      const DiscreteSpace::BoundaryConditionTypes& type = space.BoundaryConditionsType.at(mesh.Cell0DMarker(p));
      DiscreteProblemData::DOF& cell0D_DOF =  problemData.Cell0Ds_DOF[p];

      switch (type)
      {
        case DiscreteSpace::BoundaryConditionTypes::None:
        case DiscreteSpace::BoundaryConditionTypes::Weak:
          cell0D_DOF.Type = DiscreteProblemData::DOF::Types::DOF;
          cell0D_DOF.Global_Index = problemData.NumberDOFs;
          problemData.NumberDOFs++;
          break;
        case DiscreteSpace::BoundaryConditionTypes::Strong:
          cell0D_DOF.Type = DiscreteProblemData::DOF::Types::Strong;
          cell0D_DOF.Global_Index = problemData.NumberStrongs;
          problemData.NumberStrongs++;
          break;
        default:
          throw std::runtime_error("BoundaryCondition Type " +
                                   std::to_string((unsigned int)type) +
                                   " not supported");
      }
    }

    if (space.Order == 2)
    {
      for (unsigned int e = 0; e < mesh.Cell1DTotalNumber(); e++)
      {
        const DiscreteSpace::BoundaryConditionTypes& type = space.BoundaryConditionsType.at(mesh.Cell1DMarker(e));
        DiscreteProblemData::DOF& cell1D_DOF =  problemData.Cell1Ds_DOF[e];

        switch (type)
        {
          case DiscreteSpace::BoundaryConditionTypes::None:
          case DiscreteSpace::BoundaryConditionTypes::Weak:
            cell1D_DOF.Type = DiscreteProblemData::DOF::Types::DOF;
            cell1D_DOF.Global_Index = problemData.NumberDOFs;
            problemData.NumberDOFs++;
            break;
          case DiscreteSpace::BoundaryConditionTypes::Strong:
            cell1D_DOF.Type = DiscreteProblemData::DOF::Types::Strong;
            cell1D_DOF.Global_Index = problemData.NumberStrongs;
            problemData.NumberStrongs++;
            break;
          default:
            throw std::runtime_error("BoundaryCondition Type " +
                                     std::to_string((unsigned int)type) +
                                     " not supported");
        }
      }
    }

    FEM_RefElement_Langrange_PCC_Triangle_2D femRefElement;
    problemData.LocalSpace = femRefElement.Compute(space.Order);

    return problemData;
  }
  // ***************************************************************************
}
