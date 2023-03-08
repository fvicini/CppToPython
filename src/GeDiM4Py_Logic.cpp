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
  GeDiM4Py_Logic::GeDiM4Py_Logic(InterfaceData& data) :
    data(data)
  {
  }
  GeDiM4Py_Logic::~GeDiM4Py_Logic()
  {
  }
  // ***************************************************************************
  void GeDiM4Py_Logic::Initialize(const GeDiM4Py_Logic_Configuration& config)
  {
    data.p_geometryUtilitiesConfig = new Gedim::GeometryUtilitiesConfig();
    data.p_geometryUtilitiesConfig->Tolerance = config.GeometricTolerance;
    data.p_geometryUtilities = new Gedim::GeometryUtilities(*data.p_geometryUtilitiesConfig);
    data.p_meshUtilities = new Gedim::MeshUtilities();
  }
  // ***************************************************************************
  Domain2DMesh GeDiM4Py_Logic::CreateDomainMesh2D(const Domain2D& domain)
  {
    InterfaceDataDAO gedimData(data);

    Domain2DMesh mesh;
    Gedim::MeshMatricesDAO domainMesh(mesh.Mesh);

    switch (domain.DiscretizationType)
    {
      case Domain2D::DiscretizationTypes::Triangular:
      {
        gedimData.MeshUtilities().CreateTriangularMesh(domain.Vertices,
                                                       domain.MeshCellsMaximumArea,
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
}
