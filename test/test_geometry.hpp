#ifndef __TEST_GEOMETRY_H
#define __TEST_GEOMETRY_H

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>

#include "GeDiM4Py_Logic.hpp"
#include "MeshMatricesDAO.hpp"
#include "VTKUtilities.hpp"

namespace UnitTesting
{
  // ***************************************************************************
  TEST(TestGeometry, Test_Initialize)
  {
    const std::string exportFolder = "./Export";
    Gedim::Output::CreateFolder(exportFolder);

    GedimForPy::InterfaceConfiguration interfaceConfig;
    interfaceConfig.GeometricTolerance = 1.0e-8;

    GedimForPy::InterfaceData data;
    GedimForPy::InterfaceDataDAO gedimData(data);

    GedimForPy::GeDiM4Py_Logic interface;

    ASSERT_NO_THROW(interface.Initialize(interfaceConfig,
                                         data));
  }
  // ***************************************************************************
  TEST(TestGeometry, Test_CreateDomain_Square)
  {
    const std::string exportFolder = "./Export/TestGeometry/Test_CreateDomain_Square";
    Gedim::Output::CreateFolder(exportFolder);

    GedimForPy::InterfaceConfiguration interfaceConfig;
    interfaceConfig.GeometricTolerance = 1.0e-8;

    GedimForPy::InterfaceData data;
    GedimForPy::InterfaceDataDAO gedimData(data);

    GedimForPy::GeDiM4Py_Logic::Initialize(interfaceConfig,
                                           data);

    GedimForPy::Domain2D domain;
    domain.Vertices = gedimData.GeometryUtilities().CreateSquare(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                                 1.0);
    domain.VerticesBoundaryCondition = std::vector<unsigned int> { 1, 1, 1, 1 };
    domain.EdgesBoundaryCondition = std::vector<unsigned int> { 1, 1, 1, 1 };
    domain.DiscretizationType = GedimForPy::Domain2D::DiscretizationTypes::Triangular;
    domain.MeshCellsMaximumArea = 0.1;

    GedimForPy::Domain2DMesh mesh = GedimForPy::GeDiM4Py_Logic::CreateDomainMesh2D(domain,
                                                                                   gedimData);

    ASSERT_EQ(13, mesh.Mesh.NumberCell0D);
    ASSERT_EQ(28, mesh.Mesh.NumberCell1D);
    ASSERT_EQ(16, mesh.Mesh.NumberCell2D);

    // export
    {
      Gedim::MeshMatricesDAO domainMesh(mesh.Mesh);

      {
        Gedim::VTKUtilities exporter;
        exporter.AddPolygon(domain.Vertices);
        exporter.Export(exportFolder +
                        "/Domain.vtu");
      }

      {
        Gedim::VTKUtilities exporter;
        gedimData.MeshUtilities().ExportMeshToVTU(domainMesh,
                                                  exportFolder,
                                                  "Mesh");
      }
    }
  }
  // ***************************************************************************
  TEST(TestGeometry, Test_CreateDofs)
  {
    const std::string exportFolder = "./Export/TestGeometry/Test_CreateDomain_Square";
    Gedim::Output::CreateFolder(exportFolder);

    GedimForPy::InterfaceConfiguration interfaceConfig;
    interfaceConfig.GeometricTolerance = 1.0e-8;

    GedimForPy::InterfaceData data;
    GedimForPy::InterfaceDataDAO gedimData(data);

    GedimForPy::GeDiM4Py_Logic::Initialize(interfaceConfig,
                                           data);

    GedimForPy::Domain2D domain;
    domain.Vertices = gedimData.GeometryUtilities().CreateSquare(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                                 1.0);
    domain.VerticesBoundaryCondition = std::vector<unsigned int> { 1, 1, 1, 1 };
    domain.EdgesBoundaryCondition = std::vector<unsigned int> { 1, 1, 1, 1 };
    domain.DiscretizationType = GedimForPy::Domain2D::DiscretizationTypes::Triangular;
    domain.MeshCellsMaximumArea = 0.1;

    GedimForPy::Domain2DMesh mesh = GedimForPy::GeDiM4Py_Logic::CreateDomainMesh2D(domain,
                                                                                   gedimData);
    Gedim::MeshMatricesDAO meshDAO(mesh.Mesh);

    GedimForPy::DiscreteSpace discreteSpace;
    discreteSpace.Type = GedimForPy::DiscreteSpace::Types::FEM;
    discreteSpace.Order = 1;
    discreteSpace.BoundaryConditionsType = { GedimForPy::DiscreteSpace::BoundaryConditionTypes::None,
                                             GedimForPy::DiscreteSpace::BoundaryConditionTypes::Strong };

    GedimForPy::DiscreteProblemData problemData = GedimForPy::GeDiM4Py_Logic::Discretize(meshDAO,
                                                                                         discreteSpace);

  }
  // ***************************************************************************
}

#endif // __TEST_GEOMETRY_H
