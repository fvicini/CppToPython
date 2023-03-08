#ifndef __TEST_PYTHON_H
#define __TEST_PYTHON_H

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>

#include "GeDiM4Py_Interface.hpp"
#include "MeshMatricesDAO.hpp"
#include "VTKUtilities.hpp"

namespace UnitTesting
{
  // ***************************************************************************
  TEST(TestPython, Test_Initialize)
  {
    Py_Initialize();

    PyObject* dict = PyDict_New();
    PyDict_SetItemString(dict, "GeometricTolerance", Py_BuildValue("d", 1.0e-8));

    ASSERT_NO_THROW(GedimForPy_Initialize(dict));

    Py_Finalize();
  }
  // ***************************************************************************
  TEST(TestPython, Test_CreateDomainSquare)
  {
    Py_Initialize();

    const std::string exportFolder = "./Export/TestPython/Test_CreateDomainSquare";
    Gedim::Output::CreateFolder(exportFolder);

    GedimForPy::InterfaceDataDAO gedimData(GedimForPy::GeDiM4Py_Interface::InterfaceData);

    const double squareEdge = 1.0;

    GedimForPy::Domain2D domain_expected;
    domain_expected.Vertices = gedimData.GeometryUtilities().CreateSquare(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                                          squareEdge);
    domain_expected.MeshCellsMaximumArea = 0.1;
    domain_expected.DiscretizationType = GedimForPy::Domain2D::DiscretizationTypes::Triangular;
    const int discretizationType = static_cast<int>(domain_expected.DiscretizationType);

    PyObject* dict = PyDict_New();

    PyDict_SetItemString(dict, "SquareEdge", Py_BuildValue("d", squareEdge));
    PyDict_SetItemString(dict, "DiscretizationType", Py_BuildValue("i", discretizationType));
    PyDict_SetItemString(dict, "MeshCellsMaximumArea", Py_BuildValue("d", domain_expected.MeshCellsMaximumArea));

    ASSERT_NO_THROW(GedimForPy_CreateDomainSquare(dict));
    ASSERT_EQ(domain_expected.Vertices,
              GedimForPy::GeDiM4Py_Interface::Domain.Vertices);
    ASSERT_EQ(domain_expected.MeshCellsMaximumArea,
              GedimForPy::GeDiM4Py_Interface::Domain.MeshCellsMaximumArea);
    ASSERT_EQ(domain_expected.DiscretizationType,
              GedimForPy::GeDiM4Py_Interface::Domain.DiscretizationType);
    ASSERT_EQ(domain_expected.VerticesBoundaryCondition,
              GedimForPy::GeDiM4Py_Interface::Domain.VerticesBoundaryCondition);
    ASSERT_EQ(domain_expected.EdgesBoundaryCondition,
              GedimForPy::GeDiM4Py_Interface::Domain.EdgesBoundaryCondition);

    // export
    {
      Gedim::MeshMatricesDAO domainMesh(GedimForPy::GeDiM4Py_Interface::Mesh.Mesh);

      {
        Gedim::VTKUtilities exporter;
        exporter.AddPolygon(GedimForPy::GeDiM4Py_Interface::Domain.Vertices);
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

    Py_Finalize();
  }
  // ***************************************************************************
}

#endif // __TEST_PYTHON_H
