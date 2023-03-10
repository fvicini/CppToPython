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
  TEST(TestPython, Test_SquareLaplace)
  {
    const std::string exportFolder = "./Export/TestPython/Test_SquareLaplace";
    Gedim::Output::CreateFolder(exportFolder);

    Py_Initialize();

    GedimForPy::InterfaceConfiguration config_expected;
    config_expected.GeometricTolerance = 1.0e-8;

    PyObject* config = PyDict_New();
    PyDict_SetItemString(config, "GeometricTolerance", Py_BuildValue("d", config_expected.GeometricTolerance));

    ASSERT_NO_THROW(GedimForPy_Initialize(config));
    ASSERT_EQ(config_expected.GeometricTolerance,
              GedimForPy::GeDiM4Py_Interface::InterfaceConfig.GeometricTolerance);

    GedimForPy::InterfaceDataDAO gedimData(GedimForPy::GeDiM4Py_Interface::InterfaceData);

    const double squareEdge = 1.0;

    GedimForPy::Domain2D domain_expected;
    domain_expected.Vertices = gedimData.GeometryUtilities().CreateSquare(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                                          squareEdge);
    domain_expected.VerticesBoundaryCondition = { 1, 1, 1, 1 };
    domain_expected.EdgesBoundaryCondition = { 1, 1, 1, 1 };
    domain_expected.MeshCellsMaximumArea = 0.1;
    domain_expected.DiscretizationType = GedimForPy::Domain2D::DiscretizationTypes::Triangular;

    PyObject* square = PyDict_New();

    PyObject* verticesBoundaryMarker = PyList_New(4);
    PyObject* edgesBoundaryMarker = PyList_New(4);
    for (unsigned int v = 0; v < 4; v++)
    {
      PyList_SET_ITEM(verticesBoundaryMarker, v, Py_BuildValue("i", 1));
      PyList_SET_ITEM(edgesBoundaryMarker, v, Py_BuildValue("i", 1));
    }

    PyDict_SetItemString(square, "SquareEdge", Py_BuildValue("d", squareEdge));
    PyDict_SetItemString(square, "VerticesBoundaryCondition", Py_BuildValue("O", verticesBoundaryMarker));
    PyDict_SetItemString(square, "EdgesBoundaryCondition", Py_BuildValue("O", edgesBoundaryMarker));
    PyDict_SetItemString(square, "DiscretizationType", Py_BuildValue("i", static_cast<int>(domain_expected.DiscretizationType)));
    PyDict_SetItemString(square, "MeshCellsMaximumArea", Py_BuildValue("d", domain_expected.MeshCellsMaximumArea));

    double* meshCoordinates = nullptr;

    ASSERT_NO_THROW(GedimForPy_CreateDomainSquare(square,
                                                  &meshCoordinates));

    delete[] meshCoordinates;

    ASSERT_EQ(domain_expected.Vertices,
              GedimForPy::GeDiM4Py_Interface::Domain.Vertices);
    ASSERT_EQ(domain_expected.VerticesBoundaryCondition,
              GedimForPy::GeDiM4Py_Interface::Domain.VerticesBoundaryCondition);
    ASSERT_EQ(domain_expected.EdgesBoundaryCondition,
              GedimForPy::GeDiM4Py_Interface::Domain.EdgesBoundaryCondition);
    ASSERT_EQ(domain_expected.MeshCellsMaximumArea,
              GedimForPy::GeDiM4Py_Interface::Domain.MeshCellsMaximumArea);
    ASSERT_EQ(domain_expected.DiscretizationType,
              GedimForPy::GeDiM4Py_Interface::Domain.DiscretizationType);
    ASSERT_EQ(domain_expected.VerticesBoundaryCondition,
              GedimForPy::GeDiM4Py_Interface::Domain.VerticesBoundaryCondition);
    ASSERT_EQ(domain_expected.EdgesBoundaryCondition,
              GedimForPy::GeDiM4Py_Interface::Domain.EdgesBoundaryCondition);

    ASSERT_EQ(13, GedimForPy::GeDiM4Py_Interface::Mesh.Mesh.NumberCell0D);
    ASSERT_EQ(28, GedimForPy::GeDiM4Py_Interface::Mesh.Mesh.NumberCell1D);
    ASSERT_EQ(16, GedimForPy::GeDiM4Py_Interface::Mesh.Mesh.NumberCell2D);

    Gedim::MeshMatricesDAO meshDAO(GedimForPy::GeDiM4Py_Interface::Mesh.Mesh);

    // export
    {
      {
        Gedim::VTKUtilities exporter;
        exporter.AddPolygon(GedimForPy::GeDiM4Py_Interface::Domain.Vertices);
        exporter.Export(exportFolder +
                        "/Domain.vtu");
      }

      {
        Gedim::VTKUtilities exporter;
        gedimData.MeshUtilities().ExportMeshToVTU(meshDAO,
                                                  exportFolder,
                                                  "Mesh");
      }
    }

    PyObject* discreteSpace = PyDict_New();

    GedimForPy::DiscreteSpace space_expected;
    space_expected.Order = 2;
    space_expected.Type = GedimForPy::DiscreteSpace::Types::FEM;
    space_expected.BoundaryConditionsType = { GedimForPy::DiscreteSpace::BoundaryConditionTypes::None,
                                              GedimForPy::DiscreteSpace::BoundaryConditionTypes::Strong };

    PyObject* boundaryConditionTypes = PyList_New(2);
    for (unsigned int b = 0; b < space_expected.BoundaryConditionsType.size(); b++)
      PyList_SET_ITEM(boundaryConditionTypes,
                      b,
                      Py_BuildValue("i", static_cast<int>(space_expected.BoundaryConditionsType[b])));

    PyDict_SetItemString(discreteSpace, "Order", Py_BuildValue("i", space_expected.Order));
    PyDict_SetItemString(discreteSpace, "BoundaryConditionsType", Py_BuildValue("O", boundaryConditionTypes));
    PyDict_SetItemString(discreteSpace, "Type", Py_BuildValue("i", static_cast<unsigned int>(space_expected.Type)));

    double* dofsCoordinate = nullptr;
    double* strongsCoordinate = nullptr;

    ASSERT_NO_THROW(GedimForPy_Discretize(discreteSpace,
                                          &dofsCoordinate,
                                          &strongsCoordinate));
    delete[] dofsCoordinate;
    delete[] strongsCoordinate;

    ASSERT_EQ(space_expected.Order,
              GedimForPy::GeDiM4Py_Interface::Space.Order);
    ASSERT_EQ(space_expected.Type,
              GedimForPy::GeDiM4Py_Interface::Space.Type);
    ASSERT_EQ(space_expected.BoundaryConditionsType,
              GedimForPy::GeDiM4Py_Interface::Space.BoundaryConditionsType);

    GedimForPy::DiscreteProblemData& problemData = GedimForPy::GeDiM4Py_Interface::ProblemData;
    ASSERT_EQ(25, problemData.NumberDOFs);
    ASSERT_EQ(16, problemData.NumberStrongs);
    ASSERT_EQ(13, problemData.Cell0Ds_DOF.size());
    ASSERT_EQ(28, problemData.Cell1Ds_DOF.size());

    // export
    {
      {
        std::vector<double> cell0Ds_DOFType(meshDAO.Cell0DTotalNumber(), 0.0);
        std::vector<double> cell0Ds_DOFGlobalIndex(meshDAO.Cell0DTotalNumber(), 0.0);
        std::vector<double> cell1Ds_DOFType(meshDAO.Cell1DTotalNumber(), 0.0);
        std::vector<double> cell1Ds_DOFGlobalIndex(meshDAO.Cell1DTotalNumber(), 0.0);

        for (unsigned int p = 0; p < problemData.Cell0Ds_DOF.size(); p++)
        {
          const GedimForPy::DiscreteProblemData::DOF& dof = problemData.Cell0Ds_DOF[p];
          cell0Ds_DOFType[p] = (unsigned int)dof.Type;
          cell0Ds_DOFGlobalIndex[p] = dof.Global_Index;
        }

        for (unsigned int p = 0; p < problemData.Cell1Ds_DOF.size(); p++)
        {
          const GedimForPy::DiscreteProblemData::DOF& dof = problemData.Cell1Ds_DOF[p];
          cell1Ds_DOFType[p] = (unsigned int)dof.Type;
          cell1Ds_DOFGlobalIndex[p] = dof.Global_Index;
        }

        Gedim::VTKUtilities exporter;
        exporter.AddSegments(meshDAO.Cell0DsCoordinates(),
                             meshDAO.Cell1DsExtremes(),
                             {
                               {
                                 "cell0Ds_DOFType",
                                 Gedim::VTPProperty::Formats::Points,
                                 static_cast<unsigned int>(cell0Ds_DOFType.size()),
                                 cell0Ds_DOFType.data()
                               },
                               {
                                 "cell0Ds_DOFGlobalIndex",
                                 Gedim::VTPProperty::Formats::Points,
                                 static_cast<unsigned int>(cell0Ds_DOFGlobalIndex.size()),
                                 cell0Ds_DOFGlobalIndex.data()
                               },
                               {
                                 "cell1Ds_DOFType",
                                 Gedim::VTPProperty::Formats::Cells,
                                 static_cast<unsigned int>(cell1Ds_DOFType.size()),
                                 cell1Ds_DOFType.data()
                               },
                               {
                                 "cell1Ds_DOFGlobalIndex",
                                 Gedim::VTPProperty::Formats::Cells,
                                 static_cast<unsigned int>(cell1Ds_DOFGlobalIndex.size()),
                                 cell1Ds_DOFGlobalIndex.data()
                               }
                             });
        exporter.Export(exportFolder +
                        "/DOFs.vtu");
      }
    }

    Py_Finalize();
  }
  // ***************************************************************************
}

#endif // __TEST_PYTHON_H
