#include "GeDiM4Py_Logic.hpp"
#include "MeshMatricesDAO.hpp"
#include "MeshUtilities.hpp"
#include "PDE_Equation.hpp"

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

    Gedim::MapTriangle mapTriangle;
    mesh.Cell2DsMap.resize(domainMesh.Cell2DTotalNumber());
    for (unsigned int c = 0; c < domainMesh.Cell2DTotalNumber(); c++)
      mesh.Cell2DsMap.at(c) = mapTriangle.Compute(domainMesh.Cell2DVerticesCoordinates(c));

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
    problemData.Cell2Ds_DOF.resize(mesh.Cell2DTotalNumber(),
                                   std::vector<DiscreteProblemData::DOF*>(space.Order * 3, nullptr));

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

    for (unsigned int c = 0; c < mesh.Cell2DTotalNumber(); c++)
    {
      std::vector<DiscreteProblemData::DOF*>& dofs = problemData.Cell2Ds_DOF[c];
      for (unsigned int v = 0; v < 3; v++)
        dofs[v] = &problemData.Cell0Ds_DOF.at(mesh.Cell2DVertex(c, v));

      if (space.Order == 2)
      {
        for (unsigned int e = 0; e < 3; e++)
          dofs[3 + e] = &problemData.Cell1Ds_DOF.at(mesh.Cell2DEdge(c, e));
      }
    }

    FEM_RefElement_Langrange_PCC_Triangle_2D femRefElement;
    problemData.LocalSpace = femRefElement.Compute(space.Order);

    return problemData;
  }
  // ***************************************************************************
  std::list<Eigen::Triplet<double>> GeDiM4Py_Logic::AssembleStiffnessMatrix(K k,
                                                                            const Gedim::IMeshDAO& mesh,
                                                                            const std::vector<Gedim::MapTriangle::MapTriangleData>& cell2DsMap,
                                                                            const DiscreteProblemData& problemData)
  {
    std::list<Eigen::Triplet<double>> triplets;
    FEM_RefElement_Langrange_PCC_Triangle_2D femValues;
    Gedim::MapTriangle mapTriangle;
    const FEM_RefElement_Langrange_PCC_Triangle_2D::LocalSpace& localSpace = problemData.LocalSpace;
    PDE_Equation equation;

    const std::vector<Eigen::MatrixXd> referenceBasisFunctionDerivatives = femValues.Reference_BasisFunctionDerivatives(localSpace,
                                                                                                                        localSpace.ReferenceElement.InternalQuadrature.Points);
    const unsigned int numLocals = problemData.LocalSpace.NumberBasisFunctions;

    for (unsigned int cell2DIndex = 0; cell2DIndex < mesh.Cell2DTotalNumber(); cell2DIndex++)
    {
      // Compute cell geometric properties
      const Gedim::MapTriangle::MapTriangleData& cell2DMapData = cell2DsMap.at(cell2DIndex);

      const Eigen::MatrixXd cell2DQuadraturePoints = mapTriangle.F(cell2DMapData,
                                                                   localSpace.ReferenceElement.InternalQuadrature.Points);
      const Eigen::MatrixXd cell2DQuadratureWeights = localSpace.ReferenceElement.InternalQuadrature.Weights *
                                                      abs(cell2DMapData.DetB);

      const std::vector<Eigen::MatrixXd> basisFunctionDerivativeValues2D = femValues.BasisFunctionDerivatives(localSpace,
                                                                                                              cell2DMapData,
                                                                                                              referenceBasisFunctionDerivatives);

      const Eigen::VectorXd diffusioTermValues = k(cell2DQuadraturePoints);
      const Eigen::MatrixXd cellMatrixA = equation.ComputeStiffnessMatrix(numLocals,
                                                                          basisFunctionDerivativeValues2D,
                                                                          diffusioTermValues,
                                                                          cell2DQuadratureWeights);

      const std::vector<DiscreteProblemData::DOF*>& cell2D_DOF = problemData.Cell2Ds_DOF[cell2DIndex];

      for (unsigned int i = 0; i < numLocals; i++)
      {
        const DiscreteProblemData::DOF& dofI = *cell2D_DOF[i];
        if (dofI.Type != DiscreteProblemData::DOF::Types::DOF)
          continue;

        for (unsigned int j = 0; j < numLocals; j++)
        {
          const DiscreteProblemData::DOF& dofJ = *cell2D_DOF[j];

          if (dofI.Type != DiscreteProblemData::DOF::Types::DOF)
            continue;

          triplets.push_back(Eigen::Triplet<double>(dofI.Global_Index,
                                                    dofJ.Global_Index,
                                                    cellMatrixA(i, j)));
        }
      }
    }

    return triplets;
  }
  // ***************************************************************************
  Eigen::VectorXd GeDiM4Py_Logic::AssembleForcingTerm(F f,
                                                      const Gedim::IMeshDAO& mesh,
                                                      const std::vector<Gedim::MapTriangle::MapTriangleData>& cell2DsMap,
                                                      const DiscreteProblemData& problemData)
  {
    Eigen::VectorXd forcingTerm = Eigen::VectorXd::Zero(problemData.NumberDOFs);

    FEM_RefElement_Langrange_PCC_Triangle_2D femValues;
    Gedim::MapTriangle mapTriangle;
    const FEM_RefElement_Langrange_PCC_Triangle_2D::LocalSpace& localSpace = problemData.LocalSpace;
    PDE_Equation equation;

    const Eigen::MatrixXd referenceBasisFunctions = femValues.Reference_BasisFunctions(localSpace,
                                                                                       localSpace.ReferenceElement.InternalQuadrature.Points);

    const unsigned int numLocals = problemData.LocalSpace.NumberBasisFunctions;

    for (unsigned int cell2DIndex = 0; cell2DIndex < mesh.Cell2DTotalNumber(); cell2DIndex++)
    {
      // Compute cell geometric properties
      const Gedim::MapTriangle::MapTriangleData& cell2DMapData = cell2DsMap.at(cell2DIndex);

      const Eigen::MatrixXd cell2DQuadraturePoints = mapTriangle.F(cell2DMapData,
                                                                   localSpace.ReferenceElement.InternalQuadrature.Points);
      const Eigen::MatrixXd cell2DQuadratureWeights = localSpace.ReferenceElement.InternalQuadrature.Weights *
                                                      abs(cell2DMapData.DetB);

      const Eigen::MatrixXd basisFunctionValues2D = femValues.BasisFunctions(localSpace,
                                                                             cell2DMapData,
                                                                             referenceBasisFunctions);


      const Eigen::VectorXd forcingTermValues = f(cell2DQuadraturePoints);
      const Eigen::VectorXd cellForcingTerm = equation.ComputeCellForcingTerm(forcingTermValues,
                                                                              basisFunctionValues2D,
                                                                              cell2DQuadratureWeights);

      const std::vector<DiscreteProblemData::DOF*>& cell2D_DOF = problemData.Cell2Ds_DOF[cell2DIndex];

      for (unsigned int i = 0; i < numLocals; i++)
      {
        const DiscreteProblemData::DOF& dofI = *cell2D_DOF[i];
        if (dofI.Type != DiscreteProblemData::DOF::Types::DOF)
          continue;

        forcingTerm[dofI.Global_Index] = cellForcingTerm[i];
      }
    }

    return forcingTerm;
  }
  // ***************************************************************************
}
