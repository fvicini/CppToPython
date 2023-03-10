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
      const unsigned int vertexMarker = mesh.Cell0DMarker(p);
      const DiscreteSpace::BoundaryConditionTypes& type = space.BoundaryConditionsType.at(vertexMarker);
      DiscreteProblemData::DOF& cell0D_DOF =  problemData.Cell0Ds_DOF[p];

      switch (type)
      {
        case DiscreteSpace::BoundaryConditionTypes::None:
          cell0D_DOF.Type = DiscreteProblemData::DOF::Types::DOF;
          cell0D_DOF.Global_Index = problemData.NumberDOFs;
          cell0D_DOF.Boundary.Type = DiscreteProblemData::DOF::BoundaryInfo::Types::None;
          cell0D_DOF.Boundary.Marker = 0;
          problemData.NumberDOFs++;
          break;
        case DiscreteSpace::BoundaryConditionTypes::Weak:
          cell0D_DOF.Type = DiscreteProblemData::DOF::Types::DOF;
          cell0D_DOF.Global_Index = problemData.NumberDOFs;
          cell0D_DOF.Boundary.Type = DiscreteProblemData::DOF::BoundaryInfo::Types::Weak;
          cell0D_DOF.Boundary.Marker = vertexMarker;
          problemData.NumberDOFs++;
          break;
        case DiscreteSpace::BoundaryConditionTypes::Strong:
          cell0D_DOF.Type = DiscreteProblemData::DOF::Types::Strong;
          cell0D_DOF.Global_Index = problemData.NumberStrongs;
          cell0D_DOF.Boundary.Type = DiscreteProblemData::DOF::BoundaryInfo::Types::Strong;
          cell0D_DOF.Boundary.Marker = vertexMarker;
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
        const unsigned int edgeMarker = mesh.Cell1DMarker(e);
        const DiscreteSpace::BoundaryConditionTypes& type = space.BoundaryConditionsType.at(edgeMarker);
        DiscreteProblemData::DOF& cell1D_DOF =  problemData.Cell1Ds_DOF[e];

        switch (type)
        {
          case DiscreteSpace::BoundaryConditionTypes::None:
            cell1D_DOF.Type = DiscreteProblemData::DOF::Types::DOF;
            cell1D_DOF.Global_Index = problemData.NumberDOFs;
            cell1D_DOF.Boundary.Type = DiscreteProblemData::DOF::BoundaryInfo::Types::None;
            cell1D_DOF.Boundary.Marker = 0;
            problemData.NumberDOFs++;
            break;
          case DiscreteSpace::BoundaryConditionTypes::Weak:
            cell1D_DOF.Type = DiscreteProblemData::DOF::Types::DOF;
            cell1D_DOF.Global_Index = problemData.NumberDOFs;
            cell1D_DOF.Boundary.Type = DiscreteProblemData::DOF::BoundaryInfo::Types::Weak;
            cell1D_DOF.Boundary.Marker = edgeMarker;
            problemData.NumberDOFs++;
            break;
          case DiscreteSpace::BoundaryConditionTypes::Strong:
            cell1D_DOF.Type = DiscreteProblemData::DOF::Types::Strong;
            cell1D_DOF.Global_Index = problemData.NumberStrongs;
            cell1D_DOF.Boundary.Type = DiscreteProblemData::DOF::BoundaryInfo::Types::Strong;
            cell1D_DOF.Boundary.Marker = edgeMarker;
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

    problemData.DOFsCoordinate.resize(3, problemData.NumberDOFs);
    problemData.StrongsCoordinate.resize(3, problemData.NumberStrongs);

    for (unsigned int p = 0; p < mesh.Cell0DTotalNumber(); p++)
    {
      const DiscreteProblemData::DOF& cell0D_DOF = problemData.Cell0Ds_DOF[p];
      switch (cell0D_DOF.Type)
      {
        case DiscreteProblemData::DOF::Types::DOF:
          problemData.DOFsCoordinate.col(cell0D_DOF.Global_Index)<< mesh.Cell0DCoordinates(p);
          break;
        case DiscreteProblemData::DOF::Types::Strong:
          problemData.StrongsCoordinate.col(cell0D_DOF.Global_Index)<< mesh.Cell0DCoordinates(p);
          break;
        default:
          throw std::runtime_error("DOF Type " +
                                   std::to_string((unsigned int)cell0D_DOF.Type) +
                                   " not supported");
      }
    }

    if (space.Order == 2)
    {
      for (unsigned int e = 0; e < mesh.Cell1DTotalNumber(); e++)
      {
        const DiscreteProblemData::DOF& cell1D_DOF = problemData.Cell1Ds_DOF[e];
        switch (cell1D_DOF.Type)
        {
          case DiscreteProblemData::DOF::Types::DOF:
            problemData.DOFsCoordinate.col(cell1D_DOF.Global_Index)<< 0.5 * (mesh.Cell1DOriginCoordinates(e) + mesh.Cell1DEndCoordinates(e));;
            break;
          case DiscreteProblemData::DOF::Types::Strong:
            problemData.StrongsCoordinate.col(cell1D_DOF.Global_Index)<< 0.5 * (mesh.Cell1DOriginCoordinates(e) + mesh.Cell1DEndCoordinates(e));;
            break;
          default:
            throw std::runtime_error("DOF Type " +
                                     std::to_string((unsigned int)cell1D_DOF.Type) +
                                     " not supported");
        }
      }
    }

    FEM_RefElement_Langrange_PCC_Triangle_2D femRefElement;
    problemData.LocalSpace = femRefElement.Compute(space.Order);

    return problemData;
  }
  // ***************************************************************************
  void GeDiM4Py_Logic::AssembleStiffnessMatrix(K k,
                                               const Gedim::IMeshDAO& mesh,
                                               const std::vector<Gedim::MapTriangle::MapTriangleData>& cell2DsMap,
                                               const DiscreteProblemData& problemData,
                                               std::list<Eigen::Triplet<double>>& stiffnessTriplets,
                                               std::list<Eigen::Triplet<double>>& stiffnessStrongTriplets)
  {
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

      const double* diffusioTermValues = k(cell2DQuadraturePoints.cols(),
                                           cell2DQuadraturePoints.data());
      const Eigen::MatrixXd cellMatrixA = equation.ComputeStiffnessMatrix(numLocals,
                                                                          basisFunctionDerivativeValues2D,
                                                                          Eigen::Map<const Eigen::VectorXd>(diffusioTermValues,
                                                                                                            cell2DQuadraturePoints.cols()),
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

          switch (dofJ.Type)
          {
            case DiscreteProblemData::DOF::Types::DOF:
              stiffnessTriplets.push_back(Eigen::Triplet<double>(dofI.Global_Index,
                                                                 dofJ.Global_Index,
                                                                 cellMatrixA(i, j)));
              break;
            case DiscreteProblemData::DOF::Types::Strong:
              stiffnessStrongTriplets.push_back(Eigen::Triplet<double>(dofI.Global_Index,
                                                                       dofJ.Global_Index,
                                                                       cellMatrixA(i, j)));
              break;
            default:
              throw std::runtime_error("DOF Type " +
                                       std::to_string((unsigned int)dofJ.Type) +
                                       " not supported");
          }

          if (dofJ.Type != DiscreteProblemData::DOF::Types::DOF)
            continue;
        }
      }
    }
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

      const double* forcingTermValues = f(cell2DQuadraturePoints.cols(),
                                          cell2DQuadraturePoints.data());
      const Eigen::VectorXd cellForcingTerm = equation.ComputeCellForcingTerm(Eigen::Map<const Eigen::VectorXd>(forcingTermValues,
                                                                                                                cell2DQuadraturePoints.cols()),
                                                                              basisFunctionValues2D,
                                                                              cell2DQuadratureWeights);

      const std::vector<DiscreteProblemData::DOF*>& cell2D_DOF = problemData.Cell2Ds_DOF[cell2DIndex];

      for (unsigned int i = 0; i < numLocals; i++)
      {
        const DiscreteProblemData::DOF& dofI = *cell2D_DOF[i];

        if (dofI.Type != DiscreteProblemData::DOF::Types::DOF)
          continue;

        forcingTerm[dofI.Global_Index] += cellForcingTerm[i];
      }
    }

    return forcingTerm;
  }
  // ***************************************************************************
  Eigen::VectorXd GeDiM4Py_Logic::AssembleStrongSolution(Strong g,
                                                         const Gedim::IMeshDAO& mesh,
                                                         const std::vector<Gedim::MapTriangle::MapTriangleData>& cell2DsMap,
                                                         const DiscreteProblemData& problemData)
  {
    return Eigen::Map<Eigen::VectorXd>(g(problemData.NumberStrongs,
                                         problemData.StrongsCoordinate.data()),
                                       problemData.NumberStrongs);
  }
  // ***************************************************************************
  Eigen::VectorXd GeDiM4Py_Logic::AssembleWeakTerm(Weak g,
                                                   const int marker,
                                                   const Gedim::IMeshDAO& mesh,
                                                   const Eigen::VectorXd& cell2DEdgeLengths,
                                                   const Eigen::MatrixXd& cell2DEdgeTangents,
                                                   const std::vector<Gedim::MapTriangle::MapTriangleData>& cell2DsMap,
                                                   const DiscreteProblemData& problemData)
  {
    FEM_RefElement_Langrange_PCC_Triangle_2D femValues;
    Gedim::MapTriangle mapTriangle;
    const FEM_RefElement_Langrange_PCC_Triangle_2D::LocalSpace& localSpace = problemData.LocalSpace;
    PDE_Equation equation;

    const Eigen::MatrixXd referenceBasisFunctions = femValues.Reference_BasisFunctions(localSpace,
                                                                                       localSpace.ReferenceElement.InternalQuadrature.Points);

    const unsigned int numLocals = problemData.LocalSpace.NumberBasisFunctions;

    for (unsigned int cell2DIndex = 0; cell2DIndex < mesh.Cell2DTotalNumber(); cell2DIndex++)
    {
      const std::vector<DiscreteProblemData::DOF*>& cell2D_DOF = problemData.Cell2Ds_DOF[cell2DIndex];

      bool hasWeakCondition = false;
      for(unsigned int i = 0; i < numLocals; i ++)
      {
        const DiscreteProblemData::DOF& dofI = *cell2D_DOF[i];

        if (dofI.Boundary.Type != DiscreteProblemData::DOF::BoundaryInfo::Types::Weak)
          continue;

        if (dofI.Boundary.Marker != marker)
          continue;

        hasWeakCondition = true;
      }

      if (!hasWeakCondition)
        continue;

      for(unsigned int ed = 0; ed < 3; ed ++)
      {
        const unsigned int cell1DIndex = mesh.Cell2DEdge(cell2DIndex,
                                                           ed);
        if (!dofManager.IsCellWeakBoundaryCondition(cell1DIndex, 1))
          continue;

        const unsigned int edgeMarker = dofManager.CellMarker(cell1DIndex, 1);

        // map edge internal quadrature points
        const Eigen::Vector3d& edgeStart = mesh.Cell1DOriginCoordinates(cell1DIndex);
        const Eigen::Vector3d& edgeTangent = cell2DEdgeTangents.col(ed);

        const unsigned int numEdgeWeakQuadraturePoints = weakReferenceSegmentPoints.cols();
        Eigen::MatrixXd weakQuadraturePoints(3, numEdgeWeakQuadraturePoints);
        for (unsigned int q = 0; q < numEdgeWeakQuadraturePoints; q++)
        {
          weakQuadraturePoints.col(q) = edgeStart +
                                        weakReferenceSegmentPoints(0, q) *
                                        edgeTangent;
        }
        const double absMapDeterminant = std::abs(cell2DEdgeLengths[ed]);
        const Eigen::MatrixXd weakQuadratureWeights = weakReferenceSegmentWeights *
                                               absMapDeterminant;

        const Eigen::MatrixXd weakBasisFunctionsValues = femValues.ComputeBasisFunctionValues(localSpace,
                                                                                       cell2DMapData,
                                                                                       weakQuadraturePoints);

        const Eigen::VectorXd neumannValues = weakBoundaryCondition.Evaluate(edgeMarker,
                                                                      weakQuadraturePoints);

        // compute values of Neumann condition
        const Eigen::VectorXd neumannContributions = weakBasisFunctionsValues.transpose() *
                                              weakQuadratureWeights.asDiagonal() *
                                              neumannValues;

        // add contributions relative to edge extrema.
        for (unsigned int p = 0; p < 2; p++)
        {
          const unsigned int vertexLocalIndex = (ed + p) % 3;
          const unsigned int vertexGlobalIndex = mesh.Cell1DVertex(cell1DIndex,
                                                                   p);

          if (dofManager.CellMarker(vertexGlobalIndex, 0) != edgeMarker)
            continue;

          const unsigned int numCell0DLocals = dofManager.NumberLocals(vertexGlobalIndex,
                                                                       0);
          const unsigned int cell0DStartingLocalIdex = localSpace.Dof0DsIndex[vertexLocalIndex];
          const unsigned int cell0DEndingLocalIdex = localSpace.Dof0DsIndex[vertexLocalIndex + 1];
          Gedim::Output::Assert((cell0DEndingLocalIdex - cell0DStartingLocalIdex) ==
                                numCell0DLocals);

          for (unsigned int l = 0; l < numCell0DLocals; l++)
          {
            if (!dofManager.IsWeakBoundaryCondition(vertexGlobalIndex, l, 0))
              continue;

            const int globalNeumann_i = dofManager.GlobalIndex(vertexGlobalIndex,
                                                               l,
                                                               0);
            rightHandSide.AddValue(globalNeumann_i,
                                   neumannContributions(cell0DStartingLocalIdex + l));
          }
        }

        const unsigned int numCell1DLocals = dofManager.NumberLocals(cell1DIndex,
                                                                     1);
        const unsigned int cell1DStartingLocalIdex = localSpace.Dof1DsIndex[ed];
        const unsigned int cell1DEndingLocalIdex = localSpace.Dof1DsIndex[ed + 1];
        Gedim::Output::Assert((cell1DEndingLocalIdex - cell1DStartingLocalIdex) ==
                              numCell1DLocals);

        // add contributions relative to edge internal dofs
        for (unsigned int l = 0; l < numCell1DLocals; l++)
        {
          const int globalNeumann_i = dofManager.GlobalIndex(cell1DIndex, l, 1);

          rightHandSide.AddValue(globalNeumann_i,
                                 neumannContributions(cell1DStartingLocalIdex + l));
        }
      }
    }
  }
  // ***************************************************************************
}
