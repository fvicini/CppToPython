#include "GeDiM4Py_Logic.hpp"
#include "MeshMatricesDAO.hpp"
#include "MeshUtilities.hpp"
#include "PDE_Equation.hpp"
#include "FileTextReader.hpp"

namespace GedimForPy
{
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
    Gedim::MeshMatricesDAO meshDAO(mesh.Mesh);

    switch (domain.DiscretizationType)
    {
      case Domain2D::DiscretizationTypes::Triangular:
      {
        gedimData.MeshUtilities().CreateTriangularMesh(domain.Vertices,
                                                       domain.MeshCellsMaximumArea,
                                                       meshDAO);
        gedimData.MeshUtilities().ChangePolygonMeshMarkers(domain.Vertices,
                                                           domain.VerticesBoundaryCondition,
                                                           domain.EdgesBoundaryCondition,
                                                           meshDAO);
      }
        break;
      default:
        throw std::runtime_error("MeshGenerator " +
                                 std::to_string((unsigned int)domain.DiscretizationType) +
                                 " not supported");
    }

    Gedim::MapTriangle mapTriangle;
    mesh.Cell2DsMap.resize(meshDAO.Cell2DTotalNumber());
    for (unsigned int c = 0; c < meshDAO.Cell2DTotalNumber(); c++)
      mesh.Cell2DsMap.at(c) = mapTriangle.Compute(meshDAO.Cell2DVerticesCoordinates(c));

    mesh.MeshGeometricData = gedimData.MeshUtilities().FillMesh2DGeometricData(gedimData.GeometryUtilities(),
                                                                               meshDAO);

    return mesh;
  }
  // ***************************************************************************
  Domain2DMesh GeDiM4Py_Logic::ImportDomainMesh2D(const ImportMesh2D& domain,
                                                  InterfaceDataDAO& gedimData,
                                                  const bool& checkMesh)
  {
    Domain2DMesh mesh;
    Gedim::MeshMatricesDAO meshDAO(mesh.Mesh);

    Eigen::MatrixXd cell0Ds;
    vector<Eigen::VectorXi> originalCell2Ds;
    vector<Eigen::VectorXi> cell2Ds;

    {
      std::vector<std::string> cell0DsLines;
      Gedim::FileReader csvFileReader(domain.InputFolder + "/Cell0Ds.csv");

      if (!csvFileReader.Open())
        throw runtime_error("File not found in folder " + domain.InputFolder);

      csvFileReader.GetAllLines(cell0DsLines);
      csvFileReader.Close();

      unsigned int numCell0Ds = cell0DsLines.size() - 1;
      if (numCell0Ds == 0)
        throw runtime_error("File cell0Ds empty");


      cell0Ds.setZero(3, numCell0Ds);

      char temp;
      unsigned int id;
      for (unsigned int v = 0; v < numCell0Ds; v++)
      {
        istringstream converter(cell0DsLines[v + 1]);

        converter >> id;
        if (domain.Separator != ' ')
          converter >> temp;
        converter >> cell0Ds(0, v);
        if (domain.Separator != ' ')
          converter >> temp;
        converter >> cell0Ds(1, v);
      }
    }

    {
      std::vector<std::string> cell2DsLines;
      Gedim::FileReader csvFileReader(domain.InputFolder + "/Cell2Ds.csv");

      if (!csvFileReader.Open())
        throw runtime_error("File not found in folder " + domain.InputFolder);

      csvFileReader.GetAllLines(cell2DsLines);
      csvFileReader.Close();

      unsigned int numCell2Ds = cell2DsLines.size() - 1;
      if (numCell2Ds == 0)
        throw runtime_error("File cell2Ds empty");

      cell2Ds.resize(numCell2Ds, Eigen::VectorXi::Zero(3));
      originalCell2Ds.resize(numCell2Ds, Eigen::VectorXi::Zero(3));

      char temp;
      unsigned int id;
      for (unsigned int t = 0; t < numCell2Ds; t++)
      {
        istringstream converter(cell2DsLines[t + 1]);

        converter >> id;
        if (domain.Separator != ' ')
          converter >> temp;
        converter >> originalCell2Ds[t](0);
        if (domain.Separator != ' ')
          converter >> temp;
        converter >> originalCell2Ds[t](1);
        if (domain.Separator != ' ')
          converter >> temp;
        converter >> originalCell2Ds[t](2);

        Eigen::Matrix3d points;
        points.col(0)<< cell0Ds.col(originalCell2Ds.at(t)(0));
        points.col(1)<< cell0Ds.col(originalCell2Ds.at(t)(1));
        points.col(2)<< cell0Ds.col(originalCell2Ds.at(t)(2));

        std::vector<unsigned int> convexPoints = gedimData.GeometryUtilities().ConvexHull(points);
        cell2Ds.at(t)(0) = originalCell2Ds.at(t)(convexPoints.at(0));
        cell2Ds.at(t)(1) = originalCell2Ds.at(t)(convexPoints.at(1));
        cell2Ds.at(t)(2) = originalCell2Ds.at(t)(convexPoints.at(2));
      }
    }

    Gedim::MeshUtilities::ComputeMesh2DCell1DsResult cell1Ds = gedimData.MeshUtilities().ComputeMesh2DCell1Ds(cell0Ds,
                                                                                                              cell2Ds);

    gedimData.MeshUtilities().FillMesh2D(cell0Ds,
                                         cell1Ds.Cell1Ds,
                                         cell1Ds.Cell2Ds,
                                         meshDAO);

    {
      std::vector<std::string> cell2DsLines;
      Gedim::FileReader csvFileReader(domain.InputFolder + "/Cell2DsMarker.csv");

      if (!csvFileReader.Open())
        throw runtime_error("File not found in folder " + domain.InputFolder);

      csvFileReader.GetAllLines(cell2DsLines);
      csvFileReader.Close();

      unsigned int numCell2Ds = cell2DsLines.size() - 1;
      if (numCell2Ds == 0)
        throw runtime_error("File cell2DsMarker empty");

      char temp;
      unsigned int cell2DIndex, vertexIndex, marker;
      for (unsigned int t = 0; t < numCell2Ds; t++)
      {
        istringstream converter(cell2DsLines[t + 1]);

        converter >> cell2DIndex;
        if (domain.Separator != ' ')
          converter >> temp;
        converter >> vertexIndex;
        if (domain.Separator != ' ')
          converter >> temp;
        converter >> marker;


        if (marker != 0)
        {
          const unsigned int correctedOriginIndex = originalCell2Ds.at(cell2DIndex)[(vertexIndex + 1) % 3];
          const unsigned int correctedEndIndex = originalCell2Ds.at(cell2DIndex)[(vertexIndex + 2) % 3];
          const unsigned int cell1DIndex = meshDAO.Cell1DExists(correctedOriginIndex,
                                                                correctedEndIndex) ?
                                             meshDAO.Cell1DByExtremes(correctedOriginIndex, correctedEndIndex) :
                                             meshDAO.Cell1DByExtremes(correctedEndIndex, correctedOriginIndex);

          meshDAO.Cell0DSetMarker(correctedOriginIndex, marker);
          meshDAO.Cell0DSetMarker(correctedEndIndex, marker);
          meshDAO.Cell1DSetMarker(cell1DIndex, marker);
        }
      }
    }

    if (checkMesh)
    {
      Gedim::MeshUtilities::CheckMesh2DConfiguration check;
      check.Cell1D_CheckNeighbours = false;
      gedimData.MeshUtilities().CheckMesh2D(check,
                                            gedimData.GeometryUtilities(),
                                            meshDAO);
    }

    Gedim::MapTriangle mapTriangle;
    mesh.Cell2DsMap.resize(meshDAO.Cell2DTotalNumber());
    for (unsigned int c = 0; c < meshDAO.Cell2DTotalNumber(); c++)
      mesh.Cell2DsMap.at(c) = mapTriangle.Compute(meshDAO.Cell2DVerticesCoordinates(c));

    mesh.MeshGeometricData = gedimData.MeshUtilities().FillMesh2DGeometricData(gedimData.GeometryUtilities(),
                                                                               meshDAO);

    return mesh;
  }
  // ***************************************************************************
  DiscreteProblemData GeDiM4Py_Logic::Discretize(const Gedim::IMeshDAO& mesh,
                                                 const Gedim::MeshUtilities::MeshGeometricData2D& meshGeometricData,
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

    problemData.H = *max_element(std::begin(meshGeometricData.Cell2DsDiameters),
                                 std::end(meshGeometricData.Cell2DsDiameters));
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
  void GeDiM4Py_Logic::AssembleStiffnessMatrix(A a,
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
      const Gedim::MapTriangle::MapTriangleData& cell2DMapData = cell2DsMap.at(cell2DIndex);

      const Eigen::MatrixXd cell2DQuadraturePoints = mapTriangle.F(cell2DMapData,
                                                                   localSpace.ReferenceElement.InternalQuadrature.Points);
      const Eigen::VectorXd cell2DQuadratureWeights = localSpace.ReferenceElement.InternalQuadrature.Weights *
                                                      abs(cell2DMapData.DetB);

      const std::vector<Eigen::MatrixXd> basisFunctionDerivativeValues2D = femValues.BasisFunctionDerivatives(localSpace,
                                                                                                              cell2DMapData,
                                                                                                              referenceBasisFunctionDerivatives);

      const double* diffusioTermValues = a(cell2DQuadraturePoints.cols(),
                                           cell2DQuadraturePoints.data());

      Eigen::MatrixXd diffusionTerm = Eigen::MatrixXd::Zero(3, cell2DQuadraturePoints.cols());
      diffusionTerm.row(0)<< Eigen::Map<const Eigen::VectorXd>(diffusioTermValues,
                                                               cell2DQuadraturePoints.cols());
      diffusionTerm.row(2)<< Eigen::Map<const Eigen::VectorXd>(diffusioTermValues,
                                                               cell2DQuadraturePoints.cols());
      const Eigen::MatrixXd cellMatrixA = equation.ComputeStiffnessMatrix(numLocals,
                                                                          basisFunctionDerivativeValues2D,
                                                                          diffusionTerm,
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
  void GeDiM4Py_Logic::AssembleAnisotropicStiffnessMatrix(A a,
                                                          const Gedim::IMeshDAO& mesh,
                                                          const std::vector<Gedim::MapTriangle::MapTriangleData>& cell2DsMap,
                                                          const DiscreteProblemData& problemData,
                                                          std::list<Eigen::Triplet<double> >& stiffnessTriplets,
                                                          std::list<Eigen::Triplet<double> >& stiffnessStrongTriplets)
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
      const Gedim::MapTriangle::MapTriangleData& cell2DMapData = cell2DsMap.at(cell2DIndex);

      const Eigen::MatrixXd cell2DQuadraturePoints = mapTriangle.F(cell2DMapData,
                                                                   localSpace.ReferenceElement.InternalQuadrature.Points);
      const Eigen::VectorXd cell2DQuadratureWeights = localSpace.ReferenceElement.InternalQuadrature.Weights *
                                                      abs(cell2DMapData.DetB);

      const std::vector<Eigen::MatrixXd> basisFunctionDerivativeValues2D = femValues.BasisFunctionDerivatives(localSpace,
                                                                                                              cell2DMapData,
                                                                                                              referenceBasisFunctionDerivatives);

      const double* diffusioTermValues = a(cell2DQuadraturePoints.cols(),
                                           cell2DQuadraturePoints.data());
      const Eigen::MatrixXd cellMatrixA = equation.ComputeStiffnessMatrix(numLocals,
                                                                          basisFunctionDerivativeValues2D,
                                                                          Eigen::Map<const Eigen::MatrixXd>(diffusioTermValues,
                                                                                                            3,
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
  void GeDiM4Py_Logic::AssembleAdvectionMatrix(B b,
                                               const Gedim::IMeshDAO& mesh,
                                               const std::vector<Gedim::MapTriangle::MapTriangleData>& cell2DsMap,
                                               const DiscreteProblemData& trial_Functions,
                                               const DiscreteProblemData& test_Functions,
                                               std::list<Eigen::Triplet<double>>& advectionTriplets,
                                               std::list<Eigen::Triplet<double>>& advectionStrongTriplets)
  {
    FEM_RefElement_Langrange_PCC_Triangle_2D trial_femValues, test_femValues;
    Gedim::MapTriangle mapTriangle;
    const FEM_RefElement_Langrange_PCC_Triangle_2D::LocalSpace& trial_localSpace = trial_Functions.LocalSpace;
    const FEM_RefElement_Langrange_PCC_Triangle_2D::LocalSpace& test_localSpace = test_Functions.LocalSpace;
    PDE_Equation equation;

    const Eigen::MatrixXd& internal_quadraturePoints = trial_localSpace.ReferenceElement.InternalQuadrature.Points;
    const Eigen::VectorXd& internal_quadratureWeights = trial_localSpace.ReferenceElement.InternalQuadrature.Weights;

    const Eigen::MatrixXd test_referenceBasisFunctions = test_femValues.Reference_BasisFunctions(test_localSpace,
                                                                                                 internal_quadraturePoints);
    const std::vector<Eigen::MatrixXd> trial_referenceBasisFunctionDerivatives = trial_femValues.Reference_BasisFunctionDerivatives(trial_localSpace,
                                                                                                                                    internal_quadraturePoints);
    const unsigned int trial_numLocals = trial_Functions.LocalSpace.NumberBasisFunctions;
    const unsigned int test_numLocals = test_Functions.LocalSpace.NumberBasisFunctions;

    for (unsigned int cell2DIndex = 0; cell2DIndex < mesh.Cell2DTotalNumber(); cell2DIndex++)
    {
      const Gedim::MapTriangle::MapTriangleData& cell2DMapData = cell2DsMap.at(cell2DIndex);

      const Eigen::MatrixXd cell2DQuadraturePoints = mapTriangle.F(cell2DMapData,
                                                                   internal_quadraturePoints);
      const Eigen::VectorXd cell2DQuadratureWeights = internal_quadratureWeights *
                                                      abs(cell2DMapData.DetB);

      const Eigen::MatrixXd test_basisFunctionValues2D = test_femValues.BasisFunctions(test_localSpace,
                                                                                       cell2DMapData,
                                                                                       test_referenceBasisFunctions);

      const std::vector<Eigen::MatrixXd> trial_basisFunctionDerivativeValues2D = trial_femValues.BasisFunctionDerivatives(trial_localSpace,
                                                                                                                          cell2DMapData,
                                                                                                                          trial_referenceBasisFunctionDerivatives);

      const double* advectionTermValues = b(cell2DQuadraturePoints.cols(),
                                            cell2DQuadraturePoints.data());
      const Eigen::MatrixXd cellMatrixB = equation.ComputeCellAdvectionMatrix(trial_numLocals,
                                                                              test_numLocals,
                                                                              Eigen::Map<const Eigen::MatrixXd>(advectionTermValues,
                                                                                                                2,
                                                                                                                cell2DQuadraturePoints.cols()),
                                                                              test_basisFunctionValues2D,
                                                                              trial_basisFunctionDerivativeValues2D,
                                                                              cell2DQuadratureWeights);

      const std::vector<DiscreteProblemData::DOF*>& trial_cell2D_DOF = trial_Functions.Cell2Ds_DOF[cell2DIndex];
      const std::vector<DiscreteProblemData::DOF*>& test_cell2D_DOF = test_Functions.Cell2Ds_DOF[cell2DIndex];

      for (unsigned int i = 0; i < test_numLocals; i++)
      {
        const DiscreteProblemData::DOF& dofI = *test_cell2D_DOF[i];

        if (dofI.Type != DiscreteProblemData::DOF::Types::DOF)
          continue;

        for (unsigned int j = 0; j < trial_numLocals; j++)
        {
          const DiscreteProblemData::DOF& dofJ = *trial_cell2D_DOF[j];

          switch (dofJ.Type)
          {
            case DiscreteProblemData::DOF::Types::DOF:
              advectionTriplets.push_back(Eigen::Triplet<double>(dofI.Global_Index,
                                                                 dofJ.Global_Index,
                                                                 cellMatrixB(i, j)));
              break;
            case DiscreteProblemData::DOF::Types::Strong:
              advectionStrongTriplets.push_back(Eigen::Triplet<double>(dofI.Global_Index,
                                                                       dofJ.Global_Index,
                                                                       cellMatrixB(i, j)));
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
  void GeDiM4Py_Logic::AssembleReactionMatrix(C c,
                                              const Gedim::IMeshDAO& mesh,
                                              const std::vector<Gedim::MapTriangle::MapTriangleData>& cell2DsMap,
                                              const DiscreteProblemData& problemData,
                                              std::list<Eigen::Triplet<double> >& reactionTriplets,
                                              std::list<Eigen::Triplet<double> >& reactionStrongTriplets)
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
      const Gedim::MapTriangle::MapTriangleData& cell2DMapData = cell2DsMap.at(cell2DIndex);

      const Eigen::MatrixXd cell2DQuadraturePoints = mapTriangle.F(cell2DMapData,
                                                                   localSpace.ReferenceElement.InternalQuadrature.Points);
      const Eigen::VectorXd cell2DQuadratureWeights = localSpace.ReferenceElement.InternalQuadrature.Weights *
                                                      abs(cell2DMapData.DetB);

      const Eigen::MatrixXd basisFunctionValues2D = femValues.BasisFunctions(localSpace,
                                                                             cell2DMapData,
                                                                             referenceBasisFunctions);

      const double* reactionTermValues = c(cell2DQuadraturePoints.cols(),
                                           cell2DQuadraturePoints.data());
      const Eigen::MatrixXd cellMatrixC = equation.ComputeCellReactionMatrix(Eigen::Map<const Eigen::VectorXd>(reactionTermValues,
                                                                                                               cell2DQuadraturePoints.cols()),
                                                                             basisFunctionValues2D,
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
              reactionTriplets.push_back(Eigen::Triplet<double>(dofI.Global_Index,
                                                                dofJ.Global_Index,
                                                                cellMatrixC(i, j)));
              break;
            case DiscreteProblemData::DOF::Types::Strong:
              reactionStrongTriplets.push_back(Eigen::Triplet<double>(dofI.Global_Index,
                                                                      dofJ.Global_Index,
                                                                      cellMatrixC(i, j)));
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
      const Gedim::MapTriangle::MapTriangleData& cell2DMapData = cell2DsMap.at(cell2DIndex);

      const Eigen::MatrixXd cell2DQuadraturePoints = mapTriangle.F(cell2DMapData,
                                                                   localSpace.ReferenceElement.InternalQuadrature.Points);
      const Eigen::VectorXd cell2DQuadratureWeights = localSpace.ReferenceElement.InternalQuadrature.Weights *
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
                                                         const unsigned int& marker,
                                                         const Gedim::IMeshDAO& mesh,
                                                         const std::vector<Gedim::MapTriangle::MapTriangleData>& cell2DsMap,
                                                         const DiscreteProblemData& problemData)
  {
    Eigen::VectorXd strongSolution = Eigen::VectorXd::Zero(problemData.NumberStrongs);

    for (unsigned int v = 0; v < mesh.Cell0DTotalNumber(); v++)
    {
      const DiscreteProblemData::DOF& dof = problemData.Cell0Ds_DOF[v];

      if (dof.Boundary.Type != DiscreteProblemData::DOF::BoundaryInfo::Types::Strong)
        continue;

      if (dof.Boundary.Marker != marker)
        continue;

      const Eigen::Vector3d coordinate = problemData.StrongsCoordinate.col(dof.Global_Index);

      strongSolution[dof.Global_Index] = *g(1,
                                            coordinate.data());
    }

    if (problemData.LocalSpace.Order == 2)
    {
      for (unsigned int e = 0; e < mesh.Cell1DTotalNumber(); e++)
      {
        const DiscreteProblemData::DOF& dof = problemData.Cell1Ds_DOF[e];

        if (dof.Boundary.Type != DiscreteProblemData::DOF::BoundaryInfo::Types::Strong)
          continue;

        if (dof.Boundary.Marker != marker)
          continue;

        const Eigen::Vector3d coordinate = problemData.StrongsCoordinate.col(dof.Global_Index);

        strongSolution[dof.Global_Index] = *g(1,
                                              coordinate.data());
      }
    }

    return strongSolution;
  }
  // ***************************************************************************
  Eigen::VectorXd GeDiM4Py_Logic::AssembleWeakTerm(Weak g,
                                                   const unsigned int& marker,
                                                   const Gedim::IMeshDAO& mesh,
                                                   const std::vector<Eigen::MatrixXd>& cell2DsVertices,
                                                   const std::vector<Eigen::VectorXd>& cell2DsEdgeLengths,
                                                   const std::vector<Eigen::MatrixXd>& cell2DsEdgeTangents,
                                                   const std::vector<Gedim::MapTriangle::MapTriangleData>& cell2DsMap,
                                                   const DiscreteProblemData& problemData)
  {
    Eigen::VectorXd weakTerm = Eigen::VectorXd::Zero(problemData.NumberDOFs);

    FEM_RefElement_Langrange_PCC_Triangle_2D femValues;
    Gedim::MapTriangle mapTriangle;
    const FEM_RefElement_Langrange_PCC_Triangle_2D::LocalSpace& localSpace = problemData.LocalSpace;
    PDE_Equation equation;

    for (unsigned int cell2DIndex = 0; cell2DIndex < mesh.Cell2DTotalNumber(); cell2DIndex++)
    {
      const Gedim::MapTriangle::MapTriangleData& cell2DMapData = cell2DsMap.at(cell2DIndex);
      const Eigen::MatrixXd& cell2DVertices = cell2DsVertices[cell2DIndex];
      const Eigen::VectorXd& cell2DEdgeLengths = cell2DsEdgeLengths[cell2DIndex];
      const Eigen::MatrixXd& cell2DEdgeTangents = cell2DsEdgeTangents[cell2DIndex];

      for(unsigned int ed = 0; ed < 3; ed++)
      {
        const unsigned int cell1DIndex = mesh.Cell2DEdge(cell2DIndex,
                                                         ed);
        const unsigned int cell0DOriginIndex = mesh.Cell1DOrigin(cell1DIndex);
        const unsigned int cell0DEndIndex = mesh.Cell1DEnd(cell1DIndex);

        if (mesh.Cell1DMarker(cell1DIndex) != marker)
          continue;

        const DiscreteProblemData::DOF& dofOrigin = problemData.Cell0Ds_DOF[cell0DOriginIndex];
        const DiscreteProblemData::DOF& dofEnd = problemData.Cell0Ds_DOF[cell0DEndIndex];
        const DiscreteProblemData::DOF& dofEdge = problemData.Cell1Ds_DOF[cell1DIndex];

        const bool isDofOriginWeak = (dofOrigin.Boundary.Type == DiscreteProblemData::DOF::BoundaryInfo::Types::Weak &&
                                      dofOrigin.Boundary.Marker == marker);
        const bool isDofEndWeak = (dofEnd.Boundary.Type == DiscreteProblemData::DOF::BoundaryInfo::Types::Weak &&
                                   dofEnd.Boundary.Marker == marker);
        const bool isDofEdgeWeak = (dofEdge.Boundary.Type == DiscreteProblemData::DOF::BoundaryInfo::Types::Weak &&
                                    dofEdge.Boundary.Marker == marker);

        if (!isDofOriginWeak && !isDofEndWeak && !isDofEdgeWeak)
          continue;

        const Eigen::Vector3d& edgeStart = cell2DVertices.col(ed);
        const Eigen::Vector3d& edgeTangent = cell2DEdgeTangents.col(ed);

        const unsigned int numEdgeQuadraturePoints = localSpace.ReferenceElement.BorderQuadrature.Points.cols();
        Eigen::MatrixXd weakQuadraturePoints(3, numEdgeQuadraturePoints);
        for (unsigned int q = 0; q < numEdgeQuadraturePoints; q++)
        {
          weakQuadraturePoints.col(q) = edgeStart +
                                        localSpace.ReferenceElement.BorderQuadrature.Points(0, q) *
                                        edgeTangent;
        }

        const Eigen::MatrixXd weakQuadratureWeights = localSpace.ReferenceElement.BorderQuadrature.Weights *
                                                      std::abs(cell2DEdgeLengths[ed]);

        const Eigen::MatrixXd weakBasisFunctionsValues = femValues.BasisFunctionsOnPoints(localSpace,
                                                                                          cell2DMapData,
                                                                                          weakQuadraturePoints);

        const double* weakTermValues = g(numEdgeQuadraturePoints,
                                         weakQuadraturePoints.data());
        const Eigen::VectorXd neumannContributions = equation.ComputeCellWeakTerm(Eigen::Map<const Eigen::VectorXd>(weakTermValues,
                                                                                                                    numEdgeQuadraturePoints),
                                                                                  weakBasisFunctionsValues,
                                                                                  weakQuadratureWeights);

        if (isDofOriginWeak)
          weakTerm[dofOrigin.Global_Index] += neumannContributions[localSpace.Dof0DsIndex[ed]];
        if (isDofEndWeak)
          weakTerm[dofEnd.Global_Index] += neumannContributions[localSpace.Dof0DsIndex[(ed + 1) % 3]];
        if (isDofEdgeWeak)
          weakTerm[dofEdge.Global_Index] += neumannContributions[localSpace.Dof1DsIndex[ed]];
      }
    }

    return weakTerm;
  }
  // ***************************************************************************
  // ***************************************************************************
  Eigen::VectorXd GeDiM4Py_Logic::ComputeErrorL2(Exact u,
                                                 const Eigen::VectorXd& numeric,
                                                 const Eigen::VectorXd& strong,
                                                 const Gedim::IMeshDAO& mesh,
                                                 const std::vector<Gedim::MapTriangle::MapTriangleData>& cell2DsMap,
                                                 const DiscreteProblemData& problemData)
  {
    Eigen::VectorXd errorL2 = Eigen::VectorXd::Zero(mesh.Cell2DTotalNumber());

    FEM_RefElement_Langrange_PCC_Triangle_2D femValues;
    Gedim::MapTriangle mapTriangle;
    const FEM_RefElement_Langrange_PCC_Triangle_2D::LocalSpace& localSpace = problemData.LocalSpace;
    PDE_Equation equation;

    const Eigen::MatrixXd referenceBasisFunctions = femValues.Reference_BasisFunctions(localSpace,
                                                                                       localSpace.ReferenceElement.InternalQuadrature.Points);

    const unsigned int numLocals = problemData.LocalSpace.NumberBasisFunctions;

    for (unsigned int cell2DIndex = 0; cell2DIndex < mesh.Cell2DTotalNumber(); cell2DIndex++)
    {
      const Gedim::MapTriangle::MapTriangleData& cell2DMapData = cell2DsMap.at(cell2DIndex);

      const Eigen::MatrixXd cell2DQuadraturePoints = mapTriangle.F(cell2DMapData,
                                                                   localSpace.ReferenceElement.InternalQuadrature.Points);
      const Eigen::VectorXd cell2DQuadratureWeights = localSpace.ReferenceElement.InternalQuadrature.Weights *
                                                      abs(cell2DMapData.DetB);

      const Eigen::MatrixXd basisFunctionValues2D = femValues.BasisFunctions(localSpace,
                                                                             cell2DMapData,
                                                                             referenceBasisFunctions);

      const std::vector<DiscreteProblemData::DOF*>& cell2D_DOF = problemData.Cell2Ds_DOF[cell2DIndex];

      Eigen::VectorXd localNumericSolution = Eigen::VectorXd::Zero(numLocals);
      for(unsigned int i = 0; i < numLocals; ++i)
      {
        const DiscreteProblemData::DOF& dofI = *cell2D_DOF[i];

        switch (dofI.Type)
        {
          case DiscreteProblemData::DOF::Types::DOF:
            localNumericSolution[i] = numeric[dofI.Global_Index];
            break;
          case DiscreteProblemData::DOF::Types::Strong:
            localNumericSolution[i] = strong[dofI.Global_Index];
            break;
          default:
            throw std::runtime_error("DOF Type " +
                                     std::to_string((unsigned int)dofI.Type) +
                                     " not supported");
        }
      }

      const double* exactSolutionValues = u(cell2DQuadraturePoints.cols(),
                                            cell2DQuadraturePoints.data());

      const Eigen::VectorXd localError = (basisFunctionValues2D * localNumericSolution -
                                          Eigen::Map<const Eigen::VectorXd>(exactSolutionValues,
                                                                            cell2DQuadraturePoints.cols())).array().square();

      errorL2[cell2DIndex] = cell2DQuadratureWeights.transpose() * localError;
    }

    return errorL2;
  }
  // ***************************************************************************
  Eigen::VectorXd GeDiM4Py_Logic::ComputeErrorH1(ExactDerivative uDer,
                                                 const Eigen::VectorXd& numeric,
                                                 const Eigen::VectorXd& strong,
                                                 const Gedim::IMeshDAO& mesh,
                                                 const std::vector<Gedim::MapTriangle::MapTriangleData>& cell2DsMap,
                                                 const DiscreteProblemData& problemData)
  {
    Eigen::VectorXd errorH1 = Eigen::VectorXd::Zero(mesh.Cell2DTotalNumber());

    FEM_RefElement_Langrange_PCC_Triangle_2D femValues;
    Gedim::MapTriangle mapTriangle;
    const FEM_RefElement_Langrange_PCC_Triangle_2D::LocalSpace& localSpace = problemData.LocalSpace;
    PDE_Equation equation;

    const std::vector<Eigen::MatrixXd> referenceBasisFunctionDerivatives = femValues.Reference_BasisFunctionDerivatives(localSpace,
                                                                                                                        localSpace.ReferenceElement.InternalQuadrature.Points);
    const unsigned int numLocals = problemData.LocalSpace.NumberBasisFunctions;

    for (unsigned int cell2DIndex = 0; cell2DIndex < mesh.Cell2DTotalNumber(); cell2DIndex++)
    {
      const Gedim::MapTriangle::MapTriangleData& cell2DMapData = cell2DsMap.at(cell2DIndex);

      const Eigen::MatrixXd cell2DQuadraturePoints = mapTriangle.F(cell2DMapData,
                                                                   localSpace.ReferenceElement.InternalQuadrature.Points);
      const Eigen::VectorXd cell2DQuadratureWeights = localSpace.ReferenceElement.InternalQuadrature.Weights *
                                                      abs(cell2DMapData.DetB);

      const std::vector<Eigen::MatrixXd> basisFunctionDerivativeValues2D = femValues.BasisFunctionDerivatives(localSpace,
                                                                                                              cell2DMapData,
                                                                                                              referenceBasisFunctionDerivatives);

      const std::vector<DiscreteProblemData::DOF*>& cell2D_DOF = problemData.Cell2Ds_DOF[cell2DIndex];

      Eigen::VectorXd localNumericSolution = Eigen::VectorXd::Zero(numLocals);
      for(unsigned int i = 0; i < numLocals; ++i)
      {
        const DiscreteProblemData::DOF& dofI = *cell2D_DOF[i];

        switch (dofI.Type)
        {
          case DiscreteProblemData::DOF::Types::DOF:
            localNumericSolution[i] = numeric[dofI.Global_Index];
            break;
          case DiscreteProblemData::DOF::Types::Strong:
            localNumericSolution[i] = strong[dofI.Global_Index];
            break;
          default:
            throw std::runtime_error("DOF Type " +
                                     std::to_string((unsigned int)dofI.Type) +
                                     " not supported");
        }
      }

      Eigen::VectorXd localError = Eigen::VectorXd::Zero(cell2DQuadraturePoints.cols());
      for (unsigned int dim = 0; dim < 2; dim++)
      {
        double* exactSolutionDerivativeValues = uDer(dim,
                                                     cell2DQuadraturePoints.cols(),
                                                     cell2DQuadraturePoints.data());

        localError.array() += (basisFunctionDerivativeValues2D[dim] * localNumericSolution -
                               Eigen::Map<const Eigen::VectorXd>(exactSolutionDerivativeValues,
                                                                 cell2DQuadraturePoints.cols())).array().square();
      }

      errorH1[cell2DIndex] = cell2DQuadratureWeights.transpose() * localError;
    }

    return errorH1;
  }
  // ***************************************************************************
}
