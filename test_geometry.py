import numpy as np
import GeDiM4Py as gedim

def OmegaTilde1(numPoints, points):
  matPoints = gedim.make_nd_matrix(points, (3, numPoints), np.double)
  values = np.zeros(numPoints)
  for p in range(0, numPoints):
    if (matPoints[0,p] <= (1.0 + 1.0e-8)):
      values[p] = 1.
  return values.ctypes.data

def OmegaTilde2_1(numPoints, points):
  matPoints = gedim.make_nd_matrix(points, (3, numPoints), np.double)
  values = np.zeros((3, numPoints), order='F')
  for p in range(0, numPoints):
    if (matPoints[0,p] > (1.0 + 1.0e-8)):
      values[0, p] = 1.
  return values.ctypes.data

def OmegaTilde2_2(numPoints, points):
  matPoints = gedim.make_nd_matrix(points, (3, numPoints), np.double)
  values = np.zeros((3, numPoints), order='F')
  for p in range(0, numPoints):
    if (matPoints[0,p] > (1.0 + 1.0e-8)):
      values[2, p] = 1.
  return values.ctypes.data

def forcing_term(points):
  return 32.0 * (points[1,:] * (1.0 - points[1,:]) + points[0,:] * (1.0 + MU_TILDE - points[0,:]))

def forcing_term_ref_11(points):
  return 32.0 * (points[1,:] * (1.0 - points[1,:]) + points[0,:] * (1.0 - points[0,:]))
def forcing_term_ref_12(points):
  return 32.0 * points[0,:]
def forcing_term_ref_21(points):
  return 32.0 * (points[1,:] * (1.0 - points[1,:]))
def forcing_term_ref_22(points):
  return 32.0 * (points[0,:] * (2.0 - points[0,:]))
def forcing_term_ref_23(points):
  return 32.0 * (2.0 - points[0,:])
 
def exact_solution(points):
  return 16.0 * (points[1,:] * (1.0 - points[1,:]) * points[0,:] * (1.0 + MU_TILDE - points[0,:]))

def exact_derivative_solution(direction, points):
  if direction == 0:
    values = 16.0 * (1.0 + MU_TILDE - 2.0 * points[0,:]) * points[1,:] * (1.0 - points[1,:])
  elif direction == 1:
    values = 16.0 * (1.0 - 2.0 * points[1,:]) * points[0,:] * (1.0 + MU_TILDE - points[0,:])
  else:
    values = np.zeros(points.shape[1])

  return values

def ForcingTerm(numPoints, points):
  matPoints = gedim.make_nd_matrix(points, (3, numPoints), np.double)
  values = forcing_term(matPoints)
  return values.ctypes.data

def ForcingTerm11(numPoints, points):
  matPoints = gedim.make_nd_matrix(points, (3, numPoints), np.double)
  values = forcing_term_ref_11(matPoints)
  for p in range(0, numPoints):
    if (matPoints[0, p] > (1.0 + 1.0e-8)):
      values[p] = 0. 
  return values.ctypes.data
def ForcingTerm12(numPoints, points):
  matPoints = gedim.make_nd_matrix(points, (3, numPoints), np.double)
  values = forcing_term_ref_12(matPoints)
  for p in range(0, numPoints):
    if (matPoints[0, p] > (1.0 + 1.0e-8)):
      values[p] = 0. 
  return values.ctypes.data

def ForcingTerm21(numPoints, points):
  matPoints = gedim.make_nd_matrix(points, (3, numPoints), np.double)
  values = forcing_term_ref_21(matPoints)
  for p in range(0, numPoints):
    if (matPoints[0, p] <= (1.0 + 1.0e-8)):
      values[p] = 0. 
  return values.ctypes.data
def ForcingTerm22(numPoints, points):
  matPoints = gedim.make_nd_matrix(points, (3, numPoints), np.double)
  values = forcing_term_ref_22(matPoints)
  for p in range(0, numPoints):
    if (matPoints[0, p] <= (1.0 + 1.0e-8)):
      values[p] = 0. 
  return values.ctypes.data
def ForcingTerm23(numPoints, points):
  matPoints = gedim.make_nd_matrix(points, (3, numPoints), np.double)
  values = forcing_term_ref_23(matPoints)
  for p in range(0, numPoints):
    if (matPoints[0, p] <= (1.0 + 1.0e-8)):
      values[p] = 0. 
  return values.ctypes.data

def Poisson_exactSolution(numPoints, points):
  matPoints = gedim.make_nd_matrix(points, (3, numPoints), np.double)
  values = exact_solution(matPoints)
  return values.ctypes.data

def Poisson_exactDerivativeSolution(direction, numPoints, points):
  matPoints = gedim.make_nd_matrix(points, (3, numPoints), np.double)
  values = exact_derivative_solution(direction, matPoints)
  return values.ctypes.data

def Map(points, mu):
  numPoints = points.shape[1]
  mappedPoints = np.copy(points)

  for p in range(1, numPoints):
    if (points[0, p] > 1.0 + 1.0e-8):
      mappedPoints[0, p] = mu * points[0, p] + (1. - mu)
  return mappedPoints

def thetaA(mu):
  return [1.0, 1.0 / mu[0], mu[0]]
def thetaF(mu):
  return [1.0, mu[0], mu[0], mu[0] * mu[0] * mu[0], mu[0] * mu[0] * (1.0 - mu[0])]

def normX(v, X):
  return np.sqrt(np.transpose(v) @ X @ v)

def ProjectSystem(AQH, fQH, B):
    AQN = []
    fQN = []
    for AH in AQH:
        AQN.append(np.copy(np.transpose(B) @ AH @ B))
    for fH in fQH:
        fQN.append(np.copy(np.transpose(B) @ fH))
    return [AQN, fQN]

def Solve_full_order(AQH, fQH, thetaA_mu, thetaF_mu):
    A = thetaA_mu[0] * AQH[0]
    f = thetaF_mu[0] * fQH[0]
    for i in range(1, len(AQH)):
        A += thetaA_mu[i] * AQH[i]
    for i in range(1, len(fQH)):
        f += thetaF_mu[i] * fQH[i]
    return gedim.LUSolver(A, f, lib)

def Solve_reduced_order(AQN, fQN, thetaA_mu, thetaF_mu):
    A = thetaA_mu[0] * AQN[0]
    f = thetaF_mu[0] * fQN[0]
    for i in range(1, len(AQN)):
        A += thetaA_mu[i] * AQN[i]
    for i in range(1, len(fQN)):
        f += thetaF_mu[i] * fQN[i]
    return np.linalg.solve(A, f)

def POD(AQH, fQH, X, N_max, tol):
    #### snapshot matrix creation
    snapshot_matrix = []

    for mu in training_set:
        snapshot = Solve_full_order(AQH, fQH, thetaA(mu), thetaF(mu))
        snapshot_matrix.append(np.copy(snapshot))
    
    snapshot_matrix = np.array(snapshot_matrix) 

    ### covariance matrix
    C = snapshot_matrix @ X @ np.transpose(snapshot_matrix) ## metti inner product
    L_e, VM_e = np.linalg.eig(C)
    eigenvalues = []
    eigenvectors = []

    for i in range(len(L_e)):
        eig_real = L_e[i].real
        eig_complex = L_e[i].imag
        assert np.isclose(eig_complex, 0.)
        eigenvalues.append(eig_real)
        eigenvectors.append(VM_e[i].real)

    total_energy = sum(eigenvalues)
    retained_energy_vector = np.cumsum(eigenvalues)
    relative_retained_energy = retained_energy_vector/total_energy

    if all(flag==False for flag in relative_retained_energy >= (1.0 - tol)):
        N = N_max
    else:
        N = np.argmax(relative_retained_energy >= (1.0 - tol)) + 1

    # Create the basis function matrix
    basis_functions = []
    for n in range(N):
        eigenvector =  eigenvectors[n]
        # basis = (1/np.sqrt(M))*np.transpose(snapshot_matrix)@eigenvector 
        basis = np.transpose(snapshot_matrix) @ eigenvector
        norm = normX(basis, X)
        # norm = np.sqrt(np.transpose(basis)@basis)
        basis /= norm
        basis_functions.append(np.copy(basis))

    return [N, np.transpose(np.array(basis_functions))]

def TestSingleParameter(AQH, fQH, AQN, fQN, B, mu):
    reduced_solution = Solve_reduced_order(AQN, fQN, thetaA(mu), thetaF(mu))
    full_solution = Solve_full_order(AQH, fQH, thetaA(mu), thetaF(mu))

    ###### plot #######
    proj_reduced_solution = B @ reduced_solution

    ### computing error
    error_function = full_solution - proj_reduced_solution
    error_norm_squared_component = np.transpose(error_function) @ X @ error_function
    abs_err_ROM = np.sqrt(abs(error_norm_squared_component))

    full_solution_norm_squared_component = np.transpose(full_solution) @ X @ full_solution
    rel_err_ROM = abs_err_ROM / np.sqrt(abs(full_solution_norm_squared_component))

    solutionStrong = np.zeros(problemData['NumberStrongs'])

    map_mu = mu[0]
    mappedMesh = Map(mesh, map_mu)
    mappedDofs = Map(dofs, map_mu)
    mappedStrongs = Map(strongs, map_mu)

    exact_solution_Dofs = exact_solution(mappedDofs)
    exact_solution_norm_squared_component = np.transpose(exact_solution_Dofs) @ X @ exact_solution_Dofs

    if activatePlot:
      gedim.PlotSolution(mappedMesh, mappedDofs, mappedStrongs, proj_reduced_solution, solutionStrong, "MOR Solution")
      gedim.PlotSolution(mappedMesh, mappedDofs, mappedStrongs, full_solution, solutionStrong, "FOR Solution")
      gedim.PlotSolution(mappedMesh, mappedDofs, mappedStrongs, exact_solution_Dofs, solutionStrong, "EXT Solution")
    
    abs_err_L2 = gedim.ComputeErrorL2(Poisson_exactSolution, full_solution, solutionStrong, lib)
    abs_err_H1 = gedim.ComputeErrorH1(Poisson_exactDerivativeSolution, full_solution, solutionStrong, lib)
    rel_err_L2 = abs_err_L2 / np.sqrt(abs(exact_solution_norm_squared_component))
    rel_err_H1 = abs_err_H1 / np.sqrt(abs(exact_solution_norm_squared_component))
  
    return [rel_err_ROM, abs_err_ROM, rel_err_L2, abs_err_L2, rel_err_H1, abs_err_H1]

if __name__ == '__main__':

  lib = gedim.ImportLibrary("./release/GeDiM4Py.so")

  config = { 'GeometricTolerance': 1.0e-8 }
  gedim.Initialize(config, lib)
  
  ### Parameters ###
  activatePlot = True
  order = 1
  meshType = 1
  meshSizes = [0.001]
  mu_test = 2.5

  print("DOFs","N","rel_err_ROM","abs_err_ROM","rel_err_L2","abs_err_L2","rel_err_H1","abs_err_H1")

  for meshSize in meshSizes:
    if meshType == 0:
      domain = { 'RectangleBase': 2.0, 'RectangleHeight': 1.0, 'VerticesBoundaryCondition': [1,1,1,1], 'EdgesBoundaryCondition': [1,1,1,1], 'DiscretizationType': 1, 'MeshCellsMaximumArea': meshSize }
      [meshInfo, mesh] = gedim.CreateDomainRectangle(domain, lib)
    else:
      [meshInfo, mesh] = gedim.ImportDomainMesh2D(lib)

    if activatePlot:
      gedim.PlotMesh(mesh)

    discreteSpace = { 'Order': order, 'Type': 1, 'BoundaryConditionsType': [1, 2] }
    [problemData, dofs, strongs] = gedim.Discretize(discreteSpace, lib)

    if activatePlot:
      gedim.PlotDofs(mesh, dofs, strongs)

    [stiffness1, stiffnessStrong1] = gedim.AssembleStiffnessMatrix(OmegaTilde1, problemData, lib)
    [stiffness2_1, stiffnessStrong2_1] = gedim.AssembleAnisotropicStiffnessMatrix(OmegaTilde2_1, problemData, lib)
    [stiffness2_2, stiffnessStrong2_2] = gedim.AssembleAnisotropicStiffnessMatrix(OmegaTilde2_2, problemData, lib)

    forcingTerm11 = gedim.AssembleForcingTerm(ForcingTerm11, problemData, lib)
    forcingTerm12 = gedim.AssembleForcingTerm(ForcingTerm12, problemData, lib)
    forcingTerm21 = gedim.AssembleForcingTerm(ForcingTerm21, problemData, lib)
    forcingTerm22 = gedim.AssembleForcingTerm(ForcingTerm22, problemData, lib)
    forcingTerm23 = gedim.AssembleForcingTerm(ForcingTerm23, problemData, lib)

    X = stiffness1 + stiffness2_1 + stiffness2_2

    ### define the problem
    AQH = [stiffness1, stiffness2_1, stiffness2_2]
    fQH = [forcingTerm11, forcingTerm12, forcingTerm21, forcingTerm22, forcingTerm23]

    ### define the training set
    M = 100
    mu1_range = [1.0, 3.5]
    P = np.array([mu1_range])

    training_set = np.random.uniform(low=P[:, 0], high=P[:, 1], size=(M, P.shape[0]))

    tol = 1.0e-7
    N_max = 20

    ### Compute POD
    [N_POD, B_POD] = POD(AQH, fQH, X, N_max, tol)

    [AQN_POD, fQN_POD] = ProjectSystem(AQH, fQH, B_POD)

    MU_TILDE = mu_test

    [rel_err_ROM, abs_err_ROM, rel_err_L2, abs_err_L2, rel_err_H1, abs_err_H1] = TestSingleParameter(AQH, fQH, AQN_POD, fQN_POD, B_POD, [mu_test])
    
    print(problemData['NumberDOFs'],\
          N_POD,\
          '{:.4e}'.format(np.mean(rel_err_ROM)),\
          '{:.4e}'.format(np.mean(abs_err_ROM)),\
          '{:.4e}'.format(np.mean(rel_err_L2)),\
          '{:.4e}'.format(np.mean(abs_err_L2)),\
          '{:.4e}'.format(np.mean(rel_err_H1)),\
          '{:.4e}'.format(np.mean(abs_err_H1)))
  