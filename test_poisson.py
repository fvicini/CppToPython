import numpy as np
import GeDiM4Py as gedim

def Poisson_A():
	return 10.0
def Poisson_B():
	return 0.1
def Poisson_C():
	return 2.0

def Poisson_a(numPoints, points):
	values = np.ones(numPoints) * Poisson_A()
	return values.ctypes.data

def Poisson_b(numPoints, points):
	values = np.ones((2, numPoints)) * Poisson_B()
	return values.ctypes.data

def Poisson_c(numPoints, points):
	values = np.ones(numPoints) * Poisson_C()
	return values.ctypes.data

def Poisson_f(numPoints, points):
	matPoints = gedim.make_nd_matrix(points, (3, numPoints), np.double)
	values = Poisson_A() * 32.0 * (matPoints[1,:] * (1.0 - matPoints[1,:]) + matPoints[0,:] * (1.0 - matPoints[0,:])) + \
	Poisson_B() * 16.0 * (1.0 - 2.0 * matPoints[0,:]) * matPoints[1,:] * (1.0 - matPoints[1,:]) + \
	Poisson_B() * 16.0 * (1.0 - 2.0 * matPoints[1,:]) * matPoints[0,:] * (1.0 - matPoints[0,:]) + \
	Poisson_C() * 16.0 * (matPoints[1,:] * (1.0 - matPoints[1,:]) * matPoints[0,:] * (1.0 - matPoints[0,:])) + Poisson_C() * 1.1
	return values.ctypes.data

def Poisson_exactSolution(numPoints, points):
	matPoints = gedim.make_nd_matrix(points, (3, numPoints), np.double)
	values = 16.0 * (matPoints[1,:] * (1.0 - matPoints[1,:]) * matPoints[0,:] * (1.0 - matPoints[0,:])) + 1.1
	return values.ctypes.data

def Poisson_exactDerivativeSolution(direction, numPoints, points):
	matPoints = gedim.make_nd_matrix(points, (3, numPoints), np.double)

	if direction == 0:
		values = 16.0 * (1.0 - 2.0 * matPoints[0,:]) * matPoints[1,:] * (1.0 - matPoints[1,:])
	elif direction == 1:
		values = 16.0 * (1.0 - 2.0 * matPoints[1,:]) * matPoints[0,:] * (1.0 - matPoints[0,:])
	else:
		values = np.zeros(numPoints)

	return values.ctypes.data

def Poisson_strongTerm(numPoints, points):
	matPoints = gedim.make_nd_matrix(points, (3, numPoints), np.double)
	values = 16.0 * (matPoints[1,:] * (1.0 - matPoints[1,:]) * matPoints[0,:] * (1.0 - matPoints[0,:])) + 1.1
	return values.ctypes.data

def Poisson_weakTerm_right(numPoints, points):
	matPoints = gedim.make_nd_matrix(points, (3, numPoints), np.double)
	values = Poisson_A() * 16.0 * (1.0 - 2.0 * matPoints[0,:]) * matPoints[1,:] * (1.0 - matPoints[1,:])
	return values.ctypes.data
	
def Poisson_weakTerm_left(numPoints, points):
	matPoints = gedim.make_nd_matrix(points, (3, numPoints), np.double)
	values = - Poisson_A() * 16.0 * (1.0 - 2.0 * matPoints[0,:]) * matPoints[1,:] * (1.0 - matPoints[1,:])
	return values.ctypes.data

if __name__ == '__main__':

	lib = gedim.ImportLibrary("./release/GeDiM4Py.so")

	config = { 'GeometricTolerance': 1.0e-8 }
	gedim.Initialize(config, lib)
	
	meshSizes = [0.1 ]
	order = 2

	for meshSize in meshSizes:
		domain = { 'SquareEdge': 1.0, 'VerticesBoundaryCondition': [1,1,1,1], 'EdgesBoundaryCondition': [1,2,1,3], 'DiscretizationType': 1, 'MeshCellsMaximumArea': meshSize }
		[meshInfo, mesh] = gedim.CreateDomainSquare(domain, lib)

		domain = { 'RectangleBase': 1.0, 'RectangleHeight': 1.0, 'VerticesBoundaryCondition': [1,1,1,1], 'EdgesBoundaryCondition': [1,2,1,3], 'DiscretizationType': 1, 'MeshCellsMaximumArea': meshSize }
		[meshInfo, mesh] = gedim.CreateDomainRectangle(domain, lib)

		gedim.PlotMesh(mesh)

		discreteSpace = { 'Order': order, 'Type': 1, 'BoundaryConditionsType': [1, 2, 3, 3] }
		[problemData, dofs, strongs] = gedim.Discretize(discreteSpace, lib)

		gedim.PlotDofs(mesh, dofs, strongs)

		[stiffness, stiffnessStrong] = gedim.AssembleStiffnessMatrix(Poisson_a, problemData, lib)

		[advection, advectionStrong] = gedim.AssembleAdvectionMatrix(Poisson_b, problemData, lib)

		[reaction, reactionStrong] = gedim.AssembleReactionMatrix(Poisson_c, problemData, lib)

		forcingTerm = gedim.AssembleForcingTerm(Poisson_f, problemData, lib)

		solutionStrong = gedim.AssembleStrongSolution(Poisson_strongTerm, 1, problemData, lib)
		
		weakTerm_right = gedim.AssembleWeakTerm(Poisson_weakTerm_right, 2, problemData, lib)
		weakTerm_left = gedim.AssembleWeakTerm(Poisson_weakTerm_left, 3, problemData, lib)

		solution = gedim.LUSolver(stiffness + advection + reaction, \
				forcingTerm - \
				(stiffnessStrong + advectionStrong + reactionStrong) @ solutionStrong + \
				weakTerm_right + \
				weakTerm_left, lib)

		errorL2 = gedim.ComputeErrorL2(Poisson_exactSolution, solution, solutionStrong, lib)

		errorH1 = gedim.ComputeErrorH1(Poisson_exactDerivativeSolution, solution, solutionStrong, lib)

		print("dofs", "h", "errorL2", "errorH1")
		print(problemData['NumberDOFs'], '{:.16e}'.format(problemData['H']), '{:.16e}'.format(errorL2), '{:.16e}'.format(errorH1))

		gedim.PlotSolution(mesh, dofs, strongs, solution, solutionStrong)
