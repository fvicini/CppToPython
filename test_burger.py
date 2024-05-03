import numpy as np
import GeDiM4Py as gedim

def Burger_a(numPoints, points):
	values = np.ones(numPoints)
	return values.ctypes.data

def Burger_b(numPoints, points):
	values = np.ones((2, numPoints))
	return values.ctypes.data

def Burger_c(numPoints, points):
	values = np.ones(numPoints)
	return values.ctypes.data

def Burger_non_linear_b(numPoints, points, u, u_x, u_y):
	vecu = gedim.make_nd_array(u, numPoints)
	values = vecu
	return values.ctypes.data

def Burger_non_linear_c(numPoints, points, u, u_x, u_y):
	vecu_x = gedim.make_nd_array(u_x, numPoints)
	vecu_y = gedim.make_nd_array(u_y, numPoints)
	values = vecu_x + vecu_y
	return values.ctypes.data

def Burger_f(numPoints, points):
	matPoints = gedim.make_nd_matrix(points, (3, numPoints), np.double)
	values = 32.0 * (matPoints[1,:] * (1.0 - matPoints[1,:]) + matPoints[0,:] * (1.0 - matPoints[0,:])) + \
	(16.0 * (1.0 - 2.0 * matPoints[0,:]) * matPoints[1,:] * (1.0 - matPoints[1,:]) + \
	16.0 * (1.0 - 2.0 * matPoints[1,:]) * matPoints[0,:] * (1.0 - matPoints[0,:])) * \
	16.0 * (matPoints[1,:] * (1.0 - matPoints[1,:]) * matPoints[0,:] * (1.0 - matPoints[0,:]))
	return values.ctypes.data

def Burger_non_linear_f(numPoints, points, u, u_x, u_y):
	vecu = gedim.make_nd_array(u, numPoints)
	vecu_x = gedim.make_nd_array(u_x, numPoints)
	vecu_y = gedim.make_nd_array(u_y, numPoints)
	values = vecu * (vecu_x + vecu_y)
	return values.ctypes.data

def Burger_non_linear_der_f(numPoints, points, u, u_x, u_y):
	vecu = gedim.make_nd_array(u, numPoints)
	vecu_x = gedim.make_nd_array(u_x, numPoints)
	vecu_y = gedim.make_nd_array(u_y, numPoints)
	values = vecu * (vecu_x + vecu_y)
	return values.ctypes.data

def Burger_exactSolution(numPoints, points):
	matPoints = gedim.make_nd_matrix(points, (3, numPoints), np.double)
	values = 16.0 * (matPoints[1,:] * (1.0 - matPoints[1,:]) * matPoints[0,:] * (1.0 - matPoints[0,:]))
	return values.ctypes.data

def Burger_exactDerivativeSolution(direction, numPoints, points):
	matPoints = gedim.make_nd_matrix(points, (3, numPoints), np.double)

	if direction == 0:
		values = 16.0 * (1.0 - 2.0 * matPoints[0,:]) * matPoints[1,:] * (1.0 - matPoints[1,:])
	elif direction == 1:
		values = 16.0 * (1.0 - 2.0 * matPoints[1,:]) * matPoints[0,:] * (1.0 - matPoints[0,:])
	else:
		values = np.zeros(numPoints)

	return values.ctypes.data

def Ones(numPoints, points):
	values = np.ones(numPoints)
	return values.ctypes.data

def OnesDerivative(numPoints, points):
	values = np.ones(numPoints)
	return values.ctypes.data

def Zeros(numPoints, points):
	values = np.zeros(numPoints)
	return values.ctypes.data

def ZerosDerivative(direction, numPoints, points):
	values = np.zeros(numPoints)
	return values.ctypes.data

if __name__ == '__main__':

	lib = gedim.ImportLibrary("./debug/GeDiM4Py.so")

	config = { 'GeometricTolerance': 1.0e-8 }
	gedim.Initialize(config, lib)
	
	meshSizes = [0.1, 0.01]
	order = 2

	for meshSize in meshSizes:
		print("Order", order, "meshSize", meshSize)

		domain = { 'SquareEdge': 1.0, 'VerticesBoundaryCondition': [1,1,1,1], 'EdgesBoundaryCondition': [1,2,1,3], 'DiscretizationType': 1, 'MeshCellsMaximumArea': meshSize }
		[meshInfo, mesh] = gedim.CreateDomainSquare(domain, lib)

		domain = { 'RectangleBase': 1.0, 'RectangleHeight': 1.0, 'VerticesBoundaryCondition': [1,1,1,1], 'EdgesBoundaryCondition': [1,2,1,3], 'DiscretizationType': 1, 'MeshCellsMaximumArea': meshSize }
		[meshInfo, mesh] = gedim.CreateDomainRectangle(domain, lib)

		gedim.PlotMesh(mesh)

		discreteSpace = { 'Order': order, 'Type': 1, 'BoundaryConditionsType': [1, 2, 3, 3] }
		[problemData, dofs, strongs] = gedim.Discretize(discreteSpace, lib)

		gedim.PlotDofs(mesh, dofs, strongs)
		
		residual_norm = 1.0
		solution_norm = 1.0;
		newton_tol = 1.0e-6
		max_iterations = 20
		num_iteration = 1
		
		u_k = np.zeros(problemData['NumberDOFs'])
		u_strong = np.zeros(problemData['NumberStrongs'])

		[stiffness, stiffnessStrong] = gedim.AssembleStiffnessMatrix(Burger_a, problemData, lib)
		[advection, advectionStrong] = gedim.AssembleNonLinearAdvectionMatrix(Burger_b, Burger_non_linear_b, u_k, u_strong, problemData, lib)
		[reaction, reactionStrong] = gedim.AssembleNonLinearReactionMatrix(Burger_c, Burger_non_linear_c, u_k, u_strong, problemData, lib)

		forcingTerm_g = gedim.AssembleForcingTerm(Burger_f, problemData, lib)
		forcingTerm_v = gedim.AssembleNonLinearForcingTerm(Burger_f, Burger_non_linear_f, u_k, u_strong, problemData, lib)
		forcingTerm_der_v = gedim.AssembleNonLinearDerivativeForcingTerm(Burger_f, Burger_non_linear_der_f, u_k, u_strong, problemData, lib)

		du = gedim.LUSolver(stiffness + advection + reaction, \
				forcingTerm_g - forcingTerm_v - forcingTerm_der_v, \
				lib)
		
		u_k = u_k + du
		
		du_normL2 = gedim.ComputeErrorL2(Zeros, du, np.zeros(problemData['NumberStrongs']), lib)
		u_errorL2 = gedim.ComputeErrorL2(Burger_exactSolution, u_k, u_strong, lib)
		u_errorH1 = gedim.ComputeErrorH1(Burger_exactDerivativeSolution, u_k, u_strong, lib)
		u_normL2 = gedim.ComputeErrorL2(Zeros, u_k, u_strong, lib)
		u_normH1 = gedim.ComputeErrorH1(ZerosDerivative, u_k, u_strong, lib)
		
		solution_norm = u_normL2;
        residual_norm = du_normL2;

		print("dofs", "h", "errorL2", "errorH1", "residual")
		print(problemData['NumberDOFs'], '{:.16e}'.format(problemData['H']), '{:.16e}'.format(errorL2 / u_normL2), '{:.16e}'.format(errorH1 / u_normH1), '{:.16e}'.format(residual_norm / u_normL2))

		gedim.ExportSolution(Burger_exactSolution, u_k, u_strong, lib)
		gedim.PlotSolution(mesh, dofs, strongs, solution, solutionStrong)
		
		num_iteration = num_iteration + 1
