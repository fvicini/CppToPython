import numpy as np
import GeDiM4Py as gedim

def Stokes_V():
	return 1.0

def Stokes_v(numPoints, points):
	values = np.ones(numPoints) * Stokes_V()
	return values.ctypes.data

def Stokes_advection_1(numPoints, points):
	values = np.zeros((2, numPoints), order='F')
	values[0,:] = 1.0
	return values.ctypes.data

def Stokes_advection_2(numPoints, points):
	values = np.zeros((2, numPoints), order='F')
	values[1,:] = 1.0
	return values.ctypes.data

def Stokes_f_1(numPoints, points):
	matPoints = gedim.make_nd_matrix(points, (3, numPoints), np.double)
	values = - (+8.0 * np.pi * np.pi * np.cos(4.0 * np.pi * matPoints[0,:]) - 4.0 * np.pi * np.pi) * np.sin(2.0 * np.pi * matPoints[1,:]) * np.cos(2.0 * np.pi * matPoints[1,:]) + (+2.0 * np.pi * np.cos(2.0 * np.pi * matPoints[0,:]) * np.cos(2.0 * np.pi * matPoints[1,:]))
	return values.ctypes.data

def Stokes_f_2(numPoints, points):
	matPoints = gedim.make_nd_matrix(points, (3, numPoints), np.double)
	values = - (-8.0 * np.pi * np.pi * np.cos(4.0 * np.pi * matPoints[1,:]) + 4.0 * np.pi * np.pi) * np.sin(2.0 * np.pi * matPoints[0,:]) * np.cos(2.0 * np.pi * matPoints[0,:]) + (-2.0 * np.pi * np.sin(2.0 * np.pi * matPoints[0,:]) * np.sin(2.0 * np.pi * matPoints[1,:]))
	return values.ctypes.data

def Stokes_pressure_exactSolution(numPoints, points):
	matPoints = gedim.make_nd_matrix(points, (3, numPoints), np.double)
	values = np.sin(2.0 * np.pi * matPoints[0,:]) * np.cos(2.0 * np.pi * matPoints[1,:])
	return values.ctypes.data

def Stokes_speed_exactSolution_1(numPoints, points):
	matPoints = gedim.make_nd_matrix(points, (3, numPoints), np.double)
	values = +0.5 * np.sin(2.0 * np.pi * matPoints[0,:]) * np.sin(2.0 * np.pi * matPoints[0,:]) * np.sin(2.0 * np.pi * matPoints[1,:]) * np.cos(2.0 * np.pi * matPoints[1,:])
	return values.ctypes.data

def Stokes_speed_exactSolution_2(numPoints, points):
	matPoints = gedim.make_nd_matrix(points, (3, numPoints), np.double)
	values = -0.5 * np.sin(2.0 * np.pi * matPoints[1,:]) * np.sin(2.0 * np.pi * matPoints[1,:]) * np.sin(2.0 * np.pi * matPoints[0,:]) * np.cos(2.0 * np.pi * matPoints[0,:])
	return values.ctypes.data

if __name__ == '__main__':

	lib = gedim.ImportLibrary("./release/GeDiM4Py.so")

	config = { 'GeometricTolerance': 1.0e-8 }
	gedim.Initialize(config, lib)
	
	activatePlot = True
	meshSizes = [0.001]

	for meshSize in meshSizes:
		domain = { 'SquareEdge': 1.0, 'VerticesBoundaryCondition': [1,1,1,1], 'EdgesBoundaryCondition': [2,3,4,5], 'DiscretizationType': 1, 'MeshCellsMaximumArea': meshSize }
		[meshInfo, mesh] = gedim.CreateDomainSquare(domain, lib)

		if activatePlot:
			gedim.PlotMesh(mesh)

		pressure_discreteSpace = { 'Order': 1, 'Type': 1, 'BoundaryConditionsType': [1, 2, 1, 1, 1, 1] }
		speed_discreteSpace = { 'Order': 2, 'Type': 1, 'BoundaryConditionsType': [1, 2, 2, 2, 2, 2] }
		
		[pressure_problemData, pressure_dofs, pressure_strongs] = gedim.Discretize(pressure_discreteSpace, lib)
		[speed_problemData, speed_dofs, speed_strongs] = gedim.Discretize(speed_discreteSpace, lib)

		pressure_n_dofs = pressure_problemData['NumberDOFs']
		pressure_n_strongs = pressure_problemData['NumberStrongs']
		speed_n_dofs = speed_problemData['NumberDOFs']
		speed_n_strongs = speed_problemData['NumberStrongs']

		[stiffness_1, stiffnessStrong_1] = gedim.AssembleStiffnessMatrix_Shift(speed_problemData['SpaceIndex'], speed_problemData['SpaceIndex'], Stokes_v, 2 * speed_n_dofs + pressure_n_dofs, 2 * speed_n_dofs + pressure_n_dofs, 2 * speed_n_strongs + pressure_n_strongs, 0, 0, 0, lib)
		[stiffness_2, stiffnessStrong_2] = gedim.AssembleStiffnessMatrix_Shift(speed_problemData['SpaceIndex'], speed_problemData['SpaceIndex'], Stokes_v, 2 * speed_n_dofs + pressure_n_dofs, 2 * speed_n_dofs + pressure_n_dofs, 2 * speed_n_strongs + pressure_n_strongs, speed_n_dofs, speed_n_dofs, speed_n_strongs, lib)

		[advection_1, advectionStrong_1] = gedim.AssembleAdvectionMatrix_Shift(speed_problemData['SpaceIndex'], pressure_problemData['SpaceIndex'], Stokes_advection_1, 2 * speed_n_dofs + pressure_n_dofs, 2 * speed_n_dofs + pressure_n_dofs, 2 * speed_n_strongs + pressure_n_strongs, 2 * speed_n_dofs, 0, 0, lib)
		[advection_2, advectionStrong_2] = gedim.AssembleAdvectionMatrix_Shift(speed_problemData['SpaceIndex'], pressure_problemData['SpaceIndex'], Stokes_advection_2, 2 * speed_n_dofs + pressure_n_dofs, 2 * speed_n_dofs + pressure_n_dofs, 2 * speed_n_strongs + pressure_n_strongs, 2 * speed_n_dofs, speed_n_dofs, speed_n_strongs, lib)

		forcingTerm_1 = gedim.AssembleForcingTerm(Stokes_f_1, speed_problemData, lib)
		forcingTerm_2 = gedim.AssembleForcingTerm(Stokes_f_2, speed_problemData, lib)
		forcingTerm = np.concatenate([forcingTerm_1, forcingTerm_2, np.zeros(pressure_n_dofs)])

		pressure_solutionStrong = gedim.AssembleStrongSolution(Stokes_pressure_exactSolution, 1, pressure_problemData, lib)
		
		solution = gedim.LUSolver(stiffness_1 + stiffness_2 \
			     - advection_1 - advection_2 \
				 - np.transpose(advection_1) - np.transpose(advection_2), forcingTerm, lib)

		pressure_errorL2 = gedim.ComputeErrorL2(Stokes_pressure_exactSolution, solution[2 * speed_n_dofs:], pressure_solutionStrong, lib, pressure_problemData)

		print("dofs", "h", "pressure_errorL2")
		print(2 * speed_n_dofs + pressure_n_dofs, '{:.16e}'.format(pressure_problemData['H']), '{:.16e}'.format(pressure_errorL2))

		u = solution[0:2 * speed_n_dofs]
		p = solution[2 * speed_n_dofs:]

		if activatePlot:
			gedim.PlotSolution(mesh, pressure_dofs, pressure_strongs, p, pressure_solutionStrong, "Pressure")
			gedim.PlotSolution(mesh, speed_dofs, speed_strongs, u[0:speed_n_dofs], np.zeros(speed_n_strongs), "Speed X")
			gedim.PlotSolution(mesh, speed_dofs, speed_strongs, u[speed_n_dofs:], np.zeros(speed_n_strongs), "Speed Y")
			gedim.PlotSolution(mesh, speed_dofs, speed_strongs, np.sqrt(u[0:speed_n_dofs] * u[0:speed_n_dofs] + u[speed_n_dofs:] * u[speed_n_dofs:]), np.zeros(speed_n_strongs), "Speed Magnitude")

		[X_1, XStrong_1] = gedim.AssembleStiffnessMatrix_Shift(speed_problemData['SpaceIndex'], speed_problemData['SpaceIndex'], Stokes_v, 2 * speed_n_dofs, 2 * speed_n_dofs, 2 * speed_n_strongs, 0, 0, 0, lib)
		[X_2, XStrong_2] = gedim.AssembleStiffnessMatrix_Shift(speed_problemData['SpaceIndex'], speed_problemData['SpaceIndex'], Stokes_v, 2 * speed_n_dofs, 2 * speed_n_dofs, 2 * speed_n_strongs, speed_n_dofs, speed_n_dofs, speed_n_strongs, lib)

		[B_1, BStrong_1] = gedim.AssembleAdvectionMatrix_Shift(speed_problemData['SpaceIndex'], pressure_problemData['SpaceIndex'], Stokes_advection_1, pressure_n_dofs, 2 * speed_n_dofs, 2 * speed_n_strongs, 0, 0, 0, lib)
		[B_2, BStrong_2] = gedim.AssembleAdvectionMatrix_Shift(speed_problemData['SpaceIndex'], pressure_problemData['SpaceIndex'], Stokes_advection_2, pressure_n_dofs, 2 * speed_n_dofs, 2 * speed_n_strongs, 0, speed_n_dofs, speed_n_strongs, lib)
		
		supremizer = gedim.LUSolver(X_1 + X_2, np.transpose(B_1 + B_2) @ p, lib)

		if activatePlot:
			gedim.PlotSolution(mesh, speed_dofs, speed_strongs, supremizer[0:speed_n_dofs], np.zeros(speed_n_strongs), "Supremizer X")
			gedim.PlotSolution(mesh, speed_dofs, speed_strongs, supremizer[speed_n_dofs:], np.zeros(speed_n_strongs), "Supremizer Y")