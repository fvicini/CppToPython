import numpy as np
import GeDiM4Py as gedim

def Heat_R():
	return 0.5
def Heat_K():
	return 6.68
def Heat_G():
	return 0.94

def Heat_k(numPoints, points):
	matPoints = gedim.make_nd_matrix(points, (3, numPoints), np.double)
	values = np.ones(numPoints)
	for p in range(0, numPoints):
		if (matPoints[0,p] * matPoints[0,p] + matPoints[1,p] * matPoints[1,p]) <= (Heat_R() * Heat_R() + 1.0e-16):
			values[p] = Heat_K()
	return values.ctypes.data

def Heat_weakTerm_down(numPoints, points):
	values = np.ones(numPoints) * Heat_G()
	return values.ctypes.data
	
if __name__ == '__main__':

	lib = gedim.ImportLibrary("/home/geoscore/Desktop/GEO++/Courses/CppToPython/release/GeDiM4Py.so")

	config = { 'GeometricTolerance': 1.0e-8 }
	gedim.Initialize(config, lib)
	
	order = 2

	[meshInfo, mesh] = gedim.ImportDomainMesh2D(lib)

	gedim.PlotMesh(mesh)

	discreteSpace = { 'Order': order, 'Type': 1, 'BoundaryConditionsType': [1, 3, 3, 2] }
	[problemData, dofs, strongs] = gedim.Discretize(discreteSpace, lib)

	gedim.PlotDofs(mesh, dofs, strongs)

	[stiffness, stiffnessStrong] = gedim.AssembleStiffnessMatrix(Heat_k, problemData, lib)
	
	weakTerm_down = gedim.AssembleWeakTerm(Heat_weakTerm_down, 1, problemData, lib)

	solution = gedim.LUSolver(stiffness, weakTerm_down, lib)

	gedim.PlotSolution(mesh, dofs, strongs, solution, np.zeros(problemData['NumberStrongs']))