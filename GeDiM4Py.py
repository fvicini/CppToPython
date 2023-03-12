import ctypes as ct
import numpy as np
import scipy.sparse
import sys
import matplotlib.pyplot as plt
import matplotlib.tri

def make_np_sparse(nRows, nCols, c_nNonZeros, c_pointerTriplets):
	nNonZeros = c_nNonZeros.value
	triplets = make_nd_matrix(c_pointerTriplets, (3, nNonZeros))
	
	return scipy.sparse.csr_matrix((triplets[2,:], (triplets[0,:], triplets[1,:])), shape=(nRows, nCols))

def make_nd_matrix(c_pointer, shape, dtype=np.double, order='F', own_data=True):
    arr_size = np.prod(shape[:]) * np.dtype(dtype).itemsize 
    
    if sys.version_info.major >= 3:
        ct.pythonapi.PyMemoryView_FromMemory.restype = ct.py_object
        ct.pythonapi.PyMemoryView_FromMemory.argtypes = (ct.c_void_p, ct.c_int, ct.c_int)
        buffer = ct.pythonapi.PyMemoryView_FromMemory(c_pointer, arr_size, 0x100)
    else:
        ct.pythonapi.PyBuffer_FromMemory.restype = ct.py_object
        buffer = ct.pythonapi.PyBuffer_FromMemory(c_pointer, arr_size)
        
    arr = np.ndarray(tuple(shape[:]), dtype, buffer, order=order)
    if own_data and not arr.flags.owndata:
        return arr.copy()
    else:
        return arr

def make_nd_array(c_pointer, size, dtype=np.double, order='F', own_data=True):
    arr_size = np.prod(size) * np.dtype(dtype).itemsize 
    
    if sys.version_info.major >= 3:
        ct.pythonapi.PyMemoryView_FromMemory.restype = ct.py_object
        ct.pythonapi.PyMemoryView_FromMemory.argtypes = (ct.c_void_p, ct.c_int, ct.c_int)
        buffer = ct.pythonapi.PyMemoryView_FromMemory(c_pointer, arr_size, 0x100)
    else:
        ct.pythonapi.PyBuffer_FromMemory.restype = ct.py_object
        buffer = ct.pythonapi.PyBuffer_FromMemory(c_pointer, arr_size)
        
    arr = np.ndarray(size, dtype, buffer, order=order)
    if own_data and not arr.flags.owndata:
        return arr.copy()
    else:
        return arr

def Poisson_k(numPoints, points):
	values = np.ones((1, numPoints))
	return values.ctypes.data

def Poisson_f(numPoints, points):
	matPoints = make_nd_matrix(points, (3, numPoints), np.double)
	values = 32.0 * (matPoints[1,:] * (1.0 - matPoints[1,:]) + matPoints[0,:] * (1.0 - matPoints[0,:]))
	return values.ctypes.data

def Poisson_exactSolution(numPoints, points):
	matPoints = make_nd_matrix(points, (3, numPoints), np.double)
	values = 16.0 * (matPoints[1,:] * (1.0 - matPoints[1,:]) * matPoints[0,:] * (1.0 - matPoints[0,:])) + 1.1
	return values.ctypes.data
	
def Poisson_weakTerm_right(numPoints, points):
	matPoints = make_nd_matrix(points, (3, numPoints), np.double)
	values = 16.0 * (1.0 - 2.0 * matPoints[0,:]) * matPoints[1,:] * (1.0 - matPoints[1,:])
	return values.ctypes.data
	
def Poisson_weakTerm_left(numPoints, points):
	matPoints = make_nd_matrix(points, (3, numPoints), np.double)
	values = -16.0 * (1.0 - 2.0 * matPoints[0,:]) * matPoints[1,:] * (1.0 - matPoints[1,:])
	return values.ctypes.data

def ImportLibrary(path):
	return ct.cdll.LoadLibrary(path)

def Initialize(config):
	lib.GedimForPy_Initialize.argtypes = [ct.py_object]
	lib.GedimForPy_Initialize.restype = None
	lib.GedimForPy_Initialize(config)
	
def CreateDomainSquare(domain):
	lib.GedimForPy_CreateDomainSquare.argtypes = [ct.py_object, ct.POINTER(ct.POINTER(ct.c_double))]
	lib.GedimForPy_CreateDomainSquare.restype = ct.py_object

	pointerMeshCoordinates = ct.POINTER(ct.c_double)()
	meshInfo = lib.GedimForPy_CreateDomainSquare(domain, pointerMeshCoordinates)

	mesh = make_nd_matrix(pointerMeshCoordinates, (3, meshInfo['NumberCell0Ds']))
	return [meshInfo, mesh]
	
def Discretize(discreteSpace):
	lib.GedimForPy_Discretize.argtypes = [ct.py_object, ct.POINTER(ct.POINTER(ct.c_double)), ct.POINTER(ct.POINTER(ct.c_double))]
	lib.GedimForPy_Discretize.restype = ct.py_object

	pointerDofs = ct.POINTER(ct.c_double)()
	pointerStrongs = ct.POINTER(ct.c_double)()
	problemData = lib.GedimForPy_Discretize(discreteSpace, pointerDofs, pointerStrongs)

	dofs = make_nd_matrix(pointerDofs, (3, problemData['NumberDOFs']))
	strongs = make_nd_matrix(pointerStrongs, (3, problemData['NumberStrongs']))

	return [problemData, dofs, strongs]
	
def AssembleStiffnessMatrix(k, problemData):
	DiffusionFN = ct.CFUNCTYPE(np.ctypeslib.ndpointer(dtype=np.double), ct.c_int, np.ctypeslib.ndpointer(dtype=np.double))
		
	lib.GedimForPy_AssembleStiffnessMatrix.argtypes = [DiffusionFN, ct.POINTER(ct.c_int), ct.POINTER(ct.POINTER(ct.c_double)), ct.POINTER(ct.c_int), ct.POINTER(ct.POINTER(ct.c_double))]
	lib.GedimForPy_AssembleStiffnessMatrix.restype =  None
	
	pointerStifness = ct.POINTER(ct.c_double)()
	numStiffnessTriplets = ct.c_int(0)
	pointerStifnessStrong = ct.POINTER(ct.c_double)()
	numStiffnessStrongTriplets = ct.c_int(0)
	lib.GedimForPy_AssembleStiffnessMatrix(DiffusionFN(k), ct.byref(numStiffnessTriplets), ct.byref(pointerStifness), ct.byref(numStiffnessStrongTriplets), ct.byref(pointerStifnessStrong))
	
	numDofs = problemData['NumberDOFs']
	numStrongs = problemData['NumberStrongs']

	stifness = make_np_sparse(numDofs, numDofs, numStiffnessTriplets, pointerStifness)
	stifnessStrong = make_np_sparse(numDofs, numStrongs, numStiffnessStrongTriplets, pointerStifnessStrong)

	return [stifness, stifnessStrong]

def AssembleForcingTerm(f, problemData):
	ForcingTermFN = ct.CFUNCTYPE(np.ctypeslib.ndpointer(dtype=np.double), ct.c_int, np.ctypeslib.ndpointer(dtype=np.double))
		
	lib.GedimForPy_AssembleForcingTerm.argtypes = [ForcingTermFN, ct.POINTER(ct.c_int), ct.POINTER(ct.POINTER(ct.c_double))]
	lib.GedimForPy_AssembleForcingTerm.restype =  None
	
	pointerF = ct.POINTER(ct.c_double)()
	size = ct.c_int(0)
	lib.GedimForPy_AssembleForcingTerm(ForcingTermFN(f), ct.byref(size), ct.byref(pointerF))
	size = size.value
	return make_nd_array(pointerF, size)

def AssembleStrongSolution(g, marker, problemData):
	StrongFN = ct.CFUNCTYPE(np.ctypeslib.ndpointer(dtype=np.double), ct.c_int, np.ctypeslib.ndpointer(dtype=np.double))
		
	lib.GedimForPy_AssembleStrongSolution.argtypes = [StrongFN, ct.c_int, ct.POINTER(ct.c_int), ct.POINTER(ct.POINTER(ct.c_double))]
	lib.GedimForPy_AssembleStrongSolution.restype =  None
	
	pointerStrongSolution = ct.POINTER(ct.c_double)()
	size = ct.c_int(0)
	lib.GedimForPy_AssembleStrongSolution(StrongFN(g), marker, ct.byref(size), ct.byref(pointerStrongSolution))
	size = size.value
	return make_nd_array(pointerStrongSolution, size)

def AssembleWeakTerm(g, marker, problemData):
	WeakTermFN = ct.CFUNCTYPE(np.ctypeslib.ndpointer(dtype=np.double), ct.c_int, np.ctypeslib.ndpointer(dtype=np.double))
		
	lib.GedimForPy_AssembleWeakTerm.argtypes = [WeakTermFN, ct.c_int, ct.POINTER(ct.c_int), ct.POINTER(ct.POINTER(ct.c_double))]
	lib.GedimForPy_AssembleWeakTerm.restype =  None
	
	pointerWeak = ct.POINTER(ct.c_double)()
	size = ct.c_int(0)
	lib.GedimForPy_AssembleWeakTerm(WeakTermFN(g), marker, ct.byref(size), ct.byref(pointerWeak))
	size = size.value
	return make_nd_array(pointerWeak, size)

def CholeskySolver(A, f):
	[rows, cols, values] = scipy.sparse.find(A)
	nonZerosA = np.column_stack((rows, cols, values))
	lib.GedimForPy_CholeskySolver.argtypes = [ct.c_int, ct.c_int, np.ctypeslib.ndpointer(dtype=np.double), np.ctypeslib.ndpointer(dtype=np.double), ct.POINTER(ct.POINTER(ct.c_double))]
	lib.GedimForPy_CholeskySolver.restype =  None

	pointerSolution = ct.POINTER(ct.c_double)()

	lib.GedimForPy_CholeskySolver(A.shape[0], rows.shape[0], nonZerosA, f, ct.byref(pointerSolution))

	return make_nd_array(pointerSolution, A.shape[0])

def Solver(A, f):
	return scipy.sparse.linalg.spsolve(A, f)

def PlotMesh(mesh):	
	fig = plt.figure(figsize=plt.figaspect(0.5))

	ax1 = fig.add_subplot(1, 1, 1)
	ax1.set_aspect('equal')
	ax1.triplot(matplotlib.tri.Triangulation(mesh[0, :], mesh[1, :]), 'ko-', lw=1)
	ax1.grid(True)
	
	plt.show()

def PlotDofs(mesh, dofs, strongs):
	x = np.concatenate((dofs[0,:], strongs[0,:]))
	y = np.concatenate((dofs[1,:], strongs[1,:]))
	z = np.concatenate((np.arange(0, dofs.shape[1]), np.arange(0, strongs.shape[1])))

	fig = plt.figure(figsize=plt.figaspect(0.5))

	ax1 = fig.add_subplot(1, 1, 1)
	ax1.set_aspect('equal')
	ax1.triplot(matplotlib.tri.Triangulation(mesh[0, :], mesh[1, :]), 'k--', lw=1)
	ax1.scatter(x, y, c=z)
	ax1.grid(True)

	plt.show()

def PlotSolution(mesh, dofs, strongs, solutionDofs, solutionStrongs):
	x = np.concatenate((dofs[0,:], strongs[0,:]), axis=0)
	y = np.concatenate((dofs[1,:], strongs[1,:]), axis=0)
	z = np.concatenate((solutionDofs, solutionStrongs), axis=0)
	triang = matplotlib.tri.Triangulation(x, y)
	
	fig = plt.figure(figsize=plt.figaspect(0.5))

	ax1 = fig.add_subplot(1, 2, 1)
	ax1.set_aspect('equal')
	tpc = ax1.tripcolor(triang, z, shading='flat')
	ax1.triplot(matplotlib.tri.Triangulation(mesh[0, :], mesh[1, :]), 'k--', lw=1)
	fig.colorbar(tpc)

	ax2 = fig.add_subplot(1, 2, 2, projection='3d')
	ax2.plot_trisurf(x, y, z, triangles=triang.triangles, cmap=plt.cm.Spectral)

	plt.show()

if __name__ == '__main__':

	print("Importing library...")
	lib = ImportLibrary("/home/geoscore/Desktop/GEO++/Courses/CppToPython/release/GeDiM4Py.so")
	print(lib)
	print("Import library successful")

	print("Initialize...")
	config = { 'GeometricTolerance': 1.0e-8 }
	Initialize(config)
	print("Initialize successful")
	
	print("CreateDomainSquare...")
	domain = { 'SquareEdge': 1.0, 'VerticesBoundaryCondition': [1,1,1,1], 'EdgesBoundaryCondition': [1,2,1,3], 'DiscretizationType': 1, 'MeshCellsMaximumArea': 0.01 }
	[meshInfo, mesh] = CreateDomainSquare(domain)
	print("CreateDomainSquare successful")

	PlotMesh(mesh)

	print("Discretize...")
	discreteSpace = { 'Order': 2, 'Type': 1, 'BoundaryConditionsType': [1, 2, 3, 3] }
	[problemData, dofs, strongs] = Discretize(discreteSpace)
	print("Discretize successful")

	PlotDofs(mesh, dofs, strongs)

	print("AssembleStiffnessMatrix...")
	[stiffness, stiffnessStrong] = AssembleStiffnessMatrix(Poisson_k, problemData)
	print("AssembleStiffnessMatrix successful")

	print("AssembleForcingTerm...")
	forcingTerm = AssembleForcingTerm(Poisson_f, problemData)
	print("AssembleForcingTerm successful")

	print("AssembleStrongSolution...")
	solutionStrong = AssembleStrongSolution(Poisson_exactSolution, 1, problemData)
	print("AssembleStrongSolution successful")
	
	print("AssembleWeakTerm...")
	weakTerm_right = AssembleWeakTerm(Poisson_weakTerm_right, 2, problemData)
	weakTerm_left = AssembleWeakTerm(Poisson_weakTerm_left, 3, problemData)
	print("AssembleWeakTerm successful")

	print("CholeskySolver...")
	solution = CholeskySolver(stiffness, forcingTerm - stiffnessStrong @ solutionStrong + weakTerm_right + weakTerm_left)
	print("CholeskySolver successful")

	PlotSolution(mesh, dofs, strongs, solution, solutionStrong)

	print("Test successful")
