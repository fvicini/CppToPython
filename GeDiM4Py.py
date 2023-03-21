import ctypes as ct
import numpy as np
import scipy.sparse
import sys
import matplotlib.pyplot as plt
import matplotlib.tri
import os

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

def ImportLibrary(path):
	return ct.cdll.LoadLibrary(path)

def Initialize(config, lib):
	lib.GedimForPy_Initialize.argtypes = [ct.py_object]
	lib.GedimForPy_Initialize.restype = None
	lib.GedimForPy_Initialize(config)
	
def CreateDomainSquare(domain, lib):
	lib.GedimForPy_CreateDomainSquare.argtypes = [ct.py_object, ct.POINTER(ct.POINTER(ct.c_double))]
	lib.GedimForPy_CreateDomainSquare.restype = ct.py_object

	pointerMeshCoordinates = ct.POINTER(ct.c_double)()
	meshInfo = lib.GedimForPy_CreateDomainSquare(domain, pointerMeshCoordinates)

	mesh = make_nd_matrix(pointerMeshCoordinates, (3, meshInfo['NumberCell0Ds']))
	return [meshInfo, mesh]

def CreateDomainRectangle(domain, lib):
	lib.GedimForPy_CreateDomainRectangle.argtypes = [ct.py_object, ct.POINTER(ct.POINTER(ct.c_double))]
	lib.GedimForPy_CreateDomainRectangle.restype = ct.py_object

	pointerMeshCoordinates = ct.POINTER(ct.c_double)()
	meshInfo = lib.GedimForPy_CreateDomainRectangle(domain, pointerMeshCoordinates)

	mesh = make_nd_matrix(pointerMeshCoordinates, (3, meshInfo['NumberCell0Ds']))
	return [meshInfo, mesh]

def ImportDomainMesh2D(lib):
	lib.GedimForPy_ImportDomainMesh2D.argtypes = [ct.POINTER(ct.POINTER(ct.c_double))]
	lib.GedimForPy_ImportDomainMesh2D.restype = ct.py_object

	pointerMeshCoordinates = ct.POINTER(ct.c_double)()
	meshInfo = lib.GedimForPy_ImportDomainMesh2D(pointerMeshCoordinates)

	mesh = make_nd_matrix(pointerMeshCoordinates, (3, meshInfo['NumberCell0Ds']))
	return [meshInfo, mesh]

def Discretize(discreteSpace, lib):
	lib.GedimForPy_Discretize.argtypes = [ct.py_object, ct.POINTER(ct.POINTER(ct.c_double)), ct.POINTER(ct.POINTER(ct.c_double))]
	lib.GedimForPy_Discretize.restype = ct.py_object

	pointerDofs = ct.POINTER(ct.c_double)()
	pointerStrongs = ct.POINTER(ct.c_double)()
	problemData = lib.GedimForPy_Discretize(discreteSpace, pointerDofs, pointerStrongs)

	dofs = make_nd_matrix(pointerDofs, (3, problemData['NumberDOFs']))
	strongs = make_nd_matrix(pointerStrongs, (3, problemData['NumberStrongs']))

	return [problemData, dofs, strongs]
	
def AssembleStiffnessMatrix(a, problemData, lib):
	DiffusionFN = ct.CFUNCTYPE(np.ctypeslib.ndpointer(dtype=np.double), ct.c_int, np.ctypeslib.ndpointer(dtype=np.double))
		
	lib.GedimForPy_AssembleStiffnessMatrix.argtypes = [DiffusionFN, ct.POINTER(ct.c_int), ct.POINTER(ct.POINTER(ct.c_double)), ct.POINTER(ct.c_int), ct.POINTER(ct.POINTER(ct.c_double))]
	lib.GedimForPy_AssembleStiffnessMatrix.restype =  None
	
	pointerStiffness = ct.POINTER(ct.c_double)()
	numStiffnessTriplets = ct.c_int(0)
	pointerStiffnessStrong = ct.POINTER(ct.c_double)()
	numStiffnessStrongTriplets = ct.c_int(0)
	lib.GedimForPy_AssembleStiffnessMatrix(DiffusionFN(a), ct.byref(numStiffnessTriplets), ct.byref(pointerStiffness), ct.byref(numStiffnessStrongTriplets), ct.byref(pointerStiffnessStrong))
	
	numDofs = problemData['NumberDOFs']
	numStrongs = problemData['NumberStrongs']

	stiffness = make_np_sparse(numDofs, numDofs, numStiffnessTriplets, pointerStiffness)
	stiffnessStrong = make_np_sparse(numDofs, numStrongs, numStiffnessStrongTriplets, pointerStiffnessStrong)

	return [stiffness, stiffnessStrong]

def AssembleAdvectionMatrix(b, problemData, lib):
	AdvectionFN = ct.CFUNCTYPE(np.ctypeslib.ndpointer(dtype=np.double), ct.c_int, np.ctypeslib.ndpointer(dtype=np.double))
		
	lib.GedimForPy_AssembleAdvectionMatrix.argtypes = [AdvectionFN, ct.POINTER(ct.c_int), ct.POINTER(ct.POINTER(ct.c_double)), ct.POINTER(ct.c_int), ct.POINTER(ct.POINTER(ct.c_double))]
	lib.GedimForPy_AssembleAdvectionMatrix.restype =  None
	
	pointerAdvection = ct.POINTER(ct.c_double)()
	numAdvectionTriplets = ct.c_int(0)
	pointerAdvectionStrong = ct.POINTER(ct.c_double)()
	numAdvectionStrongTriplets = ct.c_int(0)
	lib.GedimForPy_AssembleAdvectionMatrix(AdvectionFN(b), ct.byref(numAdvectionTriplets), ct.byref(pointerAdvection), ct.byref(numAdvectionStrongTriplets), ct.byref(pointerAdvectionStrong))
	
	numDofs = problemData['NumberDOFs']
	numStrongs = problemData['NumberStrongs']

	advection = make_np_sparse(numDofs, numDofs, numAdvectionTriplets, pointerAdvection)
	advectionStrong = make_np_sparse(numDofs, numStrongs, numAdvectionStrongTriplets, pointerAdvectionStrong)

	return [advection, advectionStrong]

def AssembleReactionMatrix(c, problemData, lib):
	ReactionFN = ct.CFUNCTYPE(np.ctypeslib.ndpointer(dtype=np.double), ct.c_int, np.ctypeslib.ndpointer(dtype=np.double))
		
	lib.GedimForPy_AssembleReactionMatrix.argtypes = [ReactionFN, ct.POINTER(ct.c_int), ct.POINTER(ct.POINTER(ct.c_double)), ct.POINTER(ct.c_int), ct.POINTER(ct.POINTER(ct.c_double))]
	lib.GedimForPy_AssembleReactionMatrix.restype =  None
	
	pointerReaction = ct.POINTER(ct.c_double)()
	numReactionTriplets = ct.c_int(0)
	pointerReactionStrong = ct.POINTER(ct.c_double)()
	numReactionStrongTriplets = ct.c_int(0)
	lib.GedimForPy_AssembleReactionMatrix(ReactionFN(c), ct.byref(numReactionTriplets), ct.byref(pointerReaction), ct.byref(numReactionStrongTriplets), ct.byref(pointerReactionStrong))
	
	numDofs = problemData['NumberDOFs']
	numStrongs = problemData['NumberStrongs']

	reaction = make_np_sparse(numDofs, numDofs, numReactionTriplets, pointerReaction)
	reactionStrong = make_np_sparse(numDofs, numStrongs, numReactionStrongTriplets, pointerReactionStrong)

	return [reaction, reactionStrong]

def AssembleForcingTerm(f, problemData, lib):
	ForcingTermFN = ct.CFUNCTYPE(np.ctypeslib.ndpointer(dtype=np.double), ct.c_int, np.ctypeslib.ndpointer(dtype=np.double))
		
	lib.GedimForPy_AssembleForcingTerm.argtypes = [ForcingTermFN, ct.POINTER(ct.c_int), ct.POINTER(ct.POINTER(ct.c_double))]
	lib.GedimForPy_AssembleForcingTerm.restype =  None
	
	pointerF = ct.POINTER(ct.c_double)()
	size = ct.c_int(0)
	lib.GedimForPy_AssembleForcingTerm(ForcingTermFN(f), ct.byref(size), ct.byref(pointerF))
	size = size.value
	return make_nd_array(pointerF, size)

def AssembleStrongSolution(g, marker, problemData, lib):
	StrongFN = ct.CFUNCTYPE(np.ctypeslib.ndpointer(dtype=np.double), ct.c_int, np.ctypeslib.ndpointer(dtype=np.double))
		
	lib.GedimForPy_AssembleStrongSolution.argtypes = [StrongFN, ct.c_int, ct.POINTER(ct.c_int), ct.POINTER(ct.POINTER(ct.c_double))]
	lib.GedimForPy_AssembleStrongSolution.restype =  None
	
	pointerStrongSolution = ct.POINTER(ct.c_double)()
	size = ct.c_int(0)
	lib.GedimForPy_AssembleStrongSolution(StrongFN(g), marker, ct.byref(size), ct.byref(pointerStrongSolution))
	size = size.value
	return make_nd_array(pointerStrongSolution, size)

def AssembleWeakTerm(g, marker, problemData, lib):
	WeakTermFN = ct.CFUNCTYPE(np.ctypeslib.ndpointer(dtype=np.double), ct.c_int, np.ctypeslib.ndpointer(dtype=np.double))
		
	lib.GedimForPy_AssembleWeakTerm.argtypes = [WeakTermFN, ct.c_int, ct.POINTER(ct.c_int), ct.POINTER(ct.POINTER(ct.c_double))]
	lib.GedimForPy_AssembleWeakTerm.restype =  None
	
	pointerWeak = ct.POINTER(ct.c_double)()
	size = ct.c_int(0)
	lib.GedimForPy_AssembleWeakTerm(WeakTermFN(g), marker, ct.byref(size), ct.byref(pointerWeak))
	size = size.value
	return make_nd_array(pointerWeak, size)

def CholeskySolver(A, f, lib):
	[rows, cols, values] = scipy.sparse.find(A)
	nonZerosA = np.column_stack((rows, cols, values))
	lib.GedimForPy_CholeskySolver.argtypes = [ct.c_int, ct.c_int, np.ctypeslib.ndpointer(dtype=np.double), np.ctypeslib.ndpointer(dtype=np.double), ct.POINTER(ct.POINTER(ct.c_double))]
	lib.GedimForPy_CholeskySolver.restype =  None

	pointerSolution = ct.POINTER(ct.c_double)()

	lib.GedimForPy_CholeskySolver(A.shape[0], rows.shape[0], nonZerosA, f, ct.byref(pointerSolution))

	return make_nd_array(pointerSolution, A.shape[0])

def LUSolver(A, f, lib):
	[rows, cols, values] = scipy.sparse.find(A)
	nonZerosA = np.column_stack((rows, cols, values))
	lib.GedimForPy_LUSolver.argtypes = [ct.c_int, ct.c_int, np.ctypeslib.ndpointer(dtype=np.double), np.ctypeslib.ndpointer(dtype=np.double), ct.POINTER(ct.POINTER(ct.c_double))]
	lib.GedimForPy_LUSolver.restype =  None

	pointerSolution = ct.POINTER(ct.c_double)()

	lib.GedimForPy_LUSolver(A.shape[0], rows.shape[0], nonZerosA, f, ct.byref(pointerSolution))

	return make_nd_array(pointerSolution, A.shape[0])

def ComputeErrorL2(u, solution, solutionStrong, lib):
	ExactFN = ct.CFUNCTYPE(np.ctypeslib.ndpointer(dtype=np.double), ct.c_int, np.ctypeslib.ndpointer(dtype=np.double))
	lib.GedimForPy_ComputeErrorL2.argtypes = [ExactFN, np.ctypeslib.ndpointer(dtype=np.double), np.ctypeslib.ndpointer(dtype=np.double)]
	lib.GedimForPy_ComputeErrorL2.restype =  ct.c_double

	return lib.GedimForPy_ComputeErrorL2(ExactFN(u), solution, solutionStrong)

def ComputeErrorH1(uDer, solution, solutionStrong, lib):
	ExactDerivativeFN = ct.CFUNCTYPE(np.ctypeslib.ndpointer(dtype=np.double), ct.c_int, ct.c_int, np.ctypeslib.ndpointer(dtype=np.double))
	lib.GedimForPy_ComputeErrorH1.argtypes = [ExactDerivativeFN, np.ctypeslib.ndpointer(dtype=np.double), np.ctypeslib.ndpointer(dtype=np.double)]
	lib.GedimForPy_ComputeErrorH1.restype =  ct.c_double

	return lib.GedimForPy_ComputeErrorH1(ExactDerivativeFN(uDer), solution, solutionStrong)

def PythonSolver(A, f, lib):
	return scipy.sparse.linalg.spsolve(A, f)

def PlotMesh(mesh):	
	fig = plt.figure(figsize=plt.figaspect(0.5))

	ax1 = fig.add_subplot(1, 1, 1)
	ax1.set_aspect('equal')
	ax1.triplot(matplotlib.tri.Triangulation(mesh[0, :], mesh[1, :]), 'ko-', lw=1)
	ax1.grid(True)
		
	current_directory_path = os.getcwd()
	subfolder_path = os.path.join(current_directory_path, 'Images')
	if not os.path.exists(subfolder_path):
		os.makedirs(subfolder_path)
	file_name = 'Mesh.png'
	file_path = os.path.join(subfolder_path, file_name)
	print(file_path)
	plt.savefig(file_path)
	
	plt.show()

	plt.close()

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
