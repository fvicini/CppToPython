import ctypes as ct
import numpy as np
import scipy.sparse
import sys

def make_nd_array(c_pointer, shape, dtype=np.double, order='F', own_data=True):
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

def Poisson_k(numPoints, points):
	values = np.ones((1, numPoints))
	return values.ctypes.data

def ImportLibrary(path):
	return ct.cdll.LoadLibrary(path)

def Initialize(config):
	lib.GedimForPy_Initialize.argtypes = [ct.py_object]
	lib.GedimForPy_Initialize.restype = None
	lib.GedimForPy_Initialize(config)
	
def CreateDomainSquare(domain):
	lib.GedimForPy_CreateDomainSquare.argtypes = [ct.py_object]
	lib.GedimForPy_CreateDomainSquare.restype = None
	lib.GedimForPy_CreateDomainSquare(domain)
	
def Discretize(discreteSpace):
	lib.GedimForPy_Discretize.argtypes = [ct.py_object]
	lib.GedimForPy_Discretize.restype = ct.py_object
	return lib.GedimForPy_Discretize(discreteSpace)
	
def AssembleStiffnessMatrix(problemData):
	DiffusionFN = ct.CFUNCTYPE(np.ctypeslib.ndpointer(dtype=np.double), ct.c_int, np.ctypeslib.ndpointer(dtype=np.double))
		
	lib.GedimForPy_AssembleStiffnessMatrix.argtypes = [DiffusionFN, ct.POINTER(ct.c_int), ct.POINTER(ct.POINTER(ct.c_double))]
	lib.GedimForPy_AssembleStiffnessMatrix.restype =  None
	
	pointerT = ct.POINTER(ct.c_double)()
	numTriplets = ct.c_int(0)
	lib.GedimForPy_AssembleStiffnessMatrix(DiffusionFN(Poisson_k), ct.byref(numTriplets),  ct.byref(pointerT))
	numTriplets = numTriplets.value
	triplets = make_nd_array(pointerT, (3, numTriplets))
	
	numDofs = problemData['NumberDOFs']
	return scipy.sparse.csr_matrix((triplets[2,:], (triplets[0,:], triplets[1,:])), shape=(numDofs, numDofs))

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
	domain = { 'SquareEdge': 1.0, 'VerticesBoundaryCondition': [1,1,1,1], 'EdgesBoundaryCondition': [1,1,1,1], 'DiscretizationType': 1, 'MeshCellsMaximumArea': 0.1 }
	CreateDomainSquare(domain)
	print("CreateDomainSquare successful")

	print("Discretize...")
	discreteSpace = { 'Order': 2, 'Type': 1, 'BoundaryConditionsType': [1, 2] }
	problemData = Discretize(discreteSpace)
	print(problemData)
	print("Discretize successful")

	print("AssembleStiffnessMatrix...")
	stiffness = AssembleStiffnessMatrix(problemData)
	print(stiffness)
	print("AssembleStiffnessMatrix successful")

	print("Test successful")
