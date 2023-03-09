import ctypes as ct
import numpy as np
import sys

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
	ExampleFN = ct.CFUNCTYPE(np.ctypeslib.ndpointer(dtype=np.double), ct.c_int, np.ctypeslib.ndpointer(dtype=np.double))

	lib.GedimForPy_Discretize.argtypes = [ct.py_object]
	lib.GedimForPy_Discretize.restype = ct.py_object
	return lib.GedimForPy_Discretize(discreteSpace)
	
def AssembleStiffnessMatrix(discreteSpace):
	lib.GedimForPy_AssembleStiffnessMatrix.argtypes = [ct.py_object]
	lib.GedimForPy_AssembleStiffnessMatrix.restype = ct.py_object
	triplets = lib.GedimForPy_AssembleStiffnessMatrix()
	return stiffness

if __name__ == '__main__':

	print("Importing library...")
	lib = ImportLibrary("./release/GeDiM4Py.so")
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

	print("Discretize...")
	stiffness = AssembleStiffnessMatrix()
	print(stiffness)
	print("Discretize successful")

	print("Test successful")
