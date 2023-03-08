import ctypes as ct
import numpy as np
import sys

def Initialize(config):
	lib.GedimForPy_Initialize.argtypes = [ct.py_object]
	lib.GedimForPy_Initialize.restype = None
	lib.GedimForPy_Initialize(config)
	
def CreateDomainSquare(domain):
	lib.GedimForPy_CreateDomainSquare.argtypes = [ct.py_object]
	lib.GedimForPy_CreateDomainSquare.restype = None
	lib.GedimForPy_CreateDomainSquare(domain)

if __name__ == '__main__':

	print("Importing library...")
	lib = ct.cdll.LoadLibrary("./release/GeDiM4Py.so")
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

	print("Test successful")
