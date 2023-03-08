import ctypes as ct
import numpy as np
import sys

def Initialize(config):
	lib.GedimForPy_Initialize.argtypes = [ct.py_object]
	lib.GedimForPy_Initialize.restype = None
	lib.GedimForPy_Initialize(config)

if __name__ == '__main__':

	print("Importing library...")
	lib = ct.cdll.LoadLibrary("./release/GeDiM4Py.so")
	print(lib)
	print("Import library successful")

	print("Initialize...")
	config = { 'GeometricTolerance': 1.0e-8, }
	Initialize(config)
	print("Initialize successful")

	print("Test successful")
