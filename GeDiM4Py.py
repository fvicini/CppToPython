import ctypes as ct
import numpy as np
import sys

if __name__ == '__main__':

	print("Importing library...")
	lib = ct.cdll.LoadLibrary("./release/GeDiM4Py.so")
	print(lib)
	print("Import successful")



	print("Call successful")
