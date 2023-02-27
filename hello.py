import ctypes as ct
import numpy as np

def change_list(elements):
    return ct.cast((ct.c_int * len(elements))(*elements), ct.POINTER(ct.c_int))


if __name__ == '__main__':
	print("Importing library...")
	lib = ct.cdll.LoadLibrary("./debug/hello.so")
	print(lib)
	print("Import successful")

	print("Calling C++...")
	lib.hello()
	lib.testInt(102)

	lib.testString.restype = ct.c_char_p
	lib.testString.argtypes = [ct.c_char_p]
	name = ct.create_string_buffer(b"Tom")
	result = lib.testString(name)
	print(result)

	lib.testTuple.restype = ct.py_object
	result = lib.testTuple()
	print(result)

	lib.testList.restype = ct.py_object
	result = lib.testList()
	print(type(result))
	print(result)

	l = [1,2,3,4,5]
	lib.sum.restype = ct.c_longlong
	lib.sum.argtypes = [ct.c_int, ct.POINTER(ct.c_int)]
	result = lib.sum(len(l), change_list(l))
	print(result)

	lib.sum.restype = ct.c_longlong
	lib.sum.argtypes = [ct.c_int, 
                        np.ctypeslib.ndpointer(dtype=np.int32)]

	array = np.arange(0, 1000, 1, np.int32)
	result = lib.sum(len(array), array)
	print(result)

	print("Call successful")
