import ctypes as ct
import numpy as np
import sys

def change_list(elements):
    return ct.cast((ct.c_int * len(elements))(*elements), ct.POINTER(ct.c_int))

def make_nd_array(c_pointer, shape, dtype=np.double, order='C', own_data=True):
    arr_size = np.prod(shape[:]) * np.dtype(dtype).itemsize 
    print("size", arr_size)
    
    if sys.version_info.major >= 3:
        ct.pythonapi.PyMemoryView_FromMemory.restype = ct.py_object
        ct.pythonapi.PyMemoryView_FromMemory.argtypes = (ct.c_void_p, ct.c_int, ct.c_int)
        buffer = ct.pythonapi.PyMemoryView_FromMemory(c_pointer, arr_size, 0x100)
    else:
        ct.pythonapi.PyBuffer_FromMemory.restype = ct.py_object
        buffer = ct.pythonapi.PyBuffer_FromMemory(c_pointer, arr_size)
    
    print(tuple(shape[:]), buffer)
    arr = np.ndarray(tuple(shape[:]), dtype, buffer, order=order)
    if own_data and not arr.flags.owndata:
        return arr.copy()
    else:
        return arr

def example_fn(numPoints, points):        
    print('numPoints', numPoints)
    print('points', points)
    
    vec = make_nd_array(points, (3,numPoints), np.double, 'C', True)
    print('vec', vec)
    print('shape', vec.shape)
    print(type(vec))
        
    result = (vec[0,:]**2 + vec[1,:]**2) # [x^2 + y^2]
    print('result', result)
    print('shape', result.shape)
    print(type(result))
        
    return result.ctypes.data

if __name__ == '__main__':
	print("Importing library...")
	lib = ct.cdll.LoadLibrary("./release/hello.so")
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
	
	ExampleFN = ct.CFUNCTYPE(np.ctypeslib.ndpointer(dtype=np.double), ct.c_int, np.ctypeslib.ndpointer(dtype=np.double))
	lib.PassFunctionPointer.argtypes = [ExampleFN]
	lib.PassFunctionPointer.restype = None

	lib.PassFunctionPointer(ExampleFN(example_fn))
	print("Function pointer successfull")
	
	lib.CreateMatrix.argtypes = [ct.c_int, ct.c_int]
	lib.CreateMatrix.restype = np.ctypeslib.ndpointer(dtype=np.double)
	A = make_nd_array(lib.CreateMatrix(3, 4), (3, 4), np.double)
	print(A)

	print("Call successful")
