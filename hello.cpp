#include <Python.h>

#include <iostream>
#include <vector>

#include "IOUtilities.hpp"

std::string print(const std::string& what);

extern "C"
void hello() { std::cout << "Hello world!"<< std::endl; }

extern "C"
void testInt(const int test)
{
  Gedim::Output::PrintStars(std::cout);
  std::cout<< "Get "<< test<< std::endl;
  Gedim::Output::PrintStars(std::cout);
}

extern "C"
const char* testString(const char* name)
{
  std::cout<<strlen(name)<< std::endl;
  std::cout<<name<< std::endl;
  static std::string s = "Test ";
  s += name;
  return s.c_str();
}

extern "C"
PyObject* testTuple()
{
  PyObject* tupleOne = Py_BuildValue("(id)",1,2.2);
  PyObject* tupleTwo = Py_BuildValue("(ii)",3,4);
  PyObject* data = Py_BuildValue("(OO)", tupleOne, tupleTwo);

  return data;
}

extern "C"
PyObject* testList()
{
  PyObject *PList = PyList_New(0);

  std::vector<int> intVector;
  std::vector<int>::const_iterator it;

  for(int i = 0 ; i < 10 ; i++)
    intVector.push_back(i);

  for(const int& value : intVector)
    PyList_Append(PList, Py_BuildValue("i", value));

  return PList;
}

extern "C"
void testInputTuple(PyObject* args)
{
  std::cerr<< "argType->ob_refcnt "<< args->ob_refcnt<< std::endl;
  std::cerr<< "argType->tp_name "<< args->ob_type->tp_name<< std::endl;

  int test = 0;
  PyArg_ParseTuple(args, "(i)", &test);
  std::cout<< "Receive "<< test<< std::endl;
}

extern "C"
long long sum(const int n,
              const int* array)
{
  long long res = 0;
  for (int i = 0; i < n; ++i)
    res += array[i];

  return res;
}

int main() 
{
  hello();
  std::cout << print("Hello world!")<< std::endl;
  const std::string test = "test";
  return 0;
}

std::string print(const std::string& what)
{
  std::cout << what<< std::endl;
  return what + " Done";
}
