#ifndef __TEST_PYTHON_H
#define __TEST_PYTHON_H

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>

#include "GeDiM4Py_Interface.hpp"

namespace UnitTesting
{
  // ***************************************************************************
  TEST(TestPython, Test_Dict)
  {
    PyObject* dict = PyDict_New();
    PyDict_SetItemString(dict, "GeometricTolerance", Py_BuildValue("d", 1.0e-8));

    ASSERT_NO_THROW(GedimForPy_Initialize(dict));
  }
  // ***************************************************************************
}

#endif // __TEST_PYTHON_H
