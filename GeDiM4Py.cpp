#include <gtest/gtest.h>

#include "test_Poisson.hpp"
#include "test_heat_conductivity.hpp"
#include "test_Stokes.hpp"
#include "test_Burger.hpp"
#include "test_NonLinearPoisson.hpp"
#include "test_python.hpp"

int main(int argc, char *argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
