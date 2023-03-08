#include <gtest/gtest.h>

#include "test_geometry.hpp"
#include "test_python.hpp"

int main(int argc, char *argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
