#ifndef __MOCK_Cube_H
#define __MOCK_Cube_H

#include "Eigen/Eigen"
#include "GeometryUtilities.hpp"

namespace UnitTesting {

  class MockCube {
    public:
      static Gedim::GeometryUtilities::Polyhedron CubeWithVertices(const Eigen::Vector3d& origin,
                                                                   const Eigen::Vector3d& lengthVector,
                                                                   const Eigen::Vector3d& heightVector,
                                                                   const Eigen::Vector3d& widthVector)
      {
        Gedim::GeometryUtilities::Polyhedron cube;

        // create vertices
        cube.Vertices.setZero(3, 8);
        cube.Vertices.col(0)<< origin;
        cube.Vertices.col(1)<< origin + lengthVector;
        cube.Vertices.col(2)<< origin + lengthVector + widthVector;
        cube.Vertices.col(3)<< origin + widthVector;
        cube.Vertices.col(4)<< origin + heightVector;
        cube.Vertices.col(5)<< origin + heightVector + lengthVector;
        cube.Vertices.col(6)<< origin + heightVector + lengthVector + widthVector;
        cube.Vertices.col(7)<< origin + heightVector + widthVector;

        // create edges
        cube.Edges.setZero(2, 12);
        cube.Edges.col(0)<< 0, 1;
        cube.Edges.col(1)<< 1, 2;
        cube.Edges.col(2)<< 2, 3;
        cube.Edges.col(3)<< 3, 0;
        cube.Edges.col(4)<< 4, 5;
        cube.Edges.col(5)<< 5, 6;
        cube.Edges.col(6)<< 6, 7;
        cube.Edges.col(7)<< 7, 4;
        cube.Edges.col(8)<< 0, 4;
        cube.Edges.col(9)<< 1, 5;
        cube.Edges.col(10)<< 2, 6;
        cube.Edges.col(11)<< 3, 7;

        // create faces
        cube.Faces.resize(6, Eigen::MatrixXi::Zero(2, 4));
        cube.Faces[0].row(0)<< 0, 1, 2, 3;
        cube.Faces[0].row(1)<< 0, 1, 2, 3;

        cube.Faces[1].row(0)<< 4, 5, 6, 7;
        cube.Faces[1].row(1)<< 4, 5, 6, 7;

        cube.Faces[2].row(0)<< 0, 3, 7, 4;
        cube.Faces[2].row(1)<< 3, 11, 7, 8;

        cube.Faces[3].row(0)<< 1, 2, 6, 5;
        cube.Faces[3].row(1)<< 1, 10, 5, 9;

        cube.Faces[4].row(0)<< 0, 1, 5, 4;
        cube.Faces[4].row(1)<< 0, 9, 4, 8;

        cube.Faces[5].row(0)<< 3, 2, 6, 7;
        cube.Faces[5].row(1)<< 2, 10, 6, 11;

        return cube;
      }
  };
}

#endif // __MOCK_Cube_H
