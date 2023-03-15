#ifndef __GeDiM4Py_Interface_H
#define __GeDiM4Py_Interface_H

#include <Python.h>
#include <iostream>

#include "GeDiM4Py_Logic.hpp"

// ***************************************************************************
extern "C"
void GedimForPy_Initialize(PyObject* config);
extern "C"
PyObject* GedimForPy_CreateDomainSquare(PyObject* square,
                                        double** coordinates);
extern "C"
PyObject* GedimForPy_Discretize(PyObject* discreteSpace,
                                double** dofsCoordinate,
                                double** strongsCoordinate);
extern "C"
void GedimForPy_AssembleStiffnessMatrix(GedimForPy::GeDiM4Py_Logic::A a,
                                        int* numStiffnessTriplets,
                                        double** stiffnessTriplets,
                                        int* numStiffnessStrongTriplets,
                                        double** stiffnessStrongTriplets);
extern "C"
void GedimForPy_AssembleAdvectionMatrix(GedimForPy::GeDiM4Py_Logic::B b,
                                        int* numAdvectionTriplets,
                                        double** advectionTriplets,
                                        int* numAdvectionStrongTriplets,
                                        double** advectionStrongTriplets);
extern "C"
void GedimForPy_AssembleReactionMatrix(GedimForPy::GeDiM4Py_Logic::C c,
                                       int* numReactionTriplets,
                                       double** reactionTriplets,
                                       int* numReactionStrongTriplets,
                                       double** reactionStrongTriplets);
extern "C"
void GedimForPy_AssembleForcingTerm(GedimForPy::GeDiM4Py_Logic::F f,
                                    int* size,
                                    double** forcingTerm);
extern "C"
void GedimForPy_AssembleStrongSolution(GedimForPy::GeDiM4Py_Logic::Strong g,
                                       int marker,
                                       int* size,
                                       double** solutionStrong);
extern "C"
void GedimForPy_AssembleWeakTerm(GedimForPy::GeDiM4Py_Logic::Weak g,
                                 int marker,
                                 int* size,
                                 double** weakTerm);
extern "C"
void GedimForPy_CholeskySolver(const int nRows,
                               const int numNonZeros,
                               const double* pointerTriplets,
                               const double* pointerF,
                               double** solution);
extern "C"
void GedimForPy_LUSolver(const int nRows,
                         const int numNonZeros,
                         const double* pointerTriplets,
                         const double* pointerF,
                         double** solution);
extern "C"
double GedimForPy_ComputeErrorL2(GedimForPy::GeDiM4Py_Logic::Exact u,
                                 const double* pointerNumericSolution,
                                 const double* pointerStrongSolution);
extern "C"
double GedimForPy_ComputeErrorH1(GedimForPy::GeDiM4Py_Logic::ExactDerivative uDer,
                                 const double* pointerNumericSolution,
                                 const double* pointerStrongSolution);
// ***************************************************************************

namespace GedimForPy
{
  class GeDiM4Py_Interface final
  {
    public:
      static GedimForPy::InterfaceConfiguration InterfaceConfig;
      static GedimForPy::InterfaceData InterfaceData;
      static GedimForPy::Domain2D Domain;
      static GedimForPy::Domain2DMesh Mesh;
      static GedimForPy::DiscreteSpace Space;
      static GedimForPy::DiscreteProblemData ProblemData;

    private:
      template<typename T>
      static std::vector<T> ConvertToArray(PyObject* list);

    public:
      static InterfaceConfiguration ConvertConfiguration(PyObject* config);
      static Domain2D ConvertDomainSquare(InterfaceDataDAO& gedimData,
                                          PyObject* square);
      static PyObject* ConvertMesh(const Gedim::IMeshDAO& mesh,
                                   double*& coordinates);
      static DiscreteSpace ConvertDiscreteSpace(PyObject* discreteSpace);
      static PyObject* ConvertProblemData(GedimForPy::DiscreteProblemData& problemData,
                                          double*& dofsCoordinate,
                                          double*& strongsCoordinate);

      static void ConvertTriplets(const std::list<Eigen::Triplet<double>>& triplets,
                                  int& numTriplets,
                                  double*& convertedTriplets);
      static Eigen::SparseMatrix<double> ConvertSparseMatrix(const int& nRows,
                                                             const int& nCols,
                                                             const int& numNonZeros,
                                                             const double* pointerTriplets);
      static void ConvertMatrix(const Eigen::MatrixXd& matrix,
                                double*& convertedMatrix);
      static void ConvertArray(const Eigen::VectorXd& array,
                               int& size,
                               double*& convertedArray);
  };
}

#endif
