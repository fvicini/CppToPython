VectorXd localNumericVelocity = VectorXd::Zero(cellNumLocalsVelocity);
            unsigned int id = 0;
            for(unsigned int i = 0; i < cellNumLocalsVelocityX; ++i)
            {
                const unsigned int globalDofI = dofManagerVelocityX.GlobalIndex(c, i);
                if(globalDofI >= numGlobalDofVelocityX)
                    localNumericVelocity[id++] = previousIteration[globalDofI - numGlobalDofVelocityX + numGlobalDofVelocity
                                                                   + numGlobalDofPressure + 1 + cellNumPartialOffsetX];
                else
                    localNumericVelocity[id++] = previousIteration[globalDofI + cellNumGlobalOffsetX];
            }
            for(unsigned int i = 0; i < cellNumLocalsVelocityY; ++i)
            {
                const unsigned int globalDofI = dofManagerVelocityY.GlobalIndex(c, i);
                if(globalDofI >= numGlobalDofVelocityY)
                    localNumericVelocity[id++] = previousIteration[globalDofI - numGlobalDofVelocityY + numGlobalDofVelocity
                                                                   + numGlobalDofPressure + 1 + cellNumPartialOffsetY];
                else
                    localNumericVelocity[id++] = previousIteration[globalDofI + cellNumGlobalOffsetY];
            }
            for(unsigned int i = 0; i < cellNumLocalsVelocityInternal; ++i)
            {
                const unsigned int globalDofI = dofManagerVelocityInternal.GlobalIndex(c, i);
                localNumericVelocity[id++] = previousIteration[globalDofI + cellNumGlobalOffsetInternal];
            }

            vector<VectorXd> previousIterationValues(2);
            for(unsigned int d = 0; d < 2; d++)
                previousIterationValues[d] = basisFunctionValues[d] * localNumericVelocity;

            vector<VectorXd> previousIterationDerivativesValues(4);
            for(unsigned int d = 0; d < 4; d++)
                previousIterationDerivativesValues[d] = basisFunctionDerivativeValues[d] * localNumericVelocity;

            MatrixXd cellCMatrix1 = MatrixXd::Zero(cellNumLocalsVelocity, cellNumLocalsVelocity);
            MatrixXd cellCMatrix2 = MatrixXd::Zero(cellNumLocalsVelocity, cellNumLocalsVelocity);

            for(unsigned int d1 = 0; d1 < 2; d1++)
            {
                for(unsigned int d2 = 0; d2 < 2; d2++)
                    cellCMatrix1 += basisFunctionValues[d1].transpose()
                                    * internalQuadratureWeights2D.cwiseProduct(previousIterationValues[d2]).asDiagonal()
                                    * basisFunctionDerivativeValues[2 * d1 + d2];
            }

            for(unsigned int d1 = 0; d1 < 2; d1++)
            {
                for(unsigned int d2 = 0; d2 < 2; d2++)
                    cellCMatrix1 += basisFunctionValues[d1].transpose()
                                    * internalQuadratureWeights2D.cwiseProduct(previousIterationDerivativesValues[2 * d1 + d2]).asDiagonal()
                                    * basisFunctionValues[d2];
            }

            for(unsigned int d1 = 0; d1 < 2; d1++)
            {
                for(unsigned int d2 = 0; d2 < 2; d2++)
                    cellCMatrix2 += basisFunctionDerivativeValues[2 * d1 + d2].transpose()
                                    * internalQuadratureWeights2D.cwiseProduct(previousIterationValues[d2]).asDiagonal()
                                    * basisFunctionValues[d1];
            }

            for(unsigned int d1 = 0; d1 < 2; d1++)
            {
                for(unsigned int d2 = 0; d2 < 2; d2++)
                    cellCMatrix2 += basisFunctionDerivativeValues[2 * d1 + d2].transpose()
                                    * internalQuadratureWeights2D.cwiseProduct(previousIterationValues[d1]).asDiagonal()
                                    * basisFunctionValues[d2];
            }

            cellMatrixC = 0.5 * (cellCMatrix1 - cellCMatrix2);

            VectorXd cellRightHandSide1 = VectorXd::Zero(cellNumLocalsVelocity);
            VectorXd cellRightHandSide2 = VectorXd::Zero(cellNumLocalsVelocity);

            for(unsigned int d1 = 0; d1 < 2; d1++)
            {
                for(unsigned int d2 = 0; d2 < 2; d2++)
                    cellRightHandSide1 += basisFunctionValues[d1].transpose()
                                          * internalQuadratureWeights2D.asDiagonal()
                                          * previousIterationValues[d2].cwiseProduct(previousIterationDerivativesValues[2 * d1 + d2]);
            }

            for(unsigned int d1 = 0; d1 < 2; d1++)
            {
                for(unsigned int d2 = 0; d2 < 2; d2++)
                    cellRightHandSide2 += basisFunctionDerivativeValues[2 * d1 + d2].transpose()
                                          * internalQuadratureWeights2D.asDiagonal()
                                          * previousIterationValues[d2].cwiseProduct(previousIterationValues[d1]);
            }

            cellRightHandSideC = 0.5 * (cellRightHandSide1 - cellRightHandSide2);
