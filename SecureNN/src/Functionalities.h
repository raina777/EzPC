
#pragma once
#include "tools.h"
#include "connect.h"
#include "globals.h"
using namespace std;

extern void start_time();
extern void start_communication();
extern void end_time(string str);
extern void end_communication(string str);



void funcTruncate2PC(vector<myType> &a, size_t power, size_t size, size_t party_1, size_t party_2);
void funcXORModuloOdd2PC(vector<smallType> &bit, vector<myType> &shares, vector<myType> &output, size_t size); 
void funcReconstruct2PC(const vector<myType> &a, size_t size, string str, vector<myType>* b=NULL, int revealToParties = 2);
void funcReconstructBit2PC(const vector<smallType> &a, size_t size, string str);
void funcConditionalSet2PC(const vector<myType> &a, const vector<myType> &b, vector<smallType> &c, 
							vector<myType> &u, vector<myType> &v, size_t size);
void funcSecretShareConstant(const vector<myType> &cons, vector<myType> &curShare, size_t size);
void funcMatMulMPC(const vector<myType> &a, const vector<myType> &b, vector<myType> &c, 
				size_t rows, size_t common_dim, size_t columns,
			 	size_t transpose_a, size_t transpose_b, bool areInputsScaled = true, int coming_from_conv = 0);
void funcDotProductMPC(const vector<myType> &a, const vector<myType> &b, 
					   vector<myType> &c, size_t size);
void funcPrivateCompareMPC(const vector<smallType> &share_m, const vector<myType> &r, 
							  const vector<smallType> &beta, vector<smallType> &betaPrime, 
							  size_t size, size_t dim);
void funcShareConvertMPC(vector<myType> &a, size_t size);
void funcComputeMSB4PC(const vector<myType> &a, vector<smallType> &b, size_t size);
void funcComputeMSB3PC(const vector<myType> &a, vector<myType> &b, size_t size);
void funcSelectShares4PC(const vector<myType> &a, const vector<smallType> &b, vector<myType> &c, size_t size);
void funcSelectShares3PC(const vector<myType> &a, const vector<myType> &b, vector<myType> &c, size_t size);
void funcRELUPrime4PC(const vector<myType> &a, vector<smallType> &b, size_t size);
void funcRELUPrime3PC(const vector<myType> &a, vector<myType> &b, size_t size);
void funcRELUMPC(const vector<myType> &a, vector<myType> &b, size_t size);
void funcDivisionMPC(const vector<myType> &a, const vector<myType> &b, vector<myType> &quotient, 
						size_t size);
void funcMaxMPC_last(vector<myType> &a, vector<myType> &max, vector<myType> &maxIndex, 
						size_t rows, size_t columns, bool computeMaxIndex = false);
void funcMaxMPC(vector<myType> &a, vector<myType> &max, vector<myType> &maxIndex, 
						size_t rows, size_t columns, bool computeMaxIndex = false);
void funcMaxIndexMPC(vector<myType> &a, const vector<myType> &maxIndex, 
						size_t rows, size_t columns);
void aggregateCommunication();

//Wrapper calls which take integers
myType funcMult(myType a, myType b);
myType funcReluPrime(myType a);
myType funcDiv(myType a, myType b);
myType funcSSCons(myType a);
myType funcReconstruct2PCCons(myType a, int revealToParties); //revealToParties is a bitmask as to which parties should learn
																	//	the reconstructed value.


void funcConv2DCSF(int32_t N, int32_t H, int32_t W, int32_t CI, int32_t FH, int32_t FW, int32_t CO, int32_t zPadH, int32_t zPadW, int32_t strideH, int32_t strideW,
		vector< vector< vector< vector<myType> > > >& inputArr,
		vector< vector< vector< vector<myType> > > >& filterArr,
		vector< vector< vector< vector<myType> > > >& outArr,
		int64_t consSF);


//Debug
void debugDotProd();
void debugComputeMSB();
void debugPC();
void debugDivision();
void debugMax();
void debugSS();
void debugMatMul();
void debugReLUPrime();
void debugMaxIndex();

//Test
void testMatMul(size_t rows, size_t common_dim, size_t columns, size_t iter);
void testRelu(size_t r, size_t c, size_t iter);
void testReluPrime(size_t r, size_t c, size_t iter);
void testMaxPool(size_t p_range, size_t q_range, size_t px, size_t py, size_t D, size_t iter);
void testMaxPoolDerivative(size_t p_range, size_t q_range, size_t px, size_t py, size_t D, size_t iter);
void testComputeMSB(size_t size);
void testShareConvert(size_t size);

