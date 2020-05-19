#pragma once
#include "tools.h"
#include "connect.h"
#include "globals.h"
#include "Functionalities.h"
using namespace std;

void MatMulCSF2D(int32_t i, int32_t j, int32_t k, vector< vector<myType> >& A, vector< vector<myType> >& B, vector< vector<myType> >& C, int32_t consSF);

void ArgMax1(int32_t outArrS1, int32_t inArrS1, int32_t inArrS2, vector< vector<myType> >& inArr, int32_t dim, vector<myType>& outArr);
void ArgMax3(int32_t outs1, int32_t outs2, int32_t outs3, 
			   int32_t ins1, int32_t ins2, int32_t ins3, int32_t ins4,
			   vector< vector< vector< vector<myType> > > >& inArr, int32_t dim, vector< vector< vector<myType> > >& outArr);


void Relu2(int32_t s1, int32_t s2, vector< vector<myType> >& inArr, vector< vector<myType> >& outArr);
void Relu4(int32_t s1, int32_t s2, int32_t s3, int32_t s4, vector< vector< vector< vector<myType> > > >& inArr, vector< vector< vector< vector<myType> > > >& outArr);

void MaxPool44(int32_t N, int32_t H, int32_t W, int32_t C, 
			  int32_t ksizeH, int32_t ksizeW,
			  int32_t zPadH, int32_t zPadW, 
			  int32_t strideH, int32_t strideW,
			  int32_t N1, int32_t imgH, int32_t imgW, int32_t C1,
			  vector< vector< vector< vector<myType> > > >& inArr, 
			  vector< vector< vector< vector<myType> > > >& outArr);

void MaxPool4(int32_t N, int32_t H, int32_t W, int32_t C, 
			  int32_t imgH, int32_t imgW,
			  int32_t ksizeH, int32_t ksizeW,
			  int32_t zPadH, int32_t zPadW, 
			  int32_t strideH, int32_t strideW,
			  vector< vector< vector< vector<myType> > > >& inArr, 
			  vector< vector< vector< vector<myType> > > >& outArr);

void ElemWiseMul2(int32_t s1, int32_t s2, vector< vector<myType> >& arr1, vector< vector<myType> >& arr2, vector< vector<myType> >& outArr, int32_t shrout);
void ElemWiseMul4(int32_t s1, int32_t s2, int32_t s3, int32_t s4, 
				  vector< vector< vector< vector<myType> > > >& arr1, 
				  vector< vector< vector< vector<myType> > > >& arr2, 
				  vector< vector< vector< vector<myType> > > >& outArr, 
				  int32_t shrout);

void AvgPool44(int32_t N, int32_t H, int32_t W, int32_t C, 
			  int32_t ksizeH, int32_t ksizeW,
			  int32_t zPadH, int32_t zPadW, 
			  int32_t strideH, int32_t strideW,
			  int32_t N1, int32_t imgH, int32_t imgW, int32_t C1,
			  vector< vector< vector< vector<myType> > > >& inArr, 
			  vector< vector< vector< vector<myType> > > >& outArr);

void ScalarMul4(int32_t s1, int32_t s2, int32_t s3, int32_t s4, 
			   int64_t scalar, 
			   vector< vector< vector< vector<myType> > > >& inputArr, 
			   vector< vector< vector< vector<myType> > > >& outputArr, 
			   int64_t consSF);

void ScalarMul2(int32_t s1, int32_t s2,  
			   int64_t scalar, 
			   vector< vector<myType> >& inputArr, 
			   vector< vector<myType> >& outputArr, 
			   int64_t consSF);


// Temporary implementation of fused batch norm 
// TODO 
void TempFusedBatchNorm4411(int32_t s1, int32_t s2, int32_t s3, int32_t s4, 
							vector< vector< vector< vector<myType> > > >& inArr, 
							int32_t vecS1, 
							vector<myType>& multArr, 
							vector<myType>& biasArr, 
							vector< vector< vector< vector<myType> > > >& outputArr);

#ifdef CONV_OPTI
void Conv2DCSF(int32_t N, int32_t H, int32_t W, int32_t CI, int32_t FH, int32_t FW, int32_t CO, int32_t zPadH, int32_t zPadW, int32_t strideH, int32_t strideW, vector< vector< vector< vector<myType> > > >& inputArr, vector< vector< vector< vector<myType> > > >& filterArr, vector< vector< vector< vector<myType> > > >& outArr, int64_t consSF);
#endif
