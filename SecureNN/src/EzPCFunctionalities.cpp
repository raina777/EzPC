#pragma once
#include "EzPCFunctionalities.h"
#include <iostream>
#include <cassert>
#include <chrono>
#include <thread>
extern CommunicationObject commObject;

using namespace std::chrono;

/*
 * NOTE : Please keep in mind that vector should be passed by references to avoid runtime costs.
 */

void MatMulCSF2D(int32_t i, int32_t j, int32_t k, vector< vector<myType> >& A, vector< vector<myType> >& B, vector< vector<myType> >& C, int32_t consSF){
	log_print("EzPCFunctionalities : Starting MatMulCSF2D ... ");

#if (LOG_LAYERWISE)	
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	uint64_t bytesSent = commObject.totalDataSent;
	uint64_t bytesReceived = commObject.totalDataReceived;
#endif

	vector<myType> X(i*j);
	vector<myType> Y(j*k);
	vector<myType> Z(i*k);
	for (int ii=0; ii<i; ii++){
		for (int jj=0; jj<j; jj++){
			X[ii*j + jj] = A[ii][jj]; //Each row is of size j
		}
	}
	for (int ii=0; ii<j; ii++){
		for (int jj=0; jj<k; jj++){
			Y[ii*k + jj] = B[ii][jj]; //Each row is of size k
		}
	}
	funcMatMulMPC(X, Y, Z, i, j, k, 0, 0);
	for (int ii=0; ii<i; ii++){
		for (int jj=0; jj<k; jj++){
			C[ii][jj] = Z[ii*k + jj]; //Each row is of size k
		}
	}

#if (LOG_LAYERWISE)
	high_resolution_clock::time_point t2 = high_resolution_clock::now();	
	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
	auto tt = time_span.count();
	commObject.timeMatmul[1] += tt;
	commObject.dataMatmul[0] += (commObject.totalDataSent - bytesSent);
	commObject.dataMatmul[1] += (commObject.totalDataReceived - bytesReceived);
#endif
}

void empty_func1(){
	return;
}

void funcMaxMPCParallelPrototype(vector<myType> &a, vector<myType> &max, vector<myType> &maxIndex, 
							size_t rows, size_t columns, bool dec)
{
	//funcMaxMPC(a, max, maxIndex, rows, columns, true);
	//INIT_TIMER;
	//START_TIMER;
	//funcMaxMPC(a, max, maxIndex, rows, columns, true);
	//STOP_TIMER("maxMPC function");
	//Parallelize MaxMPC.
	int leaves_chunk_thread0, leaves_chunk_thread1;
	//assert(rows == 1 && "Not yet implemented parallel MaxPool for rows > 1");
	leaves_chunk_thread0 = columns/2;
	leaves_chunk_thread1 = columns - leaves_chunk_thread0;
	//TODO : Remove these vectors. They are filled with 0 for no reason. SLOW.
	
	vector<myType> max_thread0(rows);
	vector<myType> max_thread1(rows);
	vector<myType> maxIndex_thread0(rows);
	vector<myType> maxIndex_thread1(rows);
	vector<myType> a_thread0(rows*leaves_chunk_thread0);
	vector<myType> a_thread1(rows*leaves_chunk_thread1);
	//STOP_TIMER("argmax: Make vectors");
	//START_TIMER;	
	for(int i=0; i< rows*leaves_chunk_thread0; i++){
		a_thread0[i] = a[i];	
	}
	for(int i=0; i< rows*leaves_chunk_thread1; i++){
		a_thread1[i] = a[leaves_chunk_thread0+i];
	}
	//Input loaded. Create threads and call the functions.
	thread* threads	= new thread[2];
	threads[0] = thread(funcMaxMPC, ref(a_thread0), ref(max_thread0), ref(maxIndex_thread0), rows, leaves_chunk_thread0, true);
	//threads[1] = thread(funcMaxMPC, ref(a_thread1), ref(max_thread1), ref(maxIndex_thread1), rows, leaves_chunk_thread1, true);
	threads[1] = thread(empty_func1);
	for(int i=0; i<2; i++){
		threads[i].join();
	}
	delete[] threads;
	//STOP_TIMER("argmax: Leaf layer (threaded)");
	//cout<<"Done leaf layer of MaxPool using 2 threads\n";
	//START_TIMER;
	//Now call a single threaded maxpool to get the max of what both threads reported.
	vector<myType> a_squashed(2*rows);
	vector<myType> maxIndex_raw(rows);
	//TODO : change this. This is incorrect.
	for(int i=0; i<rows; i++){
		a_squashed[i*2 + 0] = max_thread0[i];
		a_squashed[i*2 + 1] = max_thread1[i];
	}
	funcMaxMPC_last(a_squashed, max, maxIndex_raw, rows, 2, dec);
	//Now we have Shares of selector bit in maxIndex_raw.
	vector<myType> diffIndex(rows);
	for(int i=0; i<rows; i++){
		diffIndex[i] = maxIndex_thread1[i] - maxIndex_thread0[i];
	}	
	funcSelectShares3PC(diffIndex, maxIndex_raw, maxIndex, rows);
	for(int i=0; i<rows; i++){
		maxIndex[i] = maxIndex[i] + maxIndex_thread0[i];
	}
	//STOP_TIMER("maxMPC function");	
}


void ArgMax1(int32_t outArrS1, int32_t inArrS1, int32_t inArrS2, vector< vector<myType> >& inArr, int32_t dim, vector<myType>& outArr){
	//TODO_nishkum : right now dim is getting ignored, but once the const optimization is done in tf compilation, dim should be a constant
	//	Then use dim to find out which axis to take max on
	
	log_print("EzPCFunctionalities : Starting ArgMax1 ... ");
	vector<myType> Arr(inArrS1*inArrS2);
	vector<myType> maxi(outArrS1);
	for(int ii=0;ii<inArrS1;ii++){
		for(int jj=0;jj<inArrS2;jj++){
			Arr[ii*inArrS2 + jj] = inArr[ii][jj]; //Each row is of size inArrS2
		}
	}
	//funcMaxMPC(Arr, maxi, outArr, inArrS1, inArrS2, true);
	funcMaxMPCParallelPrototype(Arr, maxi, outArr, inArrS1, inArrS2, true);
}

void ArgMax3(int32_t outs1, int32_t outs2, int32_t outs3, 
			   int32_t ins1, int32_t ins2, int32_t ins3, int32_t ins4,
			   vector< vector< vector< vector<myType> > > >& inArr, int32_t dim, vector< vector< vector<myType> > >& outArr){

	//TODO_nishkum : right now dim is getting ignored, but once the const optimization is done in tf compilation, dim should be a constant
	//	Then use dim to find out which axis to take max on

	//TODO : This is very temp for ccs for quick argmax
	
	log_print("EzPCFunctionalities : Starting ArgMax3 ... ");
	vector<myType> Arr(ins3*ins4);
	vector<myType> maxi(outs3);
	vector<myType> maxIndex(outs3);
	for(int ii=0;ii<ins3;ii++){
		for(int jj=0;jj<ins4;jj++){
			Arr[ii*ins4 + jj] = inArr[0][0][ii][jj]; //Each row is of size inArrS2
		}
	}
	funcMaxMPC(Arr, maxi, maxIndex, ins3, ins4, true);
	for(int ii=0;ii<outs3;ii++){
		outArr[0][0][ii] = maxIndex[ii];
	}
}

void Relu2(int32_t s1, int32_t s2, vector< vector<myType> >& inArr, vector< vector<myType> >& outArr){
	log_print("EzPCFunctionalities : Starting Relu2 ... ");

#if (LOG_LAYERWISE)
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	uint64_t bytesSent = commObject.totalDataSent;
	uint64_t bytesReceived = commObject.totalDataReceived;
#endif
	vector<myType> finArr(s1*s2, 0);
	vector<myType> foutArr(s1*s2, 0);
	for(int i=0;i<s1;i++){
		for(int j=0;j<s2;j++){
			finArr[i*s2+j] = inArr[i][j];
		}
	}
	funcRELUMPC(finArr, foutArr, s1*s2);
	for(int i=0;i<s1;i++){
		for(int j=0;j<s2;j++){
			outArr[i][j] = foutArr[i*s2+j];
		}
	}
#if (LOG_LAYERWISE)
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	
	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
	auto tt = time_span.count();
	commObject.timeRelu += tt;
	commObject.dataRelu[0] += (commObject.totalDataSent - bytesSent);
	commObject.dataRelu[1] += (commObject.totalDataReceived - bytesReceived);
#endif
}

void Relu4(int32_t s1, int32_t s2, int32_t s3, int32_t s4, vector< vector< vector< vector<myType> > > >& inArr, vector< vector< vector< vector<myType> > > >& outArr){
	log_print("EzPCFunctionalities : Starting Relu4 ... ");

#if (LOG_LAYERWISE)
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	uint64_t bytesSent = commObject.totalDataSent;
	uint64_t bytesReceived = commObject.totalDataReceived;
#endif

	vector<myType> finArr(s1*s2*s3*s4, 0);
	vector<myType> foutArr(s1*s2*s3*s4, 0);
	for(int i=0;i<s1;i++){
		for(int j=0;j<s2;j++){
			for(int k=0;k<s3;k++){
				for(int l=0;l<s4;l++){
					finArr[i*s2*s3*s4 + j*s3*s4 + k*s4 + l] = inArr[i][j][k][l];
				}
			}
		}
	}
	funcRELUMPC(finArr, foutArr, s1*s2*s3*s4);
	for(int i=0;i<s1;i++){
		for(int j=0;j<s2;j++){
			for(int k=0;k<s3;k++){
				for(int l=0;l<s4;l++){
					outArr[i][j][k][l] = foutArr[i*s2*s3*s4 + j*s3*s4 + k*s4 + l];
				}
			}
		}
	}

#if (LOG_LAYERWISE)
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	
	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
	auto tt = time_span.count();
	commObject.timeRelu += tt;
	commObject.dataRelu[0] += (commObject.totalDataSent - bytesSent);
	commObject.dataRelu[1] += (commObject.totalDataReceived - bytesReceived);
#endif
}

void ElemWiseMul2(int32_t s1, int32_t s2, vector< vector<myType> >& arr1, vector< vector<myType> >& arr2, vector< vector<myType> >& outArr, int32_t shrout)
{
	log_print("EzPCFunctionalities : Starting ElemWiseMul2 ... ");
	vector<myType> arr1Vec(s1*s2, 0);
	vector<myType> arr2Vec(s1*s2, 0);
	vector<myType> outArrVec(s1*s2, 0);
	for(int i=0;i<s1;i++){
		for(int j=0;j<s2;j++){
			arr1Vec[i*s2+j] = arr1[i][j];
		}
	}
	for(int i=0;i<s1;i++){
		for(int j=0;j<s2;j++){
			arr2Vec[i*s2+j] = arr2[i][j];
		}
	}
	funcDotProductMPC(arr1Vec,arr2Vec,outArrVec,s1*s2);
	for(int i=0;i<s1;i++){
		for(int j=0;j<s2;j++){
			outArr[i][j] = outArrVec[i*s2+j];
		}
	}
}

void ElemWiseMul4(int32_t s1, int32_t s2, int32_t s3, int32_t s4,
		vector< vector< vector< vector<myType> > > >& arr1,
		vector< vector< vector< vector<myType> > > >& arr2,
		vector< vector< vector< vector<myType> > > >& outArr,
		int32_t shrout)
{
	log_print("EzPCFunctionalities : Starting ElemWiseMul4 ... ");
	vector<myType> arr1Vec(s1*s2*s3*s4, 0);
	vector<myType> arr2Vec(s1*s2*s3*s4, 0);
	vector<myType> outArrVec(s1*s2*s3*s4, 0);
	for(int i=0;i<s1;i++){
		for(int j=0;j<s2;j++){
			for(int k=0;k<s3;k++){
				for(int l=0;l<s4;l++){
					arr1Vec[i*s2*s3*s4 + j*s3*s4 + k*s4 + l] = arr1[i][j][k][l];
				}
			}
		}
	}
	for(int i=0;i<s1;i++){
		for(int j=0;j<s2;j++){
			for(int k=0;k<s3;k++){
				for(int l=0;l<s4;l++){
					arr2Vec[i*s2*s3*s4 + j*s3*s4 + k*s4 + l] = arr2[i][j][k][l];
				}
			}
		}
	}

	funcDotProductMPC(arr1Vec, arr2Vec, outArrVec, s1*s2*s3*s4);
	for(int i=0;i<s1;i++){
		for(int j=0;j<s2;j++){
			for(int k=0;k<s3;k++){
				for(int l=0;l<s4;l++){
					outArr[i][j][k][l] = outArrVec[i*s2*s3*s4 + j*s3*s4 + k*s4 + l];
				}
			}
		}
	}
}

//////////////////////////////////////////////////
// MaxPool
/////////////////////////////////////////////////
// void MaxPool4(int32_t N, int32_t H, int32_t W, int32_t C, vector< vector< vector< vector<myType> > > >& inArr, vector< vector< vector< vector<myType> > > >& outArr){
//	//int64_al[N][2L*H][2L*W][C] inArr, int64_al[N][H][W][C] outArr
//	int32_t rows = N*H*W*C;
//	int32_t cols = 2*2; //TODO_nishkum : Hardcoding stride here

//	vector<myType> reInpArr(rows*cols, 0);
//	vector<myType> maxi(rows, 0);
//	vector<myType> maxiIdx(rows, 0);

//	for(int n=0;n<N;n++){
//		for(int h=0;h<H;h++){
//			for(int w=0;w<W;w++){
//				for(int c=0;c<C;c++){
//					int rowIdx = n*H*W*C + h*W*C + w*C + c;
//					for(int ii=0;ii<2;ii++){
//						for(int jj=0;jj<2;jj++){
//							int colIdx = ii*2 + jj;
//							int curIdx = (rowIdx * 4) + colIdx;
//							reInpArr[curIdx] = inArr[n][2*h + ii][2*w + jj][c];
//						}
//					}
//				}
//			}
//		}
//	}

//	funcMaxMPC(reInpArr, maxi, maxiIdx, rows, cols);
//	for(int n=0;n<N;n++){
//		for(int h=0;h<H;h++){
//			for(int w=0;w<W;w++){
//				for(int c=0;c<C;c++){
//					int rowIdx = n*H*W*C + h*W*C + w*C + c;
//					outArr[n][h][w][c] = maxi[rowIdx];
//				}
//			}
//		}
//	}
// }

void MaxPool4Valid(int32_t N, int32_t H, int32_t W, int32_t C,
		int32_t imgH, int32_t imgW,
		int32_t ksizeH, int32_t ksizeW,
		int32_t zPadH, int32_t zPadW,
		int32_t strideH, int32_t strideW,
		vector< vector< vector< vector<myType> > > >& inArr,
		vector< vector< vector< vector<myType> > > >& outArr)
{
	int rows = N*H*W*C;
	int cols = ksizeH*ksizeW;

	vector<myType> reInpArr(rows*cols, 0);
	vector<myType> maxi(rows, 0);
	vector<myType> maxiIdx(rows, 0);

	int rowIdx = 0;
	for(int n=0;n<N;n++){
		for(int c=0;c<C;c++){
			int32_t leftTopCornerH = 0;
			int32_t extremeRightBottomCornerH = imgH - 1;
			while((leftTopCornerH + ksizeH - 1) <= extremeRightBottomCornerH){
				int32_t leftTopCornerW = 0;
				int32_t extremeRightBottomCornerW = imgW - 1;
				while((leftTopCornerW + ksizeW - 1) <= extremeRightBottomCornerW){

					for(int fh=0;fh<ksizeH;fh++){
						for(int fw=0;fw<ksizeW;fw++){
							int32_t colIdx = fh*ksizeW + fw;
							int32_t finalIdx = rowIdx*(ksizeH*ksizeW) + colIdx;

							int32_t curPosH = leftTopCornerH + fh;
							int32_t curPosW = leftTopCornerW + fw;

							uint64_t temp = 0;
							if ((((curPosH < 0) || (curPosH >= imgH)) || ((curPosW < 0) || (curPosW >= imgW)))){
								temp = 0;
							}
							else{
								temp = inArr[n][curPosH][curPosW][c];
							}
							reInpArr[finalIdx] = temp;
						}
					}

					rowIdx += 1;
					leftTopCornerW = leftTopCornerW + strideW;
				}

				leftTopCornerH = leftTopCornerH + strideH;
			}
		}
	}

	funcMaxMPC(reInpArr, maxi, maxiIdx, rows, cols);
	for(int n=0;n<N;n++){
		for(int c=0;c<C;c++){
			for(int h=0;h<H;h++){
				for(int w=0;w<W;w++){
					int iidx = n*C*H*W + c*H*W + h*W + w;
					outArr[n][h][w][c] = maxi[iidx];
				}
			}
		}
	}
}

void MaxPool4Same(int32_t N, int32_t H, int32_t W, int32_t C,
		int32_t imgH, int32_t imgW,
		int32_t ksizeH, int32_t ksizeW,
		int32_t zPadH, int32_t zPadW,
		int32_t strideH, int32_t strideW,
		vector< vector< vector< vector<myType> > > >& inArr,
		vector< vector< vector< vector<myType> > > >& outArr)
{
	int rows = N*H*W*C;
	int cols = ksizeH*ksizeW;

	vector<myType> reInpArr(rows*cols, 0);
	vector<myType> maxi(rows, 0);
	vector<myType> maxiIdx(rows, 0);

	// vector<myType> temp1(32*32,0);
	// vector<myType> temp2(32*32,0);
	// for(int i=0;i<32;i++){
	//	for(int j=0;j<32;j++){
	//		// cout<<inArr[0][i][j][0]<< " ";
	//		temp1[i*32+j] = inArr[0][i][j][0];
	//	}
	//	// cout<<endl;
	// }
	// if (PRIMARY){
	//	funcReconstruct2PC(temp1, 32*32, "reluOutput", &temp2,2);
	//	for(int i=0;i<32;i++){
	//		for(int j=0;j<32;j++){
	//			cout<<temp2[i*32 + j]<<" ";
	//		}
	//		cout<<"\n";
	//	}
	// }

	int rowIdx = 0;
	for(int n=0;n<N;n++){
		for(int c=0;c<C;c++){
			int32_t leftTopCornerH = 0;
			int32_t extremeRightBottomCornerH = imgH - 1;
			while(leftTopCornerH <= extremeRightBottomCornerH){
				int32_t leftTopCornerW = 0;
				int32_t extremeRightBottomCornerW = imgW - 1;
				while(leftTopCornerW <= extremeRightBottomCornerW){

					for(int fh=0;fh<ksizeH;fh++){
						for(int fw=0;fw<ksizeW;fw++){
							int32_t colIdx = fh*ksizeW + fw;
							int32_t finalIdx = rowIdx*(ksizeH*ksizeW) + colIdx;

							int32_t curPosH = leftTopCornerH + fh;
							int32_t curPosW = leftTopCornerW + fw;
							uint64_t temp = 0;
							if ((((curPosH < 0) || (curPosH >= imgH)) || ((curPosW < 0) || (curPosW >= imgW)))){
								temp = 0;
							}
							else{
								temp = inArr[n][curPosH][curPosW][c];
							}
							reInpArr[finalIdx] = temp;
						}
					}

					rowIdx += 1;
					// cout<<" ** "<<rowIdx<<" "<<leftTopCornerH<<" "<<leftTopCornerW<<endl;
					leftTopCornerW = leftTopCornerW + strideW;
				}

				leftTopCornerH = leftTopCornerH + strideH;
			}

			// if (c==0){
			//	vector<myType> temp11(rowIdx*ksizeH*ksizeW,0);
			//	vector<myType> temp12(rowIdx*ksizeH*ksizeW,0);
			//	for(int i=0;i<rowIdx;i++){
			//		for(int j=0;j<ksizeH*ksizeW;j++){
			//			// cout<<reInpArr[i*(ksizeH*ksizeW) + j]<<" ";
			//			temp11[i*(ksizeH*ksizeW) + j] = reInpArr[i*(ksizeH*ksizeW) + j];
			//		}
			//		// cout<<endl;
			//	}

			//	if (PRIMARY){
			//		funcReconstruct2PC(temp11, rowIdx*ksizeH*ksizeW, "maxirearrangeOutput", &temp12,2);
			//		for(int i=0;i<rowIdx;i++){
			//			for(int j=0;j<ksizeH*ksizeW;j++){
			//				// cout<<reInpArr[i*(ksizeH*ksizeW) + j]<<" ";
			//				cout<<temp12[i*(ksizeH*ksizeW) + j]<<" ";
			//			}
			//			cout<<endl;
			//		}
			//	}
			//	onechann = rowIdx;
			// }

		}
	}

	funcMaxMPC(reInpArr, maxi, maxiIdx, rows, cols);
	// vector<myType> tempx1(onechann, 0);
	// for(int i=0;i<onechann;i++){
	//	tempx1[i] = maxi[i];
	// }

	// if (PRIMARY){
	//	vector<myType> tempx2(onechann, 0);
	//	funcReconstruct2PC(tempx1, onechann, "fooi", &tempx2, 2);
	//	cout<<"\n\n\n****************\n\n\n";
	//	for(int i=0;i<onechann;i++){
	//		cout<<tempx2[i]<<endl;
	//	}
	// }


	for(int n=0;n<N;n++){
		for(int c=0;c<C;c++){
			for(int h=0;h<H;h++){
				for(int w=0;w<W;w++){
					int iidx = n*C*H*W + c*H*W + h*W + w;
					outArr[n][h][w][c] = maxi[iidx];
				}
			}
		}
	}
}

void MaxPool44(int32_t N, int32_t H, int32_t W, int32_t C,
		int32_t ksizeH, int32_t ksizeW,
		int32_t zPadH, int32_t zPadW,
		int32_t strideH, int32_t strideW,
		int32_t N1, int32_t imgH, int32_t imgW, int32_t C1,
		vector< vector< vector< vector<myType> > > >& inArr,
		vector< vector< vector< vector<myType> > > >& outArr)
{
	log_print("EzPCFunctionalities : Starting MaxPool44 ... ");
#if (LOG_LAYERWISE)
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	uint64_t bytesSent = commObject.totalDataSent;
	uint64_t bytesReceived = commObject.totalDataReceived;
#endif
	//TODO : Unify and do properly
	if (zPadH == 0){
		MaxPool4Valid(N,H,W,C,imgH,imgW,ksizeH,ksizeW,zPadH,zPadW,strideH,strideW,inArr,outArr);
	}
	else{
		MaxPool4Same(N,H,W,C,imgH,imgW,ksizeH,ksizeW,zPadH,zPadW,strideH,strideW,inArr,outArr);
	};

#if (LOG_LAYERWISE)	
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	
	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
	auto tt = time_span.count();
	commObject.timeMaxpool += tt;
	commObject.dataMaxPool[0] += (commObject.totalDataSent - bytesSent);
	commObject.dataMaxPool[1] += (commObject.totalDataReceived - bytesReceived);
#endif
}

void MaxPool4(int32_t N, int32_t H, int32_t W, int32_t C, 
			  int32_t imgH, int32_t imgW,
			  int32_t ksizeH, int32_t ksizeW,
			  int32_t zPadH, int32_t zPadW, 
			  int32_t strideH, int32_t strideW,
			  vector< vector< vector< vector<myType> > > >& inArr, 
			  vector< vector< vector< vector<myType> > > >& outArr)
{
	//TODO : Unify and do properly
	if (zPadH == 0){
		MaxPool4Valid(N,H,W,C,imgH,imgW,ksizeH,ksizeW,zPadH,zPadW,strideH,strideW,inArr,outArr);
	}
	else{
		MaxPool4Same(N,H,W,C,imgH,imgW,ksizeH,ksizeW,zPadH,zPadW,strideH,strideW,inArr,outArr);
	};
}


//////////////////////////////////////////////////
// AvgPool
/////////////////////////////////////////////////

void AvgPool4Valid(int32_t N, int32_t H, int32_t W, int32_t C,
		int32_t imgH, int32_t imgW,
		int32_t ksizeH, int32_t ksizeW,
		int32_t zPadH, int32_t zPadW,
		int32_t strideH, int32_t strideW,
		vector< vector< vector< vector<myType> > > >& inArr,
		vector< vector< vector< vector<myType> > > >& outArr)
{
	int rows = N*H*W*C;
	vector<myType> filterAvg(rows, 0);

	int rowIdx = 0;
	for(int n=0;n<N;n++){
		for(int c=0;c<C;c++){
			int32_t leftTopCornerH = 0;
			int32_t extremeRightBottomCornerH = imgH - 1;
			while((leftTopCornerH + ksizeH - 1) <= extremeRightBottomCornerH){
				int32_t leftTopCornerW = 0;
				int32_t extremeRightBottomCornerW = imgW - 1;
				while((leftTopCornerW + ksizeW - 1) <= extremeRightBottomCornerW){

					uint64_t curFilterSum = 0;
					for(int fh=0;fh<ksizeH;fh++){
						for(int fw=0;fw<ksizeW;fw++){
							int32_t curPosH = leftTopCornerH + fh;
							int32_t curPosW = leftTopCornerW + fw;

							uint64_t temp = 0;
							if ((((curPosH < 0) || (curPosH >= imgH)) || ((curPosW < 0) || (curPosW >= imgW)))){
								temp = 0;
							}
							else{
								temp = inArr[n][curPosH][curPosW][c];
							}

							curFilterSum += temp;
						}
					}

					//IMP NOTE : The local division should always be signed division.
					//TODO : For now doing local truncation : but this will introduce error
					filterAvg[rowIdx] = static_cast<uint64_t>(static_cast<int64_t>(curFilterSum)/(ksizeH*ksizeW));

					rowIdx += 1;
					leftTopCornerW = leftTopCornerW + strideW;
				}

				leftTopCornerH = leftTopCornerH + strideH;
			}
		}
	}

	for(int n=0;n<N;n++){
		for(int c=0;c<C;c++){
			for(int h=0;h<H;h++){
				for(int w=0;w<W;w++){
					int iidx = n*C*H*W + c*H*W + h*W + w;
					outArr[n][h][w][c] = filterAvg[iidx];
				}
			}
		}
	}
}

void AvgPool4Same(int32_t N, int32_t H, int32_t W, int32_t C,
		int32_t imgH, int32_t imgW,
		int32_t ksizeH, int32_t ksizeW,
		int32_t zPadH, int32_t zPadW,
		int32_t strideH, int32_t strideW,
		vector< vector< vector< vector<myType> > > >& inArr,
		vector< vector< vector< vector<myType> > > >& outArr)
{
	int rows = N*H*W*C;
	vector<myType> filterAvg(rows, 0);

	// vector<myType> temp1(32*32,0);
	// vector<myType> temp2(32*32,0);
	// for(int i=0;i<32;i++){
	//	for(int j=0;j<32;j++){
	//		// cout<<inArr[0][i][j][0]<< " ";
	//		temp1[i*32+j] = inArr[0][i][j][0];
	//	}
	//	// cout<<endl;
	// }
	// if (PRIMARY){
	//	funcReconstruct2PC(temp1, 32*32, "reluOutput", &temp2,2);
	//	for(int i=0;i<32;i++){
	//		for(int j=0;j<32;j++){
	//			cout<<temp2[i*32 + j]<<" ";
	//		}
	//		cout<<"\n";
	//	}
	// }

	int rowIdx = 0;
	for(int n=0;n<N;n++){
		for(int c=0;c<C;c++){
			int32_t leftTopCornerH = 0;
			int32_t extremeRightBottomCornerH = imgH - 1;
			while(leftTopCornerH <= extremeRightBottomCornerH){
				int32_t leftTopCornerW = 0;
				int32_t extremeRightBottomCornerW = imgW - 1;
				while(leftTopCornerW <= extremeRightBottomCornerW){

					uint64_t curFilterSum = 0;
					for(int fh=0;fh<ksizeH;fh++){
						for(int fw=0;fw<ksizeW;fw++){
							int32_t curPosH = leftTopCornerH + fh;
							int32_t curPosW = leftTopCornerW + fw;
							uint64_t temp = 0;
							if ((((curPosH < 0) || (curPosH >= imgH)) || ((curPosW < 0) || (curPosW >= imgW)))){
								temp = 0;
							}
							else{
								temp = inArr[n][curPosH][curPosW][c];
							}
							curFilterSum += temp;
						}
					}

					//IMP NOTE : The local division should always be signed division.
					//TODO : For now doing local truncation : but this will introduce error
					filterAvg[rowIdx] = static_cast<uint64_t>(static_cast<int64_t>(curFilterSum)/(ksizeH*ksizeW));

					rowIdx += 1;
					// cout<<" ** "<<rowIdx<<" "<<leftTopCornerH<<" "<<leftTopCornerW<<endl;
					leftTopCornerW = leftTopCornerW + strideW;
				}

				leftTopCornerH = leftTopCornerH + strideH;
			}

			// if (c==0){
			//	vector<myType> temp11(rowIdx*ksizeH*ksizeW,0);
			//	vector<myType> temp12(rowIdx*ksizeH*ksizeW,0);
			//	for(int i=0;i<rowIdx;i++){
			//		for(int j=0;j<ksizeH*ksizeW;j++){
			//			// cout<<reInpArr[i*(ksizeH*ksizeW) + j]<<" ";
			//			temp11[i*(ksizeH*ksizeW) + j] = reInpArr[i*(ksizeH*ksizeW) + j];
			//		}
			//		// cout<<endl;
			//	}

			//	if (PRIMARY){
			//		funcReconstruct2PC(temp11, rowIdx*ksizeH*ksizeW, "maxirearrangeOutput", &temp12,2);
			//		for(int i=0;i<rowIdx;i++){
			//			for(int j=0;j<ksizeH*ksizeW;j++){
			//				// cout<<reInpArr[i*(ksizeH*ksizeW) + j]<<" ";
			//				cout<<temp12[i*(ksizeH*ksizeW) + j]<<" ";
			//			}
			//			cout<<endl;
			//		}
			//	}
			//	onechann = rowIdx;
			// }

		}
	}

	// vector<myType> tempx1(onechann, 0);
	// for(int i=0;i<onechann;i++){
	//	tempx1[i] = maxi[i];
	// }

	// if (PRIMARY){
	//	vector<myType> tempx2(onechann, 0);
	//	funcReconstruct2PC(tempx1, onechann, "fooi", &tempx2, 2);
	//	cout<<"\n\n\n****************\n\n\n";
	//	for(int i=0;i<onechann;i++){
	//		cout<<tempx2[i]<<endl;
	//	}
	// }


	for(int n=0;n<N;n++){
		for(int c=0;c<C;c++){
			for(int h=0;h<H;h++){
				for(int w=0;w<W;w++){
					int iidx = n*C*H*W + c*H*W + h*W + w;
					outArr[n][h][w][c] = filterAvg[iidx];
				}
			}
		}
	}
}

void AvgPool44(int32_t N, int32_t H, int32_t W, int32_t C,
		int32_t ksizeH, int32_t ksizeW,
		int32_t zPadH, int32_t zPadW,
		int32_t strideH, int32_t strideW,
		int32_t N1, int32_t imgH, int32_t imgW, int32_t C1,
		vector< vector< vector< vector<myType> > > >& inArr,
		vector< vector< vector< vector<myType> > > >& outArr)
{
	log_print("EzPCFunctionalities : Starting AvgPool44 ... ");
#if (LOG_LAYERWISE)
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	uint64_t bytesSent = commObject.totalDataSent;
	uint64_t bytesReceived = commObject.totalDataReceived;
#endif

	//TODO : Unify and do properly
	if (zPadH == 0){
		AvgPool4Valid(N,H,W,C,imgH,imgW,ksizeH,ksizeW,zPadH,zPadW,strideH,strideW,inArr,outArr);
	}
	else{
		AvgPool4Same(N,H,W,C,imgH,imgW,ksizeH,ksizeW,zPadH,zPadW,strideH,strideW,inArr,outArr);
	};
	
#if (LOG_LAYERWISE)
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	
	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
	auto tt = time_span.count();
	commObject.timeAvgPool += tt;
	commObject.dataAvgPool[0] += (commObject.totalDataSent - bytesSent);
	commObject.dataAvgPool[1] += (commObject.totalDataReceived - bytesReceived);
#endif
}

void ScalarMul4(int32_t s1, int32_t s2, int32_t s3, int32_t s4,
		int64_t scalar,
		vector< vector< vector< vector<myType> > > >& inputArr,
		vector< vector< vector< vector<myType> > > >& outputArr,
		int64_t consSF)
{
	log_print("EzPCFunctionalities : Starting ScalarMul4 ... ");
	for(int i=0;i<s1;i++){
		for(int j=0;j<s2;j++){
			for(int k=0;k<s3;k++){
				for(int l=0;l<s4;l++){
					outputArr[i][j][k][l] = inputArr[i][j][k][l]*scalar;
				}
			}
		}
	}
}

void ScalarMul2(int32_t s1, int32_t s2,
		int64_t scalar,
		vector< vector<myType> >& inputArr,
		vector< vector<myType> >& outputArr,
		int64_t consSF)
{
	log_print("EzPCFunctionalities : Starting ScalarMul2 ... ");
	for(int i=0;i<s1;i++){
		for(int j=0;j<s2;j++){
			outputArr[i][j] = inputArr[i][j]*scalar;
		}
	}
}

// Temporary implementation of fused batch norm
// TODO
void TempFusedBatchNorm4411(int32_t s1, int32_t s2, int32_t s3, int32_t s4,
		vector< vector< vector< vector<myType> > > >& inArr,
		int32_t vecS1,
		vector<myType>& multArr,
		vector<myType>& biasArr,
		vector< vector< vector< vector<myType> > > >& outputArr)
{
	log_print("EzPCFunctionalities : Starting TempFusedBatchNorm4411 ... ");
#if (LOG_LAYERWISE)
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	uint64_t bytesSent = commObject.totalDataSent;
	uint64_t bytesReceived = commObject.totalDataReceived;
#endif

	vector<myType> inArrVec(s1*s2*s3*s4, 0);
	vector<myType> multArrVec(s1*s2*s3*s4, 0);
	vector<myType> outArrVec(s1*s2*s3*s4, 0);
	for(int i=0;i<s1;i++){
		for(int j=0;j<s2;j++){
			for(int k=0;k<s3;k++){
				for(int l=0;l<s4;l++){
					inArrVec[i*s2*s3*s4 + j*s3*s4 + k*s4 + l] = inArr[i][j][k][l];
				}
			}
		}
	}
	for(int i=0;i<s1;i++){
		for(int j=0;j<s2;j++){
			for(int k=0;k<s3;k++){
				for(int l=0;l<s4;l++){
					multArrVec[i*s2*s3*s4 + j*s3*s4 + k*s4 + l] = multArr[l];
				}
			}
		}
	}

	funcDotProductMPC(inArrVec, multArrVec, outArrVec, s1*s2*s3*s4);
	for(int i=0;i<s1;i++){
		for(int j=0;j<s2;j++){
			for(int k=0;k<s3;k++){
				for(int l=0;l<s4;l++){
					outputArr[i][j][k][l] = outArrVec[i*s2*s3*s4 + j*s3*s4 + k*s4 + l] + biasArr[l];
				}
			}
		}
	}

#if (LOG_LAYERWISE)
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	
	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
	auto tt = time_span.count();
	commObject.timeBN += tt;
	commObject.dataBN[0] += (commObject.totalDataSent - bytesSent);
	commObject.dataBN[1] += (commObject.totalDataReceived - bytesReceived);
#endif
}

#ifdef CONV_OPTI
void Conv2DCSF(int32_t N, int32_t H, int32_t W, int32_t CI, int32_t FH, int32_t FW, int32_t CO, int32_t zPadH, int32_t zPadW, int32_t strideH, int32_t strideW,
		vector< vector< vector< vector<myType> > > >& inputArr,
		vector< vector< vector< vector<myType> > > >& filterArr,
		vector< vector< vector< vector<myType> > > >& outArr,
		int64_t consSF)
{
	log_print("EzPCFunctionalities : Starting Conv2DCSF ... ");
	funcConv2DCSF(N, H, W, CI, FH, FW, CO, zPadH, zPadW, strideH, strideW, inputArr, filterArr, outArr, consSF);
}
#endif

