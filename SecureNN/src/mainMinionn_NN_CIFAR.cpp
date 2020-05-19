#include "globals.h"
#ifdef F_MINIONN_CIFAR

#include<vector>
#include<math.h>
#include<cstdlib>
#include<iostream>
#include<fstream>
#include "EzPCFunctionalities.h"

using namespace std;

uint32_t public_lrshift(uint32_t x, uint32_t y){
return (x >> y);
}

int32_t public_lrshift(int32_t x, uint32_t y){
return ((int32_t)(((uint32_t)x) >> y));
}

uint64_t public_lrshift(uint64_t x, uint64_t y){
return (x >> y);
}

int64_t public_lrshift(int64_t x, uint64_t y){
return ((int64_t)(((uint64_t)x) >> y));
}

template<typename T>
vector<T> make_vector(size_t size) {
return std::vector<T>(size);
}

template <typename T, typename... Args>
auto make_vector(size_t first, Args... sizes)
{
auto inner = make_vector<T>(sizes...);
return vector<decltype(inner)>(first, inner);
}

template<typename T>
ostream& operator<< (ostream &os, const vector<T> &v)
{
for(auto it = v.begin (); it != v.end (); ++it) {
os << *it << endl;
}
return os;
}


#include "ezpc.h"

extern int partyNum;
vector<uint64_t*> toFreeMemoryLaterArr;
int NUM_OF_PARTIES;

AESObject* aes_common;
AESObject* aes_indep;
AESObject* aes_a_1;
AESObject* aes_a_2;
AESObject* aes_b_1;
AESObject* aes_b_2;
AESObject* aes_c_1;
AESObject* aes_share_conv_bit_shares_p0_p2;
AESObject* aes_share_conv_bit_shares_p1_p2;
AESObject* aes_share_conv_shares_mod_odd_p0_p2;
AESObject* aes_share_conv_shares_mod_odd_p1_p2;
AESObject* aes_comp_msb_shares_lsb_p0_p2;
AESObject* aes_comp_msb_shares_lsb_p1_p2;
AESObject* aes_comp_msb_shares_bit_vec_p0_p2;
AESObject* aes_comp_msb_shares_bit_vec_p1_p2;
AESObject* aes_conv_opti_a_1;
AESObject* aes_conv_opti_a_2;
AESObject* aes_conv_opti_b_1;
AESObject* aes_conv_opti_b_2;
AESObject* aes_conv_opti_c_1;
ParallelAESObject* aes_parallel;

output_queue out_q;











void MatAddBroadCast2(int32_t s1, int32_t s2, auto& A, auto& B, auto& outArr){
for (uint32_t i1 =  (int32_t)0; i1 < s1; i1++){
for (uint32_t i2 =  (int32_t)0; i2 < s2; i2++){
outArr[i1][i2] = A[i1][i2]+B[i2];
}
}
}

void MatAdd2(int32_t s1, int32_t s2, auto& A, auto& B, auto& outArr){
for (uint32_t i1 =  (int32_t)0; i1 < s1; i1++){
for (uint32_t i2 =  (int32_t)0; i2 < s2; i2++){
outArr[i1][i2] = A[i1][i2]+B[i1][i2];
}
}
}

void MatAddBroadCast4(int32_t s1, int32_t s2, int32_t s3, int32_t s4, auto& A, auto& B, auto& outArr){
for (uint32_t i1 =  (int32_t)0; i1 < s1; i1++){
for (uint32_t i2 =  (int32_t)0; i2 < s2; i2++){
for (uint32_t i3 =  (int32_t)0; i3 < s3; i3++){
for (uint32_t i4 =  (int32_t)0; i4 < s4; i4++){
outArr[i1][i2][i3][i4] = A[i1][i2][i3][i4]+B[i4];
}
}
}
}
}

void MatAdd4(int32_t s1, int32_t s2, int32_t s3, int32_t s4, auto& A, auto& B, auto& outArr){
for (uint32_t i1 =  (int32_t)0; i1 < s1; i1++){
for (uint32_t i2 =  (int32_t)0; i2 < s2; i2++){
for (uint32_t i3 =  (int32_t)0; i3 < s3; i3++){
for (uint32_t i4 =  (int32_t)0; i4 < s4; i4++){
outArr[i1][i2][i3][i4] = A[i1][i2][i3][i4]+B[i1][i2][i3][i4];
}
}
}
}
}

void CreateTensor1(int32_t s1, int64_t val, auto& arr){
for (uint32_t i1 =  (int32_t)0; i1 < s1; i1++){
arr[i1] = val;
}
}

void CreateTensor2(int32_t s1, int32_t s2, int64_t val, auto& arr){
for (uint32_t i1 =  (int32_t)0; i1 < s1; i1++){
for (uint32_t i2 =  (int32_t)0; i2 < s2; i2++){
arr[i1][i2] = val;
}
}
}

void CreateTensor4(int32_t s1, int32_t s2, int32_t s3, int32_t s4, int64_t val, auto& arr){
for (uint32_t i1 =  (int32_t)0; i1 < s1; i1++){
for (uint32_t i2 =  (int32_t)0; i2 < s2; i2++){
for (uint32_t i3 =  (int32_t)0; i3 < s3; i3++){
for (uint32_t i4 =  (int32_t)0; i4 < s4; i4++){
arr[i1][i2][i3][i4] = val;
}
}
}
}
}

void CopyTensor1(int32_t s1, auto& targetArr, auto& fromArr, auto& ignore){
for (uint32_t i1 =  (int32_t)0; i1 < s1; i1++){
targetArr[i1] = fromArr[i1];
}
}

void CopyTensor2(int32_t s1, int32_t s2, auto& targetArr, auto& fromArr, auto& ignore){
for (uint32_t i1 =  (int32_t)0; i1 < s1; i1++){
for (uint32_t i2 =  (int32_t)0; i2 < s2; i2++){
targetArr[i1][i2] = fromArr[i1][i2];
}
}
}

void CopyTensor4(int32_t s1, int32_t s2, int32_t s3, int32_t s4, auto& targetArr, auto& fromArr, auto& ignore){
for (uint32_t i1 =  (int32_t)0; i1 < s1; i1++){
for (uint32_t i2 =  (int32_t)0; i2 < s2; i2++){
for (uint32_t i3 =  (int32_t)0; i3 < s3; i3++){
for (uint32_t i4 =  (int32_t)0; i4 < s4; i4++){
targetArr[i1][i2][i3][i4] = fromArr[i1][i2][i3][i4];
}
}
}
}
}

void CreateIdentity11(int32_t s1, auto& fromArr, auto& newArr){
for (uint32_t i1 =  (int32_t)0; i1 < s1; i1++){
newArr[i1] = fromArr[i1];
}
}

void CreateIdentity22(int32_t s1, int32_t s2, auto& fromArr, auto& newArr){
for (uint32_t i1 =  (int32_t)0; i1 < s1; i1++){
for (uint32_t i2 =  (int32_t)0; i2 < s2; i2++){
newArr[i1][i2] = fromArr[i1][i2];
}
}
}

void CreateIdentity44(int32_t s1, int32_t s2, int32_t s3, int32_t s4, auto& fromArr, auto& newArr){
for (uint32_t i1 =  (int32_t)0; i1 < s1; i1++){
for (uint32_t i2 =  (int32_t)0; i2 < s2; i2++){
for (uint32_t i3 =  (int32_t)0; i3 < s3; i3++){
for (uint32_t i4 =  (int32_t)0; i4 < s4; i4++){
newArr[i1][i2][i3][i4] = fromArr[i1][i2][i3][i4];
}
}
}
}
}

void Concat2T444(int32_t s1, int32_t s2, int32_t s3, int32_t s4, int32_t inp1s1, int32_t inp1s2, int32_t inp1s3, int32_t inp1s4, auto& inp1, int32_t inp2s1, int32_t inp2s2, int32_t inp2s3, int32_t inp2s4, auto& inp2, int32_t axis, auto& outp){
for (uint32_t i1 =  (int32_t)0; i1 < s1; i1++){
for (uint32_t i2 =  (int32_t)0; i2 < s2; i2++){
for (uint32_t i3 =  (int32_t)0; i3 < s3; i3++){
for (uint32_t i4 =  (int32_t)0; i4 < s4; i4++){
if ((axis ==  (int32_t)0)) {
if ((i1 < inp1s1)) {
outp[i1][i2][i3][i4] = inp1[i1][i2][i3][i4];
} else {
outp[i1][i2][i3][i4] = inp2[(i1 - inp1s1)][i2][i3][i4];
}
} else {
if ((axis ==  (int32_t)1)) {
if ((i2 < inp1s2)) {
outp[i1][i2][i3][i4] = inp1[i1][i2][i3][i4];
} else {
outp[i1][i2][i3][i4] = inp2[i1][(i2 - inp1s2)][i3][i4];
}
} else {
if ((axis ==  (int32_t)2)) {
if ((i3 < inp1s3)) {
outp[i1][i2][i3][i4] = inp1[i1][i2][i3][i4];
} else {
outp[i1][i2][i3][i4] = inp2[i1][i2][(i3 - inp1s3)][i4];
}
} else {
if ((i4 < inp1s4)) {
outp[i1][i2][i3][i4] = inp1[i1][i2][i3][i4];
} else {
outp[i1][i2][i3][i4] = inp2[i1][i2][i3][(i4 - inp1s4)];
}
}
}
}
}
}
}
}
}

void RandomUniform2(int32_t s1, int32_t s2, int64_t dataType, auto& outArr){
for (uint32_t i1 =  (int32_t)0; i1 < s1; i1++){
for (uint32_t i2 =  (int32_t)0; i2 < s2; i2++){
outArr[i1][i2] = funcSSCons( (int64_t)100);
}
}
}

void Conv2DReshapeFilter(int32_t FH, int32_t FW, int32_t CI, int32_t CO, auto& inputArr, auto& outputArr){
for (uint32_t co =  (int32_t)0; co < CO; co++){
for (uint32_t fh =  (int32_t)0; fh < FH; fh++){
for (uint32_t fw =  (int32_t)0; fw < FW; fw++){
for (uint32_t ci =  (int32_t)0; ci < CI; ci++){

int32_t linIdx = ((((fh * FW) * CI) + (fw * CI)) + ci);
outputArr[co][linIdx] = inputArr[fh][fw][ci][co];
}
}
}
}
}

void Conv2DReshapeMatMulOP(int32_t N, int32_t finalH, int32_t finalW, int32_t CO, auto& inputArr, auto& outputArr){
for (uint32_t co =  (int32_t)0; co < CO; co++){
for (uint32_t n =  (int32_t)0; n < N; n++){
for (uint32_t h =  (int32_t)0; h < finalH; h++){
for (uint32_t w =  (int32_t)0; w < finalW; w++){
outputArr[n][h][w][co] = inputArr[co][((((n * finalH) * finalW) + (h * finalW)) + w)];
}
}
}
}
}

void Conv2DReshapeInput(int32_t N, int32_t H, int32_t W, int32_t CI, int32_t FH, int32_t FW, int32_t zeroPadH, int32_t zeroPadW, int32_t strideH, int32_t strideW, int32_t RRows, int32_t RCols, auto& inputArr, auto& outputArr){

int32_t linIdxFilterMult =  (int32_t)0;
for (uint32_t n =  (int32_t)0; n < N; n++){

int32_t leftTopCornerH = ( (int32_t)0 - zeroPadH);

int32_t extremeRightBottomCornerH = ((H -  (int32_t)1) + zeroPadH);
while ((((leftTopCornerH + FH) -  (int32_t)1) <= extremeRightBottomCornerH)) {

int32_t leftTopCornerW = ( (int32_t)0 - zeroPadW);

int32_t extremeRightBottomCornerW = ((W -  (int32_t)1) + zeroPadW);
while ((((leftTopCornerW + FW) -  (int32_t)1) <= extremeRightBottomCornerW)) {
for (uint32_t fh =  (int32_t)0; fh < FH; fh++){
for (uint32_t fw =  (int32_t)0; fw < FW; fw++){

int32_t curPosH = (leftTopCornerH + fh);

int32_t curPosW = (leftTopCornerW + fw);

uint64_t val = 0;
for (uint32_t ci =  (int32_t)0; ci < CI; ci++){
if ((((curPosH <  (int32_t)0) || (curPosH >= H)) || ((curPosW <  (int32_t)0) || (curPosW >= W)))) {
val = 0;
} else {
val = inputArr[n][curPosH][curPosW][ci];
}
outputArr[((((fh * FW) * CI) + (fw * CI)) + ci)][linIdxFilterMult] = val;
}
}
}
linIdxFilterMult = (linIdxFilterMult +  (int32_t)1);
leftTopCornerW = (leftTopCornerW + strideW);
}

leftTopCornerH = (leftTopCornerH + strideH);
}

}
}

#ifndef CONV_OPTI
void Conv2DCSF(int32_t N, int32_t H, int32_t W, int32_t CI, int32_t FH, int32_t FW, int32_t CO, int32_t zPadH, int32_t zPadW, int32_t strideH, int32_t strideW, auto& inputArr, auto& filterArr, auto& outArr, int64_t consSF){
int32_t reshapedFilterRows = CO;

int32_t reshapedFilterCols = ((FH * FW) * CI);

int32_t reshapedIPRows = ((FH * FW) * CI);

int32_t newH = ((((H + ( (int32_t)2 * zPadH)) - FH) / strideH) +  (int32_t)1);

int32_t newW = ((((W + ( (int32_t)2 * zPadW)) - FW) / strideW) +  (int32_t)1);

int32_t reshapedIPCols = ((N * newH) * newW);

auto filterReshaped = make_vector<uint64_t>(reshapedFilterRows, reshapedFilterCols);

auto inputReshaped = make_vector<uint64_t>(reshapedIPRows, reshapedIPCols);

auto matmulOP = make_vector<uint64_t>(reshapedFilterRows, reshapedIPCols);
Conv2DReshapeFilter(FH, FW, CI, CO, filterArr, filterReshaped);
Conv2DReshapeInput(N, H, W, CI, FH, FW, zPadH, zPadW, strideH, strideW, reshapedIPRows, reshapedIPCols, inputArr, inputReshaped);
MatMulCSF2D(reshapedFilterRows, reshapedFilterCols, reshapedIPCols, filterReshaped, inputReshaped, matmulOP, consSF);
Conv2DReshapeMatMulOP(N, newH, newW, CO, matmulOP, outArr);
}

#endif

void Conv2DCSFMain(int32_t N, int32_t H, int32_t W, int32_t CI, int32_t FH, int32_t FW, int32_t CO, int32_t zPadH, int32_t zPadW, int32_t strideH, int32_t strideW, auto& inputArr, auto& filterArr, auto& outArr, int64_t consSF)
{
#ifdef CONV_OPTI
	if ((FH>=5) || (FW>=5))
	{
		funcConv2DCSF(N, H, W, CI, FH, FW, CO, zPadH, zPadW, strideH, strideW, inputArr, filterArr, outArr, consSF);
	}
	else
	{
		Conv2DCSF(N, H, W, CI, FH, FW, CO, zPadH, zPadW, strideH, strideW, inputArr, filterArr, outArr, consSF);
	}
#else
	Conv2DCSF(N, H, W, CI, FH, FW, CO, zPadH, zPadW, strideH, strideW, inputArr, filterArr, outArr, consSF);
#endif
}


void Transpose2(int32_t s1, int32_t s2, auto& inArr, auto& outArr){
for (uint32_t i =  (int32_t)0; i < s1; i++){
for (uint32_t j =  (int32_t)0; j < s2; j++){
outArr[i][j] = inArr[j][i];
}
}
}

void Pad442(int32_t s1, int32_t s2, int32_t s3, int32_t s4, int32_t inps1, int32_t inps2, int32_t inps3, int32_t inps4, auto& inpArr, int32_t pads1, int32_t pads2, auto& paddings, auto& outArr){

int32_t lbounds1 = paddings[ (int32_t)0][ (int32_t)0];

int32_t rbounds1excl = (s1 - paddings[ (int32_t)0][ (int32_t)1]);

int32_t lbounds2 = paddings[ (int32_t)1][ (int32_t)0];

int32_t rbounds2excl = (s2 - paddings[ (int32_t)1][ (int32_t)1]);

int32_t lbounds3 = paddings[ (int32_t)2][ (int32_t)0];

int32_t rbounds3excl = (s3 - paddings[ (int32_t)2][ (int32_t)1]);

int32_t lbounds4 = paddings[ (int32_t)3][ (int32_t)0];

int32_t rbounds4excl = (s4 - paddings[ (int32_t)3][ (int32_t)1]);
for (uint32_t i =  (int32_t)0; i < s1; i++){
for (uint32_t j =  (int32_t)0; j < s2; j++){
for (uint32_t k =  (int32_t)0; k < s3; k++){
for (uint32_t l =  (int32_t)0; l < s4; l++){
if (((((((((i >= lbounds1) && (i < rbounds1excl)) && (j >= lbounds2)) && (j < rbounds2excl)) && (k >= lbounds3)) && (k < rbounds3excl)) && (l >= lbounds4)) && (l < rbounds4excl))) {
outArr[i][j][k][l] = inpArr[(i - paddings[ (int32_t)0][ (int32_t)0])][(j - paddings[ (int32_t)1][ (int32_t)0])][(k - paddings[ (int32_t)2][ (int32_t)0])][(l - paddings[ (int32_t)3][ (int32_t)0])];
} else {
outArr[i][j][k][l] = 0;
}
}
}
}
}
}

void Squeeze24(int32_t s1, int32_t s2, int32_t dim1, int32_t dim2, int32_t ins1, int32_t ins2, int32_t ins3, int32_t ins4, auto& inArr, auto& outArr){
for (uint32_t i =  (int32_t)0; i < ins1; i++){
for (uint32_t j =  (int32_t)0; j < ins2; j++){
for (uint32_t k =  (int32_t)0; k < ins3; k++){
for (uint32_t l =  (int32_t)0; l < ins4; l++){

int32_t linIdx = ((((((i * ins2) * ins3) * ins4) + ((j * ins3) * ins4)) + (k * ins4)) + l);

int32_t outIdx1 = (linIdx / s2);

int32_t outIdx2 = (linIdx % s2);
outArr[outIdx1][outIdx2] = inArr[i][j][k][l];
}
}
}
}
}


int main(int argc, char** argv)
{
parseInputs(argc, argv, false);
string whichNetwork = "No Network";

aes_indep = new AESObject(argv[4]);
aes_common = new AESObject(argv[5]);
aes_a_1 = new AESObject("files/keyD");
aes_a_2 = new AESObject("files/keyD");
aes_b_1 = new AESObject("files/keyD");
aes_b_2 = new AESObject("files/keyD");
aes_c_1 = new AESObject("files/keyD");
aes_share_conv_bit_shares_p0_p2 = new AESObject("files/keyD");
aes_share_conv_bit_shares_p1_p2 = new AESObject("files/keyD");
aes_share_conv_shares_mod_odd_p0_p2 = new AESObject("files/keyD");
aes_share_conv_shares_mod_odd_p1_p2 = new AESObject("files/keyD");
aes_comp_msb_shares_lsb_p0_p2 = new AESObject("files/keyD");
aes_comp_msb_shares_lsb_p1_p2 = new AESObject("files/keyD");
aes_comp_msb_shares_bit_vec_p0_p2 = new AESObject("files/keyD");
aes_comp_msb_shares_bit_vec_p1_p2 = new AESObject("files/keyD");
aes_conv_opti_a_1 = new AESObject("files/keyD");
aes_conv_opti_a_2 = new AESObject("files/keyD");
aes_conv_opti_b_1 = new AESObject("files/keyD");
aes_conv_opti_b_2 = new AESObject("files/keyD");
aes_conv_opti_c_1 = new AESObject("files/keyD");
aes_parallel = new ParallelAESObject(argv[5]);

if (!STANDALONE)
{
initializeMPC();
initializeCommunication(argv[3], partyNum);
synchronize(2000000); 
}

if (PARALLEL)
aes_parallel->precompute();

e_role role = partyNum;
// start_m();


auto tmp10 = make_vector<int64_t>( (int32_t)3,  (int32_t)3,  (int32_t)3,  (int32_t)64);

auto tmp11 = make_vector<uint64_t>( (int32_t)3,  (int32_t)3,  (int32_t)3,  (int32_t)64);

auto tmp12 = make_vector<uint64_t>( (int32_t)1,  (int32_t)32,  (int32_t)32,  (int32_t)64);

auto tmp13 = make_vector<uint64_t>( (int32_t)1,  (int32_t)32,  (int32_t)32,  (int32_t)64);

auto tmp14 = make_vector<int64_t>( (int32_t)3,  (int32_t)3,  (int32_t)64,  (int32_t)64);

auto tmp15 = make_vector<uint64_t>( (int32_t)3,  (int32_t)3,  (int32_t)64,  (int32_t)64);

auto tmp16 = make_vector<uint64_t>( (int32_t)1,  (int32_t)32,  (int32_t)32,  (int32_t)64);

auto tmp17 = make_vector<uint64_t>( (int32_t)1,  (int32_t)32,  (int32_t)32,  (int32_t)64);

auto tmp18 = make_vector<uint64_t>( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)64);

auto tmp19 = make_vector<int64_t>( (int32_t)3,  (int32_t)3,  (int32_t)64,  (int32_t)64);

auto tmp20 = make_vector<uint64_t>( (int32_t)3,  (int32_t)3,  (int32_t)64,  (int32_t)64);

auto tmp21 = make_vector<uint64_t>( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)64);

auto tmp22 = make_vector<uint64_t>( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)64);

auto tmp23 = make_vector<int64_t>( (int32_t)3,  (int32_t)3,  (int32_t)64,  (int32_t)64);

auto tmp24 = make_vector<uint64_t>( (int32_t)3,  (int32_t)3,  (int32_t)64,  (int32_t)64);

auto tmp25 = make_vector<uint64_t>( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)64);

auto tmp26 = make_vector<uint64_t>( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)64);

auto tmp27 = make_vector<uint64_t>( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)64);

auto tmp28 = make_vector<int64_t>( (int32_t)3,  (int32_t)3,  (int32_t)64,  (int32_t)64);

auto tmp29 = make_vector<uint64_t>( (int32_t)3,  (int32_t)3,  (int32_t)64,  (int32_t)64);

auto tmp30 = make_vector<uint64_t>( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)64);

auto tmp31 = make_vector<uint64_t>( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)64);

auto tmp32 = make_vector<int64_t>( (int32_t)1,  (int32_t)1,  (int32_t)64,  (int32_t)64);

auto tmp33 = make_vector<uint64_t>( (int32_t)1,  (int32_t)1,  (int32_t)64,  (int32_t)64);

auto tmp34 = make_vector<uint64_t>( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)64);

auto tmp35 = make_vector<uint64_t>( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)64);

auto tmp36 = make_vector<int64_t>( (int32_t)1,  (int32_t)1,  (int32_t)64,  (int32_t)16);

auto tmp37 = make_vector<uint64_t>( (int32_t)1,  (int32_t)1,  (int32_t)64,  (int32_t)16);

auto tmp38 = make_vector<uint64_t>( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)16);

auto tmp39 = make_vector<uint64_t>( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)16);

auto tmp40 = make_vector<int64_t>( (int32_t)1024,  (int32_t)10);

auto tmp41 = make_vector<uint64_t>( (int32_t)1024,  (int32_t)10);

auto tmp42 = make_vector<int64_t>( (int32_t)10);

auto tmp43 = make_vector<uint64_t>( (int32_t)10);

auto tmp44 = make_vector<int32_t>( (int32_t)2);

auto tmp45 = make_vector<uint64_t>( (int32_t)1,  (int32_t)1024);

int64_t i0;

int64_t i1;

int64_t i2;

int64_t i3;

int64_t i4;

int64_t i5;

auto tmp46 = make_vector<uint64_t>( (int32_t)1,  (int32_t)10);

auto tmp47 = make_vector<uint64_t>( (int32_t)1,  (int32_t)10);

auto tmp48 = make_vector<uint64_t>( (int32_t)1);

auto tmp0 = make_vector<uint64_t>( (int32_t)1,  (int32_t)32,  (int32_t)32,  (int32_t)3);
/* Variable to read the clear value corresponding to the input variable tmp0 at (432,1-432,45) */
uint64_t __tmp_in_tmp0;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)1; i0++){
for (uint32_t i1 =  (uint32_t)0; i1 <  (int32_t)32; i1++){
for (uint32_t i2 =  (uint32_t)0; i2 <  (int32_t)32; i2++){
for (uint32_t i3 =  (uint32_t)0; i3 <  (int32_t)3; i3++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp0;
}
tmp0[i0][i1][i2][i3] = (role == CLIENT) ? __tmp_in_tmp0 : 0;
}
}
}
}

auto tmp1 = make_vector<uint64_t>( (int32_t)3,  (int32_t)3,  (int32_t)3,  (int32_t)64);
/* Variable to read the clear value corresponding to the input variable tmp1 at (435,1-435,44) */
uint64_t __tmp_in_tmp1;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)3; i0++){
for (uint32_t i1 =  (uint32_t)0; i1 <  (int32_t)3; i1++){
for (uint32_t i2 =  (uint32_t)0; i2 <  (int32_t)3; i2++){
for (uint32_t i3 =  (uint32_t)0; i3 <  (int32_t)64; i3++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp1;
}
tmp1[i0][i1][i2][i3] = (role == CLIENT) ? __tmp_in_tmp1 : 0;
}
}
}
}

auto tmp2 = make_vector<uint64_t>( (int32_t)3,  (int32_t)3,  (int32_t)64,  (int32_t)64);
/* Variable to read the clear value corresponding to the input variable tmp2 at (438,1-438,45) */
uint64_t __tmp_in_tmp2;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)3; i0++){
for (uint32_t i1 =  (uint32_t)0; i1 <  (int32_t)3; i1++){
for (uint32_t i2 =  (uint32_t)0; i2 <  (int32_t)64; i2++){
for (uint32_t i3 =  (uint32_t)0; i3 <  (int32_t)64; i3++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp2;
}
tmp2[i0][i1][i2][i3] = (role == CLIENT) ? __tmp_in_tmp2 : 0;
}
}
}
}

auto tmp3 = make_vector<uint64_t>( (int32_t)3,  (int32_t)3,  (int32_t)64,  (int32_t)64);
/* Variable to read the clear value corresponding to the input variable tmp3 at (441,1-441,45) */
uint64_t __tmp_in_tmp3;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)3; i0++){
for (uint32_t i1 =  (uint32_t)0; i1 <  (int32_t)3; i1++){
for (uint32_t i2 =  (uint32_t)0; i2 <  (int32_t)64; i2++){
for (uint32_t i3 =  (uint32_t)0; i3 <  (int32_t)64; i3++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp3;
}
tmp3[i0][i1][i2][i3] = (role == CLIENT) ? __tmp_in_tmp3 : 0;
}
}
}
}

auto tmp4 = make_vector<uint64_t>( (int32_t)3,  (int32_t)3,  (int32_t)64,  (int32_t)64);
/* Variable to read the clear value corresponding to the input variable tmp4 at (444,1-444,45) */
uint64_t __tmp_in_tmp4;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)3; i0++){
for (uint32_t i1 =  (uint32_t)0; i1 <  (int32_t)3; i1++){
for (uint32_t i2 =  (uint32_t)0; i2 <  (int32_t)64; i2++){
for (uint32_t i3 =  (uint32_t)0; i3 <  (int32_t)64; i3++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp4;
}
tmp4[i0][i1][i2][i3] = (role == CLIENT) ? __tmp_in_tmp4 : 0;
}
}
}
}

auto tmp5 = make_vector<uint64_t>( (int32_t)3,  (int32_t)3,  (int32_t)64,  (int32_t)64);
/* Variable to read the clear value corresponding to the input variable tmp5 at (447,1-447,45) */
uint64_t __tmp_in_tmp5;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)3; i0++){
for (uint32_t i1 =  (uint32_t)0; i1 <  (int32_t)3; i1++){
for (uint32_t i2 =  (uint32_t)0; i2 <  (int32_t)64; i2++){
for (uint32_t i3 =  (uint32_t)0; i3 <  (int32_t)64; i3++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp5;
}
tmp5[i0][i1][i2][i3] = (role == CLIENT) ? __tmp_in_tmp5 : 0;
}
}
}
}

auto tmp6 = make_vector<uint64_t>( (int32_t)1,  (int32_t)1,  (int32_t)64,  (int32_t)64);
/* Variable to read the clear value corresponding to the input variable tmp6 at (450,1-450,45) */
uint64_t __tmp_in_tmp6;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)1; i0++){
for (uint32_t i1 =  (uint32_t)0; i1 <  (int32_t)1; i1++){
for (uint32_t i2 =  (uint32_t)0; i2 <  (int32_t)64; i2++){
for (uint32_t i3 =  (uint32_t)0; i3 <  (int32_t)64; i3++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp6;
}
tmp6[i0][i1][i2][i3] = (role == CLIENT) ? __tmp_in_tmp6 : 0;
}
}
}
}

auto tmp7 = make_vector<uint64_t>( (int32_t)1,  (int32_t)1,  (int32_t)64,  (int32_t)16);
/* Variable to read the clear value corresponding to the input variable tmp7 at (453,1-453,45) */
uint64_t __tmp_in_tmp7;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)1; i0++){
for (uint32_t i1 =  (uint32_t)0; i1 <  (int32_t)1; i1++){
for (uint32_t i2 =  (uint32_t)0; i2 <  (int32_t)64; i2++){
for (uint32_t i3 =  (uint32_t)0; i3 <  (int32_t)16; i3++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp7;
}
tmp7[i0][i1][i2][i3] = (role == CLIENT) ? __tmp_in_tmp7 : 0;
}
}
}
}

auto tmp8 = make_vector<uint64_t>( (int32_t)1024,  (int32_t)10);
/* Variable to read the clear value corresponding to the input variable tmp8 at (456,1-456,41) */
uint64_t __tmp_in_tmp8;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)1024; i0++){
for (uint32_t i1 =  (uint32_t)0; i1 <  (int32_t)10; i1++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp8;
}
tmp8[i0][i1] = (role == CLIENT) ? __tmp_in_tmp8 : 0;
}
}

auto tmp9 = make_vector<uint64_t>( (int32_t)10);
/* Variable to read the clear value corresponding to the input variable tmp9 at (459,1-459,35) */
uint64_t __tmp_in_tmp9;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)10; i0++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp9;
}
tmp9[i0] = (role == CLIENT) ? __tmp_in_tmp9 : 0;
}

cout<<"Starting 2nd syncronize .. "<<endl;
synchronize(2000000); 
cout<<"Syncronized .. now starting actual execution at "<<getCurrentTime()<<endl;
start_m();
#ifdef PRECOMPUTEAES
auto t1 = high_resolution_clock::now();

#ifndef RUN_SHARECONV_OPTI //If shareconv is off

#ifndef RUN_MSB_OPTI  //If both shareConv and computeMSB are off


#else //If shareConv is off, computeMSB is on


#endif

#else //If share convert opti is ON.

#ifndef RUN_MSB_OPTI //If share conv is on, and msb is off


#else //If share conv is on, and so is msb ; UPDATE : this is when all opti are on

if (partyNum == PARTY_A)
{
	aes_common->PreComputeKeys(262302 + 10, NO_CORES);
	aes_a_1->PreComputeKeys(173581 + 10, NO_CORES);
	aes_b_1->PreComputeKeys(178189 + 10, NO_CORES);
	aes_c_1->PreComputeKeys(173074 + 10, NO_CORES);
	aes_share_conv_bit_shares_p0_p2->PreComputeKeys(706034 + 10, NO_CORES);
	aes_share_conv_shares_mod_odd_p0_p2->PreComputeKeys(129796 + 10, NO_CORES);
	aes_comp_msb_shares_lsb_p0_p2->PreComputeKeys(43263 + 10, NO_CORES);
	aes_comp_msb_shares_bit_vec_p0_p2->PreComputeKeys(43263 + 10, NO_CORES);
	aes_conv_opti_a_1->PreComputeKeys(77151 + 10, NO_CORES);
	aes_conv_opti_b_1->PreComputeKeys(62975 + 10, NO_CORES);
	aes_conv_opti_c_1->PreComputeKeys(86527 + 10, NO_CORES);
}
else if (partyNum == PARTY_B)
{
	aes_common->PreComputeKeys(262302 + 10, NO_CORES);
	aes_a_2->PreComputeKeys(173581 + 10, NO_CORES);
	aes_b_2->PreComputeKeys(178189 + 10, NO_CORES);
	aes_share_conv_bit_shares_p1_p2->PreComputeKeys(706108 + 10, NO_CORES);
	aes_share_conv_shares_mod_odd_p1_p2->PreComputeKeys(129800 + 10, NO_CORES);
	aes_comp_msb_shares_lsb_p1_p2->PreComputeKeys(43268 + 10, NO_CORES);
	aes_comp_msb_shares_bit_vec_p1_p2->PreComputeKeys(43268 + 10, NO_CORES);
	aes_conv_opti_a_2->PreComputeKeys(77151 + 10, NO_CORES);
	aes_conv_opti_b_2->PreComputeKeys(62975 + 10, NO_CORES);
}
else
{
	aes_indep->PreComputeKeys(86532 + 10, NO_CORES);
	aes_a_1->PreComputeKeys(173581 + 10, NO_CORES);
	aes_a_2->PreComputeKeys(173581 + 10, NO_CORES);
	aes_b_1->PreComputeKeys(178189 + 10, NO_CORES);
	aes_b_2->PreComputeKeys(178189 + 10, NO_CORES);
	aes_c_1->PreComputeKeys(173074 + 10, NO_CORES);
	aes_share_conv_bit_shares_p0_p2->PreComputeKeys(706034 + 10, NO_CORES);
	aes_share_conv_bit_shares_p1_p2->PreComputeKeys(706108 + 10, NO_CORES);
	aes_share_conv_shares_mod_odd_p0_p2->PreComputeKeys(129796 + 10, NO_CORES);
	aes_share_conv_shares_mod_odd_p1_p2->PreComputeKeys(129800 + 10, NO_CORES);
	aes_comp_msb_shares_lsb_p0_p2->PreComputeKeys(43263 + 10, NO_CORES);
	aes_comp_msb_shares_lsb_p1_p2->PreComputeKeys(43268 + 10, NO_CORES);
	aes_comp_msb_shares_bit_vec_p0_p2->PreComputeKeys(43263 + 10, NO_CORES);
	aes_comp_msb_shares_bit_vec_p1_p2->PreComputeKeys(43268 + 10, NO_CORES);
	aes_conv_opti_a_1->PreComputeKeys(77151 + 10, NO_CORES);
	aes_conv_opti_a_2->PreComputeKeys(77151 + 10, NO_CORES);
	aes_conv_opti_b_1->PreComputeKeys(62975 + 10, NO_CORES);
	aes_conv_opti_b_2->PreComputeKeys(62975 + 10, NO_CORES);
	aes_conv_opti_c_1->PreComputeKeys(86527 + 10, NO_CORES);
}

#endif

#endif

auto t2 = high_resolution_clock::now();
auto tt = (duration_cast<duration<double>>(t2 - t1)).count();
cout<<"Time for precomputation = "<<tt<<endl;
#endif
CreateTensor4( (int32_t)3,  (int32_t)3,  (int32_t)3,  (int32_t)64,  (int64_t)327, tmp10);
CreateIdentity44( (int32_t)3,  (int32_t)3,  (int32_t)3,  (int32_t)64, tmp1, tmp11);
Conv2DCSFMain( (int32_t)1,  (int32_t)32,  (int32_t)32,  (int32_t)3,  (int32_t)3,  (int32_t)3,  (int32_t)64,  (int32_t)1,  (int32_t)1,  (int32_t)1,  (int32_t)1, tmp0, tmp11, tmp12,  (int64_t)15);
Relu4( (int32_t)1,  (int32_t)32,  (int32_t)32,  (int32_t)64, tmp12, tmp13);
CreateTensor4( (int32_t)3,  (int32_t)3,  (int32_t)64,  (int32_t)64,  (int64_t)327, tmp14);
CreateIdentity44( (int32_t)3,  (int32_t)3,  (int32_t)64,  (int32_t)64, tmp2, tmp15);
Conv2DCSFMain( (int32_t)1,  (int32_t)32,  (int32_t)32,  (int32_t)64,  (int32_t)3,  (int32_t)3,  (int32_t)64,  (int32_t)1,  (int32_t)1,  (int32_t)1,  (int32_t)1, tmp13, tmp15, tmp16,  (int64_t)15);
Relu4( (int32_t)1,  (int32_t)32,  (int32_t)32,  (int32_t)64, tmp16, tmp17);
AvgPool44( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)64,  (int32_t)2,  (int32_t)2,  (int32_t)0,  (int32_t)0,  (int32_t)2,  (int32_t)2,  (int32_t)1,  (int32_t)32,  (int32_t)32,  (int32_t)64, tmp17, tmp18);
CreateTensor4( (int32_t)3,  (int32_t)3,  (int32_t)64,  (int32_t)64,  (int64_t)327, tmp19);
CreateIdentity44( (int32_t)3,  (int32_t)3,  (int32_t)64,  (int32_t)64, tmp3, tmp20);
Conv2DCSFMain( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)64,  (int32_t)3,  (int32_t)3,  (int32_t)64,  (int32_t)1,  (int32_t)1,  (int32_t)1,  (int32_t)1, tmp18, tmp20, tmp21,  (int64_t)15);
Relu4( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)64, tmp21, tmp22);
CreateTensor4( (int32_t)3,  (int32_t)3,  (int32_t)64,  (int32_t)64,  (int64_t)327, tmp23);
CreateIdentity44( (int32_t)3,  (int32_t)3,  (int32_t)64,  (int32_t)64, tmp4, tmp24);
Conv2DCSFMain( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)64,  (int32_t)3,  (int32_t)3,  (int32_t)64,  (int32_t)1,  (int32_t)1,  (int32_t)1,  (int32_t)1, tmp22, tmp24, tmp25,  (int64_t)15);
Relu4( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)64, tmp25, tmp26);
AvgPool44( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)64,  (int32_t)2,  (int32_t)2,  (int32_t)0,  (int32_t)0,  (int32_t)1,  (int32_t)1,  (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)64, tmp26, tmp27);
CreateTensor4( (int32_t)3,  (int32_t)3,  (int32_t)64,  (int32_t)64,  (int64_t)327, tmp28);
CreateIdentity44( (int32_t)3,  (int32_t)3,  (int32_t)64,  (int32_t)64, tmp5, tmp29);
Conv2DCSFMain( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)64,  (int32_t)3,  (int32_t)3,  (int32_t)64,  (int32_t)1,  (int32_t)1,  (int32_t)2,  (int32_t)2, tmp27, tmp29, tmp30,  (int64_t)15);
Relu4( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)64, tmp30, tmp31);
CreateTensor4( (int32_t)1,  (int32_t)1,  (int32_t)64,  (int32_t)64,  (int64_t)327, tmp32);
CreateIdentity44( (int32_t)1,  (int32_t)1,  (int32_t)64,  (int32_t)64, tmp6, tmp33);
Conv2DCSFMain( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)64,  (int32_t)1,  (int32_t)1,  (int32_t)64,  (int32_t)0,  (int32_t)0,  (int32_t)1,  (int32_t)1, tmp31, tmp33, tmp34,  (int64_t)15);
Relu4( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)64, tmp34, tmp35);
CreateTensor4( (int32_t)1,  (int32_t)1,  (int32_t)64,  (int32_t)16,  (int64_t)327, tmp36);
CreateIdentity44( (int32_t)1,  (int32_t)1,  (int32_t)64,  (int32_t)16, tmp7, tmp37);
Conv2DCSFMain( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)64,  (int32_t)1,  (int32_t)1,  (int32_t)16,  (int32_t)0,  (int32_t)0,  (int32_t)1,  (int32_t)1, tmp35, tmp37, tmp38,  (int64_t)15);
Relu4( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)16, tmp38, tmp39);
CreateTensor2( (int32_t)1024,  (int32_t)10,  (int64_t)327, tmp40);
CreateIdentity22( (int32_t)1024,  (int32_t)10, tmp8, tmp41);
CreateTensor1( (int32_t)10,  (int64_t)327, tmp42);
CreateIdentity11( (int32_t)10, tmp9, tmp43);
tmp44[ (int64_t)0] =  (int32_t)-1;
tmp44[ (int64_t)1] =  (int32_t)1024;
i0 =  (int64_t)0;
i1 =  (int64_t)0;
i2 =  (int64_t)0;
i3 =  (int64_t)0;
for (uint32_t i4 =  (int32_t)0; i4 <  (int32_t)1; i4++){
for (uint32_t i5 =  (int32_t)0; i5 <  (int32_t)1024; i5++){
tmp45[i4][i5] = tmp39[i0][i1][i2][i3];
i3 = (i3 +  (int64_t)1);
if ((i3 ==  (int64_t)16)) {
i3 =  (int64_t)0;
i2 = (i2 +  (int64_t)1);
if ((i2 ==  (int64_t)8)) {
i2 =  (int64_t)0;
i1 = (i1 +  (int64_t)1);
if ((i1 ==  (int64_t)8)) {
i1 =  (int64_t)0;
i0 = (i0 +  (int64_t)1);
}
}
}
}
}
MatMulCSF2D( (int32_t)1,  (int32_t)1024,  (int32_t)10, tmp45, tmp41, tmp46,  (int64_t)15);
MatAddBroadCast2( (int32_t)1,  (int32_t)10, tmp46, tmp43, tmp47);
ArgMax1( (int32_t)1,  (int32_t)1,  (int32_t)10, tmp47,  (int32_t)1, tmp48);
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)1; i0++){
add_to_output_queue(out_q, funcReconstruct2PCCons(tmp48[i0], 1), SERVER, cout);
}
end_m(whichNetwork);
 //cout<<"aes_common rCoutner = "<<aes_common->getRCounter()<<endl;
	 //cout<<"aes_indep rCoutner = "<<aes_indep->getRCounter()<<endl;
	 //cout<<"aes_a_1 rCoutner = "<<aes_a_1->getRCounter()<<endl;
	 //cout<<"aes_a_2 rCoutner = "<<aes_a_2->getRCounter()<<endl;
	 //cout<<"aes_b_1 rCoutner = "<<aes_b_1->getRCounter()<<endl;
	 //cout<<"aes_b_2 rCoutner = "<<aes_b_2->getRCounter()<<endl;
	 //cout<<"aes_c_1 rCoutner = "<<aes_c_1->getRCounter()<<endl;
	 //cout<<"aes_share_conv_bit_shares_p0_p2 rCoutner = "<<aes_share_conv_bit_shares_p0_p2->getRCounter()<<endl;
	 //cout<<"aes_share_conv_bit_shares_p1_p2 rCoutner = "<<aes_share_conv_bit_shares_p1_p2->getRCounter()<<endl;
	 //cout<<"aes_share_conv_shares_mod_odd_p0_p2 rCoutner = "<<aes_share_conv_shares_mod_odd_p0_p2->getRCounter()<<endl;
	 //cout<<"aes_share_conv_shares_mod_odd_p1_p2 rCoutner = "<<aes_share_conv_shares_mod_odd_p1_p2->getRCounter()<<endl;
	 //cout<<"aes_comp_msb_shares_lsb_p0_p2 rCoutner = "<<aes_comp_msb_shares_lsb_p0_p2->getRCounter()<<endl;
	 //cout<<"aes_comp_msb_shares_lsb_p1_p2 rCoutner = "<<aes_comp_msb_shares_lsb_p1_p2->getRCounter()<<endl;
	 //cout<<"aes_comp_msb_shares_bit_vec_p0_p2 rCoutner = "<<aes_comp_msb_shares_bit_vec_p0_p2->getRCounter()<<endl;
	 //cout<<"aes_comp_msb_shares_bit_vec_p1_p2 rCoutner = "<<aes_comp_msb_shares_bit_vec_p1_p2->getRCounter()<<endl;
	 //cout<<"aes_conv_opti_a_1 rCoutner = "<<aes_conv_opti_a_1->getRCounter()<<endl;
	 //cout<<"aes_conv_opti_a_2 rCoutner = "<<aes_conv_opti_a_2->getRCounter()<<endl;
	 //cout<<"aes_conv_opti_b_1 rCoutner = "<<aes_conv_opti_b_1->getRCounter()<<endl;
	 //cout<<"aes_conv_opti_b_2 rCoutner = "<<aes_conv_opti_b_2->getRCounter()<<endl;
	 //cout<<"aes_conv_opti_c_1 rCoutner = "<<aes_conv_opti_c_1->getRCounter()<<endl;
cout << "----------------------------------" << endl;
cout << NUM_OF_PARTIES << "PC code, P" << partyNum << endl;
cout << NUM_ITERATIONS << " iterations, " << whichNetwork << ", batch size " << MINI_BATCH_SIZE << endl;
cout << "----------------------------------" << endl << endl;


cout<<"**************RESULTS***************"<<endl;
flush_output_queue(out_q, role);
cout<<"************************************"<<endl;

/****************************** CLEAN-UP ******************************/
delete aes_common;
delete aes_indep;
delete aes_a_1;
delete aes_a_2;
delete aes_b_1;
delete aes_b_2;
delete aes_c_1;
delete aes_share_conv_bit_shares_p0_p2;
delete aes_share_conv_bit_shares_p1_p2;
delete aes_share_conv_shares_mod_odd_p0_p2;
delete aes_share_conv_shares_mod_odd_p1_p2;
delete aes_comp_msb_shares_lsb_p0_p2;
delete aes_comp_msb_shares_lsb_p1_p2;
delete aes_comp_msb_shares_bit_vec_p0_p2;
delete aes_comp_msb_shares_bit_vec_p1_p2;
delete aes_conv_opti_a_1;
delete aes_conv_opti_a_2;
delete aes_conv_opti_b_1;
delete aes_conv_opti_b_2;
delete aes_conv_opti_c_1;
delete aes_parallel;
// delete config;
// delete l0;
// delete l1;
// delete l2;
// delete l3;
// delete network;
if (partyNum != PARTY_S)
deleteObjects();

return 0;

}

#endif
