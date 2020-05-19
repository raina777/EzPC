#include "globals.h"
#ifdef F_SQUEEZENET_CIFAR


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


auto tmp29 = make_vector<int64_t>( (int32_t)3,  (int32_t)3,  (int32_t)3,  (int32_t)64);

auto tmp30 = make_vector<uint64_t>( (int32_t)3,  (int32_t)3,  (int32_t)3,  (int32_t)64);

auto tmp31 = make_vector<int64_t>( (int32_t)64);

auto tmp32 = make_vector<uint64_t>( (int32_t)64);

auto tmp33 = make_vector<uint64_t>( (int32_t)1,  (int32_t)32,  (int32_t)32,  (int32_t)64);

auto tmp34 = make_vector<uint64_t>( (int32_t)1,  (int32_t)32,  (int32_t)32,  (int32_t)64);

auto tmp35 = make_vector<uint64_t>( (int32_t)1,  (int32_t)32,  (int32_t)32,  (int32_t)64);

auto tmp36 = make_vector<uint64_t>( (int32_t)1,  (int32_t)32,  (int32_t)32,  (int32_t)64);

auto tmp37 = make_vector<uint64_t>( (int32_t)1,  (int32_t)32,  (int32_t)32,  (int32_t)64);

auto tmp38 = make_vector<uint64_t>( (int32_t)1,  (int32_t)32,  (int32_t)32,  (int32_t)64);

auto tmp39 = make_vector<uint64_t>( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)64);

auto tmp40 = make_vector<int64_t>( (int32_t)1,  (int32_t)1,  (int32_t)64,  (int32_t)32);

auto tmp41 = make_vector<uint64_t>( (int32_t)1,  (int32_t)1,  (int32_t)64,  (int32_t)32);

auto tmp42 = make_vector<int64_t>( (int32_t)32);

auto tmp43 = make_vector<uint64_t>( (int32_t)32);

auto tmp44 = make_vector<uint64_t>( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)32);

auto tmp45 = make_vector<uint64_t>( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)32);

auto tmp46 = make_vector<uint64_t>( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)32);

auto tmp47 = make_vector<uint64_t>( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)32);

auto tmp48 = make_vector<uint64_t>( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)32);

auto tmp49 = make_vector<uint64_t>( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)32);

auto tmp50 = make_vector<int64_t>( (int32_t)1,  (int32_t)1,  (int32_t)32,  (int32_t)64);

auto tmp51 = make_vector<uint64_t>( (int32_t)1,  (int32_t)1,  (int32_t)32,  (int32_t)64);

auto tmp52 = make_vector<int64_t>( (int32_t)64);

auto tmp53 = make_vector<uint64_t>( (int32_t)64);

auto tmp54 = make_vector<uint64_t>( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)64);

auto tmp55 = make_vector<uint64_t>( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)64);

auto tmp56 = make_vector<uint64_t>( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)64);

auto tmp57 = make_vector<uint64_t>( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)64);

auto tmp58 = make_vector<uint64_t>( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)64);

auto tmp59 = make_vector<uint64_t>( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)64);

auto tmp60 = make_vector<int64_t>( (int32_t)3,  (int32_t)3,  (int32_t)32,  (int32_t)64);

auto tmp61 = make_vector<uint64_t>( (int32_t)3,  (int32_t)3,  (int32_t)32,  (int32_t)64);

auto tmp62 = make_vector<int64_t>( (int32_t)64);

auto tmp63 = make_vector<uint64_t>( (int32_t)64);

auto tmp64 = make_vector<uint64_t>( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)64);

auto tmp65 = make_vector<uint64_t>( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)64);

auto tmp66 = make_vector<uint64_t>( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)64);

auto tmp67 = make_vector<uint64_t>( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)64);

auto tmp68 = make_vector<uint64_t>( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)64);

auto tmp69 = make_vector<uint64_t>( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)64);

auto tmp70 = make_vector<uint64_t>( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)128);

auto tmp71 = make_vector<int64_t>( (int32_t)1,  (int32_t)1,  (int32_t)128,  (int32_t)32);

auto tmp72 = make_vector<uint64_t>( (int32_t)1,  (int32_t)1,  (int32_t)128,  (int32_t)32);

auto tmp73 = make_vector<int64_t>( (int32_t)32);

auto tmp74 = make_vector<uint64_t>( (int32_t)32);

auto tmp75 = make_vector<uint64_t>( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)32);

auto tmp76 = make_vector<uint64_t>( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)32);

auto tmp77 = make_vector<uint64_t>( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)32);

auto tmp78 = make_vector<uint64_t>( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)32);

auto tmp79 = make_vector<uint64_t>( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)32);

auto tmp80 = make_vector<uint64_t>( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)32);

auto tmp81 = make_vector<int64_t>( (int32_t)1,  (int32_t)1,  (int32_t)32,  (int32_t)64);

auto tmp82 = make_vector<uint64_t>( (int32_t)1,  (int32_t)1,  (int32_t)32,  (int32_t)64);

auto tmp83 = make_vector<int64_t>( (int32_t)64);

auto tmp84 = make_vector<uint64_t>( (int32_t)64);

auto tmp85 = make_vector<uint64_t>( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)64);

auto tmp86 = make_vector<uint64_t>( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)64);

auto tmp87 = make_vector<uint64_t>( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)64);

auto tmp88 = make_vector<uint64_t>( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)64);

auto tmp89 = make_vector<uint64_t>( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)64);

auto tmp90 = make_vector<uint64_t>( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)64);

auto tmp91 = make_vector<int64_t>( (int32_t)3,  (int32_t)3,  (int32_t)32,  (int32_t)64);

auto tmp92 = make_vector<uint64_t>( (int32_t)3,  (int32_t)3,  (int32_t)32,  (int32_t)64);

auto tmp93 = make_vector<int64_t>( (int32_t)64);

auto tmp94 = make_vector<uint64_t>( (int32_t)64);

auto tmp95 = make_vector<uint64_t>( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)64);

auto tmp96 = make_vector<uint64_t>( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)64);

auto tmp97 = make_vector<uint64_t>( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)64);

auto tmp98 = make_vector<uint64_t>( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)64);

auto tmp99 = make_vector<uint64_t>( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)64);

auto tmp100 = make_vector<uint64_t>( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)64);

auto tmp101 = make_vector<uint64_t>( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)128);

auto tmp102 = make_vector<uint64_t>( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)128);

auto tmp103 = make_vector<int64_t>( (int32_t)1,  (int32_t)1,  (int32_t)128,  (int32_t)32);

auto tmp104 = make_vector<uint64_t>( (int32_t)1,  (int32_t)1,  (int32_t)128,  (int32_t)32);

auto tmp105 = make_vector<int64_t>( (int32_t)32);

auto tmp106 = make_vector<uint64_t>( (int32_t)32);

auto tmp107 = make_vector<uint64_t>( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)32);

auto tmp108 = make_vector<uint64_t>( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)32);

auto tmp109 = make_vector<uint64_t>( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)32);

auto tmp110 = make_vector<uint64_t>( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)32);

auto tmp111 = make_vector<uint64_t>( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)32);

auto tmp112 = make_vector<uint64_t>( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)32);

auto tmp113 = make_vector<int64_t>( (int32_t)1,  (int32_t)1,  (int32_t)32,  (int32_t)128);

auto tmp114 = make_vector<uint64_t>( (int32_t)1,  (int32_t)1,  (int32_t)32,  (int32_t)128);

auto tmp115 = make_vector<int64_t>( (int32_t)128);

auto tmp116 = make_vector<uint64_t>( (int32_t)128);

auto tmp117 = make_vector<uint64_t>( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)128);

auto tmp118 = make_vector<uint64_t>( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)128);

auto tmp119 = make_vector<uint64_t>( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)128);

auto tmp120 = make_vector<uint64_t>( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)128);

auto tmp121 = make_vector<uint64_t>( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)128);

auto tmp122 = make_vector<uint64_t>( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)128);

auto tmp123 = make_vector<int64_t>( (int32_t)3,  (int32_t)3,  (int32_t)32,  (int32_t)128);

auto tmp124 = make_vector<uint64_t>( (int32_t)3,  (int32_t)3,  (int32_t)32,  (int32_t)128);

auto tmp125 = make_vector<int64_t>( (int32_t)128);

auto tmp126 = make_vector<uint64_t>( (int32_t)128);

auto tmp127 = make_vector<uint64_t>( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)128);

auto tmp128 = make_vector<uint64_t>( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)128);

auto tmp129 = make_vector<uint64_t>( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)128);

auto tmp130 = make_vector<uint64_t>( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)128);

auto tmp131 = make_vector<uint64_t>( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)128);

auto tmp132 = make_vector<uint64_t>( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)128);

auto tmp133 = make_vector<uint64_t>( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)256);

auto tmp134 = make_vector<int64_t>( (int32_t)1,  (int32_t)1,  (int32_t)256,  (int32_t)32);

auto tmp135 = make_vector<uint64_t>( (int32_t)1,  (int32_t)1,  (int32_t)256,  (int32_t)32);

auto tmp136 = make_vector<int64_t>( (int32_t)32);

auto tmp137 = make_vector<uint64_t>( (int32_t)32);

auto tmp138 = make_vector<uint64_t>( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)32);

auto tmp139 = make_vector<uint64_t>( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)32);

auto tmp140 = make_vector<uint64_t>( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)32);

auto tmp141 = make_vector<uint64_t>( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)32);

auto tmp142 = make_vector<uint64_t>( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)32);

auto tmp143 = make_vector<uint64_t>( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)32);

auto tmp144 = make_vector<int64_t>( (int32_t)1,  (int32_t)1,  (int32_t)32,  (int32_t)128);

auto tmp145 = make_vector<uint64_t>( (int32_t)1,  (int32_t)1,  (int32_t)32,  (int32_t)128);

auto tmp146 = make_vector<int64_t>( (int32_t)128);

auto tmp147 = make_vector<uint64_t>( (int32_t)128);

auto tmp148 = make_vector<uint64_t>( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)128);

auto tmp149 = make_vector<uint64_t>( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)128);

auto tmp150 = make_vector<uint64_t>( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)128);

auto tmp151 = make_vector<uint64_t>( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)128);

auto tmp152 = make_vector<uint64_t>( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)128);

auto tmp153 = make_vector<uint64_t>( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)128);

auto tmp154 = make_vector<int64_t>( (int32_t)3,  (int32_t)3,  (int32_t)32,  (int32_t)128);

auto tmp155 = make_vector<uint64_t>( (int32_t)3,  (int32_t)3,  (int32_t)32,  (int32_t)128);

auto tmp156 = make_vector<int64_t>( (int32_t)128);

auto tmp157 = make_vector<uint64_t>( (int32_t)128);

auto tmp158 = make_vector<uint64_t>( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)128);

auto tmp159 = make_vector<uint64_t>( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)128);

auto tmp160 = make_vector<uint64_t>( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)128);

auto tmp161 = make_vector<uint64_t>( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)128);

auto tmp162 = make_vector<uint64_t>( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)128);

auto tmp163 = make_vector<uint64_t>( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)128);

auto tmp164 = make_vector<uint64_t>( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)256);

auto tmp165 = make_vector<int64_t>( (int32_t)1,  (int32_t)1,  (int32_t)256,  (int32_t)10);

auto tmp166 = make_vector<uint64_t>( (int32_t)1,  (int32_t)1,  (int32_t)256,  (int32_t)10);

auto tmp167 = make_vector<int64_t>( (int32_t)10);

auto tmp168 = make_vector<uint64_t>( (int32_t)10);

auto tmp169 = make_vector<uint64_t>( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)10);

auto tmp170 = make_vector<uint64_t>( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)10);

auto tmp171 = make_vector<uint64_t>( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)10);

auto tmp172 = make_vector<uint64_t>( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)10);

auto tmp173 = make_vector<uint64_t>( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)10);

auto tmp174 = make_vector<uint64_t>( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)10);

auto tmp175 = make_vector<uint64_t>( (int32_t)1,  (int32_t)1,  (int32_t)1,  (int32_t)10);

auto tmp176 = make_vector<int32_t>( (int32_t)2);

auto tmp177 = make_vector<uint64_t>( (int32_t)1,  (int32_t)10);

int64_t i0;

int64_t i1;

int64_t i2;

int64_t i3;

int64_t i4;

int64_t i5;

auto tmp178 = make_vector<uint64_t>( (int32_t)1);

auto tmp0 = make_vector<uint64_t>( (int32_t)1,  (int32_t)32,  (int32_t)32,  (int32_t)3);
/* Variable to read the clear value corresponding to the input variable tmp0 at (543,1-543,45) */
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
/* Variable to read the clear value corresponding to the input variable tmp1 at (546,1-546,44) */
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

auto tmp2 = make_vector<uint64_t>( (int32_t)64);
/* Variable to read the clear value corresponding to the input variable tmp2 at (549,1-549,35) */
uint64_t __tmp_in_tmp2;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)64; i0++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp2;
}
tmp2[i0] = (role == CLIENT) ? __tmp_in_tmp2 : 0;
}

auto tmp3 = make_vector<uint64_t>( (int32_t)1,  (int32_t)1,  (int32_t)64,  (int32_t)32);
/* Variable to read the clear value corresponding to the input variable tmp3 at (552,1-552,45) */
uint64_t __tmp_in_tmp3;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)1; i0++){
for (uint32_t i1 =  (uint32_t)0; i1 <  (int32_t)1; i1++){
for (uint32_t i2 =  (uint32_t)0; i2 <  (int32_t)64; i2++){
for (uint32_t i3 =  (uint32_t)0; i3 <  (int32_t)32; i3++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp3;
}
tmp3[i0][i1][i2][i3] = (role == CLIENT) ? __tmp_in_tmp3 : 0;
}
}
}
}

auto tmp4 = make_vector<uint64_t>( (int32_t)32);
/* Variable to read the clear value corresponding to the input variable tmp4 at (555,1-555,35) */
uint64_t __tmp_in_tmp4;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)32; i0++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp4;
}
tmp4[i0] = (role == CLIENT) ? __tmp_in_tmp4 : 0;
}

auto tmp5 = make_vector<uint64_t>( (int32_t)1,  (int32_t)1,  (int32_t)32,  (int32_t)64);
/* Variable to read the clear value corresponding to the input variable tmp5 at (558,1-558,45) */
uint64_t __tmp_in_tmp5;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)1; i0++){
for (uint32_t i1 =  (uint32_t)0; i1 <  (int32_t)1; i1++){
for (uint32_t i2 =  (uint32_t)0; i2 <  (int32_t)32; i2++){
for (uint32_t i3 =  (uint32_t)0; i3 <  (int32_t)64; i3++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp5;
}
tmp5[i0][i1][i2][i3] = (role == CLIENT) ? __tmp_in_tmp5 : 0;
}
}
}
}

auto tmp6 = make_vector<uint64_t>( (int32_t)64);
/* Variable to read the clear value corresponding to the input variable tmp6 at (561,1-561,35) */
uint64_t __tmp_in_tmp6;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)64; i0++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp6;
}
tmp6[i0] = (role == CLIENT) ? __tmp_in_tmp6 : 0;
}

auto tmp7 = make_vector<uint64_t>( (int32_t)3,  (int32_t)3,  (int32_t)32,  (int32_t)64);
/* Variable to read the clear value corresponding to the input variable tmp7 at (564,1-564,45) */
uint64_t __tmp_in_tmp7;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)3; i0++){
for (uint32_t i1 =  (uint32_t)0; i1 <  (int32_t)3; i1++){
for (uint32_t i2 =  (uint32_t)0; i2 <  (int32_t)32; i2++){
for (uint32_t i3 =  (uint32_t)0; i3 <  (int32_t)64; i3++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp7;
}
tmp7[i0][i1][i2][i3] = (role == CLIENT) ? __tmp_in_tmp7 : 0;
}
}
}
}

auto tmp8 = make_vector<uint64_t>( (int32_t)64);
/* Variable to read the clear value corresponding to the input variable tmp8 at (567,1-567,35) */
uint64_t __tmp_in_tmp8;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)64; i0++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp8;
}
tmp8[i0] = (role == CLIENT) ? __tmp_in_tmp8 : 0;
}

auto tmp9 = make_vector<uint64_t>( (int32_t)1,  (int32_t)1,  (int32_t)128,  (int32_t)32);
/* Variable to read the clear value corresponding to the input variable tmp9 at (570,1-570,46) */
uint64_t __tmp_in_tmp9;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)1; i0++){
for (uint32_t i1 =  (uint32_t)0; i1 <  (int32_t)1; i1++){
for (uint32_t i2 =  (uint32_t)0; i2 <  (int32_t)128; i2++){
for (uint32_t i3 =  (uint32_t)0; i3 <  (int32_t)32; i3++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp9;
}
tmp9[i0][i1][i2][i3] = (role == CLIENT) ? __tmp_in_tmp9 : 0;
}
}
}
}

auto tmp10 = make_vector<uint64_t>( (int32_t)32);
/* Variable to read the clear value corresponding to the input variable tmp10 at (573,1-573,36) */
uint64_t __tmp_in_tmp10;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)32; i0++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp10;
}
tmp10[i0] = (role == CLIENT) ? __tmp_in_tmp10 : 0;
}

auto tmp11 = make_vector<uint64_t>( (int32_t)1,  (int32_t)1,  (int32_t)32,  (int32_t)64);
/* Variable to read the clear value corresponding to the input variable tmp11 at (576,1-576,46) */
uint64_t __tmp_in_tmp11;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)1; i0++){
for (uint32_t i1 =  (uint32_t)0; i1 <  (int32_t)1; i1++){
for (uint32_t i2 =  (uint32_t)0; i2 <  (int32_t)32; i2++){
for (uint32_t i3 =  (uint32_t)0; i3 <  (int32_t)64; i3++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp11;
}
tmp11[i0][i1][i2][i3] = (role == CLIENT) ? __tmp_in_tmp11 : 0;
}
}
}
}

auto tmp12 = make_vector<uint64_t>( (int32_t)64);
/* Variable to read the clear value corresponding to the input variable tmp12 at (579,1-579,36) */
uint64_t __tmp_in_tmp12;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)64; i0++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp12;
}
tmp12[i0] = (role == CLIENT) ? __tmp_in_tmp12 : 0;
}

auto tmp13 = make_vector<uint64_t>( (int32_t)3,  (int32_t)3,  (int32_t)32,  (int32_t)64);
/* Variable to read the clear value corresponding to the input variable tmp13 at (582,1-582,46) */
uint64_t __tmp_in_tmp13;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)3; i0++){
for (uint32_t i1 =  (uint32_t)0; i1 <  (int32_t)3; i1++){
for (uint32_t i2 =  (uint32_t)0; i2 <  (int32_t)32; i2++){
for (uint32_t i3 =  (uint32_t)0; i3 <  (int32_t)64; i3++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp13;
}
tmp13[i0][i1][i2][i3] = (role == CLIENT) ? __tmp_in_tmp13 : 0;
}
}
}
}

auto tmp14 = make_vector<uint64_t>( (int32_t)64);
/* Variable to read the clear value corresponding to the input variable tmp14 at (585,1-585,36) */
uint64_t __tmp_in_tmp14;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)64; i0++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp14;
}
tmp14[i0] = (role == CLIENT) ? __tmp_in_tmp14 : 0;
}

auto tmp15 = make_vector<uint64_t>( (int32_t)1,  (int32_t)1,  (int32_t)128,  (int32_t)32);
/* Variable to read the clear value corresponding to the input variable tmp15 at (588,1-588,47) */
uint64_t __tmp_in_tmp15;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)1; i0++){
for (uint32_t i1 =  (uint32_t)0; i1 <  (int32_t)1; i1++){
for (uint32_t i2 =  (uint32_t)0; i2 <  (int32_t)128; i2++){
for (uint32_t i3 =  (uint32_t)0; i3 <  (int32_t)32; i3++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp15;
}
tmp15[i0][i1][i2][i3] = (role == CLIENT) ? __tmp_in_tmp15 : 0;
}
}
}
}

auto tmp16 = make_vector<uint64_t>( (int32_t)32);
/* Variable to read the clear value corresponding to the input variable tmp16 at (591,1-591,36) */
uint64_t __tmp_in_tmp16;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)32; i0++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp16;
}
tmp16[i0] = (role == CLIENT) ? __tmp_in_tmp16 : 0;
}

auto tmp17 = make_vector<uint64_t>( (int32_t)1,  (int32_t)1,  (int32_t)32,  (int32_t)128);
/* Variable to read the clear value corresponding to the input variable tmp17 at (594,1-594,47) */
uint64_t __tmp_in_tmp17;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)1; i0++){
for (uint32_t i1 =  (uint32_t)0; i1 <  (int32_t)1; i1++){
for (uint32_t i2 =  (uint32_t)0; i2 <  (int32_t)32; i2++){
for (uint32_t i3 =  (uint32_t)0; i3 <  (int32_t)128; i3++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp17;
}
tmp17[i0][i1][i2][i3] = (role == CLIENT) ? __tmp_in_tmp17 : 0;
}
}
}
}

auto tmp18 = make_vector<uint64_t>( (int32_t)128);
/* Variable to read the clear value corresponding to the input variable tmp18 at (597,1-597,37) */
uint64_t __tmp_in_tmp18;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)128; i0++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp18;
}
tmp18[i0] = (role == CLIENT) ? __tmp_in_tmp18 : 0;
}

auto tmp19 = make_vector<uint64_t>( (int32_t)3,  (int32_t)3,  (int32_t)32,  (int32_t)128);
/* Variable to read the clear value corresponding to the input variable tmp19 at (600,1-600,47) */
uint64_t __tmp_in_tmp19;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)3; i0++){
for (uint32_t i1 =  (uint32_t)0; i1 <  (int32_t)3; i1++){
for (uint32_t i2 =  (uint32_t)0; i2 <  (int32_t)32; i2++){
for (uint32_t i3 =  (uint32_t)0; i3 <  (int32_t)128; i3++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp19;
}
tmp19[i0][i1][i2][i3] = (role == CLIENT) ? __tmp_in_tmp19 : 0;
}
}
}
}

auto tmp20 = make_vector<uint64_t>( (int32_t)128);
/* Variable to read the clear value corresponding to the input variable tmp20 at (603,1-603,37) */
uint64_t __tmp_in_tmp20;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)128; i0++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp20;
}
tmp20[i0] = (role == CLIENT) ? __tmp_in_tmp20 : 0;
}

auto tmp21 = make_vector<uint64_t>( (int32_t)1,  (int32_t)1,  (int32_t)256,  (int32_t)32);
/* Variable to read the clear value corresponding to the input variable tmp21 at (606,1-606,47) */
uint64_t __tmp_in_tmp21;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)1; i0++){
for (uint32_t i1 =  (uint32_t)0; i1 <  (int32_t)1; i1++){
for (uint32_t i2 =  (uint32_t)0; i2 <  (int32_t)256; i2++){
for (uint32_t i3 =  (uint32_t)0; i3 <  (int32_t)32; i3++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp21;
}
tmp21[i0][i1][i2][i3] = (role == CLIENT) ? __tmp_in_tmp21 : 0;
}
}
}
}

auto tmp22 = make_vector<uint64_t>( (int32_t)32);
/* Variable to read the clear value corresponding to the input variable tmp22 at (609,1-609,36) */
uint64_t __tmp_in_tmp22;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)32; i0++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp22;
}
tmp22[i0] = (role == CLIENT) ? __tmp_in_tmp22 : 0;
}

auto tmp23 = make_vector<uint64_t>( (int32_t)1,  (int32_t)1,  (int32_t)32,  (int32_t)128);
/* Variable to read the clear value corresponding to the input variable tmp23 at (612,1-612,47) */
uint64_t __tmp_in_tmp23;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)1; i0++){
for (uint32_t i1 =  (uint32_t)0; i1 <  (int32_t)1; i1++){
for (uint32_t i2 =  (uint32_t)0; i2 <  (int32_t)32; i2++){
for (uint32_t i3 =  (uint32_t)0; i3 <  (int32_t)128; i3++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp23;
}
tmp23[i0][i1][i2][i3] = (role == CLIENT) ? __tmp_in_tmp23 : 0;
}
}
}
}

auto tmp24 = make_vector<uint64_t>( (int32_t)128);
/* Variable to read the clear value corresponding to the input variable tmp24 at (615,1-615,37) */
uint64_t __tmp_in_tmp24;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)128; i0++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp24;
}
tmp24[i0] = (role == CLIENT) ? __tmp_in_tmp24 : 0;
}

auto tmp25 = make_vector<uint64_t>( (int32_t)3,  (int32_t)3,  (int32_t)32,  (int32_t)128);
/* Variable to read the clear value corresponding to the input variable tmp25 at (618,1-618,47) */
uint64_t __tmp_in_tmp25;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)3; i0++){
for (uint32_t i1 =  (uint32_t)0; i1 <  (int32_t)3; i1++){
for (uint32_t i2 =  (uint32_t)0; i2 <  (int32_t)32; i2++){
for (uint32_t i3 =  (uint32_t)0; i3 <  (int32_t)128; i3++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp25;
}
tmp25[i0][i1][i2][i3] = (role == CLIENT) ? __tmp_in_tmp25 : 0;
}
}
}
}

auto tmp26 = make_vector<uint64_t>( (int32_t)128);
/* Variable to read the clear value corresponding to the input variable tmp26 at (621,1-621,37) */
uint64_t __tmp_in_tmp26;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)128; i0++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp26;
}
tmp26[i0] = (role == CLIENT) ? __tmp_in_tmp26 : 0;
}

auto tmp27 = make_vector<uint64_t>( (int32_t)1,  (int32_t)1,  (int32_t)256,  (int32_t)10);
/* Variable to read the clear value corresponding to the input variable tmp27 at (624,1-624,47) */
uint64_t __tmp_in_tmp27;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)1; i0++){
for (uint32_t i1 =  (uint32_t)0; i1 <  (int32_t)1; i1++){
for (uint32_t i2 =  (uint32_t)0; i2 <  (int32_t)256; i2++){
for (uint32_t i3 =  (uint32_t)0; i3 <  (int32_t)10; i3++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp27;
}
tmp27[i0][i1][i2][i3] = (role == CLIENT) ? __tmp_in_tmp27 : 0;
}
}
}
}

auto tmp28 = make_vector<uint64_t>( (int32_t)10);
/* Variable to read the clear value corresponding to the input variable tmp28 at (627,1-627,36) */
uint64_t __tmp_in_tmp28;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)10; i0++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp28;
}
tmp28[i0] = (role == CLIENT) ? __tmp_in_tmp28 : 0;
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
	aes_common->PreComputeKeys(14 + 10, NO_CORES);
	aes_a_1->PreComputeKeys(92493 + 10, NO_CORES);
	aes_b_1->PreComputeKeys(92493 + 10, NO_CORES);
	aes_c_1->PreComputeKeys(92493 + 10, NO_CORES);
	aes_share_conv_bit_shares_p0_p2->PreComputeKeys(4 + 10, NO_CORES);
	aes_share_conv_shares_mod_odd_p0_p2->PreComputeKeys(0 + 10, NO_CORES);
	aes_comp_msb_shares_lsb_p0_p2->PreComputeKeys(0 + 10, NO_CORES);
	aes_comp_msb_shares_bit_vec_p0_p2->PreComputeKeys(0 + 10, NO_CORES);
	aes_conv_opti_a_1->PreComputeKeys(72799 + 10, NO_CORES);
	aes_conv_opti_b_1->PreComputeKeys(67071 + 10, NO_CORES);
	aes_conv_opti_c_1->PreComputeKeys(92479 + 10, NO_CORES);
}
else if (partyNum == PARTY_B)
{
	aes_common->PreComputeKeys(14 + 10, NO_CORES);
	aes_a_2->PreComputeKeys(92493 + 10, NO_CORES);
	aes_b_2->PreComputeKeys(92493 + 10, NO_CORES);
	aes_share_conv_bit_shares_p1_p2->PreComputeKeys(73 + 10, NO_CORES);
	aes_share_conv_shares_mod_odd_p1_p2->PreComputeKeys(8 + 10, NO_CORES);
	aes_comp_msb_shares_lsb_p1_p2->PreComputeKeys(4 + 10, NO_CORES);
	aes_comp_msb_shares_bit_vec_p1_p2->PreComputeKeys(4 + 10, NO_CORES);
	aes_conv_opti_a_2->PreComputeKeys(72799 + 10, NO_CORES);
	aes_conv_opti_b_2->PreComputeKeys(67071 + 10, NO_CORES);
}
else
{
	aes_indep->PreComputeKeys(4 + 10, NO_CORES);
	aes_a_1->PreComputeKeys(92493 + 10, NO_CORES);
	aes_a_2->PreComputeKeys(92493 + 10, NO_CORES);
	aes_b_1->PreComputeKeys(92493 + 10, NO_CORES);
	aes_b_2->PreComputeKeys(92493 + 10, NO_CORES);
	aes_c_1->PreComputeKeys(92493 + 10, NO_CORES);
	aes_share_conv_bit_shares_p0_p2->PreComputeKeys(0 + 10, NO_CORES);
	aes_share_conv_bit_shares_p1_p2->PreComputeKeys(73 + 10, NO_CORES);
	aes_share_conv_shares_mod_odd_p0_p2->PreComputeKeys(4 + 10, NO_CORES);
	aes_share_conv_shares_mod_odd_p1_p2->PreComputeKeys(8 + 10, NO_CORES);
	aes_comp_msb_shares_lsb_p0_p2->PreComputeKeys(0 + 10, NO_CORES);
	aes_comp_msb_shares_lsb_p1_p2->PreComputeKeys(4 + 10, NO_CORES);
	aes_comp_msb_shares_bit_vec_p0_p2->PreComputeKeys(0 + 10, NO_CORES);
	aes_comp_msb_shares_bit_vec_p1_p2->PreComputeKeys(4 + 10, NO_CORES);
	aes_conv_opti_a_1->PreComputeKeys(72799 + 10, NO_CORES);
	aes_conv_opti_a_2->PreComputeKeys(72799 + 10, NO_CORES);
	aes_conv_opti_b_1->PreComputeKeys(67071 + 10, NO_CORES);
	aes_conv_opti_b_2->PreComputeKeys(67071 + 10, NO_CORES);
	aes_conv_opti_c_1->PreComputeKeys(92479 + 10, NO_CORES);
}

#endif

#endif

auto t2 = high_resolution_clock::now();
auto tt = (duration_cast<duration<double>>(t2 - t1)).count();
cout<<"Time for precomputation = "<<tt<<endl;
#endif
CreateTensor4( (int32_t)3,  (int32_t)3,  (int32_t)3,  (int32_t)64,  (int64_t)327, tmp29);
CreateIdentity44( (int32_t)3,  (int32_t)3,  (int32_t)3,  (int32_t)64, tmp1, tmp30);
CreateTensor1( (int32_t)64,  (int64_t)0, tmp31);
CreateIdentity11( (int32_t)64, tmp2, tmp32);
Conv2DCSFMain( (int32_t)1,  (int32_t)32,  (int32_t)32,  (int32_t)3,  (int32_t)3,  (int32_t)3,  (int32_t)64,  (int32_t)1,  (int32_t)1,  (int32_t)1,  (int32_t)1, tmp0, tmp30, tmp33,  (int64_t)15);
MatAddBroadCast4( (int32_t)1,  (int32_t)32,  (int32_t)32,  (int32_t)64, tmp33, tmp32, tmp34);
ElemWiseMul4( (int32_t)1,  (int32_t)32,  (int32_t)32,  (int32_t)64, tmp34, tmp34, tmp35,  (int64_t)15);
ScalarMul4( (int32_t)1,  (int32_t)32,  (int32_t)32,  (int32_t)64,  (int64_t)327, tmp35, tmp36,  (int64_t)15);
ScalarMul4( (int32_t)1,  (int32_t)32,  (int32_t)32,  (int32_t)64,  (int64_t)3276, tmp34, tmp37,  (int64_t)15);
MatAdd4( (int32_t)1,  (int32_t)32,  (int32_t)32,  (int32_t)64, tmp36, tmp37, tmp38);
AvgPool44( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)64,  (int32_t)3,  (int32_t)3,  (int32_t)1,  (int32_t)1,  (int32_t)2,  (int32_t)2,  (int32_t)1,  (int32_t)32,  (int32_t)32,  (int32_t)64, tmp38, tmp39);
CreateTensor4( (int32_t)1,  (int32_t)1,  (int32_t)64,  (int32_t)32,  (int64_t)327, tmp40);
CreateIdentity44( (int32_t)1,  (int32_t)1,  (int32_t)64,  (int32_t)32, tmp3, tmp41);
CreateTensor1( (int32_t)32,  (int64_t)0, tmp42);
CreateIdentity11( (int32_t)32, tmp4, tmp43);
Conv2DCSFMain( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)64,  (int32_t)1,  (int32_t)1,  (int32_t)32,  (int32_t)0,  (int32_t)0,  (int32_t)1,  (int32_t)1, tmp39, tmp41, tmp44,  (int64_t)15);
MatAddBroadCast4( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)32, tmp44, tmp43, tmp45);
ElemWiseMul4( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)32, tmp45, tmp45, tmp46,  (int64_t)15);
ScalarMul4( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)32,  (int64_t)327, tmp46, tmp47,  (int64_t)15);
ScalarMul4( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)32,  (int64_t)3276, tmp45, tmp48,  (int64_t)15);
MatAdd4( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)32, tmp47, tmp48, tmp49);
CreateTensor4( (int32_t)1,  (int32_t)1,  (int32_t)32,  (int32_t)64,  (int64_t)327, tmp50);
CreateIdentity44( (int32_t)1,  (int32_t)1,  (int32_t)32,  (int32_t)64, tmp5, tmp51);
CreateTensor1( (int32_t)64,  (int64_t)0, tmp52);
CreateIdentity11( (int32_t)64, tmp6, tmp53);
Conv2DCSFMain( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)32,  (int32_t)1,  (int32_t)1,  (int32_t)64,  (int32_t)0,  (int32_t)0,  (int32_t)1,  (int32_t)1, tmp49, tmp51, tmp54,  (int64_t)15);
MatAddBroadCast4( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)64, tmp54, tmp53, tmp55);
ElemWiseMul4( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)64, tmp55, tmp55, tmp56,  (int64_t)15);
ScalarMul4( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)64,  (int64_t)327, tmp56, tmp57,  (int64_t)15);
ScalarMul4( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)64,  (int64_t)3276, tmp55, tmp58,  (int64_t)15);
MatAdd4( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)64, tmp57, tmp58, tmp59);
CreateTensor4( (int32_t)3,  (int32_t)3,  (int32_t)32,  (int32_t)64,  (int64_t)327, tmp60);
CreateIdentity44( (int32_t)3,  (int32_t)3,  (int32_t)32,  (int32_t)64, tmp7, tmp61);
CreateTensor1( (int32_t)64,  (int64_t)0, tmp62);
CreateIdentity11( (int32_t)64, tmp8, tmp63);
Conv2DCSFMain( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)32,  (int32_t)3,  (int32_t)3,  (int32_t)64,  (int32_t)1,  (int32_t)1,  (int32_t)1,  (int32_t)1, tmp49, tmp61, tmp64,  (int64_t)15);
MatAddBroadCast4( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)64, tmp64, tmp63, tmp65);
ElemWiseMul4( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)64, tmp65, tmp65, tmp66,  (int64_t)15);
ScalarMul4( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)64,  (int64_t)327, tmp66, tmp67,  (int64_t)15);
ScalarMul4( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)64,  (int64_t)3276, tmp65, tmp68,  (int64_t)15);
MatAdd4( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)64, tmp67, tmp68, tmp69);
Concat2T444( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)128,  (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)64, tmp59,  (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)64, tmp69,  (int32_t)3, tmp70);
CreateTensor4( (int32_t)1,  (int32_t)1,  (int32_t)128,  (int32_t)32,  (int64_t)327, tmp71);
CreateIdentity44( (int32_t)1,  (int32_t)1,  (int32_t)128,  (int32_t)32, tmp9, tmp72);
CreateTensor1( (int32_t)32,  (int64_t)0, tmp73);
CreateIdentity11( (int32_t)32, tmp10, tmp74);
Conv2DCSFMain( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)128,  (int32_t)1,  (int32_t)1,  (int32_t)32,  (int32_t)0,  (int32_t)0,  (int32_t)1,  (int32_t)1, tmp70, tmp72, tmp75,  (int64_t)15);
MatAddBroadCast4( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)32, tmp75, tmp74, tmp76);
ElemWiseMul4( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)32, tmp76, tmp76, tmp77,  (int64_t)15);
ScalarMul4( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)32,  (int64_t)327, tmp77, tmp78,  (int64_t)15);
ScalarMul4( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)32,  (int64_t)3276, tmp76, tmp79,  (int64_t)15);
MatAdd4( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)32, tmp78, tmp79, tmp80);
CreateTensor4( (int32_t)1,  (int32_t)1,  (int32_t)32,  (int32_t)64,  (int64_t)327, tmp81);
CreateIdentity44( (int32_t)1,  (int32_t)1,  (int32_t)32,  (int32_t)64, tmp11, tmp82);
CreateTensor1( (int32_t)64,  (int64_t)0, tmp83);
CreateIdentity11( (int32_t)64, tmp12, tmp84);
Conv2DCSFMain( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)32,  (int32_t)1,  (int32_t)1,  (int32_t)64,  (int32_t)0,  (int32_t)0,  (int32_t)1,  (int32_t)1, tmp80, tmp82, tmp85,  (int64_t)15);
MatAddBroadCast4( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)64, tmp85, tmp84, tmp86);
ElemWiseMul4( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)64, tmp86, tmp86, tmp87,  (int64_t)15);
ScalarMul4( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)64,  (int64_t)327, tmp87, tmp88,  (int64_t)15);
ScalarMul4( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)64,  (int64_t)3276, tmp86, tmp89,  (int64_t)15);
MatAdd4( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)64, tmp88, tmp89, tmp90);
CreateTensor4( (int32_t)3,  (int32_t)3,  (int32_t)32,  (int32_t)64,  (int64_t)327, tmp91);
CreateIdentity44( (int32_t)3,  (int32_t)3,  (int32_t)32,  (int32_t)64, tmp13, tmp92);
CreateTensor1( (int32_t)64,  (int64_t)0, tmp93);
CreateIdentity11( (int32_t)64, tmp14, tmp94);
Conv2DCSFMain( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)32,  (int32_t)3,  (int32_t)3,  (int32_t)64,  (int32_t)1,  (int32_t)1,  (int32_t)1,  (int32_t)1, tmp80, tmp92, tmp95,  (int64_t)15);
MatAddBroadCast4( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)64, tmp95, tmp94, tmp96);
ElemWiseMul4( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)64, tmp96, tmp96, tmp97,  (int64_t)15);
ScalarMul4( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)64,  (int64_t)327, tmp97, tmp98,  (int64_t)15);
ScalarMul4( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)64,  (int64_t)3276, tmp96, tmp99,  (int64_t)15);
MatAdd4( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)64, tmp98, tmp99, tmp100);
Concat2T444( (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)128,  (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)64, tmp90,  (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)64, tmp100,  (int32_t)3, tmp101);
AvgPool44( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)128,  (int32_t)3,  (int32_t)3,  (int32_t)1,  (int32_t)1,  (int32_t)2,  (int32_t)2,  (int32_t)1,  (int32_t)16,  (int32_t)16,  (int32_t)128, tmp101, tmp102);
CreateTensor4( (int32_t)1,  (int32_t)1,  (int32_t)128,  (int32_t)32,  (int64_t)327, tmp103);
CreateIdentity44( (int32_t)1,  (int32_t)1,  (int32_t)128,  (int32_t)32, tmp15, tmp104);
CreateTensor1( (int32_t)32,  (int64_t)0, tmp105);
CreateIdentity11( (int32_t)32, tmp16, tmp106);
Conv2DCSFMain( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)128,  (int32_t)1,  (int32_t)1,  (int32_t)32,  (int32_t)0,  (int32_t)0,  (int32_t)1,  (int32_t)1, tmp102, tmp104, tmp107,  (int64_t)15);
MatAddBroadCast4( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)32, tmp107, tmp106, tmp108);
ElemWiseMul4( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)32, tmp108, tmp108, tmp109,  (int64_t)15);
ScalarMul4( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)32,  (int64_t)327, tmp109, tmp110,  (int64_t)15);
ScalarMul4( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)32,  (int64_t)3276, tmp108, tmp111,  (int64_t)15);
MatAdd4( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)32, tmp110, tmp111, tmp112);
CreateTensor4( (int32_t)1,  (int32_t)1,  (int32_t)32,  (int32_t)128,  (int64_t)327, tmp113);
CreateIdentity44( (int32_t)1,  (int32_t)1,  (int32_t)32,  (int32_t)128, tmp17, tmp114);
CreateTensor1( (int32_t)128,  (int64_t)0, tmp115);
CreateIdentity11( (int32_t)128, tmp18, tmp116);
Conv2DCSFMain( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)32,  (int32_t)1,  (int32_t)1,  (int32_t)128,  (int32_t)0,  (int32_t)0,  (int32_t)1,  (int32_t)1, tmp112, tmp114, tmp117,  (int64_t)15);
MatAddBroadCast4( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)128, tmp117, tmp116, tmp118);
ElemWiseMul4( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)128, tmp118, tmp118, tmp119,  (int64_t)15);
ScalarMul4( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)128,  (int64_t)327, tmp119, tmp120,  (int64_t)15);
ScalarMul4( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)128,  (int64_t)3276, tmp118, tmp121,  (int64_t)15);
MatAdd4( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)128, tmp120, tmp121, tmp122);
CreateTensor4( (int32_t)3,  (int32_t)3,  (int32_t)32,  (int32_t)128,  (int64_t)327, tmp123);
CreateIdentity44( (int32_t)3,  (int32_t)3,  (int32_t)32,  (int32_t)128, tmp19, tmp124);
CreateTensor1( (int32_t)128,  (int64_t)0, tmp125);
CreateIdentity11( (int32_t)128, tmp20, tmp126);
Conv2DCSFMain( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)32,  (int32_t)3,  (int32_t)3,  (int32_t)128,  (int32_t)1,  (int32_t)1,  (int32_t)1,  (int32_t)1, tmp112, tmp124, tmp127,  (int64_t)15);
MatAddBroadCast4( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)128, tmp127, tmp126, tmp128);
ElemWiseMul4( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)128, tmp128, tmp128, tmp129,  (int64_t)15);
ScalarMul4( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)128,  (int64_t)327, tmp129, tmp130,  (int64_t)15);
ScalarMul4( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)128,  (int64_t)3276, tmp128, tmp131,  (int64_t)15);
MatAdd4( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)128, tmp130, tmp131, tmp132);
Concat2T444( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)256,  (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)128, tmp122,  (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)128, tmp132,  (int32_t)3, tmp133);
CreateTensor4( (int32_t)1,  (int32_t)1,  (int32_t)256,  (int32_t)32,  (int64_t)327, tmp134);
CreateIdentity44( (int32_t)1,  (int32_t)1,  (int32_t)256,  (int32_t)32, tmp21, tmp135);
CreateTensor1( (int32_t)32,  (int64_t)0, tmp136);
CreateIdentity11( (int32_t)32, tmp22, tmp137);
Conv2DCSFMain( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)256,  (int32_t)1,  (int32_t)1,  (int32_t)32,  (int32_t)0,  (int32_t)0,  (int32_t)1,  (int32_t)1, tmp133, tmp135, tmp138,  (int64_t)15);
MatAddBroadCast4( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)32, tmp138, tmp137, tmp139);
ElemWiseMul4( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)32, tmp139, tmp139, tmp140,  (int64_t)15);
ScalarMul4( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)32,  (int64_t)327, tmp140, tmp141,  (int64_t)15);
ScalarMul4( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)32,  (int64_t)3276, tmp139, tmp142,  (int64_t)15);
MatAdd4( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)32, tmp141, tmp142, tmp143);
CreateTensor4( (int32_t)1,  (int32_t)1,  (int32_t)32,  (int32_t)128,  (int64_t)327, tmp144);
CreateIdentity44( (int32_t)1,  (int32_t)1,  (int32_t)32,  (int32_t)128, tmp23, tmp145);
CreateTensor1( (int32_t)128,  (int64_t)0, tmp146);
CreateIdentity11( (int32_t)128, tmp24, tmp147);
Conv2DCSFMain( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)32,  (int32_t)1,  (int32_t)1,  (int32_t)128,  (int32_t)0,  (int32_t)0,  (int32_t)1,  (int32_t)1, tmp143, tmp145, tmp148,  (int64_t)15);
MatAddBroadCast4( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)128, tmp148, tmp147, tmp149);
ElemWiseMul4( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)128, tmp149, tmp149, tmp150,  (int64_t)15);
ScalarMul4( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)128,  (int64_t)327, tmp150, tmp151,  (int64_t)15);
ScalarMul4( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)128,  (int64_t)3276, tmp149, tmp152,  (int64_t)15);
MatAdd4( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)128, tmp151, tmp152, tmp153);
CreateTensor4( (int32_t)3,  (int32_t)3,  (int32_t)32,  (int32_t)128,  (int64_t)327, tmp154);
CreateIdentity44( (int32_t)3,  (int32_t)3,  (int32_t)32,  (int32_t)128, tmp25, tmp155);
CreateTensor1( (int32_t)128,  (int64_t)0, tmp156);
CreateIdentity11( (int32_t)128, tmp26, tmp157);
Conv2DCSFMain( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)32,  (int32_t)3,  (int32_t)3,  (int32_t)128,  (int32_t)1,  (int32_t)1,  (int32_t)1,  (int32_t)1, tmp143, tmp155, tmp158,  (int64_t)15);
MatAddBroadCast4( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)128, tmp158, tmp157, tmp159);
ElemWiseMul4( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)128, tmp159, tmp159, tmp160,  (int64_t)15);
ScalarMul4( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)128,  (int64_t)327, tmp160, tmp161,  (int64_t)15);
ScalarMul4( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)128,  (int64_t)3276, tmp159, tmp162,  (int64_t)15);
MatAdd4( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)128, tmp161, tmp162, tmp163);
Concat2T444( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)256,  (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)128, tmp153,  (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)128, tmp163,  (int32_t)3, tmp164);
CreateTensor4( (int32_t)1,  (int32_t)1,  (int32_t)256,  (int32_t)10,  (int64_t)327, tmp165);
CreateIdentity44( (int32_t)1,  (int32_t)1,  (int32_t)256,  (int32_t)10, tmp27, tmp166);
CreateTensor1( (int32_t)10,  (int64_t)0, tmp167);
CreateIdentity11( (int32_t)10, tmp28, tmp168);
Conv2DCSFMain( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)256,  (int32_t)1,  (int32_t)1,  (int32_t)10,  (int32_t)0,  (int32_t)0,  (int32_t)1,  (int32_t)1, tmp164, tmp166, tmp169,  (int64_t)15);
MatAddBroadCast4( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)10, tmp169, tmp168, tmp170);
ElemWiseMul4( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)10, tmp170, tmp170, tmp171,  (int64_t)15);
ScalarMul4( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)10,  (int64_t)327, tmp171, tmp172,  (int64_t)15);
ScalarMul4( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)10,  (int64_t)3276, tmp170, tmp173,  (int64_t)15);
MatAdd4( (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)10, tmp172, tmp173, tmp174);
AvgPool44( (int32_t)1,  (int32_t)1,  (int32_t)1,  (int32_t)10,  (int32_t)8,  (int32_t)8,  (int32_t)0,  (int32_t)0,  (int32_t)1,  (int32_t)1,  (int32_t)1,  (int32_t)8,  (int32_t)8,  (int32_t)10, tmp174, tmp175);
tmp176[ (int64_t)0] =  (int32_t)-1;
tmp176[ (int64_t)1] =  (int32_t)10;
i0 =  (int64_t)0;
i1 =  (int64_t)0;
i2 =  (int64_t)0;
i3 =  (int64_t)0;
for (uint32_t i4 =  (int32_t)0; i4 <  (int32_t)1; i4++){
for (uint32_t i5 =  (int32_t)0; i5 <  (int32_t)10; i5++){
tmp177[i4][i5] = tmp175[i0][i1][i2][i3];
i3 = (i3 +  (int64_t)1);
if ((i3 ==  (int64_t)10)) {
i3 =  (int64_t)0;
i2 = (i2 +  (int64_t)1);
if ((i2 ==  (int64_t)1)) {
i2 =  (int64_t)0;
i1 = (i1 +  (int64_t)1);
if ((i1 ==  (int64_t)1)) {
i1 =  (int64_t)0;
i0 = (i0 +  (int64_t)1);
}
}
}
}
}
ArgMax1( (int32_t)1,  (int32_t)1,  (int32_t)10, tmp177,  (int32_t)1, tmp178);
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)1; i0++){
add_to_output_queue(out_q, funcReconstruct2PCCons(tmp178[i0], 1), SERVER, cout);
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
