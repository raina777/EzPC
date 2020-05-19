#include "globals.h"
#ifdef F_SQUEEZENET_IMAGENET


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
//start_m();


auto tmp53 = make_vector<uint64_t>( (int32_t)1,  (int32_t)113,  (int32_t)113,  (int32_t)64);

auto tmp54 = make_vector<uint64_t>( (int32_t)1,  (int32_t)113,  (int32_t)113,  (int32_t)64);

auto tmp55 = make_vector<uint64_t>( (int32_t)1,  (int32_t)113,  (int32_t)113,  (int32_t)64);

auto tmp56 = make_vector<uint64_t>( (int32_t)1,  (int32_t)56,  (int32_t)56,  (int32_t)64);

auto tmp57 = make_vector<uint64_t>( (int32_t)1,  (int32_t)56,  (int32_t)56,  (int32_t)16);

auto tmp58 = make_vector<uint64_t>( (int32_t)1,  (int32_t)56,  (int32_t)56,  (int32_t)16);

auto tmp59 = make_vector<uint64_t>( (int32_t)1,  (int32_t)56,  (int32_t)56,  (int32_t)16);

auto tmp60 = make_vector<uint64_t>( (int32_t)1,  (int32_t)56,  (int32_t)56,  (int32_t)64);

auto tmp61 = make_vector<uint64_t>( (int32_t)1,  (int32_t)56,  (int32_t)56,  (int32_t)64);

auto tmp62 = make_vector<uint64_t>( (int32_t)1,  (int32_t)56,  (int32_t)56,  (int32_t)64);

auto tmp63 = make_vector<uint64_t>( (int32_t)1,  (int32_t)56,  (int32_t)56,  (int32_t)64);

auto tmp64 = make_vector<uint64_t>( (int32_t)1,  (int32_t)56,  (int32_t)56,  (int32_t)64);

auto tmp65 = make_vector<uint64_t>( (int32_t)1,  (int32_t)56,  (int32_t)56,  (int32_t)64);

auto tmp66 = make_vector<uint64_t>( (int32_t)1,  (int32_t)56,  (int32_t)56,  (int32_t)128);

auto tmp67 = make_vector<uint64_t>( (int32_t)1,  (int32_t)56,  (int32_t)56,  (int32_t)16);

auto tmp68 = make_vector<uint64_t>( (int32_t)1,  (int32_t)56,  (int32_t)56,  (int32_t)16);

auto tmp69 = make_vector<uint64_t>( (int32_t)1,  (int32_t)56,  (int32_t)56,  (int32_t)16);

auto tmp70 = make_vector<uint64_t>( (int32_t)1,  (int32_t)56,  (int32_t)56,  (int32_t)64);

auto tmp71 = make_vector<uint64_t>( (int32_t)1,  (int32_t)56,  (int32_t)56,  (int32_t)64);

auto tmp72 = make_vector<uint64_t>( (int32_t)1,  (int32_t)56,  (int32_t)56,  (int32_t)64);

auto tmp73 = make_vector<uint64_t>( (int32_t)1,  (int32_t)56,  (int32_t)56,  (int32_t)64);

auto tmp74 = make_vector<uint64_t>( (int32_t)1,  (int32_t)56,  (int32_t)56,  (int32_t)64);

auto tmp75 = make_vector<uint64_t>( (int32_t)1,  (int32_t)56,  (int32_t)56,  (int32_t)64);

auto tmp76 = make_vector<uint64_t>( (int32_t)1,  (int32_t)56,  (int32_t)56,  (int32_t)128);

auto tmp77 = make_vector<uint64_t>( (int32_t)1,  (int32_t)27,  (int32_t)27,  (int32_t)128);

auto tmp78 = make_vector<uint64_t>( (int32_t)1,  (int32_t)27,  (int32_t)27,  (int32_t)32);

auto tmp79 = make_vector<uint64_t>( (int32_t)1,  (int32_t)27,  (int32_t)27,  (int32_t)32);

auto tmp80 = make_vector<uint64_t>( (int32_t)1,  (int32_t)27,  (int32_t)27,  (int32_t)32);

auto tmp81 = make_vector<uint64_t>( (int32_t)1,  (int32_t)27,  (int32_t)27,  (int32_t)128);

auto tmp82 = make_vector<uint64_t>( (int32_t)1,  (int32_t)27,  (int32_t)27,  (int32_t)128);

auto tmp83 = make_vector<uint64_t>( (int32_t)1,  (int32_t)27,  (int32_t)27,  (int32_t)128);

auto tmp84 = make_vector<uint64_t>( (int32_t)1,  (int32_t)27,  (int32_t)27,  (int32_t)128);

auto tmp85 = make_vector<uint64_t>( (int32_t)1,  (int32_t)27,  (int32_t)27,  (int32_t)128);

auto tmp86 = make_vector<uint64_t>( (int32_t)1,  (int32_t)27,  (int32_t)27,  (int32_t)128);

auto tmp87 = make_vector<uint64_t>( (int32_t)1,  (int32_t)27,  (int32_t)27,  (int32_t)256);

auto tmp88 = make_vector<uint64_t>( (int32_t)1,  (int32_t)27,  (int32_t)27,  (int32_t)32);

auto tmp89 = make_vector<uint64_t>( (int32_t)1,  (int32_t)27,  (int32_t)27,  (int32_t)32);

auto tmp90 = make_vector<uint64_t>( (int32_t)1,  (int32_t)27,  (int32_t)27,  (int32_t)32);

auto tmp91 = make_vector<uint64_t>( (int32_t)1,  (int32_t)27,  (int32_t)27,  (int32_t)128);

auto tmp92 = make_vector<uint64_t>( (int32_t)1,  (int32_t)27,  (int32_t)27,  (int32_t)128);

auto tmp93 = make_vector<uint64_t>( (int32_t)1,  (int32_t)27,  (int32_t)27,  (int32_t)128);

auto tmp94 = make_vector<uint64_t>( (int32_t)1,  (int32_t)27,  (int32_t)27,  (int32_t)128);

auto tmp95 = make_vector<uint64_t>( (int32_t)1,  (int32_t)27,  (int32_t)27,  (int32_t)128);

auto tmp96 = make_vector<uint64_t>( (int32_t)1,  (int32_t)27,  (int32_t)27,  (int32_t)128);

auto tmp97 = make_vector<uint64_t>( (int32_t)1,  (int32_t)27,  (int32_t)27,  (int32_t)256);

auto tmp98 = make_vector<uint64_t>( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)256);

auto tmp99 = make_vector<uint64_t>( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)48);

auto tmp100 = make_vector<uint64_t>( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)48);

auto tmp101 = make_vector<uint64_t>( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)48);

auto tmp102 = make_vector<uint64_t>( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)192);

auto tmp103 = make_vector<uint64_t>( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)192);

auto tmp104 = make_vector<uint64_t>( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)192);

auto tmp105 = make_vector<uint64_t>( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)192);

auto tmp106 = make_vector<uint64_t>( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)192);

auto tmp107 = make_vector<uint64_t>( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)192);

auto tmp108 = make_vector<uint64_t>( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)384);

auto tmp109 = make_vector<uint64_t>( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)48);

auto tmp110 = make_vector<uint64_t>( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)48);

auto tmp111 = make_vector<uint64_t>( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)48);

auto tmp112 = make_vector<uint64_t>( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)192);

auto tmp113 = make_vector<uint64_t>( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)192);

auto tmp114 = make_vector<uint64_t>( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)192);

auto tmp115 = make_vector<uint64_t>( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)192);

auto tmp116 = make_vector<uint64_t>( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)192);

auto tmp117 = make_vector<uint64_t>( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)192);

auto tmp118 = make_vector<uint64_t>( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)384);

auto tmp119 = make_vector<uint64_t>( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)64);

auto tmp120 = make_vector<uint64_t>( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)64);

auto tmp121 = make_vector<uint64_t>( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)64);

auto tmp122 = make_vector<uint64_t>( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)256);

auto tmp123 = make_vector<uint64_t>( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)256);

auto tmp124 = make_vector<uint64_t>( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)256);

auto tmp125 = make_vector<uint64_t>( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)256);

auto tmp126 = make_vector<uint64_t>( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)256);

auto tmp127 = make_vector<uint64_t>( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)256);

auto tmp128 = make_vector<uint64_t>( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)512);

auto tmp129 = make_vector<uint64_t>( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)64);

auto tmp130 = make_vector<uint64_t>( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)64);

auto tmp131 = make_vector<uint64_t>( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)64);

auto tmp132 = make_vector<uint64_t>( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)256);

auto tmp133 = make_vector<uint64_t>( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)256);

auto tmp134 = make_vector<uint64_t>( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)256);

auto tmp135 = make_vector<uint64_t>( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)256);

auto tmp136 = make_vector<uint64_t>( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)256);

auto tmp137 = make_vector<uint64_t>( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)256);

auto tmp138 = make_vector<uint64_t>( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)512);

auto tmp139 = make_vector<uint64_t>( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)1000);

auto tmp140 = make_vector<uint64_t>( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)1000);

auto tmp141 = make_vector<uint64_t>( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)1000);

auto tmp142 = make_vector<uint64_t>( (int32_t)1,  (int32_t)1,  (int32_t)1,  (int32_t)1000);

auto tmp143 = make_vector<uint64_t>( (int32_t)1,  (int32_t)1,  (int32_t)1);

auto tmp0 = make_vector<uint64_t>( (int32_t)1,  (int32_t)227,  (int32_t)227,  (int32_t)3);
/* Variable to read the clear value corresponding to the input variable tmp0 at (473,1-473,47) */
uint64_t __tmp_in_tmp0;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)1; i0++){
for (uint32_t i1 =  (uint32_t)0; i1 <  (int32_t)227; i1++){
for (uint32_t i2 =  (uint32_t)0; i2 <  (int32_t)227; i2++){
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
/* Variable to read the clear value corresponding to the input variable tmp1 at (476,1-476,44) */
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
/* Variable to read the clear value corresponding to the input variable tmp2 at (479,1-479,35) */
uint64_t __tmp_in_tmp2;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)64; i0++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp2;
}
tmp2[i0] = (role == CLIENT) ? __tmp_in_tmp2 : 0;
}

auto tmp3 = make_vector<uint64_t>( (int32_t)1,  (int32_t)1,  (int32_t)64,  (int32_t)16);
/* Variable to read the clear value corresponding to the input variable tmp3 at (482,1-482,45) */
uint64_t __tmp_in_tmp3;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)1; i0++){
for (uint32_t i1 =  (uint32_t)0; i1 <  (int32_t)1; i1++){
for (uint32_t i2 =  (uint32_t)0; i2 <  (int32_t)64; i2++){
for (uint32_t i3 =  (uint32_t)0; i3 <  (int32_t)16; i3++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp3;
}
tmp3[i0][i1][i2][i3] = (role == CLIENT) ? __tmp_in_tmp3 : 0;
}
}
}
}

auto tmp4 = make_vector<uint64_t>( (int32_t)16);
/* Variable to read the clear value corresponding to the input variable tmp4 at (485,1-485,35) */
uint64_t __tmp_in_tmp4;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)16; i0++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp4;
}
tmp4[i0] = (role == CLIENT) ? __tmp_in_tmp4 : 0;
}

auto tmp5 = make_vector<uint64_t>( (int32_t)1,  (int32_t)1,  (int32_t)16,  (int32_t)64);
/* Variable to read the clear value corresponding to the input variable tmp5 at (488,1-488,45) */
uint64_t __tmp_in_tmp5;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)1; i0++){
for (uint32_t i1 =  (uint32_t)0; i1 <  (int32_t)1; i1++){
for (uint32_t i2 =  (uint32_t)0; i2 <  (int32_t)16; i2++){
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
/* Variable to read the clear value corresponding to the input variable tmp6 at (491,1-491,35) */
uint64_t __tmp_in_tmp6;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)64; i0++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp6;
}
tmp6[i0] = (role == CLIENT) ? __tmp_in_tmp6 : 0;
}

auto tmp7 = make_vector<uint64_t>( (int32_t)3,  (int32_t)3,  (int32_t)16,  (int32_t)64);
/* Variable to read the clear value corresponding to the input variable tmp7 at (494,1-494,45) */
uint64_t __tmp_in_tmp7;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)3; i0++){
for (uint32_t i1 =  (uint32_t)0; i1 <  (int32_t)3; i1++){
for (uint32_t i2 =  (uint32_t)0; i2 <  (int32_t)16; i2++){
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
/* Variable to read the clear value corresponding to the input variable tmp8 at (497,1-497,35) */
uint64_t __tmp_in_tmp8;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)64; i0++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp8;
}
tmp8[i0] = (role == CLIENT) ? __tmp_in_tmp8 : 0;
}

auto tmp9 = make_vector<uint64_t>( (int32_t)1,  (int32_t)1,  (int32_t)128,  (int32_t)16);
/* Variable to read the clear value corresponding to the input variable tmp9 at (500,1-500,46) */
uint64_t __tmp_in_tmp9;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)1; i0++){
for (uint32_t i1 =  (uint32_t)0; i1 <  (int32_t)1; i1++){
for (uint32_t i2 =  (uint32_t)0; i2 <  (int32_t)128; i2++){
for (uint32_t i3 =  (uint32_t)0; i3 <  (int32_t)16; i3++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp9;
}
tmp9[i0][i1][i2][i3] = (role == CLIENT) ? __tmp_in_tmp9 : 0;
}
}
}
}

auto tmp10 = make_vector<uint64_t>( (int32_t)16);
/* Variable to read the clear value corresponding to the input variable tmp10 at (503,1-503,36) */
uint64_t __tmp_in_tmp10;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)16; i0++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp10;
}
tmp10[i0] = (role == CLIENT) ? __tmp_in_tmp10 : 0;
}

auto tmp11 = make_vector<uint64_t>( (int32_t)1,  (int32_t)1,  (int32_t)16,  (int32_t)64);
/* Variable to read the clear value corresponding to the input variable tmp11 at (506,1-506,46) */
uint64_t __tmp_in_tmp11;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)1; i0++){
for (uint32_t i1 =  (uint32_t)0; i1 <  (int32_t)1; i1++){
for (uint32_t i2 =  (uint32_t)0; i2 <  (int32_t)16; i2++){
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
/* Variable to read the clear value corresponding to the input variable tmp12 at (509,1-509,36) */
uint64_t __tmp_in_tmp12;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)64; i0++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp12;
}
tmp12[i0] = (role == CLIENT) ? __tmp_in_tmp12 : 0;
}

auto tmp13 = make_vector<uint64_t>( (int32_t)3,  (int32_t)3,  (int32_t)16,  (int32_t)64);
/* Variable to read the clear value corresponding to the input variable tmp13 at (512,1-512,46) */
uint64_t __tmp_in_tmp13;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)3; i0++){
for (uint32_t i1 =  (uint32_t)0; i1 <  (int32_t)3; i1++){
for (uint32_t i2 =  (uint32_t)0; i2 <  (int32_t)16; i2++){
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
/* Variable to read the clear value corresponding to the input variable tmp14 at (515,1-515,36) */
uint64_t __tmp_in_tmp14;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)64; i0++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp14;
}
tmp14[i0] = (role == CLIENT) ? __tmp_in_tmp14 : 0;
}

auto tmp15 = make_vector<uint64_t>( (int32_t)1,  (int32_t)1,  (int32_t)128,  (int32_t)32);
/* Variable to read the clear value corresponding to the input variable tmp15 at (518,1-518,47) */
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
/* Variable to read the clear value corresponding to the input variable tmp16 at (521,1-521,36) */
uint64_t __tmp_in_tmp16;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)32; i0++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp16;
}
tmp16[i0] = (role == CLIENT) ? __tmp_in_tmp16 : 0;
}

auto tmp17 = make_vector<uint64_t>( (int32_t)1,  (int32_t)1,  (int32_t)32,  (int32_t)128);
/* Variable to read the clear value corresponding to the input variable tmp17 at (524,1-524,47) */
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
/* Variable to read the clear value corresponding to the input variable tmp18 at (527,1-527,37) */
uint64_t __tmp_in_tmp18;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)128; i0++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp18;
}
tmp18[i0] = (role == CLIENT) ? __tmp_in_tmp18 : 0;
}

auto tmp19 = make_vector<uint64_t>( (int32_t)3,  (int32_t)3,  (int32_t)32,  (int32_t)128);
/* Variable to read the clear value corresponding to the input variable tmp19 at (530,1-530,47) */
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
/* Variable to read the clear value corresponding to the input variable tmp20 at (533,1-533,37) */
uint64_t __tmp_in_tmp20;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)128; i0++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp20;
}
tmp20[i0] = (role == CLIENT) ? __tmp_in_tmp20 : 0;
}

auto tmp21 = make_vector<uint64_t>( (int32_t)1,  (int32_t)1,  (int32_t)256,  (int32_t)32);
/* Variable to read the clear value corresponding to the input variable tmp21 at (536,1-536,47) */
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
/* Variable to read the clear value corresponding to the input variable tmp22 at (539,1-539,36) */
uint64_t __tmp_in_tmp22;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)32; i0++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp22;
}
tmp22[i0] = (role == CLIENT) ? __tmp_in_tmp22 : 0;
}

auto tmp23 = make_vector<uint64_t>( (int32_t)1,  (int32_t)1,  (int32_t)32,  (int32_t)128);
/* Variable to read the clear value corresponding to the input variable tmp23 at (542,1-542,47) */
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
/* Variable to read the clear value corresponding to the input variable tmp24 at (545,1-545,37) */
uint64_t __tmp_in_tmp24;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)128; i0++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp24;
}
tmp24[i0] = (role == CLIENT) ? __tmp_in_tmp24 : 0;
}

auto tmp25 = make_vector<uint64_t>( (int32_t)3,  (int32_t)3,  (int32_t)32,  (int32_t)128);
/* Variable to read the clear value corresponding to the input variable tmp25 at (548,1-548,47) */
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
/* Variable to read the clear value corresponding to the input variable tmp26 at (551,1-551,37) */
uint64_t __tmp_in_tmp26;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)128; i0++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp26;
}
tmp26[i0] = (role == CLIENT) ? __tmp_in_tmp26 : 0;
}

auto tmp27 = make_vector<uint64_t>( (int32_t)1,  (int32_t)1,  (int32_t)256,  (int32_t)48);
/* Variable to read the clear value corresponding to the input variable tmp27 at (554,1-554,47) */
uint64_t __tmp_in_tmp27;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)1; i0++){
for (uint32_t i1 =  (uint32_t)0; i1 <  (int32_t)1; i1++){
for (uint32_t i2 =  (uint32_t)0; i2 <  (int32_t)256; i2++){
for (uint32_t i3 =  (uint32_t)0; i3 <  (int32_t)48; i3++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp27;
}
tmp27[i0][i1][i2][i3] = (role == CLIENT) ? __tmp_in_tmp27 : 0;
}
}
}
}

auto tmp28 = make_vector<uint64_t>( (int32_t)48);
/* Variable to read the clear value corresponding to the input variable tmp28 at (557,1-557,36) */
uint64_t __tmp_in_tmp28;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)48; i0++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp28;
}
tmp28[i0] = (role == CLIENT) ? __tmp_in_tmp28 : 0;
}

auto tmp29 = make_vector<uint64_t>( (int32_t)1,  (int32_t)1,  (int32_t)48,  (int32_t)192);
/* Variable to read the clear value corresponding to the input variable tmp29 at (560,1-560,47) */
uint64_t __tmp_in_tmp29;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)1; i0++){
for (uint32_t i1 =  (uint32_t)0; i1 <  (int32_t)1; i1++){
for (uint32_t i2 =  (uint32_t)0; i2 <  (int32_t)48; i2++){
for (uint32_t i3 =  (uint32_t)0; i3 <  (int32_t)192; i3++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp29;
}
tmp29[i0][i1][i2][i3] = (role == CLIENT) ? __tmp_in_tmp29 : 0;
}
}
}
}

auto tmp30 = make_vector<uint64_t>( (int32_t)192);
/* Variable to read the clear value corresponding to the input variable tmp30 at (563,1-563,37) */
uint64_t __tmp_in_tmp30;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)192; i0++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp30;
}
tmp30[i0] = (role == CLIENT) ? __tmp_in_tmp30 : 0;
}

auto tmp31 = make_vector<uint64_t>( (int32_t)3,  (int32_t)3,  (int32_t)48,  (int32_t)192);
/* Variable to read the clear value corresponding to the input variable tmp31 at (566,1-566,47) */
uint64_t __tmp_in_tmp31;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)3; i0++){
for (uint32_t i1 =  (uint32_t)0; i1 <  (int32_t)3; i1++){
for (uint32_t i2 =  (uint32_t)0; i2 <  (int32_t)48; i2++){
for (uint32_t i3 =  (uint32_t)0; i3 <  (int32_t)192; i3++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp31;
}
tmp31[i0][i1][i2][i3] = (role == CLIENT) ? __tmp_in_tmp31 : 0;
}
}
}
}

auto tmp32 = make_vector<uint64_t>( (int32_t)192);
/* Variable to read the clear value corresponding to the input variable tmp32 at (569,1-569,37) */
uint64_t __tmp_in_tmp32;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)192; i0++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp32;
}
tmp32[i0] = (role == CLIENT) ? __tmp_in_tmp32 : 0;
}

auto tmp33 = make_vector<uint64_t>( (int32_t)1,  (int32_t)1,  (int32_t)384,  (int32_t)48);
/* Variable to read the clear value corresponding to the input variable tmp33 at (572,1-572,47) */
uint64_t __tmp_in_tmp33;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)1; i0++){
for (uint32_t i1 =  (uint32_t)0; i1 <  (int32_t)1; i1++){
for (uint32_t i2 =  (uint32_t)0; i2 <  (int32_t)384; i2++){
for (uint32_t i3 =  (uint32_t)0; i3 <  (int32_t)48; i3++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp33;
}
tmp33[i0][i1][i2][i3] = (role == CLIENT) ? __tmp_in_tmp33 : 0;
}
}
}
}

auto tmp34 = make_vector<uint64_t>( (int32_t)48);
/* Variable to read the clear value corresponding to the input variable tmp34 at (575,1-575,36) */
uint64_t __tmp_in_tmp34;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)48; i0++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp34;
}
tmp34[i0] = (role == CLIENT) ? __tmp_in_tmp34 : 0;
}

auto tmp35 = make_vector<uint64_t>( (int32_t)1,  (int32_t)1,  (int32_t)48,  (int32_t)192);
/* Variable to read the clear value corresponding to the input variable tmp35 at (578,1-578,47) */
uint64_t __tmp_in_tmp35;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)1; i0++){
for (uint32_t i1 =  (uint32_t)0; i1 <  (int32_t)1; i1++){
for (uint32_t i2 =  (uint32_t)0; i2 <  (int32_t)48; i2++){
for (uint32_t i3 =  (uint32_t)0; i3 <  (int32_t)192; i3++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp35;
}
tmp35[i0][i1][i2][i3] = (role == CLIENT) ? __tmp_in_tmp35 : 0;
}
}
}
}

auto tmp36 = make_vector<uint64_t>( (int32_t)192);
/* Variable to read the clear value corresponding to the input variable tmp36 at (581,1-581,37) */
uint64_t __tmp_in_tmp36;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)192; i0++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp36;
}
tmp36[i0] = (role == CLIENT) ? __tmp_in_tmp36 : 0;
}

auto tmp37 = make_vector<uint64_t>( (int32_t)3,  (int32_t)3,  (int32_t)48,  (int32_t)192);
/* Variable to read the clear value corresponding to the input variable tmp37 at (584,1-584,47) */
uint64_t __tmp_in_tmp37;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)3; i0++){
for (uint32_t i1 =  (uint32_t)0; i1 <  (int32_t)3; i1++){
for (uint32_t i2 =  (uint32_t)0; i2 <  (int32_t)48; i2++){
for (uint32_t i3 =  (uint32_t)0; i3 <  (int32_t)192; i3++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp37;
}
tmp37[i0][i1][i2][i3] = (role == CLIENT) ? __tmp_in_tmp37 : 0;
}
}
}
}

auto tmp38 = make_vector<uint64_t>( (int32_t)192);
/* Variable to read the clear value corresponding to the input variable tmp38 at (587,1-587,37) */
uint64_t __tmp_in_tmp38;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)192; i0++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp38;
}
tmp38[i0] = (role == CLIENT) ? __tmp_in_tmp38 : 0;
}

auto tmp39 = make_vector<uint64_t>( (int32_t)1,  (int32_t)1,  (int32_t)384,  (int32_t)64);
/* Variable to read the clear value corresponding to the input variable tmp39 at (590,1-590,47) */
uint64_t __tmp_in_tmp39;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)1; i0++){
for (uint32_t i1 =  (uint32_t)0; i1 <  (int32_t)1; i1++){
for (uint32_t i2 =  (uint32_t)0; i2 <  (int32_t)384; i2++){
for (uint32_t i3 =  (uint32_t)0; i3 <  (int32_t)64; i3++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp39;
}
tmp39[i0][i1][i2][i3] = (role == CLIENT) ? __tmp_in_tmp39 : 0;
}
}
}
}

auto tmp40 = make_vector<uint64_t>( (int32_t)64);
/* Variable to read the clear value corresponding to the input variable tmp40 at (593,1-593,36) */
uint64_t __tmp_in_tmp40;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)64; i0++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp40;
}
tmp40[i0] = (role == CLIENT) ? __tmp_in_tmp40 : 0;
}

auto tmp41 = make_vector<uint64_t>( (int32_t)1,  (int32_t)1,  (int32_t)64,  (int32_t)256);
/* Variable to read the clear value corresponding to the input variable tmp41 at (596,1-596,47) */
uint64_t __tmp_in_tmp41;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)1; i0++){
for (uint32_t i1 =  (uint32_t)0; i1 <  (int32_t)1; i1++){
for (uint32_t i2 =  (uint32_t)0; i2 <  (int32_t)64; i2++){
for (uint32_t i3 =  (uint32_t)0; i3 <  (int32_t)256; i3++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp41;
}
tmp41[i0][i1][i2][i3] = (role == CLIENT) ? __tmp_in_tmp41 : 0;
}
}
}
}

auto tmp42 = make_vector<uint64_t>( (int32_t)256);
/* Variable to read the clear value corresponding to the input variable tmp42 at (599,1-599,37) */
uint64_t __tmp_in_tmp42;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)256; i0++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp42;
}
tmp42[i0] = (role == CLIENT) ? __tmp_in_tmp42 : 0;
}

auto tmp43 = make_vector<uint64_t>( (int32_t)3,  (int32_t)3,  (int32_t)64,  (int32_t)256);
/* Variable to read the clear value corresponding to the input variable tmp43 at (602,1-602,47) */
uint64_t __tmp_in_tmp43;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)3; i0++){
for (uint32_t i1 =  (uint32_t)0; i1 <  (int32_t)3; i1++){
for (uint32_t i2 =  (uint32_t)0; i2 <  (int32_t)64; i2++){
for (uint32_t i3 =  (uint32_t)0; i3 <  (int32_t)256; i3++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp43;
}
tmp43[i0][i1][i2][i3] = (role == CLIENT) ? __tmp_in_tmp43 : 0;
}
}
}
}

auto tmp44 = make_vector<uint64_t>( (int32_t)256);
/* Variable to read the clear value corresponding to the input variable tmp44 at (605,1-605,37) */
uint64_t __tmp_in_tmp44;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)256; i0++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp44;
}
tmp44[i0] = (role == CLIENT) ? __tmp_in_tmp44 : 0;
}

auto tmp45 = make_vector<uint64_t>( (int32_t)1,  (int32_t)1,  (int32_t)512,  (int32_t)64);
/* Variable to read the clear value corresponding to the input variable tmp45 at (608,1-608,47) */
uint64_t __tmp_in_tmp45;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)1; i0++){
for (uint32_t i1 =  (uint32_t)0; i1 <  (int32_t)1; i1++){
for (uint32_t i2 =  (uint32_t)0; i2 <  (int32_t)512; i2++){
for (uint32_t i3 =  (uint32_t)0; i3 <  (int32_t)64; i3++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp45;
}
tmp45[i0][i1][i2][i3] = (role == CLIENT) ? __tmp_in_tmp45 : 0;
}
}
}
}

auto tmp46 = make_vector<uint64_t>( (int32_t)64);
/* Variable to read the clear value corresponding to the input variable tmp46 at (611,1-611,36) */
uint64_t __tmp_in_tmp46;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)64; i0++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp46;
}
tmp46[i0] = (role == CLIENT) ? __tmp_in_tmp46 : 0;
}

auto tmp47 = make_vector<uint64_t>( (int32_t)1,  (int32_t)1,  (int32_t)64,  (int32_t)256);
/* Variable to read the clear value corresponding to the input variable tmp47 at (614,1-614,47) */
uint64_t __tmp_in_tmp47;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)1; i0++){
for (uint32_t i1 =  (uint32_t)0; i1 <  (int32_t)1; i1++){
for (uint32_t i2 =  (uint32_t)0; i2 <  (int32_t)64; i2++){
for (uint32_t i3 =  (uint32_t)0; i3 <  (int32_t)256; i3++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp47;
}
tmp47[i0][i1][i2][i3] = (role == CLIENT) ? __tmp_in_tmp47 : 0;
}
}
}
}

auto tmp48 = make_vector<uint64_t>( (int32_t)256);
/* Variable to read the clear value corresponding to the input variable tmp48 at (617,1-617,37) */
uint64_t __tmp_in_tmp48;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)256; i0++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp48;
}
tmp48[i0] = (role == CLIENT) ? __tmp_in_tmp48 : 0;
}

auto tmp49 = make_vector<uint64_t>( (int32_t)3,  (int32_t)3,  (int32_t)64,  (int32_t)256);
/* Variable to read the clear value corresponding to the input variable tmp49 at (620,1-620,47) */
uint64_t __tmp_in_tmp49;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)3; i0++){
for (uint32_t i1 =  (uint32_t)0; i1 <  (int32_t)3; i1++){
for (uint32_t i2 =  (uint32_t)0; i2 <  (int32_t)64; i2++){
for (uint32_t i3 =  (uint32_t)0; i3 <  (int32_t)256; i3++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp49;
}
tmp49[i0][i1][i2][i3] = (role == CLIENT) ? __tmp_in_tmp49 : 0;
}
}
}
}

auto tmp50 = make_vector<uint64_t>( (int32_t)256);
/* Variable to read the clear value corresponding to the input variable tmp50 at (623,1-623,37) */
uint64_t __tmp_in_tmp50;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)256; i0++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp50;
}
tmp50[i0] = (role == CLIENT) ? __tmp_in_tmp50 : 0;
}

auto tmp51 = make_vector<uint64_t>( (int32_t)1,  (int32_t)1,  (int32_t)512,  (int32_t)1000);
/* Variable to read the clear value corresponding to the input variable tmp51 at (626,1-626,49) */
uint64_t __tmp_in_tmp51;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)1; i0++){
for (uint32_t i1 =  (uint32_t)0; i1 <  (int32_t)1; i1++){
for (uint32_t i2 =  (uint32_t)0; i2 <  (int32_t)512; i2++){
for (uint32_t i3 =  (uint32_t)0; i3 <  (int32_t)1000; i3++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp51;
}
tmp51[i0][i1][i2][i3] = (role == CLIENT) ? __tmp_in_tmp51 : 0;
}
}
}
}

auto tmp52 = make_vector<uint64_t>( (int32_t)1000);
/* Variable to read the clear value corresponding to the input variable tmp52 at (629,1-629,38) */
uint64_t __tmp_in_tmp52;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)1000; i0++){
if ((role == CLIENT)) {
cin >> __tmp_in_tmp52;
}
tmp52[i0] = (role == CLIENT) ? __tmp_in_tmp52 : 0;
}

cout<<"Starting 2nd syncronize .. "<<endl;
synchronize(2000000); 
cout<<"Syncronized .. now starting actual execution at "<<getCurrentTime()<<endl;
start_m();

Conv2DCSFMain( (int32_t)1,  (int32_t)227,  (int32_t)227,  (int32_t)3,  (int32_t)3,  (int32_t)3,  (int32_t)64,  (int32_t)0,  (int32_t)0,  (int32_t)2,  (int32_t)2, tmp0, tmp1, tmp53,  (int64_t)15);
MatAddBroadCast4( (int32_t)1,  (int32_t)113,  (int32_t)113,  (int32_t)64, tmp53, tmp2, tmp54);
Relu4( (int32_t)1,  (int32_t)113,  (int32_t)113,  (int32_t)64, tmp54, tmp55);
MaxPool44( (int32_t)1,  (int32_t)56,  (int32_t)56,  (int32_t)64,  (int32_t)3,  (int32_t)3,  (int32_t)0,  (int32_t)0,  (int32_t)2,  (int32_t)2,  (int32_t)1,  (int32_t)113,  (int32_t)113,  (int32_t)64, tmp55, tmp56);
Conv2DCSFMain( (int32_t)1,  (int32_t)56,  (int32_t)56,  (int32_t)64,  (int32_t)1,  (int32_t)1,  (int32_t)16,  (int32_t)0,  (int32_t)0,  (int32_t)1,  (int32_t)1, tmp56, tmp3, tmp57,  (int64_t)15);
MatAddBroadCast4( (int32_t)1,  (int32_t)56,  (int32_t)56,  (int32_t)16, tmp57, tmp4, tmp58);
Relu4( (int32_t)1,  (int32_t)56,  (int32_t)56,  (int32_t)16, tmp58, tmp59);
Conv2DCSFMain( (int32_t)1,  (int32_t)56,  (int32_t)56,  (int32_t)16,  (int32_t)1,  (int32_t)1,  (int32_t)64,  (int32_t)0,  (int32_t)0,  (int32_t)1,  (int32_t)1, tmp59, tmp5, tmp60,  (int64_t)15);
MatAddBroadCast4( (int32_t)1,  (int32_t)56,  (int32_t)56,  (int32_t)64, tmp60, tmp6, tmp61);
Relu4( (int32_t)1,  (int32_t)56,  (int32_t)56,  (int32_t)64, tmp61, tmp62);
Conv2DCSFMain( (int32_t)1,  (int32_t)56,  (int32_t)56,  (int32_t)16,  (int32_t)3,  (int32_t)3,  (int32_t)64,  (int32_t)1,  (int32_t)1,  (int32_t)1,  (int32_t)1, tmp59, tmp7, tmp63,  (int64_t)15);
MatAddBroadCast4( (int32_t)1,  (int32_t)56,  (int32_t)56,  (int32_t)64, tmp63, tmp8, tmp64);
Relu4( (int32_t)1,  (int32_t)56,  (int32_t)56,  (int32_t)64, tmp64, tmp65);
Concat2T444( (int32_t)1,  (int32_t)56,  (int32_t)56,  (int32_t)128,  (int32_t)1,  (int32_t)56,  (int32_t)56,  (int32_t)64, tmp62,  (int32_t)1,  (int32_t)56,  (int32_t)56,  (int32_t)64, tmp65,  (int32_t)3, tmp66);
Conv2DCSFMain( (int32_t)1,  (int32_t)56,  (int32_t)56,  (int32_t)128,  (int32_t)1,  (int32_t)1,  (int32_t)16,  (int32_t)0,  (int32_t)0,  (int32_t)1,  (int32_t)1, tmp66, tmp9, tmp67,  (int64_t)15);
MatAddBroadCast4( (int32_t)1,  (int32_t)56,  (int32_t)56,  (int32_t)16, tmp67, tmp10, tmp68);
Relu4( (int32_t)1,  (int32_t)56,  (int32_t)56,  (int32_t)16, tmp68, tmp69);
Conv2DCSFMain( (int32_t)1,  (int32_t)56,  (int32_t)56,  (int32_t)16,  (int32_t)1,  (int32_t)1,  (int32_t)64,  (int32_t)0,  (int32_t)0,  (int32_t)1,  (int32_t)1, tmp69, tmp11, tmp70,  (int64_t)15);
MatAddBroadCast4( (int32_t)1,  (int32_t)56,  (int32_t)56,  (int32_t)64, tmp70, tmp12, tmp71);
Relu4( (int32_t)1,  (int32_t)56,  (int32_t)56,  (int32_t)64, tmp71, tmp72);
Conv2DCSFMain( (int32_t)1,  (int32_t)56,  (int32_t)56,  (int32_t)16,  (int32_t)3,  (int32_t)3,  (int32_t)64,  (int32_t)1,  (int32_t)1,  (int32_t)1,  (int32_t)1, tmp69, tmp13, tmp73,  (int64_t)15);
MatAddBroadCast4( (int32_t)1,  (int32_t)56,  (int32_t)56,  (int32_t)64, tmp73, tmp14, tmp74);
Relu4( (int32_t)1,  (int32_t)56,  (int32_t)56,  (int32_t)64, tmp74, tmp75);
Concat2T444( (int32_t)1,  (int32_t)56,  (int32_t)56,  (int32_t)128,  (int32_t)1,  (int32_t)56,  (int32_t)56,  (int32_t)64, tmp72,  (int32_t)1,  (int32_t)56,  (int32_t)56,  (int32_t)64, tmp75,  (int32_t)3, tmp76);
MaxPool44( (int32_t)1,  (int32_t)27,  (int32_t)27,  (int32_t)128,  (int32_t)3,  (int32_t)3,  (int32_t)0,  (int32_t)0,  (int32_t)2,  (int32_t)2,  (int32_t)1,  (int32_t)56,  (int32_t)56,  (int32_t)128, tmp76, tmp77);
Conv2DCSFMain( (int32_t)1,  (int32_t)27,  (int32_t)27,  (int32_t)128,  (int32_t)1,  (int32_t)1,  (int32_t)32,  (int32_t)0,  (int32_t)0,  (int32_t)1,  (int32_t)1, tmp77, tmp15, tmp78,  (int64_t)15);
MatAddBroadCast4( (int32_t)1,  (int32_t)27,  (int32_t)27,  (int32_t)32, tmp78, tmp16, tmp79);
Relu4( (int32_t)1,  (int32_t)27,  (int32_t)27,  (int32_t)32, tmp79, tmp80);
Conv2DCSFMain( (int32_t)1,  (int32_t)27,  (int32_t)27,  (int32_t)32,  (int32_t)1,  (int32_t)1,  (int32_t)128,  (int32_t)0,  (int32_t)0,  (int32_t)1,  (int32_t)1, tmp80, tmp17, tmp81,  (int64_t)15);
MatAddBroadCast4( (int32_t)1,  (int32_t)27,  (int32_t)27,  (int32_t)128, tmp81, tmp18, tmp82);
Relu4( (int32_t)1,  (int32_t)27,  (int32_t)27,  (int32_t)128, tmp82, tmp83);
Conv2DCSFMain( (int32_t)1,  (int32_t)27,  (int32_t)27,  (int32_t)32,  (int32_t)3,  (int32_t)3,  (int32_t)128,  (int32_t)1,  (int32_t)1,  (int32_t)1,  (int32_t)1, tmp80, tmp19, tmp84,  (int64_t)15);
MatAddBroadCast4( (int32_t)1,  (int32_t)27,  (int32_t)27,  (int32_t)128, tmp84, tmp20, tmp85);
Relu4( (int32_t)1,  (int32_t)27,  (int32_t)27,  (int32_t)128, tmp85, tmp86);
Concat2T444( (int32_t)1,  (int32_t)27,  (int32_t)27,  (int32_t)256,  (int32_t)1,  (int32_t)27,  (int32_t)27,  (int32_t)128, tmp83,  (int32_t)1,  (int32_t)27,  (int32_t)27,  (int32_t)128, tmp86,  (int32_t)3, tmp87);
Conv2DCSFMain( (int32_t)1,  (int32_t)27,  (int32_t)27,  (int32_t)256,  (int32_t)1,  (int32_t)1,  (int32_t)32,  (int32_t)0,  (int32_t)0,  (int32_t)1,  (int32_t)1, tmp87, tmp21, tmp88,  (int64_t)15);
MatAddBroadCast4( (int32_t)1,  (int32_t)27,  (int32_t)27,  (int32_t)32, tmp88, tmp22, tmp89);
Relu4( (int32_t)1,  (int32_t)27,  (int32_t)27,  (int32_t)32, tmp89, tmp90);
Conv2DCSFMain( (int32_t)1,  (int32_t)27,  (int32_t)27,  (int32_t)32,  (int32_t)1,  (int32_t)1,  (int32_t)128,  (int32_t)0,  (int32_t)0,  (int32_t)1,  (int32_t)1, tmp90, tmp23, tmp91,  (int64_t)15);
MatAddBroadCast4( (int32_t)1,  (int32_t)27,  (int32_t)27,  (int32_t)128, tmp91, tmp24, tmp92);
Relu4( (int32_t)1,  (int32_t)27,  (int32_t)27,  (int32_t)128, tmp92, tmp93);
Conv2DCSFMain( (int32_t)1,  (int32_t)27,  (int32_t)27,  (int32_t)32,  (int32_t)3,  (int32_t)3,  (int32_t)128,  (int32_t)1,  (int32_t)1,  (int32_t)1,  (int32_t)1, tmp90, tmp25, tmp94,  (int64_t)15);
MatAddBroadCast4( (int32_t)1,  (int32_t)27,  (int32_t)27,  (int32_t)128, tmp94, tmp26, tmp95);
Relu4( (int32_t)1,  (int32_t)27,  (int32_t)27,  (int32_t)128, tmp95, tmp96);
Concat2T444( (int32_t)1,  (int32_t)27,  (int32_t)27,  (int32_t)256,  (int32_t)1,  (int32_t)27,  (int32_t)27,  (int32_t)128, tmp93,  (int32_t)1,  (int32_t)27,  (int32_t)27,  (int32_t)128, tmp96,  (int32_t)3, tmp97);
MaxPool44( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)256,  (int32_t)3,  (int32_t)3,  (int32_t)0,  (int32_t)0,  (int32_t)2,  (int32_t)2,  (int32_t)1,  (int32_t)27,  (int32_t)27,  (int32_t)256, tmp97, tmp98);
Conv2DCSFMain( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)256,  (int32_t)1,  (int32_t)1,  (int32_t)48,  (int32_t)0,  (int32_t)0,  (int32_t)1,  (int32_t)1, tmp98, tmp27, tmp99,  (int64_t)15);
MatAddBroadCast4( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)48, tmp99, tmp28, tmp100);
Relu4( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)48, tmp100, tmp101);
Conv2DCSFMain( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)48,  (int32_t)1,  (int32_t)1,  (int32_t)192,  (int32_t)0,  (int32_t)0,  (int32_t)1,  (int32_t)1, tmp101, tmp29, tmp102,  (int64_t)15);
MatAddBroadCast4( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)192, tmp102, tmp30, tmp103);
Relu4( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)192, tmp103, tmp104);
Conv2DCSFMain( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)48,  (int32_t)3,  (int32_t)3,  (int32_t)192,  (int32_t)1,  (int32_t)1,  (int32_t)1,  (int32_t)1, tmp101, tmp31, tmp105,  (int64_t)15);
MatAddBroadCast4( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)192, tmp105, tmp32, tmp106);
Relu4( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)192, tmp106, tmp107);
Concat2T444( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)384,  (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)192, tmp104,  (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)192, tmp107,  (int32_t)3, tmp108);
Conv2DCSFMain( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)384,  (int32_t)1,  (int32_t)1,  (int32_t)48,  (int32_t)0,  (int32_t)0,  (int32_t)1,  (int32_t)1, tmp108, tmp33, tmp109,  (int64_t)15);
MatAddBroadCast4( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)48, tmp109, tmp34, tmp110);
Relu4( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)48, tmp110, tmp111);
Conv2DCSFMain( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)48,  (int32_t)1,  (int32_t)1,  (int32_t)192,  (int32_t)0,  (int32_t)0,  (int32_t)1,  (int32_t)1, tmp111, tmp35, tmp112,  (int64_t)15);
MatAddBroadCast4( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)192, tmp112, tmp36, tmp113);
Relu4( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)192, tmp113, tmp114);
Conv2DCSFMain( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)48,  (int32_t)3,  (int32_t)3,  (int32_t)192,  (int32_t)1,  (int32_t)1,  (int32_t)1,  (int32_t)1, tmp111, tmp37, tmp115,  (int64_t)15);
MatAddBroadCast4( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)192, tmp115, tmp38, tmp116);
Relu4( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)192, tmp116, tmp117);
Concat2T444( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)384,  (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)192, tmp114,  (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)192, tmp117,  (int32_t)3, tmp118);
Conv2DCSFMain( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)384,  (int32_t)1,  (int32_t)1,  (int32_t)64,  (int32_t)0,  (int32_t)0,  (int32_t)1,  (int32_t)1, tmp118, tmp39, tmp119,  (int64_t)15);
MatAddBroadCast4( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)64, tmp119, tmp40, tmp120);
Relu4( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)64, tmp120, tmp121);
Conv2DCSFMain( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)64,  (int32_t)1,  (int32_t)1,  (int32_t)256,  (int32_t)0,  (int32_t)0,  (int32_t)1,  (int32_t)1, tmp121, tmp41, tmp122,  (int64_t)15);
MatAddBroadCast4( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)256, tmp122, tmp42, tmp123);
Relu4( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)256, tmp123, tmp124);
Conv2DCSFMain( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)64,  (int32_t)3,  (int32_t)3,  (int32_t)256,  (int32_t)1,  (int32_t)1,  (int32_t)1,  (int32_t)1, tmp121, tmp43, tmp125,  (int64_t)15);
MatAddBroadCast4( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)256, tmp125, tmp44, tmp126);
Relu4( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)256, tmp126, tmp127);
Concat2T444( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)512,  (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)256, tmp124,  (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)256, tmp127,  (int32_t)3, tmp128);
Conv2DCSFMain( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)512,  (int32_t)1,  (int32_t)1,  (int32_t)64,  (int32_t)0,  (int32_t)0,  (int32_t)1,  (int32_t)1, tmp128, tmp45, tmp129,  (int64_t)15);
MatAddBroadCast4( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)64, tmp129, tmp46, tmp130);
Relu4( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)64, tmp130, tmp131);
Conv2DCSFMain( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)64,  (int32_t)1,  (int32_t)1,  (int32_t)256,  (int32_t)0,  (int32_t)0,  (int32_t)1,  (int32_t)1, tmp131, tmp47, tmp132,  (int64_t)15);
MatAddBroadCast4( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)256, tmp132, tmp48, tmp133);
Relu4( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)256, tmp133, tmp134);
Conv2DCSFMain( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)64,  (int32_t)3,  (int32_t)3,  (int32_t)256,  (int32_t)1,  (int32_t)1,  (int32_t)1,  (int32_t)1, tmp131, tmp49, tmp135,  (int64_t)15);
MatAddBroadCast4( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)256, tmp135, tmp50, tmp136);
Relu4( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)256, tmp136, tmp137);
Concat2T444( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)512,  (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)256, tmp134,  (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)256, tmp137,  (int32_t)3, tmp138);
Conv2DCSFMain( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)512,  (int32_t)1,  (int32_t)1,  (int32_t)1000,  (int32_t)0,  (int32_t)0,  (int32_t)1,  (int32_t)1, tmp138, tmp51, tmp139,  (int64_t)15);
MatAddBroadCast4( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)1000, tmp139, tmp52, tmp140);
Relu4( (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)1000, tmp140, tmp141);
AvgPool44( (int32_t)1,  (int32_t)1,  (int32_t)1,  (int32_t)1000,  (int32_t)13,  (int32_t)13,  (int32_t)0,  (int32_t)0,  (int32_t)1,  (int32_t)1,  (int32_t)1,  (int32_t)13,  (int32_t)13,  (int32_t)1000, tmp141, tmp142);
ArgMax3( (int32_t)1,  (int32_t)1,  (int32_t)1,  (int32_t)1,  (int32_t)1,  (int32_t)1,  (int32_t)1000, tmp142,  (int32_t)3, tmp143);
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)1; i0++){
for (uint32_t i1 =  (uint32_t)0; i1 <  (int32_t)1; i1++){
for (uint32_t i2 =  (uint32_t)0; i2 <  (int32_t)1; i2++){
add_to_output_queue(out_q, funcReconstruct2PCCons(tmp143[i0][i1][i2], 1), SERVER, cout);
}
}
}
end_m(whichNetwork);

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
