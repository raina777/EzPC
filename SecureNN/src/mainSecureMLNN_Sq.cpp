#include "globals.h"
#ifdef F_SECUREMLNN_SQ
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


auto tmp7 = make_vector<int64_t>( (int32_t)784,  (int32_t)128);

auto tmp8 = make_vector<uint64_t>( (int32_t)784,  (int32_t)128);

auto tmp9 = make_vector<int64_t>( (int32_t)128);

auto tmp10 = make_vector<uint64_t>( (int32_t)128);

auto tmp11 = make_vector<uint64_t>( (int32_t)1,  (int32_t)128);

auto tmp12 = make_vector<uint64_t>( (int32_t)1,  (int32_t)128);

auto tmp13 = make_vector<uint64_t>( (int32_t)1,  (int32_t)128);

auto tmp14 = make_vector<int64_t>( (int32_t)128,  (int32_t)128);

auto tmp15 = make_vector<uint64_t>( (int32_t)128,  (int32_t)128);

auto tmp16 = make_vector<int64_t>( (int32_t)128);

auto tmp17 = make_vector<uint64_t>( (int32_t)128);

auto tmp18 = make_vector<uint64_t>( (int32_t)1,  (int32_t)128);

auto tmp19 = make_vector<uint64_t>( (int32_t)1,  (int32_t)128);

auto tmp20 = make_vector<uint64_t>( (int32_t)1,  (int32_t)128);

auto tmp21 = make_vector<int64_t>( (int32_t)128,  (int32_t)10);

auto tmp22 = make_vector<uint64_t>( (int32_t)128,  (int32_t)10);

auto tmp23 = make_vector<int64_t>( (int32_t)10);

auto tmp24 = make_vector<uint64_t>( (int32_t)10);

auto tmp25 = make_vector<uint64_t>( (int32_t)1,  (int32_t)10);

auto tmp26 = make_vector<uint64_t>( (int32_t)1,  (int32_t)10);

auto tmp27 = make_vector<uint64_t>( (int32_t)1);

auto tmp0 = make_vector<uint64_t>( (int32_t)1,  (int32_t)784);
/* Variable to read the clear value corresponding to the input variable tmp0 at (408,1-408,39) */
uint64_t __tmp_in_tmp0;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)1; i0++){
for (uint32_t i1 =  (uint32_t)0; i1 <  (int32_t)784; i1++){
if ((role == CLIENT)) {
//cin >> __tmp_in_tmp0;
}
tmp0[i0][i1] = (role == CLIENT) ? __tmp_in_tmp0 : 0;
}
}

auto tmp1 = make_vector<uint64_t>( (int32_t)784,  (int32_t)128);
/* Variable to read the clear value corresponding to the input variable tmp1 at (411,1-411,41) */
uint64_t __tmp_in_tmp1;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)784; i0++){
for (uint32_t i1 =  (uint32_t)0; i1 <  (int32_t)128; i1++){
if ((role == CLIENT)) {
//cin >> __tmp_in_tmp1;
}
tmp1[i0][i1] = (role == CLIENT) ? __tmp_in_tmp1 : 0;
}
}

auto tmp2 = make_vector<uint64_t>( (int32_t)128);
/* Variable to read the clear value corresponding to the input variable tmp2 at (414,1-414,36) */
uint64_t __tmp_in_tmp2;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)128; i0++){
if ((role == CLIENT)) {
//cin >> __tmp_in_tmp2;
}
tmp2[i0] = (role == CLIENT) ? __tmp_in_tmp2 : 0;
}

auto tmp3 = make_vector<uint64_t>( (int32_t)128,  (int32_t)128);
/* Variable to read the clear value corresponding to the input variable tmp3 at (417,1-417,41) */
uint64_t __tmp_in_tmp3;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)128; i0++){
for (uint32_t i1 =  (uint32_t)0; i1 <  (int32_t)128; i1++){
if ((role == CLIENT)) {
//cin >> __tmp_in_tmp3;
}
tmp3[i0][i1] = (role == CLIENT) ? __tmp_in_tmp3 : 0;
}
}

auto tmp4 = make_vector<uint64_t>( (int32_t)128);
/* Variable to read the clear value corresponding to the input variable tmp4 at (420,1-420,36) */
uint64_t __tmp_in_tmp4;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)128; i0++){
if ((role == CLIENT)) {
//cin >> __tmp_in_tmp4;
}
tmp4[i0] = (role == CLIENT) ? __tmp_in_tmp4 : 0;
}

auto tmp5 = make_vector<uint64_t>( (int32_t)128,  (int32_t)10);
/* Variable to read the clear value corresponding to the input variable tmp5 at (423,1-423,40) */
uint64_t __tmp_in_tmp5;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)128; i0++){
for (uint32_t i1 =  (uint32_t)0; i1 <  (int32_t)10; i1++){
if ((role == CLIENT)) {
//cin >> __tmp_in_tmp5;
}
tmp5[i0][i1] = (role == CLIENT) ? __tmp_in_tmp5 : 0;
}
}

auto tmp6 = make_vector<uint64_t>( (int32_t)10);
/* Variable to read the clear value corresponding to the input variable tmp6 at (426,1-426,35) */
uint64_t __tmp_in_tmp6;
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)10; i0++){
if ((role == CLIENT)) {
//cin >> __tmp_in_tmp6;
}
tmp6[i0] = (role == CLIENT) ? __tmp_in_tmp6 : 0;
}

cout<<"Starting 2nd syncronize .. "<<endl;
synchronize(2000000); 
cout<<"Syncronized .. now starting actual execution at "<<getCurrentTime()<<endl;
start_m();

CreateTensor2( (int32_t)784,  (int32_t)128,  (int64_t)3276, tmp7);
CreateIdentity22( (int32_t)784,  (int32_t)128, tmp1, tmp8);
CreateTensor1( (int32_t)128,  (int64_t)3276, tmp9);
CreateIdentity11( (int32_t)128, tmp2, tmp10);
MatMulCSF2D( (int32_t)1,  (int32_t)784,  (int32_t)128, tmp0, tmp8, tmp11,  (int64_t)15);
MatAddBroadCast2( (int32_t)1,  (int32_t)128, tmp11, tmp10, tmp12);
ElemWiseMul2( (int32_t)1,  (int32_t)128, tmp12, tmp12, tmp13,  (int64_t)15);
CreateTensor2( (int32_t)128,  (int32_t)128,  (int64_t)3276, tmp14);
CreateIdentity22( (int32_t)128,  (int32_t)128, tmp3, tmp15);
CreateTensor1( (int32_t)128,  (int64_t)3276, tmp16);
CreateIdentity11( (int32_t)128, tmp4, tmp17);
MatMulCSF2D( (int32_t)1,  (int32_t)128,  (int32_t)128, tmp13, tmp15, tmp18,  (int64_t)15);
MatAddBroadCast2( (int32_t)1,  (int32_t)128, tmp18, tmp17, tmp19);
ElemWiseMul2( (int32_t)1,  (int32_t)128, tmp19, tmp19, tmp20,  (int64_t)15);
CreateTensor2( (int32_t)128,  (int32_t)10,  (int64_t)3276, tmp21);
CreateIdentity22( (int32_t)128,  (int32_t)10, tmp5, tmp22);
CreateTensor1( (int32_t)10,  (int64_t)3276, tmp23);
CreateIdentity11( (int32_t)10, tmp6, tmp24);
MatMulCSF2D( (int32_t)1,  (int32_t)128,  (int32_t)10, tmp20, tmp22, tmp25,  (int64_t)15);
MatAddBroadCast2( (int32_t)1,  (int32_t)10, tmp25, tmp24, tmp26);
ArgMax1( (int32_t)1,  (int32_t)1,  (int32_t)10, tmp26,  (int32_t)1, tmp27);
for (uint32_t i0 =  (uint32_t)0; i0 <  (int32_t)1; i0++){
add_to_output_queue(out_q, funcReconstruct2PCCons(tmp27[i0], 1), SERVER, cout);
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
