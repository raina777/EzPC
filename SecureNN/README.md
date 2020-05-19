Packages required 
(Ubuntu): 
* sudo apt-get install libssl-dev
* sudo apt-get install g++
* sudo apt-get install make

(Mac, not working):
* brew install libssl
* cd /usr/local/include
* ln -s ../opt/openssl/include/openssl . 
Also, use the second build command for Mac to be safe. Read the Caveats part of brew install instructions, since Apple has deprecated the use of OpenSSL.

To-Do: AESObjects are used as global variables. They need to be passed through all the functions involved.  

Note: libmiracl.a is compiled locally, if it does not work try the libmiracl_old.a or download the source files from https://github.com/miracl/MIRACL.git and compile miracl.a yourself (and rename to libmriacl.a when copying into this repo)

Note: when extracting bits out of \_\_m128i, val[0] corresponds to the LSB. The set functions of \_\_m128i take the first argument as val[0].

Note: If myType != uint64_t, myType multiplication function won't work.

Note: Matrix multiplication assembly code only works for Intel C/C++ compiler. The non-assembly code has correctness issues when both the multiplicands are large uint64_t's

Note: Check the security of funcPrivateCompare4PC, specifically if(r[index2] == MINUS_ONE)

Note: Valgrind usage
* Install valgrind `sudo apt-get install valgrind`
* Edit make file flags to -g -O0 instead of -O3
* Run `make clean; make`
* Run standalone code `valgrind --tool=memcheck --leak-check=full --track-origins=yes --dsymutil=yes ./BMRPassive.out STANDALONE 4 files/parties_localhost files/keyA files/keyAB files/data/mnist_data_8_samples files/data/mnist_labels_8_samples files/data/mnist_data_8_samples files/data/mnist_labels_8_samples`

### version1.0: Oakland submission code

### version2.0: Oakland submission code 
* bug fixes in standalone code
* modified makefile
* Small testing parameter values
* Including small dataset of 8 samples

### version3.0: 4PC code with new implemetation
* Working testing and debugging functions
* AES precompute not done in this version

### version4.0: Complete 3PC code with new implemetation
* AES precompute not done in this version 

### version5.0: CNN implemetation
* Layer Types added using object polymorphism
* CNN implemented
* AES precompute done
* Non blocking sockets implemented
* Gazelle, MiniONN, SecureML, Chameleon networks added

### version6.0: Beaver Optimization
* Beaver Triplet optimization
* ReLU, MaxPool optimization