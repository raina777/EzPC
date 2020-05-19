
#ifndef AESOBJECT_H
#define AESOBJECT_H

#pragma once
#include <algorithm>
#include "globals.h"


/*
* NOTE : This class is not thread safe. 
		 Donot allow mutiple threads to read from this class at once.
		 TODO_niskum : Do we need to fix this ?
*/
class AESObject
{
private:
	//AES variables
	__m128i pseudoRandomString[RANDOM_COMPUTE];
	__m128i tempSecComp[RANDOM_COMPUTE];
	unsigned long rCounter = -1;
	AES_KEY_TED aes_key;

	//Extraction variables
	__m128i randomBitNumber {0};
	uint8_t randomBitCounter = 0;
	__m128i random8BitNumber {0};
	uint8_t random8BitCounter = 0; 
	__m128i random64BitNumber {0};
	bool fetch64New = true;

	//Private extraction functions
	__m128i newRandomNumber();

	//Private helper functions
	smallType AES_random(int i);

#ifdef PRECOMPUTEAES
	__m128i* preComputedKeys = NULL;
	__m128i* tempKeyArray = NULL;
	void PreComputeKeysFunc(uint64_t startKeyNum, uint64_t numKeys);
#endif

public:
	//Constructor
	AESObject(char* filename);
	~AESObject();

#ifdef PRECOMPUTEAES
	void PreComputeKeys(uint64_t numKeysToPrecompute, int32_t numThreads);
#endif
	
	//Randomness functions
	myType get64Bits();
	smallType get8Bits();
	smallType getBit();
	
	//Other randomness functions
	smallType randModPrime();
	smallType randNonZeroModPrime();
	myType randModuloOdd();
	void AES_random_shuffle(vector<smallType> &vec, size_t begin_offset, size_t end_offset);

	unsigned long getRCounter()
	{
		return rCounter;
	}

#ifdef PRECOMPUTEAES

	__m128i* getPreComputedKeysPtr()
	{
		return preComputedKeys;
	}
	
	void fillWithRandomBits64(uint64_t* arr, size_t size);
	void fillWithRandomBits8(uint8_t* arr, size_t size);
	void fillWithRandomModuloPrimeBits(uint8_t* arr, size_t size);
	void fillWithRandomModuloOddBits(uint64_t* arr, size_t size);
	
#endif
};



#endif