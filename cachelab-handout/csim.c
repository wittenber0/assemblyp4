// anrus-rwwittenberg


#include "cachelab.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <unistd.h>
#include <math.h>

// State information for a cache line.
struct line {
	unsigned int tag;
	int valid;
	unsigned int timeStamp;
};

// Split the memory address into the cache fields.
struct cacheParam {
	unsigned int s;
	unsigned int b;
	unsigned int tag;
};

// Global variables.
int  numSetIndexBits;
int  numSets;
int  numLines;
int  blockSize;
int  verbosityFlag = 0;
char *trace;
int  hits = 0;
int  misses = 0;
int  globalTimeStamp = 0;
int  evictions = 0;
char line[256];

struct line **cachePtr;

void printUsage(char *argv[]) {
	printf("Usage: %s [-hv] -s <num> -E <num> -b <num> -t <file>", argv[0]);
	printf("Options:\n"
			"  -h         Print this help message.\n"
			"  -v         Optional verbose flag.\n"
			"  -s <num>   Number of set index bits.\n"
			"  -E <num>   Number of lines per set.\n"
			"  -b <num>   Number of block offset bits.\n"
			"  -t <file>  Trace file.\n\n");
	printf("Examples:\n"
			"  linux>  ./csim-ref -s 4 -E 1 -b 4 -t traces/yi.trace\n"
			"  linux>  ./csim-ref -v -s 8 -E 2 -b 4 -t traces/yi.trace\n");
}

// Allocate and initialize the cache.
void createCache() {
	cachePtr = malloc(numSets * sizeof(struct line *));

	for (int setIndex = 0; setIndex < numSets; setIndex ++) {
		cachePtr[setIndex] = malloc(numLines * sizeof(struct line));

		for (int lineIndex = 0; lineIndex < numLines; lineIndex++) {
			cachePtr[setIndex][lineIndex].valid = 0;
		}
	}
}

// When the trace is prefixed by an "L", then that means to try and load the memory value into the cache.
void loadOperation(struct cacheParam *cacheParamPtr) {
	struct line *setPtr = cachePtr[cacheParamPtr->s];

	// Check to see if the value is already in the cache.
	for (int lines = 0; lines < numLines; lines++) {
		if (setPtr[lines].valid && setPtr[lines].tag == cacheParamPtr->tag) {
			hits++;
			if (verbosityFlag) {
				printf("%s hit\n", &line[1]);
			}
			setPtr[lines].timeStamp = globalTimeStamp;
			return;
		}
	}
	misses++;
	if (verbosityFlag) {
		printf("%s miss\n", &line[1]);
	}

	// See if there is an unused line for our value.
	for (int lines = 0; lines < numLines; lines++) {
		if (!setPtr[lines].valid) {
			setPtr[lines].tag = cacheParamPtr->tag;
			setPtr[lines].timeStamp = globalTimeStamp;
			setPtr[lines].valid = 1;
			return;
		}
	}

	struct line *LRU = setPtr;
	// Find the oldest line and evict it to be replaced with our value.
	for (int lines = 0; lines < numLines; lines++) {
		if (setPtr[lines].timeStamp < LRU->timeStamp) {
			LRU = &setPtr[lines];
		}
	}
	LRU->tag = cacheParamPtr->tag;
	LRU->timeStamp = globalTimeStamp;
	LRU->valid = 1;
	evictions++;
	if (verbosityFlag) {
		printf("%s eviction\n", &line[1]);
	}

}

// In this assignment, load and store behave essentially the same, so we just call load operation, but we are really attempting
// to store the value in the cache.
void storeOperation(struct cacheParam *cacheParamPtr) {
	loadOperation(cacheParamPtr);
}

// In this assignment, we are calling both a load and store operation together, so we just call them individually.
void modifyOperation(struct cacheParam *cacheParamPtr) {
	loadOperation(cacheParamPtr);
	storeOperation(cacheParamPtr);
}


// In order to properly parse the memory address, we need to create a mask in which we specify the starting bit and the ending bit.
unsigned int getField(unsigned lowBit, unsigned highBit, unsigned int memAddr) {

	unsigned r = 0;
	for (unsigned i = lowBit; i <= highBit; i++)
		r |= 1 << i;

	unsigned result = r & memAddr;
	result >>= lowBit;
	return result;
}

// Parse each address and store the correct values in the parameter fields so that we may compare them with those in the cache.
void parseAddress(unsigned int memAddr, struct cacheParam *cacheParamPtr) {
	cacheParamPtr->tag = memAddr >> (blockSize + numSetIndexBits);
	cacheParamPtr->s = getField(blockSize, (numSetIndexBits + blockSize -1), memAddr);
	cacheParamPtr->b = getField(0, (blockSize -1), memAddr);
}

//  Read the traces from the specified file and run the corresponding operation.
void runSimulation() {
	FILE* file;
	file = fopen(trace, "r");

	if (file == NULL) {
		printf("Error could not open file.\n");
		exit(EXIT_FAILURE);
	}

	char operation;
	unsigned int memAddr;
	struct cacheParam cacheParam1;

	while (fgets(line, sizeof(line), file)) {
		if (line[0] == 'I') {
			continue;
		}

		strtok(line, "\n");
		sscanf(line, " %c %X", &operation, &memAddr);
		parseAddress(memAddr, &cacheParam1);

		switch(operation) {

		case 'L':
			loadOperation(&cacheParam1);
			break;
		case 'S':
			storeOperation(&cacheParam1);
			break;
		case 'M':
			modifyOperation(&cacheParam1);
			break;
		default:
			printf("Unknown Operation");
			break;
		}

		globalTimeStamp++;

	}

	fclose(file);
}


struct memBlock {
	int timeStamp;
	int startEndAddr;
};

// Read the command line arguments and assign the appropriate variables.
void getArgs(int argc, char *argv[]) {
	int opt;
	if (argc == 1) {
		printf("Missing required command line argument\n ");
		printUsage(argv);
		exit(1);
	}

	while ((opt = getopt (argc, argv, "s:E:b:t:vh")) != -1) {

		switch (opt) {

		case 's':
			numSetIndexBits = atoi(optarg);
			numSets = (int)pow((double)2, numSetIndexBits);
			break;
		case 'E':
			numLines = atoi(optarg);
			break;
		case 'b':
			blockSize = atoi(optarg);
			break;
		case 't':
			trace = optarg;
			break;
		case 'v':
			verbosityFlag = 1;
			break;
		case 'h':
		default:
			printUsage(argv);
			exit(0);
		}
	}
}

// Call our functions to read the arguments, create the cache, and run the simulation using the designated trace file.
int main(int argc, char *argv[]) {

	getArgs(argc, argv);

	createCache(numSetIndexBits, numLines, blockSize);

	runSimulation();

	// We could free the cache memory here, but exiting will free all the memory anyway.

	printSummary(hits, misses, evictions);
	return 0;
}
