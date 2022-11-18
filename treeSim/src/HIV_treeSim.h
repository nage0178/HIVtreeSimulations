#ifndef HIV_HEADER
#define HIV_HEADER
#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <assert.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <regex.h>
#include <getopt.h>

#define CalcProb(rates, probs, totRate) { \
  totRate  = rates[0] +  rates[1] + rates[2] + rates[3] + rates[4] + rates[5] + rates[6] + rates[7] + rates[8] + rates[9] + rates[10];\
  probs[0] = rates[0] / totRate; /*uninfectedBirth*/\
  probs[1] = rates[1] / totRate; /*uninfect to infect*/\
  probs[2] = rates[2] / totRate; /*uninfectDeath*/\
  probs[3] = rates[3] / totRate; /*infect death*/\
  probs[4] = rates[4] / totRate; /*virus birth*/\
  probs[5] = rates[5] / totRate; /*virus death*/\
  probs[6] = rates[6] / totRate; /*uninfectTolatIncomp*/\
  probs[7] = rates[7] / totRate; /*uninfect to lat comp*/\
  probs[8] = rates[8] / totRate; /*incomp death*/\
  probs[9] = rates[9] / totRate; /*comp death*/\
  probs[10] = rates[10]/ totRate; /*reactivate*/\
  }


/* Bifurcating tree */
struct Node {
  float birth_time;  /* time of birth event at which virus arose */
  float sample_time; /* 0 if not yet sampled otherwise time of sampling */
  char isSample;       /* 2 is default, 1 if yes */
  struct Node* left;  /* NULL if no children*/
  struct Node* right; /* only viruses with no children are replicating */
  struct Node* parent;/* Pointer to the parent node, null for stem */

  char isLatent;         /* 0 if active, 1 if latent */
  float timeToLatent;  /* Time the sequence became latent, -l if not latent*/
  float totTimeLatent; /* Total time in a latent state on a branch*/
  long unsigned arrayNum; /* Used for saving to array*/

} ;

struct cellCountTimes {
  long int numVirus;
  long int numCellInfect;
  long int numCellUninfect;
  long int numLatentComp;
  long int numLatentIncomp;

  int numVirusSample;
  int numLatentSample;

  double totTime;
  double waitTime;
  long int sizeTree;
  long int sizeLatent;

  double mLBlood;
  double lambda;
  double kappa;
  double d;
  double delta;
  double p;
  double c;

  //MAKE SURE YOU DON"T DIVIDE THESE BY ML BLOOD
  double probLatent;       /* eta, probLatent */
  double reactLatent;       /* alpha, reactLatent */
  double probDefect;       /* gamma, probDefect */
  double latIncompDeath;    /* tau, latIncompDeath  */
  double latCompDeath; /* sigma, latCompDeath */

  long int maxActive;
  long int maxVirus;
  long int maxIncomp;
  long int maxComp;

  int numSampleTimes;

  double timeLastSample;
  int sampleCounter;

  int totSampleVirus;
  int totSampleLatent;

  unsigned long int numEvents;

  double ARTstart;

} ;

struct sampleTimesLine {
  double time;
  int numActive;
  int numLatent;
};


bool FloatEquals(double a, double b, double threshold);
struct Node* newNode(double birth_add, int latentState, struct Node* parent_ptr);
void removeVirus(struct Node* array[], long int arrayIndex, long int lastArrayIndex, long int* count);
void reactivateLatent(struct Node* latentArray[], struct Node* activeArray[], long int latentIndex, long int* latentLast, long int* activeIndex, double totTime);
void birthEvent(struct Node* parentArray[], struct Node* daughterArray[], long int parentIndex, long int* daughterIndex, double time, int latentState, long int* cellDecrease, unsigned long int* totMem);
void clearTreeMemory(struct Node* node_ptr, unsigned long int* mem);
void setArrayNull(struct Node** array, long int size);
void setVirusArrays(struct Node** activeArray, struct Node** virusArray, struct Node** latentIncompArray, struct Node** latentCompArray,
struct Node** virusSampleArray, struct Node** latentSampleArray, long int maxActive, long int maxVirus, long int maxIncomp,
long int maxComp, int totSampleVirus, int totSampleLatent);
void ClearMemory(struct Node* root_ptr, struct Node* stem, struct Node* active[], struct Node* latent[],
struct Node* latentDefect[], struct Node* virus[], struct Node* sampleVirus[], struct Node* sampleLatent[], unsigned long int* totMem);
void pruneTip(struct Node* tip, unsigned long int* totMem);
void printTree(struct Node* node);
void sample(struct Node* arrayToSample[], struct Node* sampleArray[], int pickVirus, int* sampleIndex, double totTime, long int* numVirus);

void setRates(
  double* prodInfectionRates, double* parameters, double reactLatent, double latIncompDeath, double latCompDeath,
  double* rate, long int numVirus, long int numCellInfect, long int numCellUninfect, long int numLatentComp, long int numLatentIncomp) ;


void setInfectionRates (double* rate, double* prodInfectionRates, long int numCellUninfect, long int numVirus);

void sampleEvent(gsl_rng *r, int* sampleCounter, int numSampleVirus[], int numSampleLatent[],
  struct Node* virusArray[], struct Node* virusSampleArray[], struct Node* latentSampleArray[],
  struct Node* latentIncompArray[], struct Node* latentCompArray[],
  int* numVirusSample, int* numLatentSample, long int* numVirus, long int* numLatentIncomp, long int* numLatentComp, double totTime);

struct Node* reset(struct Node* stem, struct Node* root_ptr, int* sampleCounter, unsigned long int* numEvents, unsigned long int* totMem,
struct Node** activeArray, struct Node** virusArray, struct Node** latentIncompArray, struct Node** latentCompArray,
struct Node** virusSampleArray, struct Node** latentSampleArray, long int maxActive, long int maxVirus, long int maxIncomp,
long int maxComp, int totSampleVirus, int totSampleLatent);
void resetCounts(long int* numVirus, long int numVirusInit, long int* numCellUninfect, long int numCellUninfectInit,
    long int* numCellInfect, long int numCellInfectInit, long int* numLatentComp, long int numLatentCompInit,
    long int* numLatentIncomp, long int numLatentIncompInit, int* numLatentSample, int* numVirusSample, int* sampleCounter, double* totTime);


int findFileLength(char* filename);
void readSampleTimes (double sampleTimes[], int numSampleVirus[], int numSampleLatent[], char* filename, int fileLength);

void writeTxt(char* fileName, char* fileContents);
char* NewickSampled(struct Node* node, char* treeString, double previousTime, long int* nodeName) ;

void findSampledParent(struct Node* node, struct Node* stem);
void findSampledAll(struct Node** sampleActiveArray, struct Node** sampleLatentArray, int sizeActive, int sizeLatent, struct Node* stem);

char* NewickAll(struct Node* node, char* treeString, long int* nodeName);
void writeSeed(unsigned int RGSeed, char* seedFile);
char* LatentSampled(struct Node* node, char* latentString, double timeLatent, long int* nodeName) ;
void checkValidSampleTimes (double sampleTimes[], int numSampleVirus[], int numSampleLatent[], int fileLength);

struct Node** reallocArray(long int* sizeArrayOld, struct Node** array, unsigned long int* totMem);
struct Node** reallocArrayDecrease(long int* sizeArrayOld, struct Node** array, unsigned long int* totMem);

void parseControlFile(char* fileName, char* sampleTime, double* mLBlood,
double* parameters, double* probLatent, double* reactLatent, double* probDefect,
double* latIncompDeath, double* latCompDeath, double* ARTstart, unsigned int* seed, int* seedChange, int* volChange, int* sampleFileChange);


void findCellTypes(struct Node** activeArray, struct Node** virusArray, struct Node** latentIncompArray, struct Node** latentCompArray,
  struct Node** virusSampleArray, struct Node** latentSampleArray,
  long int numInfect, long int numVirus, long int numLatentIncomp, long int numLatentComp, int numVirusSample, int numLatentSample);
char* checkPointTree(struct Node* node, char* treeString, long int* nodeName) ;
char* checkPointLatent(struct Node* node, char* latentString, long int* nodeName);
int writeBinary(char* treeString, char * fileName);
int writeCountsBinary(struct cellCountTimes* counts, char * fileName);
struct Node* newNodeReload(double birthTime, int num, struct Node* parent) ;
int readInBinaries(char* countsFile, char* gslFile, char* latentFile, char* treeFile,
  struct cellCountTimes * counts, char ** treeFilePointer, char ** latentFilePointer, gsl_rng * r) ;
int makeTree(struct Node* node_ptr, char* treeString, int currChar, unsigned long int *totMem);
struct Node* newNodeSampledTree(struct Node* parent, unsigned long int* totMem);
void readLatent (struct Node* node, char* latentFile, struct Node** activeArray, struct Node** virusArray, struct Node** latentIncompArray, struct Node** latentCompArray,
struct Node** virusSampleArray, struct Node** latentSampleArray, char** ptr);

int writeSampleFileBinary(int sampleFileLength, double* sampleTimes, int* numSampleVirus, int* numSampleLatent, char * fileName);
int readSampleFileBinary(int sampleFileLength, double* sampleTimes, int* numSampleVirus, int* numSampleLatent,
char * fileName);
void writeGslState(gsl_rng *r, char* outfile);
void writeCounts(long int numVirus, long int numCellInfect, long int numCellUninfect, long int numLatentComp,
long int numLatentIncomp, int numVirusSample, int numLatentSample, double totTime, double waitTime, char* treeString, char* latentString,
double mLBlood, double* parameters, double probLatent, double reactLatent, double probDefect, double latIncompDeath,
double latCompDeath, int sampleFileLength, int sampleCounter, long int maxActive, long int maxVirus, long int maxIncomp,
long int maxComp, int totSampleVirus, int totSampleLatent, int timeLastSample, unsigned long int numEvents, double ARTstart, char* outfile);

unsigned long int findTreeStringMem(long int numVirus, long int numCellInfect, long int numLatentComp, long int numLatentIncomp, int numVirusSample, int numLatentSample, double timeLastSample) ;
unsigned long int findLatentStringMem(long int numVirus, long int numCellInfect, long int numLatentComp, long int numLatentIncomp, int numVirusSample, int numLatentSample, double timeLastSample);

void writeLog(char* outfile, double * parameters, int sampleFileLength, long int numVirusInit,
                long int numCellUninfectInit, long int numCellInfectInit, long int numLatentCompInit,
               long int numLatentIncompInit, double * sampleTimes, int * numSampleVirus,
               int * numSampleLatent, double mLBlood, double probLatent, double reactLatent, 
               double probDefect, double latIncompDeath, double latCompDeath, double lambdaInput, 
	       double kappaInput);
#endif
