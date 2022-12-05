/*
    Copyright (C) 2022 Anna Nagel

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include "HIV_treeSim.h"
#include <ctype.h>

#define BUF_LEN 200

/* From 21st Century C by Ben Klemens
Used to concatenate strings of unknown length */
#define Sasprintf(write_to,  ...) {                     \
    char *tmp_string_for_extend = (write_to);           \
    int returnval = asprintf(&(write_to), __VA_ARGS__); \
    if (returnval < 0) {                                \
      fprintf(stderr, "Problem with sasprintf\n");      \
      exit(1);                                          \
    }                                                   \
    free(tmp_string_for_extend);                        \
}


/* Checks if two doubles are with a threshold value of each other */
bool FloatEquals(double a, double b, double threshold) {
  return fabs(a-b) < threshold;
}


void printHelp() {
        printf("Program options:\n\n");
        printf("c: Name of control file. Does not have an effect if loading from a checkpoint. \n\n");
        printf("h: Prints the program options and exits.\n\n");
        printf("i: Input file with sampling times. In csv file format without a header. First column is the time, second is number of active sequences to sample, third is the number of latent sequences to sample. Does not have an effect if loading from a checkpoint.\n\n");
        printf("l: Loads checkpoint files and restarts the run. Requires the prefix of the checkpoint files.\n\n");
        printf("o: Output file prefix. \n\n");
        printf("p: Frequency of printing counts and memory usage. If 0, counts are not printed. Otherwise, every p events the counts are printed.\n\n");
        printf("r: Restart if population goes extinct. 0 if no, otherwise yes. \n\n");
        printf("s: Starting seed. Does not have an effect if loading from a checkpoint.. \n\n");
        printf("t: Time to checkpoint the program. The program cannot load a checkpoint and have another checkpoint later. \n\n");
        printf("v: Volume of the simulation in mL. Does not have an effect if loading from a checkpoint.\n\n");

	printf("For a more detail description, see https://github.com/nage0178/HIVtreeSimulations\n");

}

/* Create a new new node. Allocates memory for the node.
Assignments values to the variables in the node structure */
struct Node* newNode(double birth_add, int latentState, struct Node* parent_ptr) {
  /* Allocates memory for a new node */
  struct Node* node_ptr;
  node_ptr  = malloc(sizeof(struct Node));
  if (node_ptr == NULL) {
    fprintf(stderr, "Insufficient Memory\n");
    exit(1);
  }

  node_ptr->birth_time = birth_add;
  node_ptr->sample_time = 0;
  node_ptr->isSample = 2;

  node_ptr->isLatent = latentState;
  /* If latent, became latent at the same time the node was born */
  if (latentState == 1) {
    node_ptr->timeToLatent = birth_add;
  } else {
    node_ptr->timeToLatent = -1;
  }
  /* totTimeLatent is adjusted when a sequence becomes active */
  node_ptr->totTimeLatent = 0;

  /* Sets pointers to daughter nodes and parent node*/
  node_ptr->left = NULL;
  node_ptr->right = NULL;
  node_ptr->parent = parent_ptr;
  node_ptr->arrayNum = 0;

  return(node_ptr);
}


/* Removes a virus from the array */
/* arrayIndex is the element that is being removed,
arrayIndex is the last occupied element in the array,
count is the number of viruses in the given class of viruses*/
void removeVirus(struct Node* array[], long int arrayIndex, long int lastArrayIndex, long int* count) {
  /* Puts the last element in the array in the location of the element that is
  being removed. Sets the last occupied element in the array to null */
    array[arrayIndex] = array[lastArrayIndex];
    array[lastArrayIndex] = NULL;

    --(*count);
  return;
}

/* Takes a virus from the latent array and moves it to the active array */
void reactivateLatent(struct Node* latentArray[], struct Node* activeArray[], long int latentIndex,
  long int* latentLast, long int* activeIndex, double totTime){
     --(*latentLast);

    /* Sets node variables */
    latentArray[latentIndex]->isLatent = 0;
    latentArray[latentIndex]->totTimeLatent = latentArray[latentIndex]->totTimeLatent + (totTime - latentArray[latentIndex]->timeToLatent);
    latentArray[latentIndex]->timeToLatent = -1;

    /* Moves the virus to the other array */
    activeArray[*activeIndex] = latentArray[latentIndex];
    latentArray[latentIndex] = latentArray[*latentLast];
    latentArray[*latentLast] = NULL;

    ++(*activeIndex);

    return;
}

/* Creates a two new daughter nodes for a given parent */
void birthEvent(struct Node* parentArray[], struct Node* daughterArray[], long int parentIndex, long int* daughterIndex, double time, int latentState, long int* cellDecrease, unsigned long int* totMem){
  /* Makes new nodes for the daughter viruses*/
  /* The left daughter is the same cell/virus as the parent at time present (so it's in the same array
  at the parents). */
  parentArray[parentIndex]->left = newNode(time, 0, parentArray[parentIndex]); //Latent viruses can't be parents
  /* Daughter that is a different cell/virus type than the parent */
  parentArray[parentIndex]->right = newNode(time, latentState, parentArray[parentIndex]); //O if active, 1 if latent

  /* Right daughter is added to the daughter array */
  daughterArray[*daughterIndex] = parentArray[parentIndex]->right;
  /*The parent is replaced in the array by its left daughter */
  parentArray[parentIndex] = parentArray[parentIndex]->left;

  ++(*daughterIndex);
  if (cellDecrease) {
    --(*cellDecrease);
  }
  *totMem = *totMem + 2 * sizeof(struct Node); /* Memory allocation occurs in the newNode function*/

  return;
}

/* Clears memory for the tree. Stem must be cleared separately */
void clearTreeMemory(struct Node* node_ptr, unsigned long int* mem) {
  if (node_ptr->left != NULL) {
    clearTreeMemory(node_ptr->left, mem);
    clearTreeMemory(node_ptr->right, mem);
  }

  /* Clear node */
  free(node_ptr);
  *mem = *mem - sizeof(struct Node);
  return;
}

// Set all elements in an array of node pointers of size "size" to be NULL
void setArrayNull(struct Node** array, long int size) {
  for (int i = 0; i < size; i ++) {
    array[i] = NULL;
  }
  return;
}

/* Sets alls arrays used to store virus pointers to NULL for all elements */
void setVirusArrays(struct Node** activeArray, struct Node** virusArray, struct Node** latentIncompArray, struct Node** latentCompArray,
struct Node** virusSampleArray, struct Node** latentSampleArray, long int maxActive, long int maxVirus, long int maxIncomp,
long int maxComp, int totSampleVirus, int totSampleLatent) {

  setArrayNull(activeArray, maxActive);
  setArrayNull(virusArray, maxVirus);
  setArrayNull(latentIncompArray, maxIncomp);
  setArrayNull(latentCompArray, maxComp );
  setArrayNull(virusSampleArray, totSampleVirus);
  setArrayNull(latentSampleArray, totSampleLatent);

}

/* Clears all memory in program */
void ClearMemory(struct Node* root_ptr, struct Node* stem, struct Node* active[], struct Node* latent[],
struct Node* latentDefect[], struct Node* virus[], struct Node* sampleVirus[], struct Node* sampleLatent[],
unsigned long int* totMem) {

  /* Clear arrays */
  free(active);
  free(latent);
  free(latentDefect);
  free(virus);
  free(sampleVirus);
  free(sampleLatent);

  /* Clear tree */
  free(stem);
  clearTreeMemory(root_ptr, totMem);

  return;
}

/* Prunes a tip and the parent of that tip off of the tree */
void pruneTip(struct Node* tip, unsigned long int* totMem) {
  /* Finds the pointer to the sister*/
  struct Node* sister;
  if (tip->parent->right == tip) {
    sister = tip->parent->left;

  } else {
    sister = tip->parent->right;

  }

  /* Set sister to have parent birth time. Add latent time together*/
  sister->birth_time = tip->parent->birth_time;
  sister->totTimeLatent = sister->totTimeLatent + tip->parent->totTimeLatent;
  /* Set sister's parents to be the grandparents */
  sister->parent = tip->parent->parent;

  /* Set grandparents to point to sister */
  if (tip->parent->parent->left == tip->parent) {
    tip->parent->parent->left = sister;

  } else {
    tip->parent->parent->right = sister;

  }
  *totMem = *totMem - 2* sizeof(struct Node);
  free(tip->parent);
  free(tip);
  return;
}

/*Prints only the sampled branches in Newick format */
char* NewickAll(struct Node* node, char* treeString, long int* nodeName) {
  if (node->left == NULL && node->right == NULL) {
    Sasprintf(treeString, "%s%ld", treeString, *nodeName);
    *nodeName = *nodeName +  1;
    return (treeString);

  } else {

      Sasprintf(treeString, "%s(", treeString);
      treeString = NewickAll(node->left, treeString, nodeName);
      Sasprintf(treeString, "%s,", treeString);
      treeString = NewickAll(node->right, treeString, nodeName);
      Sasprintf(treeString, "%s)%ld", treeString, *nodeName);
      *nodeName = *nodeName +  1;
  }
  return(treeString);
}

// /* Prints a Newick format tree to the screen*/
// void printTree(struct Node* node) {
//
//   if (node->left == NULL && node->right == NULL) {
//     printf("%d", node->name);
//     return;
//
//   } else {
//     printf("(");
//     printTree(node->left);
//     printf(",");
//     printTree(node->right);
//     printf(")");
//     printf("%d", node->name);
//
//   }
// return;
// }

/* Samples a virus by moving it from an arrayTosample to the sampleArray
pickVirus is the index of the virus being sampled. sampleIndex is the index of the virus in the sampled array */
void sample(struct Node* arrayToSample[], struct Node* sampleArray[], int pickVirus, int* sampleIndex, double totTime, long int* numVirus) {
  long int lastIndex = *numVirus;
  /* Changes the variables related to sampling */
  arrayToSample[pickVirus]->sample_time = totTime;
  arrayToSample[pickVirus]->isSample = 1;

  /* If the sampled virus is latent, find the total time it was latent */
  if (arrayToSample[pickVirus]->isLatent == 1) {
    arrayToSample[pickVirus]->totTimeLatent = arrayToSample[pickVirus]->totTimeLatent + (totTime - arrayToSample[pickVirus]->timeToLatent);
  }

  /* Moves the virus to the sampled array */
  sampleArray[*sampleIndex] = arrayToSample[pickVirus];
  /* Take the virus out of the array to sample */
  removeVirus(arrayToSample, pickVirus, lastIndex - 1, numVirus);
  *sampleIndex = *sampleIndex + 1;
}

/* These set rates that change when the number ov viruses or the number of uninfected cells change */
void setInfectionRates (double* rate,
double* prodInfectionRates, long int numCellUninfect, long int numVirus) {
  long int virusXuninfected = numVirus * numCellUninfect;
  rate[1] = prodInfectionRates[0] * virusXuninfected;
  rate[6] = prodInfectionRates[1] * virusXuninfected;
  rate[7] = prodInfectionRates[2] * virusXuninfected;
}

/* Calculates all rates in the model */
void setRates(
  double* prodInfectionRates, double* parameters, double reactLatent, double latIncompDeath, double latCompDeath,
  double* rate, long int numVirus, long int numCellInfect, long int numCellUninfect, long int numLatentComp, long int numLatentIncomp) {

    rate[0] = parameters[0];
    rate[1] = prodInfectionRates[0] * numVirus * numCellUninfect;
    rate[2] = parameters[2] * numCellUninfect;
    rate[3] = parameters[3] * numCellInfect;
    rate[4] = parameters[4] * numCellInfect;
    rate[5] = parameters[5] * numVirus;

    rate[6] = prodInfectionRates[1] * numCellUninfect * numVirus;
    rate[7] = prodInfectionRates[2] * numCellUninfect * numVirus;
    rate[8] =  latIncompDeath * numLatentIncomp;
    rate[9] = latCompDeath * numLatentComp;
    rate[10] = reactLatent * numLatentComp;
  }


/* When a sampling time is reached, this function moves randomly choosen free viruses and latent viruses
to their respective sampled virus arrays. */
void sampleEvent(gsl_rng *r, int* sampleCounter, int numSampleVirus[], int numSampleLatent[],
  struct Node* virusArray[], struct Node* virusSampleArray[], struct Node* latentSampleArray[],
  struct Node* latentIncompArray[], struct Node* latentCompArray[],
  int* numVirusSample, int* numLatentSample, long int* numVirus, long int* numLatentIncomp, long int* numLatentComp, double totTime) {

  /* Used to store the randomly choosen viruses that are sampled */
  long int pickVirus;

  /* Sample Active */
  for (int count = 0; count < numSampleVirus[*sampleCounter]; count++) {
    pickVirus = gsl_rng_uniform_int (r, *numVirus);
    sample(virusArray, virusSampleArray, pickVirus, numVirusSample, totTime, numVirus);
  }

  /* Sample Latent */
  /* Confirms there are enough latent viruses to sample */
  if (! (numSampleLatent[*sampleCounter] < *numLatentIncomp + *numLatentComp || numSampleLatent[*sampleCounter] == 0)) {
    printf("Time latent incomp latent comp %f %ld %ld\n", totTime, *numLatentIncomp , *numLatentComp );
  }
  assert(numSampleLatent[*sampleCounter] < *numLatentIncomp + *numLatentComp || numSampleLatent[*sampleCounter] == 0);

  for (int count = 0; count < numSampleLatent[*sampleCounter]; count++) {
    pickVirus = gsl_rng_uniform_int (r, *numLatentIncomp + *numLatentComp);

    /* Sample latent incompetent virus */
    if (pickVirus < *numLatentIncomp) {
      sample(latentIncompArray, latentSampleArray, pickVirus, numLatentSample, totTime, numLatentIncomp);

    /* Sample latent competent virus */
    } else {
      pickVirus = pickVirus - *numLatentIncomp;
      sample(latentCompArray, latentSampleArray, pickVirus, numLatentSample, totTime, numLatentComp);

    }
  }
  /* Keeps track of which sampling event is next */
  ++(*sampleCounter);

  return;
}

/* Resets variables and data structures to their original values */
struct Node* reset(struct Node* stem, struct Node* root_ptr, int* sampleCounter, unsigned long int* numEvents, unsigned long int* totMem,
struct Node** activeArray, struct Node** virusArray, struct Node** latentIncompArray, struct Node** latentCompArray,
struct Node** virusSampleArray, struct Node** latentSampleArray, long int maxActive, long int maxVirus, long int maxIncomp,
long int maxComp, int totSampleVirus, int totSampleLatent) {
  clearTreeMemory(stem->left, totMem);
  setVirusArrays(activeArray, virusArray, latentIncompArray, latentCompArray, virusSampleArray, latentSampleArray,
    maxActive, maxVirus, maxIncomp, maxComp, totSampleVirus, totSampleLatent);

  root_ptr = newNode(0, 0, stem); /* Root of the tree */
  stem->left = root_ptr;
  activeArray[0] = root_ptr;
  *sampleCounter = 0;
  *numEvents = 0 ;
  *totMem = *totMem +  sizeof(struct Node); //for root

  return root_ptr;
}

/* Writes a file filename with contents fileContents */
void writeTxt(char* fileName, char* fileContents) {
  FILE *file;
  file = fopen(fileName, "w");
  if (! file) {
	  fprintf(stderr,"%s failed to open. Exiting the program.\n", fileName);
	  exit(1);
  }
  fputs(fileContents, file);
  fclose(file);
  free(fileContents);
  return;
}

/*Prints only the sampled branches in Newick format */
char* NewickSampled(struct Node* node, char* treeString, double previousTime, long int* nodeName) {
  long int tmp_nodeName;
  if (node->left == NULL && node->right == NULL) {
    Sasprintf(treeString, "%s%ld:%f", treeString, *nodeName, node->sample_time - previousTime);
    *nodeName = *nodeName +  1;
    return (treeString);

  } else {
    /* Both daughter nodes are sampled, new node */
    if (node->left->isSample == 1 && node->right->isSample == 1){
      tmp_nodeName = *nodeName;
      *nodeName = *nodeName +  1;
      Sasprintf(treeString, "%s(", treeString);
      treeString = NewickSampled(node->left, treeString, node->left->birth_time, nodeName);
      Sasprintf(treeString, "%s,", treeString);
      treeString = NewickSampled(node->right, treeString, node->right->birth_time, nodeName);
      Sasprintf(treeString, "%s)%ld:%f", treeString, tmp_nodeName, node->left->birth_time - previousTime);


    } else if (node->left->isSample == 1) {
      treeString = NewickSampled(node->left, treeString, previousTime, nodeName);

    } else if (node->right->isSample == 1) {
      treeString = NewickSampled(node->right, treeString, previousTime, nodeName);
    } else {
      fprintf(stderr, "There are no sampled viruses\n");
      exit(1);
    }

  }
  return(treeString);
}

/* Used to determine which nodes are part of the sampled tree. Traces a sampled daughter lineage backward
until it reaches the stem or a lineage that has already been marked as sampled */
void findSampledParent(struct Node* node, struct Node* stem) {

  // If the lineage has already been marked as sampled, no need to check again
  if (node->isSample == 1) {
    return;

  // Mark the lineage as sampled
  } else {
    node->isSample = 1;
  }

  // If the lineage isn't the stem, follow the lineage backward
  if (node == stem){
    return;
  } else {
    findSampledParent(node->parent, stem);
  }
  return;
}

/* Marks all of the nodes as sampled that are ancestral to any sampled lineages */
void findSampledAll(struct Node** sampleActiveArray, struct Node** sampleLatentArray, int sizeActive, int sizeLatent, struct Node* stem) {
  /* Note: tips are already marked as sampled, so start with the parents of the sampled tips */
  /* For every free virus sampled, mark all of its ancestors as sampled  */
  for(int i = 0; i < sizeActive; i++) {
    findSampledParent(sampleActiveArray[i]->parent, stem);
  }

  /* For every latent virus sampled, mark all of its ancestors as sampled  */
  for(int i = 0; i < sizeLatent; i++) {
    findSampledParent(sampleLatentArray[i]->parent, stem);

  }
  return;
}

/* Counts the number of lines in the sample times file */
int findFileLength(char* filename) {
  FILE *fp;
  fp = fopen(filename, "r");
  if (!fp) {
    fprintf(stderr, "Cannot open sample times file.\nExiting the program.\n");
    exit(1);
  }
  int count = 0;
  char buff[200];

  int checkScan = fscanf(fp, "%s", buff);

  while(checkScan != EOF) {
    checkScan = fscanf(fp, "%s", buff);
    count++;
  }
  fclose(fp);
  return count;
}

/* Reads in a csv formal file with the samples times, number of active viruses to sample, number of latent viruses to sample */
void readSampleTimes (double sampleTimes[], int numSampleVirus[], int numSampleLatent[], char* filename, int fileLength) {

  regex_t regDecimal, regInt;
  int regCompiled, regMatchDec, regMatchInt;
  //regCompiled = regcomp(&regDecimal, "^([0-9]+)((.)([0-9]+))$" , REG_EXTENDED);
  /* Regular expression for decimal number, integer, integer */
  regCompiled = regcomp(&regDecimal, "^([0-9]*)((\\.)([0-9]+))(,([0-9]+))(,([0-9]+))$" , REG_EXTENDED);
  if (regCompiled == 1) {
    fprintf(stderr, "Regular expression did not compile.\n");
    exit(1);
  }
  /* Regular expression for integer, integer, integer */
  regCompiled = regcomp(&regInt, "^([0-9]+),([0-9]+),([0-9]+)$" , REG_EXTENDED);
  if (regCompiled == 1) {
    fprintf(stderr, "Regular expression did not compile.\n");
    exit(1);
  }

  FILE *sampleTimesFile;
  sampleTimesFile = fopen(filename, "r");
  if (! sampleTimesFile) {
	  fprintf(stderr, "%s failed to open. Exiting the program.\n", filename);
	  exit(1);
  }

  char buff[200];
  char* ptr;
  int checkScan;
  int i = 0;
  char *token;

  checkScan = fscanf(sampleTimesFile, "%s", buff);
  regMatchDec = regexec(&regDecimal, buff, 0, NULL, 0);
  regMatchInt = regexec(&regInt, buff, 0, NULL, 0);

  if (regMatchDec != 0 && regMatchInt != 0 ) {
    fprintf(stderr, "Sample times file is in the incorrect format. \nExiting the program. \n");
    exit(1);
  }

  while(checkScan != EOF) {
    /* Storing the samples times and number of viruses to sample in the appropriate arrays*/
    sampleTimes[i] = strtod(buff, &ptr);
    token = strtok(ptr + 1, ",");
    numSampleVirus[i] = atoi(token);
    token = strtok(ptr + strlen(token) + 2, "\n");
    numSampleLatent[i] = atoi(token);

    /* Reads in next line and checks the format */
    checkScan = fscanf(sampleTimesFile, "%s", buff);
    regMatchDec = regexec(&regDecimal, buff, 0, NULL, 0);
    regMatchInt = regexec(&regInt, buff, 0, NULL, 0);

    if (checkScan != EOF && regMatchDec != 0 && regMatchInt != 0 ) {
      fprintf(stderr, "Sample times file is in the incorrect format. \nExiting the program. \n");
      exit(1);
    }

    i++;
  }

  if (i != fileLength) {
    fprintf(stderr, "Problem reading in sample time file. \n");
    exit(1);
  }
  regfree(&regDecimal);
  regfree(&regInt);
  fclose(sampleTimesFile);
  return;
}

/* Checks sampling input meets basic requirements - times are sorted, the number of viruses to sample is non negative*/
void checkValidSampleTimes(double sampleTimes[], int numSampleVirus[], int numSampleLatent[], int fileLength) {

  if(numSampleVirus[0] < 0 || numSampleLatent[0] < 0) {
    fprintf(stderr, "The number of sampled viruses must be non-negative. \n");
    exit(1);
  }

  for(int i = 1; i < fileLength; i++) {
    if (sampleTimes[i-1] >=sampleTimes[i]) {
      fprintf(stderr, "Sample times must be strictly increasing. \n");
      exit(1);
    }

    if(numSampleVirus[i] < 0 || numSampleLatent[i] < 0) {
      fprintf(stderr, "The number of sampled viruses must be non-negative. \n");
      exit(1);
    }
  }
  return;
}

/* Resets variables to their initial conditions */
void resetCounts(long int* numVirus, long int numVirusInit, long int* numCellUninfect, long int numCellUninfectInit,
    long int* numCellInfect, long int numCellInfectInit, long int* numLatentComp, long int numLatentCompInit,
    long int* numLatentIncomp, long int numLatentIncompInit, int* numLatentSample, int* numVirusSample, int* sampleCounter, double* totTime) {
  *numVirus        = numVirusInit;
  *numCellUninfect = numCellUninfectInit;
  *numCellInfect   = numCellInfectInit;
  *numLatentComp   = numLatentCompInit;
  *numLatentIncomp = numLatentIncompInit;
  *numLatentSample = 0;
  *numVirusSample  = 0;
  *sampleCounter   = 0;

  *totTime = 0;

  return;
}

/*Prints only the sampled branches in Newick format */
char* LatentSampled(struct Node* node, char* latentString, double timeLatent, long int* nodeName) {
  if (node->left == NULL && node->right == NULL) {
    Sasprintf(latentString, "%s%ld:%f,%f,%d\n", latentString, *nodeName,
    timeLatent + node->totTimeLatent, node->sample_time, node->isLatent);

    *nodeName = *nodeName + 1;
    return (latentString);

  } else {
    /* Both daughter nodes are sampled, new node */
    if (node->left->isSample == 1 && node->right->isSample == 1) {
      Sasprintf(latentString, "%s%ld:%f,%f,%d\n", latentString, *nodeName, timeLatent + node->totTimeLatent, node->sample_time, node->isLatent); /* -1 indicates the internal node is not sampled */
      *nodeName = *nodeName + 1;

      latentString = LatentSampled(node->left, latentString, 0, nodeName);
      latentString = LatentSampled(node->right, latentString, 0, nodeName);

    } else if (node->left->isSample == 1) {
      latentString = LatentSampled(node->left, latentString, timeLatent + node->totTimeLatent, nodeName);

    } else if (node->right->isSample == 1) {
      latentString = LatentSampled(node->right, latentString, timeLatent + node->totTimeLatent, nodeName);
    } else {
      fprintf(stderr, "There are no sampled viruses\n");
      exit(1);
    }

  }
  return(latentString);
}

/* Writes the starting seed to an output file*/
void writeSeed(unsigned int RGSeed, char* outfile) {
  char seed[15];
  sprintf(seed, "%u", RGSeed);

  char seedFile[50] = "";
  strcat(seedFile, outfile);
  strcat(seedFile, "seed.txt");
  FILE *file;
  file = fopen(seedFile, "w");
  if (! file) {
	  fprintf("%s failed to open. Exiting the program.\n", seedFile);
	  exit(1);
  }
  fputs(seed, file);
  fclose(file);

  return;
}

/* Increases the memory allocated to an array of nodes by 20 %*/
struct Node** reallocArray(long int* sizeArrayOld, struct Node** array, unsigned long int* totMem) {
  *totMem = *totMem + ((long int)(*sizeArrayOld * 0.2)) * sizeof(struct Node *);
  *sizeArrayOld = *sizeArrayOld * 1.2;
  struct Node** nodeArray = realloc(array, *sizeArrayOld * sizeof(struct Node *));
  if(!nodeArray) {
    fprintf(stderr, "Problem reallocating array memory\n");
    exit(1);
  }
  return nodeArray;
}

/* Decreases the memory allocated to an array of nodes by 1/6 */
struct Node** reallocArrayDecrease(long int* sizeArrayOld, struct Node** array, unsigned long int* totMem) {
  long int amountDecrease = ((long int)(*sizeArrayOld / 6 ))* sizeof(struct Node *);
  *totMem = *totMem - amountDecrease;
  *sizeArrayOld = *sizeArrayOld -  ((long int)(*sizeArrayOld / 6 ));
  struct Node** nodeArray = realloc(array, *sizeArrayOld * sizeof(struct Node *));

  if(!nodeArray) {
    fprintf(stderr, "Problem reallocating array memory\n");
    exit(1);
  }
  return nodeArray;
}

/* Used to parse the control file. Assigns variables from control file*/
void parseControlFile(char* fileName, char* sampleTime, double* mLBlood,
double* parameters, double* probLatent, double* reactLatent, double* probDefect,
double* latIncompDeath, double* latCompDeath, double* ARTstart, unsigned int* seed,
int* seedChange, int* volChange, int* sampleFileChange) {

  FILE *file;
  file = fopen(fileName, "r");

  if (!file) {
    fprintf(stderr, "Problem opening control file. Exiting.\n");
    exit(1);
  }
  char buf[BUF_LEN];
  char removeWhite[BUF_LEN];
  int charIndex = 0;
  int noWhiteIndex = 0;
  char *ptr;

  regex_t regDecimal, regInt, regString;
  int regCompiled, regMatchDec, regMatchInt, regMatchString;
  /* Regex for word = decimal number */
  regCompiled = regcomp(&regDecimal, "^([a-zA-Z]+)(=)([0-9]*)((\\.)([0-9]+))$" , REG_EXTENDED);
  if (regCompiled == 1) {
    fprintf(stderr, "Regular expression did not compile.\n");
    exit(1);
  }

  /* Regex for word=integer*/
  regCompiled = regcomp(&regInt, "^([a-zA-Z]+)(=)([0-9]+)$" , REG_EXTENDED);
  if (regCompiled == 1) {
    fprintf(stderr, "Regular expression did not compile.\n");
    exit(1);
  }
  /*Regex for sample_times=word, word has letters, numbers, period, underscore, or dashes */
   regCompiled = regcomp(&regString, "^(sample_times)(=)([[:alnum:]._/-]+)$" , REG_EXTENDED);
  if (regCompiled == 1) {
    fprintf(stderr, "Regular expression did not compile.\n");
    exit(1);
  }

  /* Reads in next line*/
  while(fgets(buf, BUF_LEN, file)) {
    noWhiteIndex = 0;
    charIndex = 0;
    /* Strips whitespace from line*/
    while(buf[charIndex]!= '\0') {
      if(!isspace(buf[charIndex])) {
        removeWhite[noWhiteIndex] = buf[charIndex];
        noWhiteIndex++;
      }
      charIndex++;
    }
    removeWhite[noWhiteIndex] = '\0';

    /* Checks if line matches regular expressions */
    regMatchDec = regexec(&regDecimal, removeWhite, 0, NULL, 0);
    regMatchInt = regexec(&regInt, removeWhite, 0, NULL, 0);
    regMatchString = regexec(&regString, removeWhite, 0, NULL, 0);

    if(regMatchDec && regMatchInt && regMatchString) {
      fprintf(stderr, "Control file line is in the wrong format. \n %s \n Exiting.\n", buf);
      exit(1);
    }

    strtok(removeWhite, "=");

    char* input = removeWhite + strlen(removeWhite) + 1;

    /* Checks if the regular expression matched was for the string */
    if ( regMatchString == 0) {
	(*sampleFileChange)++;
      strcpy(sampleTime, input);
    } else if (strcmp(removeWhite, "volume") == 0) {
	(*volChange)++;
        *mLBlood = strtod(input, &ptr);

    } else if (strcmp(removeWhite, "lamdba") == 0) {
        parameters[0] = strtod(input, &ptr);

    } else if (strcmp(removeWhite, "kappa") == 0) {
        parameters[1] = strtod(input, &ptr);

    } else if (strcmp(removeWhite, "d") == 0) {
        parameters[2] = strtod(input, &ptr);

    } else if (strcmp(removeWhite, "delta") == 0) {
        parameters[3] = strtod(input, &ptr);

    } else if (strcmp(removeWhite, "pi") == 0) {
        parameters[4] = strtod(input, &ptr);

    } else if (strcmp(removeWhite, "c") == 0) {
        parameters[5] = strtod(input, &ptr);

    } else if (strcmp(removeWhite, "eta") == 0) {
        *probLatent = strtod(input, &ptr);

    } else if (strcmp(removeWhite, "alpha") == 0) {
        *reactLatent = strtod(input, &ptr);

    } else if (strcmp(removeWhite, "gamma") == 0) {
        *probDefect = strtod(input, &ptr);

    } else if (strcmp(removeWhite, "tau") == 0) {
        *latIncompDeath = strtod(input, &ptr);

    } else if (strcmp(removeWhite, "seed") == 0) {
        *seed = strtoul(input, &ptr, 10);
	(*seedChange)++;

    } else if (strcmp(removeWhite, "sigma") == 0) {
        *latCompDeath = strtod(input, &ptr);

    } else if (strcmp(removeWhite, "ARTstart") == 0) {
        *ARTstart = strtod(input, &ptr);

    } else {
      fprintf(stderr, "Control file keyword not recognized. Exiting.\n");
      exit(1);
    }
  }

  regfree(&regDecimal);
  regfree(&regInt);
  regfree(&regString);

}

unsigned long int findTreeStringMem(long int numVirus, long int numCellInfect, long int numLatentComp, long int numLatentIncomp, int numVirusSample, int numLatentSample, double timeLastSample) {
  unsigned long int numTips = numVirus + numCellInfect + numLatentComp + numLatentIncomp + numVirusSample + numLatentSample;
  unsigned long int internalNode = numTips - 1;
  unsigned long int numBranches = numTips * 2 - 1;
  unsigned long int numNodes = 2 * numTips - 1;

  unsigned long int maxLength = 0;
   /* Number of commas and paraentheses*/
   maxLength = maxLength + 3 * internalNode;

   /* Length of node names*/
   maxLength = maxLength + (log10(numNodes) + 1) * numNodes;

   /* length of branch lengths  (digit before decimal + 8 decimals + (decimal  + colon))*/
   maxLength = maxLength + ((log10(timeLastSample) + 1) + 8 + 2) * numBranches;

   /* Number of commas and paraentheses*/
   maxLength = maxLength + 3 * internalNode;

   return(maxLength);
}

unsigned long int findLatentStringMem(long int numVirus, long int numCellInfect, long int numLatentComp, long int numLatentIncomp, int numVirusSample, int numLatentSample, double timeLastSample) {
  unsigned long int numTips = numVirus + numCellInfect + numLatentComp + numLatentIncomp + numVirusSample + numLatentSample;
  unsigned long int numNodes = 2 * numTips - 1;

  unsigned long int maxLength = 0;
   /* Each double- 8 decimal places + 1 decimal (9), log10 + 1 digits before the decimal
     1 -is latent, int between 0 and 8
     10 is long int for the array index
     6 is the 4 commas, colon, and newline
     */
   maxLength = (((log10(timeLastSample) + 1) + 9) * 3 + 1 + 10 + 6)* numNodes;


   return(maxLength);
}

/// New checkpoint
/* Saves the entire tree structure as a newick format tree. The birth times for
each node are stored as the branch lengths are typically in a newick format. */
char* checkPointTree(struct Node* node, char* treeString, long int* nodeName){
  long int tmp_nodeName;
  if (node->left == NULL && node->right == NULL) {
    sprintf(treeString, "%ld:%.8f", *nodeName, node->birth_time);
    *nodeName = *nodeName +  1;
    return (treeString + strlen(treeString));

  } else {
      tmp_nodeName = *nodeName;
      *nodeName = *nodeName +  1;
      sprintf(treeString, "(");
      treeString = treeString + 1;
      treeString = checkPointTree(node->left, treeString, nodeName);
      sprintf(treeString, ",");
      treeString = treeString + 1;
      treeString = checkPointTree(node->right, treeString, nodeName);
      sprintf(treeString, ")%ld:%.8f", tmp_nodeName, node->birth_time);
      treeString = treeString + strlen(treeString);
  }
  return(treeString);
}

/* Saves the latent history for the entire tree structure. Each line corresponds
to one node. The nodes are in the order of a preorder tree tranversal */
char* checkPointLatent(struct Node* node, char* latentString, long int* nodeName) {
  if (node->left == NULL && node->right == NULL) {
    sprintf(latentString, "%ld:%.8f,%.8f,%.8f,%d,%ld\n", *nodeName, node->timeToLatent,
    node->totTimeLatent, node->sample_time, node->isLatent, node->arrayNum);

    *nodeName = *nodeName + 1;
    return (latentString + strlen(latentString));

  } else {
      sprintf(latentString, "%ld:%.8f,%.8f,%.8f,%d,%ld\n", *nodeName, node->timeToLatent, node->totTimeLatent, node->sample_time, node->isLatent, node->arrayNum); /* -1 indicates the internal node is not sampled */
      *nodeName = *nodeName + 1;

      latentString = latentString + strlen(latentString);
      latentString = checkPointLatent(node->left, latentString,  nodeName);
      latentString = checkPointLatent(node->right, latentString, nodeName);

  }
  return(latentString  + strlen(latentString));
}

/* Write a string at a binary file */
int writeBinary(char* treeString, char * fileName) {
  long long int size = strlen(treeString) + 1;

  int status = EXIT_SUCCESS;
  FILE *fp;
  if ((fp = fopen(fileName, "wb")) == NULL) {
      fputs("Cannot open binary file\n", stderr);
      return EXIT_FAILURE;
    }

  if (fwrite(treeString, size, 1, fp) != 1) {
      fputs("Cannot write binary file\n", stderr);
      status = EXIT_FAILURE;
    }

    if (fclose(fp) == EOF) {
      fputs("Failed to close file\n", stderr);
      status = EXIT_FAILURE;
    }
    return status;
}

/* Writes a struct to a binary file */
int writeCountsBinary(struct cellCountTimes* counts, char * fileName) {
  int size = sizeof(struct cellCountTimes);
  int status = EXIT_SUCCESS;
  FILE *fp;
  if ((fp = fopen(fileName, "wb")) == NULL) {
      fputs("Cannot open binary file\n", stderr);
      return EXIT_FAILURE;
    }

  if (fwrite(counts, size, 1, fp) != 1) {
      fputs("Cannot write binary file\n", stderr);
      status = EXIT_FAILURE;
    }

    if (fclose(fp) == EOF) {
      fputs("Failed to close file\n", stderr);
      status = EXIT_FAILURE;
    }
    return status;
}

/* Saves the sampling times, number of active and latent samples as a binary file
sampleFileLength is the number of sampling events*/
int writeSampleFileBinary(int sampleFileLength, double* sampleTimes, int* numSampleVirus, int* numSampleLatent,
char * fileName) {
  int status = EXIT_SUCCESS;
  struct sampleTimesLine line;
  int size = sizeof(struct sampleTimesLine);
  FILE *fp;

  if ((fp = fopen(fileName, "wb")) == NULL) {
      fputs("Cannot open binary file for sampling times\n", stderr);
      return EXIT_FAILURE;
    }

  for(int i = 0; i < sampleFileLength; i++) {
    line.time = sampleTimes[i];
    line.numActive = numSampleVirus[i];
    line.numLatent = numSampleLatent[i];

    if (fwrite(&line, size, 1, fp) != 1) {
        fputs("Cannot write binary file for sample times\n", stderr);
        status = EXIT_FAILURE;
      }
  }

  if (fclose(fp) == EOF) {
    fputs("Failed to close file\n", stderr);
    status = EXIT_FAILURE;
  }
  return status;

}

/* Reads in the binary file with the sample times and number of samples to be taken  */
int readSampleFileBinary(int sampleFileLength, double* sampleTimes, int* numSampleVirus, int* numSampleLatent,
char * outfile) {

  char fileName[50] = "\0";
  strcat(fileName, outfile);
  strcat(fileName, "sampleTimes.chkpt");

  int status = EXIT_SUCCESS;
  FILE *fp;
  int size = sizeof(struct sampleTimesLine);
  struct sampleTimesLine line;

  if ((fp = fopen(fileName, "rb")) == NULL) {
       fputs("Cannot open counts file\n", stderr);
       return EXIT_FAILURE;
    }

  for (int i = 0; i < sampleFileLength; i++) {
    if (fread(&line, size, 1, fp) != 1) {
         fputs("Cannot read from counts file\n", stderr);
         status = EXIT_FAILURE;
      }
    sampleTimes[i] = line.time;
    numSampleVirus[i] = line.numActive;
    numSampleLatent[i] = line.numLatent;
  }

    fclose(fp);
    return status;
}


/* Reads in binary files with the state of the random number generator, strings
of the tree in newick format, the latent history, and a struct with parameters, counts,
and other relevant states of the simulation */
int readInBinaries(char* countsFile, char* gslFile, char* latentFile, char* treeFile,
  struct cellCountTimes * counts, char ** treeFilePointer, char ** latentFilePointer, gsl_rng * r) {

  FILE *f_rng;
  if ((f_rng = fopen(gslFile, "rb")) == NULL) {
       fputs("Cannot open gsl file\n", stderr);
       return EXIT_FAILURE;
    }
  int exit = gsl_rng_fread(f_rng, r);

  if (exit == GSL_EFAILED) {
    fputs("Cannot read random generator binary file\n", stderr);
    return EXIT_FAILURE;
  }
  fclose(f_rng);

  int status = EXIT_SUCCESS;
    FILE *fp;
    int size = sizeof(struct cellCountTimes);

  if ((fp = fopen(countsFile, "rb")) == NULL) {
       fputs("Cannot open counts file\n", stderr);
       return EXIT_FAILURE;
    }
  if (fread(counts, size, 1, fp) != 1) {
       fputs("Cannot read from counts file\n", stderr);
       status = EXIT_FAILURE;
    }
    fclose(fp);

    /* Allocates memory to read in the tree file and latent history file */
    char* treeString = malloc(counts->sizeTree * sizeof(char));
    char* latentString = malloc(counts->sizeLatent * sizeof(char));
    *treeFilePointer = treeString;
    *latentFilePointer = latentString;

    size_t sizeTree = counts->sizeTree * sizeof(char);
    size_t sizeLatent = counts->sizeLatent * sizeof(char);

    /* Reads in string with tree in newick format*/
    FILE * fTree;
    if ((fTree = fopen(treeFile, "rb")) == NULL) {
         fputs("Cannot open tree binary file\n", stderr);
         return EXIT_FAILURE;
      }
    if (fread(treeString, sizeTree, 1, fTree) != 1) {
         fputs("Cannot read from tree binary  file\n", stderr);
         status = EXIT_FAILURE;
      }
    fclose(fTree);

    /* Reads in string with latent history information */
    FILE * fLatent;
    if ((fLatent = fopen(latentFile, "rb")) == NULL) {
         fputs("Cannot open latent binary file\n", stderr);
         return EXIT_FAILURE;
      }

    if (fread(latentString, sizeLatent, 1, fLatent) != 1) {
         fputs("Cannot read from latent binary  file\n", stderr);
         status = EXIT_FAILURE;
      }
    fclose(fLatent);
    return status;
}


// Read in files
/* Create a new new node. Allocates memory for the node and the latency states
structure. Assignment values to the variables in the node structure */
struct Node* newNodeSampledTree(struct Node* parent, unsigned long int* totMem) {
  /* Allocates memory for a new node */
  struct Node* node_ptr;
  node_ptr  = malloc(sizeof(struct Node));
  *totMem = *totMem + sizeof(struct Node);

  if (node_ptr == NULL) {
    fprintf(stderr, "Insufficient Memory\n");
    exit(1);
  }

  node_ptr->birth_time = -1; /* Sets to impossible value. Must get reset */
  node_ptr->sample_time = -1;
  node_ptr->isSample = 3;

  /* Sets pointers to daughter nodes */
  node_ptr->left = NULL;
  node_ptr->right = NULL;
  node_ptr->parent = parent;

  node_ptr->isLatent = 8; /* Sets to value that is not used. Must get reset */
  node_ptr->timeToLatent = -1;
  node_ptr->totTimeLatent = -1;
  node_ptr->arrayNum = 0;
  return(node_ptr);
}

/* The function parses a string into a tree structure. Every time there's an open
parenthesis, two new nodes are created and the function is called twice. The
function keeps track of which character in the treeString (currChar) it is looking
at. */
int makeTree(struct Node* node_ptr, char* treeString, int currChar, unsigned long int *totMem) {

  /* Note: You should never actually reach the end of string character. The last
  time you return is after the last closed parenthesis */
  while( treeString[currChar] != '\0') {

    /* If the current character is an open parenthesis, make new node, call
    makeTree recursively. currChar is updated along the way */
    if (treeString[currChar] == '(') {
      node_ptr->left = newNodeSampledTree(node_ptr, totMem);
      currChar = makeTree(node_ptr->left, treeString, currChar + 1, totMem);
      node_ptr->right = newNodeSampledTree(node_ptr, totMem);
      currChar = makeTree(node_ptr->right, treeString, currChar, totMem);


    } else {
      /* First, find the node name. Then, find the branch length. Return the
      current character */
      char* ptr; /* Needed for the strtod function */
      //char* residualTreeString = malloc(strlen(treeString) + 1); /* Used so the original string isn't modified */
      char* residualTreeString = treeString + currChar;
      //strcpy(residualTreeString, treeString + currChar);

      regex_t regDecimal, regInt;
      int regCompiled, regMatch;
      /* Regex of decimal number */
      regCompiled = regcomp(&regDecimal, "^([0-9]+)((\\.)([0-9]+))$" , REG_EXTENDED);
      if (regCompiled == 1) {
        fprintf(stderr, "Regular expression did not compile.\n");
        exit(1);
      }
      /* Regex of integer number */
      regCompiled = regcomp(&regInt, "^([0-9]+)$" , REG_EXTENDED);
      if (regCompiled == 1) {
        fprintf(stderr, "Regular expression did not compile.\n");
        exit(1);
      }

      /* Finds the nodeName by looking for the next ":". Convert it into an
      integer, save it in the node structure. Update currChar. */
      char* nodeName = strtok(residualTreeString, ":"); /*Note this function modified residualTreeString */

      regMatch = regexec(&regInt, nodeName, 0, NULL, 0);
      if (regMatch != 0) {
        fprintf(stderr, "Problem reading in tree file. Regular expression does not match a integers.\n");
        exit(1);
      }

      currChar = currChar + strlen(nodeName) + 1;

      //residualTreeString = strcpy(residualTreeString, treeString + currChar);
      residualTreeString = treeString + currChar;

      /* Finds the "branch length", converts it to a double, saves it in the node
      structure. What is usually the branch length in newick format is actually
      the birth time  */
      char* branchLength = strtok(residualTreeString, ",)");

      regMatch = regexec(&regDecimal, branchLength, 0, NULL, 0);
      if (regMatch != 0) {
        fprintf(stderr, "Problem reading in tree file. Regular expression does not match a decimal number.\n");
        exit(1);
      }
      node_ptr->birth_time = strtod(branchLength, &ptr);
      currChar = currChar + strlen(branchLength) + 1 ;

      //free(residualTreeString);
      regfree(&regDecimal);
      regfree(&regInt);

      /* Returns the updated current character */
      return(currChar);
    }
  }
  return(currChar);
}

/* Reads in a file with the amount of time spent in a latent state */
/* Each row in the file is in the format:
nodeName (number), node->timeToLatent, node->totTimeLatent, node->sample_time, node->isLatent, node->arrayNum
isLatent is coresponds to the cell classes.
Places tips in the appropriate arrays */
void readLatent (struct Node* node, char* latentFile, struct Node** activeArray, struct Node** virusArray, struct Node** latentIncompArray, struct Node** latentCompArray,
struct Node** virusSampleArray, struct Node** latentSampleArray, char** ptr) {

  double latentTime, sampleTime, timeToLatent;
  char* endPtr; /* Used in strtoul function */
  char* token = strtok(*ptr, ":");

  timeToLatent = strtod(*ptr + strlen(token) + 1, ptr);
  node->timeToLatent = timeToLatent;

  latentTime = strtod(*ptr+1, ptr);
  node->totTimeLatent = latentTime;

  sampleTime = strtod(*ptr+1, ptr);
  node->sample_time = sampleTime;

  node->isLatent = atoi(*ptr + 1);
  if (node->isLatent < 0 && node->isLatent  > 7) {
    fprintf(stderr, "Problem reading in latent states. The latent state is not between 2 and 7. \n");
    exit(1);
  }

  node->arrayNum = strtoul(*ptr + 3, &endPtr, 10);

  if (node->isLatent > 1 && node->isLatent  < 8 ) {
    if (node->isLatent == 2) {
      activeArray[node->arrayNum] = node;
      node->isLatent = 0;
      node->isSample = 2;

    } else if (node->isLatent == 3) {
      virusArray[node->arrayNum] = node;
      node->isLatent = 0;
      node->isSample = 2;

    } else if (node->isLatent == 4) {
      latentIncompArray[node->arrayNum] = node;
      node->isLatent = 1;
      node->isSample = 2;

    } else if (node->isLatent == 5) {
      latentCompArray[node->arrayNum] = node;
      node->isLatent = 1;
      node->isSample = 2;

    } else if (node->isLatent == 6) {
      virusSampleArray[node->arrayNum] = node;
      node->isLatent = 0;
      node->isSample = 1;

    } else if (node->isLatent == 7) {
      latentSampleArray[node->arrayNum] = node;
      node->isLatent = 1;
      node->isSample = 1;

    } else{
      printf("Problem with loading arrays for checkpoint.\n");
      assert(0==1);
    }

  } else {
    node->isSample = 2;
  }

  *ptr = endPtr;

  if (node->left == NULL && node->right == NULL) {
    return ;

  } else {
    readLatent(node->left, latentFile, activeArray, virusArray, latentIncompArray, latentCompArray,
    virusSampleArray, latentSampleArray, ptr);

    readLatent(node->right, latentFile, activeArray, virusArray, latentIncompArray, latentCompArray,
    virusSampleArray, latentSampleArray, ptr);

  }
  return;
}

/* For all of the tips in the tree, change the latent state in the node struct
based off of what cell class the node is. */
void findCellTypes(struct Node** activeArray, struct Node** virusArray, struct Node** latentIncompArray, struct Node** latentCompArray,
  struct Node** virusSampleArray, struct Node** latentSampleArray,
  long int numInfect, long int numVirus, long int numLatentIncomp, long int numLatentComp, int numVirusSample, int numLatentSample) {
    for (unsigned long int i = 0; i < numInfect; i++) {
      activeArray[i]->isLatent = 2;
      activeArray[i]->arrayNum = i;
    }

    for (unsigned long int i = 0; i < numVirus; i++) {
      virusArray[i]->isLatent = 3;
      virusArray[i]->arrayNum = i;
    }

    for (unsigned long int i = 0; i < numLatentIncomp; i++) {
      latentIncompArray[i]->isLatent = 4;
      latentIncompArray[i]->arrayNum = i;
    }

    for (unsigned long int i = 0; i < numLatentComp; i++) {
      latentCompArray[i]->isLatent = 5;
      latentCompArray[i]->arrayNum = i;
    }

    for (unsigned long int i = 0; i < numVirusSample; i++) {
      virusSampleArray[i]->isLatent = 6;
      virusSampleArray[i]->arrayNum = i;
    }

    for (unsigned long int i = 0; i < numLatentSample; i++) {
      latentSampleArray[i]->isLatent = 7;
      latentSampleArray[i]->arrayNum = i;
    }
}


/* Saves the state of the random number generator as a binary file */
void writeGslState(gsl_rng *r, char* outfile) {

  char gslBiFile[100] = "\0";
  strcat(gslBiFile, outfile);
  strcat(gslBiFile, "gsl.chkpt");

  FILE *gslFile;
  gslFile = fopen(gslBiFile, "wb");
  if (gslFile == NULL) {
      fputs("Cannot open file to write random number generator \n", stderr);
    }
  if(gsl_rng_fwrite(gslFile, r) == GSL_EFAILED) {
    fputs("Cannot write random number generator to binary file\n", stderr);
  }
  if (fclose(gslFile) == EOF) {
    fputs("Failed to close random number generator file\n", stderr);
  }
}

/* Places information that will be needed to restart the simulation into a structure.
Saves the structure as a binary file */
void writeCounts(long int numVirus, long int numCellInfect, long int numCellUninfect, long int numLatentComp,
long int numLatentIncomp, int numVirusSample, int numLatentSample, double totTime, double waitTime, char* treeString, char* latentString,
double mLBlood, double* parameters, double probLatent, double reactLatent, double probDefect, double latIncompDeath,
double latCompDeath, int sampleFileLength, int sampleCounter, long int maxActive, long int maxVirus, long int maxIncomp,
long int maxComp, int totSampleVirus, int totSampleLatent, int timeLastSample, unsigned long int numEvents, double ARTstart, char* outfile) {
  struct cellCountTimes counts;
  counts.numVirus = numVirus;
  counts.numCellInfect = numCellInfect;
  counts.numCellUninfect = numCellUninfect;
  counts.numLatentComp = numLatentComp;
  counts.numLatentIncomp = numLatentIncomp;

  counts.numVirusSample = numVirusSample;
  counts.numLatentSample = numLatentSample;
  counts.totTime = totTime;
  counts.waitTime = waitTime;
  counts.sizeTree = strlen(treeString) + 1;
  counts.sizeLatent = strlen(latentString) + 1;

  counts.mLBlood = mLBlood;
  counts.lambda = parameters[0];
  counts.kappa = parameters[1];
  counts.d = parameters[2];
  counts.delta = parameters[3];
  counts.p = parameters[4];
  counts.c = parameters[5];

  counts.probLatent = probLatent;       /* eta, probLatent */
  counts.reactLatent = reactLatent;       /* alpha, reactLatent */
  counts.probDefect = probDefect;       /* gamma, probDefect */
  counts.latIncompDeath = latIncompDeath;    /* tau, latIncompDeath  */
  counts.latCompDeath = latCompDeath; /* sigma, latCompDeath */

  counts.numSampleTimes = sampleFileLength;
  counts.sampleCounter = sampleCounter;

  counts.maxActive = maxActive;
  counts.maxVirus = maxVirus;
  counts.maxIncomp = maxIncomp;
  counts.maxComp = maxComp;

  counts.totSampleVirus = totSampleVirus;
  counts.totSampleLatent = totSampleLatent;
  counts.timeLastSample = timeLastSample;
  counts.numEvents = numEvents;
  counts.ARTstart = ARTstart;

  char countsBiFile[50] = "\0";
  strcat(countsBiFile, outfile);
  strcat(countsBiFile, "counts.chkpt");

  writeCountsBinary(&counts, countsBiFile);
}

void writeLog(char* outfile, double * parameters, int sampleFileLength, long int numVirusInit,
		long int numCellUninfectInit, long int numCellInfectInit, long int numLatentCompInit,
	       long int	numLatentIncompInit, double * sampleTimes, int * numSampleVirus, 
	       int * numSampleLatent, double mLBlood, double probLatent, double reactLatent, 
	       double probDefect, double latIncompDeath, double latCompDeath, double lambdaInput, 
	       double kappaInput, double ARTStart) {
	
  char logFile[50] = "";
  strcat(logFile, outfile);
  strcat(logFile, "log.txt");
  FILE *file;
  file = fopen(logFile, "w");
  if (!file) {
	  fprintf("%s failed to open. Exiting the program.\n", logFile);
	  exit(1);
  }
  char line[100];

  sprintf(line, "mL Blood: \t%f\n", mLBlood);
  fputs(line, file);

  sprintf(line, "lambda: \t%f\n", lambdaInput);
  fputs(line, file);

  sprintf(line, "kappa: \t\t%.10f\n", kappaInput);
  fputs(line, file);

  sprintf(line, "d: \t\t%f\n", parameters[2]);
  fputs(line, file);

  sprintf(line, "delta: \t\t%f\n", parameters[3]);
  fputs(line, file);

  sprintf(line, "pi: \t\t%f\n", parameters[4]);
  fputs(line, file);

  sprintf(line, "c: \t\t%f\n", parameters[5]);
  fputs(line, file);

  sprintf(line, "eta: \t\t%f\n", probLatent);
  fputs(line, file);

  sprintf(line, "alpha: \t\t%f\n", reactLatent);
  fputs(line, file);
  
  sprintf(line, "gamma: \t\t%f\n", probDefect);
  fputs(line, file);

  sprintf(line, "sigma: \t\t%f\n", latCompDeath);
  fputs(line, file);

  sprintf(line, "tau: \t\t%f\n", latIncompDeath);
  fputs(line, file);



  //Sample times info, initial cell concentrations
  //
  sprintf(line, "\nInitial Conditions:\n");
  fputs(line, file);


  sprintf(line, "Initial number of viruses: %ld\n", numVirusInit);
  fputs(line, file);
  
  sprintf(line, "Initial number of uninfected cells: %ld\n", numCellUninfectInit);
  fputs(line, file);

  sprintf(line, "Initial number of infected cells: %ld\n", numCellInfectInit);
  fputs(line, file);

  sprintf(line, "Initial number of latent replication competent cells: %ld\n", numLatentCompInit);
  fputs(line, file);

  sprintf(line, "Initial number of latent replication incompetent cells: %ld\n", numLatentIncompInit);
  fputs(line, file);

  if (ARTStart)
  	sprintf(line, "\nART Initiation: %f\n", ARTStart);
  else 
  	sprintf(line, "\nNo ART\n");
  fputs(line,file);

  sprintf(line, "\nSample Times:\n");
  fputs(line, file);

  for(int i = 0; i < sampleFileLength; i++ ) {
  sprintf(line, "%f,%d,%d\n", sampleTimes[i], numSampleVirus[i], numSampleLatent[i]);
  fputs(line, file);

  }
  fclose(file);

  return;
}
