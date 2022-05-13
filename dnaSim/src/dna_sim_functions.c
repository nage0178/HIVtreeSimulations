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

#include "dna_sim.h"
#define _GNU_SOURCE
#include <ctype.h>
#include <math.h>

#define Sasprintf(write_to,  ...) {                     \
    char *tmp_string_for_extend = (write_to);           \
    int returnval = asprintf(&(write_to), __VA_ARGS__); \
    if (returnval < 0) {                                \
      fprintf(stderr, "Problem with sasprintf\n");      \
      exit(1);                                          \
    }                                                   \
    free(tmp_string_for_extend);                        \
}


/* Checks if two doubles are within a threshold value of each other */
bool FloatEquals(double a, double b, double threshold) {
  return fabs(a-b) < threshold;
}

/* Checks the input format for the parameters in the rate matrix */
void checkModelInput(int regMatch) {
  if (regMatch != 0) {
    fprintf(stderr, "Model parameters are not in the correct format.\nThe correct format is statFreq_A:statFreq_C:statFreq_G:statFreq-rateParameter_1:rateParameter_2:rateParameter_3:rateParameter_4:rateParameter_5:rateParameter_6\nNumbers must be in decimal format with a digit after the decimal. \nExiting the program. \n");
    exit(1);
  }
}

/* Create a new new node. Allocates memory for the node and the latency states
structure. Assignment values to the variables in the node structure */
struct NodeSampledTree* newNodeSampledTree(double branchLength, int num) {
  /* Allocates memory for a new node */
  struct NodeSampledTree* node_ptr;
  node_ptr  = (struct NodeSampledTree*)malloc(sizeof(struct NodeSampledTree));

  if (node_ptr == NULL) {
    fprintf(stderr, "Insufficient Memory\n");
    exit(1);
  }

  node_ptr->branchLength = branchLength;
  node_ptr->sample_time = 0;
  node_ptr->name = num;

  node_ptr->isLatent = 0;
  node_ptr->totTimeLatent = -1;

  /* Sets pointers to daughter nodes */
  node_ptr->left = NULL;
  node_ptr->right = NULL;
  return(node_ptr);
}

// Read in files

/* The function parses a string into a tree structure. Every time there's an open
parenthesis, two new nodes are created and the function is called twice. The
function keeps track of which character in the treeString (currChar) it is looking
at. */
int makeTree(struct NodeSampledTree* node_ptr, char* treeString, int currChar) {

  /* Note: You should never actually reach the end of string character. The last
  time you return is after the last closed parenthesis */
  while( treeString[currChar] != '\0') {

    /* If the current character is an open parenthesis, make new node, call
    makeTree recursively. currChar is updated along the way */
    if (treeString[currChar] == '(') {
      node_ptr->left = newNodeSampledTree(-1, -1);
      currChar = makeTree(node_ptr->left, treeString, currChar + 1);
      node_ptr->right = newNodeSampledTree(-1, -1);
      currChar = makeTree(node_ptr->right, treeString, currChar);


    } else {
      /* First, find the node name. Then, find the branch length. Return the
      current character */
      char* ptr; /* Needed for the strtod function */
      char* residualTreeString = malloc(strlen(treeString) + 1); /* Used so the original string isn't modified */
      strcpy(residualTreeString, treeString + currChar);

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

      node_ptr->name = atoi(nodeName);
      currChar = currChar + strlen(nodeName) + 1;

      residualTreeString = strcpy(residualTreeString, treeString + currChar);

      /* Finds the branch length, converts it to a double, saves it in the node
      structure */
      char* branchLength = strtok(residualTreeString, ",)");

      regMatch = regexec(&regDecimal, branchLength, 0, NULL, 0);
      if (regMatch != 0) {
        fprintf(stderr, "Problem reading in tree file. Regular expression does not match a decimal number.\n");
        exit(1);
      }
      node_ptr->branchLength = strtod(branchLength, &ptr);
      currChar = currChar + strlen(branchLength) + 1 ;

      free(residualTreeString);
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
nodeName(integer):totalTimeLatent,sampleTime,isLatent
isLatent is 1 if the tip is latent when sampled and 0 otherwise*/
void readLatent (struct NodeSampledTree* node, FILE* latentFile) {
  char buff[200];
  char* nodeName;
  char* ptr;
  double latentTime, sampleTime;
  int checkScan;

  checkScan = fscanf(latentFile, "%s", buff);
  if (checkScan == EOF) {
    fprintf(stderr, "Problem reading in latent states. fscanf returns EOF. Not enough lines in the file. \n");
    exit(1);
  }
  nodeName = strtok(buff, ":");

  /* Checks the node in the input file matches the current node */
  if (node->name != atoi(nodeName)) {
    fprintf(stderr, "Problem reading in latent states. States are not in the correct order or does not contain the correct nodes. \n");
    exit(1);
  }

  latentTime = strtod(buff + strlen(buff) + 1, &ptr);
  node->totTimeLatent = latentTime;

  sampleTime = strtod(ptr+1, &ptr);
  node->sample_time = sampleTime;

  node->isLatent = atoi(ptr + 1);
  if (node->isLatent != 0 && node->isLatent != 1) {
    fprintf(stderr, "Problem reading in latent states. The latent state is not a 0 or a 1. \n");
    exit(1);
  }

  if (node->left == NULL && node->right == NULL) {
    return ;

  } else {
    readLatent(node->left, latentFile);
    readLatent(node->right, latentFile);

  }
  return;
}

// Simulation Functions

/* Simulates substitions along a tree for a single position in a dna seqeunce.
Stores the base at tips in the alignment.
node: Pointer to the current node
base: Base at the node as an integer (0=A, 1=C, 2=G, 3=T)
instRate: Instantaneous rate matrix
r: random number generator
alignment: string array storing the alignment
baseIndex: Index of the position in the alignment that is being simulated
nodeIndex: Index for the node in the alignment. Ordering is based off the order
you each node is visited in a preorder traversal where the left node is visited
before the right node
nodeNames: pointer to an array to store the node names. names are stored in the
array in this function
*/
int simulateDnaTree(struct NodeSampledTree* node, int base, double instRate[4][4], const gsl_rng * r,
  char** alignment, int baseIndex, int nodeIndex, int* nodeNames, int* nodeLatentstate, double* nodeSampleTime, double *nodeLatentTime, double gammaDraw) {
    assert(node->totTimeLatent >= 0);
  /* If the node is a tip */
  if (node->left == NULL && node->right == NULL) {
    //printf("If %d\n", node->name);
    nodeNames[nodeIndex] = node->name;
    nodeLatentstate[nodeIndex] = node->isLatent;
    nodeSampleTime[nodeIndex] = node->sample_time;
    nodeLatentTime[nodeIndex] = node->totTimeLatent;

    base = simulateDNABranch(base, instRate, (node->branchLength - node->totTimeLatent) * gammaDraw,  r);

    /* Converts from integer base to the corresponding letter and saves in the
    alignment */
    if (base == 0) {
      alignment[nodeIndex][baseIndex] = 'A';
    } else if (base == 1) {
      alignment[nodeIndex][baseIndex] = 'C';
    } else if (base == 2) {
      alignment[nodeIndex][baseIndex] = 'G';
    } else if (base == 3 ){
       alignment[nodeIndex][baseIndex] = 'T';
    } else {
      fprintf(stderr, "Problem with simulating DNA \n");
      exit(1);
    }

    return (nodeIndex + 1);

  } else {
    int prevBase = base;
    /* Simulates substitions along the branch of the current node */
    base = simulateDNABranch(base, instRate, (node->branchLength - node->totTimeLatent) * gammaDraw,  r);
    if (node->branchLength - node->totTimeLatent == 0) {
      assert(prevBase ==base);
    }
    /* Recursively calls simulateDnaTree to simulate substitions along subtrees */
    /* nodeIndex is updated, everytime a tip is visited */
    nodeIndex = simulateDnaTree(node->left, base, instRate, r, alignment, baseIndex,
      nodeIndex, nodeNames, nodeLatentstate, nodeSampleTime, nodeLatentTime, gammaDraw); //I think the branch lengths are correct
    nodeIndex = simulateDnaTree(node->right, base, instRate, r, alignment, baseIndex,
      nodeIndex, nodeNames, nodeLatentstate, nodeSampleTime, nodeLatentTime, gammaDraw);

  }
return (nodeIndex);
}


/* Used to simulate DNA substitions over a single branch with a single base.
base: starting base (0=A, 1=C, 2=G, 3=T, used for indexing convenience)
instRate: Instantaneous rate matrix, this is completely general. It doesn't need
modifications for different substition models.
branchLength: branch length
r: random number generator
Returns: the base at the end of the branch */
int simulateDNABranch(int base, double instRate[4][4], double branchLength, const gsl_rng * r) {
  double time = 0;        /* The time that has passed in the simulation */
  double totRate, unif;   /* Total rate of substitions, depends on starting state.
                          unif is used with the random number generator to pick the next base */
  double probBase[4];     /* Given an substition occured, the probability the new base is 0-3 (A-T)*/

  /* Until the end of the branch is reached */
  while (time <= branchLength) {
    /* Rate depends on starting state */
    totRate = -instRate[base][base];

    /* Time of next substition*/
    time = time + gsl_ran_exponential (r, 1 / totRate);

    if (time > branchLength) {
      return base;
    }

    /* Finds the probability of each substition */
    for (int i = 0; i < 4; i++) {
      if (i == base) {
        probBase[i] = 0.0;
      } else {
        probBase[i] = instRate[base][i] / totRate;
      }
    }
    assert(FloatEquals(probBase[0] + probBase[1] + probBase[2] + probBase[3], 1, 1e-4));

    /* Updates the current base */
    unif = gsl_rng_uniform (r);
    if (unif <= probBase[0]) {
      base = 0;
    } else if (unif <= probBase[0] + probBase[1]) {
      base = 1;
    } else if (unif <= probBase[0] + probBase[1] + probBase[2]) {
      base = 2;
    } else {
      base = 3;
    }
  }
  return base;
}

/* Counts the number of tips on the tree structure.*/
int countTips(struct NodeSampledTree* node, int count) {
  if (node->left == NULL && node->right == NULL) {
    return(count + 1);

  } else {
    count = countTips(node->left, count);
    count = countTips(node->right, count);
    return count;
  }
}

/* Prints the simulated alignment in fastsa format */
void printAlignment(char** alignment, int numSeq, int* nodeName, int* nodeLatentstate, double* sampleTime, double* nodeLatentTime, char* outfile) {
  char outFilename[120] = "\0";
  strcat(outFilename, outfile);
  strcat(outFilename, ".fa");

  FILE *file;
  file = fopen(outFilename, "w");
  for (int i = 0; i < numSeq; i++) {
    if (nodeName[i] == -2) {
      fprintf(file, ">outgroup\n%s\n",alignment[i]);

    } else {
      fprintf(file, ">Node_%d_%d_%0.0f_%0.0f\n%s\n", nodeLatentstate[i], nodeName[i], sampleTime[i] - nodeLatentTime[i],sampleTime[i], alignment[i]);

    }

  }
  fclose(file);
}


/* This function prints the tree structure in newick format.
This can be used to check that the file is identical to the file that was
read in. A different function is required than the function originally used to
print the tree since the birth age is actually a branch length when it is read
n.*/
void printTree(struct NodeSampledTree* root_ptr, char* fileName) {
  char *treeStringCheck = strdup("");
  treeStringCheck = NewickStringBranchLength(root_ptr, treeStringCheck, 0);

  FILE *treeFile;
  treeFile = fopen(fileName, "w");
  fputs(treeStringCheck, treeFile);
  fclose(treeFile);
  free(treeStringCheck);
}

/* Used to print the tree to test that it matches the tree that was read in. See
printTree for more detail. */
char* NewickStringBranchLength(struct NodeSampledTree* node, char* treeString, double previousTime) {

  if (node->left == NULL && node->right == NULL) {
    Sasprintf(treeString, "%s%d:%f", treeString, node->name, node->branchLength);
    return (treeString);

  } else {
    Sasprintf(treeString, "%s(", treeString);
    treeString = NewickStringBranchLength(node->left, treeString, node->branchLength);
    Sasprintf(treeString, "%s,", treeString);
    treeString = NewickStringBranchLength(node->right, treeString, node->branchLength);
    Sasprintf(treeString, "%s)%d:%f", treeString, node->name, node->branchLength);

  }
  return (treeString);
}


/* This function prints the tree structure in newick format.
The names match those in the alignment. A different function
is required than the function originally used to print the tree
since the birth age is actually a branch length when it is read
in.*/
void printTreeNewNames(struct NodeSampledTree* root_ptr, char *outfile, int outgroup) {
  /* Writes tree file without outgroup */
  char outFilename[120] = "\0";
  strcat(outFilename, outfile);
  strcat(outFilename, "Tree.txt");

  char *treeStringCheck = strdup("");

  if (outgroup) {
    treeStringCheck = treeNewNames(root_ptr->left, treeStringCheck, 0);

  } else {
    treeStringCheck = treeNewNames(root_ptr, treeStringCheck, 0);
  }

  Sasprintf(treeStringCheck , "%s;", treeStringCheck);

  FILE *treeFile;
  treeFile = fopen(outFilename, "w");
  fputs(treeStringCheck, treeFile);
  fclose(treeFile);
  free(treeStringCheck);

 /* Writes tree file with outgroup */
  if (outgroup) {
    outFilename[0] = '\0';
    strcat(outFilename, outfile);
    strcat(outFilename, "TreeOG.txt");

    char *treeStringOG = strdup("");

    treeStringOG = treeNewNames(root_ptr, treeStringOG, 0);

    Sasprintf(treeStringOG , "%s;", treeStringOG);

    FILE *treeFileOG;
    treeFileOG = fopen(outFilename, "w");
    fputs(treeStringOG, treeFileOG);
    fclose(treeFileOG);
    free(treeStringOG);

  }
}

/* Find the newick format of the tree with names that match the alignment */
char* treeNewNames(struct NodeSampledTree* node, char* treeString, double previousTime) {

  if (node->name == -2) {
    Sasprintf(treeString, "%soutgroup:%f", treeString, node->branchLength);
    return (treeString);
  } else if (node->left == NULL && node->right == NULL) {
    Sasprintf(treeString, "%sNode_%d_%d_%0.0f_%0.0f:%f", treeString, node->isLatent, node->name, node->sample_time - node->totTimeLatent, node->sample_time, node->branchLength);
    return (treeString);

  } else {
    Sasprintf(treeString, "%s(", treeString);
    treeString = treeNewNames(node->left, treeString, node->branchLength);
    Sasprintf(treeString, "%s,", treeString);
    treeString = treeNewNames(node->right, treeString, node->branchLength);
    Sasprintf(treeString, "%s)Node_%d_%d_%0.0f_%0.0f:%f", treeString, node->isLatent, node->name, node->sample_time - node->totTimeLatent, node->sample_time, node->branchLength);

  }
  return (treeString);
}

/* Find the newick format of the tree with names that match the alignment */
void findLastSample(struct NodeSampledTree* node, double* lastSample) {

  if (node->left == NULL && node->right == NULL) {
    if (node->sample_time > *lastSample) {
      *lastSample = node->sample_time;
    }

  } else {
    findLastSample(node->left, lastSample);
    findLastSample(node->right, lastSample);

  }
  return;
}

/* Clears memory for the tree. Stem must be cleared separately */
void clearTree(struct NodeSampledTree* node_ptr) {
  if (node_ptr->left != NULL) {
    clearTree(node_ptr->left);
    clearTree(node_ptr->right);
  }

  /* Clear node */
  free(node_ptr);
  return;
}

/* Reads in a fasta file and sets it as the starting sequence*/
int* readFasta(char* fastaSeq, int* numBases) {
  int count = 0;

  /* Keep going until you reach the end of the line with the ">" */
  while(*fastaSeq!= '\n') {
    fastaSeq++;
  }
  fastaSeq++;

  /* Count the number of bases in the file. Checks that the characters are DNA bases*/
  char* tmp = fastaSeq;
  while(*tmp != '\0'){
    if (*tmp == 'A' || *tmp == 'C'|| *tmp == 'G' || *tmp == 'T'
        || *tmp == 'a' || *tmp == 'c' || *tmp == 'g' || *tmp == 't') {
        count++;

    }  else if (!isspace(*tmp)) {
      fprintf(stderr, "Problem with fasta file. There are characters other than bases and whitespace in the sequence. Ambiguity codes are not accomodated. \n");
      exit(1);
    }
    tmp++;
  }

  if (*numBases != count && *numBases != 0) {
    fprintf(stderr, "The number of bases in the alignment does not match the input number. The whole alignment will be used.\n");
  }

  *numBases = count;

  /* Sets the starting sequence to be the same as the input file */
  int i = 0;
  int* startSeq = malloc(sizeof(int) * count);
  while(i < count) {

    if (*fastaSeq == 'A' || *fastaSeq == 'a'){
      startSeq[i] = 0;
    } else if (*fastaSeq == 'C' || *fastaSeq == 'c'){
      startSeq[i] = 1;
    } else if (*fastaSeq == 'G' || *fastaSeq == 'g'){
      startSeq[i] = 2;
    } else if (*fastaSeq == 'T' || *fastaSeq == 't'){
      startSeq[i] = 3;
    } else {
      i--;
    }
    i++;
    fastaSeq++;
  }

  return startSeq;
}

/* Writes the starting seed to a file*/
void writeSeed(unsigned int RGSeed, char* outfile) {
  char filename[120] = "\0";
  strcat(filename, outfile);
  strcat(filename, "seed.txt");

  char seed[15];
  sprintf(seed, "%u", RGSeed);

  FILE *file;
  file = fopen(filename, "w");
  fputs(seed, file);
  fclose(file);

  return;
}

/* Sets the starting sequences. If randomStart is true, then the bases are drawn from the
stationary frequencies. Otherwise, a starting sequence is read in from a fasta file.*/
int* setStartSeq (int* numBases, int randomStart, double statFreq[], char* fastaSeq, gsl_rng *r) {
  double unif;
  int* startSeq;
  if(randomStart) {
    //If the user did not provide the length of the sequence
    if (*numBases == 0) {
      *numBases = 500; //Default number of bases
    }
    startSeq = malloc(sizeof(int) * *numBases);

    // Assigns starting bases from the stationary distribution
    for (int indBase = 0; indBase < *numBases; indBase ++) {
      unif = gsl_rng_uniform (r);
      if (unif <= statFreq[0]) {
        startSeq[indBase] = 0;
      } else if (unif <= statFreq[0] + statFreq[1]) {
        startSeq[indBase] = 1;
      } else if (unif <= statFreq[0] + statFreq[1] + statFreq[2]) {
        startSeq[indBase] = 2;
      } else {
        startSeq[indBase] = 3;
      }
    }

  } else {
    // If the user provided a starting sequence, assign it as the starting sequence
    startSeq = readFasta(fastaSeq, numBases);
  }

  return startSeq;
}

//Sets the seed based on the system time if the user did not provide a seed.
// Writes the seed to a file
void setSeed(int RGSeed, gsl_rng *r, char* outfile) {
  if (RGSeed == 0) {
    RGSeed = time(0);
  }
  gsl_rng_set(r, RGSeed);
  writeSeed(RGSeed, outfile);

  return;
}

/* Makes an alignment based on the tree, starting sequence, and instantaneous rate matrix.
Puts the node names, latent state, and sample times into arrays. This is so this information can be printed in the sequence name.
Each site is simulated independently. */
void makeAlignment(struct NodeSampledTree* root_ptr, int startSeq[], double instRate[4][4], char** alignment, int nodeNames[], int nodeLatentstate[], double nodeSampleTime[], double nodelatentTime[], int numBases, int numTips, gsl_rng *r, double alpha, int gammaRateVariation) {
  int checkSimOutput; /* Checks the dna simulator return the correct nodeIndex.*/

  double gammaDraw = 1;

  for (int indBase = 0; indBase < numBases; indBase ++) {
    if (gammaRateVariation) {
      gammaDraw = gsl_ran_gamma(r, alpha, 1/alpha);
    }

    checkSimOutput = simulateDnaTree(root_ptr, startSeq[indBase], instRate, r, alignment, indBase, 0, nodeNames, nodeLatentstate, nodeSampleTime, nodelatentTime, gammaDraw);
    assert(checkSimOutput == numTips);
  }
  return;
}

/* Frees some of the memory in the program. (There's still a leak)*/
void clearMemory(struct NodeSampledTree* root_ptr, int numTips, int* nodeLatentstate, int* nodeNames, double* nodeSampleTime, double* nodeLatentTime, char** alignment, char* treeString, gsl_rng *r) {
  for (int i = 0; i < numTips; i++ ) {
    free(alignment[i]);
  }

  clearTree(root_ptr);
  free(alignment);
  free(nodeNames);
  free(nodeLatentstate);
  free(nodeLatentTime);

  free(treeString);
  gsl_rng_free (r);

  return;
}

/* Allocates memory for the alignment. The alignment is an array of strings.
Each string in the number of bases (plus the null byte)*/
char** allocateAlignmentMem(int numTips, int numBases) {
  char** alignment = malloc(sizeof(char*) * numTips);
  if (alignment == NULL) {
    fprintf(stderr, "Insufficient Memory\n");
    exit(1);
  }

  /* Allocates memory for every sequence in the alignment.
  Sets the end of the string character for every sequence in the alignment */
  for (int i = 0; i < numTips; i++ ) {
    alignment[i] = malloc(sizeof(char) * numBases + 1);

    if (alignment[i] == NULL) {
      fprintf(stderr, "Insufficient Memory\n");
      exit(1);
    }
    alignment[i][numBases] = '\0';
  }

  return alignment;
}

/* Makes the instantaneous rate matrix given the stationary frequencies, substition rate, and rate parameters (in instRate)*/
void makeInstantaneousRate(double instRate[4][4], double statFreq[], double mu) {

  double diag,  weightSum  = 0;
  /* Finds the diagonal elements and the total substition rate*/
  for (int i = 0; i < 4; i++) {
    instRate[i][i] = 0;
    diag = 0;
    for (int j = 0; j < 4; j ++) {
      diag = diag + instRate[i][j];
    }
    weightSum = diag * statFreq[i] + weightSum;
    instRate[i][i] = -diag;

  }
  if (weightSum != 1) {
    printf("Rescaling instantaneous rate matrix so the average substition rate is one.\nThen multiplying by the substituion rate. \n");
  }
  if (weightSum == 0) {
    fprintf(stderr, "The average substition rate is zero. Exiting. \n");
    exit(1);
  }

  /* Rescales the instantaneous rate matrix so that the average substitution rate is equal to mu */
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j ++) {
      instRate[i][j] = instRate[i][j] / weightSum *  mu;
    }
  }

  return;
}

/* Checks the input parameters are valid */
void checkInput(double statFreq[], double rateParams[], double mu) {
  if (!FloatEquals(statFreq[0] + statFreq[1] + statFreq[2] + statFreq[3], 1, 1e-4)) {
    fprintf(stderr, "Stationary frequencies do not sum to 1. Check -f inputs. \n");
    exit(1);
  }

  double sum = 0;
  for (int i = 0; i < 6; i++) {
    sum = sum + rateParams[i];
    if (rateParams[i] < 0) {
      fprintf(stderr, " Transitions rates must be non-negative. Check -f inputs. \n");
      exit(1);
    }
  }

  if (sum <= 0) {
    fprintf(stderr, "At least one of the transition rates must be positive. Check -f inputs. \n");
    exit(1);
  }

  if (mu <= 0) {
    fprintf(stderr, "The mutation rate must be positive. Check -u input.  \n");
    exit(1);
  }
}
