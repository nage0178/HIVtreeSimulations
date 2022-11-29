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
#include "string_utilities.h"

gsl_rng * r;  /* global generator */

int main(int argc, char **argv) {

  /* Prints the command that the program is running with */
  printf("Command: \n");
  for (int i =0; i < argc; i ++ ) {
    printf("%s ", argv[i]);

  }
  printf("\n\n");

  /* Set up random number generator */
  const gsl_rng_type * T;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);

  int c; /* Used for switch statement for getopt */
  char* endPtr; /* Needed for strtoul function*/
  char* token; /* Used to tokenize input strings */

  unsigned int RGSeed = 0;   /* Seed for random number generator */
  int randomStart = 1;       /* Determines if the sequence at the root
                             is given or choosen from the stationary distribution*/
  int numBases = 1000;       /* Number of bases in the sequence */
  double mu = 3.6*1e-5;     /* Per base mutation rate, empirically this is ballpark 1.8-3.6 *10^-5 per day */
  double statFreq[4] = {.25, .25, .25, .25}; /* Stationary frequencies, of A, C, G, T */
  double rateParams[6] = {1, 1, 1, 1, 1, 1}; /* Parameters in the instantaneous rate matrix */
  double alpha = 1;
  int gammaRateVariation = 0;
  int outgroup = 0;
  int outFileNameGiven = 0;

  /* Used to confirm the required input files were given */
  int treeFileEntered = 0;
  int latentHistoryFileEntered = 0;

  char* outfile  = malloc(sizeof(char)*100);          /* Prefix of output file */
  char* treeFile = malloc(sizeof(char)*100);          /* Name of tree file, output from HIV_final_sim.c*/
  char* latentHistoryFile = malloc(sizeof(char)*100); /* Name of latent history file, output from HIV_final_sim.c*/
  char* inputFasta = malloc(sizeof(char)*30);         /* Name of file containing the ancestral sequence, optional*/
  char *fastaSeq =  NULL;                             /* Starting dna seqeunce */

  regex_t regDecimal, regInt;
  int regCompiled, regMatch;
  /* Regex for decimal number */
  regCompiled = regcomp(&regDecimal, "^([0-9]*)((.)([0-9]+))$" , REG_EXTENDED);
  if (regCompiled == 1) {
    fprintf(stderr, "Regular expression did not compile.\n");
    exit(1);
  }
  /* Regex for integer number */
  regCompiled = regcomp(&regInt, "^([0-9]+)$" , REG_EXTENDED);
  if (regCompiled == 1) {
    fprintf(stderr, "Regular expression did not compile.\n");
    exit(1);
  }

  /* Process command line arguements */
  while((c = getopt(argc, argv, "b:f:i:s:u:t:l:a:o:hr")) != -1){
    switch(c) {
      /* Number of bases */
      case 'b':
        token = strtok(optarg, " ");
        regMatch = regexec(&regInt, token, 0, NULL, 0);
        if (regMatch != 0) {
          fprintf(stderr, "Number of bases is in the incorrect format. The number must be an integer. \nExiting the program. \n");
          exit(1);
        }
        numBases = atoi(optarg);
        break;

      /* Stationary frequencies followed by rate matrix parameters */
      case 'f':
        token = strtok(optarg, ":");
        regMatch = regexec(&regDecimal, token, 0, NULL, 0);
        checkModelInput(regMatch);
        statFreq[0] = atof(token);

        for (int i = 1; i < 3; i++) {
          token = strtok(NULL, ":");
          regMatch = regexec(&regDecimal, token, 0, NULL, 0);
          checkModelInput(regMatch);
          statFreq[i] = atof(token);
        }

        token = strtok(NULL, "-");
        regMatch = regexec(&regDecimal, token, 0, NULL, 0);
        checkModelInput(regMatch);
        statFreq[3] = atof(token);

        for (int i = 0; i < 6; i++) {
          token = strtok(NULL, ":");
          regMatch = regexec(&regDecimal, token, 0, NULL, 0);
          checkModelInput(regMatch);
          rateParams[i] = atof(token);
        }
        break;

      /* File for ancestral DNA sequence */
      case 'i':
        strncpy(inputFasta, optarg, 30);
        fastaSeq = string_from_file(inputFasta);
        randomStart = 0;
        break;

      /* Tree file from HIV_final_sim.c  */
      case 't':
        strncpy(treeFile, optarg, 100);
        treeFileEntered = 1;
        break;

      /* Latent history file from HIV_final_sim.c  */
      case 'l':
        strncpy(latentHistoryFile, optarg, 100);
        latentHistoryFileEntered = 1;
        break;

      /* Starting seed */
      case 's':
        token = strtok(optarg, " ");
        regMatch = regexec(&regInt, token, 0, NULL, 0);
        if (regMatch != 0) {
          fprintf(stderr, "The seed is in the incorrect format. The number must be an integer. \nExiting the program. \n");
          exit(1);
        }
        RGSeed = strtoul(optarg, &endPtr,10);
        break;

      /* Substitution rate */
      case 'u':
        token = strtok(optarg, " ");
        regMatch = regexec(&regDecimal, token, 0, NULL, 0);
        if (regMatch != 0) {
          fprintf(stderr, "Substitution rate is not in the correct format. \nIt must be a decimal number with a digit after the decimal.\nExiting the program.\n");
          exit(1);
        }
        mu = atof(token);

        break;
      case 'a':
        token = strtok(optarg, " ");
        regMatch = regexec(&regDecimal, token, 0, NULL, 0);
        if (regMatch != 0) {
          fprintf(stderr, "alpha parameter of gamma distribution is not in the correct format. \nIt must be a decimal number with a digit after the decimal.\nExiting the program.\n");
          exit(1);
        }
        alpha = atof(token);
        gammaRateVariation = 1 ;
        break;

      case 'o':
        outFileNameGiven = 1;
        strncpy(outfile, optarg, 100);
        break;

      case 'r':
        outgroup = 1;
        break;

     case 'h':
      printHelp();
      exit(0);
      break;

    default:
      abort();
    }
  }

  if (!latentHistoryFileEntered || !treeFileEntered) {
    fprintf(stderr, "Did not provide the two required input files, the tree file and the latent history file.\n");
    printHelp();
    fprintf(stderr, "Exiting the program.\n");
    exit(1);
  }

  if (!outFileNameGiven) {
    fprintf(stderr, "Did not provide the output file name.\n");
    printHelp();
    fprintf(stderr, "Exiting the program.\n");
    exit(1);
  }

  checkInput(statFreq, rateParams, mu); /* Checks input parameters are valid */
  setSeed(RGSeed, r, outfile);
  int* startSeq = setStartSeq (&numBases, randomStart, statFreq, fastaSeq, r);

  /* Read in newick format file. File needs to have internal node labels and
  branch lengths */
  char *treeString = string_from_file(treeFile);

  /* Make a tree from the input tree file */
  struct NodeSampledTree* root_ptr = newNodeSampledTree(0, 1);
  int endChar = makeTree(root_ptr, treeString, 0);
  if (treeString[endChar] != '\0') {
    fprintf(stderr, "Problem reading in tree. Did not reach the end of the file \n");
    exit(1);
  }

  /* Add the latent history to the tree from the input file */
  FILE *latentFile;
  latentFile = fopen(latentHistoryFile, "r");
  readLatent(root_ptr, latentFile);
  fclose(latentFile);

  if (outgroup) {
    struct NodeSampledTree* new_root_ptr = newNodeSampledTree(0, -1);
    new_root_ptr->totTimeLatent = 0;
    double lastSample = 0;
    findLastSample(root_ptr, &lastSample);

    struct NodeSampledTree* outgroup = newNodeSampledTree(lastSample, -2);
    outgroup->totTimeLatent = 0;
    outgroup-> sample_time = lastSample;

    new_root_ptr->left = root_ptr;
    new_root_ptr->right = outgroup;
    root_ptr = new_root_ptr;
  }

  /* Set the values of the instantaneous rate matrix. Sets diagonal and rescales to makes the mutation rate in the function. */
  double instRate[4][4] = { {                          0, rateParams[0] * statFreq[1],    rateParams[1] * statFreq[2],  rateParams[2] * statFreq[3]},
                {rateParams[0] * statFreq[0],                           0,    rateParams[3] * statFreq[2],  rateParams[4] * statFreq[3]},
                {rateParams[1] * statFreq[0], rateParams[3] * statFreq[1],                              0,  rateParams[5] * statFreq[3]},
                {rateParams[2] * statFreq[0], rateParams[4] * statFreq[1],    rateParams[5] * statFreq[2],                           0}};
  makeInstantaneousRate(instRate, statFreq, mu);

  int numTips = countTips(root_ptr, 0); /* Number of tips in the tree */

  /* Needed for naming the sequences in the alignment. */
  int* nodeNames         = malloc(sizeof(int) * numTips);
  int* nodeLatentstate   = malloc(sizeof(int) * numTips);
  double* nodeSampleTime = malloc(sizeof(double) * numTips);
  double* nodeLatentTime = malloc(sizeof(double) * numTips);
  /* Allocates the memory required for the alignment */
  char** alignment = allocateAlignmentMem(numTips, numBases);

  /* Generates and prints an alignment */
  makeAlignment(root_ptr, startSeq, instRate, alignment, nodeNames, nodeLatentstate, nodeSampleTime, nodeLatentTime, numBases, numTips, r, alpha, gammaRateVariation);
  printAlignment(alignment, numTips, nodeNames, nodeLatentstate, nodeSampleTime, nodeLatentTime, outfile);
  printTreeNewNames(root_ptr, outfile, outgroup);

  clearMemory(root_ptr, numTips, nodeLatentstate, nodeNames, nodeSampleTime, nodeLatentTime, alignment, treeString, r);
  free(treeFile);
  free(latentHistoryFile);
  return(0);
}
