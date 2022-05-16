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

gsl_rng * r;  /* global generator */

int main(int argc, char **argv) {
  /* Prints the command that the program is running with */
  printf("Command: \n");
  for (int i =0; i < argc; i ++ ) {
    printf("%s ", argv[i]);

  }
  printf("\n\n");

  double checkPointTime = -1;
  char* sampleTimesFile = malloc(sizeof(char)*100);
  int sampleFileLength;
  double timeLastSample;
  char* controlFile     = malloc(sizeof(char)*100);
  char* loadName        = malloc(sizeof(char)*100);

  char* outfile           = malloc(sizeof(char)*100);
  outfile[0] = '\0';
  char countsFile[100]     = "\0";
  char gslFile[100]        = "\0";
  char latentFileRead[100] = "\0";
  char treeFileRead[100]   = "\0";

  unsigned long int totMem = 0; /* Memory currently allocated in the simulation */
  totMem = sizeof(char) * 100 * 3;
  unsigned int RGSeed = 0;

  int printFreq = 100000000;
  int print = 1;   /* If true (nonzero), prints population size every printFreq events*/
  int restart = 1; /* If nonzero, restart the simulation if the population goes extinct. Used for debugging. */

  char* endPtr; /* Used in strtoul function */
  int c;        /* Used in switch statement for command line arguements */

  /*Set up random number generator */
  const gsl_rng_type * T;
  gsl_rng_env_setup();
  T = gsl_rng_taus2;  /*gsl_rng_default*/
  r = gsl_rng_alloc(T);

  /* Regular expression parts*/
  regex_t regDecimal, regInt;
  int regCompiled, regMatch, regMatch2;
  /* Regular expression for decimal number */
  regCompiled = regcomp(&regDecimal, "^([0-9]*)((\\.)([0-9]+))$" , REG_EXTENDED);
  if (regCompiled == 1) {
    fprintf(stderr, "Regular expression did not compile.\n");
    exit(1);
  }
  /* Regular expression for integer */
  regCompiled = regcomp(&regInt, "^([0-9]+)$" , REG_EXTENDED);
  if (regCompiled == 1) {
    fprintf(stderr, "Regular expression did not compile.\n");
    exit(1);
  }

  /* Simulation parameters */
  double mLBlood = 5; /* mL of blood in the simulation.*/
  double parameters[6];
  parameters[0] = 1e2;                  /* lambda, uninfectedBirth */
  parameters[1] = 4.6 * 1e-7;           /* kappa, uninfectToInfect */
  parameters[2] = 0.01;                 /* d, uninfectDeath*/
  parameters[3] = .4;                    /* delta, infectedDeath */
  parameters[4] = 6500;                 /* p (or pi), virusBirth */
  parameters[5] = 23;                   /* c, virusDeath */

  double probLatent       = .000045;    /* eta, probLatent */
  double reactLatent      = .000057;    /* alpha, reactLatent */
  double probDefect       = 0.98;       /* gamma, probDefect */
  double latIncompDeath   = 0.00011;    /* tau, latIncompDeath  */
  double latCompDeath     = 0.00052;    /* sigma, latCompDeath */

  /* Used to ensure the command line arguement is used if values are given in the
  command line and in the control file*/
  int seedChange = 0, volChange = 0, sampleFileChange = 0;

  /* Checkpoint related */
  /* Used to save all variables and parameters not related to the tree for a checkpoint */
  struct cellCountTimes* counts = malloc(sizeof(struct cellCountTimes));
  char* treeFilePointer = NULL;   /*Pointer to binary tree file from checkpoint  */
  char* latentFilePointer = NULL; /*Pointer to binary latent history file from checkpoint  */
  int loadCheckPoint = 0;         /* Was the checkpoint loaded from previous time? Yes = 1, no = 0 */
  int startReload = 0;            /* Used so that an extra random number isn't drawn.
                                     Use the number from before the checkpoint. */

  /* Stem will point to the root, needed for pruning the tree to work*/
  struct Node* stem = newNode(0, 0, NULL);
  totMem = totMem + sizeof(struct Node);
  struct Node* root_ptr;

  /* Current length of each array*/
  long int maxActive;
  long int maxVirus;
  long int maxIncomp;
  long int maxComp;

  /* Total number of virions and latent viruses to sample
  Also the length of the arrays for sampling */
  int totSampleVirus;
  int totSampleLatent;
  /* Number of viruses that have bene sampled. Used as an
  index to place sampled viruses into the sampled array */
  int numLatentSample;
  int numVirusSample;

  /* Store viruses that are tips  of each class */
  struct Node** activeArray;
  struct Node** virusArray;
  struct Node** latentIncompArray;
  struct Node** latentCompArray;
  struct Node** virusSampleArray;
  struct Node** latentSampleArray;

  while((c = getopt(argc, argv, "s:i:o:r:p:v:c:l:t:h")) != -1) {
    switch(c) {
      /* Control file */
      case 'c':
        strncpy(controlFile, optarg, 100);
        parseControlFile(controlFile, sampleTimesFile, &mLBlood, parameters,
          &probLatent, &reactLatent, &probDefect, &latIncompDeath, &latCompDeath,
          &RGSeed, &seedChange, &volChange, &sampleFileChange);
        break;
      /* Seed */
      case 's':
        regMatch = regexec(&regInt, optarg, 0, NULL, 0);
        if (regMatch != 0) {
          fprintf(stderr, "The seed is in the incorrect format. The number must be an integer. \nExiting the program. \n");
          exit(1);
        }
        RGSeed = strtoul(optarg, &endPtr, 10);
        seedChange++;
        break;
      /* Input file */
      case 'i':
        strncpy(sampleTimesFile, optarg, 100);
        sampleFileChange++;
        break;
      /* Output file */
      case 'o':
        strncpy(outfile, optarg, 100);
        break;
      /* Print frequency */
      case 'p':
        regMatch = regexec(&regInt, optarg, 0, NULL, 0);
        if (regMatch != 0) {
          fprintf(stderr, "The print frequency must be an integer. \nExiting the program. \n");
          exit(1);
        }
        printFreq = atoi(optarg);
        if(printFreq <= 0 ) {
          print = 0;
        }
        break;
      /* Restarting if extinction */
      case 'r':
        regMatch = regexec(&regInt, optarg, 0, NULL, 0);
        if (regMatch != 0) {
          fprintf(stderr, "The restart flag must be an integer. \nExiting the program. \n");
          exit(1);
        }
        restart = atoi(optarg);
        break;
      /* Volume of simulation */
      case 'v':
        regMatch = regexec(&regDecimal, optarg, 0, NULL, 0);
        regMatch2 = regexec(&regInt, optarg, 0, NULL, 0);
        if (regMatch != 0 && regMatch2 != 0) {
          fprintf(stderr, "The volume is not in the correct format. \nIt must be a decimal number with a digit after the decimal or an integer.\nExiting the program.\n");
          exit(1);
        }
        mLBlood = atof(optarg);
        volChange++;
        break;
      /* Changes the time of the checkpoint from the default time */
      case 't':
        regMatch = regexec(&regDecimal, optarg, 0, NULL, 0);
        regMatch2 = regexec(&regInt, optarg, 0, NULL, 0);
        if (regMatch != 0 && regMatch2 != 0) {
          fprintf(stderr, "The checkpoint time is not in the correct format. \nIt must be a decimal number with a digit after the decimal or an integer.\nExiting the program.\n");
          exit(1);
        }
        checkPointTime = atof(optarg);
        break;
      /* Load files from previous run */
      case 'l':
        loadCheckPoint = 1;
        startReload = 1;

        /* loadName is the prefix for all files from the checkpoint */
        strncpy(loadName, optarg, 100);

        strcat(countsFile, loadName);
        strcat(countsFile, "counts.chkpt");

        strcat(gslFile, loadName);
        strcat(gslFile, "gsl.chkpt");

        strcat(latentFileRead, loadName);
        strcat(latentFileRead, "latent.chkpt");

        strcat(treeFileRead, loadName);
        strcat(treeFileRead, "tree.chkpt");

        /* Reads in all of the files from the checkpoint */
        if(readInBinaries(countsFile, gslFile, latentFileRead, treeFileRead,
          counts, &treeFilePointer, &latentFilePointer, r) == EXIT_FAILURE) {
            fputs("Problem loading checkpoint files\nExiting\n", stderr);
            exit(1);
        }

        /* Sets the size of the viral arrays */
        maxActive = counts->maxActive;
        maxVirus = counts->maxVirus;
        maxIncomp = counts->maxIncomp;
        maxComp = counts->maxComp;
        totSampleVirus  = counts->totSampleVirus;
        totSampleLatent = counts->totSampleLatent;
        numVirusSample = counts->numVirusSample;
        numLatentSample = counts->numLatentSample;

        /* Makes the viral arrays to store the tips of the tree*/
        activeArray       = malloc(sizeof(struct Node*) * maxActive);
        virusArray        = malloc(sizeof(struct Node*) * maxVirus);
        latentIncompArray = malloc(sizeof(struct Node*) * maxIncomp);
        latentCompArray   = malloc(sizeof(struct Node*) * maxComp);
        virusSampleArray  = malloc(sizeof(struct Node*) * totSampleVirus);
        latentSampleArray = malloc(sizeof(struct Node*) * totSampleLatent);

        if (activeArray == NULL || virusArray == NULL || latentIncompArray == NULL || latentCompArray == NULL || virusSampleArray == NULL || latentSampleArray == NULL) {
          fprintf(stderr, "Insufficient Memory\n");
          exit(1);
        }

        /* Sets the array elements to null*/
        setVirusArrays(activeArray, virusArray, latentIncompArray, latentCompArray, virusSampleArray, latentSampleArray,
          maxActive, maxVirus, maxIncomp, maxComp, totSampleVirus, totSampleLatent);

        /* Creates the tree from before the checkpoint */
        root_ptr = newNodeSampledTree(stem, &totMem);
        int endChar = makeTree(root_ptr, treeFilePointer, 0, &totMem);
        if (treeFilePointer[endChar] != '\0') {
          printf("%s\n", treeFilePointer +endChar);
          fprintf(stderr, "Problem reading in tree. Did not reach the end of the file \n");
          exit(1);
        }

        /* Sets all of the variables in the node struct for the tree from before the checkpoint */
        char* ptr = latentFilePointer;
        readLatent(root_ptr, latentFilePointer, activeArray, virusArray, latentIncompArray, latentCompArray,
        virusSampleArray, latentSampleArray, &ptr);

        free(treeFilePointer);
        free(latentFilePointer);

        break;
      case 'h':
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
        exit(1);
        break;
      default:
        abort();
    }
  }

/* If the run is starting fresh (no checkpoint)*/
if (loadCheckPoint == 0) {
    if (seedChange > 1 || volChange > 1 || sampleFileChange > 1) {
	fprintf(stderr, "Options cannot be specified multiple times. Check no options are specified both in the command line and in the control file. Exiting.\n");

	exit(1);
    }

  /* Ensuring units match the volume of blood*/
  parameters[0] = parameters[0] * mLBlood;
  parameters[1] = parameters[1] / mLBlood;

  /* Use system time for seed if no seed is given */
  if (RGSeed == 0) {
    RGSeed = time(0);
  }
  gsl_rng_set(r,RGSeed);
  writeSeed(RGSeed, outfile);

  /* Empirical measurements of kappa are for productively infected cells.
  We need to take this into account */
  parameters[1] = parameters[1] / (1 - probDefect) / (1 - probLatent);

 // Sampling
  sampleFileLength = findFileLength(sampleTimesFile);

} else {


  // Sets parameters to be equal to those used before checkpoint
    mLBlood = counts->mLBlood;
    parameters[0] = counts->lambda;
    parameters[1] = counts->kappa;
    parameters[2] = counts->d;
    parameters[3] = counts->delta;
    parameters[4] = counts->p;
    parameters[5] = counts->c;

    probLatent = counts->probLatent;
    reactLatent = counts->reactLatent;
    probDefect = counts->probDefect;
    latIncompDeath = counts->latIncompDeath;
    latCompDeath = counts->latCompDeath;

    //Sampling
    sampleFileLength = counts->numSampleTimes;
  }

  /* Counts of State Variables */
  /* Initial values */
  long int numVirusInit = 0;
  long int numCellUninfectInit = 1e4 * mLBlood;
  long int numCellInfectInit = 1;
  long int numLatentCompInit = 0;
  long int numLatentIncompInit = 0;

  double prodInfectionRates[3]; /* These are used frequently, so they are only calculated once here */
  prodInfectionRates[0] = (1 - probDefect) * (1 - probLatent) * parameters[1]; // UninfectToInfect
  prodInfectionRates[1] = probDefect * probLatent * parameters[1];             // UninfectToLatIncomp
  prodInfectionRates[2] = (1 - probDefect) * probLatent * parameters[1];       // UninfectToLatComp

  /* Probabilities of each of the 11 events */
  // 0 UninfectedBirth;
  // 1 UninfectToInfect;
  // 2 UninfectDeath;
  // 3 InfectedDeath;
  // 4 VirusBirth;
  // 5 VirusDeath;
  //
  // 6 UninfectToLatIncomp;
  // 7 UninfectToLatComp;
  // 8 LatIncompDeath;
  // 9 LatCompDeath;
  // 10 LatReact;

  double rates[11]; /* Rate of each of the 11 event types */
  double probs[11]; /* Probability of each of the 11 event types */
  double totRate;   /* Total rate of events */

  /* Updated throughout the program */
  long int numVirus;
  long int numCellInfect;
  long int numCellUninfect;
  long int numLatentComp;
  long int numLatentIncomp;

  /* Sampling */
  double sampleTimes[sampleFileLength];
  int numSampleVirus[sampleFileLength];
  int numSampleLatent[sampleFileLength];
  int sampleCounter = 0;

  /* Times */
  double totTime = 0;
  unsigned long int numEvents = 0; /* Total number of events that have occured */

  /* Random number */
  double waitTime = -1;    /* Exponential waiting time */
  double unif;        /* Uniform random number */
  long int pickVirus; /* Uniform random number (discrete), used to choose which virus the event changes */

 /* If not loading from a checkpoint */
  if (loadCheckPoint == 0) {
    /* Find the sample times, check they are valid sample times*/
    readSampleTimes(sampleTimes, numSampleVirus, numSampleLatent, sampleTimesFile, sampleFileLength);
    checkValidSampleTimes(sampleTimes, numSampleVirus, numSampleLatent, sampleFileLength);

    timeLastSample = sampleTimes[sampleFileLength - 1];
    totSampleVirus = 0;
    totSampleLatent = 0;

    /* Counts the number of samples to be taken  */
    for (int i = 0; i < sampleFileLength; i++) {
      totSampleVirus = totSampleVirus + numSampleVirus[i];
      totSampleLatent = totSampleLatent + numSampleLatent[i];
    }

    /* Create arrays for the tips to be stored in */
    maxActive =  10000 * mLBlood;
    maxVirus = 1000000 * mLBlood;
    maxIncomp = 5000 * mLBlood;
    maxComp = 100 * mLBlood;

    activeArray       = malloc(sizeof(struct Node*) * maxActive);
    virusArray        = malloc(sizeof(struct Node*) * maxVirus);
    latentIncompArray = malloc(sizeof(struct Node*) * maxIncomp);
    latentCompArray   = malloc(sizeof(struct Node*) * maxComp);
    virusSampleArray  = malloc(sizeof(struct Node*) * totSampleVirus);
    latentSampleArray = malloc(sizeof(struct Node*) * totSampleLatent);

    if (activeArray == NULL || virusArray == NULL || latentIncompArray == NULL || latentCompArray == NULL || virusSampleArray == NULL || latentSampleArray == NULL) {
      fprintf(stderr, "Insufficient Memory\n");
      exit(1);
    }

    /* Sets the array elements to null */
    setVirusArrays(activeArray, virusArray, latentIncompArray, latentCompArray, virusSampleArray, latentSampleArray,
      maxActive, maxVirus, maxIncomp, maxComp, totSampleVirus, totSampleLatent);

    root_ptr = newNode(0, 0, stem); /* Root of the tree */
    activeArray[0] = root_ptr; /* Adds the initially infected cell */
    totMem = totMem + sizeof(struct Node);

  } else {
    /* Reads in the sample file */
    if (readSampleFileBinary(sampleFileLength, sampleTimes, numSampleVirus, numSampleLatent, loadName) == EXIT_FAILURE) {
        fputs("Problem loading checkpoint sample times\nExiting\n", stderr);
        exit(1);
    }

    /* Sets the counts, times, and indexes from before the checkpoint */
    numVirus = counts->numVirus;
    numCellInfect = counts->numCellInfect;
    numCellUninfect = counts->numCellUninfect;
    numLatentComp = counts->numLatentComp;
    numLatentIncomp = counts->numLatentIncomp;
    totTime  = counts->totTime;
    waitTime = counts->waitTime;

    timeLastSample = counts->timeLastSample;
    sampleCounter = counts->sampleCounter;
    numEvents = counts->numEvents;

  }
  free(sampleTimesFile);

  if (checkPointTime == -1) {
    checkPointTime = timeLastSample + 1; /* If no checkpoint time is given, there is no checkpoint  */
  }

  totMem = totMem - sizeof(char) * 100;
  /* Adds total memory allocated thus far */
  totMem = totMem + sizeof(struct Node*) * (maxActive + maxVirus + maxIncomp + maxComp + totSampleVirus + totSampleLatent);
  stem->left = root_ptr; /* Stem points to the root */

  int runFinish = 0; /* Keeps track of if the last sampling time was reached */
  unsigned long int maxMem = 0;

  int timesReset = 0; /* Used to determine whether a checkpoint was reloaded and then the virus went extinct */
  struct Node* tmp_Node; /* Used when the tree needs to be pruned */

  if (print) {
    printf("totTime, numVirus, numCellInfect, numCellUninfect, numLatentComp, numLatentIncomp, memory\n");
  }

  printf("checkpoint time %f\n", checkPointTime);
  /* Repeat the simulation numRep times */
  while (runFinish == 0){
    /* Reset the initial conditions for every simulation */

    if (loadCheckPoint == 0) {
      resetCounts(&numVirus, numVirusInit, &numCellUninfect, numCellUninfectInit,
          &numCellInfect, numCellInfectInit, &numLatentComp, numLatentCompInit,
          &numLatentIncomp, numLatentIncompInit, &numLatentSample, &numVirusSample, &sampleCounter, &totTime);
    } else {
      if(timesReset != 0 ) {
        fprintf(stderr, "Extinction in simulation loaded from checkpoint. Exiting.\n");
        exit(1);
      }
      timesReset = timesReset + 1;
    }

    setRates(prodInfectionRates, parameters, reactLatent,  latIncompDeath, latCompDeath,
      rates, numVirus, numCellInfect, numCellUninfect, numLatentComp, numLatentIncomp);

    /* Continue simulating until the last time is reached or there are no viruses/infectee cells or no cells */
    while (totTime < timeLastSample && numVirus + numCellInfect > 0 && numCellInfect + numCellUninfect > 0) {

      /* Check the arrays are not too small */
      if (maxMem < totMem) {
        maxMem = totMem;
      }

      CalcProb(rates, probs, totRate);
      assert(totMem  / 25 < 1073741824 );

      /* Prints population sizes */
      if (print && numEvents % printFreq == 0 ) {
        printf("%f, %ld, %ld, %ld, %ld, %ld, %f\n", totTime, numVirus, numCellInfect, numCellUninfect, numLatentComp, numLatentIncomp, totMem/1073741824.0);
        fflush(stdout);
      }

      /* Draw time until next event */
      if (startReload) {
        /* Do not draw a new waiting time if the first time through loop after
        loading the checkpoint */
        startReload = 0;
      } else {
        waitTime = gsl_ran_exponential (r, 1 / totRate);
      }
      if (totTime + waitTime > checkPointTime) {

        // Record cell types in tips
        findCellTypes(activeArray, virusArray, latentIncompArray, latentCompArray, virusSampleArray, latentSampleArray, numCellInfect,
          numVirus, numLatentIncomp, numLatentComp, numVirusSample, numLatentSample);

        unsigned long treeMem = findTreeStringMem(numVirus, numCellInfect, numLatentComp, numLatentIncomp, numVirusSample, numLatentSample, timeLastSample);
        char* treeString = malloc(treeMem * sizeof(char));


        // Save tree structure
        long int nodeNum = 0;

        checkPointTree(stem->left, treeString, &nodeNum);

        char treeBiFile[150] = "\0";
        strcat(treeBiFile, outfile);
        strcat(treeBiFile, "tree.chkpt");

        writeBinary(treeString, treeBiFile);

        // save latent history/states/types
        nodeNum = 0;
        unsigned long latentMem = findLatentStringMem(numVirus, numCellInfect, numLatentComp, numLatentIncomp, numVirusSample, numLatentSample, timeLastSample);

        char * latentString = malloc(latentMem * sizeof(char));
        checkPointLatent(stem->left, latentString, &nodeNum);

        nodeNum = 0;

        char latentBiFile[150] = "\0";
        strcat(latentBiFile, outfile);
        strcat(latentBiFile, "latent.chkpt");

        writeBinary(latentString, latentBiFile);

        /* Saves the sample times and number of samples to be taken */
        char sampleBiFile[150] = "\0";
        strcat(sampleBiFile, outfile);
        strcat(sampleBiFile, "sampleTimes.chkpt");

        writeSampleFileBinary(sampleFileLength, sampleTimes, numSampleVirus, numSampleLatent, sampleBiFile);

        /* Saves all the parameters, counts, indeces, etc.*/
        writeCounts(numVirus, numCellInfect, numCellUninfect, numLatentComp,
        numLatentIncomp, numVirusSample, numLatentSample, totTime, waitTime, treeString, latentString,
        mLBlood, parameters, probLatent, reactLatent, probDefect, latIncompDeath,
        latCompDeath, sampleFileLength, sampleCounter, maxActive, maxVirus, maxIncomp,
        maxComp, totSampleVirus, totSampleLatent, timeLastSample, numEvents, outfile);

        /* Saves the state of the random number generator*/
        writeGslState(r, outfile);

        // stop program/exit gracefully
        runFinish = 2; /* memory should get freed  */
        printf("%f, %ld, %ld, %ld, %ld, %ld, %f\n", totTime, numVirus, numCellInfect, numCellUninfect, numLatentComp, numLatentIncomp, totMem/1073741824.0);
        totTime = timeLastSample + 1; /* Ensures the simulation is not restarted */
        fflush(stdout);

      }
      /* Checks if passed the sample time */
       else if (totTime + waitTime > sampleTimes[sampleCounter]) {
        totTime = sampleTimes[sampleCounter];

        /* If there are not enough viruses to sample, the virus should be on its way to going extinct. Stops the program */
        if (numSampleVirus[sampleCounter] > numVirus) {
          printf("Not enough viruses to sample %f, %ld, %ld, %ld, %ld, %ld, %f\n", totTime, numVirus, numCellInfect, numCellUninfect, numLatentComp, numLatentIncomp, totRate);
          break;
        }
          sampleEvent(r, &sampleCounter, numSampleVirus, numSampleLatent,
          virusArray, virusSampleArray, latentSampleArray, latentIncompArray, latentCompArray,
          &numVirusSample, &numLatentSample, &numVirus, &numLatentIncomp, &numLatentComp, totTime);

        setRates(prodInfectionRates, parameters, reactLatent, latIncompDeath, latCompDeath,
          rates, numVirus, numCellInfect, numCellUninfect, numLatentComp, numLatentIncomp);

      } else {
        totTime = totTime + waitTime;

        /* Pick an event */
        unif = gsl_rng_uniform (r);
        ++numEvents;

        /* Uninfected birth */
        if (unif < probs[0]) {
          // printf("Uninfected birth\n");
          ++numCellUninfect;
          setInfectionRates(rates, prodInfectionRates, numCellUninfect, numVirus);
          rates[2] = rates[2] + parameters[2];

          /* Uninfected cell becomes infected */
        } else if (unif < probs[0] + probs[1]) {
          // printf("Uninfected to infected");
          pickVirus = gsl_rng_uniform_int (r, numVirus);
          birthEvent(virusArray, activeArray, pickVirus, &numCellInfect, totTime, 1, &numCellUninfect, &totMem);
          setInfectionRates(rates, prodInfectionRates, numCellUninfect, numVirus);

          rates[2] = rates[2] - parameters[2];
          rates[3] = rates[3] + parameters[3];
          rates[4] = rates[4] + parameters[4];

          if (numCellInfect == maxActive - 1) {
            activeArray = reallocArray(&maxActive, activeArray, &totMem);
          }

          // assert(numCellInfect < maxActive- 1);
          /* Uninfected cell dies */
        } else if (unif <  probs[0] + probs[1] + probs[2]) {
          // printf("Uninfected death\n");
          --numCellUninfect;
          setInfectionRates(rates, prodInfectionRates, numCellUninfect, numVirus);
          rates[2] = rates[2] - parameters[2];

        /* Infected cell dies */
      } else if (unif < probs[0] + probs[1] + probs[2] + probs[3]) {
        // printf("Infected death\n");
          if (numCellInfect - 1 + numVirus > 0) {
            pickVirus = gsl_rng_uniform_int (r, numCellInfect);
            tmp_Node = activeArray[pickVirus];
            removeVirus(activeArray, pickVirus, numCellInfect - 1, &numCellInfect);
            pruneTip(tmp_Node, &totMem);

          } else {
            --numCellInfect;
          }

          /* Decrease memory allocation to array*/
          if((numCellInfect < maxActive * 5 /6  - 500) && totTime > 50) {
            activeArray = reallocArrayDecrease(&maxActive, activeArray, &totMem);
          }

          rates[3] = rates[3] - parameters[3];
          rates[4] = rates[4] - parameters[4];

          /* Virus is born */
        } else if (unif < probs[0] + probs[1] + probs[2] + probs[3] + probs[4]) {
          // printf("virus birth\n");
          pickVirus = gsl_rng_uniform_int (r, numCellInfect);
          birthEvent(activeArray, virusArray, pickVirus, &numVirus, totTime, 0, NULL, &totMem);

          setInfectionRates(rates, prodInfectionRates, numCellUninfect, numVirus);
          rates[5] = rates[5] + parameters[5];
          if (numVirus ==  maxVirus - 1) {
            virusArray = reallocArray(&maxVirus, virusArray, &totMem);
          }
          // assert(numVirus <  maxVirus - 1);
          /* Virus dies */
        } else if (unif < probs[0] + probs[1] + probs[2] + probs[3] + probs[4] + probs[5]) {
          // printf("Virus death \n");
            if (numCellInfect - 1 + numVirus > 0) {
              pickVirus = gsl_rng_uniform_int (r, numVirus);
              tmp_Node = virusArray[pickVirus];
              removeVirus(virusArray, pickVirus, numVirus - 1, &numVirus);
              pruneTip(tmp_Node, &totMem);
              } else {
              --numVirus;
            }

            /* Decrease memory allocation to array*/
            if ((numVirus <  maxVirus * 5/6 - 500) && totTime > 50) {
              virusArray = reallocArrayDecrease(&maxVirus, virusArray, &totMem);
            }

            setInfectionRates(rates, prodInfectionRates, numCellUninfect, numVirus);
            rates[5] = rates[5] - parameters[5];

          /* Uninfected cell becomes latent replication incompetent cell*/
        } else if (unif < probs[0] + probs[1] + probs[2] + probs[3] + probs[4] + probs[5] + probs[6]) {
          // printf("Uninfected to latent incompetent\n");

            pickVirus = gsl_rng_uniform_int (r, numVirus);
            birthEvent(virusArray, latentIncompArray, pickVirus, &numLatentIncomp, totTime, 1, &numCellUninfect, &totMem);

            setInfectionRates(rates, prodInfectionRates, numCellUninfect, numVirus);
            rates[2] = rates[2] - parameters[2];
            rates[8] = rates[8] + latIncompDeath;

            if (numLatentIncomp == maxIncomp - 1) {
              latentIncompArray = reallocArray(&maxIncomp, latentIncompArray, &totMem);
            }

              // assert(numLatentIncomp < maxIncomp  - 1);
          /* Uninfected cell becomes latent replication competent cell*/
        } else if (unif < probs[0] + probs[1] + probs[2] + probs[3] + probs[4] + probs[5] + probs[6] + probs[7]) {
            // printf("Uninfected to latent competent\n");

            pickVirus = gsl_rng_uniform_int (r, numVirus);
            birthEvent(virusArray, latentCompArray, pickVirus, &numLatentComp, totTime, 1, &numCellUninfect, &totMem );

            setInfectionRates(rates, prodInfectionRates, numCellUninfect, numVirus);
            rates[2] = rates[2] - parameters[2];
            rates[9] = rates[9] + latCompDeath;
            rates[10] = rates[10] + reactLatent;

            if (numLatentComp == maxComp - 1) {
              latentCompArray = reallocArray(&maxComp, latentCompArray, &totMem);
            }

            // assert(numLatentComp < maxComp - 1);
          /* Latent replication incompetent cell death */
        } else if (unif < probs[0] + probs[1] + probs[2] + probs[3] + probs[4] + probs[5] + probs[6] + probs[7] + probs[8]) {
          // printf("Latent incompetent death\n");

            pickVirus = gsl_rng_uniform_int (r, numLatentIncomp);
            tmp_Node = latentIncompArray[pickVirus];
            removeVirus(latentIncompArray, pickVirus, numLatentIncomp - 1, &numLatentIncomp);
            pruneTip(tmp_Node, &totMem);

            /* Decrease memory allocation to array*/
            if ((numLatentIncomp < maxIncomp * 5 /6 - 500) && totTime > 50) {
              latentIncompArray = reallocArrayDecrease(&maxIncomp, latentIncompArray, &totMem);
            }

            rates[8] = rates[8] - latIncompDeath;

        /* Latent replication competent cell death */
      } else if (unif < probs[0] + probs[1] + probs[2] + probs[3] + probs[4] + probs[5] + probs[6] + probs[7] + probs[8] + probs[9] ) {
          // printf("latent competent death\n");

            pickVirus = gsl_rng_uniform_int (r, numLatentComp);
            tmp_Node = latentCompArray[pickVirus];
            removeVirus(latentCompArray, pickVirus, numLatentComp - 1, &numLatentComp);
            pruneTip(tmp_Node, &totMem);

            /* Decrease memory allocation to array*/
            if ((numLatentComp < maxComp * 5/6 - 500) && totTime > 50) {
              latentCompArray = reallocArrayDecrease(&maxComp, latentCompArray, &totMem);
            }

            rates[9]  = rates[9] - latCompDeath;
            rates[10] = rates[10] - reactLatent;

        /* Latent virus reactivates */
        } else {
          // printf("reactivate\n");
          pickVirus = gsl_rng_uniform_int (r, numLatentComp);
          reactivateLatent(latentCompArray, activeArray, pickVirus, &numLatentComp, &numCellInfect, totTime);

          rates[3] = rates[3] + parameters[3];
          rates[4] = rates[4] + parameters[4];
          rates[9] = rates[9] - latCompDeath;
          rates[10] = rates[10] - reactLatent;

          if (numCellInfect == maxActive - 1) {
            activeArray = reallocArray(&maxActive, activeArray, &totMem);
          }

          /* Decrease memory allocation to array */
          if ((numLatentComp < maxComp * 5/6 - 500) && totTime > 50) {
            latentCompArray = reallocArrayDecrease(&maxComp, latentCompArray, &totMem);
          }
          // assert(numCellInfect < maxActive - 1);
        }
      }
    }

    /* If the simulation didn't finish, start over */
   if (totTime < timeLastSample) {
      printf("Did not finish %ld %f %ld %ld %ld %ld %ld\n", numEvents, totTime, numCellUninfect, numCellInfect, numVirus, numLatentComp, numLatentIncomp);
      if (restart) {
        printf("Starting again\n");

        root_ptr = reset(stem, root_ptr, &sampleCounter, &numEvents, &totMem,
        activeArray, virusArray, latentIncompArray, latentCompArray, virusSampleArray, latentSampleArray,
        maxActive, maxVirus, maxIncomp, maxComp, totSampleVirus, totSampleLatent); //not going back to starting size-- if you got that large before you will probably get that large again for arrays

      } else {
        runFinish = -1;
      }

    } else {
      if (runFinish != 2) {
        runFinish = 1;
      }
    }
  }

  if (runFinish == 1) {

    findSampledAll(virusSampleArray, latentSampleArray, totSampleVirus, totSampleLatent, stem);
    long int nodeName = 0;
    char *sampleTreeString = strdup("");
    sampleTreeString = NewickSampled(stem->left, sampleTreeString, 0, &nodeName);

    char treeFile[150] = "\0";
    strcat(treeFile, outfile);
    strcat(treeFile, "SampledTree.txt");
    writeTxt(treeFile, sampleTreeString);

    nodeName = 0;
    char* latentString = strdup("");
    latentString = LatentSampled(stem->left, latentString, 0, &nodeName);
    char latentFileOut[150] = "\0";
    strcat(latentFileOut, outfile);
    strcat(latentFileOut, "LatentStates.txt");
    writeTxt(latentFileOut, latentString);

    printf("Max memory used %f\n", maxMem / 1073741824.0 );

  }
  free(controlFile);
  free(counts);
  regfree(&regDecimal);
  regfree(&regInt);
  ClearMemory(stem->left, stem, activeArray, latentCompArray, latentIncompArray, virusArray, virusSampleArray, latentSampleArray, &totMem);

  gsl_rng_free (r);
  return(0);
}
