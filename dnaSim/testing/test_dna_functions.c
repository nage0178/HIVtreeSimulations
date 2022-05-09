# include "../src/dna_sim.h"
# include "../src/string_utilities.h"
# include <math.h>

gsl_rng * r;  /* global generator */

// Checks the FloatEquals functions
void testFloatEquals() {
  assert(FloatEquals(1.1, 1.1001, 1e-4));
  assert(!(FloatEquals(1.1, 1.01, 1e-4)));
  assert(FloatEquals(2.0, 2.0 + 1e-15, 1e-14));
  assert(FloatEquals(-2.3, -2.30006, 1e-3));
  assert(!FloatEquals(-2.3, -2.1, 1e-3));
  return;
}

void checkMakeTree() {
  char *treeString = string_from_file("testing/testingSampledTree.txt");

  /*Make a tree */
  struct NodeSampledTree* root_ptr = newNodeSampledTree(0, 1);
  makeTree(root_ptr, treeString, 0);
  printTree(root_ptr, "testTreeRead.txt");

  char *newTreeString = string_from_file("testTreeRead.txt");
  assert(strcmp(treeString, newTreeString) == 0);

  int out = system("rm testTreeRead.txt");

  if (out == -1) {
    fprintf(stderr, "Problem in checkMakeTree. \n");
    exit(1);
  }

  free(treeString);
  free(newTreeString);
  clearTree(root_ptr);
}

void checkCountTips() {
  char *treeString = string_from_file("testing/testingSampledTree.txt");

  /*Make a tree */
  struct NodeSampledTree* root_ptr = newNodeSampledTree(0, 1);
  makeTree(root_ptr, treeString, 0);
  int numTips = countTips(root_ptr, 0);
  assert(numTips == 11);

  free(treeString);
  clearTree(root_ptr);
}

char*  printLatent(struct NodeSampledTree* node, char* latentString ) {
  Sasprintf(latentString, "%s%d:%f,%f,%d\n", latentString, node->name, node->totTimeLatent, node->sample_time, node->isLatent); /* -1 indicates the internal node is not sampled */
  if (node->left) {
    latentString = printLatent(node->left, latentString);
    latentString = printLatent(node->right, latentString);
  }
  return latentString;
}

void checkReadLatent() {
  char *treeString = string_from_file("testing/testSampledTree.txt");

  /*Make a tree */
  /* Note that  birth times are actually branch lengths-figure out if they are where you want them to be */
  struct NodeSampledTree* root_ptr = newNodeSampledTree(0, 1);
  int endChar = makeTree(root_ptr, treeString, 0);
  if (treeString[endChar] != '\0') {
    fprintf(stderr, "Problem reading in tree. Did not reach the end of the file \n");
    exit(1);
  }

  FILE *latentFile;
  latentFile = fopen("testing/testLatentStates.txt", "r");
  readLatent(root_ptr, latentFile);
  fclose(latentFile);
  free(treeString);

  char* latentString = strdup("");
  latentString = printLatent(root_ptr, latentString);

  FILE* file = fopen("checkLatentStates.txt", "w");
  fputs(latentString, file);
  fclose(file);
  free(latentString);

  int one = system("diff testing/testLatentStates.txt checkLatentStates.txt > compareLatent");
  char* compare = string_from_file("compareLatent");
  assert(strlen(compare) == 0);
  int two = system("rm checkLatentStates.txt  compareLatent");

  if (one == -1 || two == -1) {
    fprintf(stderr, "Problem checking Latent file. \n");
    exit(1);
  }
  free(compare);
  clearTree(root_ptr);
}

void checkSimulateDNABranch(const gsl_rng * r) {
  int base = 0;
  double statFreq[4] = {.25, .25, .25, .25};
  double rateParams[6] = {1, 1, 1, 1, 1, 1};
  double instRate[4][4] = { {                          0, rateParams[0] * statFreq[1],    rateParams[1] * statFreq[2],  rateParams[2] * statFreq[3]},
                {rateParams[0] * statFreq[0],                           0,    rateParams[3] * statFreq[2],  rateParams[4] * statFreq[3]},
                {rateParams[1] * statFreq[0], rateParams[3] * statFreq[1],                              0,  rateParams[5] * statFreq[3]},
                {rateParams[2] * statFreq[0], rateParams[4] * statFreq[1],    rateParams[5] * statFreq[2],                           0}};
  double mu = .5;
  makeInstantaneousRate(instRate, statFreq, mu);

  double branchLength = 1.5;
  int returnBase;
  double baseFreq[] = {0, 0, 0, 0};
  int numRep = 10000000;

  for (int i = 0; i < numRep; i ++) {
    returnBase = simulateDNABranch(base, instRate, branchLength, r);

    if (returnBase == 0) {
      baseFreq[0]++;
    } else if (returnBase == 1) {
      baseFreq[1]++;
    } else if (returnBase == 2) {
      baseFreq[2]++;
    } else {
      baseFreq[3]++;
    }
  }

  for (int i = 0; i < 4; i++) {
    baseFreq[i] = baseFreq[i] /  numRep;
  }

  double probSame = 1/4.0 + 3/4.0 * exp(-4 * (mu/3) * branchLength);
  double probOther = 1/4.0 - 1/4.0 * exp (-4 * (mu/3) * branchLength);
  //printf("%f %f %f %f %f %f \n", baseFreq[0], baseFreq[1], baseFreq[2], baseFreq[3], probSame, probOther);

  assert(FloatEquals(baseFreq[0], probSame, 0.005));
  assert(FloatEquals(baseFreq[1], probOther, 0.005));
  assert(FloatEquals(baseFreq[2], probOther, 0.005));
  assert(FloatEquals(baseFreq[3], probOther, 0.005));

  double a = 1.0, b = 3.5;
  double k = a / b;
  instRate[0][1] = a;
  instRate[0][2] = b;
  instRate[0][3] = b;

  instRate[1][0] = a;
  instRate[1][2] = b;
  instRate[1][3] = b;

  instRate[2][0] = b;
  instRate[2][1] = b;
  instRate[2][3] = a;

  instRate[3][0] = b;
  instRate[3][1] = b;
  instRate[3][2] = a;

  mu = 1.2;

  makeInstantaneousRate(instRate, statFreq, mu);

  for (int i = 0; i < 4; i++) {
    baseFreq[i] = 0;
  }

  base = 1;
  branchLength = 2.3;

  for (int i = 0; i < numRep; i ++) {
    returnBase = simulateDNABranch(base, instRate, branchLength, r);

    if (returnBase == 0) {
      baseFreq[0]++;
    } else if (returnBase == 1) {
      baseFreq[1]++;
    } else if (returnBase == 2) {
      baseFreq[2]++;
    } else {
      baseFreq[3]++;
    }
  }

  for (int i = 0; i < 4; i++) {
    baseFreq[i] = baseFreq[i] /  numRep;
  }

  double d = mu * branchLength;
  double probZero = 1/4.0 + 1/4.0 * exp(- 4.0 * d / (k + 2)) + 1/2.0  * exp(-2 * d * (k + 1) / (k + 2));
  double probOne = 1/4.0  + 1/4.0  * exp(- 4.0 *  d / (k + 2)) - 1/2.0  * exp(-2 * d * (k + 1) / (k + 2));
  double probTwo = 1/4.0  - 1/4.0  * exp(- 4.0 * d / (k + 2));

  assert(FloatEquals(baseFreq[0], probOne, 0.005));
  assert(FloatEquals(baseFreq[1], probZero, 0.005));
  assert(FloatEquals(baseFreq[2], probTwo, 0.005));
  assert(FloatEquals(baseFreq[3], probTwo, 0.005));

  // TN93 model
  double a1 = 1.2, a2 = .75;
  b = .3;

  for (int i = 0; i < 4; i ++) {
    statFreq[i] = (i + 1)  * 0.1;
  }

  instRate[0][1] = b  * statFreq[1];
  instRate[0][2] = a2 * statFreq[2];
  instRate[0][3] = b  * statFreq[3];

  instRate[1][0] = b  * statFreq[0];
  instRate[1][2] = b  * statFreq[2];
  instRate[1][3] = a1 * statFreq[3];

  instRate[2][0] = a2 * statFreq[0];
  instRate[2][1] = b  * statFreq[1];
  instRate[2][3] = b  * statFreq[3];

  instRate[3][0] = b  * statFreq[0];
  instRate[3][1] = a1 * statFreq[1];
  instRate[3][2] = b  * statFreq[2];


  base = 2;
  branchLength = 1.9;
  makeInstantaneousRate(instRate, statFreq, mu);

  for (int i = 0; i < numRep; i ++) {
    returnBase = simulateDNABranch(base, instRate, branchLength, r);

    if (returnBase == 0) {
      baseFreq[0]++;
    } else if (returnBase == 1) {
      baseFreq[1]++;
    } else if (returnBase == 2) {
      baseFreq[2]++;
    } else {
      baseFreq[3]++;
    }
  }

  for (int i = 0; i < 4; i++) {
    baseFreq[i] = baseFreq[i] /  numRep;
  }

  double probGA, probGC, probGG, probGT;
  double e2, e3, piR, piY;
  piR = statFreq[0] + statFreq[2];
  piY = statFreq[1] + statFreq[3];

  a2 = instRate[0][2] / statFreq[2];
  b = instRate[0][1] / statFreq[1];

  e2 = exp(- b * branchLength);
  e3 = exp(- branchLength * (piR * a2 +  piY * b));
  probGA = statFreq[0] + statFreq[0] * piY / piR * e2 - statFreq[0] / piR * e3;
  probGC = statFreq[1] * (1 - e2);
  probGG = statFreq[2] + statFreq[2] * piY / piR * e2 + statFreq[0] / piR * e3;
  probGT = statFreq[3] * (1 - e2);

  assert(FloatEquals(baseFreq[0], probGA, 0.005));
  assert(FloatEquals(baseFreq[1], probGC, 0.005));
  assert(FloatEquals(baseFreq[2], probGG, 0.005));
  assert(FloatEquals(baseFreq[3], probGT, 0.005));
  //printf("%f %f %f %f \n%f %f %f %f\n", baseFreq[0], baseFreq[1], baseFreq[2], baseFreq[3], probGA, probGC, probGG, probGT);

}

void k80prob (double d, double k, double *zero, double* one, double* two) {
  *zero  = 1/4.0 + 1/4.0 * exp(- 4.0 * d / (k + 2)) + 1/2.0  * exp(-2 * d * (k + 1) / (k + 2));
  *one = 1/4.0  + 1/4.0  * exp(- 4.0 * d / (k + 2)) - 1/2.0  * exp(-2 * d * (k + 1) / (k + 2));
  *two = 1/4.0  - 1/4.0  * exp(- 4.0 * d / (k + 2));
}

void probTip (double probTrans[4], int tipBase, double probZero, double probOne, double probTwo) {

  int transitionMat[4][4] = {{0, 1, 2, 2},
                             {1, 0, 2, 2},
                             {2, 2, 0, 1},
                             {2, 2, 1, 0}};

  for (int i = 0; i < 4; i ++) {
    if (transitionMat[i][tipBase] == 0) {
      probTrans[i] = probZero;

    } else if (transitionMat[i][tipBase] == 1) {
      probTrans[i] = probOne;

    } else if (transitionMat[i][tipBase] == 2) {
      probTrans[i] = probTwo;
    }
  }

  return;
}

void probInternal (double probToRoot[4], double probToTips[4], double probZero, double probOne, double probTwo) {

int transitionMat[4][4] = {{0, 1, 2, 2},
                           {1, 0, 2, 2},
                           {2, 2, 0, 1},
                           {2, 2, 1, 0}};
double probTrans;
  for (int i = 0; i < 4; i ++) {
    probToRoot[i] = 0;

    for (int tipBase = 0; tipBase < 4; tipBase ++) {
      if (transitionMat[i][tipBase] == 0) {
        probTrans = probZero;

      } else if (transitionMat[i][tipBase] == 1) {
        probTrans = probOne;

      } else if (transitionMat[i][tipBase] == 2) {
        probTrans = probTwo;
      }

      probToRoot[i] = probToRoot[i] + probToTips[tipBase] * probTrans;
    }
  }

}

double poorlyWritingCheckSim (struct NodeSampledTree* root_ptr, double mu, double k, int site1, int site2, int site3, int site4) {
  double d = mu * (root_ptr->right->right->left->branchLength);
  double probZero, probOne, probTwo;
  k80prob(d, k, &probZero, &probOne, &probTwo);

  double d2 = mu * (root_ptr->right->right->right->branchLength - root_ptr->right->right->right->totTimeLatent);
  double probZero2, probOne2, probTwo2;
  k80prob(d2, k, &probZero2, &probOne2, &probTwo2);

  double probGTg[4];
  double probGg[4], probTg[4];
  probTip(probGg, site3, probZero, probOne, probTwo);
  probTip(probTg, site4, probZero2, probOne2, probTwo2);

  for (int i = 0; i < 4; i++) {
    probGTg[i] = probGg[i] * probTg[i];
  }

  d = mu * (root_ptr->right->left->branchLength);
  k80prob(d, k, &probZero, &probOne, &probTwo);
  d2 = mu * (root_ptr->right->right->branchLength);
  k80prob(d2, k, &probZero2, &probOne2, &probTwo2);

  double probCGTg[4];
  double probCg[4];
  probTip(probCg, site2, probZero, probOne, probTwo);
  double probGTbranch[4];
  probInternal (probGTbranch, probGTg, probZero2, probOne2, probTwo2);

  for (int i = 0; i < 4; i++) {
    probCGTg[i] = probGTbranch[i] * probCg[i];
  }

  d = mu * (root_ptr->left->branchLength);
  k80prob(d, k, &probZero, &probOne, &probTwo);
  d2 = mu * (root_ptr->right->branchLength);
  k80prob(d2, k, &probZero2, &probOne2, &probTwo2);

  probTip(probTg, site1, probZero, probOne, probTwo);
  double probCGTbranch[4];
  probInternal(probCGTbranch, probCGTg, probZero2, probOne2, probTwo2);
  double probFinal[4];
  for(int i = 0; i < 4; i++) {
    probFinal[i] = probTg[i]*probCGTbranch[i];
  }

return(probFinal[0]);
}

double poorlyWritingCheckSimBalanced (struct NodeSampledTree* root_ptr, double mu, double k, int site1, int site2, int site3, int site4) {
  double d = mu * (root_ptr->left->left->branchLength);
  double probZero, probOne, probTwo;
  k80prob(d, k, &probZero, &probOne, &probTwo);

  double d2 = mu * (root_ptr->left->right->branchLength);
  double probZero2, probOne2, probTwo2;
  k80prob(d2, k, &probZero2, &probOne2, &probTwo2);

  double probLeftg[4];
  double probSite1g[4], probSite2g[4];
  probTip(probSite1g, site1, probZero, probOne, probTwo);
  probTip(probSite2g, site2, probZero2, probOne2, probTwo2);

  for (int i = 0; i < 4; i++) {
    probLeftg[i] = probSite1g[i] * probSite2g[i];
  }

  d = mu * (root_ptr->right->left->branchLength);
  k80prob(d, k, &probZero, &probOne, &probTwo);

  d2 = mu * (root_ptr->right->right->branchLength);
  k80prob(d2, k, &probZero2, &probOne2, &probTwo2);

  double probRightg[4];
  double probSite3g[4], probSite4g[4];
  probTip(probSite3g, site3, probZero, probOne, probTwo);
  probTip(probSite4g, site4, probZero2, probOne2, probTwo2);

  for (int i = 0; i < 4; i++) {
    probRightg[i] = probSite3g[i] * probSite4g[i];
  }

  d = mu * (root_ptr->left->branchLength);
  k80prob(d, k, &probZero, &probOne, &probTwo);
  d2 = mu * (root_ptr->right->branchLength - root_ptr->right->totTimeLatent);
  k80prob(d2, k, &probZero2, &probOne2, &probTwo2);

  double probrootLeftg[4];
  probInternal (probrootLeftg, probLeftg, probZero, probOne, probTwo);
  double probrootRightg[4];
  probInternal (probrootRightg, probRightg, probZero2, probOne2, probTwo2);

  double probRoot[4];
  for (int i = 0; i < 4; i++) {
    probRoot[i] = probrootLeftg[i] * probrootRightg[i];
  }

return(probRoot[2]);
}


void checksimulateDNATree(const gsl_rng * r) {
  struct NodeSampledTree* root_ptr = newNodeSampledTree(0.0, 1);
  root_ptr->totTimeLatent = 0;
  root_ptr->left = newNodeSampledTree(4.5, 2);
  root_ptr->left->totTimeLatent = 0;
  root_ptr->right = newNodeSampledTree(.3, 3);
  root_ptr->right->totTimeLatent = 0;

  root_ptr->right->left = newNodeSampledTree(4.2, 4);
  root_ptr->right->left->totTimeLatent = 0;
  root_ptr->right->right = newNodeSampledTree(3.5, 5);
  root_ptr->right->right->totTimeLatent = 0;

  root_ptr->right->right->left = newNodeSampledTree(.7, 6);
  root_ptr->right->right->left->totTimeLatent = 0;
  root_ptr->right->right->right = newNodeSampledTree(.7, 7);
  root_ptr->right->right->right->totTimeLatent = .3;
  root_ptr->right->right->right->isLatent = 1;


  int numTips = 4;
  int numBases = 1;
  double a = 1.5;
  double b = 1.2;
  double k = a/b;

  double statFreq[4] = {.25, .25, .25, .25};
  double instRate[4][4] = {{0, a, b, b},
                           {a, 0, b, b},
                           {b, b, 0, a},
                           {b, b, a, 0}};
  double mu = .2;
  makeInstantaneousRate(instRate, statFreq, mu);

  int base = 0;
  int* nodeNames         = malloc(sizeof(int) * numTips);
  int* nodeLatentstate   = malloc(sizeof(int) * numTips);
  double* nodeSampleTime = malloc(sizeof(double) * numTips);
  double* nodeLatentTime = malloc(sizeof(double) * numTips);
  char** alignment = allocateAlignmentMem(numTips, numBases);
  double probAAAC = 0;
  double probCAAA = 0;
  double probAAAA = 0;
  long int numRep = 1000000;
  int baseIndex = 0;
  double gammaDraw = .65;

  for (long int i= 0; i < numRep; i++) {
    simulateDnaTree(root_ptr, base, instRate, r, alignment, baseIndex, 0, nodeNames, nodeLatentstate, nodeSampleTime, nodeLatentTime, gammaDraw);
     // printf("%c %c %c %c\n", alignment[0][0], alignment[1][0], alignment[2][0], alignment[3][0]);
    if (alignment[0][0] == 'A' && alignment[1][0] == 'A' && alignment[2][0] == 'A' && alignment[3][0] == 'C') {
      probAAAC++;
    }

    if (alignment[0][0] == 'C' && alignment[1][0] == 'A' && alignment[2][0] == 'A' && alignment[3][0] == 'A') {
      probCAAA++;
    }
    if (alignment[0][0] == 'A' && alignment[1][0] == 'A' && alignment[2][0] == 'A' && alignment[3][0] == 'A') {
      probAAAA++;
    }


  }
  probAAAC = probAAAC / numRep;
  probCAAA = probCAAA / numRep;
  probAAAA = probAAAA / numRep;

  assert(FloatEquals(probAAAA, poorlyWritingCheckSim (root_ptr, mu * gammaDraw,  k, 0, 0, 0, 0), .0006));
  assert(FloatEquals(probCAAA, poorlyWritingCheckSim (root_ptr, mu * gammaDraw,  k, 1, 0, 0, 0), .0005));
  assert(FloatEquals(probAAAC, poorlyWritingCheckSim (root_ptr, mu * gammaDraw,  k, 0, 0, 0, 1), .0005));

  assert(nodeNames[0] == 2);
  assert(nodeNames[1] == 4);
  assert(nodeNames[2] == 6);
  assert(nodeNames[3] == 7);

  assert(nodeLatentstate[0] == 0);
  assert(nodeLatentstate[1] == 0);
  assert(nodeLatentstate[2] == 0);
  assert(nodeLatentstate[3] == 1);

  clearTree(root_ptr);

  struct NodeSampledTree* root_ptr2 = newNodeSampledTree(0.0, 1);
  root_ptr2->totTimeLatent = 0;
  root_ptr2->left = newNodeSampledTree(2, 2);
  root_ptr2->left->totTimeLatent = 0;
  root_ptr2->right = newNodeSampledTree(2.5, 3);
  root_ptr2->right->totTimeLatent = 1.5;

  root_ptr2->left->left = newNodeSampledTree(1, 4);
  root_ptr2->left->left->totTimeLatent = 0;
  root_ptr2->left->right = newNodeSampledTree(1, 5);
  root_ptr2->left->right->totTimeLatent = 0;

  root_ptr2->right->left = newNodeSampledTree(.5, 6);
  root_ptr2->right->left->totTimeLatent = 0;
  root_ptr2->right->right = newNodeSampledTree(.5, 7);
  root_ptr2->right->right->totTimeLatent = 0;

  base = 2;
  double probGGGG= 0;
  gammaDraw = 1.2;

  for (long int i= 0; i < numRep; i++) {
    simulateDnaTree(root_ptr2, base, instRate, r, alignment, baseIndex, 0, nodeNames, nodeLatentstate, nodeSampleTime, nodeLatentTime, gammaDraw);
     // printf("%c %c %c %c\n", alignment[0][0], alignment[1][0], alignment[2][0], alignment[3][0]);
    if (alignment[0][0] == 'G' && alignment[1][0] == 'G' && alignment[2][0] == 'G' && alignment[3][0] == 'G') {
      probGGGG++;
    }
  }
  probGGGG = probGGGG / numRep;

  assert(FloatEquals(probGGGG, poorlyWritingCheckSimBalanced (root_ptr2, mu * gammaDraw, k, 2, 2, 2, 2), .0005));
  free(nodeNames);
  free(nodeLatentstate);
  free(nodeSampleTime);
  for (int i = 0; i < numTips; i++ ) {
    free(alignment[i]);
  }
  free(alignment);
  clearTree(root_ptr2);
}

void checkMakeInstantaneousRate(){
  double instRate[4][4] = {{1, 2, 3, 4},
                    {1.1, 2.2, 3.3, 4.4},
                    {1.1, 2.1, 3.1, 4.1},
                    {4, 3, 2, 2}};
  double statFreq[4] = {.2, .3, .35, .15};
  double mu = .3;
  makeInstantaneousRate(instRate, statFreq, mu);
  assert(FloatEquals(instRate[0][0], - (instRate[0][1] + instRate[0][2] + instRate[0][3]), 1e-14));
  assert(FloatEquals(instRate[1][1], - (instRate[1][0] + instRate[1][2] + instRate[1][3]), 1e-14));
  assert(FloatEquals(instRate[2][2], - (instRate[2][0] + instRate[2][1] + instRate[2][3]), 1e-14));
  assert(FloatEquals(instRate[3][3], - (instRate[3][0] + instRate[3][1] + instRate[3][2]), 1e-14));

  assert(FloatEquals(instRate[0][1], 2 / 8.345 * mu, 1e-14));
  assert(FloatEquals(instRate[0][2], 3 / 8.345 * mu, 1e-14));
  assert(FloatEquals(instRate[0][3], 4 / 8.345 * mu, 1e-14));

  assert(FloatEquals(instRate[1][0], 1.1 / 8.345 * mu, 1e-14));
  assert(FloatEquals(instRate[1][2], 3.3 / 8.345 * mu, 1e-14));
  assert(FloatEquals(instRate[1][3], 4.4 / 8.345 * mu, 1e-14));

  assert(FloatEquals(instRate[2][0], 1.1 / 8.345 * mu, 1e-14));
  assert(FloatEquals(instRate[2][1], 2.1 / 8.345 * mu, 1e-14));
  assert(FloatEquals(instRate[2][3], 4.1 / 8.345 * mu, 1e-14));

  assert(FloatEquals(instRate[3][0], 4 / 8.345 * mu, 1e-14));
  assert(FloatEquals(instRate[3][1], 3 / 8.345 * mu, 1e-14));
  assert(FloatEquals(instRate[3][2], 2 / 8.345 * mu, 1e-14));
}

void checkNewNodeSampleTree() {
  struct NodeSampledTree* root = newNodeSampledTree(5, 1);
  assert(root->branchLength == 5);
  assert(root->sample_time == 0);
  assert(root->name == 1);
  assert(root->isLatent == 0);
  assert(root->totTimeLatent == -1);
  assert(root->left == NULL);
  assert(root->right == NULL);

  root->right = newNodeSampledTree(10, 2);
  assert(root->right->branchLength == 10);
  assert(root->right->sample_time == 0);
  assert(root->right->name == 2);
  assert(root->right->isLatent == 0);
  assert(root->right->totTimeLatent == -1);
  assert(root->right->left == NULL);
  assert(root->right->right == NULL);

  free(root->right);
  free(root);
}

void checkReadFasta() {
  int numBases = 0;
  char *fastaSeq;
  fastaSeq = string_from_file("testing/HIV1_ref_env.fasta");
  int* startSeq = readFasta(fastaSeq, &numBases);
  assert(numBases == 2571);

  FILE *file;
  file = fopen("testReadAlignment.fa", "w");
  for (int i = 0; i < numBases; i++) {
    fprintf(file, "%d", startSeq[i]);
  }

  fclose(file);


  int one = system("sed 's/0/A/g' testReadAlignment.fa |sed 's/1/C/g'|sed 's/2/G/g' | sed 's/3/T/g' >compareFasta1");
  int two = system("grep -v \">\" testing/HIV1_ref_env.fasta |tr -d \" \t\n\r\" >compareFasta2");
  int three = system("diff compareFasta1 compareFasta2 > compareFasta");
  char* compare = string_from_file("compareFasta");
  assert(strlen(compare) == 0);
  int four = system("rm compareFasta compareFasta1 compareFasta2 testReadAlignment.fa");

  if (one == -1 || two == -1 || three == -1 || four == -1) {
    fprintf(stderr, "Problem checking fasta file. \n");
    exit(1);
  }
  free(fastaSeq);
  free(compare);
  free(startSeq);
  return;
}

void checkSetStartSeq (gsl_rng *r) {
  int numBases = 100000;
  int randomStart = 1;
  double statFreq[] = {.1, .2, .3, .4};
  char* fastaSeq = NULL;
  int* startSeq = setStartSeq (&numBases, randomStart, statFreq,  fastaSeq, r);

  double numA = 0;
  double numC = 0;
  double numG = 0;
  double numT = 0;

  for (int i = 0; i < numBases; i++) {
    if (startSeq[i] == 0) {
      numA++;
    } else if (startSeq[i] == 1) {
      numC++;
    } else if (startSeq[i] == 2) {
      numG++;
    } else if (startSeq[i] == 3) {
      numT++;
    } else {
      fprintf(stderr, "Issue with creating the starting sequence. Exiting.\n");
      exit(1);
    }
  }
  numA = numA / numBases;
  numC = numC / numBases;
  numG = numG / numBases;
  numT = numT / numBases;
  assert(FloatEquals(numA, statFreq[0], .005));
  assert(FloatEquals(numC, statFreq[1], .005));
  assert(FloatEquals(numG, statFreq[2], .005));
  assert(FloatEquals(numT, statFreq[3], .005));
  free(startSeq);
}

void checkMakeAlignment() {
  struct NodeSampledTree* root_ptr = newNodeSampledTree(0.0, 1);
  root_ptr->totTimeLatent = 0;
  root_ptr->left = newNodeSampledTree(4.5, 2);
  root_ptr->left->totTimeLatent = 0;
  root_ptr->right = newNodeSampledTree(.3, 3);
  root_ptr->right->totTimeLatent = 0;

  root_ptr->right->left = newNodeSampledTree(4.2, 4);
  root_ptr->right->left->totTimeLatent = 0;
  root_ptr->right->right = newNodeSampledTree(3.5, 5);
  root_ptr->right->right->totTimeLatent = 0;

  root_ptr->right->right->left = newNodeSampledTree(.7, 6);
  root_ptr->right->right->left->totTimeLatent = 0;
  root_ptr->right->right->right = newNodeSampledTree(.7, 7);
  root_ptr->right->right->right->totTimeLatent = .3;
  root_ptr->right->right->right->isLatent = 1;

  int numTips = 4;
  int numBases = 1000000;
  double a = 1.5;
  double b = 1.2;
  double k = a/b;

  double statFreq[4] = {.25, .25, .25, .25};
  double instRate[4][4] = {{0, a, b, b},
                           {a, 0, b, b},
                           {b, b, 0, a},
                           {b, b, a, 0}};
  double mu = .2;
  makeInstantaneousRate(instRate, statFreq, mu);


  int* nodeNames         = malloc(sizeof(int) * numTips);
  int* nodeLatentstate   = malloc(sizeof(int) * numTips);
  double* nodeSampleTime = malloc(sizeof(double) * numTips);
  double* nodeLatentTime = malloc(sizeof(double) * numTips);
  char** alignment = allocateAlignmentMem(numTips, numBases);


  int* startSeq = malloc(sizeof(int) * numBases);
  for(int i = 0; i < numBases; i++) {
    startSeq[i] = 0;
  }
  makeAlignment(root_ptr, startSeq, instRate, alignment, nodeNames, nodeLatentstate, nodeSampleTime, nodeLatentTime, numBases, numTips, r, 1, 0);

  double freqBases[4][4] = {{0, 0, 0, 0},
                          {0, 0, 0, 0},
                          {0, 0, 0, 0},
                          {0, 0, 0, 0}};
  for(int j = 0; j < 4; j++) {
    for(int i = 0; i < numBases; i++) {
      if (alignment[j][i] == 'A') {
        freqBases[j][0]++;
      } else if (alignment[j][i] == 'C') {
        freqBases[j][1]++;
      } else if (alignment[j][i] == 'G') {
        freqBases[j][2]++;
      } else if (alignment[j][i] == 'T') {
        freqBases[j][3]++;
      }

    }
  }

  for(int j = 0; j < 4; j++) {
    for(int i = 0; i < 4; i++) {
      freqBases[j][i] = freqBases[j][i] / numBases;
    }
  }

  double timeTip1 = root_ptr->branchLength + root_ptr->left->branchLength;
  //Tip2 and tip3 have the same length as tip1
  double timeTip4 = root_ptr->branchLength + root_ptr->right->branchLength + root_ptr->right->right->branchLength + root_ptr->right->right->right->branchLength - root_ptr->right->right->right->totTimeLatent;

  double d;
  d = mu * timeTip1;
  double probZero = 1/4.0 + 1/4.0 * exp(- 4.0 * d / (k + 2)) + 1/2.0  * exp(-2 * d * (k + 1) / (k + 2));
  double probOne = 1/4.0  + 1/4.0  * exp(- 4.0 *  d / (k + 2)) - 1/2.0  * exp(-2 * d * (k + 1) / (k + 2));
  double probTwo = 1/4.0  - 1/4.0  * exp(- 4.0 * d / (k + 2));

  for (int i = 0; i < 3; i ++) {
    assert(FloatEquals(probZero, freqBases[i][0], .005));
    assert(FloatEquals(probOne, freqBases[i][1], .005));
    assert(FloatEquals(probTwo, freqBases[i][2], .005));
    assert(FloatEquals(probTwo, freqBases[i][3], .005));
  }


  d = mu * timeTip4;
  probZero = 1/4.0 + 1/4.0 * exp(- 4.0 * d / (k + 2)) + 1/2.0  * exp(-2 * d * (k + 1) / (k + 2));
  probOne = 1/4.0  + 1/4.0  * exp(- 4.0 *  d / (k + 2)) - 1/2.0  * exp(-2 * d * (k + 1) / (k + 2));
  probTwo = 1/4.0  - 1/4.0  * exp(- 4.0 * d / (k + 2));

  assert(FloatEquals(probZero, freqBases[3][0], .005));
  assert(FloatEquals(probOne, freqBases[3][1], .005));
  assert(FloatEquals(probTwo, freqBases[3][2], .005));
  assert(FloatEquals(probTwo, freqBases[3][3], .005));

  // All starting with G

  for(int i = 0; i < numBases; i++) {
    startSeq[i] = 2;
  }
  makeAlignment(root_ptr, startSeq, instRate, alignment, nodeNames, nodeLatentstate, nodeSampleTime, nodeLatentTime, numBases, numTips, r, 1, 0);

  for(int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j ++) {
      freqBases[i][j] = 0;
    }
  }

  for(int j = 0; j < 4; j++) {
    for(int i = 0; i < numBases; i++) {
      if (alignment[j][i] == 'A') {
        freqBases[j][0]++;
      } else if (alignment[j][i] == 'C') {
        freqBases[j][1]++;
      } else if (alignment[j][i] == 'G') {
        freqBases[j][2]++;
      } else if (alignment[j][i] == 'T') {
        freqBases[j][3]++;
      }

    }
  }

  for(int j = 0; j < 4; j++) {
    for(int i = 0; i < 4; i++) {
      freqBases[j][i] = freqBases[j][i] / numBases;
    }
  }


  d = mu * timeTip1;
  probZero = 1/4.0 + 1/4.0 * exp(- 4.0 * d / (k + 2)) + 1/2.0  * exp(-2 * d * (k + 1) / (k + 2));
  probOne = 1/4.0  + 1/4.0  * exp(- 4.0 *  d / (k + 2)) - 1/2.0  * exp(-2 * d * (k + 1) / (k + 2));
  probTwo = 1/4.0  - 1/4.0  * exp(- 4.0 * d / (k + 2));

  for (int i = 0; i < 3; i ++) {
    assert(FloatEquals(probZero, freqBases[i][2], .005));
    assert(FloatEquals(probOne, freqBases[i][3], .005));
    assert(FloatEquals(probTwo, freqBases[i][0], .005));
    assert(FloatEquals(probTwo, freqBases[i][1], .005));
  }


  d = mu * timeTip4;
  probZero = 1/4.0 + 1/4.0 * exp(- 4.0 * d / (k + 2)) + 1/2.0  * exp(-2 * d * (k + 1) / (k + 2));
  probOne = 1/4.0  + 1/4.0  * exp(- 4.0 *  d / (k + 2)) - 1/2.0  * exp(-2 * d * (k + 1) / (k + 2));
  probTwo = 1/4.0  - 1/4.0  * exp(- 4.0 * d / (k + 2));

  assert(FloatEquals(probZero, freqBases[3][2], .005));
  assert(FloatEquals(probOne, freqBases[3][3], .005));
  assert(FloatEquals(probTwo, freqBases[3][0], .005));
  assert(FloatEquals(probTwo, freqBases[3][1], .005));

  free(nodeNames);
  free(nodeLatentstate);
  free(nodeSampleTime);
  for (int i = 0; i < numTips; i++ ) {
    free(alignment[i]);
  }
  free(alignment);
  free(startSeq);
  clearTree(root_ptr);
}


int main () {
  const gsl_rng_type * T;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);

  testFloatEquals();
  checkMakeTree();// Check by seeing if you can print out the same tree with printTree
  checkCountTips();
  checkSimulateDNABranch(r);
  checkMakeInstantaneousRate();
  checkNewNodeSampleTree();
  checksimulateDNATree(r);
  checkReadFasta();
  checkSetStartSeq(r);
  checkMakeAlignment(r);
  checkReadLatent();


   // Check with valgrind
   // clearMemory();
   // clearTree();
   // NewickStringBranchLength(); // just for printing the tree- not used except for checking
   // checkInput(); // Don't want to check because it would exit the program
   printf("\nFinished\n");
   gsl_rng_free (r);

}
