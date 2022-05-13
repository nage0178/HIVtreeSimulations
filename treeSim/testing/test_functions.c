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

#include "../src/HIV_treeSim.h"
gsl_rng * r;  /* global generator */

/* Creates an array of node pointers of size arraySize. Creates new nodes for
lastIndex nodes and sets the rest of the array equal to null. Returns the pointer
to the array */
struct Node** setArrays (int lastIndex, int arraySize) {
  struct Node** array = malloc(sizeof(struct Node*) * arraySize);

  for (int i = 0; i < lastIndex; i++) {
    array[i] = newNode(i , 0, NULL);
  }

  for(int i = lastIndex; i < arraySize; i++) {
    array[i] = NULL;
  }

  return array;
}

/* Frees the memory in an array of pointers to node structs for both the node
structs and the memory of the arrray of pointers. lastIndex is the number of
nodes in the array. They must be the first lastIndex - 1 elements of the array */
void freeArrays(struct Node* array[], int lastIndex) {
  for (int i = 0; i < lastIndex; i++) {
    free(array[i]);
  }
  free(array);
  return;
}

// Checks the FloatEquals functions
void testFloatEquals() {
  assert(FloatEquals(1.1, 1.1001, 1e-4));
  assert(!(FloatEquals(1.1, 1.01, 1e-4)));
  assert(FloatEquals(2.0, 2.0 + 1e-15, 1e-14));
  assert(FloatEquals(-2.3, -2.30006, 1e-3));
  assert(!FloatEquals(-2.3, -2.1, 1e-3));
  return;
}

/* Checks whether all of the elements in an array of size are null.
Returns 0 if they are and 1 if they are not */
int checkNull(struct Node** array, int size) {
  int fail = 0;
  for(int i = 0; i < size; i++) {
    if (array[i] != NULL) {
      fail = 1;
      break;
    }
  }
  return fail;
}

/* Checks that setArrayNull actually set the array to Null*/
void testSetArrayNull() {
  //Makes two arrays, sets elements to NULL, checks all elements are NULL
  int arraySize1 = 100;
  long int arraySize2 = 1000;
  struct Node** array1 = malloc(sizeof(struct Node*) * arraySize1);
  struct Node** array2 = malloc(sizeof(struct Node*) * arraySize2);
  setArrayNull(array1, arraySize1);
  setArrayNull(array2, arraySize2);

  assert(checkNull(array1, arraySize1) == 0);
  assert(checkNull(array2, arraySize2) == 0);

  // Changes elements in the arrays
  struct Node* node1 = newNode(0, 0, NULL);
  struct Node* node2 = newNode(0, 0, NULL);
  array1[31] = node1;
  array1[0] = node1;
  array2[30] = node2;
  array2[50] = node1;
  array2[arraySize2 -1] = node2;

  // CheckNull should fail
  assert(checkNull(array1, arraySize1) == 1);
  assert(checkNull(array2, arraySize2) == 1);

  //Sets to null
  setArrayNull(array1, arraySize1);
  setArrayNull(array2, arraySize2);

  // checkNull should pass after having the elements set to zero
  assert(checkNull(array1, arraySize1) == 0);
  assert(checkNull(array2, arraySize2) == 0);

  free(node1);
  free(node2);
  free(array1);
  free(array2);
  return;
}

/* Checks that setVirusArrays sets all of the relevant arrays to NULL*/
void testSetVirusArrays() {
  long int maxActive = 1200000; /* per mL */
  long int maxVirus = 10000000;
  long int maxIncomp = 5000;
  long int maxComp = 5000;
  int mLBlood = 50;
  int totSampleVirus = 100;
  int totSampleLatent = 200;
  struct Node** activeArray       = malloc(sizeof(struct Node*) * maxActive * mLBlood);
  struct Node** virusArray        = malloc(sizeof(struct Node*) * maxVirus * mLBlood);
  struct Node** latentIncompArray = malloc(sizeof(struct Node*) * maxIncomp * mLBlood);
  struct Node** latentCompArray   = malloc(sizeof(struct Node*) * maxComp * mLBlood);
  struct Node** virusSampleArray  = malloc(sizeof(struct Node*) * totSampleVirus);
  struct Node** latentSampleArray = malloc(sizeof(struct Node*) * totSampleLatent);

  /* sets all of the elements in the arrays to NULL */
  setVirusArrays(activeArray, virusArray, latentIncompArray, latentCompArray,
  virusSampleArray, latentSampleArray, maxActive * mLBlood, maxVirus * mLBlood, maxIncomp * mLBlood,
  maxComp * mLBlood, totSampleVirus, totSampleLatent);

  /* Checks all of the elements in the arays are NULL */
  assert(checkNull(activeArray, maxActive * mLBlood) == 0);
  assert(checkNull(virusArray, maxVirus * mLBlood) == 0);
  assert(checkNull(latentIncompArray, maxIncomp * mLBlood) == 0);
  assert(checkNull(latentCompArray, maxComp * mLBlood) == 0);
  assert(checkNull(virusSampleArray, totSampleVirus) == 0);
  assert(checkNull(latentSampleArray, totSampleLatent) == 0);

  // Adds some elements to the arrays
  struct Node* node1 = newNode(0, 0, NULL);
  struct Node* node2 = newNode(0, 0, NULL);
  virusSampleArray[50] = node1;
  activeArray[300] = node2;
  activeArray[0] = node1;
  latentSampleArray[totSampleLatent - 1] = node2;
  /* The arrays should not all be NULL*/
  assert(checkNull(virusSampleArray, totSampleVirus) == 1);
  assert(checkNull(activeArray, maxActive * mLBlood) == 1);

  /* Sets all of the elements in the arrays to NULL */
  setVirusArrays(activeArray, virusArray, latentIncompArray, latentCompArray,
  virusSampleArray, latentSampleArray, maxActive, maxVirus, maxIncomp,
  maxComp, totSampleVirus, totSampleLatent);

  /* Checks all of the elements in the arays are NULL */
  assert(checkNull(activeArray, maxActive * mLBlood) == 0);
  assert(checkNull(virusArray, maxVirus * mLBlood) == 0);
  assert(checkNull(latentIncompArray, maxIncomp * mLBlood) == 0);
  assert(checkNull(latentCompArray, maxComp * mLBlood) == 0);
  assert(checkNull(virusSampleArray, totSampleVirus) == 0);
  assert(checkNull(latentSampleArray, totSampleLatent) == 0);

  free(activeArray);
  free(virusArray);
  free(latentIncompArray);
  free(latentCompArray);
  free(virusSampleArray);
  free(latentSampleArray);
  free(node1);
  free(node2);
  return;
}


void testNewNode () {
  /* Makes a tree. Makes the tree doublely linked */
  struct Node* node1 = newNode(0, 0, NULL);
  struct Node* node2 = newNode(.5, 1, node1);
  struct Node* node3 = newNode(.5, 0, node1);
  node1->left = node2;
  node1->right = node3;
  node1->left->left = newNode(1, 0, node1->left);

  /* Checks all variables in node 1 */
  assert(node1->birth_time == 0);
  assert(node1->sample_time == 0);
  assert(node1->isLatent == 0);
  assert(node1->timeToLatent == -1);
  assert(node1->totTimeLatent == 0);
  assert(node1->left == node2);
  assert(node1->right == node3);
  assert(node1->parent == NULL);

  /* Checks all variables in node 2 */
  assert(node2->birth_time == .5);
  assert(node2->sample_time == 0);
  assert(node2->isLatent == 1);
  assert(node2->timeToLatent == .5);
  assert(node2->totTimeLatent == 0);
  assert(node2->right == NULL);
  assert(node2->parent == node1);

  /* Checks all variables in node 3 */
  assert(node3->birth_time == .5);
  assert(node3->sample_time == 0);
  assert(node3->isLatent == 0);
  assert(node3->timeToLatent == -1);
  assert(node3->totTimeLatent == 0);
  assert(node3->left == NULL);
  assert(node3->right == NULL);
  assert(node3->parent == node1);

  /* Checks all variables in node 4 */
  assert(node2->left->birth_time == 1);
  assert(node2->left->sample_time == 0);
  assert(node2->left->isLatent == 0);
  assert(node2->left->timeToLatent == -1);
  assert(node2->left->totTimeLatent == 0);
  assert(node2->left->left == NULL);
  assert(node2->left->right == NULL);
  assert(node2->left->parent == node2);

  free(node2->left);
  free(node2);
  free(node1);
  free(node3);

  return;
}

/* Checks the first part of an array up to lastIndex
is not NULL and the second part is*/
void checkArrayEndNull(struct Node* testArray[], int lastIndex, int sizeArray) {
  for (int i = 0; i < lastIndex; i++) {
    assert(testArray[i] != NULL);
  }
  for (int i = lastIndex; i < sizeArray; i++){
    assert(testArray[i] == NULL);
  }
}

void testRemoveVirus () {
  /* Makes an array with node pointers */
  long int count = 20;
  long int tmpCount = count;
  int lastIndex = 8;
  int arraySize = 10;
  struct Node** testArray = setArrays(lastIndex, arraySize);
  struct Node* tmp_Node;
  struct Node* rm_Node;

  /* Removes a virus from the array */
  tmp_Node = testArray[lastIndex - 1];
  rm_Node = testArray[5];
  removeVirus(testArray, 5, lastIndex - 1, &count);
  assert(count == tmpCount - 1);
  lastIndex--;
  /* Checks the appropriate elements in the array are null */
  checkArrayEndNull(testArray, lastIndex, arraySize);
  /* Checks the last pointer in the array got moved to where the element that
  got removed was */
  assert(tmp_Node == testArray[5]);
  free(rm_Node);

  /* Removes a virus from the array */
  tmp_Node = testArray[lastIndex - 1];
  rm_Node = testArray[4];
  tmpCount = count;
  removeVirus(testArray, 4, lastIndex - 1, &count);
  assert(count == tmpCount - 1);
  lastIndex--;
  /* Checks the appropriate elements in the array are null */
  checkArrayEndNull(testArray, lastIndex, arraySize);
  /* Checks the last pointer in the array got moved to where the element that
  got removed was */
  assert(tmp_Node == testArray[4]);
  free(rm_Node);
  freeArrays(testArray, lastIndex);

  /* Makes an array with node pointers */
  lastIndex = 10;
  struct Node** testArray2 = setArrays(lastIndex, arraySize);

  /* Removes the last element from the array */
  rm_Node = testArray2[9];
  tmpCount = count;
  removeVirus(testArray2, 9, lastIndex - 1, &count);
  assert(count == tmpCount - 1);
  lastIndex--;
  /* Checks the appropriate elements in the array are null */
  checkArrayEndNull(testArray2, lastIndex, arraySize);
  free(rm_Node);

  freeArrays(testArray2, lastIndex);

  /* Makes an array with node pointers */
  lastIndex = 1;
  arraySize = 20;
  struct Node** testArray3 = setArrays(lastIndex, arraySize);

  /* Removes the first element in the array */
  rm_Node = testArray3[0];
  tmpCount = count;
  removeVirus(testArray3, 0, lastIndex - 1, &count);
  assert(count == tmpCount - 1);
  lastIndex--;
  /* Checks the whole array is null */
  checkArrayEndNull(testArray3, lastIndex, arraySize);
  free(rm_Node);
  free(testArray3);
  return;
}

void testReactivateLatent () {
  /* Makes arrays with node pointers */
  long int lastIndexLatent = 9;
  long int lastIndexActive = 20;
  int arraySize1 = 20;
  int arraySize2 = 40;
  struct Node** testLatentArray = setArrays(lastIndexLatent, arraySize1);
  struct Node** testActiveArray = setArrays(lastIndexActive, arraySize2);
  long int latentBefore, activeBefore;

  int latentIndex = 7;
  double totTime = 5.0;
  struct Node* toReactivate = testLatentArray[latentIndex];
  struct Node* lastLatentBefore = testLatentArray[lastIndexLatent - 1];

  /* Sets the timeToLatent as 0 for all viruses in the Latent array*/
  for (int i = 0; i < lastIndexLatent; i++) {
    testLatentArray[i]->timeToLatent = 0;
  }

  latentBefore = lastIndexLatent;
  activeBefore = lastIndexActive;
  /* Reactivate a latent virus */
  reactivateLatent(testLatentArray, testActiveArray, latentIndex, &lastIndexLatent, &lastIndexActive, totTime);

  /* Check index is changed */
  assert(latentBefore - 1 == lastIndexLatent);
  assert(activeBefore + 1 == lastIndexActive);

  /* Check the pointers in the arrays have been moved */
  assert(toReactivate == testActiveArray[lastIndexActive - 1]);
  assert(lastLatentBefore == testLatentArray[latentIndex]);
  assert(testLatentArray[lastIndexLatent] == NULL);
  /* Check the variables in the reactived node struct have been changed correctly */
  assert(testActiveArray[lastIndexActive - 1]->isLatent == 0);
  assert(testActiveArray[lastIndexActive - 1]->totTimeLatent == totTime);
  assert(testActiveArray[lastIndexActive - 1]->timeToLatent == -1);

  /* Checks the correct elements in the array are null */
  checkArrayEndNull(testLatentArray, lastIndexLatent, arraySize1);
  checkArrayEndNull(testActiveArray, lastIndexActive, arraySize2);

  /* Reactivates the last virus in the latent array */
  toReactivate = testLatentArray[latentIndex];
  totTime = 6;
  latentBefore = lastIndexLatent;
  activeBefore = lastIndexActive;

  reactivateLatent(testLatentArray, testActiveArray, latentIndex, &lastIndexLatent, &lastIndexActive, totTime);

  /* Check index is changed */
  assert(latentBefore - 1 == lastIndexLatent);
  assert(activeBefore + 1 == lastIndexActive);
  /* Check the pointers in the arrays have been moved */
  assert(toReactivate == testActiveArray[lastIndexActive - 1]);
  assert(testLatentArray[lastIndexLatent] == NULL);
  /* Check the variables in the reactived node struct have been changed correctly */
  assert(testActiveArray[lastIndexActive - 1]->isLatent == 0);
  assert(testActiveArray[lastIndexActive - 1]->totTimeLatent == totTime);
  assert(testActiveArray[lastIndexActive - 1]->timeToLatent == -1);

  /* Checks the correct elements in the array are null */
  checkArrayEndNull(testLatentArray, lastIndexLatent, arraySize1);
  checkArrayEndNull(testActiveArray, lastIndexActive, arraySize2);

  freeArrays(testLatentArray, lastIndexLatent);
  freeArrays(testActiveArray, lastIndexActive);

  /* Makes arrays of node pointers */
  lastIndexLatent = 1;
  lastIndexActive = 10;
  arraySize1 = 10;
  arraySize2 = 20;
  struct Node** testLatentArray2 = setArrays(lastIndexLatent, arraySize1);
  struct Node** testActiveArray2 = setArrays(lastIndexActive, arraySize2);
  latentIndex = 0;
  totTime = 7.0;
  testLatentArray2[0]->timeToLatent = 0;

  toReactivate = testLatentArray2[latentIndex];
  latentBefore = lastIndexLatent;
  activeBefore = lastIndexActive;

  /* Reactivates the only latent virus in the array */
  reactivateLatent(testLatentArray2, testActiveArray2, latentIndex, &lastIndexLatent, &lastIndexActive, totTime);

  /* Check index is changed */
  assert(latentBefore - 1 == lastIndexLatent);
  assert(activeBefore + 1 == lastIndexActive);
  /* Check the pointers in the arrays have been moved */
  assert(toReactivate == testActiveArray2[lastIndexActive - 1]);
  assert(testLatentArray2[lastIndexLatent] == NULL);
  /* Check the variables in the reactived node struct have been changed correctly */
  assert(testActiveArray2[lastIndexActive - 1 ]->isLatent == 0);
  assert(testActiveArray2[lastIndexActive - 1]->totTimeLatent == totTime);
  assert(testActiveArray2[lastIndexActive - 1]->timeToLatent == -1);

  /* Checks the correct elements in the array are null */
  checkArrayEndNull(testLatentArray2, lastIndexLatent, arraySize1);
  checkArrayEndNull(testActiveArray2, lastIndexActive, arraySize2);
  freeArrays(testActiveArray2, lastIndexActive);
  free(testLatentArray2);

  return;
}

void testBirthEvent () {
  double time = 3.3;
  int latentState = 1;

  long int lastIndexParent = 9;
  long int lastIndexDaughter = 20;
  int parentIndex = 5;
  unsigned long int totMem, tmpMem;
  long int tmp_lastIndexParent = lastIndexParent;
  long int tmp_lastIndexDaughter = lastIndexDaughter;


  /* Makes two arrays of node pointers */
  struct Node** testParentArray = setArrays(lastIndexParent, 20);
  struct Node** testDaughterArray = setArrays(lastIndexDaughter, 40);
  struct Node* parent = testParentArray[parentIndex];

  totMem = sizeof(testParentArray) + sizeof(testDaughterArray) + sizeof(struct Node)*(lastIndexParent +lastIndexDaughter);
  tmpMem = totMem;
  birthEvent(testParentArray, testDaughterArray, parentIndex, &lastIndexDaughter, time, latentState, &lastIndexParent, &totMem);
  assert(tmpMem + 2* (sizeof(struct Node)) == totMem);

  /* Check the parent is replaced by the left daughter in the array and all
  variables are asigned properly */
  assert(parent->left == testParentArray[parentIndex]);
  assert(FloatEquals(testParentArray[parentIndex]->birth_time, time, 1e-7));
  assert(testParentArray[parentIndex]->isLatent == 0);
  assert(testParentArray[parentIndex]->parent == parent);

  /* Check the right daughter is placed in the daughter array and all Variables
  are assinged properly */
  assert(parent->right == testDaughterArray[lastIndexDaughter - 1]);
  assert(FloatEquals(testDaughterArray[lastIndexDaughter - 1]->birth_time, time, 1e-7));
  assert(testDaughterArray[lastIndexDaughter - 1]->isLatent == latentState);
  assert(testDaughterArray[lastIndexDaughter - 1]->parent == parent);

  assert(tmp_lastIndexParent == lastIndexParent + 1);
  assert(tmp_lastIndexDaughter == lastIndexDaughter - 1);

  free(parent);
  freeArrays(testParentArray, tmp_lastIndexParent);
  freeArrays(testDaughterArray, lastIndexDaughter);
  return;
}

void testPruneTip () {

  unsigned long int totMem = sizeof(struct Node) * 8;

  /* Creates a tree ((4,5)2,(6,7)3)1*/
  struct Node* stem = newNode(0, 0, NULL);
  struct Node* root_ptr = newNode(0, 0, stem); /* Root of the tree */
  stem->left = root_ptr; /* Stem points to the root */

  root_ptr->totTimeLatent = .06;
  root_ptr->left = newNode(.5, 0, root_ptr);
  root_ptr->left->totTimeLatent = .05;
  root_ptr->right = newNode(.5, 0, root_ptr);

  root_ptr->left->left = newNode(.7, 0, root_ptr->left);
  root_ptr->left->right = newNode(.7, 0, root_ptr->left);

  root_ptr->right->left = newNode(.8, 0, root_ptr->right);
  root_ptr->right->right = newNode(.8, 0, root_ptr->right);

  /* Checks the tree structure is as expected */
  char *sampleTreeString = strdup("");
  long int nodeName = 0;
  sampleTreeString = NewickAll(stem->left, sampleTreeString, &nodeName);
  assert(strcmp(sampleTreeString, "((0,1)2,(3,4)5)6") == 0);
  free(sampleTreeString);

  /* Checks the tree structure after pruning. Checks varaibles have been changed
  correctly*/
  pruneTip(stem->left->left->left, &totMem); // Removes node 4 and parent node 2

  sampleTreeString = strdup("");
  nodeName = 0;
  sampleTreeString = NewickAll(stem->left, sampleTreeString, &nodeName);
  assert(strcmp(sampleTreeString, "(0,(1,2)3)4") == 0);
  free(sampleTreeString);

  assert(stem->left->left->birth_time == .5); // Birth time of node 2
  assert(FloatEquals(stem->left->left->totTimeLatent, .05, 1e-7)); // Time latent of node 2
  assert(stem->left->left->parent == stem->left);

  pruneTip(stem->left->left, &totMem); // Removes node 5 and parent node 1

  sampleTreeString = strdup("");
  nodeName = 0;
  sampleTreeString = NewickAll(stem->left, sampleTreeString, &nodeName);
  assert(strcmp(sampleTreeString, "(0,1)2") == 0);
  free(sampleTreeString);

  assert(FloatEquals(stem->left->totTimeLatent, .06, 1e-7));
  assert(stem->left->birth_time == 0);
  assert(stem->left->left->parent == stem->left);

  pruneTip(stem->left->right, &totMem); // Removes node 7 and parent node 3
  sampleTreeString = strdup("");
  nodeName = 0;
  sampleTreeString = NewickAll(stem->left, sampleTreeString, &nodeName);
  assert(strcmp(sampleTreeString, "0") == 0);
  assert(stem == stem->left->parent);

  free(stem->left);
  free(stem);
  free(sampleTreeString);

  return;
}


void testSample(){
  // Creates arrays with node structs
  long int lastIndex = 8;
  struct Node** testArray = setArrays (lastIndex, 10);
  struct Node** sampleArray = setArrays(0, 5);
  struct Node* tmp1;
  struct Node* tmp2;
  double totTime = 2.0;
  int sampleIndex = 0;

  tmp1 = testArray[5];
  tmp2 = testArray[lastIndex - 1];
  sample(testArray, sampleArray, 5, &sampleIndex, totTime, &lastIndex);

  // Checks pointers are moved/changed correctly
  assert(tmp1  == sampleArray[0]);
  assert(tmp2 == testArray[5]);
  assert(testArray[7] == NULL);
  //Checks other sampling related variables
  assert(sampleArray[0]->sample_time == totTime);
  assert(sampleArray[0]->isSample == 1);


  testArray[lastIndex - 1]->isLatent = 1;
  testArray[lastIndex - 1]->totTimeLatent = 5.0;
  testArray[lastIndex - 1]-> timeToLatent = 7.5;
  totTime = 15;

  sampleIndex = 1;
  tmp1 = testArray[lastIndex - 1];

  sample(testArray, sampleArray, lastIndex - 1, &sampleIndex, totTime, &lastIndex);

  // Checks pointers are moved/changed correctly
  assert(tmp1 == sampleArray[1]);
  assert(sampleArray[1]->isLatent == 1);
  assert(testArray[lastIndex] == NULL);
  //Checks other sampling related variables
  assert(sampleArray[1]->sample_time == totTime);
  assert(sampleArray[1]->isSample == 1);
  assert(sampleArray[1]->totTimeLatent == 7.5 + 5);

  freeArrays(testArray, lastIndex);
  freeArrays(sampleArray, 2);

  return;
}

void testCalcProb() {
  double rates[11];
  double totRate;
  double probs[11];
  double calcTotRate = 0;
  for(int i =0; i < 11; i++) {
    rates[i] = (i+1) * .1;
    calcTotRate = calcTotRate + rates[i];
  }

  CalcProb(rates, probs,totRate);
  assert(FloatEquals(totRate, calcTotRate, 1e-14));

  for (int i = 0; i < 11; i++) {
    assert(FloatEquals(probs[i], rates[i]/calcTotRate, 1e-14));
  }

  return;
}

void testSampleEvent() {

  /*Set up random number generator */
  const gsl_rng_type * T;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);

  long int numVirus = 95;
  long int numLatentComp = 9;
  long int numLatentIncomp = 18;
  int numVirusSample = 0;
  int numLatentSample = 0;
  double totTime = 4.0;

  long int tmp_numVirus = numVirus;
  long int tmp_numLatentComp = numLatentComp;
  long int tmp_numLatentIncomp = numLatentIncomp;
  int tmp_numVirusSample = numVirusSample;
  int tmp_numLatentSample = numLatentSample;

  int numSampleVirus[3] = {5, 10, 10};
  int numSampleLatent[3] = {0, 0, 5};

  int arraySize1 = 100;
  int arraySize2 = 10;
  int arraySize3 = 20;
  int arraySize4 = 25;
  int arraySize5 = 5;

  int sampleCounter = 0;
  int tmp_sampleCounter = sampleCounter;

  struct Node** virusArray = setArrays(numVirus, arraySize1);
  struct Node** latentCompArray = setArrays(numLatentComp, arraySize2);
  struct Node** latentIncompArray = setArrays(numLatentIncomp, arraySize3);
  struct Node** sampleActiveArray = setArrays(numVirusSample, arraySize4);
  struct Node** sampleLatentArray = setArrays(numLatentSample, arraySize5);
  //Need to set latent times

  for(int i = 0; i < numLatentComp; i++) {
    latentCompArray[i]->timeToLatent = 3;
  }

  for(int i = 0; i < numLatentIncomp; i++) {
    latentIncompArray[i]->timeToLatent = 2;
  }
  //SampleEvent 1

  sampleEvent(r, &sampleCounter, numSampleVirus, numSampleLatent,
  virusArray, sampleActiveArray, sampleLatentArray, latentIncompArray, latentCompArray,
  &numVirusSample, &numLatentSample, &numVirus, &numLatentIncomp, &numLatentComp, totTime);


  assert(tmp_numVirus - numSampleVirus[0] == numVirus);
  assert(tmp_numLatentComp == numLatentComp);
  assert(tmp_numLatentIncomp == numLatentIncomp);
  assert(tmp_numVirusSample + numSampleVirus[0] == numVirusSample);
  assert(tmp_numLatentSample == numLatentSample);
  assert(tmp_sampleCounter + 1 == sampleCounter);

  checkArrayEndNull(virusArray, numVirus, arraySize1);
  checkArrayEndNull(latentCompArray, numLatentComp, arraySize2);
  checkArrayEndNull(latentIncompArray, numLatentIncomp, arraySize3);
  checkArrayEndNull(sampleActiveArray, numVirusSample, arraySize4);
  checkArrayEndNull(sampleLatentArray, numLatentSample, arraySize5);

  //Sample Event 2
  tmp_numVirus = numVirus;
  tmp_numLatentComp = numLatentComp;
  tmp_numLatentIncomp = numLatentIncomp;
  tmp_numVirusSample = numVirusSample;
  tmp_numLatentSample = numLatentSample;
  tmp_sampleCounter = sampleCounter;

  sampleEvent(r, &sampleCounter, numSampleVirus, numSampleLatent,
  virusArray, sampleActiveArray, sampleLatentArray, latentIncompArray, latentCompArray,
  &numVirusSample, &numLatentSample, &numVirus, &numLatentIncomp, &numLatentComp, totTime);


  assert(tmp_numVirus - numSampleVirus[1] == numVirus);
  assert(tmp_numLatentComp == numLatentComp);
  assert(tmp_numLatentIncomp == numLatentIncomp);
  assert(tmp_numVirusSample + numSampleVirus[1] == numVirusSample);
  assert(tmp_numLatentSample == numLatentSample);
  assert(tmp_sampleCounter + 1 == sampleCounter);

  checkArrayEndNull(virusArray, numVirus, arraySize1);
  checkArrayEndNull(latentCompArray, numLatentComp, arraySize2);
  checkArrayEndNull(latentIncompArray, numLatentIncomp, arraySize3);
  checkArrayEndNull(sampleActiveArray, numVirusSample, arraySize4);
  checkArrayEndNull(sampleLatentArray, numLatentSample, arraySize5);

  //Sample event 3
  tmp_numVirus = numVirus;
  tmp_numLatentComp = numLatentComp;
  tmp_numLatentIncomp = numLatentIncomp;
  tmp_numVirusSample = numVirusSample;
  tmp_numLatentSample = numLatentSample;
  tmp_sampleCounter = sampleCounter;

  sampleEvent(r, &sampleCounter, numSampleVirus, numSampleLatent,
  virusArray, sampleActiveArray, sampleLatentArray, latentIncompArray, latentCompArray,
  &numVirusSample, &numLatentSample, &numVirus, &numLatentIncomp, &numLatentComp, totTime);


  assert(tmp_numVirus - numSampleVirus[2] == numVirus);
  assert(tmp_numLatentComp + tmp_numLatentIncomp - numSampleLatent[2] == numLatentComp + numLatentIncomp);
  assert(tmp_numVirusSample + numSampleVirus[2] == numVirusSample);
  assert(tmp_numLatentSample + numSampleLatent[2] == numLatentSample);
  assert(tmp_sampleCounter + 1 == sampleCounter);

  checkArrayEndNull(virusArray, numVirus, arraySize1);
  checkArrayEndNull(latentCompArray, numLatentComp, arraySize2);
  checkArrayEndNull(latentIncompArray, numLatentIncomp, arraySize3);
  checkArrayEndNull(sampleActiveArray, numVirusSample, arraySize4);
  checkArrayEndNull(sampleLatentArray, numLatentSample, arraySize5);


  freeArrays(virusArray, numVirus);
  freeArrays(latentCompArray, numLatentComp);
  freeArrays(latentIncompArray, numLatentIncomp);
  freeArrays(sampleActiveArray, numVirusSample);
  freeArrays(sampleLatentArray, numLatentSample);
  gsl_rng_free (r);

  return;
}

void testReset() {

  // Makes tree
  /* Creates a tree ((4,5)2,(6,7)3)1*/
  struct Node* stem = newNode(0, 0, NULL);
  struct Node* root_ptr = newNode(0, 0, stem); /* Root of the tree */
  stem->left = root_ptr; /* Stem points to the root */
  root_ptr->totTimeLatent = .06;
  root_ptr->left = newNode(.5, 0, root_ptr);
  root_ptr->left->totTimeLatent = .05;

  root_ptr->right = newNode(.5, 0, root_ptr);

  root_ptr->left->left = newNode(.7, 0, root_ptr->left);
  root_ptr->left->right = newNode(.7, 0, root_ptr->left);

  root_ptr->right->left = newNode(.8, 0, root_ptr->right);
  root_ptr->right->right = newNode(.8, 0, root_ptr->right);

  //
  int sampleCounter = 5;
  unsigned long int numEvents = 15;

  int mLBlood = 50;
  long int maxActive = 300000 * mLBlood; /* per mL */
  long int maxVirus = 2100000 * mLBlood;
  long int maxIncomp = 4000 * mLBlood;
  long int maxComp = 5000* mLBlood;

  int totSampleVirus = 100;
  int totSampleLatent = 200;
  struct Node** activeArray       = malloc(sizeof(struct Node*) * maxActive);
  struct Node** virusArray        = malloc(sizeof(struct Node*) * maxVirus);
  struct Node** latentIncompArray = malloc(sizeof(struct Node*) * maxIncomp);
  struct Node** latentCompArray   = malloc(sizeof(struct Node*) * maxComp);
  struct Node** virusSampleArray  = malloc(sizeof(struct Node*) * totSampleVirus);
  struct Node** latentSampleArray = malloc(sizeof(struct Node*) * totSampleLatent);

  long unsigned int totMem = sizeof(activeArray) + sizeof(virusArray) + sizeof(latentIncompArray) + sizeof(latentCompArray)
  + sizeof(virusSampleArray) + sizeof(latentSampleArray) + sizeof(struct Node) * 8;

  activeArray[10] = root_ptr->left;
  virusArray[15] = root_ptr->right;
  latentIncompArray[20] = root_ptr->left->left;
  latentCompArray[25] = root_ptr->left->right;
  virusSampleArray[2] = root_ptr->right->left;
  latentSampleArray[1] = root_ptr->right->right;


  root_ptr = reset(stem, root_ptr, &sampleCounter, &numEvents, &totMem,
  activeArray, virusArray, latentIncompArray, latentCompArray,
  virusSampleArray, latentSampleArray, maxActive, maxVirus, maxIncomp,
  maxComp, totSampleVirus, totSampleLatent);


  checkArrayEndNull(activeArray, 1, maxActive);
  assert(checkNull(virusArray,        maxVirus)  == 0);
  assert(checkNull(latentIncompArray, maxIncomp) == 0);
  assert(checkNull(latentCompArray,   maxComp)   == 0);
  assert(checkNull(virusSampleArray,  totSampleVirus)      == 0);
  assert(checkNull(latentSampleArray, totSampleLatent)     == 0);

  assert(sampleCounter == 0);
  assert(numEvents == 0);
  assert(stem->left == root_ptr);
  assert(stem == root_ptr->parent);
  assert(totMem == sizeof(activeArray) + sizeof(virusArray) + sizeof(latentIncompArray) + sizeof(latentCompArray)
  + sizeof(virusSampleArray) + sizeof(latentSampleArray) + sizeof(struct Node) * 2);

  clearTreeMemory(root_ptr, &totMem);
  free(stem);
  free(activeArray);
  free(virusArray);
  free(latentIncompArray);
  free(latentCompArray);
  free(virusSampleArray);
  free(latentSampleArray);

  return;
}

void testNewickSampled() {
  struct Node** sampleArray = malloc(sizeof(struct Node*) * 2);
  struct Node** sampleArrayLatent = malloc(sizeof(struct Node*) * 1);
  struct Node* stem = newNode(0, 0, NULL);
  struct Node* root_ptr = newNode(0, 0, stem); /* Root of the tree */
  stem->left = root_ptr; /* Stem points to the root */
  unsigned long int totMem = sizeof(newNode) * 10;

  root_ptr->left = newNode(.5, 0, root_ptr);
  root_ptr->right = newNode(.5, 0, root_ptr);

  root_ptr->left->left = newNode(.6, 0, root_ptr->left);
  root_ptr->left->right = newNode(.6, 0, root_ptr->left);

  root_ptr->right->left = newNode(.75, 0, root_ptr->right);
  root_ptr->right->right = newNode(.75, 0, root_ptr->right);

  root_ptr->right->left->left = newNode(.8, 0, root_ptr->right->left);
  root_ptr->right->left->right = newNode(.8, 0, root_ptr->right->left);

  root_ptr->right->left->right->isSample = 1;
  root_ptr->left->left->isSample = 1;
  root_ptr->right->left->left->isSample = 1;

  root_ptr->right->left->right->sample_time = 1;
  root_ptr->left->left->sample_time = 1;
  root_ptr->right->left->left->sample_time = 1.5;

  sampleArray[0] = root_ptr->right->left->right;
  sampleArray[1] = root_ptr->left->left;
  sampleArrayLatent[0] = root_ptr->right->left->left;

  findSampledAll(sampleArray, sampleArrayLatent, 2, 1, stem);

  long int nodeName = 0;
  char *sampleTreeString = strdup("");
  sampleTreeString = NewickSampled(stem->left, sampleTreeString, 0, &nodeName);
  assert(strcmp("(1:0.500000,(3:0.700000,4:0.200000)2:0.300000)0:0.500000", sampleTreeString) == 0);
  free(sampleTreeString);
  clearTreeMemory(stem->left, &totMem);
  free(stem);
  free(sampleArray);
  free(sampleArrayLatent);

  return;
}

void testResetCounts() {
  long int numVirus = 500;
  long int numVirusInit = 0;
  long int numCellUninfect = 400000;
  long int numCellUninfectInit = 500000;
  long int numCellInfect = 3000;
  long int numCellInfectInit = 1;
  long int numLatentComp = 200;
  long int numLatentCompInit = 0;
  long int numLatentIncomp = 600;
  long int numLatentIncompInit = 0;
  int numLatentSample = 20;
  int numVirusSample = 100;
  int sampleCounter = 4;
  double totTime = 50.2;
  resetCounts(&numVirus, numVirusInit, &numCellUninfect, numCellUninfectInit,
      &numCellInfect, numCellInfectInit, &numLatentComp, numLatentCompInit,
      &numLatentIncomp, numLatentIncompInit, &numLatentSample, &numVirusSample, &sampleCounter, &totTime);

  assert(numVirus == 0);
  assert(numCellUninfect == 500000);
  assert(numCellInfect == 1);
  assert(numLatentComp == 0);
  assert(numLatentIncomp == 0);
  assert(numLatentSample == 0);
  assert(numVirusSample == 0);
  assert(sampleCounter == 0);
  assert(totTime == 0);

  return;
}

void testFindSampledParent () {

  struct Node* stem = newNode(0, 0, NULL);
  struct Node* root_ptr = newNode(0, 0, stem); /* Root of the tree */
  stem->left = root_ptr; /* Stem points to the root */
  unsigned long int totMem = sizeof(newNode) * 10;

  root_ptr->left = newNode(.5, 0, root_ptr);
  root_ptr->right = newNode(.5, 0, root_ptr);

  root_ptr->left->left = newNode(.5, 0, root_ptr->left);
  root_ptr->left->right = newNode(.5, 0, root_ptr->left);

  root_ptr->right->left = newNode(.5, 0, root_ptr->right);
  root_ptr->right->right = newNode(.5, 0, root_ptr->right);

  root_ptr->right->left->left = newNode(.5, 0, root_ptr->right->left);
  root_ptr->right->left->right = newNode(.5, 0, root_ptr->right->left);

  root_ptr->right->left->right->isSample = 1;

  findSampledParent(root_ptr->right->left->right->parent, stem);

  assert(root_ptr->right->left->right->isSample == 1);
  assert(root_ptr->right->left->isSample == 1);
  assert(root_ptr->right->isSample == 1);
  assert(root_ptr->isSample == 1);
  assert(stem->isSample == 1);

  clearTreeMemory(root_ptr, &totMem);
  free(stem);
}

void testFindSampledAll () {
  struct Node** sampleArray = malloc(sizeof(struct Node*) * 2);
  struct Node** sampleArrayLatent = malloc(sizeof(struct Node*) * 1);
  struct Node* stem = newNode(0, 0, NULL);
  struct Node* root_ptr = newNode(0, 0, stem); /* Root of the tree */
  stem->left = root_ptr; /* Stem points to the root */
  unsigned long int totMem = sizeof(newNode) * 10;

  root_ptr->left = newNode(.5, 0, root_ptr);
  root_ptr->right = newNode(.5, 0, root_ptr);

  root_ptr->left->left = newNode(.5, 0, root_ptr->left);
  root_ptr->left->right = newNode(.5, 0, root_ptr->left);

  root_ptr->right->left = newNode(.5, 0, root_ptr->right);
  root_ptr->right->right = newNode(.5, 0, root_ptr->right);

  root_ptr->right->left->left = newNode(.5, 0, root_ptr->right->left);
  root_ptr->right->left->right = newNode(.5, 0, root_ptr->right->left);

  root_ptr->right->left->right->isSample = 1;
  root_ptr->left->left->isSample = 1;
  root_ptr->right->left->left->isSample = 1;

  sampleArray[0] = root_ptr->right->left->right;
  sampleArray[1] = root_ptr->left->left;
  sampleArrayLatent[0] = root_ptr->right->left->left;

  findSampledAll(sampleArray, sampleArrayLatent, 2, 1, stem);

  assert(root_ptr->right->left->right->isSample == 1);
  assert(root_ptr->right->left->isSample == 1);
  assert(root_ptr->right->isSample == 1);
  assert(root_ptr->isSample == 1);


  assert(root_ptr->right->left->left->isSample == 1);
  assert(root_ptr->left->left->isSample == 1);
  assert(root_ptr->left->isSample == 1);

  assert(root_ptr->left->right->isSample == 2);

  assert(root_ptr->right->right->isSample == 2);
  assert(stem->isSample == 1);


  clearTreeMemory(root_ptr, &totMem);
  free(stem);
  free(sampleArray);
  free(sampleArrayLatent);
}

void testReadSampleTimes() {
  char sampleTimesFile[] = "testing/testReadSampleTimes.csv";
  int sampleFileLength = findFileLength(sampleTimesFile);
  assert(sampleFileLength == 5);

  double sampleTimes[sampleFileLength];
  int numSampleVirus[sampleFileLength];
  int numSampleLatent[sampleFileLength];

  readSampleTimes(sampleTimes, numSampleVirus, numSampleLatent, sampleTimesFile, sampleFileLength);

  assert(sampleTimes[0] == 5);
  assert(sampleTimes[1] == 10);
  assert(sampleTimes[2] == 25);
  assert(sampleTimes[3] == 150);
  assert(sampleTimes[4] == 200);

  assert(numSampleVirus[0] == 20);
  assert(numSampleVirus[1] == 11);
  assert(numSampleVirus[2] == 90);
  assert(numSampleVirus[3] == 0);
  assert(numSampleVirus[4] == 10);

  assert(numSampleLatent[0] == 0);
  assert(numSampleLatent[1] == 0);
  assert(numSampleLatent[2] == 12);
  assert(numSampleLatent[3] == 11);
  assert(numSampleLatent[4] == 40);
}

void testLatentSampled(){
  struct Node** sampleArray = malloc(sizeof(struct Node*) * 3);
  struct Node** sampleArrayLatent = malloc(sizeof(struct Node*) * 1);
  struct Node* stem = newNode(0, 0, NULL);
  struct Node* root_ptr = newNode(0, 0, stem); /* Root of the tree */
  stem->left = root_ptr; /* Stem points to the root */
  unsigned long int totMem = sizeof(newNode) * 10;

  root_ptr->left = newNode(.5, 0, root_ptr);
  root_ptr->right = newNode(.5, 0, root_ptr);

  root_ptr->left->left = newNode(.6, 0, root_ptr->left);
  root_ptr->left->right = newNode(.6, 0, root_ptr->left);

  root_ptr->right->left = newNode(.75, 0, root_ptr->right);
  root_ptr->right->right = newNode(.75, 0, root_ptr->right);

  root_ptr->right->left->left = newNode(.8, 0, root_ptr->right->left);
  root_ptr->right->left->right = newNode(.8, 0, root_ptr->right->left);

  root_ptr->right->left->right->isSample = 1;
  root_ptr->left->left->isSample = 1;
  root_ptr->right->left->left->isSample = 1;

  root_ptr->right->left->right->sample_time = 1;
  root_ptr->left->left->sample_time = 1;
  root_ptr->right->left->left->sample_time = 1.5;

  root_ptr->left->totTimeLatent = .3;
  root_ptr->right->left->totTimeLatent= .01;
  root_ptr->right->left->left->totTimeLatent = .1;
  root_ptr->right->left->left->isLatent = 1;

  root_ptr->right->right->totTimeLatent = .03;
  root_ptr->right->right->left = newNode(1, 0, root_ptr->right->right);
  root_ptr->right->right->right = newNode(1, 0, root_ptr->right->right);
  root_ptr->right->right->right->totTimeLatent = .02;

  root_ptr->right->right->right->left = newNode(1.2, 0, root_ptr->right->right->right);
  root_ptr->right->right->right->right = newNode(1.2, 0, root_ptr->right->right->right);
  root_ptr->right->right->right->right->isSample = 1;

  sampleArray[0] = root_ptr->right->left->right;
  sampleArray[1] = root_ptr->left->left;
  sampleArrayLatent[0] = root_ptr->right->left->left;
  sampleArray[2] = root_ptr->right->right->right->right;

  findSampledAll(sampleArray, sampleArrayLatent, 3, 1, stem);

  long int nodeName = 0;
  char *latentString = strdup("");
  latentString = LatentSampled(stem->left, latentString, 0, &nodeName) ;
  assert(strcmp("0:0.000000,0.000000,0\n1:0.300000,1.000000,0\n2:0.000000,0.000000,0\n3:0.010000,0.000000,0\n4:0.100000,1.500000,1\n5:0.000000,1.000000,0\n6:0.050000,0.000000,0\n",
  latentString) == 0);

  free(latentString);
  clearTreeMemory(stem->left, &totMem);
  free(stem);
  free(sampleArray);
  free(sampleArrayLatent);

}

void testSetInfectionRates() {
  double rates[11];
  long int numInfected, numVirus;
  numVirus = 123456;
  numInfected = 987;

  double prodInfectionRates[3];
  prodInfectionRates[0] = 0.0001234;
  prodInfectionRates[1] = 0.0000044;
  prodInfectionRates[2] = 0.000053;
  setInfectionRates(rates, prodInfectionRates, numInfected, numVirus);

  assert(FloatEquals(rates[1], numVirus * numInfected * prodInfectionRates[0], 1e-14));
  assert(FloatEquals(rates[6], numVirus * numInfected * prodInfectionRates[1], 1e-14));
  assert(FloatEquals(rates[7], numVirus * numInfected * prodInfectionRates[2], 1e-14));
}

void testReallocArray() {
  long int sizeArray = 100;
  long int oldSizeArray = sizeArray;
  unsigned long int totMem = 34 + sizeof(struct Node*) * sizeArray;
  unsigned long int oldTotmem = totMem;
  struct Node** myArray = setArrays (sizeArray, sizeArray);

  struct Node* index50 = myArray[50];
  struct Node* index30 = myArray[30];
  struct Node* index0 = myArray[0];
  struct Node* index99 = myArray[99];

  myArray = reallocArray(&sizeArray, myArray, &totMem);
  assert(sizeArray == oldSizeArray *  1.2 );
  assert(index50 == myArray[50]);
  assert(index30 == myArray[30]);
  assert(index99 == myArray[99]);
  assert(index0 == myArray[0]);
  assert(totMem == oldTotmem + oldSizeArray * sizeof(struct Node*) * 0.2 );

  for(int i = 0; i < oldSizeArray; i++) {
    free(myArray[i]);
  }
  free(myArray);
}

void testReallocArrayDecrease() {
  long int sizeArray = 120;
  long int oldSizeArray = sizeArray;
  unsigned long int totMem = 34 + sizeof(struct Node*) * sizeArray;
  unsigned long int oldTotmem = totMem;
  struct Node** myArray = setArrays (100, sizeArray);

  struct Node* index50 = myArray[50];
  struct Node* index30 = myArray[30];
  struct Node* index0 = myArray[0];
  struct Node* index99 = myArray[99];

  myArray = reallocArrayDecrease(&sizeArray, myArray, &totMem);
  assert(sizeArray == oldSizeArray *  5/6);
  assert(index50 == myArray[50]);
  assert(index30 == myArray[30]);
  assert(index99 == myArray[99]);
  assert(index0 == myArray[0]);
  assert(totMem == oldTotmem - oldSizeArray * sizeof(struct Node*) * 1/6 );

  for(int i = 0; i < 100; i++) {
    free(myArray[i]);
  }
  free(myArray);
}

void testSetRates() {
  double prodInfectinRates[3];
  prodInfectinRates[0] = 1.1;
  prodInfectinRates[1] = 2.2;
  prodInfectinRates[2] = 3.3;

  double rates[11];
  double parameters[6];
  for(int i = 0; i < 6; i++) {
    parameters[i] = 0.01 * (i + 1);
  }

  double reactLatent = .7;
  double latIncompDeath = .8;
  double latCompDeath = .9;
  long int numVirus = 123456;
  long int numCellInfect = 111;
  long int numCellUninfect = 1123456;
  long int numLatentComp = 2;
  long int numLatentIncomp = 10;

  setRates(prodInfectinRates, parameters, reactLatent, latIncompDeath, latCompDeath, rates, numVirus, numCellInfect, numCellUninfect, numLatentComp, numLatentIncomp);

  assert(FloatEquals(rates[0], parameters[0], 1e-14));
  assert(FloatEquals(rates[1], prodInfectinRates[0] * numVirus * numCellUninfect, 1e-14));
  assert(FloatEquals(rates[2], parameters[2] * numCellUninfect, 1e-14));
  assert(FloatEquals(rates[3], parameters[3] * numCellInfect, 1e-14));
  assert(FloatEquals(rates[4], parameters[4] * numCellInfect, 1e-14));
  assert(FloatEquals(rates[5], parameters[5] * numVirus, 1e-14));
  assert(FloatEquals(rates[6], prodInfectinRates[1] * numCellUninfect * numVirus, 1e-14));
  assert(FloatEquals(rates[7], prodInfectinRates[2] * numCellUninfect * numVirus, 1e-14));
  assert(FloatEquals(rates[8], latIncompDeath * numLatentIncomp, 1e-14));
  assert(FloatEquals(rates[9], latCompDeath * numLatentComp, 1e-14));
  assert(FloatEquals(rates[10], reactLatent * numLatentComp, 1e-14));

}
int main(){

  testSample();
  testFloatEquals();
  testSetArrayNull();
  testNewNode();
  testRemoveVirus();
  testReactivateLatent();
  testReset();
  testSetInfectionRates();
  testResetCounts();
  testNewickSampled();
  testReadSampleTimes();
  testLatentSampled();
  testSampleEvent();
  testFindSampledParent();
  testSetVirusArrays();
  testBirthEvent();
  testFindSampledAll();
  testPruneTip();
  testReallocArray();
  testReallocArrayDecrease();
  testSetRates();
  testCalcProb();


  printf("Finished\n");
  //Testing printTree - not actually going to use
  //Testing ClearMemory- run valgrind to test
  //Testing clearTreeMemory - not necessary- run valgrind

  return 0;
}
