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


#include <stdio.h>
#include <stdlib.h>
#include "dna_sim.h"
#include "string_utilities.h"
void findLastSampleTime(struct NodeSampledTree* node, int* nodeName, double* maxSampleTime);
int findHeight(struct NodeSampledTree* node, int nodeName, double* height);
double findRootAge (struct NodeSampledTree* node) ;
void treeLength(struct NodeSampledTree* node, double* length);
double percentTimeLatent(struct NodeSampledTree* root);
void totLatentTimeTree(struct NodeSampledTree* node, double* timeLatent);

/* Finds the youngest sample time */
void findLastSampleTime(struct NodeSampledTree* node, int* nodeName, double* maxSampleTime) {
  if(node->left == NULL) {
    if(node->sample_time > *maxSampleTime) {
      *nodeName = node->name;
      *maxSampleTime = node->sample_time;
    }
  } else {
    findLastSampleTime(node->left, nodeName, maxSampleTime);
    findLastSampleTime(node->right, nodeName, maxSampleTime);
  }
  return;
}

/* Finds the height of the tree by looking for the node with the youngest sampling time
and adding the branch lengths back from it */
int findHeight(struct NodeSampledTree* node, int nodeName, double* height) {
  int left, right;
  if(node->left == NULL) {
    if(node->name == nodeName) {
      *height = node->branchLength;
      return 1;
    } else {
      return 0;
    }
  } else {
    left = findHeight(node->left, nodeName, height);
    right = findHeight(node->right, nodeName, height);
  }
  if (left || right) {
    *height = *height + node->branchLength;
    return 1;
  } else {
    return 0;
  }
}

/* Finds the root age */ 
double findRootAge (struct NodeSampledTree* node) {
  int nodeName;
  double maxSampleTime, height;
  findLastSampleTime(node, &nodeName, &maxSampleTime);
  findHeight(node, nodeName, &height);
  /* Have to subtract the branch length from the root, or all tree with the same sampling times
  will have the same ages. We want the age at the root node*/
  return height - node->branchLength;

}

/* Calculates the percent of time spend in a latent state in the tree. Not used right now */
double percentTimeLatent(struct NodeSampledTree* root) {
  double length = 0;
  treeLength(root, &length);
  double timeLatent = 0;
  totLatentTimeTree(root, &timeLatent);
  return timeLatent / length;
}

/* Calculates the tree length (in days)*/
void treeLength(struct NodeSampledTree* node, double* length) {
  if(node->left == NULL) {
    *length = *length + node->branchLength;
  } else {
    treeLength(node->left, length);
    treeLength(node->right, length );
    *length = *length + node->branchLength;

  }

}

/* Calculates the total amount of time latent in the tree (in days)*/
void totLatentTimeTree(struct NodeSampledTree* node, double* timeLatent) {
  if(node->left == NULL) {
    *timeLatent = *timeLatent +  node->totTimeLatent;
  } else {
    totLatentTimeTree(node->left, timeLatent);
    totLatentTimeTree(node->right, timeLatent);
    *timeLatent = *timeLatent + node->totTimeLatent;
  }

}

int main (int argc, char **argv) {
  char* treeFile = argv[1];
  char* latentHistoryFile = argv[2];
  char *treeString = string_from_file(treeFile);

  /*Make a tree */
  /* Note that  birth times are actually branch lengths-figure out if they are where you want them to be */
  struct NodeSampledTree* root_ptr = newNodeSampledTree(0, 1);
  int endChar = makeTree(root_ptr, treeString, 0);
  if (treeString[endChar] != '\0') {
    fprintf(stderr, "Problem reading in tree. Did not reach the end of the file \n");
    exit(1);
  }

  FILE *latentFile;
  latentFile = fopen(latentHistoryFile, "r");
  readLatent(root_ptr, latentFile);
  fclose(latentFile);

  double length = 0 ;
  treeLength(root_ptr, &length);
  printf("Tree Length: %f\n", length);

  double timeLatent = 0 ;
  totLatentTimeTree(root_ptr ,&timeLatent);
  printf("Total time Latent: %f\n", timeLatent);

  printf("Percent time latent: %f\n", percentTimeLatent(root_ptr));

  double rootage = findRootAge(root_ptr);
  printf("Root age: %f\n", rootage);


}
