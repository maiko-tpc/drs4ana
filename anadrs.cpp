#include <stdio.h>
#include <unistd.h>
#include <fcntl.h>
#include <fstream>

#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>

int main(int iarg, char *argv[]) {
  if (iarg != 3) {
    printf("usage: ./anadrs <input> <output>\n");
    printf("eg.)\n");
    printf("$ ./anadrs 0123.rdf 0123.root\n");
    exit(0);
  }

  return 0;
}
