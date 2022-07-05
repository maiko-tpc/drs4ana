#include <stdio.h>
#include <unistd.h>
#include <fcntl.h>
#include <iostream>
#include <fstream>
#include <arpa/inet.h>

#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>

#define N_ADC 16
#define N_CLK 1024
#define N_CLK_USE 1023

//////////// Analysis Parameters /////////////////
const int  N_ANA_ADC = 1;
const int  BASE_WID = 50;
const int  DECAY_FIT_MIN = 200;
const int  DECAY_FIT_MAX = 500;
//////////////////////////////////////////////////

TFile *fileout;
TTree *tree;

// Tree branch
unsigned int eve;
int adc[N_ADC][N_CLK_USE];
float adc_cor[N_ADC][N_CLK_USE];
float baseline[N_ADC];
float integ[N_ADC];
float max[N_ADC];
float decay_time[N_ADC];

void analysis(char *filename);
void anaevt();

int main(int iarg, char *argv[]) {
  if (iarg != 3) {
    printf("usage: ./anadrs <input> <output>\n");
    printf("eg.)\n");
    printf("$ ./anadrs 0123.bin 0123.root\n");
    exit(0);
  }

  fileout = new TFile(argv[2], "recreate");
  tree = new TTree("tree", "tree");

  tree->Branch("eve", &eve, "eve/I");
  tree->Branch("adc", adc, Form("adc[%d][%d]/I",N_ADC, N_CLK_USE));
  tree->Branch("baseline", baseline, Form("baseline[%d]/F",N_ADC));  
  tree->Branch("adc_cor", adc_cor, Form("adc_cor[%d][%d]/F",N_ADC, N_CLK_USE));
  tree->Branch("integ", integ, Form("integ[%d]/F",N_ADC));
  tree->Branch("max", max, Form("max[%d]/F",N_ADC));
  tree->Branch("decay_time", decay_time, Form("decay_time[%d]/F",N_ADC));      
  
  analysis(argv[1]);
  
  tree->Write();
  fileout->Close();

  return 0;
}

void analysis(char *filename) {
  unsigned int readcnt=0;
  short tmpbuf[2]={0,0};
  short tmpadc;
  unsigned short oldbuf=0;
  unsigned short buf;
  unsigned int cnt_header=0;

  eve=0;
  
  /* open data file */
  std::ifstream fdin;
  fdin.open(filename, std::ios_base::in|std::ios::binary);
  if (!fdin){
    printf("ERROR: Cannot open %s\n", filename);
    exit(0);
  }

  while(!fdin.eof()){  
    oldbuf = tmpbuf[0];
    fdin.read((char*)&tmpbuf, sizeof(short));
    readcnt ++;
    buf = tmpbuf[0];
    if(oldbuf==0x5555 && buf==0xaaaa){ //header
      cnt_header++;

      for(int i=0; i<N_ADC; i++){
	for(int j=0; j<N_CLK; j++){
	  oldbuf = tmpbuf[0];
	  fdin.read((char*)&tmpbuf, sizeof(short));
	  readcnt ++;
	  buf = tmpbuf[0];
	  if(j>0 && j<N_CLK_USE){
	    tmpadc = (short)htons(buf);
	    adc[i][j-1] = tmpadc;	    	    
	  }
	} // end of j loop
      } // end of i loop

      // detailed analysis
      anaevt();
      
      tree->Fill();
      eve++;
 
     if(eve%100==0) printf("Analyzed %d events.\n", eve);
    }

  }
  
  fdin.close();

  printf("Analyzed %d events.\n", eve);
  
}

void anaevt(){
  // Initialize parameters
  for(int i=0; i<N_ANA_ADC; i++){
    baseline[i]=0;
    integ[i]=0;
    max[i]=0;
    for(int j=0; j<N_CLK_USE; j++){
      adc_cor[i][j]=0;
    }
  }

  // baseline calculation
  for(int i=0; i<N_ANA_ADC; i++){
    for(int j=0; j<BASE_WID; j++){
      baseline[i] += adc[i][j];
    }
    baseline[i] = baseline[i]/((float)(BASE_WID));
  }  

  // baseline subtraction
  for(int i=0; i<N_ANA_ADC; i++){
    for(int j=0; j<N_CLK_USE; j++){
      adc_cor[i][j] = baseline[i] - adc[i][j];
      integ[i] += adc_cor[i][j];
    }
  }

  // pulse height search
  float tmp_max;
  for(int i=0; i<N_ANA_ADC; i++){
    tmp_max=0;
    for(int j=0; j<N_CLK_USE; j++){
      if(adc_cor[i][j] > tmp_max) tmp_max = adc_cor[i][j];
    }
    max[i] = tmp_max;
  }

  // define TGraph for fitting the decay time
  TGraph* gr_adc_cor[N_ADC];
  float clk[N_CLK_USE];
  for(int i=0; i<N_CLK; i++) clk[i] = i;

  for(int i=0; i<N_ADC; i++){
    //    gr_adc_cor[i] = new TGraph(Form("gr %d", i), N_CLK_USE, clk, adc_cor[i]);
    gr_adc_cor[i] = new TGraph(N_CLK_USE, clk, adc_cor[i]);    
  }
}
