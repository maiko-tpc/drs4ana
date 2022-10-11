#include <stdio.h>
#include <unistd.h>
#include <fcntl.h>
#include <iostream>
#include <fstream>
#include <arpa/inet.h>

#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TROOT.h>
#include <TCanvas.h>
#include "TApplication.h"

#define N_ADC 16
#define N_CLK 1024
#define N_CLK_USE 1023

//////////// Analysis Parameters /////////////////
const int   N_ANA_ADC = 1;
const int   BASE_WID = 10;
const int   DECAY_FIT_MIN = 180;
const int   DECAY_FIT_MAX = 600;
const float INTEG_TH = 100;
//////////////////////////////////////////////////

TFile *fileout;
TTree *tree;

// Tree branch
unsigned int eve;
int adc[N_ADC][N_CLK_USE];
float adc_cor[N_ADC][N_CLK_USE];
float baseline[N_ADC];
float integ[N_ADC];
float integ_gate[N_ADC];
int   pulse_wid[N_ADC];
float max[N_ADC];
float decay_time[N_ADC];

TCanvas* c1;

TH2F *h_sum[N_ADC];


void analysis(char *filename);
void anaevt();


int main(int iarg, char *argv[]) {

  for(int i=0; i<N_ADC; i++){
    h_sum[i] = new TH2F(Form("h%d", i+10), Form("sum of ch%d", i),
			N_CLK_USE, 0, N_CLK_USE, 4096, -2048, 2048);
  }

  if (iarg != 3) {
    printf("iarg=%d\n", iarg);
    printf("usage: ./anadrs <input> <output>\n");
    printf("eg.)\n");
    printf("$ ./anadrs 0123.bin 0123.root\n");
    exit(0);
  }

  
  fileout = new TFile(argv[2], "recreate");
  tree = new TTree("tree", "tree");

  //  TApplication app("app", &iarg, argv);
  
  tree->Branch("eve", &eve, "eve/I");
  tree->Branch("adc", adc, Form("adc[%d][%d]/I",N_ADC, N_CLK_USE));
  tree->Branch("baseline", baseline, Form("baseline[%d]/F",N_ADC));  
  tree->Branch("adc_cor", adc_cor, Form("adc_cor[%d][%d]/F",N_ADC, N_CLK_USE));
  tree->Branch("integ", integ, Form("integ[%d]/F",N_ADC));
  tree->Branch("integ_gate", integ_gate, Form("integ_gate[%d]/F",N_ADC));
  tree->Branch("pulse_wid", pulse_wid, Form("pulse_wid[%d]/I",N_ADC));    
  tree->Branch("max", max, Form("max[%d]/F",N_ADC));
  tree->Branch("decay_time", decay_time, Form("decay_time[%d]/F",N_ADC));      
  
  c1 = new TCanvas();
  
  analysis(argv[1]);

  //  app.Run();
  
  tree->Write();

  for(int i=0; i<N_ADC; i++){
    h_sum[i]->Write();
  }
  
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
	  if(j>0 && j<=N_CLK_USE){
	    tmpadc = (short)htons(buf);
	    adc[i][j-1] = tmpadc;	    	    
	  }
	} // end of j loop
      } // end of i loop

      // detailed analysis
      anaevt();
      
      tree->Fill();
      eve++;

      //      int a = getchar();
      
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
    integ_gate[i]=0;
    pulse_wid[i] =0;
    max[i]=0;
    decay_time[i]=0;
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
      h_sum[i]->Fill(j, adc_cor[i][j], 1.0);
      integ[i] += adc_cor[i][j];
      if(adc_cor[i][j]>INTEG_TH){
	integ_gate[i] += adc_cor[i][j];
	pulse_wid[i]++;
      }
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
  TH1F* h_adc_cor[N_ADC];  
  float clk[N_CLK_USE];
  for(int i=0; i<N_CLK; i++) clk[i] = i;

  for(int i=0; i<N_ADC; i++){
    //    gr_adc_cor[i] = new TGraph(Form("gr %d", i), N_CLK_USE, clk, adc_cor[i]);
    gr_adc_cor[i] = new TGraph(N_CLK_USE, clk, adc_cor[i]);
    gr_adc_cor[i]->SetMarkerStyle(20);
    h_adc_cor[i] = new TH1F(Form("h%d", i), "hist", 1024, 0, 1024);

    for(int j=0; j<N_CLK_USE; j++){
      h_adc_cor[i]->SetBinContent(j+1, adc_cor[i][j]);
    }
  }

  // define fit function for 
  TF1 *func = new TF1("func", "[0]+[1]*exp(-1.0*(x-[2])/[3])", 0, 1024);

  // fit routine
  for(int i=0; i<N_ANA_ADC; i++){
    func->SetParameter(0, 0);
    func->SetParameter(1, 1000);
    func->SetParameter(2, 0);
    func->SetParameter(3, 100);    
    func->SetParLimits(0, -100, 100);
    func->SetParLimits(1, 0, 10000);
    func->SetParLimits(3, 0, 10000);    
    
    //    gr_adc_cor[i]->Fit("func", "Q", "", DECAY_FIT_MIN, DECAY_FIT_MAX);
    h_adc_cor[i]->Fit("func", "Q", "", DECAY_FIT_MIN, DECAY_FIT_MAX);    
    decay_time[i] = func->GetParameter(3);
    if(decay_time[i]>20000) decay_time[i]=0;

  for(int i=0; i<N_ADC; i++){
    gr_adc_cor[i]->Delete();
    h_adc_cor[i]->Delete();
  }
//    if(i==0){
//      printf("decay time=%f\n", decay_time[0]);
//      c1->Clear();
//      c1->cd();
//      //      gr_adc_cor[0]->Draw("AL");
//      h_adc_cor[0]->Draw("");      
//      func->Draw("LSAME");
//      //c1->Modified(); c1->Update();
//      c1->Print(Form("eve%d.pdf", eve));
//    }
    
  }
}
