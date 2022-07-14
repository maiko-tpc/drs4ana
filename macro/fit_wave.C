const int N_ADC=16;
const int N_CLK_USE=1023;

const double fit_min = 180;
const double fit_max = 600;

void fit_wave(int eve=0, int ch=0){
  printf("event=%d, ch=%d\n", eve, ch);

  float adc_cor[N_ADC][N_CLK_USE];
  tree->SetBranchAddress("adc_cor", adc_cor);
  tree->GetEntry(eve); 

  TH1F *hist = new TH1F("hist",
			Form("waveform eve=%d ch=%d", eve, ch),
			1024, 0, 1024);

  hist->GetXaxis()->SetTitle("Time (ns)");
  hist->GetYaxis()->SetTitle("Pulse height (ch)");  
  
  for(int i=0; i<N_CLK_USE; i++){
    hist->SetBinContent(i+1, adc_cor[ch][i]);
  }

  TCanvas *c = new TCanvas("c", "canvas", 600, 400);
  hist->Draw();

  TF1 *func = new TF1("func", "[0]+[1]*exp(-1.0*(x-[2])/[3])", 0, 1024);

  func->SetParameter(0, 0);
  func->SetParameter(1, 1000);
  func->SetParameter(2, 0);
  func->SetParameter(3, 100);    

  func->SetParLimits(0, -100, 100);
  func->SetParLimits(1, 0, 10000);

  
  func->SetLineWidth(3);
  
  hist->Fit("func", "", "", fit_min, fit_max);

  printf("decay constant=%f\n", func->GetParameter(3));
}
