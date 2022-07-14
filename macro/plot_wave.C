const int N_ADC=16;
const int N_CLK_USE=1023;

void plot_wave(int eve=0, int ch=0){
  printf("event=%d, ch=%d\n", eve, ch);

  float adc_cor[N_ADC][N_CLK_USE];
  tree->SetBranchAddress("adc_cor", adc_cor);
  tree->GetEntry(eve); 

  TH1F *hist = new TH1F("hist",
			Form("waveform eve=%d ch=%d", eve, ch),
			N_CLK_USE, 0, N_CLK_USE);
  
  hist->GetXaxis()->SetTitle("Clock");
  hist->GetYaxis()->SetTitle("ADC");  

  for(int i=0; i<N_CLK_USE; i++){
    hist->SetBinContent(i, adc_cor[ch][i]);
  }

  TCanvas *c = new TCanvas("c", "canvas", 600, 400);
  hist->Draw();
}
