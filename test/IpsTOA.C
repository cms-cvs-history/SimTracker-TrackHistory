{

   gROOT->SetStyle("Plain");

   TFile file("test.root");

   TCanvas * c1 = new TCanvas;

   TH1D * hist = (TH1D*) gDirectory->Get("ips3d_0_35");
   hist->SetFillColor(49);
   hist->Draw("hbar2");

   TCanvas * c2 = new TCanvas;
   hist = (TH1D*) gDirectory->Get("ips3d_4_35");
   hist->SetFillColor(49);
   hist->Draw("hbar2"); 

}
