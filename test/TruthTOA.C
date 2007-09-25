{

   gROOT->SetStyle("Plain");

   TFile file("test.root");

   TCanvas * c1 = new TCanvas;

   TH1D * hist = (TH1D*) gDirectory->Get("TrackingParticleCollection");
   hist->SetFillColor(49);
   hist->Draw("hbar2");

   TCanvas * c2 = new TCanvas;
   hist = (TH1D*) gDirectory->Get("TrackingParticleCollection_pi+");
   hist->SetFillColor(49);
   hist->Draw("hbar2"); 

}
