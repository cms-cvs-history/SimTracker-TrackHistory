
{

    gROOT->SetStyle("Plain");

    TFile file("test.root");

    TCanvas * c1 = new TCanvas;

    gPad->SetLogx();

    TH1D * hist = (TH1D*) gDirectory->Get("vertexTrackHistory");
    hist->SetFillColor(kRed);
    hist->Draw("hbar1");

    TH1D * histplus = (TH1D*) gDirectory->Get("vertexTrackHistory_plus_error");
    histplus->SetFillColor(kBlack);
    histplus->SetFillStyle(3005);
    histplus->Draw("hbar same");

    TH1D * histminus = (TH1D*) gDirectory->Get("vertexTrackHistory_minus_error");
    histminus->SetFillColor(kRed);
    histminus->Draw("hbar1 same");

}
