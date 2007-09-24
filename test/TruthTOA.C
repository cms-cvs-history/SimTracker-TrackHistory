{

TFile file("test.root");

TCanvas c1;
TPie * pie = (TPie*) gDirectory->Get("TrackingParticleCollection");
pie->Draw();

TCanvas c2;
pie = (TPie*) gDirectory->Get("TrackingParticleCollection_pi+");
pie->Draw();

}
