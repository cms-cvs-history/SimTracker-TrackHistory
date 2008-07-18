
{

    TFile file("test.root");

    TCanvas c1;
    TPie * pie = (TPie*) gDirectory->Get("recoTrackCollection");
    pie->Draw();

    TCanvas c2;
    pie = (TPie*) gDirectory->Get("recoTrackCollection_pi+");
    pie->Draw();

}
