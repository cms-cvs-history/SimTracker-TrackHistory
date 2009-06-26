#include "FWCore/Framework/interface/ESHandle.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "SimTracker/TrackHistory/test/testVertexAssociationAnalyzer.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimTracker/Records/interface/VertexAssociatorRecord.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorBase.h"
#include "SimTracker/VertexAssociation/interface/VertexAssociatorBase.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h" 

#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "DataFormats/Math/interface/Vector.h"

#include <memory>
#include <iostream>
#include <string>

class TrackAssociatorByHits;
class TrackerHitAssociator;

using namespace reco;
using namespace std;
using namespace edm;


testVertexAssociationAnalyzer::testVertexAssociationAnalyzer(edm::ParameterSet const& conf) : theConfig_(conf) {

  vertexCollection_ = conf.getUntrackedParameter<edm::InputTag>( "vertexCollection" );

  n_event_ = 0;
  n_rs_vertices_ = 0;
  n_rs_vtxassocs_ = 0;
  n_sr_vertices_ = 0;
  n_sr_vtxassocs_ = 0;

}

testVertexAssociationAnalyzer::~testVertexAssociationAnalyzer() {


}

void testVertexAssociationAnalyzer::beginJob(const EventSetup & setup) {

  edm::ESHandle<MagneticField> theMF;
  setup.get<IdealMagneticFieldRecord>().get(theMF);

  edm::ESHandle<TrackAssociatorBase> theChiAssociator;
  setup.get<TrackAssociatorRecord>().get("TrackAssociatorByChi2",theChiAssociator);
  associatorByChi2 = (TrackAssociatorBase *) theChiAssociator.product();

  edm::ESHandle<TrackAssociatorBase> theHitsAssociator;
  setup.get<TrackAssociatorRecord>().get("TrackAssociatorByHits",theHitsAssociator);
  associatorByHits = (TrackAssociatorBase *) theHitsAssociator.product();

  edm::ESHandle<VertexAssociatorBase> theTracksAssociator;
  setup.get<VertexAssociatorRecord>().get("VertexAssociatorByTracks",theTracksAssociator);
  associatorByTracks = (VertexAssociatorBase *) theTracksAssociator.product();

   edm::Service<TFileService> fs; 

   //RecoToSim Histos

  rs_nvtx = fs->make<TH1F>("rs_nvtx","# of RecoToSim vertices per event",11,-0.5,10.5);
  rs_dist = fs->make<TH1F>("rs_dist","r Miss Distance (cm)",100,0,0.1);
  rs_recz = fs->make<TH1F>("rs_recz","z, Reconstructed Vertex (cm)", 200, -25.0,25.0);
  rs_simz = fs->make<TH1F>("rs_simz","z, Simulated Vertex (cm)",     200, -25.0,25.0);
  rs_nsimtrk = fs->make<TH1F>("rs_nsimtrk","# of tracks, Simulated",    501,-0.5,500.5);  
  rs_nrectrk = fs->make<TH1F>("rs_nrectrk","# of tracks, Reconstructed",501,-0.5,500.5);
  rs_qual = fs->make<TH1F>("rs_qual","Quality of Match",51,-0.01,1.01);
  rs_chi2norm = fs->make<TH1F>("rs_chi2norm","Vertex Normalized Chi2",100, 0, 10.);
  rs_chi2prob = fs->make<TH1F>("rs_chi2prob","Vertex Chi2 Probability",100, 0, 1.);
  rs_resx = fs->make<TH1F>("rs_resx","rec - sim Distance (cm)",100,-0.05,0.05);  
  rs_resy = fs->make<TH1F>("rs_resy","rec - sim Distance (cm)",100,-0.05,0.05);  
  rs_resz = fs->make<TH1F>("rs_resz","rec - sim Distance (cm)",100,-0.05,0.05);  
  rs_pullx = fs->make<TH1F>("rs_pullx","(rec - sim)/err_rec ",100,-10.,10.);  
  rs_pully = fs->make<TH1F>("rs_pully","(rec - sim)/err_rec ",100,-10.,10.); 
  rs_pullz = fs->make<TH1F>("rs_pullz","(rec - sim)/err_rec ",100,-10.,10.); 
  rs_pt = fs->make<TH1F>("rs_pt","Vertex Transverse Momentum (GeV)",2000, 0, 1000.);
  rs_eta = fs->make<TH1F>("rs_eta","Vertex Eta", 200, -3., 3.);
  rs_mass = fs->make<TH1F>("rs_mass","Vertex Mass",65, 0, 6.5);
  rs_charge = fs->make<TH1F>("rs_charge","charge of reconstructed vertex",201,-100.5,100.5);
  rs_averageweight = fs->make<TH1F>("rs_averageweight","Average sum of track weights in vertex", 100, -0.1, 1.1);
  rs_simpt = fs->make<TH1F>("rs_simpt","Vertex Transverse Momentum (GeV)",2000, 0, 1000.);
  rs_simeta = fs->make<TH1F>("rs_simeta","Vertex Eta", 200, -3., 3.);
  rs_simcharge = fs->make<TH1F>("rs_simcharge","charge of reconstructed vertex",201,-100.5,100.5);

  //SimToReco Histos

  sr_nvtx = fs->make<TH1F>("sr_nvtx","# of SimToReco vertices per event",11,-0.5,10.5);
  sr_dist = fs->make<TH1F>("sr_dist","r Miss Distance (cm)",100,0,0.1);
  sr_recz = fs->make<TH1F>("sr_recz","z, Reconstructed Vertex (cm)", 200, -25.0,25.0);
  sr_simz = fs->make<TH1F>("sr_simz","z, Simulated Vertex (cm)",     200, -25.0,25.0);
  sr_nsimtrk = fs->make<TH1F>("sr_nsimtrk","# of tracks, Simulated",    501,-0.5,500.5);  
  sr_nrectrk = fs->make<TH1F>("sr_nrectrk","# of tracks, Reconstructed",501,-0.5,500.5);
  sr_qual = fs->make<TH1F>("sr_qual","Quality of Match",51,-0.01,1.01);
  sr_chi2norm = fs->make<TH1F>("sr_chi2norm","Vertex Normalized Chi2",100, 0, 10.);
  sr_chi2prob = fs->make<TH1F>("sr_chi2prob","Vertex Chi2 Probability",100, 0, 1.);
  sr_resx = fs->make<TH1F>("sr_resx","rec - sim Distance (cm)",100,-0.05,0.05);  
  sr_resy = fs->make<TH1F>("sr_resy","rec - sim Distance (cm)",100,-0.05,0.05);  
  sr_resz = fs->make<TH1F>("sr_resz","rec - sim Distance (cm)",100,-0.05,0.05);  
  sr_pullx = fs->make<TH1F>("sr_pullx","(rec - sim)/err_rec ",100,-10.,10.);  
  sr_pully = fs->make<TH1F>("sr_pully","(rec - sim)/err_rec ",100,-10.,10.); 
  sr_pullz = fs->make<TH1F>("sr_pullz","(rec - sim)/err_rec ",100,-10.,10.); 
  sr_pt = fs->make<TH1F>("sr_pt","Vertex Transverse Momentum (GeV)",2000, 0, 1000.);
  sr_eta = fs->make<TH1F>("sr_eta","Vertex Eta", 200, -3., 3.);
  sr_mass = fs->make<TH1F>("sr_mass","Vertex Mass",65, 0, 6.5);
  sr_charge = fs->make<TH1F>("sr_charge","charge of reconstructed vertex",201,-100.5,100.5);
  sr_averageweight = fs->make<TH1F>("sr_averageweight","Average sum of track weights in vertex", 100, -0.1, 1.1);
  sr_simpt = fs->make<TH1F>("sr_simpt","Vertex Transverse Momentum (GeV)",2000, 0, 1000.);
  sr_simeta = fs->make<TH1F>("sr_simeta","Vertex Eta", 200, -3., 3.);
  sr_simcharge = fs->make<TH1F>("sr_simcharge","charge of reconstructed vertex",201,-100.5,100.5);


}

void testVertexAssociationAnalyzer::endJob() {
 
  std::cout << std::endl;
  std::cout << " ====== Total Number of analyzed events: " << n_event_        << " ====== " << std::endl;
  std::cout << " ====== Total Number of R2S vertices:    " << n_rs_vertices_  << " ====== " << std::endl;
  std::cout << " ====== Total Number of R2S vtx assocs:  " << n_rs_vtxassocs_ << " ====== " << std::endl;
  std::cout << " ====== Total Number of S2R vertices:    " << n_sr_vertices_  << " ====== " << std::endl;
  std::cout << " ====== Total Number of S2R vtx assocs:  " << n_sr_vtxassocs_ << " ====== " << std::endl;
}

void testVertexAssociationAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& setup)
{

  using namespace edm;
  using namespace reco;

  ++n_event_;

  std::cout << "*** Analyzing " << event.id() << " n_event = " << n_event_ << std::endl << std::endl;


  Handle<edm::View<reco::Track> > trackCollection;
  event.getByLabel("generalTracks",trackCollection);
  const  edm::View<reco::Track>   tC = *(trackCollection.product());

  edm::Handle<TrackingParticleCollection>  TPCollection ;
  event.getByLabel("mergedtruth","MergedTrackTruth",TPCollection);
  const TrackingParticleCollection tPC   = *(TPCollection.product());

  edm::Handle<TrackingVertexCollection>  TVCollection;
  event.getByLabel("mergedtruth","MergedTrackTruth",TVCollection);
  const TrackingVertexCollection tVC   = *(TVCollection.product());

  //Vertex Collection
  edm::Handle<reco::VertexCollection>  vertexCollection;
  event.getByLabel(vertexCollection_, vertexCollection);
  const reco::VertexCollection vC   = *(vertexCollection.product());

  cout << endl;
  cout << "                      ****************** Before Assocs ****************** " << endl << endl;  

  cout << "trackCollection.size()  = " << tC.size() << endl;
  cout << "TPCollection.size()     = " << tPC.size() << endl;
  cout << "vertexCollection.size() = " << vC.size() << endl;
  cout << "TVCollection.size()     = " << tVC.size() << endl;
  
  cout << endl;
  cout << "                      ****************** Reco To Sim ****************** " << endl << endl;

  //cout << "-- Associator by hits --" << endl;
  reco::RecoToSimCollection r2sTracks = associatorByHits->associateRecoToSim (trackCollection,TPCollection,&event );

  reco::SimToRecoCollection s2rTracks = associatorByHits->associateSimToReco (trackCollection,TPCollection,&event );
  //associatorByChi2->associateRecoToSim (trackCollection,TPCollection,&event );

  reco::VertexRecoToSimCollection r2sVertices = associatorByTracks->associateRecoToSim(vertexCollection,TVCollection,event,r2sTracks);

  reco::VertexSimToRecoCollection s2rVertices = associatorByTracks->associateSimToReco(vertexCollection,TVCollection,event,s2rTracks);

  cout << endl;
  cout << "associateRecoToSim tracks size = " <<  r2sTracks.size() << " ; associateSimToReco tracks size = " << s2rTracks.size() << endl;
  cout << "VertexRecoToSim size           = " <<  r2sVertices.size() << " ; VertexSimToReco size           = " << r2sVertices.size() << " " << endl;
 
 
  cout << endl << " [testVertexAssociationAnalyzer] Analyzing Reco To Sim" << endl;

  double thePiMass = 0.13957;


  rs_nvtx->Fill( r2sVertices.size() );

  int cont_recvR2S = 0;
  for (reco::VertexRecoToSimCollection::const_iterator iR2S = r2sVertices.begin(); iR2S != r2sVertices.end(); ++iR2S, ++cont_recvR2S) {

    ++n_rs_vertices_;

    reco::VertexRef recVertex = iR2S->key;
    math::XYZPoint recPos = recVertex->position();

    double nrectrk = recVertex->tracksSize();

    double sum_weight = 0.;
    double weight = 0.;
    math::XYZVector momentum;
    math::XYZTLorentzVector sum;
    int charge = 0;
    for ( reco::Vertex::trackRef_iterator recDaughter = recVertex->tracks_begin() ; recDaughter != recVertex->tracks_end(); ++recDaughter) {

      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double> > vec;

      vec.SetPx( (**recDaughter).px() );
      vec.SetPy( (**recDaughter).py() );
      vec.SetPz( (**recDaughter).pz() );
      vec.SetM(thePiMass); 

      sum += vec;
  
      weight = recVertex->trackWeight(*recDaughter); 
      sum_weight += weight;

      math::XYZVector p;
      p = (**recDaughter).momentum();
      momentum += p; 

      charge += (*recDaughter)->charge();

    }// end loop recDaughter

    double mass = sum.M();
    double pt = sqrt( momentum.Perp2() );
    double eta = momentum.Eta();

    std::vector<std::pair<TrackingVertexRef, double> > simVertices = iR2S->val;

    int cont_simvR2S = 0;
    for (std::vector<std::pair<TrackingVertexRef, double> >::const_iterator iMatch = simVertices.begin(); iMatch != simVertices.end(); ++iMatch, ++cont_simvR2S) {

      TrackingVertexRef simVertex =  iMatch->first;
      math::XYZTLorentzVectorD simVec = (iMatch->first)->position();
      math::XYZPoint simPos = math::XYZPoint(simVec.x(),simVec.y(),simVec.z());

      ++n_rs_vtxassocs_;

      cout << "rec vertex " << cont_recvR2S << " has associated sim vertex " << cont_simvR2S << endl;
 
      double nsimtrk = simVertex->daughterTracks().size();
      double qual  = iMatch->second;

      double chi2norm = recVertex->normalizedChi2();
      double chi2prob = ChiSquaredProbability( recVertex->chi2(), recVertex->ndof() ); 
        
      double resx = recVertex->x() - simVertex->position().x();
      double resy = recVertex->y() - simVertex->position().y();
      double resz = recVertex->z() - simVertex->position().z();	
      double pullx = (recVertex->x() - simVertex->position().x())/recVertex->xError();
      double pully = (recVertex->y() - simVertex->position().y())/recVertex->yError();
      double pullz = (recVertex->z() - simVertex->position().z())/recVertex->zError();
      double dist = sqrt(resx*resx+resy*resy+resz*resz);

      cout << "            R2S: recPos = " << recPos << " ; simPos = " << simPos << endl;

      math::XYZVector simmomentum;
      int simcharge = 0;
      for ( TrackingVertex::tp_iterator simDaughter = simVertex->daughterTracks_begin(); simDaughter != simVertex->daughterTracks_end(); ++simDaughter ){

	math::XYZVector p;
	p = (**simDaughter).momentum();
	simmomentum += p; 

	simcharge += (*simDaughter)->charge();

      }//end loop simDaughter

      double simpt = sqrt( simmomentum.Perp2() );
      double simeta = simmomentum.Eta();


      rs_resx->Fill(resx);
      rs_resy->Fill(resy);
      rs_resz->Fill(resz);
      rs_pullx->Fill(pullx);
      rs_pully->Fill(pully);
      rs_pullz->Fill(pullz);
      rs_dist->Fill(dist);
      rs_simz->Fill(simPos.Z());
      rs_recz->Fill(recPos.Z());
      rs_nsimtrk->Fill(nsimtrk);
      rs_nrectrk->Fill(nrectrk);
      rs_qual->Fill(qual);
      rs_chi2norm->Fill(chi2norm);
      rs_chi2prob->Fill(chi2prob);
      rs_pt->Fill(pt);
      rs_eta->Fill(eta);
      rs_charge->Fill(charge);
      rs_mass->Fill(mass);
      rs_averageweight->Fill(sum_weight/nrectrk);
      rs_simpt->Fill(simpt);
      rs_simeta->Fill(simeta);
      rs_simcharge->Fill(simcharge);

    }//end simVertices
    
  }//end iR2S


  cout << endl << "                      ****************** Sim To Reco ****************** " << endl << endl;

  cout << endl << "[testVertexAssociationAnalyzer] Analyzing Sim To Reco" << endl;

  sr_nvtx->Fill( s2rVertices.size() );

  int cont_simvS2R = 0;
  for (reco::VertexSimToRecoCollection::const_iterator iS2R = s2rVertices.begin(); iS2R != s2rVertices.end(); ++iS2R, ++cont_simvS2R) {

    ++n_sr_vertices_;

    TrackingVertexRef simVertex = iS2R->key;
    math::XYZTLorentzVectorD simVec = simVertex->position();
    math::XYZPoint   simPos = math::XYZPoint(simVec.x(),simVec.y(),simVec.z());

    double nsimtrk = simVertex->daughterTracks().size();

    std::vector<std::pair<VertexRef, double> > recoVertices = iS2R->val;

    math::XYZVector simmomentum;
    int simcharge = 0;
    for ( TrackingVertex::tp_iterator simDaughter = simVertex->daughterTracks_begin(); simDaughter != simVertex->daughterTracks_end(); simDaughter++ ) {

      math::XYZVector p;
      p = (**simDaughter).momentum();
      simmomentum += p; 

      simcharge += (*simDaughter)->charge();

    }//end loop a simDaughters

    double simpt = sqrt( simmomentum.Perp2() );
    double simeta = simmomentum.Eta();

    int cont_recvS2R = 0;
    for (std::vector<std::pair<VertexRef, double> >::const_iterator iMatch = recoVertices.begin(); iMatch != recoVertices.end(); ++iMatch, ++cont_recvS2R) {

      VertexRef recVertex = iMatch->first;
      math::XYZPoint recPos = recVertex->position();

      ++n_sr_vtxassocs_;

      cout << "sim vertex " << cont_simvS2R << " has associated rec vertex " << cont_recvS2R << endl;

      double nrectrk = recVertex->tracksSize();
      double qual  = iMatch->second;

      double chi2norm = recVertex->normalizedChi2();
      double chi2prob = ChiSquaredProbability( recVertex->chi2(), recVertex->ndof() ); 
         
      double resx = recVertex->x() - simVertex->position().x();
      double resy = recVertex->y() - simVertex->position().y(); 
      double resz = recVertex->z() - simVertex->position().z(); 	
      double pullx = ( recVertex->x() - simVertex->position().x() )/recVertex->xError();
      double pully = ( recVertex->y() - simVertex->position().y() )/recVertex->yError();
      double pullz = ( recVertex->z() - simVertex->position().z() )/recVertex->zError();
      double dist = sqrt(resx*resx+resy*resy+resz*resz);

      cout << "            S2R: simPos = " << simPos << " ; recPos = " << recPos << endl;

      double sum_weight = 0.;
      double weight = 0.;
      math::XYZVector momentum;
      math::XYZTLorentzVector sum;
      int charge = 0;
      for ( reco::Vertex::trackRef_iterator recDaughter = recVertex->tracks_begin() ; recDaughter != recVertex->tracks_end(); ++recDaughter ) {

	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double> > vec;

	vec.SetPx( (**recDaughter).px() );
	vec.SetPy( (**recDaughter).py() );
	vec.SetPz( (**recDaughter).pz() );
	vec.SetM(thePiMass); 

	sum += vec;
  
	weight = recVertex->trackWeight(*recDaughter); 
	sum_weight += weight;

	math::XYZVector p;
	p = (**recDaughter).momentum();
	momentum += p; 

	charge += (*recDaughter)->charge();

      }//end loop recDaughter  

      double mass = sum.M();
      double pt = sqrt( momentum.Perp2() );
      double eta = momentum.Eta();


      sr_resx->Fill(resx);
      sr_resy->Fill(resy);
      sr_resz->Fill(resz);
      sr_pullx->Fill(pullx);
      sr_pully->Fill(pully);
      sr_pullz->Fill(pullz);
      sr_dist->Fill(dist);
      sr_simz->Fill(simPos.Z());
      sr_recz->Fill(recPos.Z());
      sr_nsimtrk->Fill(nsimtrk);
      sr_nrectrk->Fill(nrectrk);
      sr_qual->Fill(qual);
      sr_chi2norm->Fill(chi2norm);
      sr_chi2prob->Fill(chi2prob);	
      sr_pt->Fill(pt);
      sr_eta->Fill(eta);
      sr_charge->Fill(charge);
      sr_mass->Fill(mass);
      sr_averageweight->Fill(sum_weight/nrectrk);
      sr_simpt->Fill(simpt);
      sr_simeta->Fill(simeta);
      sr_simcharge->Fill(simcharge);  

    }//end recoVertices

  }//end iS2R

  cout << endl;

}

#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(testVertexAssociationAnalyzer);

