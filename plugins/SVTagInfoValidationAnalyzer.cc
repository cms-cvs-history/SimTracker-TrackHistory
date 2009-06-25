
#include <map>
#include <string>

#include "TH1F.h"

#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

#include "SimTracker/TrackHistory/interface/VertexClassifierByProxy.h"

//
// class decleration
//

class SVTagInfoValidationAnalyzer : public edm::EDAnalyzer
{

public:

    explicit SVTagInfoValidationAnalyzer(const edm::ParameterSet&);

private:

    virtual void analyze(const edm::Event&, const edm::EventSetup&);

    // Member data

    VertexClassifierByProxy<reco::SecondaryVertexTagInfoCollection> classifier_;

    Int_t numberVertexClassifier_;

    edm::InputTag trackingTruth_;

    edm::InputTag svTagInfoProducer_;

    // Get the file service
    edm::Service<TFileService> fs_;

    // Bookeeping of all the histograms per category
    void bookRecoToSim(std::string const &);
    void bookSimToReco(std::string const &);

    // Fill all histogram per category
    void fillRecoToSim(std::string const &, reco::Vertex const &, TrackingVertexRef const &);
    void fillSimToReco(std::string const &, reco::VertexRef const &, TrackingVertexRef const &);

    // Histogram handlers
    std::map<std::string, TH1D *> TH1Index_;

};


SVTagInfoValidationAnalyzer::SVTagInfoValidationAnalyzer(const edm::ParameterSet& config) : classifier_(config)
{
    // Get the track collection
    svTagInfoProducer_ = config.getUntrackedParameter<edm::InputTag> ( "svTagInfoProducer" );

    // Name of the traking pariticle collection
    trackingTruth_ = config.getUntrackedParameter<edm::InputTag> ( "trackingTruth" );

    // Number of track categories
    numberVertexClassifier_ = VertexCategories::Unknown+1;

    // Define histogram for counting categories
    TH1Index_["VertexClassifier"] = fs_->make<TH1D>(
                                        "VertexClassifier",
                                        "Frequency for the different track categories",
                                        numberVertexClassifier_,
                                        -0.5,
                                        numberVertexClassifier_ - 0.5
                                    );

    // Set the proper categories names
    for (Int_t i = 0; i < numberVertexClassifier_; ++i)
        TH1Index_["VertexClassifier"]->GetXaxis()->SetBinLabel(i+1, VertexCategories::Names[i]);

    // book histograms
    bookRecoToSim("All");
    bookRecoToSim("BWeakDecay");
    bookRecoToSim("CWeakDecay");
    bookRecoToSim("Light");
    bookSimToReco("ReconstructedTV");    
}


void SVTagInfoValidationAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& setup)
{
    // Set the classifier for a new event
    classifier_.newEvent(event, setup);

    // Vertex collection
    edm::Handle<reco::SecondaryVertexTagInfoCollection> svTagInfoCollection;
    event.getByLabel(svTagInfoProducer_, svTagInfoCollection);

    // Get a constant reference to the track history associated to the classifier
    VertexHistory const & tracer = classifier_.history();

    // Loop over the svTagInfo collection.
    for (std::size_t index = 0; index < svTagInfoCollection->size(); ++index)
    {
        reco::SecondaryVertexTagInfoRef svTagInfo(svTagInfoCollection, index);

        // Loop over the vertexes in svTagInfo
        for (std::size_t vindex = 0; vindex < svTagInfo->nVertices(); ++vindex)
        {

            // Classify the tracks
            classifier_.evaluate(svTagInfo, vindex);

            // Fill the histogram with the categories
            for (Int_t i = 0; i != numberVertexClassifier_; ++i)
                if (
                    classifier_.is( (VertexCategories::Category) i )
                )
                    TH1Index_["VertexClassifier"]->Fill(i);

            if ( !classifier_.is(VertexCategories::Fake) )
            {
                // Fill histograms
                fillRecoToSim("All", svTagInfo->secondaryVertex(vindex), tracer.simVertex());
                if ( classifier_.is(VertexCategories::BWeakDecay) )
                    fillRecoToSim("BWeakDecay", svTagInfo->secondaryVertex(vindex), tracer.simVertex());
                else if ( classifier_.is(VertexCategories::CWeakDecay) )
                    fillRecoToSim("CWeakDecay", svTagInfo->secondaryVertex(vindex), tracer.simVertex());
                else
                    fillRecoToSim("Light", svTagInfo->secondaryVertex(vindex), tracer.simVertex());
            }

        }
    }
    
    // Vertex collection
    edm::Handle<TrackingVertexCollection> trackingTruth;
    event.getByLabel(trackingTruth_, trackingTruth);

    // Tracking particle information
    edm::Handle<TrackingVertexCollection>  TVCollection;
    event.getByLabel(trackingTruth_, TVCollection);
    
    // Loop over the TV collection.
    for (std::size_t index = 0; index < TVCollection->size(); ++index)
    {        
    	TrackingVertexRef trackingVertex(TVCollection, index);

    	classifier_.evaluate(trackingVertex);
    	
    	if ( classifier_.is(VertexCategories::Reconstructed) )
            fillSimToReco("ReconstructedTV", tracer.recoVertex(), trackingVertex);
    }
    
}


void SVTagInfoValidationAnalyzer::bookRecoToSim(std::string const & prefix)
{
    // Book pull histograms

    std::string name = prefix + "VertexPullx";
    TH1Index_[name] = fs_->make<TH1D>(name.c_str(), name.c_str(), 50, -10., 10.);
    name = prefix + "VertexPully";
    TH1Index_[name] = fs_->make<TH1D>(name.c_str(), name.c_str(), 50, -10., 10.);
    name = prefix + "VertexPullz";
    TH1Index_[name] = fs_->make<TH1D>(name.c_str(), name.c_str(), 50, -10., 10.);
}


void SVTagInfoValidationAnalyzer::bookSimToReco(std::string const & prefix)
{
    // Book pull histograms

    std::string name = prefix + "VertexPullx";
    TH1Index_[name] = fs_->make<TH1D>(name.c_str(), name.c_str(), 50, -10., 10.);
    name = prefix + "VertexPully";
    TH1Index_[name] = fs_->make<TH1D>(name.c_str(), name.c_str(), 50, -10., 10.);
    name = prefix + "VertexPullz";
    TH1Index_[name] = fs_->make<TH1D>(name.c_str(), name.c_str(), 50, -10., 10.);
}


void SVTagInfoValidationAnalyzer::fillRecoToSim(std::string const & prefix, reco::Vertex const & vertex, TrackingVertexRef const & simVertex)
{
    // Histograms for efficiency
    double pullx = (vertex.x() - simVertex->position().x())/vertex.xError();
    double pully = (vertex.y() - simVertex->position().y())/vertex.yError();
    double pullz = (vertex.z() - simVertex->position().z())/vertex.zError();

    TH1Index_[prefix + "VertexPullx"]->Fill(pullx);
    TH1Index_[prefix + "VertexPully"]->Fill(pully);
    TH1Index_[prefix + "VertexPullz"]->Fill(pullz);
}


void SVTagInfoValidationAnalyzer::fillSimToReco(std::string const & prefix, reco::VertexRef const & vertex, TrackingVertexRef const & simVertex)
{
    // Fill pull histograms

    double pullx = (vertex->x() - simVertex->position().x())/vertex->xError();
    double pully = (vertex->y() - simVertex->position().y())/vertex->yError();
    double pullz = (vertex->z() - simVertex->position().z())/vertex->zError();

    TH1Index_[prefix + "VertexPullx"]->Fill(pullx);
    TH1Index_[prefix + "VertexPully"]->Fill(pully);
    TH1Index_[prefix + "VertexPullz"]->Fill(pullz);
}

DEFINE_ANOTHER_FWK_MODULE(SVTagInfoValidationAnalyzer);