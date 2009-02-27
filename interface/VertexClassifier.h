
#ifndef VertexClassifier_h
#define VertexClassifier_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"

#include "SimTracker/TrackHistory/interface/VertexCategories.h"
#include "SimTracker/TrackHistory/interface/VertexHistory.h"

//! Get track history and classify it in function of their .
class VertexClassifier
{

public:

    //! Constructor by ParameterSet
    VertexClassifier(edm::ParameterSet const & pset);

    //! Pre-process event information (for accessing reconstraction information)
    void newEvent(edm::Event const &, edm::EventSetup const &);

    //! Classify the RecoVertex in categories.
    VertexClassifier const & evaluate (reco::VertexRef const &);

    //! Classify the TrackingVertex in categories.
    VertexClassifier const & evaluate (TrackingVertexRef const &);

    //! Returns track flag for a given category.
    bool is(VertexCategories::Category category) const
    {
        return flags_[category];
    }

    //! Returns vertex flags with the categories description.
    const VertexCategories::Flags & flags() const
    {
        return flags_;
    }

    //! Returns a reference to the vertex history used in the classification.
    VertexHistory const & history() const
    {
        return tracer_;
    }

private:

    VertexHistory tracer_;

    const edm::InputTag hepMCLabel_;

    double longLivedDecayLength_;
    double vertexClusteringDistance_;

    struct G4
    {
        enum Process
        {
            Undefined = 0,
            Unknown,
            Primary,
            Hadronic,
            Decay,
            Compton,
            Annihilation,
            EIoni,
            HIoni,
            MuIoni,
            Photon,
            MuPairProd,
            Conversions,
            EBrem,
            SynchrotronRadiation,
            MuBrem,
            MuNucl
        };
    };

    VertexCategories::Flags flags_;

    edm::Handle<edm::HepMCProduct> mcInformation_;

    edm::ESHandle<ParticleDataTable> particleDataTable_;

    //! Reset the categories flags.
    void reset()
    {
        flags_ = VertexCategories::Flags(VertexCategories::Unknown + 1, false);
    }

    //! Get reconstruction information
    void reconstructionInformation(reco::TrackBaseRef const &);

    //! Get all the information related to the simulation details
    void simulationInformation();

    //! Get hadron flavor of the initial hadron
    // void hadronFlavor();

    //! Get all the information related to decay process
    void processesAtGenerator();

    //! Get information about conversion and other interactions
    void processesAtSimulation();

    //! Get geometrical information about the vertices
    void vertexInformation();

    // Check for unkown classification
    void unknownVertex();

    //! Auxiliary class holding simulated primary vertices
    struct GeneratedPrimaryVertex
    {
        GeneratedPrimaryVertex(double x1,double y1,double z1): x(x1), y(y1), z(z1), ptsq(0), nGenTrk(0) {}

        bool operator< ( GeneratedPrimaryVertex const & reference) const
        {
            return ptsq < reference.ptsq;
        }

        double x, y, z;
        double ptsq;
        int nGenTrk;

        HepMC::FourVector ptot;

        std::vector<int> finalstateParticles;
        std::vector<int> simTrackIndex;
        std::vector<int> genVertex;
    };

    std::vector<GeneratedPrimaryVertex> genpvs_;

    // Auxiliary function to get the generated primary vertex
    bool isFinalstateParticle(const HepMC::GenParticle *);
    bool isCharged(const HepMC::GenParticle *);
    void genPrimaryVertices();

};

// Operation overload for printing the categories
std::ostream & operator<< (std::ostream &, VertexClassifier const &);

#endif