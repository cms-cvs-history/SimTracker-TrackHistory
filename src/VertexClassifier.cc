/*
 *  VertexClassifier.C
 */

#include <math.h>
#include <cstdlib>
#include <iostream>

#include "HepPDT/ParticleID.hh"

#include "SimTracker/TrackHistory/interface/VertexClassifier.h"


#define update(a, b) do { (a) = (a) | (b); } while(0)


VertexClassifier::VertexClassifier(edm::ParameterSet const & config) : VertexCategories(), tracer_(config),
        hepMCLabel_( config.getUntrackedParameter<edm::InputTag>("hepMC") )
{
    // Set the history depth after hadronization
    tracer_.depth(-2);

    // Set the minimum decay length for detecting long decays
    longLivedDecayLength_ = config.getUntrackedParameter<double>("longLivedDecayLength");

    // Set the distance for clustering vertices
    vertexClusteringDistance_ = config.getUntrackedParameter<double>("vertexClusteringDistance");
}


void VertexClassifier::newEvent(edm::Event const & event, edm::EventSetup const & setup)
{
    // Get the new event information for the tracer
    tracer_.newEvent(event, setup);

    // Get hepmc of the event
    event.getByLabel(hepMCLabel_, mcInformation_);

    // Get the partivle data table
    setup.getData(particleDataTable_);

    // Create the list of primary vertices associated to the event
    genPrimaryVertices();
}


VertexClassifier const & VertexClassifier::evaluate (reco::VertexRef const & vertex)
{
    // Initializing the category vector
    reset();

    // Associate and evaluate the vertex history (check for fakes)
    if ( tracer_.evaluate(vertex) )
    {
        // Get all the information related to the simulation details
        simulationInformation();

        // Get all the information related to decay process
        processesAtGenerator();

        // Get information about conversion and other interactions
        processesAtSimulation();

        // Get geometrical information about the vertices
        vertexInformation();

        // Check for unkown classification
        unknownVertex();
    }
    else
        flags_[Fake] = true;

    return *this;
}


VertexClassifier const & VertexClassifier::evaluate (TrackingVertexRef const & vertex)
{
    // Initializing the category vector
    reset();

    // Trace the history for the given TP
    tracer_.evaluate(vertex);

    // Check for a reconstructed track
    if ( tracer_.recoVertex().isNonnull() )
        flags_[Reconstructed] = true;
    else
        flags_[Reconstructed] = false;

    // Get all the information related to the simulation details
    simulationInformation();

    // Get all the information related to decay process
    processesAtGenerator();

    // Get information about conversion and other interactions
    processesAtSimulation();

    // Get geometrical information about the vertices
    vertexInformation();

    // Check for unkown classification
    unknownVertex();

    return *this;
}


void VertexClassifier::simulationInformation()
{
    // Get the event id for the initial TP.
    EncodedEventId eventId = tracer_.simVertex()->eventId();
    // Check for signal events
    flags_[SignalEvent] = !eventId.bunchCrossing() && !eventId.event();
}


void VertexClassifier::processesAtGenerator()
{
    // Get the generated vetices from track history
    VertexHistory::GenVertexTrail const & genVertexTrail = tracer_.genVertexTrail();

    // Loop over the generated vertices
    for (
        VertexHistory::GenVertexTrail::const_iterator ivertex = genVertexTrail.begin();
        ivertex != genVertexTrail.end();
        ++ivertex
    )
    {
        // Get the pointer to the vertex by removing the const-ness (no const methos in HepMC::GenVertex)
        HepMC::GenVertex * vertex = const_cast<HepMC::GenVertex *>(*ivertex);

        // Loop over the sources looking for specific decays
        for (
            HepMC::GenVertex::particle_iterator iparent = vertex->particles_begin(HepMC::parents);
            iparent != vertex->particles_end(HepMC::parents);
            ++iparent
        )
        {
            // Collect the pdgid of the parent
            int pdgid = std::abs((*iparent)->pdg_id());
            // Get particle type
            HepPDT::ParticleID particleID(pdgid);

            // Check if the particle type is valid one
            if (particleID.isValid())
            {
                // Get particle data
                ParticleData const * particleData = particleDataTable_->particle(particleID);
                // Check if the particle exist in the table
                if (particleData)
                {
                    // Check if their life time is bigger than longLivedDecayLength_
                    if ( particleData->lifetime() > longLivedDecayLength_ )
                    {
                        // Check for B, C weak decays and long lived decays
                        update(flags_[BWeakDecay], particleID.hasBottom());
                        update(flags_[CWeakDecay], particleID.hasCharm());
                        update(flags_[LongLivedDecay], true);
                    }
                    // Check Tau, Ks and Lambda decay
                    update(flags_[TauDecay], pdgid == 15);
                    update(flags_[KsDecay], pdgid == 310);
                    update(flags_[LambdaDecay], pdgid == 3122);
                    update(flags_[JpsiDecay], pdgid == 443);
                    update(flags_[XiDecay], pdgid == 3312);
                    update(flags_[OmegaDecay], pdgid == 3334);
                    update(flags_[SigmaPlusDecay], pdgid == 3222);
                    update(flags_[SigmaMinusDecay], pdgid == 3112);
                }
            }
        }
    }
}


void VertexClassifier::processesAtSimulation()
{
    VertexHistory::SimVertexTrail const & simVertexTrail = tracer_.simVertexTrail();

    for (
        VertexHistory::SimVertexTrail::const_iterator ivertex = simVertexTrail.begin();
        ivertex != simVertexTrail.end();
        ++ivertex
    )
    {
        // pdgid of the real source parent vertex
        int pdgid = 0;

        // select the original source in case of combined vertices
        bool flag = false;
        TrackingVertex::tp_iterator itd, its;

        for (its = (*ivertex)->sourceTracks_begin(); its != (*ivertex)->sourceTracks_end(); ++its)
        {
            for (itd = (*ivertex)->daughterTracks_begin(); itd != (*ivertex)->daughterTracks_end(); ++itd)
                if (itd != its)
                {
                    flag = true;
                    break;
                }
            if (flag)
                break;
        }
        // Collect the pdgid of the original source track
        if ( its != (*ivertex)->sourceTracks_end() )
            pdgid = std::abs((*its)->pdgId());
        else
            pdgid = 0;

        // Loop over the simulated particles
        for (
            TrackingVertex::tp_iterator iparticle = (*ivertex)->daughterTracks_begin();
            iparticle != (*ivertex)->daughterTracks_end();
            ++iparticle
        )
        {

            if ( (*iparticle)->matchedHit() )
            {
                // Collect the G4 process of the first psimhit (it should be the same for all of them)
                unsigned short process = (*iparticle)->pSimHit_begin()->processType();

                // Look for conversion process
                flags_[Conversion] = (process == G4::Conversions);

                // Special treatment for decays
                if (process == G4::Decay)
                {
                    // Get particle type
                    HepPDT::ParticleID particleID(pdgid);
                    // Check if the particle type is valid one
                    if (particleID.isValid())
                    {
                        // Get particle data
                        ParticleData const * particleData = particleDataTable_->particle(particleID);
                        // Check if the particle exist in the table
                        if (particleData)
                        {
                            // Check if their life time is bigger than 1e-14
                            if ( particleDataTable_->particle(particleID)->lifetime() > longLivedDecayLength_ )
                            {
                                // Check for B, C weak decays and long lived decays
                                update(flags_[BWeakDecay], particleID.hasBottom());
                                update(flags_[CWeakDecay], particleID.hasCharm());
                                update(flags_[LongLivedDecay], true);
                           }
                            // Check Tau, Ks and Lambda decay
                            update(flags_[TauDecay], pdgid == 15);
                            update(flags_[KsDecay], pdgid == 310);
                            update(flags_[LambdaDecay], pdgid == 3122);
                            update(flags_[JpsiDecay], pdgid == 443);
                            update(flags_[XiDecay], pdgid == 3312);
                            update(flags_[OmegaDecay], pdgid == 3334);
                            update(flags_[SigmaPlusDecay], pdgid == 3222);
                            update(flags_[SigmaMinusDecay], pdgid == 3112);
                        }
                    }
                    update(
                        flags_[Interaction],
                        !flags_[BWeakDecay] &&
                        !flags_[CWeakDecay] &&
                        !flags_[LongLivedDecay] &&
                        !flags_[TauDecay] &&
                        !flags_[KsDecay] &&
                        !flags_[LambdaDecay] &&
                        !flags_[JpsiDecay] &&
                        !flags_[XiDecay] &&
                        !flags_[OmegaDecay] &&
                        !flags_[SigmaPlusDecay] &&
                        !flags_[SigmaMinusDecay]
                    );
                }
                else
                {
                    update(
                        flags_[Interaction],
                        process != G4::Undefined &&
                        process != G4::Unknown &&
                        process != G4::Primary &&
                        process != G4::Hadronic &&
                        process != G4::Conversions
                    );
                }
            }
        }
    }
}


void VertexClassifier::vertexInformation()
{
    // Helper class for clusterization
    typedef std::multimap<double, HepMC::ThreeVector> Clusters;
    typedef std::pair<double, HepMC::ThreeVector> ClusterPair;

    Clusters clusters;

    // Get the main primary vertex from the list
    GeneratedPrimaryVertex const & genpv = genpvs_.back();

    // Get the generated history of the tracks
    VertexHistory::GenVertexTrail & genVertexTrail = const_cast<VertexHistory::GenVertexTrail &> (tracer_.genVertexTrail());

    // Unit transformation from mm to cm
    double const mm = 0.1;

    // Loop over the generated vertexes
    for (
        VertexHistory::GenVertexTrail::const_iterator ivertex = genVertexTrail.begin();
        ivertex != genVertexTrail.end();
        ++ivertex
    )
    {
        // Check vertex exist
        if (*ivertex)
        {
            // Measure the distance2 respecto the primary vertex
            HepMC::ThreeVector p = (*ivertex)->point3d();
            double distance = sqrt( pow(p.x() * mm - genpv.x, 2) + pow(p.y() * mm - genpv.y, 2) + pow(p.z() * mm - genpv.z, 2) );

            // If there is not any clusters add the first vertex.
            if ( clusters.empty() )
            {
                clusters.insert( ClusterPair(distance, HepMC::ThreeVector(p.x() * mm, p.y() * mm, p.z() * mm)) );
                continue;
            }

            // Check if there is already a cluster in the given distance from primary vertex
            Clusters::const_iterator icluster = clusters.lower_bound(distance - vertexClusteringDistance_);

            if ( icluster == clusters.upper_bound(distance + vertexClusteringDistance_) )
            {
                clusters.insert ( ClusterPair(distance, HepMC::ThreeVector(p.x() * mm, p.y() * mm, p.z() * mm)) );
                break;
            }

            // Looping over the vertex clusters of a given distance from primary vertex
            for (;
                    icluster != clusters.upper_bound(distance + vertexClusteringDistance_);
                    ++icluster
                )
            {
                double difference = sqrt (
                                        pow(p.x() * mm - icluster->second.x(), 2) +
                                        pow(p.y() * mm - icluster->second.y(), 2) +
                                        pow(p.z() * mm - icluster->second.z(), 2)
                                    );


                if ( difference > vertexClusteringDistance_ )
                {
                    clusters.insert ( ClusterPair(distance, HepMC::ThreeVector(p.x() * mm, p.y() * mm, p.z() * mm)) );
                    break;
                }
            }
        }
    }

    VertexHistory::SimVertexTrail & simVertexTrail = const_cast<VertexHistory::SimVertexTrail &> (tracer_.simVertexTrail());

    // Loop over the generated particles
    for (
        VertexHistory::SimVertexTrail::reverse_iterator ivertex = simVertexTrail.rbegin();
        ivertex != simVertexTrail.rend();
        ++ivertex
    )
    {
        // Look for those with production vertex
        TrackingVertex::LorentzVector p = (*ivertex)->position();

        double distance = sqrt( pow(p.x() - genpv.x, 2) + pow(p.y() - genpv.y, 2) + pow(p.z() - genpv.z, 2) );

        // If there is not any clusters add the first vertex.
        if ( clusters.empty() )
        {
            clusters.insert( ClusterPair(distance, HepMC::ThreeVector(p.x(), p.y(), p.z())) );
            continue;
        }

        // Check if there is already a cluster in the given distance from primary vertex
        Clusters::const_iterator icluster = clusters.lower_bound(distance - vertexClusteringDistance_);

        if ( icluster == clusters.upper_bound(distance + vertexClusteringDistance_) )
        {
            clusters.insert ( ClusterPair(distance, HepMC::ThreeVector(p.x(), p.y(), p.z())) );
            break;
        }

        // Looping over the vertex clusters of a given distance from primary vertex
        for (;
                icluster != clusters.upper_bound(distance + vertexClusteringDistance_);
                ++icluster
            )
        {
            double difference = sqrt (
                                    pow(p.x() - icluster->second.x(), 2) +
                                    pow(p.y() - icluster->second.y(), 2) +
                                    pow(p.z() - icluster->second.z(), 2)
                                );

            if ( difference > vertexClusteringDistance_ )
            {
                clusters.insert ( ClusterPair(distance, HepMC::ThreeVector(p.x(), p.y(), p.z())) );
                break;
            }
        }
    }

    if ( clusters.size() == 1 )
        flags_[PrimaryVertex] = true;
    else if ( clusters.size() == 2 )
        flags_[SecondaryVertex] = true;
    else
        flags_[TertiaryVertex] = true;
}


bool VertexClassifier::isFinalstateParticle(const HepMC::GenParticle * p)
{
    return !p->end_vertex() && p->status() == 1;
}


bool VertexClassifier::isCharged(const HepMC::GenParticle * p)
{
    const ParticleData * part = particleDataTable_->particle( p->pdg_id() );
    if (part)
        return part->charge()!=0;
    else
    {
        // the new/improved particle table doesn't know anti-particles
        return  particleDataTable_->particle( -p->pdg_id() ) != 0;
    }
}


void VertexClassifier::genPrimaryVertices()
{
    genpvs_.clear();

    const HepMC::GenEvent * event = mcInformation_->GetEvent();

    if (event)
    {
        int idx = 0;

        // Loop over the different GenVertex
        for ( HepMC::GenEvent::vertex_const_iterator ivertex = event->vertices_begin(); ivertex != event->vertices_end(); ++ivertex )
        {
            bool hasParentVertex = false;

            // Loop over the parents looking to see if they are coming from a production vertex
            for (
                HepMC::GenVertex::particle_iterator iparent = (*ivertex)->particles_begin(HepMC::parents);
                iparent != (*ivertex)->particles_end(HepMC::parents);
                ++iparent
            )
                if ( (*iparent)->production_vertex() )
                {
                    hasParentVertex = true;
                    break;
                }

            // Reject those vertices with parent vertices
            if (hasParentVertex) continue;

            // Get the position of the vertex
            HepMC::FourVector pos = (*ivertex)->position();

            double const mm = 0.1;

            GeneratedPrimaryVertex pv(pos.x()*mm, pos.y()*mm, pos.z()*mm);

            std::vector<GeneratedPrimaryVertex>::iterator ientry = genpvs_.begin();

            // Search for a VERY close vertex in the list
            for (; ientry != genpvs_.end(); ++ientry)
            {
                double distance = sqrt( pow(pv.x - ientry->x, 2) + pow(pv.y - ientry->y, 2) + pow(pv.z - ientry->z, 2) );
                if ( distance < vertexClusteringDistance_ )
                    break;
            }

            // Check if there is not a VERY close vertex and added to the list
            if (ientry == genpvs_.end())
                ientry = genpvs_.insert(ientry,pv);

            // Add the vertex barcodes to the new or existent vertices
            ientry->genVertex.push_back((*ivertex)->barcode());

            // Collect final state descendants
            for (
                HepMC::GenVertex::particle_iterator idecendants  = (*ivertex)->particles_begin(HepMC::descendants);
                idecendants != (*ivertex)->particles_end(HepMC::descendants);
                ++idecendants
            )
            {
                if (isFinalstateParticle(*idecendants))
                    if ( find(ientry->finalstateParticles.begin(), ientry->finalstateParticles.end(), (*idecendants)->barcode()) == ientry->finalstateParticles.end() )
                    {
                        ientry->finalstateParticles.push_back((*idecendants)->barcode());
                        HepMC::FourVector m = (*idecendants)->momentum();

                        ientry->ptot.setPx(ientry->ptot.px() + m.px());
                        ientry->ptot.setPy(ientry->ptot.py() + m.py());
                        ientry->ptot.setPz(ientry->ptot.pz() + m.pz());
                        ientry->ptot.setE(ientry->ptot.e() + m.e());
                        ientry->ptsq += m.perp() * m.perp();

                        if ( m.perp() > 0.8 && std::abs(m.pseudoRapidity()) < 2.5 && isCharged(*idecendants) ) ientry->nGenTrk++;
                    }
            }
            idx++;
        }
    }

    std::sort(genpvs_.begin(), genpvs_.end());
}
