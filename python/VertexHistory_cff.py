import FWCore.ParameterSet.Config as cms

# Magnetic Field, Geometry, TransientTracks
from Configuration.StandardSequences.MagneticField_cff import *

# Track Associators
from SimTracker.TrackAssociation.TrackAssociatorByChi2_cfi import *
from SimTracker.TrackAssociation.TrackAssociatorByHits_cfi import *
from SimTracker.VertexAssociation.VertexAssociatorByTracks_cfi.py import *


# Track history parameters
vertexHistory = cms.PSet(
    bestMatchByMaxValue = cms.untracked.bool(True),
    trackingTruth = cms.untracked.InputTag("mergedtruth","MergedTrackTruth"),
    trackAssociator = cms.untracked.string('TrackAssociatorByHits'),
    trackProducer = cms.untracked.InputTag("generalTracks"),
    vertexAssciator = cms.untracked.string('VertexAssociatorByTracks'),
    vertexProducer = cms.untracked.InputTag("?")
)


