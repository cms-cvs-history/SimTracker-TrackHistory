import FWCore.ParameterSet.Config as cms

from SimGeneral.HepPDTESSource.pythiapdt_cfi import *

from TrackingTools.TransientTrack.TransientTrackBuilder_cfi import *

from SimTracker.TrackHistory.TrackHistory_cff import *
from SimTracker.TrackHistory.TrackQuality_cff import *

trackClassifier = cms.PSet(
    trackHistory,
    trackQuality,
    hepMC = cms.untracked.InputTag("source"),
    beamSpot = cms.untracked.InputTag("offlineBeamSpot"),
    badD0Pull = cms.untracked.double(3.0),
    longLivedDecayLength = cms.untracked.double(1e-14),
    vertexClusteringDistance = cms.untracked.double(0.0001),
    numberOfInnerLayers = cms.untracked.uint32(2)
)