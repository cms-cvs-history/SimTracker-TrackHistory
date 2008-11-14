# Playback file
import FWCore.ParameterSet.Config as cms

# import of standard configurations
from Configuration.StandardSequences.Services_cff import *
from FWCore.MessageService.MessageLogger_cfi import *
from Configuration.StandardSequences.MixingNoPileUp_cff import *
from Configuration.StandardSequences.GeometryPilot2_cff import *
from Configuration.StandardSequences.MagneticField_38T_cff import *
from Configuration.StandardSequences.RawToDigi_cff import *
from Configuration.StandardSequences.Reconstruction_cff import *
from Configuration.StandardSequences.FrontierConditions_GlobalTag_cff import *
from Configuration.EventContent.EventContent_cff import *
from Configuration.StandardSequences.Digi_cff import *

# playback
del RandomNumberGeneratorService.theSource
RandomNumberGeneratorService.restoreStateLabel = cms.untracked.string('randomEngineStateProducer')
from SimGeneral.MixingModule.mixNoPU_cfi import *
mix.playback = cms.untracked.bool(True)

# trackingTruth
from SimGeneral.TrackingAnalysis.trackingParticles_cfi import *

# tracking truth and SiStrip(Pixel)DigiSimLinks
trackingTruth = cms.Sequence(mix * doAllDigi * trackingParticles)

# reconstruction
playback = cms.Sequence(RawToDigi + trackingTruth + reconstruction)