import FWCore.ParameterSet.Config as cms

process = cms.Process("TrackOriginAnalyzerTest")

# Message logger
process.load("FWCore.MessageLogger.MessageLogger_cfi")

# Conditions
process.load("Configuration.StandardSequences.FakeConditions_cff")

# Geometry
process.load("Configuration.StandardSequences.GeometryIdeal_cff")

# TrackHistory setup
process.load("SimTracker.TrackHistory.TrackClassifier_cff")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('')
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

process.trackHistoryAnalyzer = cms.EDFilter("SimpleTHA",
    process.trackClassifier
)

process.p = cms.Path(process.trackHistoryAnalyzer)

import SimTracker.TrackHistory.DBSPlugin as DBSPlugin

process.PoolSource.fileNames = DBSPlugin.get(
    dataset = "/RelValQCD_Pt_80_120/CMSSW_2_1_2_STARTUP_V5_10TeV_v1/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO",
    site = "srm.cern.ch"
)

