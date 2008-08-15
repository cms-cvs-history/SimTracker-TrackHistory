import FWCore.ParameterSet.Config as cms

process = cms.Process("TrackOriginAnalyzerTest")

# Message logger
process.load("FWCore.MessageLogger.MessageLogger_cfi")

# Conditions
process.load("Configuration.StandardSequences.FakeConditions_cff")

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
    dataset = "/RelValTTbar/CMSSW_2_1_0_STARTUP_V4_v1/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO",
    site = "cmssrm.fnal.gov"
)

