import FWCore.ParameterSet.Config as cms

process = cms.Process("TrackOriginAnalyzerTest")

# Message logger
process.load("FWCore.MessageLogger.MessageLogger_cfi")

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

