import FWCore.ParameterSet.Config as cms

process = cms.Process("TrackOriginAnalyzerTest")

process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.load("Configuration.StandardSequences.FakeConditions_cff")

process.load("Configuration.StandardSequences.GeometryIdeal_cff")

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.load("SimTracker.TrackHistory.TrackHistory_cff")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('')
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)
process.trackHistoryAnalyzer = cms.EDFilter("VertexTHA",
    process.trackHistory,
    sourceCut = cms.untracked.double(0.001),
    rootFile = cms.untracked.string('file:test.root')
)

process.p = cms.Path(process.trackHistoryAnalyzer)

