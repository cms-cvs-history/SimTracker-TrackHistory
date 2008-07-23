import FWCore.ParameterSet.Config as cms

process = cms.Process("TrackOriginAnalyzerTest")

# Message logger
process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

# TrackHistory setup
process.load("SimTracker.TrackHistory.TrackHistory_cff")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('')
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

process.trackOriginAnalyzer = cms.EDFilter("TruthTOA",
    process.trackHistory,

    # output file 
    rootFile = cms.untracked.string('file:test.root'),

    # false = count particle and antiparticle as the same type
    # true = caount particle and antiparticle as different type
    antiparticles = cms.untracked.bool(False),

    status2 = cms.untracked.bool(False),

    veto = cms.untracked.PSet(
        list1 = cms.untracked.vstring(),
        list2 = cms.untracked.vstring('pi+')
    ),

    # run the analysis per each jet(no jet overlap) 
    insideJet = cms.untracked.bool(True),

    # Minimum number of matched hits
    matchedHits = cms.untracked.int32(8)
)

process.p = cms.Path(process.trackOriginAnalyzer)

