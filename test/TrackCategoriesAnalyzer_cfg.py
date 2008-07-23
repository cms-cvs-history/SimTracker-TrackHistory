import FWCore.ParameterSet.Config as cms

process = cms.Process("TrackOriginAnalyzerTest")

process.load("FWCore.MessageLogger.MessageLogger_cfi");

process.load("SimTracker.TrackHistory.TrackClassifier_cff")

process.load("Configuration.StandardSequences.FakeConditions_cff")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('')
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)

process.add_( 
  cms.Service("TFileService",
      fileName = cms.string("test.root")
  )
)

process.trackCategoriesAnalyzer = cms.EDFilter("TrackCategoriesAnalyzer",
    process.trackClassifier,
    minimumNumberOfHits = cms.untracked.int32(8),
    minimumTransverseMomentum = cms.untracked.double(1.),
    minimumNumberOfPixelHits = cms.untracked.int32(2),
    maximumChiSquared = cms.untracked.double(5.),
    trackQualityClass = cms.untracked.string('loose')
)

process.p = cms.Path(process.trackCategoriesAnalyzer)


