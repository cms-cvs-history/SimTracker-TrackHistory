# Playback test file

import FWCore.ParameterSet.Config as cms

process = cms.Process('TrackCategoriesWithPlayback')

process.load("SimTracker.TrackHistory.Playback_cff")
process.load("SimTracker.TrackHistory.TrackClassifier_cff")

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

# Other statements
process.GlobalTag.globaltag = 'IDEAL_V9::All'

# Path
process.path = cms.Path(process.playback * process.trackCategoriesAnalyzer)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
process.source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( ( 
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/00208556-CA8C-DD11-A1BC-000E0C4D3F34.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/0099EC84-CF8C-DD11-9C30-00E08134B780.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/009BE5C8-FE8C-DD11-9350-001C23C0F175.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/009C43B6-898D-DD11-9959-001A92544590.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/02136B1F-CE8F-DD11-BB66-001E6849D24C.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/02140F50-F68C-DD11-8C67-00304834BB58.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/022A3A5B-B48C-DD11-935A-00163691DEF2.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/0254FA7B-BD8C-DD11-9873-001A9227D34F.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/0261C762-C58C-DD11-ADFA-0016368E0C2C.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/0273F700-C78C-DD11-A079-001EC94BF93F.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/02987061-F08C-DD11-9D0C-001A9227D387.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/02CBBC29-B38C-DD11-BD08-001E8C7A922C.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/02CC0F9D-048D-DD11-B9EE-0002B3E92671.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/02FDCB6B-BC8C-DD11-9CF6-00163691DEF6.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/043FC27B-BD8C-DD11-9583-001731DEBFF8.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/04751647-D18C-DD11-91D0-0016368E0AE8.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/049DA585-4E8D-DD11-964D-003048335548.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/06193997-2A8E-DD11-A608-00E0813006E6.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/065F2760-8D8D-DD11-88F0-001E8CC04088.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/067F964B-C88C-DD11-8704-0016368E0820.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/06E5FCEC-958D-DD11-9299-001A9227D32F.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/08229C3C-BE8C-DD11-B49C-00163691D196.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/0877A8BB-7D8D-DD11-8BE4-001A92971CCE.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/08929760-E88C-DD11-A690-00163691DA92.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/0892EA4B-EB8C-DD11-A4A0-001A9243D516.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/08BF6F1D-BD8C-DD11-9B7A-0002B3E92671.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/0A5440A3-838D-DD11-B5E0-001A92971C6A.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/0A837F75-438D-DD11-A851-000E0C4D3F34.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/0A894285-CF8C-DD11-869B-00E081339578.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/0AB5332F-C28C-DD11-81F9-0016368E0DE0.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/0AE1EC02-138D-DD11-96F4-0002B3E92671.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/0C009ABC-ED8E-DD11-9ACF-00E0813006E6.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/0C01A50A-968C-DD11-AFEA-0002B3E92671.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/0C2116F0-BA8C-DD11-9332-0016368E0820.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/0C34B6F6-B58C-DD11-A402-001E8CCCE140.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/0C5DF1CB-9A8C-DD11-8E52-0002B3E92671.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/0C6A7C72-D38D-DD11-AF9E-001AA033166A.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/0C791DB0-728D-DD11-9C4F-00304865C49C.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/0C7C4BA9-B28C-DD11-89EC-00163691DCC6.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/0C8D7768-C58C-DD11-BBA5-00163691DAE6.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/0C9F71AA-C78C-DD11-8791-00E081339574.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/0CC4D7A2-7F8D-DD11-98D3-001A9227D20D.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/0CC67F47-BE8C-DD11-A888-00163691D196.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/0CE03147-AC8C-DD11-828B-0016368E0A9C.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/0CF0BDF8-EA8C-DD11-902F-000E0C4D3F34.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/0E0BEEE2-718D-DD11-8DBE-00304865C28A.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/0E39018D-CD8C-DD11-8BD4-00E08133E4B2.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/0E6D2DAC-B28C-DD11-AA4C-00163691DCC6.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/0EAAB9A3-C98C-DD11-9640-00E08133D50E.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/0EBBC760-E88C-DD11-9007-00163691DA92.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/0EDFE904-A48C-DD11-9D04-0016368E0A9C.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/100F304E-AB8C-DD11-B475-001E681EA272.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/101FA157-B78C-DD11-A6EE-003048334511.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/10354024-C28C-DD11-9BB4-00163691DC86.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/1043837F-388D-DD11-B11E-0002B3E92671.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/104595FB-C68C-DD11-A9E9-001EC94BA187.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/104827AD-C68C-DD11-8B0A-001EC94BFB3E.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/10A370D8-D38C-DD11-B8E6-0016368E0DE0.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/1229384F-D18C-DD11-9FD9-001A9227D3CB.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/123DC0B1-778D-DD11-8F25-001A9227D34F.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/127D11B2-F08C-DD11-A00F-0030487D07BA.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/12D8EF5D-E88C-DD11-8CAF-00163691DC0A.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/12EC462C-CD8C-DD11-A814-003048554726.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/1407A0B7-B28C-DD11-A561-0016368E0A9C.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/143792B4-C58C-DD11-8323-00E081339574.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/14593B58-CF8C-DD11-9500-0002B3E92671.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/145A4DF1-B98C-DD11-A0BC-001A9243D668.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/147BBC70-098D-DD11-8AC6-000E0C4D3F34.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/1492AF64-E88C-DD11-93EE-0016368E0A9C.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/14D7451E-C28C-DD11-94F9-00163691DEF2.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/14F1C291-D08C-DD11-B4E8-00163691DA92.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/14FD53CB-808D-DD11-B78E-001A92971CCE.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/1641A9AB-508D-DD11-868A-003048335548.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/1671F3A4-C98C-DD11-8C89-00E08133E4B2.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/16A023B3-C58C-DD11-A81D-00E08133E4B2.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/16AC0341-B08C-DD11-BE70-000E0C4D3F34.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/16B83310-A58C-DD11-AB47-0002B3E92671.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/16F3F94E-968D-DD11-86C0-001A9227D383.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/182C2652-1D8D-DD11-AF66-000E0C4D3F34.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/182C295B-938D-DD11-975D-001A9227D3C1.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/182FBE61-CD8C-DD11-B7C3-000E0C4D3F34.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/184C8567-D48C-DD11-B316-001EC94BF6F7.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/185507A2-CD8C-DD11-8B3F-00E08133CD36.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/18ABED42-BE8C-DD11-B7BF-00163691D992.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/18AEC3E6-D78C-DD11-9AE2-0002B3E92671.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/18CECE7B-BD8C-DD11-BA4F-001A9227D371.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/1A0E706E-EC8C-DD11-B310-001C23C102BB.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/1A343A59-5C8D-DD11-982D-0002B3E92671.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/1A49E6BE-ED8E-DD11-93A5-00E0813006E6.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/1A57FAF5-B98C-DD11-8149-0016368E0DE0.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/1A6CF41D-A68C-DD11-9859-0016368E0A9C.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/1AAA9620-968C-DD11-8833-000E0C4D3F34.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/1AC820EF-198D-DD11-9F3A-0002B3E92671.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/1ACE5ADF-AD8C-DD11-93C7-0030482EE630.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/1AEB3863-E88C-DD11-B900-00163691D762.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/1C007185-768D-DD11-A500-001A9227D3C1.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/1C053291-E18C-DD11-BB92-001A9227D3A9.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/1C298BF5-B98C-DD11-969C-0016368E0DE0.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/1C3B16EB-218D-DD11-855C-003048322BF6.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/1C9C66A2-CD8C-DD11-A2B6-00E08133CDA0.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/1CACF859-AB8C-DD11-A358-003048344A8C.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/1CB44869-078D-DD11-BAA3-000E0C4D3F34.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/1CC6BE61-CD8C-DD11-A8B2-000E0C4D3F34.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/1CCEB850-B48C-DD11-BD7C-00163691DE22.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/1E5D81BF-D28C-DD11-A16D-001C23C105F2.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/1E5F20DB-D48C-DD11-BC9F-001EC94BF999.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/1E70A26F-BC8C-DD11-A66F-00163691DC3E.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/1E71232A-2F8D-DD11-929A-0002B3E92671.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/1E7FC95A-908D-DD11-B691-001A9227D3BD.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/1EB492DF-538D-DD11-A4C2-0002B3E92671.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/1EDD120F-FE8C-DD11-89A9-00E0814296A8.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/2024D2E6-C78D-DD11-9824-001A9243D528.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/20D64F57-6E8D-DD11-B73F-0002B3E92671.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/20E12066-D48C-DD11-AA29-001EC94BA187.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/2209A4EB-C48C-DD11-9E3E-001A9227D3A9.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/2263B103-908D-DD11-AFFE-00304858A541.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/226E453C-B48C-DD11-BBCC-003048322A68.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/229DCAD9-DA8C-DD11-BD97-00304858A6B3.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/240C8456-CA8C-DD11-8A63-000E0C4D3F34.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/242DF4E0-178D-DD11-9721-00304833457E.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/244BBA99-CB8C-DD11-8BD3-00E08133E4B2.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/2469503E-BE8C-DD11-9B68-00163691DA0A.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/2471E96A-458D-DD11-B86E-001A92BCDC7A.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/249F168B-CA8C-DD11-98DB-0002B3E92671.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/24B954E8-D38C-DD11-94C7-0016368E0C2C.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/24B9FC4A-B58C-DD11-B2A6-00304834BB00.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/24FF292E-B48C-DD11-97E3-0030483229E0.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/2676D40A-868D-DD11-8FBF-001A9243D528.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/267DDE56-108D-DD11-9AAA-001A9227D3D1.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/267FA678-D38D-DD11-AF9B-00E08130071E.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/268ADE4B-C88C-DD11-B241-0016368E0A9C.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/26C5167A-DB8C-DD11-ABD9-000E0C4D3F34.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/28044D49-E58C-DD11-AF77-000E0C4D3F34.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/285EBDA4-C88D-DD11-8CD6-001A92971C6A.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/2872C006-4B8D-DD11-B14D-0002B3E92671.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/288CE55D-EC8C-DD11-AFC5-0030485546C2.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/2890ABB2-F08C-DD11-BC1D-0030487D676A.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/28A204A5-C58C-DD11-B499-00163691DC0E.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/28B90B0A-3D8D-DD11-BA81-0002B3E92671.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/28B96F4B-7E8D-DD11-89C9-001A9227D3DD.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/28F93C40-BE8C-DD11-8D4A-00163691DEEE.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/28FEB145-EE8C-DD11-8597-003048335594.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/2A33EFA3-C08D-DD11-AF1A-001E8C7A922C.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/2A370F4B-AB8C-DD11-9783-001E682F2052.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/2A3952DD-C98C-DD11-98D2-00163691DCC6.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/2A45ED7F-618D-DD11-B839-0002B3E92671.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/2A7FD170-BC8C-DD11-A338-001B243DEF3F.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/2A845B46-578D-DD11-ABCB-000E0C4D3F34.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/2A8F2E80-3F8D-DD11-AB87-001EC94B4F63.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/2AC61645-AC8C-DD11-A8CD-00163691DDA2.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/2AC91BAA-C18C-DD11-B10D-001A9227D3BD.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/2C93B6D2-B98C-DD11-91AA-00163691DD9A.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/2CB0416B-008D-DD11-B670-000E0C4D3F34.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/2CB30D0A-3D8D-DD11-A31F-0002B3E92671.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/2E535FB7-F68C-DD11-A3E1-0002B3E92671.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/2E8FD016-B88C-DD11-9F1C-00304834BB00.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/2EA1A26E-AE8C-DD11-BD7D-0016368E0AE8.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/2EB5E1A7-C78C-DD11-A6CA-00E08133D50E.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/2EC248AC-088D-DD11-A9CC-0002B3E92671.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/2EEE0B30-908D-DD11-A814-00304858A663.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/2EF750DF-388D-DD11-BA5D-001E8C7A91E4.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/305025ED-C98C-DD11-867D-0002B3E92671.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/30A4D5DB-948D-DD11-8B4D-001A9227D359.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/30C56E26-B48C-DD11-8495-00163691DA92.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/30F024D5-F68C-DD11-A8DF-000E0C4D3F34.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/32148721-968C-DD11-B6B8-000E0C4D3F34.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/3221828C-CD8C-DD11-8B73-00E08133D50E.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/323781D7-948D-DD11-A0A6-00E081403270.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/323F027E-428D-DD11-B5A4-001E8CCCE148.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/32AE97B1-C98C-DD11-924A-00163691D762.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/32D1D743-AC8C-DD11-911E-00163691D1AE.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/32D67D78-AE8C-DD11-A4C3-00163691DC0E.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/32D8C42C-518D-DD11-86FC-003048335548.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/341564D2-B98C-DD11-854C-00163691DD9A.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/341C57BD-298D-DD11-A73E-0002B3E92671.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/3420B34F-5A8D-DD11-A9DE-001A9227D363.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/347E89F7-B58C-DD11-AD4B-001A9227D3ED.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/349957BD-748D-DD11-8F66-001A9227D3BD.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/34A289A9-2F8D-DD11-B45F-000E0C4D3F34.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/34D38406-B68C-DD11-ADFE-001E8CCCE148.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/34D644B2-C98C-DD11-B28F-0016368E0A9C.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/34F6ACCC-B88C-DD11-AB8E-0002B3E92671.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/36327F7B-AE8C-DD11-84AB-0016368E0A9C.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/364275C4-188D-DD11-A272-003048322BF6.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/3677A26C-228D-DD11-8EB6-000E0C4D3F34.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/36793F12-B68C-DD11-ABF3-001E8CC0413A.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/368CC256-128D-DD11-BDC6-003048322C3E.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/36C66E5F-EC8C-DD11-BD0E-001C23C105CF.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/36D37A22-158D-DD11-B4FB-003048322C3E.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/380711A0-C18C-DD11-A7FD-000E0C4D3F34.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/380A8231-CB8C-DD11-8CE9-00304833529A.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/381916C1-8D8D-DD11-BAEA-001A9227D1EB.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/383D4C28-4A8D-DD11-9058-002215031696.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/384B5246-BE8C-DD11-AFA9-00163691D18E.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/38587158-438D-DD11-AC9A-001A9243D516.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/3888731A-968C-DD11-8496-000E0C4D3F34.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/388BC8A2-D38C-DD11-B91C-00163691DEEA.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/38AEB8D5-D38C-DD11-8EE6-00163691DC86.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/38B5F426-D18C-DD11-8831-0016368E0AE8.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/38DF5819-008D-DD11-A10F-0002B3E92671.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/38F01141-B48C-DD11-A15D-0030483355D4.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/38FD92FD-B98C-DD11-8D28-00163691DEF2.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/3A08CE47-4F8D-DD11-BEAC-0030483344E2.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/3A33C2F0-B98C-DD11-A570-00163691DC86.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/3A33CF1B-C88C-DD11-9887-00163691DEEA.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/3ABDBD6A-468D-DD11-9056-0002B3E92671.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/3AE4010B-968C-DD11-A0E1-0002B3E92671.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/3AFB2299-CF8C-DD11-84B0-00E08134420C.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/3C0ED4BF-BF8C-DD11-A8C5-00E08133E4B2.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/3C5F8DDD-008D-DD11-80AD-0019B9CABFAC.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/3C9980F6-B58C-DD11-8076-001A9243D528.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/3CA679C8-D38C-DD11-88F5-00163691DEEA.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/3CCD9E29-B38C-DD11-B69C-001E8CCCE140.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/3CE44D27-CC8C-DD11-B8B1-001A9227D32F.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/3CFC0F95-CF8C-DD11-B799-00E08134B780.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/3E5FA34A-E48D-DD11-9E35-00E08130071E.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/3E7DEA44-B38C-DD11-B3A0-00163691DC62.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/3ECA463C-C28C-DD11-AFB9-00163691DCC6.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/3EF08E1D-C88C-DD11-A5E4-0016368E0C2C.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/401A40DA-C98C-DD11-92CF-001B243DE10F.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/4081E37C-3F8D-DD11-9719-001EC94BA146.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/408426C3-948D-DD11-8B98-00E081403270.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/408575BA-438D-DD11-BD95-0030483345FC.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/40A8AB8B-AB8C-DD11-9828-000E0C4DE422.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/4206B0D7-728D-DD11-91C0-0002B3E92671.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/423931C4-F08C-DD11-946E-003048553D02.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/424EB333-FE8C-DD11-BBE5-000E0C4D3F34.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/425DDD4C-428D-DD11-8BD5-001A92971CDC.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/4272F0B4-C68C-DD11-B803-00E081237B25.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/42791F85-CF8C-DD11-B5C7-00E08133CD36.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/42943AC4-958D-DD11-8FB6-001A9227D1EB.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/42943F1F-B48C-DD11-88CC-00163691DE22.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/42964470-D48C-DD11-A094-001C23C0B673.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/42E08831-488D-DD11-8297-00304833442E.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/42E399EE-178D-DD11-83D9-003048322C3E.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/42F27E78-AE8C-DD11-A25B-00163691DC0E.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/440E7147-018D-DD11-B9F2-001EC94BF6F7.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/440F5AF2-7F8D-DD11-AC30-001A9227D1EB.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/4415801F-968C-DD11-82CF-000E0C4D3F34.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/4422B5D2-B48C-DD11-A69B-00163691DE22.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/44355BD8-B48C-DD11-9451-0002B3E92671.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/44376B4C-388D-DD11-AB0E-0002B3E92671.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/44A50BF9-908D-DD11-8F92-001731DEC014.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/462FD356-128D-DD11-90B0-003048322CAA.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/4669ECEA-758D-DD11-9A05-0018F387307A.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/46790276-748D-DD11-834A-001A9227D3BD.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/46DE6AF7-C68C-DD11-A411-001EC94BA3B8.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/46E13D3E-BE8C-DD11-8047-00163691D996.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/48060967-AD8C-DD11-A40F-001E682F1FF4.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/4826B947-AC8C-DD11-8E1D-00163691DC0E.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/48480584-CF8C-DD11-B8E8-00E08133CDA0.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/485E467C-BC8C-DD11-B136-00163691D196.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/485F9322-C88C-DD11-AA85-00163691DEEA.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/4877E447-BC8C-DD11-A609-00163691DA8E.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/48AE2A9A-D38C-DD11-9223-001EC94BF09F.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/48AF1DA2-CD8C-DD11-865B-00E08133CDA0.root',
       '/store/mc/Summer08/InclusiveMu5Pt50/GEN-SIM-RAW/IDEAL_V9_v1/0000/48B4B22C-C28C-DD11-86F9-00163691DC22.root') );

