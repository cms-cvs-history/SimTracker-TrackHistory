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
process.trackHistoryAnalyzer = cms.EDFilter("VertexTHA",
    process.trackHistory,
    sourceCut = cms.untracked.double(0.001),
    rootFile = cms.untracked.string('file:test.root')
)

process.p = cms.Path(process.trackHistoryAnalyzer)
process.PoolSource.fileNames = ['/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0003/1C57E5B2-C040-DD11-AE7A-000423D98804.root', '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0003/84384243-C140-DD11-AC97-001617DBCF1E.root', '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0003/966EF25F-BB40-DD11-A8A7-000423D986A8.root', '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0003/C2EA6F89-C140-DD11-9A3E-001617DBD5AC.root', '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0003/DA99824A-BD40-DD11-9ACC-001617E30E2C.root', 
    '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0004/02367270-E840-DD11-A5F6-000423D6B358.root', '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0004/0A8354F2-0141-DD11-8823-001617E30E28.root', '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0004/10AE569A-DB40-DD11-ADF4-000423D98DB4.root', '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0004/1615B267-2C41-DD11-80CA-000423D6CA6E.root', '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0004/180E6089-0241-DD11-B2F2-001617DBCF1E.root', 
    '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0004/18854915-DA40-DD11-9027-000423D6A6F4.root', '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0004/2A9CFBA4-0F41-DD11-95C9-000423D6CA72.root', '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0004/3AF44B7C-0241-DD11-8439-001617DC1F70.root', '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0004/425C94E0-F440-DD11-8315-001617DF785A.root', '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0004/4E54E2C7-0141-DD11-AF11-001617C3B706.root', 
    '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0004/4ED5F79D-0D41-DD11-A707-001617E30CA4.root', '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0004/5297C0CF-0141-DD11-99F2-001617C3B6DE.root', '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0004/56959392-CE40-DD11-BF8B-0019DB2F3F9B.root', '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0004/581AE09F-1A41-DD11-B1EC-000423D6CA6E.root', '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0004/5E577CB3-1841-DD11-90B7-001617DF785A.root', 
    '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0004/608CEC77-EE40-DD11-A28E-001617DC1F70.root', '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0004/609AFAC5-0141-DD11-B850-001617C3B76E.root', '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0004/64B30AFD-FC40-DD11-BD9B-000423D6BA18.root', '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0004/6823588D-1941-DD11-A9B7-000423D6B358.root', '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0004/682F1CA0-FF40-DD11-B921-000423D6C8EE.root', 
    '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0004/68CDBFCF-0141-DD11-8A02-001617C3B70E.root', '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0004/789261F9-0141-DD11-911A-001617E30CE8.root', '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0004/80E4E1D8-C540-DD11-8448-001617E30F58.root', '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0004/8A9065E5-ED40-DD11-B89A-000423D6CA6E.root', '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0004/8EB0B875-FF40-DD11-A5F9-000423DD2F34.root', 
    '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0004/90CCB731-C740-DD11-B56A-000423D98804.root', '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0004/B4DFA622-E240-DD11-A6CD-000423D6CA72.root', '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0004/CCF2C7F7-0141-DD11-B834-001617DBD472.root', '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0004/D2E409DF-F440-DD11-8B3F-001617C3B64C.root', '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0004/E8EE6A3B-D140-DD11-8710-000423D6CA02.root', 
    '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0004/EA3F37EB-0141-DD11-93A7-000423D991D4.root', '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0004/EA7D0C5D-C340-DD11-B799-001617C3B79A.root', '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0004/ECD138E7-ED40-DD11-AED0-000423D6B42C.root', '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0004/FEF75486-C840-DD11-AC7D-000423D6CA02.root', '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0005/124C388D-1641-DD11-A791-001617DBD224.root', 
    '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0005/126C1DC7-1341-DD11-B350-000423D9890C.root', '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0005/14EB564B-1641-DD11-91C6-001617C3B706.root', '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0005/181AFF1C-1941-DD11-92B7-001617DBCF90.root', '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0005/1E187CD7-0541-DD11-9A09-001617E30E2C.root', '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0005/1EEF0BD9-1A41-DD11-A67C-000423D6006E.root', 
    '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0005/2A6FDDEA-0541-DD11-952B-000423D94990.root', '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0005/34601CF4-0841-DD11-AC15-001D09F2983F.root', '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0005/3AC2B264-0341-DD11-806A-001617E30F4C.root', '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0005/52C2F593-0541-DD11-92E9-001D09F23A84.root', '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0005/58995EF3-1E41-DD11-921E-000423D6CA6E.root', 
    '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0005/5AE251C6-0841-DD11-B618-001D09F2503C.root', '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0005/5E6F381A-1A41-DD11-887E-000423D99896.root', '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0005/5E71DD4C-0E41-DD11-A162-000423D98804.root', '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0005/64CC1485-1F41-DD11-A514-000423D9939C.root', '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0005/663AC722-0541-DD11-B683-001617E30D06.root', 
    '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0005/6E81C3EB-0441-DD11-B32B-001D09F2447F.root', '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0005/743D9D37-1D41-DD11-A6EB-000423D9870C.root', '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0005/7AECEC36-0A41-DD11-A416-001D09F241F0.root', '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0005/88E5FBD6-0841-DD11-9092-0030487C5CFA.root', '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0005/8A8F338D-0541-DD11-B4F1-000423D99996.root', 
    '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0005/9080DB9B-1B41-DD11-8428-000423D98EC8.root', '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0005/B08E97C7-1941-DD11-89C1-001617C3B6DC.root', '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0005/B452397E-0341-DD11-848F-0030487A322E.root', '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0005/BC0E7058-0B41-DD11-B5A3-001D09F244BB.root', '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0005/BC8F6B08-1B41-DD11-9260-001617DBCF1E.root', 
    '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0005/D0570008-1941-DD11-8883-000423D6C8EE.root', '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0005/D4D1500B-0B41-DD11-839D-0019B9F709A4.root', '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0005/DAC316EC-0741-DD11-912D-001617DBD556.root', '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0005/E4733E33-0741-DD11-A2D7-001D09F23F2A.root', '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0005/F45994AC-0541-DD11-972F-000423D98B08.root', 
    '/store/relval/2008/6/22/RelVal-RelValTTbar-1214048167-IDEAL_V2-2nd/0005/FC6700F1-1641-DD11-B606-000423D98BC4.root']


