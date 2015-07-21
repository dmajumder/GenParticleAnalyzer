import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/user/d/devdatta/afswork/CMSREL/CMSSW_7_1_14/src/Configuration/Generator/test/TpbTotH_topdecay_RH_narrow_M1200_TuneCUETP8M1_13TeV_GEN.root') 
    )

process.TFileService = cms.Service("TFileService",
       fileName = cms.string("TpbTotH_topdecay_RH_narrow_M1200_GENHists.root") 
       )

process.gen = cms.EDAnalyzer('GenParticleAnalyzer'
)


process.p = cms.Path(process.gen)
