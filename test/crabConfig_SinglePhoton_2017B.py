from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'MVA_monopho17BUL'
config.General.workArea = 'crab_projects'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'run_monopho17.py'
config.JobType.allowUndistributedCMSSW = True

config.section_("Data")
config.Data.inputDataset = '/SinglePhoton/Run2017B-09Aug2019_UL2017-v1/MINIAOD'
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 50
config.Data.lumiMask = 'Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt'
#config.Data.runRange = '275776-275782'
config.Data.outLFNDirBase = '/store/user/xdelin/monopho17/MVA'


config.section_("Site")
config.Site.storageSite = 'T3_US_FNALLPC'
