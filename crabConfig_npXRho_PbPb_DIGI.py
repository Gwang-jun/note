#

from CRABClient.UserUtilities import config
config = config()
config.General.transferOutputs = True
config.General.requestName = 'MC_20190516_BsToJpsiPhi_Bpt5p0_Pthat5_TuneCP5_DIGI_RAW_Drum5F'
config.General.workArea = 'crab_projects'
config.JobType.psetName = 'step1_digi.py'
config.JobType.pluginName = 'PrivateMC'
config.JobType.maxMemoryMB = 2500
config.JobType.pyCfgParams = ['noprint']
config.JobType.numCores = 4
# config.JobType.allowUndistributedCMSSW = True

##
# config.Data.inputDataset = '/Hydjet_Pythia8_X3872ToJpsiRho_nonprompt_Pthat5_TuneCP5_5020GeV_Drum5F/wangj-MC_20190516_Xpt0p0_5020GeV_GEN_SIM-c6f4f87e1b47e594872980bfb96a5efc/USER'
# config.Data.inputDataset = '/Hydjet_Pythia8_X3872ToJpsiRho_nonprompt_Pthat10_TuneCP5_5020GeV_Drum5F/wangj-MC_20190516_Xpt0p0_5020GeV_GEN_SIM-032242104c7c9b007039250e8074bbb7/USER'
# config.Data.inputDataset = '/Hydjet_Pythia8_X3872ToJpsiRho_nonprompt_Pthat15_TuneCP5_5020GeV_Drum5F/wangj-MC_20190516_Xpt0p0_5020GeV_GEN_SIM-c4286cd8b1f1e070c5095886808cb75e/USER'
# config.Data.inputDataset = '/Hydjet_Pythia8_X3872ToJpsiRho_nonprompt_Pthat30_TuneCP5_5020GeV_Drum5F/wangj-MC_20190516_Xpt0p0_5020GeV_GEN_SIM-2d38c01c52e179527f7585006907a509/USER'
config.Data.inputDataset = '/Hydjet_Pythia8_BsToJpsiPhi_Bpt5p0_Pthat5_TuneCP5_5020GeV_Drum5F/gwkim-MC_20190516_Bpt5p0_5020GeV_GEN_SIM-a1870d500cc892054235382816bf3d3d/USER'
config.Data.inputDBS = 'phys03'
config.Data.unitsPerJob = 50000 ##
config.Data.totalUnits = -1
config.Data.splitting = 'FileBased'
config.Data.publication = True
config.Data.outputDatasetTag = 'MC_20190516_Bpt5p0_5020GeV_DIGI_RAW'

##
config.Site.storageSite = 'T2_US_MIT'
config.Site.blacklist = ['T2_US_Nebraska','T2_US_UCSD','T2_US_Wisconsin','T2_FR_GRIF_IRFU','T3_US_UCR','T3_US_Rutgers','T3_BG_UNI_SOFIA','T3_IT_Perugia']
# config.Site.ignoreGlobalBlacklist = True
# config.section_("Debug")
# config.Debug.extraJDL = ['+CMS_ALLOW_OVERFLOW=False']

