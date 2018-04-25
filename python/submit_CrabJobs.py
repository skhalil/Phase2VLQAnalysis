#
#
# example configuration for multicrab submission. Needs to be adapted to analysis needs
# (thanks to Sandhya for providing this)
#
#

from CRABClient.UserUtilities import config, getUsernameFromSiteDB
import os
from CRABAPI.RawCommand import crabCommand
from CRABClient.ClientExceptions import ClientException
from httplib import HTTPException
from multiprocessing import Process

#batch = "signal"
#batch = "data"
batch = "background"

if "signal" in batch:
	mydatasets = 'ana_crab_input_files_signal.txt'
	filesplit = 'FileBased'
	units = 1
if "data" in batch:
        mydatasets = 'ana_crab_input_files_data.txt'
        filesplit = 'LumiBased'
	units = 10	
if "background" in batch:
        mydatasets = 'ana_crab_input_files_background.txt'
        filesplit = 'FileBased'
	units = 1

config = config()

config.General.requestName = None
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.workArea = 'crab_files_test4_25April18'

config.JobType.psetName = '/home/t3-ku/jaking/CMSSW_9_3_2/src/Upgrades/VLQAnalyzer/python/ConfFile_cfg.py'
config.JobType.pluginName = 'Analysis'

config.JobType.maxJobRuntimeMin = 2000
config.JobType.maxMemoryMB = 2500
config.Data.inputDataset = None
config.Data.inputDBS = 'global'
#config.Data.inputDBS = 'phys03'
config.Data.splitting = filesplit
config.Data.unitsPerJob = units
config.Data.outLFNDirBase = '/store/user/jaking/datasets/semi_lep/'
config.Data.allowNonValidInputDataset = True
config.Data.publication = False
config.Data.secondaryInputDataset = None

config.section_("Site")
#config.Site.whitelist = [ 'T2_US_*' ]
config.Site.storageSite = 'T2_US_Nebraska'

def submit(config):
    try:
        crabCommand('submit', config = config)
    except HTTPException as hte:
        print "Failed submitting task: %s" % (hte.headers)
    except ClientException as cle:
        print "Failed submitting task: %s" % (cle)

datasetsFile = open( mydatasets )
jobsLines = datasetsFile.readlines()
for ijob in jobsLines :
    s = ijob.rstrip()
    if (len(s)==0 or s[0][0]=='#'): continue
    cdi = s + 'MINIAODSIM'
    cgr = s.split('/')[1]
    r = s.replace("MiniAOD","DR" )
    cds = r + 'GEN-SIM-RECO'
    config.Data.inputDataset = cdi
    config.General.requestName = cgr
    config.Data.secondaryInputDataset = cds
    print "Submitting to Crab:"
    print "Inputdataset: ",cdi
    print "requestName: ",cgr
    print "SecondaryInputdataset: ",cds
    print 
    submit(config)



#config.Data.inputDataset = '/TprimeBToTH_M-2000_Width-10p_LH_TuneCUETP8M2T4_14TeV-madgraph-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2-v1/MINIAODSIM'
#config.General.requestName = 'TprimeBToTH_M-2000_Width-10p_LH_TuneCUETP8M2T4_14TeV-madgraph-pythia8'
#config.Data.secondaryInputDataset = '/TprimeBToTH_M-2000_Width-10p_LH_TuneCUETP8M2T4_14TeV-madgraph-pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2-v1/GEN-SIM-RECO'
#submit(config)



