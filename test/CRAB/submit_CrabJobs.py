from CRABClient.UserUtilities import config, getUsernameFromSiteDB
import os
from CRABAPI.RawCommand import crabCommand
from CRABClient.ClientExceptions import ClientException
from httplib import HTTPException
from multiprocessing import Process

PPNdir = os.environ['CMSSW_BASE']+'/src/Upgrades/VLQAnalyzer/'

#batch = "signal"
batch = "background"

if "signal" in batch:
	mydatasets = 'ana_crab_input_files_signal.txt'
	filesplit = 'FileBased'
	units = 1	
if "background" in batch:
        mydatasets = 'ana_crab_input_files_background.txt'
        filesplit = 'FileBased'
	units = 1

config = config()

config.section_('General')
config.General.requestName = None
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.workArea = batch

config.section_('JobType')
config.JobType.psetName = PPNdir+'test/VLQAnalyzer_cfg.py'
config.JobType.pluginName = 'Analysis'
config.JobType.maxJobRuntimeMin = 2500#2000
config.JobType.maxMemoryMB = 4000 #2500
#config.JobType.allowUndistributedCMSSW = True
#config.JobType.disableAutomaticOutputCollection = True

config.section_('Data')
config.Data.inputDataset = None
config.Data.secondaryInputDataset = None
config.Data.inputDBS = 'global'
config.Data.splitting = filesplit
config.Data.unitsPerJob = units
#config.Data.outLFNDirBase = '/store/group/lpcbprime/noreplica/skhalil/Upgrade'
config.Data.outLFNDirBase = '/store/user/skhalil/upgradeNtuples/'
#physical location:  /mnt/hadoop/user/uscms01/pnfs/unl.edu/data4/cms/store/group/skhalil/upgradeNtuples/
config.Data.allowNonValidInputDataset = True
#config.Data.ignoreLocality = False
config.Data.publication = False

config.section_('Site')
config.Site.storageSite = 'T2_US_Nebraska' #'T3_US_FNALLPC'

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
    print 's: ', s
    cdi = s + 'MINIAODSIM'
    cgr = s.split('/')[1]
    #print 'cgr: ', cgr  
    if 'ext1' in s: cgr_v1 = cgr+'_ext1' 
    else: cgr_v1 = cgr 
    r = s.replace("MiniAOD","DR" )
    cds = r + 'GEN-SIM-RECO'
    if 'TT' in cgr:
            cds_v1 = cds.replace("v3", "v1")
    else:   cds_v1 = cds
    config.Data.inputDataset = cdi
    config.General.requestName = cgr_v1
    config.Data.secondaryInputDataset = cds_v1
    print "Submitting to Crab:"
    print "Inputdataset: ",cdi
    print "requestName: ",cgr_v1
    print "SecondaryInputdataset: ",cds_v1
    print 
    submit(config)



#config.Data.inputDataset = '/TprimeBToTH_M-2000_Width-10p_LH_TuneCUETP8M2T4_14TeV-madgraph-pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2-v1/MINIAODSIM'
#config.General.requestName = 'TprimeBToTH_M-2000_Width-10p_LH_TuneCUETP8M2T4_14TeV-madgraph-pythia8'
#config.Data.secondaryInputDataset = '/TprimeBToTH_M-2000_Width-10p_LH_TuneCUETP8M2T4_14TeV-madgraph-pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2-v1/GEN-SIM-RECO'
#submit(config)
