from WMCore.Configuration import Configuration
config = Configuration()
config.section_('General')
config.General.requestName = 'VLQAnalysis_histos'
config.section_('JobType')
config.JobType.psetName = '/home/t3-ku/z184o935/CMSSW_9_3_2/src/Upgrades/VLQAnalyzer/test/do_nothing_cfg.py'
config.JobType.pluginName = 'PrivateMC'
config.JobType.inputFiles = ['/home/t3-ku/z184o935/CMSSW_9_3_2/src/Upgrades/VLQAnalyzer/test/FrameworkJobReport.xml','/home/t3-ku/z184o935/CMSSW_9_3_2/src/Upgrades/VLQAnalyzer/test/VLQAnalyzer.py','/home/t3-ku/z184o935/CMSSW_9_3_2/src/Upgrades/VLQAnalyzer/test/VLQAnalysis_Ntuples.txt']
config.JobType.outputFiles = ['VLQAnalysis_output.tar']
config.JobType.scriptExe = 'myscript.sh'
config.section_('Data')
config.Data.outputDatasetTag = 'VLQAnalysis_histos'
config.Data.publication = False
config.Data.unitsPerJob = 1
config.Data.splitting = 'EventBased'
config.Data.outLFNDirBase = '/store/user/zolson/VLQ/'
config.Data.outputPrimaryDataset = 'Combine'
config.Data.totalUnits = 29
config.section_('User')
config.section_('Site')
config.Site.blacklist = ['T3_IT_Bologna', 'T3_US_UMiss', 'T2_RU_PNPI']
config.Site.storageSite = 'T2_US_Nebraska'
