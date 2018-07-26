# VLQAnaAnalyzer
Instructions to check out code:
Set up CMSSW area:
```
setenv SCRAM_ARCH slc6_amd64_gcc630
cmsrel CMSSW_9_3_2
cd CMSSW_9_3_2/src
cmsenv
git cms-init
git cms-merge-topic -u nsmith-:EgammaFromMultiCl_932v2
mkdir -p RecoEgamma && pushd RecoEgamma
git clone -b integrated https://github.com/nsmith-/Phase2InterimID.git
popd
git clone git@github.com:skhail/VLQAnalyzer.git Upgrades/VLQAnalyzer
scram b -j4
```
To run:
cd test
python testVLQAnalyzer.py -f "NAME-OF-ROOT-FILE-YOU-WOULD-LIKE-TO-RUN-OVER.root"



Note: To run over individual ntuple files, you must use the testVLQAnalyzer.py file and not the VLQAnalyzer.py file, because the VLQAnalyzer.py file is set up to be used in CRAB jobs.
