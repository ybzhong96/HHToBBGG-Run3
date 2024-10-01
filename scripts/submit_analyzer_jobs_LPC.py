#!/usr/bin/python

import os
import datetime
import time
import subprocess
import glob
import sys
from collections import OrderedDict

option = 0
label = "option"+str(option)

analysis = "HHTo2B2GNtupler"
outputfile = "HHTo2B2GNtuple" + "_" + label

#analysis = "JetHTTriggerNtupler"
#outputfile = "JetHTTriggerNtuple" + "_" + label

#analysis = "MakeMCPileupDistribution"
#outputfile = "MCPileupDistribution" + "_" + label

cmsswReleaseVersion = "CMSSW_14_0_7"



outputDirectoryBase = "/store/user/yzhong/Run3/HiggsDNA_selection/analyzer/"+analysis+"/"+label+"/"
filesPerJob = 1

datasetList = OrderedDict()

#2022 ntuples
#datasetList["nano/run3/2022/EGamma_2022A.list"] = [1, 1, "2022", "", 1]
#datasetList["nano/run3/2022/EGamma_2022B.list"] = [1, 1, "2022", "", 1]
#datasetList["nano/run3/2022/EGamma_2022C.list"] = [1, 1, "2022", "", 1]
#datasetList["nano/run3/2022/EGamma_2022D.list"] = [1, 1, "2022", "", 1]

datasetList["nano/run3/2022/Run2022B-22Sep2023-v2.list"] = [1, 1, "2022", "", 1]
datasetList["nano/run3/2022/Run2022C-22Sep2023-v1.list"] = [1, 1, "2022", "", 1]
datasetList["nano/run3/2022/Run2022D-22Sep2023-v1.list"] = [1, 1, "2022", "", 1]

datasetList["nano/run3/2022/Run2022E-22Sep2023-v1.list"] = [1, 1, "2022", "", 1]
datasetList["nano/run3/2022/Run2022F-22Sep2023-v1.list"] = [1, 1, "2022", "", 1]
datasetList["nano/run3/2022/Run2022G-22Sep2023-v2.list"] = [1, 1, "2022", "", 1]

#datasetList['nano/run3/2022/JetMET_2022C.list'] = [1, 1, "2022", "", 1]
#datasetList['nano/run3/2022/JetMET_2022D.list'] = [1, 1, "2022", "", 1]
#datasetList['nano/run3/2022/JetMET_2022D-v2.list'] = [1, 1, "2022", "", 10]
#datasetList['nano/run3/2022/JetMET_2022E.list'] = [1, 1, "2022", "", 1]
#datasetList['nano/run3/2022/JetMET_2022F.list'] = [1, 1, "2022", "", 1]
#datasetList['nano/run3/2022/JetMET_2022G.list'] = [1, 1, "2022", "", 1]

datasetList['nano/run3/2022/ggHHTo2B2G_cHHH_1_PrivProd.list'] = [2, 1, "2022", "", 1]
datasetList['nano/run3/2022/GG-Box-3Jets_MGG-0to40_13p6TeV_sherpa_2022.list'] = [2, 1, "2022", "", 1]
datasetList['nano/run3/2022/GG-Box-3Jets_MGG-40to80_13p6TeV_sherpa_2022.list'] = [2, 1, "2022", "", 1]
datasetList['nano/run3/2022/GG-Box-3Jets_MGG-80_13p6TeV_sherpa_2022.list'] = [2, 1, "2022", "", 1]
datasetList['nano/run3/2022/GG-Box-3Jets_MGG-0to40_13p6TeV_sherpa_2022EE.list'] = [2, 1, "2022", "", 1]
datasetList['nano/run3/2022/GG-Box-3Jets_MGG-40to80_13p6TeV_sherpa_2022EE.list'] = [2, 1, "2022", "", 1]
datasetList['nano/run3/2022/GG-Box-3Jets_MGG-80_13p6TeV_sherpa_2022EE.list'] = [2, 1, "2022", "", 1]
datasetList['nano/run3/2022/GJet_PT-10to40_DoubleEMEnriched_TuneCP5_13p6TeV_pythia8_2022.list']=[2,1,"2022","",1]
datasetList['nano/run3/2022/GJet_PT-10to40_DoubleEMEnriched_TuneCP5_13p6TeV_pythia8_2022EE.list']=[2,1,"2022","",1]
datasetList['nano/run3/2022/GJet_PT-40_DoubleEMEnriched_TuneCP5_13p6TeV_pythia8_2022.list']=[2,1,"2022","",1]
datasetList['nano/run3/2022/GJet_PT-40_DoubleEMEnriched_TuneCP5_13p6TeV_pythia8_2022EE.list']=[2,1,"2022","",1]
datasetList['nano/run3/2022/GJet_PT-20to40_DoubleEMEnriched_MGG-80_TuneCP5_13p6TeV_pythia8_2022.list']=[2,1,"2022","",1]
datasetList['nano/run3/2022/GJet_PT-20to40_DoubleEMEnriched_MGG-80_TuneCP5_13p6TeV_pythia8_2022EE.list']=[2,1,"2022","",1]
datasetList['nano/run3/2022/GJet_PT-40_DoubleEMEnriched_MGG-80_TuneCP5_13p6TeV_pythia8_2022.list']=[2,1,"2022","",1]
datasetList['nano/run3/2022/GJet_PT-40_DoubleEMEnriched_MGG-80_TuneCP5_13p6TeV_pythia8_2022EE.list']=[2,1,"2022","",1]

datasetList['nano/run3/2022/GluGluHtoGG_M-120_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_2022EE.list']=[2,1,"2022","",1]
datasetList['nano/run3/2022/GluGluHtoGG_M-125_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_2022EE.list']=[2,1,"2022","",1]
datasetList['nano/run3/2022/GluGluHtoGG_M-130_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_2022EE.list']=[2,1,"2022","",1]
datasetList['nano/run3/2022/VBFHtoGG_M-120_TuneCP5_13p6TeV_amcatnlo-pythia8_2022EE.list']=[2,1,"2022","",1]
datasetList['nano/run3/2022/VBFHtoGG_M-125_TuneCP5_13p6TeV_amcatnlo-pythia8_2022EE.list']=[2,1,"2022","",1]
datasetList['nano/run3/2022/VBFHtoGG_M-130_TuneCP5_13p6TeV_amcatnlo-pythia8_2022EE.list']=[2,1,"2022","",1]
datasetList['nano/run3/2022/VHtoGG_M-120_TuneCP5_13p6TeV_amcatnloFXFX-madspin-pythia8_2022EE.list']=[2,1,"2022","",1]
datasetList['nano/run3/2022/VHtoGG_M-125_TuneCP5_13p6TeV_amcatnloFXFX-madspin-pythia8_2022EE.list']=[2,1,"2022","",1]
datasetList['nano/run3/2022/VHtoGG_M-130_TuneCP5_13p6TeV_amcatnloFXFX-madspin-pythia8_2022EE.list']=[2,1,"2022","",1]
datasetList['nano/run3/2022/ttHtoGG_M-120_TuneCP5_13p6TeV_amcatnloFXFX-madspin-pythia8_2022EE.list']=[2,1,"2022","",1]
datasetList['nano/run3/2022/ttHtoGG_M-125_TuneCP5_13p6TeV_amcatnloFXFX-madspin-pythia8_2022EE.list']=[2,1,"2022","",1]
datasetList['nano/run3/2022/ttHtoGG_M-130_TuneCP5_13p6TeV_amcatnloFXFX-madspin-pythia8_2022EE.list']=[2,1,"2022","",1]
datasetList['nano/run3/2022/GluGlutoHHto2B2G_kl-1p00_kt-1p00_c2-0p00_TuneCP5_13p6TeV_powheg-pythia8_2022EE.list']=[2,1,"2022","",1]
#datasetList['nano/run3/2022/QCDPt80to120_postEE.list']=[2,1,"2022","",1]
#datasetList['nano/run3/2022/QCDPt120to170_postEE.list']=[2,1,"2022","",1]
#datasetList['nano/run3/2022/QCDPt170to300_postEE.list']=[2,1,"2022","",1]
#datasetList['nano/run3/2022/QCDPt300to470_postEE.list']=[2,1,"2022","",1]
#datasetList['nano/run3/2022/QCDPt470to600_postEE.list']=[2,1,"2022","",1]
#datasetList['nano/run3/2022/QCDPt600to800_postEE.list']=[2,1,"2022","",1]
#datasetList['nano/run3/2022/QCDPt800to1000_postEE.list']=[2,1,"2022","",1]
#datasetList['nano/run3/2022/QCDPt1000to1400_postEE.list']=[2,1,"2022","",1]
#datasetList['nano/run3/2022/QCDPt1400to1800_postEE.list']=[2,1,"2022","",1]
#datasetList['nano/run3/2022/QCDPt1800to2400_postEE.list']=[2,1,"2022","",1]
#datasetList['nano/run3/2022/QCDPt2400to3200_postEE.list']=[2,1,"2022","",1]

CMSSW_BASE_DIR = os.getenv('CMSSW_BASE')
Analyzer_DIR = CMSSW_BASE_DIR+"/src/HHToBBGG-Run3/"
print(Analyzer_DIR)
#create directory for condor jobs

for listfile in datasetList.keys():

    datasetName = listfile.replace(".list","")
    print ("Preparing analyzer workflow for dataset :" + datasetName + "\n")
    if not os.path.exists(Analyzer_DIR+"/list/" + listfile):
        print ("listfile: " + listfile + " does not exist. skipping.")
        continue

    outputDirectory = outputDirectoryBase + datasetName + "/"
    tmpListFile = open(Analyzer_DIR + "/list/" + listfile,"r")

    year = datasetList[listfile][2]
    sampleName = datasetList[listfile][3]
    numberOfJobsPerFile = datasetList[listfile][4]

    #####################################
    #Job Splitting
    #####################################
    isData = "no"
    if (datasetList[listfile][0] == 1): 
        isData = "yes"
    filesPerJob = datasetList[listfile][1]
    tmpJobFileCount = 0
    nJobs = numberOfJobsPerFile

    if os.path.exists(Analyzer_DIR+"/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/" + "input_list_" + str(nJobs) + ".txt"):
        print ("Warning: condor directory " + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + " is not empty. Skipping.")
        continue
        
    #create condor directories
    os.system("mkdir -p " + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName )
    os.system("mkdir -p " + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/log/")
    os.system("mkdir -p " + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/out/")
    os.system("mkdir -p " + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/err/")

    tmpOutputListFile = open( Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/" + "input_list_" + str(nJobs/numberOfJobsPerFile) + ".txt","w")
    for line in tmpListFile:
                
        #open list file for new job
        if tmpJobFileCount >= filesPerJob:
            tmpOutputListFile.close()
            tmpJobFileCount = 0
            nJobs = nJobs + numberOfJobsPerFile
            tmpOutputListFile = open( Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/" + "input_list_" + str(nJobs/numberOfJobsPerFile) + ".txt","w")
          
        #write input file into job list file
        tmpOutputListFile.write(line)
        tmpJobFileCount += 1

    tmpOutputListFile.close()    
    os.system("cd " + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/; tar czf input_list.tgz input_list_*.txt")

    ###################################################
    #Copy code tarball, run script, and executable
    ###################################################
    #os.system("cd " + Analyzer_DIR + "; tar vczf " + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/code.tgz " + " include/ src/ app/ Makefile" )
    os.system("cp " + Analyzer_DIR + "/scripts/run_job_LPC.sh " + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/")
    os.system("cp " + Analyzer_DIR + "Run" + analysis + " " + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/")
    os.system("mkdir -p " + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/HHToBBGG-Run3/data/")
    os.system("cp " + Analyzer_DIR + "/data/JetHTTriggerEfficiency_2016" + ".root " + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/HHToBBGG-Run3/data/")
    os.system("cp " + Analyzer_DIR + "/data/JetHTTriggerEfficiency_2017" + ".root " + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/HHToBBGG-Run3/data/")
    os.system("cp " + Analyzer_DIR + "/data/JetHTTriggerEfficiency_2018" + ".root " + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/HHToBBGG-Run3/data/")
    #os.system("cp " + Analyzer_DIR + "/data/JetHTTriggerEfficiency_" + year + ".root " + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/HHToBBGG-Run3/data/")
    os.system("cp " + Analyzer_DIR + "/data/JetHTTriggerEfficiency_Summer16.root " + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/HHToBBGG-Run3/data/")
    os.system("cp " + Analyzer_DIR + "/data/JetHTTriggerEfficiency_Fall17.root " + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/HHToBBGG-Run3/data/")
    os.system("cp " + Analyzer_DIR + "/data/JetHTTriggerEfficiency_Fall18.root " + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/HHToBBGG-Run3/data/")

    os.system("cp " + Analyzer_DIR + "/Run3_2022_2023_Golden.json " + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/HHToBBGG-Run3/data/")

    os.system("mkdir -p " + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/HHToBBGG-Run3/data/PileupWeights/")
    os.system("cp " + Analyzer_DIR + "/data/PileupWeights/PileupWeights.root " + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/HHToBBGG-Run3/data/PileupWeights/")
    os.system("mkdir -p " + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/HHToBBGG-Run3/data/JEC/Summer16_07Aug2017_V11_MC/")
    os.system("cp " + Analyzer_DIR + "/data/JEC/Summer16_07Aug2017_V11_MC/Summer16_07Aug2017_V11_MC_Uncertainty_AK8PFPuppi.txt " + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/HHToBBGG-Run3/data/JEC/Summer16_07Aug2017_V11_MC/")
    os.system("mkdir -p " + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/HHToBBGG-Run3/data/JEC/Fall17_17Nov2017_V32_MC/")
    os.system("cp " + Analyzer_DIR + "/data/JEC/Fall17_17Nov2017_V32_MC/Fall17_17Nov2017_V32_MC_Uncertainty_AK8PFPuppi.txt " + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/HHToBBGG-Run3/data/JEC/Fall17_17Nov2017_V32_MC/")
    os.system("mkdir -p " + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/HHToBBGG-Run3/data/JEC/Autumn18_V19_MC/")
    os.system("cp " + Analyzer_DIR + "/data/JEC/Autumn18_V19_MC/Autumn18_V19_MC_Uncertainty_AK8PFPuppi.txt " + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/HHToBBGG-Run3/data/JEC/Autumn18_V19_MC/")
   os.system("mkdir -p " + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/HHToBBGG-Run3/data/JEC/Summer22_22Sep2023_RunCD_V2_DATA/")
    os.system("cp " + Analyzer_DIR + "/data/JEC/Summer22_22Sep2023_RunCD_V2_DATA/Summer22_22Sep2023_RunCD_V2_DATA_L1FastJet_AK4PFchs.txt " + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/HHToBBGG-Run3/data/JEC/Summer22_22Sep2023_RunCD_V2_DATA/")


    #####################################
    #Create Condor JDL file
    #####################################
    tmpCondorJDLFile = open(Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/task.jdl","w+")
    tmpCondorJDLFileTemplate = """
Universe  = vanilla
Executable = ./run_job_LPC.sh
"""
    tmpCondorJDLFile.write(tmpCondorJDLFileTemplate)
    tmpCondorJDLFile.write("Arguments = " + analysis + " " + str(isData) + " " + str(option) + " " + "$(I) " + str(numberOfJobsPerFile) + " " + outputfile + " " + outputDirectory + " " + cmsswReleaseVersion + " " + year + " " + sampleName + "\n")

    tmpCondorJDLFileTemplate = """
Log = log/job.$(Cluster).$(Process).log
Output = out/job.$(Cluster).$(Process).out
Error = err/job.$(Cluster).$(Process).err
x509userproxy = $ENV(X509_USER_PROXY)
"""
    tmpCondorJDLFile.write(tmpCondorJDLFileTemplate)
    tmpCondorJDLFile.write(
        "transfer_input_files = " 
        + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/run_job_LPC.sh, " 
        + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/input_list.tgz, " 
        + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/Run" + analysis + ", "
        #+ Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/code.tgz, " 
        + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/HHToBBGG-Run3/data/JetHTTriggerEfficiency_2016" + ".root, "
        + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/HHToBBGG-Run3/data/JetHTTriggerEfficiency_2017" + ".root, "
        + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/HHToBBGG-Run3/data/JetHTTriggerEfficiency_2018" + ".root, "
        + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/HHToBBGG-Run3/data/Run3_2022_2023_Golden.json, "  # json file
        + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/HHToBBGG-Run3/data/JetHTTriggerEfficiency_Summer16.root, " 
        + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/HHToBBGG-Run3/data/JetHTTriggerEfficiency_Fall17.root, " 
        + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/HHToBBGG-Run3/data/JetHTTriggerEfficiency_Fall18.root, " 
        + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/HHToBBGG-Run3/data/PileupWeights/PileupWeights.root, " 
        + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/HHToBBGG-Run3/data/JEC/Summer16_07Aug2017_V11_MC/Summer16_07Aug2017_V11_MC_Uncertainty_AK8PFPuppi.txt, " 
        + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/HHToBBGG-Run3/data/JEC/Fall17_17Nov2017_V32_MC/Fall17_17Nov2017_V32_MC_Uncertainty_AK8PFPuppi.txt, " 
        + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/HHToBBGG-Run3/data/JEC/Autumn18_V19_MC/Autumn18_V19_MC_Uncertainty_AK8PFPuppi.txt " + "\n"
        + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/HHToBBGG-Run3/data/JEC/Summer22_22Sep2023_RunCD_V2_DATA/Summer22_22Sep2023_RunCD_V2_DATA_L1FastJet_AK4PFchs.txt " + "\n"
    )

    tmpCondorJDLFileTemplate = """
should_transfer_files = YES
when_to_transfer_output = ON_EXIT

# Resources request
"""
    tmpCondorJDLFile.write(tmpCondorJDLFileTemplate)
    tmpCondorJDLFile.write("RequestMemory = 2000 \n")

    tmpCondorJDLFileTemplate = """

# Jobs selection
Queue I from (
"""

    tmpCondorJDLFile.write(tmpCondorJDLFileTemplate)
    for i in range(1,nJobs+1):
        tmpCondorJDLFile.write(str(i)+"\n")
    tmpCondorJDLFile.write(")\n")
    tmpCondorJDLFile.close()





