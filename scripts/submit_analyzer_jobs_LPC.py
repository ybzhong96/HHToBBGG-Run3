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


outputDirectoryBase = "/store/group/lpcdihiggsboost/sixie/analyzer/"+analysis+"/"+label+"/"
filesPerJob = 1

datasetList = OrderedDict()

#Make line from list of inputs
#ls | awk '{print "datasetList[\"nano/run3/2022EE/"$1"\"] = [1, 1, \"2022EE\", \"\", 1]"}'

########################################################
#2022 ntuples
########################################################






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
    os.system("cp " + Analyzer_DIR + "/data/Run3_2022_2023_Golden.json " + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/HHToBBGG-Run3/data/")
    os.system("mkdir -p " + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/HHToBBGG-Run3/data/PileupWeights/")
    os.system("cp " + Analyzer_DIR + "/data/PileupWeights/PileupReweight_Summer22.root " + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/HHToBBGG-Run3/data/PileupWeights/")
    os.system("cp " + Analyzer_DIR + "/data/PileupWeights/PileupReweight_Summer22EE.root " + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/HHToBBGG-Run3/data/PileupWeights/")
    os.system("cp " + Analyzer_DIR + "/data/PileupWeights/PileupReweight_Summer23.root " + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/HHToBBGG-Run3/data/PileupWeights/")
    os.system("cp " + Analyzer_DIR + "/data/PileupWeights/PileupReweight_Summer23BPix.root " + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/HHToBBGG-Run3/data/PileupWeights/")
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
        + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/HHToBBGG-Run3/data/Run3_2022_2023_Golden.json, "  # json file       
        + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/HHToBBGG-Run3/data/JEC/Summer16_07Aug2017_V11_MC/Summer16_07Aug2017_V11_MC_Uncertainty_AK8PFPuppi.txt, " 
        + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/HHToBBGG-Run3/data/JEC/Fall17_17Nov2017_V32_MC/Fall17_17Nov2017_V32_MC_Uncertainty_AK8PFPuppi.txt, " 
        + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/HHToBBGG-Run3/data/JEC/Autumn18_V19_MC/Autumn18_V19_MC_Uncertainty_AK8PFPuppi.txt, " 
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





