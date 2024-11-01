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

########################################################
#2022 ntuples
########################################################
datasetList["nano/run3/2022/EGamma_2022C.list"] = [1, 1, "2022", "", 1]
datasetList["nano/run3/2022/EGamma_2022D.list"] = [1, 1, "2022", "", 1]
datasetList['nano/run3/2022/GluGlutoHHto2B2G_kl-1p00_kt-1p00_c2-0p00_TuneCP5_13p6TeV_powheg-pythia8.list'] = [2, 1, "2022", "", 1]
datasetList['nano/run3/2022/GluGlutoHHto2B2G_kl-0p00_kt-1p00_c2-0p00_TuneCP5_13p6TeV_powheg-pythia8.list'] = [2, 1, "2022", "", 1]
datasetList['nano/run3/2022/GluGlutoHHto2B2G_kl-5p00_kt-1p00_c2-0p00_TuneCP5_13p6TeV_powheg-pythia8.list'] = [2, 1, "2022", "", 1]
datasetList['nano/run3/2022/GluGlutoHHto2B2G_kl-2p45_kt-1p00_c2-0p00_LHEweights_TuneCP5_13p6TeV_powheg-pythia.list'] = [2, 1, "2022", "", 1]
datasetList['nano/run3/2022/VBFHHto2B2G_CV-1p74_C2V-1p37_C3-14p4_TuneCP5_13p6TeV_madgraph-pythia8.list'] = [2, 1, "2022", "", 1]
datasetList['nano/run3/2022/VBFHHto2B2G_CV-m0p758_C2V-1p44_C3-m19p3_TuneCP5_13p6TeV_madgraph-pythia8.list'] = [2, 1, "2022", "", 1]
datasetList['nano/run3/2022/VBFHHto2B2G_CV-m0p012_C2V-0p030_C3-10p2_TuneCP5_13p6TeV_madgraph-pythia8.list'] = [2, 1, "2022", "", 1]
datasetList['nano/run3/2022/VBFHHto2B2G_CV-m0p962_C2V-0p959_C3-m1p43_TuneCP5_13p6TeV_madgraph-pythia8.list'] = [2, 1, "2022", "", 1]
datasetList['nano/run3/2022/VBFHHto2B2G_CV-m1p60_C2V-2p72_C3-m1p36_TuneCP5_13p6TeV_madgraph-pythia8.list'] = [2, 1, "2022", "", 1]
datasetList['nano/run3/2022/VBFHHto2B2G_CV-m1p21_C2V-1p94_C3-m0p94_TuneCP5_13p6TeV_madgraph-pythia8.list'] = [2, 1, "2022", "", 1]
datasetList['nano/run3/2022/VBFHHto2B2G_CV_1_C2V_0_C3_1_TuneCP5_13p6TeV_madgraph-pythia8.list'] = [2, 1, "2022", "", 1]
datasetList['nano/run3/2022/VBFHHto2B2G_CV-m2p12_C2V-3p87_C3-m5p96_TuneCP5_13p6TeV_madgraph-pythia8.list'] = [2, 1, "2022", "", 1]
datasetList['nano/run3/2022/VBFHHto2B2G_CV_1_C2V_1_C3_1_TuneCP5_13p6TeV_madgraph-pythia8.list'] = [2, 1, "2022", "", 1]
datasetList['nano/run3/2022/WHH_HHto2B2G_CV-1p0_C2V-0p0_C3-1p0_TuneCP5_13p6TeV_madgraph-pythia8.list'] = [2, 1, "2022", "", 1]
datasetList['nano/run3/2022/WHH_HHto2B2G_CV-0p5_C2V-1p0_C3-1p0_TuneCP5_13p6TeV_madgraph-pythia8.list'] = [2, 1, "2022", "", 1]
datasetList['nano/run3/2022/WHH_HHto2B2G_CV-1p0_C2V-1p0_C3-0p0_TuneCP5_13p6TeV_madgraph-pythia8.list'] = [2, 1, "2022", "", 1]
datasetList['nano/run3/2022/WHH_HHto2B2G_CV-1p0_C2V-10p0_C3-1p0_TuneCP5_13p6TeV_madgraph-pythia8.list'] = [2, 1, "2022", "", 1]
datasetList['nano/run3/2022/WHH_HHto2B2G_CV-1p0_C2V-1p0_C3-1p0_TuneCP5_13p6TeV_madgraph-pythia8.list'] = [2, 1, "2022", "", 1]
datasetList['nano/run3/2022/WHH_HHto2B2G_CV-1p0_C2V-1p0_C3-2p0_TuneCP5_13p6TeV_madgraph-pythia8.list'] = [2, 1, "2022", "", 1]
datasetList['nano/run3/2022/WHH_HHto2B2G_CV-1p0_C2V-1p0_C3-20p0_TuneCP5_13p6TeV_madgraph-pythia8.list'] = [2, 1, "2022", "", 1]
datasetList['nano/run3/2022/WHH_HHto2B2G_CV-1p0_C2V-2p0_C3-1p0_TuneCP5_13p6TeV_madgraph-pythia8.list'] = [2, 1, "2022", "", 1]
datasetList['nano/run3/2022/ZHH_HHto2B2G_CV-0p5_C2V-1p0_C3-1p0_TuneCP5_13p6TeV_madgraph-pythia8.list'] = [2, 1, "2022", "", 1]
datasetList['nano/run3/2022/WHH_HHto2B2G_CV-1p5_C2V-1p0_C3-1p0_TuneCP5_13p6TeV_madgraph-pythia8.list'] = [2, 1, "2022", "", 1]
datasetList['nano/run3/2022/ZHH_HHto2B2G_CV-1p0_C2V-10p0_C3-1p0_TuneCP5_13p6TeV_madgraph-pythia8.list'] = [2, 1, "2022", "", 1]
datasetList['nano/run3/2022/ZHH_HHto2B2G_CV-1p0_C2V-0p0_C3-1p0_TuneCP5_13p6TeV_madgraph-pythia8.list'] = [2, 1, "2022", "", 1]
datasetList['nano/run3/2022/ZHH_HHto2B2G_CV-1p0_C2V-1p0_C3-0p0_TuneCP5_13p6TeV_madgraph-pythia8.list'] = [2, 1, "2022", "", 1]
datasetList['nano/run3/2022/ZHH_HHto2B2G_CV-1p0_C2V-1p0_C3-20p0_TuneCP5_13p6TeV_madgraph-pythia8.list'] = [2, 1, "2022", "", 1]
datasetList['nano/run3/2022/ZHH_HHto2B2G_CV-1p0_C2V-1p0_C3-1p0_TuneCP5_13p6TeV_madgraph-pythia8.list'] = [2, 1, "2022", "", 1]
datasetList['nano/run3/2022/ZHH_HHto2B2G_CV-1p0_C2V-1p0_C3-2p0_TuneCP5_13p6TeV_madgraph-pythia8.list'] = [2, 1, "2022", "", 1]
datasetList['nano/run3/2022/ZHH_HHto2B2G_CV-1p5_C2V-1p0_C3-1p0_TuneCP5_13p6TeV_madgraph-pythia8.list'] = [2, 1, "2022", "", 1]
datasetList['nano/run3/2022/ZHH_HHto2B2G_CV-1p0_C2V-2p0_C3-1p0_TuneCP5_13p6TeV_madgraph-pythia8.list'] = [2, 1, "2022", "", 1]

datasetList['nano/run3/2022/GluGluHtoGG_M-125_TuneCP5_13p6TeV_amcatnloFXFX-pythia8.list'] = [2, 1, "2022", "", 1]
datasetList['nano/run3/2022/GluGluHToGG_M-125_TuneCP5_13p6TeV_powheg-pythia8.list'] = [2, 1, "2022", "", 1]
datasetList['nano/run3/2022/GG-Box-3Jets_MGG-0to40_13p6TeV_sherpa_2022.list'] = [2, 1, "2022", "", 1]
datasetList['nano/run3/2022/GG-Box-3Jets_MGG-40to80_13p6TeV_sherpa_2022.list'] = [2, 1, "2022", "", 1]
datasetList['nano/run3/2022/GG-Box-3Jets_MGG-80_13p6TeV_sherpa_2022.list'] = [2, 1, "2022", "", 1]
datasetList['nano/run3/2022/GJet_PT-10to40_DoubleEMEnriched_TuneCP5_13p6TeV_pythia8_2022.list']=[2,1,"2022","",1]
datasetList['nano/run3/2022/GJet_PT-40_DoubleEMEnriched_TuneCP5_13p6TeV_pythia8_2022.list']=[2,1,"2022","",1]
datasetList['nano/run3/2022/GJet_PT-20to40_DoubleEMEnriched_MGG-80_TuneCP5_13p6TeV_pythia8_2022.list']=[2,1,"2022","",1]
datasetList['nano/run3/2022/GJet_PT-40_DoubleEMEnriched_MGG-80_TuneCP5_13p6TeV_pythia8_2022.list']=[2,1,"2022","",1]
datasetList['nano/run3/2022/QCD_PT-30to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13p6TeV_pythia8.list'] = [2, 1, "2022", "", 1]
datasetList['nano/run3/2022/QCD_PT-30toInf_DoubleEMEnriched_MGG-40to80_TuneCP5_13p6TeV_pythia8.list'] = [2, 1, "2022", "", 1]
datasetList['nano/run3/2022/QCD_PT-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13p6TeV_pythia8.list'] = [2, 1, "2022", "", 1]
datasetList['nano/run3/2022/TTG-1Jets_PTG-100to200_TuneCP5_13p6TeV_amcatnloFXFXold-pythia8.list'] = [2, 1, "2022", "", 1]
datasetList['nano/run3/2022/TJGG_TuneCP5_13p6TeV_amcatnlo-madspin-pythia8.list'] = [2, 1, "2022", "", 1]
datasetList['nano/run3/2022/TTG-1Jets_PTG-200_TuneCP5_13p6TeV_amcatnloFXFXold-pythia8.list'] = [2, 1, "2022", "", 1]
datasetList['nano/run3/2022/TTG-1Jets_PTG-10to100_TuneCP5_13p6TeV_amcatnloFXFXold-pythia8.list'] = [2, 1, "2022", "", 1]
datasetList['nano/run3/2022/TTGG_TuneCP5_13p6TeV_madgraph-madspin-pythia8.list'] = [2, 1, "2022", "", 1]


########################################################
#2022EE ntuples
########################################################
datasetList["nano/run3/2022/EGamma_2022E.list"] = [1, 1, "2022EE", "", 1]
datasetList["nano/run3/2022/EGamma_2022F.list"] = [1, 1, "2022EE", "", 1]
datasetList["nano/run3/2022/EGamma_2022G.list"] = [1, 1, "2022EE", "", 1]

datasetList["nano/run3/2022EE/GluGlutoHHto2B2G_kl-0p00_kt-1p00_c2-0p00_TuneCP5_13p6TeV_powheg-pythia8.list"] = [1, 1, "2022EE", "", 1]
datasetList["nano/run3/2022EE/GluGlutoHHto2B2G_kl-1p00_kt-1p00_c2-0p00_TuneCP5_13p6TeV_powheg-pythia8.list"] = [1, 1, "2022EE", "", 1]
datasetList["nano/run3/2022EE/GluGlutoHHto2B2G_kl-2p45_kt-1p00_c2-0p00_LHEweights_TuneCP5_13p6TeV_powheg-pythia.list"] = [1, 1, "2022EE", "", 1]
datasetList["nano/run3/2022EE/GluGlutoHHto2B2G_kl-5p00_kt-1p00_c2-0p00_TuneCP5_13p6TeV_powheg-pythia8.list"] = [1, 1, "2022EE", "", 1]
datasetList["nano/run3/2022EE/VBFHHto2B2G_CV-1p74_C2V-1p37_C3-14p4_TuneCP5_13p6TeV_madgraph-pythia8.list"] = [1, 1, "2022EE", "", 1]
datasetList["nano/run3/2022EE/VBFHHto2B2G_CV-m0p012_C2V-0p030_C3-10p2_TuneCP5_13p6TeV_madgraph-pythia8.list"] = [1, 1, "2022EE", "", 1]
datasetList["nano/run3/2022EE/VBFHHto2B2G_CV-m0p758_C2V-1p44_C3-m19p3_TuneCP5_13p6TeV_madgraph-pythia8.list"] = [1, 1, "2022EE", "", 1]
datasetList["nano/run3/2022EE/VBFHHto2B2G_CV-m0p962_C2V-0p959_C3-m1p43_TuneCP5_13p6TeV_madgraph-pythia8.list"] = [1, 1, "2022EE", "", 1]
datasetList["nano/run3/2022EE/VBFHHto2B2G_CV-m1p21_C2V-1p94_C3-m0p94_TuneCP5_13p6TeV_madgraph-pythia8.list"] = [1, 1, "2022EE", "", 1]
datasetList["nano/run3/2022EE/VBFHHto2B2G_CV-m1p60_C2V-2p72_C3-m1p36_TuneCP5_13p6TeV_madgraph-pythia8.list"] = [1, 1, "2022EE", "", 1]
datasetList["nano/run3/2022EE/VBFHHto2B2G_CV-m1p83_C2V-3p57_C3-m3p39_TuneCP5_13p6TeV_madgraph-pythia8.list"] = [1, 1, "2022EE", "", 1]
datasetList["nano/run3/2022EE/VBFHHto2B2G_CV-m2p12_C2V-3p87_C3-m5p96_TuneCP5_13p6TeV_madgraph-pythia8.list"] = [1, 1, "2022EE", "", 1]
datasetList["nano/run3/2022EE/VBFHHto2B2G_CV_1_C2V_0_C3_1_TuneCP5_13p6TeV_madgraph-pythia8.list"] = [1, 1, "2022EE", "", 1]
datasetList["nano/run3/2022EE/VBFHHto2B2G_CV_1_C2V_1_C3_1_TuneCP5_13p6TeV_madgraph-pythia8.list"] = [1, 1, "2022EE", "", 1]
datasetList["nano/run3/2022EE/WHH_HHto2B2G_CV-0p5_C2V-1p0_C3-1p0_TuneCP5_13p6TeV_madgraph-pythia8.list"] = [1, 1, "2022EE", "", 1]
datasetList["nano/run3/2022EE/WHH_HHto2B2G_CV-1p0_C2V-0p0_C3-1p0_TuneCP5_13p6TeV_madgraph-pythia8.list"] = [1, 1, "2022EE", "", 1]
datasetList["nano/run3/2022EE/WHH_HHto2B2G_CV-1p0_C2V-10p0_C3-1p0_TuneCP5_13p6TeV_madgraph-pythia8.list"] = [1, 1, "2022EE", "", 1]
datasetList["nano/run3/2022EE/WHH_HHto2B2G_CV-1p0_C2V-1p0_C3-0p0_TuneCP5_13p6TeV_madgraph-pythia8.list"] = [1, 1, "2022EE", "", 1]
datasetList["nano/run3/2022EE/WHH_HHto2B2G_CV-1p0_C2V-1p0_C3-1p0_TuneCP5_13p6TeV_madgraph-pythia8.list"] = [1, 1, "2022EE", "", 1]
datasetList["nano/run3/2022EE/WHH_HHto2B2G_CV-1p0_C2V-1p0_C3-20p0_TuneCP5_13p6TeV_madgraph-pythia8.list"] = [1, 1, "2022EE", "", 1]
datasetList["nano/run3/2022EE/WHH_HHto2B2G_CV-1p0_C2V-1p0_C3-2p0_TuneCP5_13p6TeV_madgraph-pythia8.list"] = [1, 1, "2022EE", "", 1]
datasetList["nano/run3/2022EE/WHH_HHto2B2G_CV-1p0_C2V-2p0_C3-1p0_TuneCP5_13p6TeV_madgraph-pythia8.list"] = [1, 1, "2022EE", "", 1]
datasetList["nano/run3/2022EE/WHH_HHto2B2G_CV-1p5_C2V-1p0_C3-1p0_TuneCP5_13p6TeV_madgraph-pythia8.list"] = [1, 1, "2022EE", "", 1]
datasetList["nano/run3/2022EE/ZHH_HHto2B2G_CV-0p5_C2V-1p0_C3-1p0_TuneCP5_13p6TeV_madgraph-pythia8.list"] = [1, 1, "2022EE", "", 1]
datasetList["nano/run3/2022EE/ZHH_HHto2B2G_CV-1p0_C2V-0p0_C3-1p0_TuneCP5_13p6TeV_madgraph-pythia8.list"] = [1, 1, "2022EE", "", 1]
datasetList["nano/run3/2022EE/ZHH_HHto2B2G_CV-1p0_C2V-10p0_C3-1p0_TuneCP5_13p6TeV_madgraph-pythia8.list"] = [1, 1, "2022EE", "", 1]
datasetList["nano/run3/2022EE/ZHH_HHto2B2G_CV-1p0_C2V-1p0_C3-0p0_TuneCP5_13p6TeV_madgraph-pythia8.list"] = [1, 1, "2022EE", "", 1]
datasetList["nano/run3/2022EE/ZHH_HHto2B2G_CV-1p0_C2V-1p0_C3-1p0_TuneCP5_13p6TeV_madgraph-pythia8.list"] = [1, 1, "2022EE", "", 1]
datasetList["nano/run3/2022EE/ZHH_HHto2B2G_CV-1p0_C2V-1p0_C3-20p0_TuneCP5_13p6TeV_madgraph-pythia8.list"] = [1, 1, "2022EE", "", 1]
datasetList["nano/run3/2022EE/ZHH_HHto2B2G_CV-1p0_C2V-1p0_C3-2p0_TuneCP5_13p6TeV_madgraph-pythia8.list"] = [1, 1, "2022EE", "", 1]
datasetList["nano/run3/2022EE/ZHH_HHto2B2G_CV-1p0_C2V-2p0_C3-1p0_TuneCP5_13p6TeV_madgraph-pythia8.list"] = [1, 1, "2022EE", "", 1]
datasetList["nano/run3/2022EE/ZHH_HHto2B2G_CV-1p5_C2V-1p0_C3-1p0_TuneCP5_13p6TeV_madgraph-pythia8.list"] = [1, 1, "2022EE", "", 1]

datasetList["nano/run3/2022EE/GluGluHToGG_M-125_TuneCP5_13p6TeV_powheg-pythia8.list"] = [1, 1, "2022EE", "", 1]
datasetList["nano/run3/2022EE/GluGluHtoGG_M-125_TuneCP5_13p6TeV_amcatnloFXFX-pythia8.list"] = [1, 1, "2022EE", "", 1]
datasetList["nano/run3/2022EE/VBFHtoGG_M-125_TuneCP5_13p6TeV_amcatnlo-pythia8.list"] = [1, 1, "2022EE", "", 1]
datasetList["nano/run3/2022EE/VHtoGG_M-125_TuneCP5_13p6TeV_amcatnloFXFX-madspin-pythia8.list"] = [1, 1, "2022EE", "", 1]
datasetList["nano/run3/2022EE/ttHtoGG_M-125_TuneCP5_13p6TeV_amcatnloFXFX-madspin-pythia8.list"] = [1, 1, "2022EE", "", 1]
datasetList["nano/run3/2022EE/GG-Box-3Jets_MGG-0to40_13p6TeV_sherpa.list"] = [1, 1, "2022EE", "", 1]
datasetList["nano/run3/2022EE/GG-Box-3Jets_MGG-40to80_13p6TeV_sherpa.list"] = [1, 1, "2022EE", "", 1]
datasetList["nano/run3/2022EE/GG-Box-3Jets_MGG-80_13p6TeV_sherpa.list"] = [1, 1, "2022EE", "", 1]
datasetList["nano/run3/2022EE/GJet_PT-10to40_DoubleEMEnriched_TuneCP5_13p6TeV_pythia8.list"] = [1, 1, "2022EE", "", 1]
datasetList["nano/run3/2022EE/GJet_PT-20_DoubleEMEnriched_MGG-40to80_TuneCP5_13p6TeV_pythia8.list"] = [1, 1, "2022EE", "", 1]
datasetList["nano/run3/2022EE/GJet_PT-20to40_DoubleEMEnriched_MGG-80_TuneCP5_13p6TeV_pythia8.list"] = [1, 1, "2022EE", "", 1]
datasetList["nano/run3/2022EE/GJet_PT-40_DoubleEMEnriched_MGG-80_TuneCP5_13p6TeV_pythia8.list"] = [1, 1, "2022EE", "", 1]
datasetList["nano/run3/2022EE/GJet_PT-40_DoubleEMEnriched_TuneCP5_13p6TeV_pythia8.list"] = [1, 1, "2022EE", "", 1]
datasetList["nano/run3/2022EE/QCD_PT-30to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13p6TeV_pythia8.list"] = [1, 1, "2022EE", "", 1]
datasetList["nano/run3/2022EE/QCD_PT-30toInf_DoubleEMEnriched_MGG-40to80_TuneCP5_13p6TeV_pythia8.list"] = [1, 1, "2022EE", "", 1]
datasetList["nano/run3/2022EE/QCD_PT-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13p6TeV_pythia8.list"] = [1, 1, "2022EE", "", 1]
datasetList["nano/run3/2022EE/TJGG_TuneCP5_13p6TeV_amcatnlo-madspin-pythia8.list"] = [1, 1, "2022EE", "", 1]
datasetList["nano/run3/2022EE/TTG-1Jets_PTG-100to200_TuneCP5_13p6TeV_amcatnloFXFXold-pythia8.list"] = [1, 1, "2022EE", "", 1]
datasetList["nano/run3/2022EE/TTG-1Jets_PTG-10to100_TuneCP5_13p6TeV_amcatnloFXFXold-pythia8.list"] = [1, 1, "2022EE", "", 1]
datasetList["nano/run3/2022EE/TTG-1Jets_PTG-200_TuneCP5_13p6TeV_amcatnloFXFXold-pythia8.list"] = [1, 1, "2022EE", "", 1]
datasetList["nano/run3/2022EE/TTGG_TuneCP5_13p6TeV_madgraph-madspin-pythia8.list"] = [1, 1, "2022EE", "", 1]



########################################################
#2023 ntuples
########################################################
datasetList['nano/run3/2023/EGamma_2023C_v1.list'] = [1, 1, "2023", "", 1]
datasetList['nano/run3/2023/EGamma_2023C_v2.list'] = [1, 1, "2023", "", 1]
datasetList['nano/run3/2023/EGamma_2023C_v3.list'] = [1, 1, "2023", "", 1]
datasetList['nano/run3/2023/EGamma_2023C_v4.list'] = [1, 1, "2023", "", 1]



########################################################
#2023BPix ntuples
########################################################
datasetList['nano/run3/2023/EGamma_2023D_v1.list'] = [1, 1, "2023BPix", "", 1]
datasetList['nano/run3/2023/EGamma_2023D_v2.list'] = [1, 1, "2023BPix", "", 1]


########################################################
#2024 ntuples
########################################################
datasetList['nano/run3/2024/EGamma_2024B.list'] = [1, 1, "2024", "", 1]
datasetList['nano/run3/2024/EGamma_2024C.list'] = [1, 1, "2024", "", 1]
datasetList['nano/run3/2024/EGamma_2024D.list'] = [1, 1, "2024", "", 1]
datasetList['nano/run3/2024/EGamma_2024E-v1.list'] = [1, 1, "2024", "", 1]
datasetList['nano/run3/2024/EGamma_2024E-v2.list'] = [1, 1, "2024", "", 1]
datasetList['nano/run3/2024/EGamma_2024F.list'] = [1, 1, "2024", "", 1]
datasetList['nano/run3/2024/EGamma_2024G.list'] = [1, 1, "2024", "", 1]
datasetList['nano/run3/2024/EGamma_2024H.list'] = [1, 1, "2024", "", 1]



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





