#define Events_cxx
#include <iostream>
#include <vector>
//LOCAL
#include "Events.hh"
//ROOT
#include <TFile.h>
#include <TH1F.h>
#include <TLorentzVector.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>


void Events::DeleteHeap()
{
  DeleteHeapFatJets();
  DeleteHeapSubJets();
  DeleteHeapCorrT1METJet();
  DeleteHeapElectrons();
  DeleteHeapFsrPhotons();
  DeleteHeapGenVariables();
  DeleteHeapIsoTracks();
  DeleteHeapPhotons();
  DeleteHeapJets();
  DeleteHeapMuons();
  DeleteHeapSoftActivityJet();
  DeleteHeapTaus();
  DeleteHeapTriggerObject();
  DeleteHeapSVandOtherPV();
};

//clean fatjets form d-memory
void Events::DeleteHeapFatJets()
{
 // delete[] FatJet_particleNet_massCorr;
 // delete[] FatJet_particleNetWithMass_HbbvsQCD;
 // delete[] FatJet_particleNetWithMass_HccvsQCD;
 // delete[] FatJet_particleNetWithMass_QCD;
  delete[] FatJet_particleNet_QCD;
  delete[] FatJet_particleNet_QCD0HF;
  delete[] FatJet_particleNet_QCD1HF;
  delete[] FatJet_particleNet_QCD2HF;
  delete[] FatJet_particleNet_XbbVsQCD;
  delete[] FatJet_particleNet_XccVsQCD;
  delete[] FatJet_particleNet_XggVsQCD;
  delete[] FatJet_particleNet_XqqVsQCD; 
  
  delete[] FatJet_area;   //[nFatJet]
  delete[] FatJet_btagDDBvLV2;   //[nFatJet]
  delete[] FatJet_btagDDCvBV2;   //[nFatJet]
  delete[] FatJet_btagDDCvLV2;   //[nFatJet]
  delete[] FatJet_btagHbb;   //[nFatJet]
  
 // delete[] FatJet_deepTagMD_WvsQCD;
 // delete[] FatJet_deepTagMD_ZvsQCD;
 // delete[] FatJet_deepTag_WvsQCD;
 // delete[] FatJet_deepTag_ZvsQCD;
  
 //delete[] FatJet_dRLep;   //[nFatJet]
  //delete[] FatJet_bTagHbb;   //[nFatJet]
  //delete[] FatJet_bTagHcc;   //[nFatJet]
  //delete[] FatJet_deepTagHqqqq;   //[nFatJet]
  //delete[] FatJet_deepTagMDHbb;   //[nFatJet]
  //delete[] FatJet_deepTagMDHcc;   //[nFatJet]
  //delete[] FatJet_deepTagMDHqqqq;   //[nFatJet]
  //delete[] FatJet_deepTagMDQCDbb;   //[nFatJet]
 // delete[] FatJet_deepTagMDQCDcc;   //[nFatJet]
  //delete[] FatJet_deepTagMDWcq;   //[nFatJet]
  //delete[] FatJet_deepTagMDWqq;   //[nFatJet]
  //delete[] FatJet_deepTagMDZbb;   //[nFatJet]
  //delete[] FatJet_deepTagMDZcc;   //[nFatJet]
  //delete[] FatJet_deepTagMDZqq;   //[nFatJet]
  //delete[] FatJet_deepTagQCDbb;   //[nFatJet]
  //delete[] FatJet_deepTagQCDcc;   //[nFatJet]
  //delete[] FatJet_deepTagWcq;   //[nFatJet]
  //delete[] FatJet_deepTagWqq;   //[nFatJet]
  //delete[] FatJet_deepTagZbb;   //[nFatJet]
  //delete[] FatJet_deepTagZcc;   //[nFatJet]
  //delete[] FatJet_deepTagZqq;   //[nFatJet]

  delete[] FatJet_eta;   //[nFatJet]
  delete[] FatJet_lsf3;   //[nFatJet]
  delete[] FatJet_mass;   //[nFatJet]
  delete[] FatJet_msoftdrop;   //[nFatJet]
  delete[] FatJet_n2b1;   //[nFatJet]
  delete[] FatJet_n3b1;   //[nFatJet]
  delete[] FatJet_phi;   //[nFatJet]
  delete[] FatJet_pt;   //[nFatJet]
  delete[] FatJet_rawFactor;   //[nFatJet]
  
//delete[] FatJet_rawmsoftdrop;   //[nFatJet]


  delete[] FatJet_tau1;   //[nFatJet]
  delete[] FatJet_tau2;   //[nFatJet]
  delete[] FatJet_tau3;   //[nFatJet]
  delete[] FatJet_tau4;   //[nFatJet]
  delete[] FatJet_electronIdx3SJ;   //[nFatJet]
  
//delete[]   FatJet_idLep;   //[nFatJet]

  delete[] FatJet_jetId;   //[nFatJet]
  delete[] FatJet_muonIdx3SJ;   //[nFatJet]
  delete[] FatJet_nBHadrons;   //[nFatJet]
  delete[] FatJet_nCHadrons;   //[nFatJet]
  delete[] FatJet_nConstituents;   //[nFatJet]
  delete[] FatJet_subJetIdx1;   //[nFatJet]
  delete[] FatJet_subJetIdx2;   //[nFatJet]
};

void Events::DeleteHeapCorrT1METJet()
{
  delete[] CorrT1METJet_area; //[nCorrT1METJet]
  delete[] CorrT1METJet_eta; //[nCorrT1METJet]
  delete[] CorrT1METJet_muonSubtrFactor; //[nCorrT1METJet]
  delete[] CorrT1METJet_phi; //[nCorrT1METJet]
  delete[] CorrT1METJet_rawPt; //[nCorrT1METJet]
};

void Events::DeleteHeapSubJets()
{
  delete[] SubJet_area;  //[nSubJet]
  delete[] SubJet_btagCSVV2;  //[nSubJet]
  delete[] SubJet_btagDeepB;  //[nSubJet]
  delete[] SubJet_eta;  //[nSubJet]
  delete[] SubJet_mass;  //[nSubJet]
  delete[] SubJet_phi;  //[nSubJet]
  delete[] SubJet_pt;  //[nSubJet]
  delete[] SubJet_rawFactor;  //[nSubJet]
  delete[] SubJet_nBHadrons;   //[nSubJet]
  delete[] SubJet_nCHadrons;   //[nSubJet]
};

void Events::DeleteHeapElectrons()
{
  delete[] Electron_deltaEtaSC;;   //[nElectron]
  delete[] Electron_dr03EcalRecHitSumEt;;   //[nElectron]
  delete[] Electron_dr03HcalDepth1TowerSumEt;;   //[nElectron]
  delete[] Electron_dr03TkSumPt;;   //[nElectron]
  delete[] Electron_dr03TkSumPtHEEP;;   //[nElectron]
  delete[] Electron_dxy;;   //[nElectron]
  delete[] Electron_dxyErr;;   //[nElectron]
  delete[] Electron_dz;;   //[nElectron]
  delete[] Electron_dzErr;;   //[nElectron]
 // delete[] Electron_eCorr;;   //[nElectron]
  delete[] Electron_eInvMinusPInv;;   //[nElectron]
  delete[] Electron_energyErr;;   //[nElectron]
  delete[] Electron_eta;;   //[nElectron]
  delete[] Electron_hoe;;   //[nElectron]
  delete[] Electron_ip3d;;   //[nElectron]
  delete[] Electron_jetPtRelv2;;   //[nElectron]
  delete[] Electron_jetRelIso;;   //[nElectron]
  delete[] Electron_mass;;   //[nElectron]
  delete[] Electron_miniPFRelIso_all;;   //[nElectron]
  delete[] Electron_miniPFRelIso_chg;;   //[nElectron]
 // delete[] Electron_mvaFall17V1Iso;;   //[nElectron]
 // delete[] Electron_mvaFall17V1noIso;;   //[nElectron]
  //delete[] Electron_mvaFall17V2Iso;;   //[nElectron]
  //delete[] Electron_mvaFall17V2noIso;;   //[nElectron]
  delete[] Electron_pfRelIso03_all;;   //[nElectron]
  delete[] Electron_pfRelIso03_chg;;   //[nElectron]
  delete[] Electron_phi;;   //[nElectron]
  delete[] Electron_pt;;   //[nElectron]
  delete[] Electron_r9;;   //[nElectron]
  delete[] Electron_sieie;;   //[nElectron]
  delete[] Electron_sip3d;;   //[nElectron]
  delete[] Electron_mvaTTH;;   //[nElectron]
  delete[]   Electron_charge;;   //[nElectron]
  delete[]   Electron_cutBased;;   //[nElectron]
  //delete[]   Electron_cutBased_Fall17_V1;;   //[nElectron]
  delete[]   Electron_jetIdx;;   //[nElectron]
  delete[]   Electron_pdgId;;   //[nElectron]
  delete[]   Electron_photonIdx;;   //[nElectron]
  delete[]   Electron_tightCharge;;   //[nElectron]
  delete[]   Electron_vidNestedWPBitmap;;   //[nElectron]
  delete[]   Electron_vidNestedWPBitmapHEEP;;   //[nElectron]
  delete[]  Electron_convVeto;;   //[nElectron]
  delete[]  Electron_cutBased_HEEP;;   //[nElectron]
  delete[]  Electron_isPFcand;;   //[nElectron]
  delete[]  Electron_lostHits;;   //[nElectron]
  //delete[]  Electron_mvaFall17V1Iso_WP80;;   //[nElectron]
  //delete[]  Electron_mvaFall17V1Iso_WP90;;   //[nElectron]
  //delete[]  Electron_mvaFall17V1Iso_WPL;;   //[nElectron]
  //delete[]  Electron_mvaFall17V1noIso_WP80;;   //[nElectron]
  //delete[]  Electron_mvaFall17V1noIso_WP90;;   //[nElectron]
  //delete[]  Electron_mvaFall17V1noIso_WPL;;   //[nElectron]
  //delete[]  Electron_mvaFall17V2Iso_WP80;;   //[nElectron]
  //delete[]  Electron_mvaFall17V2Iso_WP90;;   //[nElectron]
  //delete[]  Electron_mvaFall17V2Iso_WPL;;   //[nElectron]
  //delete[]  Electron_mvaFall17V2noIso_WP80;;   //[nElectron]
  //delete[]  Electron_mvaFall17V2noIso_WP90;;   //[nElectron]
  //delete[]  Electron_mvaFall17V2noIso_WPL;;   //[nElectron]
  delete[]  Electron_seedGain;;   //[nElectron]
  delete[]   Electron_genPartIdx;;   //[nElectron]
  delete[]  Electron_genPartFlav;;   //[nElectron]
  //delete[]  Electron_cleanmask;;   //[nElectron]
};

void Events::DeleteHeapFsrPhotons()
{
  delete[] FsrPhoton_dROverEt2;   //[nFsrPhoton]
  delete[] FsrPhoton_eta;   //[nFsrPhoton]
  delete[] FsrPhoton_phi;   //[nFsrPhoton]
  delete[] FsrPhoton_pt;   //[nFsrPhoton]
  delete[] FsrPhoton_relIso03;   //[nFsrPhoton]
  delete[] FsrPhoton_muonIdx;   //[nFsrPhoton]
};

void Events::DeleteHeapGenVariables()
{
  delete[] GenJetAK8_eta;   //[nGenJetAK8]
  delete[] GenJetAK8_mass;   //[nGenJetAK8]
  delete[] GenJetAK8_phi;   //[nGenJetAK8]
  delete[] GenJetAK8_pt;   //[nGenJetAK8]
  delete[] GenJetAK8_partonFlavour;   //[nGenJetAK8]
  delete[] GenJetAK8_hadronFlavour;   //[nGenJetAK8]

  delete[] GenJet_eta;   //[nGenJet]
  delete[] GenJet_mass;   //[nGenJet]
  delete[] GenJet_phi;   //[nGenJet]
  delete[] GenJet_pt;   //[nGenJet]
  delete[] GenJet_partonFlavour;   //[nGenJet]
  delete[] GenJet_hadronFlavour;   //[nGenJet]

  delete[] GenPart_eta;   //[nGenPart]
  delete[] GenPart_mass;   //[nGenPart]
  delete[] GenPart_phi;   //[nGenPart]
  delete[] GenPart_pt;   //[nGenPart]
  delete[]  GenPart_genPartIdxMother;   //[nGenPart]
  delete[]  GenPart_pdgId;   //[nGenPart]
  delete[]  GenPart_status;   //[nGenPart]
  delete[]  GenPart_statusFlags;   //[nGenPart]

  delete[] SubGenJetAK8_eta;   //[nSubGenJetAK8]
  delete[] SubGenJetAK8_mass;   //[nSubGenJetAK8]
  delete[] SubGenJetAK8_phi;   //[nSubGenJetAK8]
  delete[] SubGenJetAK8_pt;   //[nSubGenJetAK8]

  delete[] GenVisTau_eta;   //[nGenVisTau]
  delete[] GenVisTau_mass;   //[nGenVisTau]
  delete[] GenVisTau_phi;   //[nGenVisTau]
  delete[] GenVisTau_pt;   //[nGenVisTau]
  delete[]   GenVisTau_charge;   //[nGenVisTau]
  delete[]   GenVisTau_genPartIdxMother;   //[nGenVisTau]
  delete[]   GenVisTau_status;   //[nGenVisTau]

  delete[] LHEPdfWeight;   //[nLHEPdfWeight]

  delete[] LHEReweightingWeight;   //[nLHEReweightingWeight]

  delete[] PSWeight;   //[nPSWeight]

  delete[] LHEPart_pt;   //[nLHEPart]
  delete[] LHEPart_eta;   //[nLHEPart]
  delete[] LHEPart_phi;   //[nLHEPart]
  delete[] LHEPart_mass;   //[nLHEPart]
  delete[]   LHEPart_pdgId;   //[nLHEPart]

  delete[] GenDressedLepton_eta;   //[nGenDressedLepton]
  delete[] GenDressedLepton_mass;   //[nGenDressedLepton]
  delete[] GenDressedLepton_phi;   //[nGenDressedLepton]
  delete[] GenDressedLepton_pt;   //[nGenDressedLepton]
  delete[] GenDressedLepton_pdgId;   //[nGenDressedLepton]
  delete[] GenDressedLepton_hasTauAnc;   //[nGenDressedLepton]
};

void Events::DeleteHeapIsoTracks()
{
  delete[] IsoTrack_dxy;   //[nIsoTrack]
  delete[] IsoTrack_dz;   //[nIsoTrack]
  delete[] IsoTrack_eta;   //[nIsoTrack]
  delete[] IsoTrack_pfRelIso03_all;   //[nIsoTrack]
  delete[] IsoTrack_pfRelIso03_chg;   //[nIsoTrack]
  delete[] IsoTrack_phi;   //[nIsoTrack]
  delete[] IsoTrack_pt;   //[nIsoTrack]
  delete[] IsoTrack_miniPFRelIso_all;   //[nIsoTrack]
  delete[] IsoTrack_miniPFRelIso_chg;   //[nIsoTrack]
  delete[]   IsoTrack_fromPV;   //[nIsoTrack]
  delete[]   IsoTrack_pdgId;   //[nIsoTrack]
  delete[]  IsoTrack_isHighPurityTrack;   //[nIsoTrack]
  delete[]  IsoTrack_isPFcand;   //[nIsoTrack]
  delete[]  IsoTrack_isFromLostTrack;   //[nIsoTrack]
};

void Events::DeleteHeapPhotons()
{
 // delete[] Photon_eCorr;   //[nPhoton]
  delete[] Photon_energyErr;   //[nPhoton]
  delete[] Photon_eta;   //[nPhoton]
  delete[] Photon_hoe;   //[nPhoton]
  //delete[] Photon_mass;   //[nPhoton]
  delete[] Photon_mvaID;   //[nPhoton]
  //delete[] Photon_mvaIDV1;   //[nPhoton]
  //delete[] Photon_pfRelIso03_all;   //[nPhoton]
  //delete[] Photon_pfRelIso03_chg;   //[nPhoton]
  delete[] Photon_phi;   //[nPhoton]
  delete[] Photon_pt;   //[nPhoton]
  delete[] Photon_r9;   //[nPhoton]
  delete[] Photon_sieie;   //[nPhoton]
  //delete[]   Photon_charge;   //[nPhoton]
  //delete[]   Photon_cutBasedBitmap;   //[nPhoton]
  //delete[]   Photon_cutBasedV1Bitmap;   //[nPhoton]
  delete[]   Photon_electronIdx;   //[nPhoton]
  delete[]   Photon_jetIdx;   //[nPhoton]
 // delete[]   Photon_pdgId;   //[nPhoton]
  delete[]   Photon_vidNestedWPBitmap;   //[nPhoton]
  delete[]   Photon_genPartIdx;   //[nPhoton]
  delete[]  Photon_electronVeto;   //[nPhoton]
  delete[]  Photon_isScEtaEB;   //[nPhoton]
  delete[]  Photon_isScEtaEE;   //[nPhoton]
  delete[]  Photon_mvaID_WP80;   //[nPhoton]
  delete[]  Photon_mvaID_WP90;   //[nPhoton]
  delete[]  Photon_pixelSeed;   //[nPhoton]
  delete[] Photon_seedGain;   //[nPhoton]
  delete[] Photon_genPartFlav;   //[nPhoton]
  //delete[] Photon_cleanmask;   //[nPhoton]
};

void Events::DeleteHeapJets()
{
  delete[] Jet_area;   //[nJet]
 // delete[] Jet_btagCMVA;   //[nJet]
  //delete[] Jet_btagCSVV2;   //[nJet]
 // delete[] Jet_btagDeepB;   //[nJet]
 // delete[] Jet_btagDeepC;   //[nJet]
  delete[] Jet_btagDeepFlavB;   //[nJet]
  //delete[] Jet_btagDeepFlavC;   //[nJet]
  delete[] Jet_chEmEF;   //[nJet]
  delete[] Jet_chHEF;   //[nJet]
  delete[] Jet_eta;   //[nJet]
  //delete[] Jet_jercCHF;   //[nJet]
  //delete[] Jet_jercCHPUF;   //[nJet]
  delete[] Jet_mass;   //[nJet]
  delete[] Jet_muEF;   //[nJet]
  delete[] Jet_muonSubtrFactor;   //[nJet]
  delete[] Jet_neEmEF;   //[nJet]
  delete[] Jet_neHEF;   //[nJet]
  delete[] Jet_phi;   //[nJet]
  delete[] Jet_pt;   //[nJet]
  //delete[] Jet_qgl;   //[nJet]
  delete[] Jet_rawFactor;   //[nJet]
  //delete[] Jet_bRegCorr;   //[nJet]
  //delete[] Jet_bRegRes;   //[nJet]
  delete[] Jet_electronIdx1;   //[nJet]
  delete[] Jet_electronIdx2;   //[nJet]
  delete[] Jet_jetId;   //[nJet]
  delete[] Jet_muonIdx1;   //[nJet]
  delete[] Jet_muonIdx2;   //[nJet]
  delete[] Jet_nConstituents;   //[nJet]
  delete[] Jet_nElectrons;   //[nJet]
  delete[] Jet_nMuons;   //[nJet]
  //delete[] Jet_puId;   //[nJet]
  delete[] Jet_genJetIdx;   //[nJet]
  delete[] Jet_hadronFlavour;   //[nJet]
  delete[] Jet_partonFlavour;   //[nJet]
  //delete[] Jet_cleanmask;   //[nJet]
};

void Events::DeleteHeapMuons()
{
  delete[] Muon_dxy;   //[nMuon]
  delete[] Muon_dxyErr;   //[nMuon]
  delete[] Muon_dz;   //[nMuon]
  delete[] Muon_dzErr;   //[nMuon]
  delete[] Muon_eta;   //[nMuon]
  delete[] Muon_ip3d;   //[nMuon]
  delete[] Muon_jetPtRelv2;   //[nMuon]
  delete[] Muon_jetRelIso;   //[nMuon]
  delete[] Muon_mass;   //[nMuon]
  delete[] Muon_miniPFRelIso_all;   //[nMuon]
  delete[] Muon_miniPFRelIso_chg;   //[nMuon]
  delete[] Muon_pfRelIso03_all;   //[nMuon]
  delete[] Muon_pfRelIso03_chg;   //[nMuon]
  delete[] Muon_pfRelIso04_all;   //[nMuon]
  delete[] Muon_phi;   //[nMuon]
  delete[] Muon_pt;   //[nMuon]
  delete[] Muon_ptErr;   //[nMuon]
  delete[] Muon_segmentComp;   //[nMuon]
  delete[] Muon_sip3d;   //[nMuon]
  delete[] Muon_softMva;   //[nMuon]
  delete[] Muon_tkRelIso;   //[nMuon]
  delete[] Muon_tunepRelPt;   //[nMuon]
  delete[] Muon_mvaLowPt;   //[nMuon]
  delete[] Muon_mvaTTH;   //[nMuon]
  delete[]   Muon_charge;   //[nMuon]
  delete[]   Muon_jetIdx;   //[nMuon]
  delete[]   Muon_nStations;   //[nMuon]
  delete[]   Muon_nTrackerLayers;   //[nMuon]
  delete[]   Muon_pdgId;   //[nMuon]
  delete[]   Muon_tightCharge;   //[nMuon]
  delete[]   Muon_fsrPhotonIdx;   //[nMuon]
  delete[]   Muon_genPartIdx;   //[nMuon]
  delete[]  Muon_inTimeMuon;   //[nMuon]
  delete[]  Muon_isGlobal;   //[nMuon]
  delete[]  Muon_isPFcand;   //[nMuon]
  delete[]  Muon_isTracker;   //[nMuon]
  delete[]  Muon_looseId;   //[nMuon]
  delete[]  Muon_mediumId;   //[nMuon]
  delete[]  Muon_mediumPromptId;   //[nMuon]
  delete[]  Muon_softId;   //[nMuon]
  delete[]  Muon_softMvaId;   //[nMuon]
  delete[]  Muon_tightId;   //[nMuon]
  delete[]  Muon_triggerIdLoose;   //[nMuon]
  delete[] Muon_miniIsoId;   //[nMuon]
  delete[] Muon_multiIsoId;   //[nMuon]
  delete[] Muon_mvaMuID;   //[nMuon]
  delete[] Muon_pfIsoId;   //[nMuon]
  delete[] Muon_highPtId;   //[nMuon]
 // delete[] Muon_cleanmask;   //[nMuon]
  delete[] Muon_tkIsoId;   //[nMuon]
  delete[] Muon_genPartFlav;   //[nMuon]
};

void Events::DeleteHeapSoftActivityJet()
{
  delete[] SoftActivityJet_eta;   //[nSoftActivityJet]
  delete[] SoftActivityJet_phi;   //[nSoftActivityJet]
  delete[] SoftActivityJet_pt;   //[nSoftActivityJet]
};

void Events::DeleteHeapTaus()
{
  delete[] Tau_chargedIso;   //[nTau]
  delete[] Tau_dxy;   //[nTau]
  delete[] Tau_dz;   //[nTau]
  delete[] Tau_eta;   //[nTau]
  delete[] Tau_leadTkDeltaEta;   //[nTau]
  delete[] Tau_leadTkDeltaPhi;   //[nTau]
  delete[] Tau_leadTkPtOverTauPt;   //[nTau]
  delete[] Tau_mass;   //[nTau]
  delete[] Tau_neutralIso;   //[nTau]
  delete[] Tau_phi;   //[nTau]
  delete[] Tau_photonsOutsideSignalCone;   //[nTau]
  delete[] Tau_pt;   //[nTau]
  delete[] Tau_puCorr;   //[nTau]
 // delete[] Tau_rawAntiEle;   //[nTau]
  //delete[] Tau_rawAntiEle2018;   //[nTau]
  delete[] Tau_rawDeepTau2017v2p1VSe;   //[nTau]
  delete[] Tau_rawDeepTau2017v2p1VSjet;   //[nTau]
  delete[] Tau_rawDeepTau2017v2p1VSmu;   //[nTau]
  delete[] Tau_rawIso;   //[nTau]
  delete[] Tau_rawIsodR03;   //[nTau]
  //delete[] Tau_rawMVAnewDM2017v2;   //[nTau]
  //delete[] Tau_rawMVAoldDM;   //[nTau]
  //delete[] Tau_rawMVAoldDM2017v1;   //[nTau]
  //delete[] Tau_rawMVAoldDM2017v2;   //[nTau]
  //delete[] Tau_rawMVAoldDMdR032017v2;   //[nTau]
  delete[]   Tau_charge;   //[nTau]
  delete[]   Tau_decayMode;   //[nTau]
  delete[]   Tau_jetIdx;   //[nTau]
  //delete[]   Tau_rawAntiEleCat;   //[nTau]
  //delete[]   Tau_rawAntiEleCat2018;   //[nTau]
  delete[]   Tau_genPartIdx;   //[nTau]
  //delete[]  Tau_idDecayMode;   //[nTau]
  delete[]  Tau_idDecayModeNewDMs;   //[nTau]
  //delete[] Tau_idAntiEle;   //[nTau]
  //delete[] Tau_idAntiEle2018;   //[nTau]
  delete[] Tau_idAntiMu;   //[nTau]
  delete[] Tau_idDeepTau2017v2p1VSe;   //[nTau]
  delete[] Tau_idDeepTau2017v2p1VSjet;   //[nTau]
  delete[] Tau_idDeepTau2017v2p1VSmu;   //[nTau]
  //delete[] Tau_idMVAnewDM2017v2;   //[nTau]
  //delete[] Tau_idMVAoldDM;   //[nTau]
  //delete[] Tau_idMVAoldDM2017v1;   //[nTau]
 // delete[] Tau_idMVAoldDM2017v2;   //[nTau]
 // delete[] Tau_idMVAoldDMdR032017v2;   //[nTau]
  delete[] Tau_genPartFlav;   //[nTau]
  //delete[] Tau_cleanmask;   //[nTau]
};

void Events::DeleteHeapTriggerObject()
{
  delete[] TrigObj_pt;   //[nTrigObj]
  delete[] TrigObj_eta;   //[nTrigObj]
  delete[] TrigObj_phi;   //[nTrigObj]
  delete[] TrigObj_l1pt;   //[nTrigObj]
  delete[] TrigObj_l1pt_2;   //[nTrigObj]
  delete[] TrigObj_l2pt;   //[nTrigObj]
  delete[] TrigObj_id;   //[nTrigObj]
  delete[] TrigObj_l1iso;   //[nTrigObj]
  delete[] TrigObj_l1charge;   //[nTrigObj]
  delete[] TrigObj_filterBits;   //[nTrigObj]
};

void Events::DeleteHeapSVandOtherPV()
{
  delete[] OtherPV_z;   //[nOtherPV]

  delete[] SV_dlen;   //[nSV]
  delete[] SV_dlenSig;   //[nSV]
  delete[] SV_dxy;   //[nSV]
  delete[] SV_dxySig;   //[nSV]
  delete[] SV_pAngle;   //[nSV]
  delete[] SV_chi2;   //[nSV]
  delete[] SV_eta;   //[nSV]
  delete[] SV_mass;   //[nSV]
  delete[] SV_ndof;   //[nSV]
  delete[] SV_phi;   //[nSV]
  delete[] SV_pt;   //[nSV]
  delete[] SV_x;   //[nSV]
  delete[] SV_y;   //[nSV]
  delete[] SV_z;   //[nSV]
};

void Events::CreateOutputTree()
{
  tree_out = new TTree("hh", "hh");

  //define event output branches;
  tree_out->Branch("run",               &run,     "run/i");      //run number
  tree_out->Branch("luminosityBlock",   &luminosityBlock,     "luminosityBlock/i");      //lumi
  tree_out->Branch("event",             &event,     "event/l");      //event number
  //gen event-genWeight
  tree_out->Branch("genWeight",         &genWeight,     "genWeight/F");      //event number
  //MET
  tree_out->Branch("ChsMET_phi",        &ChsMET_phi,     "ChsMET_phi/F");      //
  tree_out->Branch("ChsMET_pt",         &ChsMET_pt,     "ChsMET_pt/F");      //
  tree_out->Branch("ChsMET_sumEt",      &ChsMET_sumEt,     "ChsMET_sumEt/F");      //
  //tree_out->Branch("",      &,     "[nFatJet]/I");      //
  //tree_out->Branch("",      &,     "/");      //

  //---------------------test Gen
  //tree_out->Branch("nGenJet", &nGenJet, "nGenJet/i" );

  //define hh-gen-level
  tree_out->Branch("h_gen_pt",       h_gen_pt,      "h_gen_pt[2]/F");      //
  tree_out->Branch("h_gen_eta",      h_gen_eta,     "h_gen_eta[2]/F");      //
  tree_out->Branch("h_gen_phi",      h_gen_phi,     "h_gen_phi[2]/F");      //

  //index to fatjet in hh candidate
  tree_out->Branch("hh_fatjet_idx",   hh_fatjet_idx, "hh_fatjet_idx[2]/i");

  //define fat-jet variables
  tree_out->Branch("nFatJet",      &nFatJet,     "nFatJet/i");      //
  tree_out->Branch("FatJet_msoftdrop",      FatJet_msoftdrop,     "FatJet_msoftdrop[nFatJet]/F");      //
  tree_out->Branch("FatJet_n2b1",      FatJet_n2b1,     "FatJet_n2b1[nFatJet]/F");      //
  tree_out->Branch("FatJet_n3b1",      FatJet_n3b1,     "FatJet_n3b1[nFatJet]/F");      //
  tree_out->Branch("FatJet_pt",      FatJet_pt,     "FatJet_pt[nFatJet]/F");      //
  //tree_out->Branch("FatJet_LSrawmsoftdrop",      FatJet_LSrawmsoftdrop,     "FatJet_LSrawmsoftdrop[nFatJet]/F");      //
  //tree_out->Branch("FatJet_LSsubJet1btagDeepB",      FatJet_LSsubJet1btagDeepB,     "FatJet_LSsubJet1btagDeepB[nFatJet]/F");      //
  //tree_out->Branch("FatJet_LSsubJet2btagDeepB",      FatJet_LSsubJet2btagDeepB,     "FatJet_LSsubJet2btagDeepB[nFatJet]/F");      //
  tree_out->Branch("FatJet_tau1",      FatJet_tau1,     "FatJet_tau1[nFatJet]/F");      //
  tree_out->Branch("FatJet_tau2",      FatJet_tau2,     "FatJet_tau2[nFatJet]/F");      //
  tree_out->Branch("FatJet_tau3",      FatJet_tau3,     "FatJet_tau3[nFatJet]/F");      //
  tree_out->Branch("FatJet_tau4",      FatJet_tau4,     "FatJet_tau4[nFatJet]/F");      //
  tree_out->Branch("FatJet_area",      FatJet_area,     "FatJet_area[nFatJet]/F");      //
 // tree_out->Branch("FatJet_btagDDBvLV2",      FatJet_btagDDBvLV2,     "FatJet_btagDDBvLV2[nFatJet]/F");      //
 // tree_out->Branch("FatJet_btagDDCvBV2",      FatJet_btagDDCvBV2,     "FatJet_btagDDCvBV2[nFatJet]/F");      //
 // tree_out->Branch("FatJet_btagDDCvLV2",      FatJet_btagDDCvLV2,     "FatJet_btagDDCvLV2[nFatJet]/F");      //
 // tree_out->Branch("FatJet_btagHbb",      FatJet_btagHbb,     "FatJet_btagHbb[nFatJet]/F");      //
  //tree_out->Branch("FatJet_dRLep",      FatJet_btagHbb,     "FatJet_btagHbb[nFatJet]/F");      //
  //tree_out->Branch("FatJet_deepTagHbb",      FatJet_deepTagHbb,     "FatJet_deepTagHbb[nFatJet]/F");      //
  //tree_out->Branch("FatJet_deepTagHcc",      FatJet_deepTagHcc,     "FatJet_deepTagHcc[nFatJet]/F");      //
  //tree_out->Branch("FatJet_deepTagHqqqq",      FatJet_deepTagHqqqq,     "FatJet_deepTagHqqqq[nFatJet]/F");      //
  //tree_out->Branch("FatJet_deepTagMDHbb",      FatJet_deepTagMDHbb,     "FatJet_deepTagMDHbb[nFatJet]/F");      //
  //tree_out->Branch("FatJet_deepTagMDHcc",      FatJet_deepTagMDHcc,     "FatJet_deepTagMDHcc[nFatJet]/F");      //
  //tree_out->Branch("FatJet_deepTagMDHqqqq",      FatJet_deepTagMDHqqqq,     "FatJet_deepTagMDHqqqq[nFatJet]/F");      //
  //tree_out->Branch("FatJet_deepTagMDQCDbb",      FatJet_deepTagMDQCDbb,     "FatJet_deepTagMDQCDbb[nFatJet]/F");      //
  //tree_out->Branch("FatJet_deepTagMDQCDcc",      FatJet_deepTagMDQCDcc,     "FatJet_deepTagMDQCDcc[nFatJet]/F");      //
  //tree_out->Branch("FatJet_deepTagMDWcq",      FatJet_deepTagMDWcq,     "FatJet_deepTagMDWcq[nFatJet]/F");      //
  //tree_out->Branch("FatJet_deepTagMDWqq",      FatJet_deepTagMDWqq,     "FatJet_deepTagMDWqq[nFatJet]/F");      //
  //tree_out->Branch("FatJet_deepTagMDZbb",      FatJet_deepTagMDZbb,     "FatJet_deepTagMDZbb[nFatJet]/F");      //
  //tree_out->Branch("FatJet_deepTagMDZcc",      FatJet_deepTagMDZcc,     "FatJet_deepTagMDZcc[nFatJet]/F");      //
  //tree_out->Branch("FatJet_deepTagMDZqq",      FatJet_deepTagMDZqq,     "FatJet_deepTagMDZqq[nFatJet]/F");      //
  //tree_out->Branch("FatJet_deepTagQCDbb",      FatJet_deepTagQCDbb,     "FatJet_deepTagQCDbb[nFatJet]/F");      //
  //tree_out->Branch("FatJet_deepTagQCDcc",      FatJet_deepTagQCDcc,     "FatJet_deepTagQCDcc[nFatJet]/F");      //
  //tree_out->Branch("FatJet_deepTagWcq",      FatJet_deepTagWcq,     "FatJet_deepTagWcq[nFatJet]/F");      //
  //tree_out->Branch("FatJet_deepTagWqq",      FatJet_deepTagWqq,     "FatJet_deepTagWqq[nFatJet]/F");      //
  //tree_out->Branch("FatJet_deepTagZbb",      FatJet_deepTagZbb,     "FatJet_deepTagZbb[nFatJet]/F");      //
  //tree_out->Branch("FatJet_deepTagZcc",      FatJet_deepTagZcc,     "FatJet_deepTagZcc[nFatJet]/F");      //
  //tree_out->Branch("FatJet_deepTagZqq",      FatJet_deepTagZqq,     "FatJet_deepTagZqq[nFatJet]/F");      //
  tree_out->Branch("FatJet_eta",      FatJet_eta,     "FatJet_eta[nFatJet]/F");      //
  tree_out->Branch("FatJet_lsf3",      FatJet_lsf3,     "FatJet_lsf3[nFatJet]/F");      //
  tree_out->Branch("FatJet_mass",      FatJet_mass,     "FatJet_mass[nFatJet]/F");      //
  //tree_out->Branch("FatJet_msoftdrop",      FatJet_msoftdrop,     "FatJet_msoftdrop[nFatJet]/F");      //
 // tree_out->Branch("FatJet_n2b1",      FatJet_n2b1,     "FatJet_n2b1[nFatJet]/F");      //
 // tree_out->Branch("FatJet_n3b1",      FatJet_n3b1,     "FatJet_n3b1[nFatJet]/F");      //
//  tree_out->Branch("FatJet_phi",      FatJet_phi,     "FatJet_phi[nFatJet]/F");      //
//  tree_out->Branch("FatJet_pt",      FatJet_pt,     "FatJet_pt[nFatJet]/F");      //
  tree_out->Branch("FatJet_rawFactor",      FatJet_rawFactor,     "FatJet_rawFactor[nFatJet]/F");      //
  //tree_out->Branch("FatJet_rawmsoftdrop",      FatJet_rawmsoftdrop,     "FatJet_rawmsoftdrop[nFatJet]/F");      //
  //tree_out->Branch("FatJet_tau1",      FatJet_tau1,     "FatJet_tau1[nFatJet]/F");      //
  //tree_out->Branch("FatJet_tau2",      FatJet_tau2,     "FatJet_tau2[nFatJet]/F");      //
  //tree_out->Branch("FatJet_tau3",      FatJet_tau3,     "FatJet_tau3[nFatJet]/F");      //
  //tree_out->Branch("FatJet_tau4",      FatJet_tau4,     "FatJet_tau4[nFatJet]/F");      //
  tree_out->Branch("FatJet_electronIdx3SJ",      FatJet_electronIdx3SJ,     "FatJet_electronIdx3SJ[nFatJet]/I");      //
  //tree_out->Branch("FatJet_idLep",      FatJet_idLep,     "FatJet_idLep[nFatJet]/I");      //
//  tree_out->Branch("FatJet_jetId",      FatJet_jetId,     "FatJet_jetId[nFatJet]/I");      //
  tree_out->Branch("FatJet_muonIdx3SJ",      FatJet_muonIdx3SJ,     "FatJet_muonIdx3SJ[nFatJet]/I");      //
  tree_out->Branch("FatJet_nBHadrons",      FatJet_nBHadrons,     "FatJet_nBHadrons[nFatJet]/I");      //
  tree_out->Branch("FatJet_nCHadrons",      FatJet_nCHadrons,     "FatJet_nCHadrons[nFatJet]/I");      //
  tree_out->Branch("FatJet_nConstituents",      FatJet_nConstituents,     "FatJet_nConstituents[nFatJet]/I");      //
  tree_out->Branch("FatJet_subJetIdx1",      FatJet_subJetIdx1,     "FatJet_subJetIdx1[nFatJet]/I");      //
  tree_out->Branch("FatJet_subJetIdx2",      FatJet_subJetIdx2,     "FatJet_subJetIdx2[nFatJet]/I");      //
  tree_out->Branch("FatJet_Hmatch",      FatJet_Hmatch,     "FatJet_Hmatch[nFatJet]/O");
  tree_out->Branch("FatJet_HgenIdx",      FatJet_HgenIdx,     "FatJet_HgenIdx[nFatJet]/I");
  tree_out->Branch("FatJet_HminDR",      FatJet_HminDR,     "FatJet_HminDR[nFatJet]/F");

  //triggers -- directly from Events TTree
  tree_out->Branch("HLT_PFHT1050",                                        &HLT_PFHT1050,                                       "HLT_PFHT1050/O");
  tree_out->Branch("HLT_AK8PFJet360_TrimMass30",                          &HLT_AK8PFJet360_TrimMass30,                         "HLT_AK8PFJet360_TrimMass30/O");
  tree_out->Branch("HLT_AK8PFJet380_TrimMass30",                          &HLT_AK8PFJet380_TrimMass30,                         "HLT_AK8PFJet380_TrimMass30/O");
  tree_out->Branch("HLT_AK8PFJet400_TrimMass30",                          &HLT_AK8PFJet400_TrimMass30,                         "HLT_AK8PFJet400_TrimMass30/O");
  tree_out->Branch("HLT_AK8PFJet420_TrimMass30",                          &HLT_AK8PFJet420_TrimMass30,                         "HLT_AK8PFJet420_TrimMass30/O");
  tree_out->Branch("HLT_AK8PFHT800_TrimMass50",                           &HLT_AK8PFHT800_TrimMass50,                          "HLT_AK8PFHT800_TrimMass50/O");
  tree_out->Branch("HLT_PFJet500",                                        &HLT_PFJet500,                                       "HLT_PFJet500/O");
  tree_out->Branch("HLT_AK8PFJet500",                                     &HLT_AK8PFJet500,                                    "HLT_AK8PFJet500/O");
  tree_out->Branch("HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p17",     &HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p17,    "HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p17/O");
  tree_out->Branch("HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p1",      &HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p1,     "HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p1/O");
  tree_out->Branch("HLT_AK8PFJet330_PFAK8BTagCSV_p17",                    &HLT_AK8PFJet330_PFAK8BTagCSV_p17,                   "HLT_AK8PFJet330_PFAK8BTagCSV_p17/O");
  tree_out->Branch("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",                                        &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ,                                       "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ/O");
  tree_out->Branch("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",                                        &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL,                                       "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL/O");
  tree_out->Branch("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL",                                        &HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL,                                       "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL/O");
  tree_out->Branch("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",                                        &HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ,                                       "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ/O");
  tree_out->Branch("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",                                        &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ,                                       "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ/O");
  tree_out->Branch("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL",                                        &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL,                                       "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL/O");
  //tree_out->Branch("",      &,     "/O");
  //tree_out->Branch("",      &,     "/O");
  //tree_out->Branch("",      &,     "/O");

  //hh candidate info
  tree_out->Branch("hh_pt",      &hh_pt,     "hh_pt/F");
  tree_out->Branch("hh_eta",     &hh_eta,    "hh_eta/F");
  tree_out->Branch("hh_phi",     &hh_phi,    "hh_phi/F");
  tree_out->Branch("hh_mass",    &hh_mass,   "hh_mass/F");
  //
  tree_out->Branch("hh_gen_pt",      &hh_gen_pt,     "hh_gen_pt/F");
  tree_out->Branch("hh_gen_eta",     &hh_gen_eta,    "hh_gen_eta/F");
  tree_out->Branch("hh_gen_phi",     &hh_gen_phi,    "hh_gen_phi/F");
  tree_out->Branch("hh_gen_mass",    &hh_gen_mass,   "hh_gen_mass/F");	
  //lepton candidate info
  tree_out->Branch("nElectron",      &nElectron,     "nElectron/i");      //
  tree_out->Branch("Electron_pt",      Electron_pt,     "Electron_pt[nElectron]/F");      //
  tree_out->Branch("Electron_eta",      Electron_eta,     "Electron_eta[nElectron]/F");      //
  tree_out->Branch("Electron_phi",      Electron_phi,     "Electron_phi[nElectron]/F");      //
  tree_out->Branch("Electron_charge",      Electron_charge,     "Electron_charge[nElectron]/I");      //
 // tree_out->Branch("Electron_mvaFall17V2Iso_WP80",      Electron_mvaFall17V2Iso_WP80,     "Electron_mvaFall17V2Iso_WP80[nElectron]/O");      //
 // tree_out->Branch("Electron_mvaFall17V2Iso_WP90",      Electron_mvaFall17V2Iso_WP90,     "Electron_mvaFall17V2Iso_WP90[nElectron]/O");      //
  //tree_out->Branch("Electron_mvaFall17V2Iso_WPL",      Electron_mvaFall17V2Iso_WPL,     "Electron_mvaFall17V2Iso_WPL[nElectron]/O");      //
  tree_out->Branch("nMuon",      &nMuon,     "nMuon/i");      //
  tree_out->Branch("Muon_pt",      Muon_pt,     "Muon_pt[nMuon]/F");      //
  tree_out->Branch("Muon_eta",      Muon_eta,     "Muon_eta[nMuon]/F");      //
  tree_out->Branch("Muon_phi",      Muon_phi,     "Muon_phi[nMuon]/F");      //
  tree_out->Branch("Muon_charge",      Muon_charge,     "Muon_charge[nMuon]/I");      //
  tree_out->Branch("Muon_looseId",      Muon_looseId,     "Muon_looseId[nMuon]/O");      //
  tree_out->Branch("Muon_mediumId",      Muon_mediumId,     "Muon_mediumId[nMuon]/O");      //
  tree_out->Branch("Muon_tightId",      Muon_tightId,     "Muon_tightId[nMuon]/O");      //
    
  //jets
  tree_out->Branch("nJet",      &nJet,     "nJet/i");      //
  tree_out->Branch("Jet_pt",      Jet_pt,     "Jet_pt[nJet]/F");      //
  tree_out->Branch("Jet_eta",      Jet_eta,     "Jet_eta[nJet]/F");      //
  tree_out->Branch("Jet_phi",      Jet_phi,     "Jet_phi[nJet]/F");      //
  //tree_out->Branch("Jet_btagCSVV2",      Jet_btagCSVV2,     "Jet_btagCSVV2[nJet]/F");      //

  

};

void Events::ResetOutputTreeVariables()
{
  hh_pt   = -999.;
  hh_eta  = -999.;
  hh_phi  = -999.;
  hh_mass = -999.;
  //
  hh_gen_pt   = -999.;
  hh_gen_eta  = -999.;
  hh_gen_phi  = -999.;
  hh_gen_mass = -999.;
  //
  for(int i = 0; i < 2; i++)
  {
    h_gen_pt[i]  = -999.;
    h_gen_eta[i] = -999.;
    h_gen_phi[i] = -999.;
    hh_fatjet_idx[i] = 999;
  }

  for(int i = 0; i < 100; i++)
  {
    FatJet_HminDR[i] = -999.;
    FatJet_Hmatch[i] = false;
    FatJet_HgenIdx[i] = -1;
  }

}
void Events::Loop()
{
  //   In a ROOT session, you can do:
  //      root> .L Events.C
  //      root> Events t
  //      root> t.GetEntry(12); // Fill t data members with entry number 12
  //      root> t.Show();       // Show values of entry 12
  //      root> t.Show(16);     // Read and show values of entry 16
  //      root> t.Loop();       // Loop on all entries
  //

  //     This is the loop skeleton where:
  //    jentry is the global entry number in the chain
  //    ientry is the entry number in the current Tree
  //  Note that the argument to GetEntry must be:
  //    jentry for TChain::GetEntry
  //    ientry for TTree::GetEntry and TBranch::GetEntry
  //
  //       To read only selected branches, Insert statements like:
  // METHOD1:
  //    fChain->SetBranchStatus("*",0);  // disable all branches
  //    fChain->SetBranchStatus("branchname",1);  // activate branchname
  // METHOD2: replace line
  //    fChain->GetEntry(jentry);       //read all branches
  //by  b_branchname->GetEntry(ientry); //read only this branch
  this->CreateOutputTree();
  TH1F* h_h1_pt = new TH1F("h1_pt", "h1_pt", 100, 0, 1000 );
  TH1F* h_h2_pt = new TH1F("h2_pt", "h2_pt", 100, 0, 1000 );
  TH2F* h_h1_h2_pt = new TH2F("h1_h2_pt","h1_h2_pt", 100, 0, 1000, 100, 0, 1000);

  //-------------------------------------
  //Normalization histograms
  //------------------------------------
  TH1F* NEvents = new TH1F("NEvents", "NEvents", 1, 0, 1);
  TH1F* NEvents_genweight = new TH1F("NEvents_genweight", "NEvents_genweight", 1, 0, 1);

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++)
  {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    NEvents->Fill( 0.5, 1.0 );
    NEvents_genweight->Fill(0.5, genWeight);


    // if (Cut(ientry) < 0) continue;
    //std::cout << ientry << std::endl;
    this->ResetOutputTreeVariables();
    //std::cout << "============H_GEN===============" << std::endl;
    //------------------------------
    //----find gen-higgs------------
    //------------------------------
    int current_mIndex = -1;
    std::vector< TLorentzVector > h_vector;

  for(int i = 0; i < nGenPart; i++)
    {
      if( abs(GenPart_pdgId[i]) == 5  && GenPart_pdgId[GenPart_genPartIdxMother[i]] == 25 && current_mIndex != GenPart_genPartIdxMother[i] )
      {
        //std::cout << GenPart_genPartIdxMother[i] << std::endl;
        // std::cout << "mother: " << GenPart_pdgId[GenPart_genPartIdxMother[i]]
        // << " PT: " << GenPart_pt[GenPart_genPartIdxMother[i]]
        // << " eta: " << GenPart_eta[GenPart_genPartIdxMother[i]]
        // << " phi: " << GenPart_phi[GenPart_genPartIdxMother[i]] << std::endl;
        TLorentzVector h;
        h.SetPtEtaPhiM( GenPart_pt[GenPart_genPartIdxMother[i]], GenPart_eta[GenPart_genPartIdxMother[i]], GenPart_phi[GenPart_genPartIdxMother[i]], GenPart_mass[GenPart_genPartIdxMother[i]] );
        h_vector.push_back(h);
        current_mIndex = GenPart_genPartIdxMother[i];
      }
    }

    if(h_vector.size() > 1)
    {
      h_h1_pt->Fill(h_vector.at(0).Pt());
      h_h2_pt->Fill(h_vector.at(1).Pt());
      h_h1_h2_pt->Fill(h_vector.at(0).Pt(), h_vector.at(1).Pt());
    }
    //------------------------------
    //-------find fatJet------------
    //------------------------------
    //std::cout << "======================" << std::endl;
    for(unsigned int i = 0; i < nFatJet; i++ )
    {
      TLorentzVector tmp_fatJet;
      tmp_fatJet.SetPtEtaPhiM(FatJet_pt[i],FatJet_eta[i],FatJet_phi[i],FatJet_msoftdrop[i]);
      float minDR = 999.;
      int match_idx = -1;
      for( int j = 0; j < h_vector.size(); j++)
      {
        if(tmp_fatJet.DeltaR(h_vector.at(j)) < minDR)
        {
          minDR = tmp_fatJet.DeltaR(h_vector.at(j));
          match_idx = j;
        }
      }

      FatJet_HminDR[i] = minDR;
      if(FatJet_HminDR[i] < 0.4)
      {
        FatJet_Hmatch[i] = true;
        FatJet_HgenIdx[i] = match_idx;
      }
      //FatJet_Hmatch[i] = true;
      // std::cout << "FatJet_deepTagHbb: " << FatJet_deepTagHbb[i] <<
      // "; pT: " << FatJet_pt[i] << " eta: " <<  FatJet_eta[i]
      // << " phi: " << FatJet_phi[i]
      // << "; FatJet_nBHadrons: " << FatJet_nBHadrons[i] <<
      // "; mass: " << FatJet_mass[i] << "; minDR:" << minDR << std::endl;
    }

    //-----------------------------
    //----------get hh cand -------
    //-----------------------------
    double sum_hh_pt = 0;
    TLorentzVector hh_candidate(0,0,0,0);
    unsigned int hh_candidate_idx1 = -1;
    unsigned int hh_candidate_idx2 = -1;
    for( int i  = 0; i < nFatJet; i++ )
    {
      for(int j = i+1; j < nFatJet; j++)
      {
        if ( FatJet_pt[i] + FatJet_pt[j] > sum_hh_pt)//pick hh candite with largest scalar PT sum
        {
          sum_hh_pt = FatJet_pt[i] + FatJet_pt[j];
          //h1
          TLorentzVector h1;
          h1.SetPtEtaPhiM(FatJet_pt[i],FatJet_eta[i],FatJet_phi[i],FatJet_msoftdrop[i]);
          //h2
          TLorentzVector h2;
          h2.SetPtEtaPhiM(FatJet_pt[j],FatJet_eta[j],FatJet_phi[j],FatJet_msoftdrop[j]);
          hh_candidate = h1+h2;
	  hh_candidate_idx1 = i;
	  hh_candidate_idx2 = j;
        }

      }
    }

    if(h_vector.size() >= 1)
      {
	//filling tree_out variables
	this->h_gen_pt[0] = h_vector.at(0).Pt();
	this->h_gen_eta[0] = h_vector.at(0).Eta();
	this->h_gen_phi[0] = h_vector.at(0).Phi();
	//
	if(h_vector.size() >= 2)
	  {
	    this->h_gen_pt[1] = h_vector.at(1).Pt();
	    this->h_gen_eta[1] = h_vector.at(1).Eta();
	    this->h_gen_phi[1] = h_vector.at(1).Phi();
	  }
      }
    //filling hh candidate variable
    this->hh_pt   = hh_candidate.Pt();
    this->hh_eta  = hh_candidate.Eta();
    this->hh_phi  = hh_candidate.Phi();
    this->hh_mass = hh_candidate.M();
    this->hh_fatjet_idx[0] = hh_candidate_idx1;
    this->hh_fatjet_idx[1] = hh_candidate_idx2;
    //gen level
    if(h_vector.size() > 1)
    {
      this->hh_gen_pt   = (h_vector.at(0)+h_vector.at(1)).Pt();
      this->hh_gen_eta  = (h_vector.at(0)+h_vector.at(1)).Eta();
      this->hh_gen_phi  = (h_vector.at(0)+h_vector.at(1)).Phi();
      this->hh_gen_mass = (h_vector.at(0)+h_vector.at(1)).M();
    }


    //------------------------------
    //apply skimming
    //-----------------------------
    //if( !(HLT_AK8PFJet360_TrimMass30 || HLT_AK8PFJet400_TrimMass30 || HLT_AK8PFJet420_TrimMass30) ) continue;

    if ( !(this->hh_pt > 0.0 && FatJet_pt[this->hh_fatjet_idx[0]] > this->fatjet_pt_trh && FatJet_pt[this->hh_fatjet_idx[1]] > this->fatjet_pt_trh) ) continue;

    this->tree_out->Fill();
  }

  TFile* fout = new TFile( fout_name.c_str(), "recreate");
  //h_h1_pt->Write();
  //h_h2_pt->Write();
  //h_h1_h2_pt->Write();
  NEvents->Write();
  NEvents_genweight->Write();
  tree_out->Write();
  delete fout;
}
