#include "genAna.h"

struct elecTriggers{
    bool eleTrigMatch_0_HLT_e120_lhloose; 
    bool eleTrigMatch_0_HLT_e140_lhloose_nod0; 
    bool eleTrigMatch_0_HLT_e24_lhmedium_L1EM20VH; 
    bool eleTrigMatch_0_HLT_e26_lhtight_nod0_ivarloose; 
    bool eleTrigMatch_0_HLT_e60_lhmedium; 
    bool eleTrigMatch_0_HLT_e60_lhmedium_nod0; 
};

struct muonTriggers{
    bool muTrigMatch_0_HLT_mu20_iloose_L1MU15; 
    bool muTrigMatch_0_HLT_mu26_ivarmedium; 
    bool muTrigMatch_0_HLT_mu50; 
};

struct elecSF{
    double elec_0_NOMINAL_EleEffSF_Isolation_TightLLH_d0z0_v13_FCLoose;
    double elec_0_NOMINAL_EleEffSF_Isolation_TightLLH_d0z0_v13_FCTight;
    double elec_0_NOMINAL_EleEffSF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_TightLLH_d0z0_v13_isolFCTight; 
    double elec_0_NOMINAL_EleEffSF_offline_RecoTrk;
    double elec_0_NOMINAL_EleEffSF_offline_TightLLH_d0z0_v13;
    double elec_0_NOMINAL_efficiency_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_TightLLH_d0z0_v13_isolFCTight;
};

struct muonSF{
    double muon_0_NOMINAL_MuEffSF_HLT_mu20_iloose_L1MU15_OR_HLT_mu50_QualMedium; 
    double muon_0_NOMINAL_MuEffSF_HLT_mu26_ivarmedium_OR_HLT_mu50_QualMedium; 
    double muon_0_NOMINAL_MuEffSF_IsoFCLoose; 
    double muon_0_NOMINAL_MuEffSF_IsoFCLoose_FixedRad; 
    double muon_0_NOMINAL_MuEffSF_IsoFCTight; 
    double muon_0_NOMINAL_MuEffSF_IsoFCTightTrackOnly; 
    double muon_0_NOMINAL_MuEffSF_IsoFCTightTrackOnly_FixedRad; 
    double muon_0_NOMINAL_MuEffSF_IsoFCTight_FixedRad; 
    double muon_0_NOMINAL_MuEffSF_IsoFixedCutHighPtTrackOnly; 
    double muon_0_NOMINAL_MuEffSF_Reco_QualMedium; 
    double muon_0_NOMINAL_MuEffSF_TTVA;
};


struct tauHadSF{
    double tau_0_NOMINAL_TauEffSF_JetRNNloose; 
    double tau_0_NOMINAL_TauEffSF_JetRNNmedium; 
    double tau_0_NOMINAL_TauEffSF_JetRNNtight; 
    double tau_0_NOMINAL_TauEffSF_LooseEleBDT_electron; 
    double tau_0_NOMINAL_TauEffSF_MediumEleBDT_electron; 
    double tau_0_NOMINAL_TauEffSF_reco;
};

struct jetSF{
    double jet_NOMINAL_central_jets_global_effSF_JVT; 
    double jet_NOMINAL_central_jets_global_ineffSF_JVT; 
    double jet_NOMINAL_forward_jets_global_effSF_JVT; 
    double jet_NOMINAL_forward_jets_global_ineffSF_JVT; 
    double jet_NOMINAL_global_effSF_MV2c10_FixedCutBEff_85; 
    double jet_NOMINAL_global_ineffSF_MV2c10_FixedCutBEff_85;
};

struct elecTriggers getElecTriggers(uint eleTrigMatch_0_HLT_e120_lhloose, uint eleTrigMatch_0_HLT_e140_lhloose_nod0, uint eleTrigMatch_0_HLT_e24_lhmedium_L1EM20VH, uint eleTrigMatch_0_HLT_e26_lhtight_nod0_ivarloose, uint eleTrigMatch_0_HLT_e60_lhmedium, uint eleTrigMatch_0_HLT_e60_lhmedium_nod0)
{
    elecTriggers ret;
    ret.eleTrigMatch_0_HLT_e120_lhloose                    = (bool)eleTrigMatch_0_HLT_e120_lhloose;
    ret.eleTrigMatch_0_HLT_e140_lhloose_nod0               = (bool)eleTrigMatch_0_HLT_e140_lhloose_nod0;
    ret.eleTrigMatch_0_HLT_e24_lhmedium_L1EM20VH           = (bool)eleTrigMatch_0_HLT_e24_lhmedium_L1EM20VH; 
    ret.eleTrigMatch_0_HLT_e26_lhtight_nod0_ivarloose      = (bool)eleTrigMatch_0_HLT_e26_lhtight_nod0_ivarloose; 
    ret.eleTrigMatch_0_HLT_e60_lhmedium                    = (bool)eleTrigMatch_0_HLT_e60_lhmedium;
    ret.eleTrigMatch_0_HLT_e60_lhmedium_nod0               = (bool)eleTrigMatch_0_HLT_e60_lhmedium_nod0; 

    return ret;
}

struct muonTriggers getMuonTriggers(uint muTrigMatch_0_HLT_mu20_iloose_L1MU15, uint muTrigMatch_0_HLT_mu26_ivarmedium, uint muTrigMatch_0_HLT_mu50){
    muonTriggers ret;
    ret.muTrigMatch_0_HLT_mu20_iloose_L1MU15              = (bool)muTrigMatch_0_HLT_mu20_iloose_L1MU15; 
    ret.muTrigMatch_0_HLT_mu26_ivarmedium                 = (bool)muTrigMatch_0_HLT_mu26_ivarmedium;
    ret.muTrigMatch_0_HLT_mu50                            = (bool)muTrigMatch_0_HLT_mu50;
    return ret;
}

struct elecSF getElecSF(float elec_0_NOMINAL_EleEffSF_Isolation_TightLLH_d0z0_v13_FCLoose,
                        float elec_0_NOMINAL_EleEffSF_Isolation_TightLLH_d0z0_v13_FCTight,
                        float elec_0_NOMINAL_EleEffSF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_TightLLH_d0z0_v13_isolFCTight, 
                        float elec_0_NOMINAL_EleEffSF_offline_RecoTrk,
                        float elec_0_NOMINAL_EleEffSF_offline_TightLLH_d0z0_v13,
                        float elec_0_NOMINAL_efficiency_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_TightLLH_d0z0_v13_isolFCTight)
{
    elecSF ret;
    ret.elec_0_NOMINAL_EleEffSF_Isolation_TightLLH_d0z0_v13_FCLoose =elec_0_NOMINAL_EleEffSF_Isolation_TightLLH_d0z0_v13_FCLoose;
    ret.elec_0_NOMINAL_EleEffSF_Isolation_TightLLH_d0z0_v13_FCTight =elec_0_NOMINAL_EleEffSF_Isolation_TightLLH_d0z0_v13_FCTight;
    ret.elec_0_NOMINAL_EleEffSF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_TightLLH_d0z0_v13_isolFCTight=elec_0_NOMINAL_EleEffSF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_TightLLH_d0z0_v13_isolFCTight;
    ret.elec_0_NOMINAL_EleEffSF_offline_RecoTrk=elec_0_NOMINAL_EleEffSF_offline_RecoTrk;
    ret.elec_0_NOMINAL_EleEffSF_offline_TightLLH_d0z0_v13=elec_0_NOMINAL_EleEffSF_offline_TightLLH_d0z0_v13;
    ret.elec_0_NOMINAL_efficiency_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_TightLLH_d0z0_v13_isolFCTight=elec_0_NOMINAL_efficiency_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_TightLLH_d0z0_v13_isolFCTight;
    return ret;
} 

struct muonSF getMuonSF(float muon_0_NOMINAL_MuEffSF_HLT_mu20_iloose_L1MU15_OR_HLT_mu50_QualMedium, 
                        float muon_0_NOMINAL_MuEffSF_HLT_mu26_ivarmedium_OR_HLT_mu50_QualMedium, 
                        float muon_0_NOMINAL_MuEffSF_IsoFCLoose, 
                        float muon_0_NOMINAL_MuEffSF_IsoFCLoose_FixedRad, 
                        float muon_0_NOMINAL_MuEffSF_IsoFCTight, 
                        float muon_0_NOMINAL_MuEffSF_IsoFCTightTrackOnly, 
                        float muon_0_NOMINAL_MuEffSF_IsoFCTightTrackOnly_FixedRad, 
                        float muon_0_NOMINAL_MuEffSF_IsoFCTight_FixedRad, 
                        float muon_0_NOMINAL_MuEffSF_IsoFixedCutHighPtTrackOnly, 
                        float muon_0_NOMINAL_MuEffSF_Reco_QualMedium, 
                        float muon_0_NOMINAL_MuEffSF_TTVA){
    muonSF ret;
    ret.muon_0_NOMINAL_MuEffSF_HLT_mu20_iloose_L1MU15_OR_HLT_mu50_QualMedium = muon_0_NOMINAL_MuEffSF_HLT_mu20_iloose_L1MU15_OR_HLT_mu50_QualMedium;
    ret.muon_0_NOMINAL_MuEffSF_HLT_mu26_ivarmedium_OR_HLT_mu50_QualMedium = muon_0_NOMINAL_MuEffSF_HLT_mu26_ivarmedium_OR_HLT_mu50_QualMedium;
    ret.muon_0_NOMINAL_MuEffSF_IsoFCLoose = muon_0_NOMINAL_MuEffSF_IsoFCLoose;
    ret.muon_0_NOMINAL_MuEffSF_IsoFCLoose_FixedRad = muon_0_NOMINAL_MuEffSF_IsoFCLoose_FixedRad;
    ret.muon_0_NOMINAL_MuEffSF_IsoFCTight = muon_0_NOMINAL_MuEffSF_IsoFCTight;
    ret.muon_0_NOMINAL_MuEffSF_IsoFCTightTrackOnly = muon_0_NOMINAL_MuEffSF_IsoFCTightTrackOnly;
    ret.muon_0_NOMINAL_MuEffSF_IsoFCTightTrackOnly_FixedRad = muon_0_NOMINAL_MuEffSF_IsoFCTightTrackOnly_FixedRad;
    ret.muon_0_NOMINAL_MuEffSF_IsoFCTight_FixedRad = muon_0_NOMINAL_MuEffSF_IsoFCTight_FixedRad;
    ret.muon_0_NOMINAL_MuEffSF_IsoFixedCutHighPtTrackOnly = muon_0_NOMINAL_MuEffSF_IsoFixedCutHighPtTrackOnly;
    ret.muon_0_NOMINAL_MuEffSF_Reco_QualMedium = muon_0_NOMINAL_MuEffSF_Reco_QualMedium;
    ret.muon_0_NOMINAL_MuEffSF_TTVA = muon_0_NOMINAL_MuEffSF_TTVA;
    return ret;
}

    
struct tauHadSF getTauHadSF(float tau_0_NOMINAL_TauEffSF_JetRNNloose, 
                            float tau_0_NOMINAL_TauEffSF_JetRNNmedium, 
                            float tau_0_NOMINAL_TauEffSF_JetRNNtight, 
                            float tau_0_NOMINAL_TauEffSF_LooseEleBDT_electron, 
                            float tau_0_NOMINAL_TauEffSF_MediumEleBDT_electron, 
                            float tau_0_NOMINAL_TauEffSF_reco)
{
    tauHadSF ret;
    ret.tau_0_NOMINAL_TauEffSF_JetRNNloose= tau_0_NOMINAL_TauEffSF_JetRNNloose;
    ret.tau_0_NOMINAL_TauEffSF_JetRNNmedium= tau_0_NOMINAL_TauEffSF_JetRNNmedium;
    ret.tau_0_NOMINAL_TauEffSF_JetRNNtight= tau_0_NOMINAL_TauEffSF_JetRNNtight;
    ret.tau_0_NOMINAL_TauEffSF_LooseEleBDT_electron= tau_0_NOMINAL_TauEffSF_LooseEleBDT_electron;
    ret.tau_0_NOMINAL_TauEffSF_MediumEleBDT_electron= tau_0_NOMINAL_TauEffSF_MediumEleBDT_electron;
    ret.tau_0_NOMINAL_TauEffSF_reco= tau_0_NOMINAL_TauEffSF_reco;
    return ret;
}

struct jetSF getJetSF(  float jet_NOMINAL_central_jets_global_effSF_JVT,
                        float jet_NOMINAL_central_jets_global_ineffSF_JVT,
                        float jet_NOMINAL_forward_jets_global_effSF_JVT,
                        float jet_NOMINAL_forward_jets_global_ineffSF_JVT,
                        float jet_NOMINAL_global_effSF_MV2c10_FixedCutBEff_85,
                        float jet_NOMINAL_global_ineffSF_MV2c10_FixedCutBEff_85)
{
    jetSF ret;
    ret.jet_NOMINAL_central_jets_global_effSF_JVT = jet_NOMINAL_central_jets_global_effSF_JVT;
    ret.jet_NOMINAL_central_jets_global_ineffSF_JVT = jet_NOMINAL_central_jets_global_ineffSF_JVT;
    ret.jet_NOMINAL_forward_jets_global_effSF_JVT = jet_NOMINAL_forward_jets_global_effSF_JVT;
    ret.jet_NOMINAL_forward_jets_global_ineffSF_JVT = jet_NOMINAL_forward_jets_global_ineffSF_JVT;
    ret.jet_NOMINAL_global_effSF_MV2c10_FixedCutBEff_85 = jet_NOMINAL_global_effSF_MV2c10_FixedCutBEff_85;
    ret.jet_NOMINAL_global_ineffSF_MV2c10_FixedCutBEff_85 = jet_NOMINAL_global_ineffSF_MV2c10_FixedCutBEff_85;
    return ret;
}
   
bool elecSel(ROOT::Math::PtEtaPhiMVector vec, uint isolation, int id, elecTriggers trigs)
{
    return  vec.Pt() > 30 && 
            (bool)isolation && 
            (bool)id && 
            (trigs.eleTrigMatch_0_HLT_e24_lhmedium_L1EM20VH || 
            trigs.eleTrigMatch_0_HLT_e26_lhtight_nod0_ivarloose ||
            trigs.eleTrigMatch_0_HLT_e60_lhmedium || 
            trigs.eleTrigMatch_0_HLT_e60_lhmedium_nod0 ||
            trigs.eleTrigMatch_0_HLT_e120_lhloose ||
            trigs.eleTrigMatch_0_HLT_e140_lhloose_nod0);
}


bool muonSel(ROOT::Math::PtEtaPhiMVector vec, uint isolation, int id, muonTriggers trigs)
{
    return  vec.Pt()>30 && 
            (bool)isolation && 
            (bool)id && 
            (trigs.muTrigMatch_0_HLT_mu20_iloose_L1MU15 || 
            trigs.muTrigMatch_0_HLT_mu26_ivarmedium || 
            trigs.muTrigMatch_0_HLT_mu50);
}

bool bjetSel(ROOT::Math::PtEtaPhiMVector vec, int btag)
{
    return vec.Pt() > 5 && (bool)btag;
}

bool tauHadSel(ROOT::Math::PtEtaPhiMVector vec, uint rnn){
    return vec.Pt()>5 && rnn;
}

double getElecSFNum(elecSF sf)
{
    return sf.elec_0_NOMINAL_EleEffSF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_TightLLH_d0z0_v13_isolFCTight;
}

double getMuonSFNum(muonSF sf)
{
    return  sf.muon_0_NOMINAL_MuEffSF_HLT_mu26_ivarmedium_OR_HLT_mu50_QualMedium * 
            sf.muon_0_NOMINAL_MuEffSF_HLT_mu20_iloose_L1MU15_OR_HLT_mu50_QualMedium * 
            sf.muon_0_NOMINAL_MuEffSF_IsoFCTightTrackOnly_FixedRad;
}

double getJetSFNum(jetSF sf)
{
    return sf.jet_NOMINAL_global_effSF_MV2c10_FixedCutBEff_85;
}

double getTauHadSFNum(tauHadSF sf)
{
    return sf.tau_0_NOMINAL_TauEffSF_JetRNNtight;
}

double getTotalSF(double lepsf, double tausf, double jetsf)
{
    if (lepsf == 0 || tausf==0 || jetsf == 0) std::cout<<"error: zero scale fac\n";
    return lepsf * tausf * jetsf;
}