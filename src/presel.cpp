//#define DEBUG
#include "presel.h"

void presel(const char* infile, const char* outfile)
{
    #ifndef DEBUG
    ROOT::EnableImplicitMT();
    #endif
    ROOT::RDataFrame df("NOMINAL", infile);


    std::vector<std::string> elec_triggers_name = { "eleTrigMatch_0_HLT_e120_lhloose", 
                                                    "eleTrigMatch_0_HLT_e140_lhloose_nod0", 
                                                    "eleTrigMatch_0_HLT_e24_lhmedium_L1EM20VH", 
                                                    "eleTrigMatch_0_HLT_e26_lhtight_nod0_ivarloose", 
                                                    "eleTrigMatch_0_HLT_e60_lhmedium", 
                                                    "eleTrigMatch_0_HLT_e60_lhmedium_nod0" };

    std::vector<std::string> muon_triggers_name = { "muTrigMatch_0_HLT_mu20_iloose_L1MU15", 
                                                    "muTrigMatch_0_HLT_mu26_ivarmedium", 
                                                    "muTrigMatch_0_HLT_mu50" };

    std::vector<std::string> elec_sf_name       = { "elec_0_NOMINAL_EleEffSF_Isolation_TightLLH_d0z0_v13_FCLoose", 
                                                    "elec_0_NOMINAL_EleEffSF_Isolation_TightLLH_d0z0_v13_FCTight", 
                                                    "elec_0_NOMINAL_EleEffSF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_TightLLH_d0z0_v13_isolFCTight", 
                                                    "elec_0_NOMINAL_EleEffSF_offline_RecoTrk", 
                                                    "elec_0_NOMINAL_EleEffSF_offline_TightLLH_d0z0_v13", 
                                                    "elec_0_NOMINAL_efficiency_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_TightLLH_d0z0_v13_isolFCTight"};
    
    std::vector<std::string> muon_sf_name       = { "muon_0_NOMINAL_MuEffSF_HLT_mu20_iloose_L1MU15_OR_HLT_mu50_QualMedium", 
                                                    "muon_0_NOMINAL_MuEffSF_HLT_mu26_ivarmedium_OR_HLT_mu50_QualMedium", 
                                                    "muon_0_NOMINAL_MuEffSF_IsoFCLoose", 
                                                    "muon_0_NOMINAL_MuEffSF_IsoFCLoose_FixedRad", 
                                                    "muon_0_NOMINAL_MuEffSF_IsoFCTight", 
                                                    "muon_0_NOMINAL_MuEffSF_IsoFCTightTrackOnly", 
                                                    "muon_0_NOMINAL_MuEffSF_IsoFCTightTrackOnly_FixedRad", 
                                                    "muon_0_NOMINAL_MuEffSF_IsoFCTight_FixedRad", 
                                                    "muon_0_NOMINAL_MuEffSF_IsoFixedCutHighPtTrackOnly", 
                                                    "muon_0_NOMINAL_MuEffSF_Reco_QualMedium", 
                                                    "muon_0_NOMINAL_MuEffSF_TTVA" };

    std::vector<std::string> tau_had_sf_name    = { "tau_0_NOMINAL_TauEffSF_JetRNNloose", 
                                                    "tau_0_NOMINAL_TauEffSF_JetRNNmedium", 
                                                    "tau_0_NOMINAL_TauEffSF_JetRNNtight", 
                                                    "tau_0_NOMINAL_TauEffSF_LooseEleBDT_electron", 
                                                    "tau_0_NOMINAL_TauEffSF_MediumEleBDT_electron", 
                                                    "tau_0_NOMINAL_TauEffSF_reco" };

    std::vector<std::string> jet_sf_name        = { "jet_NOMINAL_central_jets_global_effSF_JVT", 
                                                    "jet_NOMINAL_central_jets_global_ineffSF_JVT", 
                                                    "jet_NOMINAL_forward_jets_global_effSF_JVT", 
                                                    "jet_NOMINAL_forward_jets_global_ineffSF_JVT", 
                                                    "jet_NOMINAL_global_effSF_MV2c10_FixedCutBEff_85", 
                                                    "jet_NOMINAL_global_ineffSF_MV2c10_FixedCutBEff_85" };

    
    if (df.HasColumn("weight_mc")){
        auto ef = df
                    .Define("elec_trigs",   getElecTriggers,    elec_triggers_name)
                    .Define("muon_trigs",   getMuonTriggers,    muon_triggers_name)

                    .Define("elec_pass",    elecSel,            {"elec_0_p4_n", "elec_0_iso_FCTight", "elec_0_id_tight", "elec_trigs"})
                    .Define("muon_pass",    muonSel,            {"muon_0_p4_n", "muon_0_iso_FCTightTrackOnly_FixedRad", "muon_0_id_tight", "muon_trigs"})
                    .Define("tau_had_pass", tauHadSel,          {"tau_0_p4_n", "tau_0_jet_rnn_tight"})
                    .Define("bjet_0_pass",  bjetSel,            {"bjet_0_p4_n", "bjet_0_b_tagged_MV2c10_FixedCutBEff_85"})
                    .Define("bjet_1_pass",  bjetSel,            {"bjet_1_p4_n", "bjet_1_b_tagged_MV2c10_FixedCutBEff_85"})

                    .Filter([](bool elec, bool muon){return elec != muon;}, {"elec_pass", "muon_pass"})
                    .Filter([](bool bjet_0, bool bjet_1){return bjet_0 && bjet_1;}, {"bjet_0_pass", "bjet_1_pass"})
                    .Filter([](bool tau_had){return tau_had;}, {"tau_had_pass"})

                    .Define("elec_sfs",     getElecSF,          elec_sf_name)
                    .Define("muon_sfs",     getMuonSF,          muon_sf_name)
                    .Define("tau_had_sfs",  getTauHadSF,        tau_had_sf_name)
                    .Define("jet_sfs",      getJetSF,           jet_sf_name)

                    .Define("elec_sf",      getElecSFNum,       {"elec_sfs"})
                    .Define("muon_sf",      getMuonSFNum,       {"muon_sfs"})
                    .Define("tau_had_sf",   getTauHadSFNum,     {"tau_had_sfs"})
                    .Define("jet_sf",       getJetSFNum,        {"jet_sfs"})

                    .Define("lep_0_p4_n",   [](ROOT::Math::PtEtaPhiMVector ele, ROOT::Math::PtEtaPhiMVector mu, bool ele_pass){return ele_pass ? ele : mu;}, {"elec_0_p4_n", "muon_0_p4_n", "elec_pass"})
                    .Define("lep_sf",       [](double ele, double mu, bool ele_pass){return ele_pass ? ele : mu;},  {"elec_sf", "muon_sf", "elec_pass"})
                    .Define("total_sf",     getTotalSF, {"lep_sf", "tau_had_sf", "jet_sf"})
                    ;

        all_branch_name_mc.push_back("total_sf");
        all_branch_name_mc.push_back("lep_0_p4_n");
        all_branch_name_mc.push_back("muon_pass");
        all_branch_name_mc.push_back("elec_pass");
        all_branch_name_mc.push_back("elec_sf");
        all_branch_name_mc.push_back("muon_sf" );
        all_branch_name_mc.push_back("tau_had_sf");
        all_branch_name_mc.push_back("jet_sf");

        

        ef.Snapshot("NOMINAL", outfile, all_branch_name_mc);
    }
    else 
    {
        auto ef = df
                    .Define("elec_trigs",   getElecTriggers,    elec_triggers_name)
                    .Define("muon_trigs",   getMuonTriggers,    muon_triggers_name)

                    .Define("elec_pass",    elecSel,            {"elec_0_p4_n", "elec_0_iso_FCTight", "elec_0_id_tight", "elec_trigs"})
                    .Define("muon_pass",    muonSel,            {"muon_0_p4_n", "muon_0_iso_FCTightTrackOnly_FixedRad", "muon_0_id_tight", "muon_trigs"})
                    .Define("tau_had_pass", tauHadSel,          {"tau_0_p4_n", "tau_0_jet_rnn_tight"})
                    .Define("bjet_0_pass",  bjetSel,            {"bjet_0_p4_n", "bjet_0_b_tagged_MV2c10_FixedCutBEff_85"})
                    .Define("bjet_1_pass",  bjetSel,            {"bjet_1_p4_n", "bjet_1_b_tagged_MV2c10_FixedCutBEff_85"})

                    .Filter([](bool elec, bool muon){return elec != muon;}, {"elec_pass", "muon_pass"})
                    .Filter([](bool bjet_0, bool bjet_1){return bjet_0 && bjet_1;}, {"bjet_0_pass", "bjet_1_pass"})
                    .Filter([](bool tau_had){return tau_had;}, {"tau_had_pass"})
                    .Define("lep_0_p4_n",   [](ROOT::Math::PtEtaPhiMVector ele, ROOT::Math::PtEtaPhiMVector mu, bool ele_pass){return ele_pass ? ele : mu;}, {"elec_0_p4_n", "muon_0_p4_n", "elec_pass"})
                    ;

        all_branch_name_data.push_back("lep_0_p4_n");
        all_branch_name_mc.push_back("muon_pass");
        all_branch_name_mc.push_back("elec_pass");

        ef.Snapshot("NOMINAL", outfile, all_branch_name_data);
    }
}
#ifdef DEBUG

int main()
{
    presel("output/01_convert_out/410470.txt.root","debug.root");
    return 0;
}

#else

int main(int argc, char** argv)
{
    presel(argv[1], argv[2]);
    copyMeta(argv[1], argv[2]);
    return 0;
}

#endif