
#include "genAna.h"

// #define MT
#define DEBUG
ROOT::Math::PtEtaPhiMVector TLV2MLV(TLorentzVector v)
{
    ROOT::Math::PtEtaPhiMVector ret(v.Pt(), v.Eta(), v.Phi(), v.M());
    return ret;
}

std::vector<std::string> getFileName(std::string path){
    std::vector<std::string> ret;
    for (const auto & entry : std::filesystem::directory_iterator(path))
    {
        if(entry.path() != path+"/meta.txt")
            ret.push_back(entry.path());
    }
    return ret;
}

void convert(std::vector<std::string> in_filelist, char *out_file)
{
#ifdef MT
    ROOT::EnableImplicitMT(4);
#endif
    ROOT::RDataFrame df("NOMINAL", in_filelist);
    if(df.HasColumn("weight_mc"))
    {
        auto ef = df.Filter([](TLorentzVector v){return v.Pt() > 0.1;}, {"bjet_0_p4"})
                    .Filter([](TLorentzVector v){return v.Pt() > 0.1;}, {"bjet_1_p4"})
                    .Define("bjet_0_p4_n", TLV2MLV, {"bjet_0_p4"})
                    .Define("bjet_0_matched_p4_n", TLV2MLV, {"bjet_0_matched_p4"})
                    .Define("bjet_1_p4_n", TLV2MLV, {"bjet_1_p4"})
                    .Define("bjet_1_matched_p4_n", TLV2MLV, {"bjet_1_matched_p4"})
                    .Define("elec_0_p4_n", TLV2MLV, {"elec_0_p4"})
                    .Define("elec_0_matched_p4_n", TLV2MLV, {"elec_0_matched_p4"})
                    .Define("ljet_0_p4_n", TLV2MLV, {"ljet_0_p4"})
                    .Define("ljet_0_matched_p4_n", TLV2MLV, {"ljet_0_matched_p4"})
                    .Define("ljet_1_p4_n", TLV2MLV, {"ljet_1_p4"})
                    .Define("ljet_1_matched_p4_n", TLV2MLV, {"ljet_1_matched_p4"})
                    .Define("ljet_2_p4_n", TLV2MLV, {"ljet_2_p4"})
                    .Define("ljet_2_matched_p4_n", TLV2MLV, {"ljet_2_matched_p4"})
                    .Define("ljet_3_p4_n", TLV2MLV, {"ljet_3_p4"})
                    .Define("ljet_3_matched_p4_n", TLV2MLV, {"ljet_3_matched_p4"})
                    .Define("ljet_4_p4_n", TLV2MLV, {"ljet_4_p4"})
                    .Define("ljet_4_matched_p4_n", TLV2MLV, {"ljet_4_matched_p4"})
                    .Define("met_reco_p4_n", TLV2MLV, {"met_reco_p4"})
                    .Define("muon_0_p4_n", TLV2MLV, {"muon_0_p4"})
                    .Define("muon_0_matched_p4_n", TLV2MLV, {"muon_0_matched_p4"})
                    .Define("tau_0_p4_n", TLV2MLV, {"tau_0_p4"})
                    .Define("tau_0_truth_p4_n", TLV2MLV, {"tau_0_truth_p4"})
                    .Define("tau_0_truth_vis_charged_p4_n", TLV2MLV, {"tau_0_truth_vis_charged_p4"})
                    .Define("tau_0_truth_vis_neutral_others_p4_n", TLV2MLV, {"tau_0_truth_vis_neutral_others_p4"})
                    .Define("tau_0_truth_vis_neutral_p4_n", TLV2MLV, {"tau_0_truth_vis_neutral_p4"})
                    .Define("tau_0_truth_vis_neutral_pions_p4_n", TLV2MLV, {"tau_0_truth_vis_neutral_pions_p4"});

        std::vector<std::string> snap_list
                        {"HLT_e120_lhloose", "HLT_e140_lhloose_nod0", "HLT_e24_lhmedium_L1EM20VH", "HLT_e26_lhtight_nod0_ivarloose", 
                        "HLT_e60_lhmedium", "HLT_e60_lhmedium_nod0", "HLT_mu20_iloose_L1MU15", "HLT_mu26_ivarmedium", "HLT_mu50", 

                        "NOMINAL_pileup_combined_weight", "NOMINAL_pileup_random_run_number", 
                        
                        "bjet_0", "bjet_0_b_tagged_MV2c10_FixedCutBEff_85", 
                        "bjet_0_matched", "bjet_0_matched_classifierParticleOrigin", "bjet_0_matched_classifierParticleType", "bjet_0_matched_origin", 
                        "bjet_0_matched_p4_n", "bjet_0_matched_pdgId", "bjet_0_matched_pz", "bjet_0_matched_q", "bjet_0_matched_status", "bjet_0_matched_type", 
                        "bjet_0_origin", "bjet_0_p4_n", "bjet_0_type", 
                        
                        "bjet_1", "bjet_1_b_tagged_MV2c10_FixedCutBEff_85", "bjet_1_matched", 
                        "bjet_1_matched_classifierParticleOrigin", "bjet_1_matched_classifierParticleType", "bjet_1_matched_origin", 
                        "bjet_1_matched_p4_n", "bjet_1_matched_pdgId", "bjet_1_matched_pz", "bjet_1_matched_q", "bjet_1_matched_status", 
                        "bjet_1_matched_type", "bjet_1_origin", "bjet_1_p4_n", "bjet_1_type", 
                        
                        "eleTrigMatch_0_HLT_e120_lhloose", 
                        "eleTrigMatch_0_HLT_e140_lhloose_nod0", "eleTrigMatch_0_HLT_e24_lhmedium_L1EM20VH", "eleTrigMatch_0_HLT_e26_lhtight_nod0_ivarloose", 
                        "eleTrigMatch_0_HLT_e60_lhmedium", "eleTrigMatch_0_HLT_e60_lhmedium_nod0", "eleTrigMatch_0_trigger_matched", 
                        
                        "elec_0", 
                        "elec_0_NOMINAL_EleEffSF_Isolation_TightLLH_d0z0_v13_FCLoose", "elec_0_NOMINAL_EleEffSF_Isolation_TightLLH_d0z0_v13_FCTight", 
                        "elec_0_NOMINAL_EleEffSF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_TightLLH_d0z0_v13_isolFCTight", 
                        "elec_0_NOMINAL_EleEffSF_offline_RecoTrk", "elec_0_NOMINAL_EleEffSF_offline_TightLLH_d0z0_v13", 
                        "elec_0_NOMINAL_efficiency_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_TightLLH_d0z0_v13_isolFCTight", 
                        
                        "elec_0_cluster_eta", "elec_0_cluster_eta_be2", "elec_0_id_medium", "elec_0_id_tight", "elec_0_iso_FCLoose", "elec_0_iso_FCLoose_FixedRad", "elec_0_iso_FCTight", 
                        "elec_0_iso_FCTightTrackOnly_FixedRad", "elec_0_iso_FixedCutLoose", "elec_0_iso_FixedCutTight", "elec_0_iso_FixedCutTightCaloOnly", 
                        
                        "elec_0_matched_classifierParticleOrigin", 
                        "elec_0_matched_classifierParticleType", "elec_0_matched_origin", "elec_0_matched_p4_n", "elec_0_matched_pdgId", "elec_0_matched_pz", "elec_0_matched_q", "elec_0_matched_type", 
                        "elec_0_p4_n", "elec_0_q", "elec_0_trk_d0_sig", "elec_0_trk_pvx_z0_sintheta", "elec_0_trk_z0_sintheta", 
                        
                        "event_is_bad_batman", "event_number", 
                        
                        "jet_NOMINAL_central_jets_global_effSF_JVT", 
                        "jet_NOMINAL_central_jets_global_ineffSF_JVT", "jet_NOMINAL_forward_jets_global_effSF_JVT", "jet_NOMINAL_forward_jets_global_ineffSF_JVT", "jet_NOMINAL_global_effSF_MV2c10_FixedCutBEff_85", 
                        "jet_NOMINAL_global_ineffSF_MV2c10_FixedCutBEff_85", 
                        
                        "ljet_0", "ljet_0_b_tagged_MV2c10_FixedCutBEff_85", "ljet_0_flavorlabel", "ljet_0_flavorlabel_cone", "ljet_0_flavorlabel_part", 
                        "ljet_0_matched", "ljet_0_matched_classifierParticleOrigin", "ljet_0_matched_classifierParticleType", "ljet_0_matched_origin", "ljet_0_matched_p4_n", "ljet_0_matched_pdgId", 
                        "ljet_0_matched_pz", "ljet_0_matched_q", "ljet_0_matched_status", "ljet_0_matched_type", "ljet_0_origin", "ljet_0_p4_n", "ljet_0_type", 
                        
                        "ljet_1", "ljet_1_b_tagged_MV2c10_FixedCutBEff_85", 
                        "ljet_1_flavorlabel", "ljet_1_flavorlabel_cone", "ljet_1_flavorlabel_part", "ljet_1_matched", "ljet_1_matched_classifierParticleOrigin", "ljet_1_matched_classifierParticleType",
                        "ljet_1_matched_origin", "ljet_1_matched_p4_n", "ljet_1_matched_pdgId", "ljet_1_matched_pz", "ljet_1_matched_q", "ljet_1_matched_status", "ljet_1_matched_type", "ljet_1_origin", 
                        "ljet_1_p4_n", "ljet_1_type", 
                        
                        "ljet_2", "ljet_2_b_tag_quantile", "ljet_2_b_tag_score", "ljet_2_b_tagged_MV2c10_FixedCutBEff_85", "ljet_2_cleanJet_EC_LooseBad", "ljet_2_fjvt", 
                        "ljet_2_flavorlabel", "ljet_2_flavorlabel_cone", "ljet_2_flavorlabel_part", "ljet_2_is_Jvt_HS", "ljet_2_jvt", "ljet_2_matched", "ljet_2_matched_classifierParticleOrigin", 
                        "ljet_2_matched_classifierParticleType", "ljet_2_matched_mother_pdgId", "ljet_2_matched_mother_status", "ljet_2_matched_origin", "ljet_2_matched_p4_n", "ljet_2_matched_pdgId", 
                        "ljet_2_matched_pz", "ljet_2_matched_q", "ljet_2_matched_status", "ljet_2_matched_type", "ljet_2_origin", "ljet_2_p4_n", "ljet_2_q", "ljet_2_type", "ljet_2_width", 
                        
                        "ljet_3", 
                        "ljet_3_b_tag_quantile", "ljet_3_b_tag_score", "ljet_3_b_tagged_MV2c10_FixedCutBEff_85", "ljet_3_cleanJet_EC_LooseBad", "ljet_3_fjvt", "ljet_3_flavorlabel", "ljet_3_flavorlabel_cone", 
                        "ljet_3_flavorlabel_part", "ljet_3_is_Jvt_HS", "ljet_3_jvt", "ljet_3_matched", "ljet_3_matched_classifierParticleOrigin", "ljet_3_matched_classifierParticleType", "ljet_3_matched_mother_pdgId", 
                        "ljet_3_matched_mother_status", "ljet_3_matched_origin", "ljet_3_matched_p4_n", "ljet_3_matched_pdgId", "ljet_3_matched_pz", "ljet_3_matched_q", "ljet_3_matched_status", 
                        "ljet_3_matched_type", "ljet_3_origin", "ljet_3_p4_n", "ljet_3_q", "ljet_3_type", "ljet_3_width", 
                        
                        "ljet_4", "ljet_4_b_tag_quantile", "ljet_4_b_tag_score", "ljet_4_b_tagged_MV2c10_FixedCutBEff_85", 
                        "ljet_4_cleanJet_EC_LooseBad", "ljet_4_fjvt", "ljet_4_flavorlabel", "ljet_4_flavorlabel_cone", "ljet_4_flavorlabel_part", "ljet_4_is_Jvt_HS", "ljet_4_jvt", "ljet_4_matched", 
                        "ljet_4_matched_classifierParticleOrigin", "ljet_4_matched_classifierParticleType", "ljet_4_matched_mother_pdgId", "ljet_4_matched_mother_status", "ljet_4_matched_origin", 
                        "ljet_4_matched_p4_n", "ljet_4_matched_pdgId", "ljet_4_matched_pz", "ljet_4_matched_q", "ljet_4_matched_status", "ljet_4_matched_type", "ljet_4_origin", "ljet_4_p4_n", "ljet_4_q", 
                        "ljet_4_type", "ljet_4_width", 
                        
                        "met_reco_p4_n", 
                        
                        "muTrigMatch_0_HLT_mu20_iloose_L1MU15", "muTrigMatch_0_HLT_mu26_ivarmedium", "muTrigMatch_0_HLT_mu50", "muTrigMatch_0_trigger_matched", 
                        
                        "muon_0", "muon_0_NOMINAL_MuEffSF_HLT_mu20_iloose_L1MU15_OR_HLT_mu50_QualMedium", "muon_0_NOMINAL_MuEffSF_HLT_mu26_ivarmedium_OR_HLT_mu50_QualMedium", "muon_0_NOMINAL_MuEffSF_IsoFCLoose", 
                        "muon_0_NOMINAL_MuEffSF_IsoFCLoose_FixedRad", "muon_0_NOMINAL_MuEffSF_IsoFCTight", "muon_0_NOMINAL_MuEffSF_IsoFCTightTrackOnly", "muon_0_NOMINAL_MuEffSF_IsoFCTightTrackOnly_FixedRad", 
                        "muon_0_NOMINAL_MuEffSF_IsoFCTight_FixedRad", "muon_0_NOMINAL_MuEffSF_IsoFixedCutHighPtTrackOnly", "muon_0_NOMINAL_MuEffSF_Reco_QualMedium", "muon_0_NOMINAL_MuEffSF_TTVA", 
                        
                        "muon_0_id_medium", "muon_0_id_tight", "muon_0_iso_FCLoose", "muon_0_iso_FCLoose_FixedRad", "muon_0_iso_FCTight", "muon_0_iso_FCTightTrackOnly_FixedRad", "muon_0_iso_FixedCutLoose", "muon_0_iso_FixedCutTight", 
                        "muon_0_iso_FixedCutTightCaloOnly", 
                        
                        "muon_0_matched_classifierParticleOrigin", "muon_0_matched_classifierParticleType", "muon_0_matched_p4_n", "muon_0_matched_pdgId", "muon_0_matched_pz", 
                        "muon_0_matched_q", "muon_0_matched_type", "muon_0_p4_n", "muon_0_q", "muon_0_trk_d0_sig", "muon_0_trk_pvx_z0_sig", "muon_0_trk_pvx_z0_sintheta", "muon_0_trk_z0_sintheta", 
                        
                        "n_actual_int", 
                        "n_actual_int_cor", "n_avg_int", "n_avg_int_cor", "n_bjets_MV2c10_FixedCutBEff_85", "n_electrons", "n_electrons_olr", "n_jets", "n_muons", "n_pvx", "n_taus", "n_taus_met", "n_taus_olr", 
                        "n_taus_rnn_loose", "n_taus_rnn_medium", "n_taus_rnn_tight", "n_taus_rnn_veryloose", "n_vx", 
                        
                        "run_number", 
                        
                        "tau_0", "tau_0_NOMINAL_TauEffSF_JetRNNloose", "tau_0_NOMINAL_TauEffSF_JetRNNmedium", 
                        "tau_0_NOMINAL_TauEffSF_JetRNNtight", "tau_0_NOMINAL_TauEffSF_LooseEleBDT_electron", "tau_0_NOMINAL_TauEffSF_MediumEleBDT_electron", "tau_0_NOMINAL_TauEffSF_reco", "tau_0_decay_mode", 
                        "tau_0_ele_bdt_eff_sf", "tau_0_ele_bdt_loose", "tau_0_ele_bdt_loose_retuned", "tau_0_ele_bdt_medium", "tau_0_ele_bdt_medium_retuned", "tau_0_ele_bdt_score", "tau_0_ele_bdt_score_retuned", 
                        "tau_0_ele_bdt_score_trans", "tau_0_ele_bdt_score_trans_retuned", "tau_0_ele_bdt_tight", "tau_0_ele_bdt_tight_retuned", "tau_0_ele_match_lhscore", "tau_0_ele_olr_pass", "tau_0_jetFakeFakeFlavour", 
                        "tau_0_jet_bdt_loose", "tau_0_jet_bdt_medium", "tau_0_jet_bdt_score", "tau_0_jet_bdt_score_trans", "tau_0_jet_bdt_tight", "tau_0_jet_bdt_veryloose", "tau_0_jet_rnn_loose", "tau_0_jet_rnn_medium", 
                        "tau_0_jet_rnn_score", "tau_0_jet_rnn_score_trans", "tau_0_jet_rnn_tight", "tau_0_n_charged_tracks", "tau_0_p4_n", "tau_0_q", 
                        
                        "tau_0_truth", "tau_0_truth_classifierParticleOrigin", 
                        "tau_0_truth_classifierParticleType", "tau_0_truth_isEle", "tau_0_truth_isHadTau", "tau_0_truth_isJet", "tau_0_truth_isMuon", "tau_0_truth_isTau", "tau_0_truth_isTruthMatch", 
                        "tau_0_truth_p4_n", "tau_0_truth_pdgId", "tau_0_truth_pz", "tau_0_truth_vis_charged_p4_n", "tau_0_truth_vis_neutral_others_p4_n", "tau_0_truth_vis_neutral_p4_n", "tau_0_truth_vis_neutral_pions_p4_n", 
                        
                        "triggerSF_em_NOMINAL", 
                        
                        "weight_mc", "weight_mc_v", "weight_total" };
        ef.Snapshot("NOMINAL", out_file, snap_list);
    }
    else
    {
        auto ef = df.Filter([](TLorentzVector v){return v.Pt() > 0.1;}, {"bjet_0_p4"})
                    .Filter([](TLorentzVector v){return v.Pt() > 0.1;}, {"bjet_1_p4"})
                    .Define("bjet_0_p4_n", TLV2MLV, {"bjet_0_p4"})
                    .Define("bjet_1_p4_n", TLV2MLV, {"bjet_1_p4"})
                    .Define("elec_0_p4_n", TLV2MLV, {"elec_0_p4"})
                    .Define("ljet_0_p4_n", TLV2MLV, {"ljet_0_p4"})
                    .Define("ljet_1_p4_n", TLV2MLV, {"ljet_1_p4"})
                    .Define("ljet_2_p4_n", TLV2MLV, {"ljet_2_p4"})
                    .Define("ljet_3_p4_n", TLV2MLV, {"ljet_3_p4"})
                    .Define("ljet_4_p4_n", TLV2MLV, {"ljet_4_p4"})
                    .Define("met_reco_p4_n", TLV2MLV, {"met_reco_p4"})
                    .Define("muon_0_p4_n", TLV2MLV, {"muon_0_p4"})
                    .Define("tau_0_p4_n", TLV2MLV, {"tau_0_p4"});
                    
        std::vector<std::string> snap_list
                    {"HLT_e120_lhloose", "HLT_e140_lhloose_nod0", "HLT_e24_lhmedium_L1EM20VH", "HLT_e26_lhtight_nod0_ivarloose", "HLT_e60_lhmedium", "HLT_e60_lhmedium_nod0", 
                    "HLT_mu20_iloose_L1MU15", "HLT_mu26_ivarmedium", "HLT_mu50", 
                    
                    "bjet_0", "bjet_0_b_tagged_MV2c10_FixedCutBEff_85", "bjet_0_origin", "bjet_0_p4_n", "bjet_0_type", 
                    
                    "bjet_1", "bjet_1_b_tagged_MV2c10_FixedCutBEff_85", "bjet_1_origin", "bjet_1_p4_n", "bjet_1_type", 
                    
                    "eleTrigMatch_0_HLT_e120_lhloose", "eleTrigMatch_0_HLT_e140_lhloose_nod0", 
                    "eleTrigMatch_0_HLT_e24_lhmedium_L1EM20VH", "eleTrigMatch_0_HLT_e26_lhtight_nod0_ivarloose", "eleTrigMatch_0_HLT_e60_lhmedium", "eleTrigMatch_0_HLT_e60_lhmedium_nod0", 
                    "eleTrigMatch_0_trigger_matched", 
                    
                    "elec_0", "elec_0_cluster_eta", "elec_0_cluster_eta_be2", "elec_0_id_medium", "elec_0_id_tight", "elec_0_iso_FCLoose", 
                    "elec_0_iso_FCLoose_FixedRad", "elec_0_iso_FCTight", "elec_0_iso_FCTightTrackOnly_FixedRad", "elec_0_iso_FixedCutLoose", "elec_0_iso_FixedCutTight", 
                    "elec_0_iso_FixedCutTightCaloOnly", "elec_0_p4_n", "elec_0_q", "elec_0_trk_d0_sig", "elec_0_trk_pvx_z0_sintheta", "elec_0_trk_z0_sintheta", 
                    
                    "event_is_bad_batman", "event_number", "lb_number", 
                    
                    "ljet_0", "ljet_0_b_tagged_MV2c10_FixedCutBEff_85", "ljet_0_flavorlabel", "ljet_0_flavorlabel_cone", "ljet_0_flavorlabel_part", "ljet_0_origin", "ljet_0_p4_n", "ljet_0_type", 
                    
                    "ljet_1", "ljet_1_b_tagged_MV2c10_FixedCutBEff_85", "ljet_1_flavorlabel", "ljet_1_flavorlabel_cone", "ljet_1_flavorlabel_part", "ljet_1_origin", "ljet_1_p4_n", "ljet_1_type", 
                    
                    "ljet_2", "ljet_2_b_tag_quantile", "ljet_2_b_tag_score", "ljet_2_b_tagged_MV2c10_FixedCutBEff_85", "ljet_2_cleanJet_EC_LooseBad", "ljet_2_fjvt", 
                    "ljet_2_flavorlabel", "ljet_2_flavorlabel_cone", "ljet_2_flavorlabel_part", "ljet_2_is_Jvt_HS", "ljet_2_jvt", "ljet_2_origin", "ljet_2_p4_n", "ljet_2_q", "ljet_2_type", 
                    "ljet_2_width", 
                    
                    "ljet_3", "ljet_3_b_tag_quantile", "ljet_3_b_tag_score", "ljet_3_b_tagged_MV2c10_FixedCutBEff_85", "ljet_3_cleanJet_EC_LooseBad", "ljet_3_fjvt", 
                    "ljet_3_flavorlabel", "ljet_3_flavorlabel_cone", "ljet_3_flavorlabel_part", "ljet_3_is_Jvt_HS", "ljet_3_jvt", "ljet_3_origin", "ljet_3_p4_n", "ljet_3_q", "ljet_3_type", 
                    "ljet_3_width", 
                    
                    "ljet_4", "ljet_4_b_tag_quantile", "ljet_4_b_tag_score", "ljet_4_b_tagged_MV2c10_FixedCutBEff_85", "ljet_4_cleanJet_EC_LooseBad", "ljet_4_fjvt", 
                    "ljet_4_flavorlabel", "ljet_4_flavorlabel_cone", "ljet_4_flavorlabel_part", "ljet_4_is_Jvt_HS", "ljet_4_jvt", "ljet_4_origin", "ljet_4_p4_n", "ljet_4_q", "ljet_4_type", 
                    "ljet_4_width", 
                    
                    "met_reco_p4_n", 
                    
                    "muTrigMatch_0_HLT_mu20_iloose_L1MU15", "muTrigMatch_0_HLT_mu26_ivarmedium", "muTrigMatch_0_HLT_mu50", "muTrigMatch_0_trigger_matched", 

                    "muon_0", "muon_0_id_medium", "muon_0_id_tight", "muon_0_iso_FCLoose", "muon_0_iso_FCLoose_FixedRad", "muon_0_iso_FCTight", "muon_0_iso_FCTightTrackOnly_FixedRad", 
                    "muon_0_iso_FixedCutLoose", "muon_0_iso_FixedCutTight", "muon_0_iso_FixedCutTightCaloOnly", "muon_0_p4_n", "muon_0_q", "muon_0_trk_d0_sig", "muon_0_trk_pvx_z0_sig", 
                    "muon_0_trk_pvx_z0_sintheta", "muon_0_trk_z0_sintheta", 
                    
                    "n_actual_int", "n_actual_int_cor", "n_avg_int", "n_avg_int_cor", "n_bjets_MV2c10_FixedCutBEff_85", 
                    "n_electrons", "n_electrons_olr", "n_jets", "n_muons", "n_pvx", "n_taus", "n_taus_met", "n_taus_olr", "n_taus_rnn_loose", "n_taus_rnn_medium", "n_taus_rnn_tight", 
                    "n_taus_rnn_veryloose", "n_vx", 
                    
                    "run_number", 
                    
                    "tau_0", "tau_0_decay_mode", "tau_0_ele_bdt_eff_sf", "tau_0_ele_bdt_loose", "tau_0_ele_bdt_loose_retuned", "tau_0_ele_bdt_medium", 
                    "tau_0_ele_bdt_medium_retuned", "tau_0_ele_bdt_score", "tau_0_ele_bdt_score_retuned", "tau_0_ele_bdt_score_trans", "tau_0_ele_bdt_score_trans_retuned", "tau_0_ele_bdt_tight", 
                    "tau_0_ele_bdt_tight_retuned", "tau_0_ele_match_lhscore", "tau_0_ele_olr_pass", "tau_0_jet_bdt_loose", "tau_0_jet_bdt_medium", "tau_0_jet_bdt_score", "tau_0_jet_bdt_score_trans", 
                    "tau_0_jet_bdt_tight", "tau_0_jet_bdt_veryloose", "tau_0_jet_rnn_loose", "tau_0_jet_rnn_medium", "tau_0_jet_rnn_score", "tau_0_jet_rnn_score_trans", "tau_0_jet_rnn_tight", 
                    "tau_0_n_charged_tracks", "tau_0_p4_n", "tau_0_q" };
        ef.Snapshot("NOMINAL", out_file, snap_list);
    }
    
    
}

void combineMeta(std::vector<std::string> in_filelist, const char *out_file)
{
    std::vector<TH1D> metas;
    for (auto entry:in_filelist){
        TFile* in_file = TFile::Open(&entry[0], "READ");
        TH1D tmp = *((TH1D*)in_file->Get("h_metadata"));
        metas.push_back(tmp);
        in_file->Close();
    }
    for(uint i=1; i < metas.size(); i++)
    {
        metas[0].Add(&metas[i]);
        std::cout<<metas[i].GetBinContent(6)<<"\n";
    }
    TFile* out = TFile::Open(out_file, "UPDATE");
    std::cout<<metas[0].GetBinContent(6)<<"\n";
    metas[0].Write();

    out->Close();
}

#ifndef DEBUG
int main(int argc, char **argv)
{
    auto fl = getFileName(argv[1]);
    convert(fl, argv[2]);
    combineMeta(fl, argv[2]);
}
#else
int main(int argc, char **argv)
{
    auto fl =getFileName("/mnt/NVME/HH/data/qidong.v24.mc16_13TeV.363356.Sherpa_221_NNPDF30NNLO_ZqqZll.deriv.DAOD_STDM4.e5525_e5984_s3126_r10724_r10726_p3975");
    convert(fl, "debug.root");
    combineMeta(fl, "debug.root");
}
#endif
