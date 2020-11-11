#include "genAna.h"
#include "GPUScaleB.cuh"

struct weight_rela
{
    double xspb;
    double kfac;
    double fileff;
    double sumofw;
};


struct weight_rela getWeightComps(const char* infile_name){
    std::string sentence(infile_name);
    std::istringstream iss(sentence);
    std::vector<std::string> tokens;
    std::string token;
    while (std::getline(iss, token, '.')) 
    {
        if (!token.empty())
            tokens.push_back(token);
    }
    std::string sentence2(tokens[0]);
    std::istringstream iss2(sentence2);
    std::vector<std::string> tokens2;
    std::string token2;
    while (std::getline(iss2, token2, '/')) 
    {
        if (!token2.empty())
            tokens2.push_back(token2);
    }
    ROOT::RDataFrame df("hehe", "/mnt/NVME/HH/output/weight_related.root");
    auto ef = df.Filter([token2](std::string s){return s.compare(token2) == 0;}, {"SampleID"});
    auto xspb = (ef.Take<double>("xspb").GetValue())[0];
    auto kfac = (ef.Take<double>("kfac").GetValue())[0];
    auto fileff = (ef.Take<double>("fileff").GetValue())[0];

    TFile* in = TFile::Open(infile_name, "READ");
    TH1D* tmp = (TH1D*)in->Get("h_metadata");
    double sumofw = tmp->GetBinContent(8);
    in->Close();
    
    return {xspb, kfac, fileff, sumofw};
}

enum lepType {MUON, ELEC, NO_LEP, MULTI_LEP};

enum vPos {INSIDE_HL_CLOSE_TO_H, INSIDE_HL_CLOSE_TO_L, OUTSIDE_HL_CLOSE_TO_H, OUTSIDE_HL_CLOSE_TO_L, FARAWAY_FROM_HL};

enum hypoTopCons {HAD_B0_LEP_B1, HAD_B1_LEP_B0};

double getMTauTauVis(ROOT::Math::PtEtaPhiMVector tau_had, ROOT::Math::PtEtaPhiMVector ele, ROOT::Math::PtEtaPhiMVector muon, uint is_ele, uint is_mu){
    double ret;
    if (is_mu && is_ele){
        //std::cout<<"ele and mu presented"<<std::endl;
        ret = -999;
    }
    else if(is_mu){
        ret =  (tau_had+muon).M();
    }
    else if(is_ele){
        ret =  (tau_had+ele).M();
    }
    return ret;
}

double getDeltaR(ROOT::Math::PtEtaPhiMVector a, ROOT::Math::PtEtaPhiMVector b)
{
    return std::sqrt((a.eta()-b.eta())*(a.eta()-b.eta()) + (a.phi()-b.phi())*(a.phi()-b.phi()));
}

double getdphi(double v1phi,  double v2phi)
{
    double ret;
    if(std::abs(v1phi - v2phi) > M_PI)
        ret = ((double)(v1phi > 0) -  (double)(v1phi <= 0)) * (2 * M_PI - std::abs(v1phi - v2phi));
    else
        ret = v2phi - v1phi;
    return ret;
}

uint getLepType(uint is_ele, uint is_mu)
{
    lepType ret;
    if (is_mu && is_ele ) ret = MULTI_LEP;
    else if(is_mu) ret = MUON;
    else if(is_ele) ret = ELEC;
    else ret = NO_LEP;
    return ret;
}

ROOT::Math::PtEtaPhiMVector getLepP4(uint type_uint, ROOT::Math::PtEtaPhiMVector elec, ROOT::Math::PtEtaPhiMVector muon)
{   
    lepType type = (lepType)type_uint;
    ROOT::Math::PtEtaPhiMVector lep(0,0,0,0);
    if(type == MUON) lep = muon;
    else if (type == ELEC) lep = elec;
    
    return lep;
}

uint getVPos(ROOT::Math::PtEtaPhiMVector met, ROOT::Math::PtEtaPhiMVector tau_had, ROOT::Math::PtEtaPhiMVector lep)
{
    vPos ret;
    double dphi_hl = getdphi(tau_had.phi(), lep.phi());
    double dphi_hv = getdphi(tau_had.phi(), met.phi());
    double dphi_lv = getdphi(lep.phi(),     met.phi());
    bool inside_hl = (dphi_hl * dphi_hv > 0) && (std::abs(dphi_hl) > std::abs(dphi_hv));
    bool close_to_h = std::abs(dphi_hv) < std::abs(dphi_lv);
    bool far_from_hl = (!inside_hl) && (std::abs(dphi_hv) > M_PI_2 || std::abs(dphi_lv) > M_PI_2);
    if (inside_hl)
    {
        if (close_to_h)         ret = INSIDE_HL_CLOSE_TO_H;
        else                    ret = INSIDE_HL_CLOSE_TO_L;
    }
    else
    {
        if(far_from_hl)         ret = FARAWAY_FROM_HL;
        else if (close_to_h)    ret = OUTSIDE_HL_CLOSE_TO_H;
        else                    ret = OUTSIDE_HL_CLOSE_TO_L;
    }
    return ret;
}

ROOT::Math::PtEtaPhiMVector getVHP4(uint v_pos_uint, ROOT::Math::PtEtaPhiMVector met, ROOT::Math::PtEtaPhiMVector tau_had, ROOT::Math::PtEtaPhiMVector lep)
{
    vPos v_pos = (vPos)v_pos_uint;
    ROOT::Math::PtEtaPhiMVector vh(0,0,0,0);
    double dphi_hv = getdphi(tau_had.phi(), met.phi());
    double dphi_hl = getdphi(tau_had.phi(), lep.phi());
    if (v_pos == OUTSIDE_HL_CLOSE_TO_H) 
    {
        vh = {met.Pt() * std::cos(std::abs(dphi_hv)), tau_had.Eta(), tau_had.Phi(), 0};
    }
    else if (v_pos == INSIDE_HL_CLOSE_TO_H || v_pos == INSIDE_HL_CLOSE_TO_L)
    {
        double vhpt = met.Pt() * std::cos(std::abs(dphi_hv)) - met.Pt() * std::sin(std::abs(dphi_hv)) * (1/std::tan(std::abs(dphi_hl)));
        vh = {vhpt, tau_had.Eta(), tau_had.Phi(), 0};
    }
    return vh;
}

ROOT::Math::PtEtaPhiMVector getVLP4(uint v_pos_uint, ROOT::Math::PtEtaPhiMVector met, ROOT::Math::PtEtaPhiMVector tau_had, ROOT::Math::PtEtaPhiMVector lep)
{
    vPos v_pos = (vPos)v_pos_uint;
    ROOT::Math::PtEtaPhiMVector vl(0,0,0,0);
    double dphi_lv = getdphi(lep.phi(),     met.phi());
    double dphi_hl = getdphi(tau_had.phi(), lep.phi());
    double dphi_hv = getdphi(tau_had.phi(), met.phi());
    if (v_pos == OUTSIDE_HL_CLOSE_TO_L) 
    {
        vl = {met.Pt() * std::cos(std::abs(dphi_lv)), lep.Eta(), lep.Phi(), 0};
    }
    else if (v_pos == INSIDE_HL_CLOSE_TO_H || v_pos == INSIDE_HL_CLOSE_TO_L)
    {
        double vlpt = met.Pt() * std::sin(std::abs(dphi_hv)) / std::sin(std::abs(dphi_hl));
        vl = {vlpt, lep.Eta(), lep.Phi(), 0};
    }
    return vl;
}

double getMTauTauColin(ROOT::Math::PtEtaPhiMVector tau_had, ROOT::Math::PtEtaPhiMVector lep, ROOT::Math::PtEtaPhiMVector vh, ROOT::Math::PtEtaPhiMVector vl, uint v_pos_uint)
{
    double ret;
    vPos v_pos = (vPos)v_pos_uint;
    if (v_pos==FARAWAY_FROM_HL) ret = -999;
    else ret = (vh+vl+tau_had+lep).M();
    return ret;
}

double getM4Body(ROOT::Math::PtEtaPhiMVector tau_had, ROOT::Math::PtEtaPhiMVector lep, ROOT::Math::PtEtaPhiMVector vh, ROOT::Math::PtEtaPhiMVector vl, ROOT::Math::PtEtaPhiMVector bjet0, ROOT::Math::PtEtaPhiMVector bjet1)
{
    double ret;
    ret = (vh+vl+tau_had+lep+bjet0+bjet1).M();
    return ret;
}

uint getTopType(ROOT::Math::PtEtaPhiMVector tau_had, ROOT::Math::PtEtaPhiMVector lep, ROOT::Math::PtEtaPhiMVector vh, ROOT::Math::PtEtaPhiMVector vl, ROOT::Math::PtEtaPhiMVector bjet0, ROOT::Math::PtEtaPhiMVector bjet1)
{
    hypoTopCons top_cons;
    double m_diff_tau_had_b0 = std::abs((vh+tau_had+bjet0).M() - TOP_MASS);
    double m_diff_tau_had_b1 = std::abs((vh+tau_had+bjet1).M() - TOP_MASS);
    double m_diff_tau_lep_b0 = std::abs((vl+lep+bjet0).M() - TOP_MASS);
    double m_diff_tau_lep_b1 = std::abs((vl+lep+bjet1).M() - TOP_MASS);

    double m_diff_had_0_lep_1 = m_diff_tau_had_b0 + m_diff_tau_lep_b1;
    double m_diff_had_1_lep_0 = m_diff_tau_had_b1 + m_diff_tau_lep_b0;

    if (m_diff_had_0_lep_1 > m_diff_had_1_lep_0)    top_cons = HAD_B1_LEP_B0;
    else                                            top_cons = HAD_B0_LEP_B1;
    return top_cons;
}

ROOT::Math::PtEtaPhiMVector getTop0P4(uint top_cons_uint, ROOT::Math::PtEtaPhiMVector tau_had, ROOT::Math::PtEtaPhiMVector lep, ROOT::Math::PtEtaPhiMVector vh, ROOT::Math::PtEtaPhiMVector vl, ROOT::Math::PtEtaPhiMVector bjet0, ROOT::Math::PtEtaPhiMVector bjet1)
{
    hypoTopCons top_cons = (hypoTopCons) top_cons_uint;
    ROOT::Math::PtEtaPhiMVector ret;

    if(top_cons == HAD_B1_LEP_B0)
        ret = vh+tau_had+bjet1;
    else
        ret = vh+tau_had+bjet0;

    return ret;
}

ROOT::Math::PtEtaPhiMVector getTop1P4(uint top_cons_uint, ROOT::Math::PtEtaPhiMVector tau_had, ROOT::Math::PtEtaPhiMVector lep, ROOT::Math::PtEtaPhiMVector vh, ROOT::Math::PtEtaPhiMVector vl, ROOT::Math::PtEtaPhiMVector bjet0, ROOT::Math::PtEtaPhiMVector bjet1)
{
    hypoTopCons top_cons = (hypoTopCons) top_cons_uint;
    ROOT::Math::PtEtaPhiMVector ret;

    if(top_cons == HAD_B1_LEP_B0)
        ret = vl+lep+bjet0;
    else
        ret = vl+lep+bjet1;

    return ret;
}

double getTop0DeltaRBTau(uint top_cons_uint, ROOT::Math::PtEtaPhiMVector tau_had, ROOT::Math::PtEtaPhiMVector bjet0, ROOT::Math::PtEtaPhiMVector bjet1)
{
    hypoTopCons top_cons = (hypoTopCons) top_cons_uint;
    double ret;
    if(top_cons == HAD_B1_LEP_B0) 
        ret = getDeltaR(tau_had, bjet1);
    else 
        ret = getDeltaR(tau_had, bjet0);
    return ret;
}

double getTop1DeltaRBTau(uint top_cons_uint, ROOT::Math::PtEtaPhiMVector lep, ROOT::Math::PtEtaPhiMVector bjet0, ROOT::Math::PtEtaPhiMVector bjet1)
{
    hypoTopCons top_cons = (hypoTopCons) top_cons_uint;
    double ret;
    if(top_cons == HAD_B1_LEP_B0) 
        ret = getDeltaR(lep, bjet0);
    else 
        ret = getDeltaR(lep, bjet1);
    return ret;
}

struct event{
    //4-vec
    ROOT::Math::PtEtaPhiMVector bjet0;
    ROOT::Math::PtEtaPhiMVector bjet1;
    ROOT::Math::PtEtaPhiMVector lep0;
    ROOT::Math::PtEtaPhiMVector tau_had0;
    ROOT::Math::PtEtaPhiMVector met;

    double bjet0_pt_scale_fac;
    double bjet1_pt_scale_fac;

    double score;
    double score_m_bb;
    double score_m_tt;
    double score_met;
    double score_chi;

    bool veto;
};



struct event getEvent(ROOT::Math::PtEtaPhiMVector bjet0, ROOT::Math::PtEtaPhiMVector bjet1, ROOT::Math::PtEtaPhiMVector lep0, ROOT::Math::PtEtaPhiMVector tau_had0, ROOT::Math::PtEtaPhiMVector met)
{

    event evt;
    evt.bjet0 = bjet0;
    evt.bjet1 = bjet1;
    evt.lep0  = lep0;
    evt.met   = met;
    evt.tau_had0 = tau_had0;
    evt.veto = 1;
    
    std::vector<double> gpu_pass_sfs = GPUScaleB(  bjet0.Pt(), bjet0.Eta(), bjet0.Phi(), bjet0.M(),
                                                    bjet1.Pt(), bjet1.Eta(), bjet1.Phi(), bjet1.M(),
                                                    lep0.Pt(), lep0.Eta(), lep0.Phi(), lep0.M(),
                                                    tau_had0.Pt(), tau_had0.Eta(), tau_had0.Phi(), tau_had0.M(),
                                                    met.Pt(), met.Eta(), met.Phi(), met.M());
    if(gpu_pass_sfs.size() > 0)
    {
        evt.veto=0;
        evt.bjet0_pt_scale_fac = gpu_pass_sfs[0];
        evt.bjet1_pt_scale_fac = gpu_pass_sfs[1];
        
        evt.score              = gpu_pass_sfs[2];

        evt.score_m_bb         = gpu_pass_sfs[3];
        evt.score_m_tt         = gpu_pass_sfs[4];
        evt.score_met          = gpu_pass_sfs[5];
        evt.score_chi          = gpu_pass_sfs[6];
    }
    
    //CPU INPLIMENTATION
    // double min_m_tautau_diff_scaled = 1000;
    // ROOT::Math::PtEtaPhiMVector bjet0_scaled;
    // ROOT::Math::PtEtaPhiMVector bjet1_scaled;
    // ROOT::Math::PxPyPzEVector met_scaled;
    // for (double i = 0.5; i < 2.0; i += 0.01)
    // {
    //     for (double j = 0.5; j < 2.0; j+=0.01)
    //     {
    //         bjet0_scaled = bjet0;
    //         bjet1_scaled = bjet1;
    //         met_scaled   = met;
    //         bjet0_scaled.SetPt(i*bjet0.Pt());
    //         bjet1_scaled.SetPt(j*bjet1.Pt());
    //         //std::cout<<"Pt_ori\t"<<met.Pt()<<"Phi_ori\t"<<met.Phi()<<"\n";
    //         met_scaled.SetPx(met.Px() - (bjet0_scaled.Px() - bjet0.Px()) - (bjet1_scaled.Px() - bjet1.Px()));
    //         met_scaled.SetPy(met.Py() - (bjet0_scaled.Py() - bjet0.Py()) - (bjet1_scaled.Py() - bjet1.Py()));
    //         //std::cout<<"Pt_aft\t"<<met_scaled.Pt()<<"Phi_aft\t"<<met_scaled.Phi()<<"\n\n";
    //         double m_bb_scaled = (bjet0_scaled + bjet1_scaled).M();
    //         if(std::abs(m_bb_scaled - Z_MASS) < 1.0)
    //         {
    //             double dphi_hl = getdphi(tau_had0.phi(), lep0.phi());
    //             double dphi_hv_scaled = getdphi(tau_had0.phi(), met_scaled.phi());
    //             double dphi_lv_scaled = getdphi(lep0.phi(),     met_scaled.phi());
    //             bool inside_hl_scaled = (dphi_hl * dphi_hv_scaled > 0) && (std::abs(dphi_hl) > std::abs(dphi_hv_scaled));
    //             bool close_to_h_scaled = std::abs(dphi_hv_scaled) < 0.17453292;// 10 degree
    //             bool close_to_l_scaled = std::abs(dphi_lv_scaled) < 0.17453292;
    //             bool v_pos_pass_scaled = inside_hl_scaled || close_to_h_scaled || close_to_l_scaled;
    //             if(v_pos_pass_scaled)
    //             {
    //                 ROOT::Math::PtEtaPhiMVector vh_scaled(0,0,0,0);
    //                 ROOT::Math::PtEtaPhiMVector vl_scaled(0,0,0,0);
    //                 if (!inside_hl_scaled && close_to_h_scaled) 
    //                 {
    //                     vh_scaled = {met_scaled.Pt() * std::cos(std::abs(dphi_hv_scaled)), tau_had0.Eta(), tau_had0.Phi(), 0};
    //                 }
    //                 else if (!inside_hl_scaled && close_to_l_scaled)
    //                 {
    //                     vl_scaled = {met_scaled.Pt() * std::cos(std::abs(dphi_lv_scaled)), lep0.Eta(), lep0.Phi(), 0};
    //                 }
    //                 else if (inside_hl_scaled)
    //                 {
    //                     double vhpt_scaled = met_scaled.Pt() * std::cos(std::abs(dphi_hv_scaled)) - met_scaled.Pt() * std::sin(std::abs(dphi_hv_scaled)) * (1/std::tan(std::abs(dphi_hl)));
    //                     vh_scaled = {vhpt_scaled, tau_had0.Eta(), tau_had0.Phi(), 0};
    //                     double vlpt_scaled = met_scaled.Pt() * std::sin(std::abs(dphi_hv_scaled)) / std::sin(std::abs(dphi_hl));
    //                     vl_scaled = {vlpt_scaled, lep0.Eta(), lep0.Phi(), 0};
    //                 }
    //                 double m_tautau_scaled = (vh_scaled+vl_scaled+tau_had0+lep0).M();
    //                 double m_tautau_diff_scaled = std::abs(m_tautau_scaled - Z_MASS);

    //                 if (m_tautau_diff_scaled < 3)
    //                 {
    //                     if(min_m_tautau_diff_scaled > m_tautau_diff_scaled)
    //                     {
    //                         min_m_tautau_diff_scaled = m_tautau_diff_scaled;
    //                         evt.bjet0_pt_scale_fac = i;
    //                         evt.bjet1_pt_scale_fac = j;
    //                         evt.bjet0_scaled = bjet0_scaled;
    //                         evt.bjet1_scaled = bjet1_scaled;
    //                         evt.met_scaled   = met_scaled;
    //                         evt.vh_scaled    = vh_scaled;
    //                         evt.vl_scaled    = vl_scaled;
    //                         evt.m_bb_scaled  = m_bb_scaled;
    //                         evt.m_tt_scaled  = m_tautau_scaled;
    //                         evt.veto = false;
    //                     }
    //                 }
    //             }
    //         }
    //     }
    // }
    return evt;
}



