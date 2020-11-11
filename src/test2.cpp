#include "genAna.h"
#define CUDART_PI_F 3.1415926f
#define CUDART_PIO4_F (CUDART_PI_F / 4.0f)
#define Z_MASS 91.1876

#include <cstdio>
#include <vector>

#define Z_MASS 91.1876

#define M_BB_SC     62.5
#define M_TT_SC     62.5
#define MET_SC_1    0.08
#define MET_SC_2    0.005
#define MET_SC_3    0.00125
#define CHI_SC      0.1

 float getPx(float pt, float phi){
    return pt * cosf(phi);
}

 float getPy(float pt, float phi){
    return pt * sinf(phi);
}

 float getPz(float pt, float eta){
    return pt * sinhf(eta);
}

 float getPScaler(float pt, float eta){
    return pt * coshf(eta);
}

 float getE(float pt, float eta, float m){
    float p = getPScaler(pt, eta);
    return sqrtf(m * m + p * p);
}

 float getMass( float energy, float px, float py, float pz){

    return sqrtf(energy * energy - px * px - py * py - pz * pz);
}

 float getdphi(float v1phi,  float v2phi)
{
    float ret;
    if(fabsf(v1phi - v2phi) > CUDART_PI_F)
        ret = ((float)(v1phi > 0) -  (float)(v1phi <= 0)) * (2 * CUDART_PI_F - fabsf(v1phi - v2phi));
    else
        ret = v2phi - v1phi;
    return ret;
}

 float getPhiFromPxPy(float px, float py)
{
    return  (float)(px>0) * atanf(py / px) + 
            (float)(px<=0 && py > 0) * (atanf(py / px) + CUDART_PI_F) + 
            (float)(px<=0 && py <= 0) * (atanf(py / px) - CUDART_PI_F);
}

 float getPsiMbb(float m_bb, float m_baseline, float scale_const_bb)
{
    return (m_bb - m_baseline) * (m_bb - m_baseline) / scale_const_bb;
}

 float getPsiMtt(float m_tt, float m_baseline, float scale_const_tt)
{
    return (m_tt - m_baseline) * (m_tt - m_baseline) / scale_const_tt;
}

 float getPsiMET(float vl_pt, float vh_pt, float tau_pt, float lep_pt, float omega,float sc_1, float sc_2, float sc_3)
{
    float f_tau_l = vl_pt / (vl_pt + lep_pt);
    float f_tau_h = vh_pt / (vh_pt + tau_pt);
    float A_1 = 0.5 * ((f_tau_l - f_tau_h) / (f_tau_l + f_tau_h) + 1);
    return (float)(A_1 > 0 && A_1 <= 0.8) * (A_1 - 0.8) * (A_1 - 0.8) / sc_1 + 
           (float)(A_1 > 0.8 && A_1 <= 1) * 0 + 
           (float)(A_1 > 1)               * (omega - 1) * (omega - 1) / sc_2 + 
           (float)(A_1 <= 0)               * (8 + omega * omega / sc_3);
}

 float getPsiChi(float chi0, float chi1, float sc_chi)
{
    float log2_chi0 = log2f(chi0);
    float log2_chi1 = log2f(chi1);
    float r2 = log2_chi0 * log2_chi0 + log2_chi1 * log2_chi1;
    return r2 / sc_chi;
}

std::vector<double> gpu(   float bjet0_pt, float bjet0_eta,    float bjet0_phi,    float bjet0_m,
                            float bjet1_pt, float bjet1_eta,    float bjet1_phi,    float bjet1_m,
                            float lep0_pt,  float lep0_eta,     float lep0_phi,     float lep0_m,
                            float tau0_pt,  float tau0_eta,     float tau0_phi,     float tau0_m,
                            float met_pt,   float met_eta,      float met_phi,      float met_m,
                            float chi0, float chi1)
{
    std::vector<double> ret;
    //mbb calculation
    float bjet0_scaled_pt = bjet0_pt * chi0;
    float bjet1_scaled_pt = bjet1_pt * chi1;

    float met_scaled_px = getPx(met_pt, met_phi) - (getPx(bjet0_scaled_pt, bjet0_phi) - getPx(bjet0_pt, bjet0_phi)) - (getPx(bjet1_scaled_pt, bjet1_phi) - getPx(bjet1_pt, bjet1_phi));
    float met_scaled_py = getPy(met_pt, met_phi) - (getPy(bjet0_scaled_pt, bjet0_phi) - getPy(bjet0_pt, bjet0_phi)) - (getPy(bjet1_scaled_pt, bjet1_phi) - getPy(bjet1_pt, bjet1_phi));
    
    float bjet0_scaled_E = getE(bjet0_scaled_pt, bjet0_eta, bjet0_m);
    float bjet0_scaled_px = getPx(bjet0_scaled_pt, bjet0_phi);
    float bjet0_scaled_py = getPy(bjet0_scaled_pt, bjet0_phi);
    float bjet0_scaled_pz = getPz(bjet0_scaled_pt, bjet0_eta);

    float bjet1_scaled_E =  getE(bjet1_scaled_pt, bjet1_eta, bjet1_m);
    float bjet1_scaled_px = getPx(bjet1_scaled_pt, bjet1_phi);
    float bjet1_scaled_py = getPy(bjet1_scaled_pt, bjet1_phi);
    float bjet1_scaled_pz = getPz(bjet1_scaled_pt, bjet1_eta);


    float m_bb_scaled = getMass(bjet0_scaled_E + bjet1_scaled_E, 
                                bjet0_scaled_px + bjet1_scaled_px,
                                bjet0_scaled_py + bjet1_scaled_py,
                                bjet0_scaled_pz + bjet1_scaled_pz);   
    
    //omega and other
    float dphi_hl = getdphi(tau0_phi, lep0_phi);
    float met_scaled_phi = getPhiFromPxPy(met_scaled_px, met_scaled_py);
    float met_scaled_pt  = sqrtf(met_scaled_py * met_scaled_py + met_scaled_px * met_scaled_px);
    float dphi_hv_scaled = getdphi(tau0_phi, met_scaled_phi);
    float dphi_lv_scaled = getdphi(lep0_phi, met_scaled_phi);
    bool inside_hl_scaled = (dphi_hl * dphi_hv_scaled > 0) && (fabsf(dphi_hl) > fabsf(dphi_hv_scaled));
    bool close_to_h_scaled = fabsf(dphi_hv_scaled) < CUDART_PIO4_F;// 10 degree
    bool close_to_l_scaled = fabsf(dphi_lv_scaled) < CUDART_PIO4_F;
    bool v_pos_pass_scaled = inside_hl_scaled || close_to_h_scaled || close_to_l_scaled;
    float omega = dphi_hv_scaled / dphi_hl;

    //p4 of vl vh
    float vh_scaled_pt =0;
    float vh_scaled_eta=0;
    float vh_scaled_phi=0;
    float vh_scaled_m  =0;

    float vl_scaled_pt =0;
    float vl_scaled_eta=0;
    float vl_scaled_phi=0;
    float vl_scaled_m  =0;

    float  m_tautau_scaled = -999;
    if(v_pos_pass_scaled)
    {
        if (!inside_hl_scaled && close_to_h_scaled) 
        {
            vh_scaled_pt  = met_scaled_pt * cosf(fabsf(dphi_hv_scaled));
            vh_scaled_eta = tau0_eta;
            vh_scaled_phi = tau0_phi; 
            vh_scaled_m   = 0;
        }
        else if (!inside_hl_scaled && close_to_l_scaled)
        {
            vl_scaled_pt = met_scaled_pt * cosf(fabsf(dphi_lv_scaled));
            vl_scaled_eta = lep0_eta;
            vl_scaled_phi = lep0_phi;
            vl_scaled_m   =  0;
        }
        else if (inside_hl_scaled)
        {
            vh_scaled_pt    = met_scaled_pt * cosf(fabsf(dphi_hv_scaled)) - met_scaled_pt * sinf(fabsf(dphi_hv_scaled)) * (1/tanf(fabsf(dphi_hl)));
            vh_scaled_eta   = tau0_eta;
            vh_scaled_phi   = tau0_phi; 
            vh_scaled_m     = 0;
            vl_scaled_pt    = met_scaled_pt * sinf(fabsf(dphi_hv_scaled)) / sinf(fabsf(dphi_hl));
            vl_scaled_eta   = lep0_eta;
            vl_scaled_phi   = lep0_phi;
            vl_scaled_m     = 0;
        }
        float  tautau_px = (getPx(vh_scaled_pt, vh_scaled_phi)+getPx(vl_scaled_pt, vl_scaled_phi)+getPx(tau0_pt, tau0_phi)+getPx(lep0_pt, lep0_phi));
        float  tautau_py = (getPy(vh_scaled_pt, vh_scaled_phi)+getPy(vl_scaled_pt, vl_scaled_phi)+getPy(tau0_pt, tau0_phi)+getPy(lep0_pt, lep0_phi));
        float  tautau_pz = (getPz(vh_scaled_pt, vh_scaled_eta)+getPz(vl_scaled_pt, vl_scaled_eta)+getPz(tau0_pt, tau0_eta)+getPz(lep0_pt, lep0_eta));
        float  tautau_e  = (getE(vh_scaled_pt, vh_scaled_eta, vh_scaled_m)+getE(vl_scaled_pt, vl_scaled_eta, vl_scaled_m)+getE(tau0_pt, tau0_eta, tau0_m)+getE(lep0_pt, lep0_eta, lep0_m));
        m_tautau_scaled = getMass(tautau_e, tautau_px, tautau_py, tautau_pz);
    }
    if(m_tautau_scaled > 0)
    {
        float score_m_bb = getPsiMbb(m_bb_scaled, Z_MASS, M_BB_SC);
        float score_m_tt = getPsiMtt(m_tautau_scaled, Z_MASS, M_TT_SC);
        float score_met  = getPsiMET(vl_scaled_pt, vh_scaled_pt, tau0_pt, lep0_pt, omega, MET_SC_1, MET_SC_2, MET_SC_3);
        float score_chi  = getPsiChi(chi0, chi1, CHI_SC);
        ret.push_back(score_m_bb);
        ret.push_back(score_m_tt);
        ret.push_back(score_met);
        ret.push_back(score_chi);
        ret.push_back(sqrtf(score_chi * score_chi + score_m_bb * score_m_bb + score_m_tt * score_m_tt + score_met * score_met));
    }
    return ret;
}

int main(){
    ROOT::RDataFrame df("NOMINAL", "debug.root");
    auto bjet0 = (df.Take<ROOT::Math::PtEtaPhiMVector>("bjet_0_p4_n").GetValue())[0];
    auto bjet1 = (df.Take<ROOT::Math::PtEtaPhiMVector>("bjet_1_p4_n").GetValue())[0];
    auto lep0 =  (df.Take<ROOT::Math::PtEtaPhiMVector>("lep_0_p4_n").GetValue())[0];
    auto tau0 =  (df.Take<ROOT::Math::PtEtaPhiMVector>("tau_0_p4_n").GetValue())[0];
    auto met =   (df.Take<ROOT::Math::PtEtaPhiMVector>("met_reco_p4_n").GetValue())[0];
    double chi0 = (df.Take<double>("bjet_0_pt_scale_fac").GetValue())[0];
    double chi1 = (df.Take<double>("bjet_1_pt_scale_fac").GetValue())[0];
    std::vector<double> gpu_res = gpu(bjet0.Pt(), bjet0.Eta(), bjet0.Phi(), bjet0.M(),
                    bjet1.Pt(), bjet1.Eta(), bjet1.Phi(), bjet1.M(),
                    lep0.Pt(), lep0.Eta(), lep0.Phi(), lep0.M(),
                    tau0.Pt(), tau0.Eta(), tau0.Phi(), tau0.M(),
                    met.Pt(), met.Eta(), met.Phi(), met.M(),
                    chi0, chi1);
}