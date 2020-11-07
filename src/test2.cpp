#include "genAna.h"
#define CUDART_PI_F 3.1415926f
#define Z_MASS 91.1876

double getdphi(double v1phi,  double v2phi)
{
    double ret;
    if(std::abs(v1phi - v2phi) > M_PI)
        ret = ((double)(v1phi > 0) -  (double)(v1phi <= 0)) * (2 * M_PI - std::abs(v1phi - v2phi));
    else
        ret = v2phi - v1phi;
    return ret;
}

float getPx(float pt,float phi)
{
    return pt * cosf(phi);
}

float getPy(float pt,float phi){
    return pt * sinf(phi);
}

float getPz(float pt,float eta){
    return pt * sinhf(eta);
}

float getPScaler(float pt,float eta){
    return pt * coshf(eta);
}

float getE(float pt,float eta,float m){
   float p = getPScaler(pt, eta);
    return sqrtf(m * m + p * p);
}

float getMass(float energy,float px,float py,float pz){

    return sqrtf(energy * energy - px * px - py * py - pz * pz);
}

float getdphi(float v1phi, float v2phi)
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

bool cpu(ROOT::Math::PtEtaPhiMVector bjet0, 
            ROOT::Math::PtEtaPhiMVector bjet1, 
            ROOT::Math::PtEtaPhiMVector lep0, 
            ROOT::Math::PtEtaPhiMVector tau_had0, 
            ROOT::Math::PtEtaPhiMVector met,
            double chi0, double chi1)
{
    bool pass=0;
    double min_m_tautau_diff_scaled = 1000;
    ROOT::Math::PtEtaPhiMVector bjet0_scaled;
    ROOT::Math::PtEtaPhiMVector bjet1_scaled;
    ROOT::Math::PxPyPzEVector met_scaled;
    double i = chi0;
    double j = chi1;
    bjet0_scaled = bjet0;
    bjet1_scaled = bjet1;
    met_scaled   = met;
    bjet0_scaled.SetPt(i*bjet0.Pt());
    bjet1_scaled.SetPt(j*bjet1.Pt());
    //std::cout<<"Pt_ori\t"<<met.Pt()<<"Phi_ori\t"<<met.Phi()<<"\n";
    met_scaled.SetPx(met.Px() - (bjet0_scaled.Px() - bjet0.Px()) - (bjet1_scaled.Px() - bjet1.Px()));
    met_scaled.SetPy(met.Py() - (bjet0_scaled.Py() - bjet0.Py()) - (bjet1_scaled.Py() - bjet1.Py()));
    double met_scaled_px = met_scaled.Px();
    double met_scaled_py = met_scaled.Py();
    //std::cout<<"Pt_aft\t"<<met_scaled.Pt()<<"Phi_aft\t"<<met_scaled.Phi()<<"\n\n";
    double m_bb_scaled = (bjet0_scaled + bjet1_scaled).M();
    if(std::abs(m_bb_scaled - Z_MASS) < 1.0)
    {
        double dphi_hl = getdphi(tau_had0.phi(), lep0.phi());
        double dphi_hv_scaled = getdphi(tau_had0.phi(), met_scaled.phi());
        double dphi_lv_scaled = getdphi(lep0.phi(),     met_scaled.phi());
        bool inside_hl_scaled = (dphi_hl * dphi_hv_scaled > 0) && (std::abs(dphi_hl) > std::abs(dphi_hv_scaled));
        bool close_to_h_scaled = std::abs(dphi_hv_scaled) < 0.17453292;// 10 degree
        bool close_to_l_scaled = std::abs(dphi_lv_scaled) < 0.17453292;
        bool v_pos_pass_scaled = inside_hl_scaled || close_to_h_scaled || close_to_l_scaled;
        if(v_pos_pass_scaled)
        {
            ROOT::Math::PtEtaPhiMVector vh_scaled(0,0,0,0);
            ROOT::Math::PtEtaPhiMVector vl_scaled(0,0,0,0);
            if (!inside_hl_scaled && close_to_h_scaled) 
            {
                vh_scaled = {met_scaled.Pt() * std::cos(std::abs(dphi_hv_scaled)), tau_had0.Eta(), tau_had0.Phi(), 0};
            }
            else if (!inside_hl_scaled && close_to_l_scaled)
            {
                vl_scaled = {met_scaled.Pt() * std::cos(std::abs(dphi_lv_scaled)), lep0.Eta(), lep0.Phi(), 0};
            }
            else if (inside_hl_scaled)
            {
                double vhpt_scaled = met_scaled.Pt() * std::cos(std::abs(dphi_hv_scaled)) - met_scaled.Pt() * std::sin(std::abs(dphi_hv_scaled)) * (1/std::tan(std::abs(dphi_hl)));
                vh_scaled = {vhpt_scaled, tau_had0.Eta(), tau_had0.Phi(), 0};
                double vlpt_scaled = met_scaled.Pt() * std::sin(std::abs(dphi_hv_scaled)) / std::sin(std::abs(dphi_hl));
                vl_scaled = {vlpt_scaled, lep0.Eta(), lep0.Phi(), 0};
            }
            double m_tautau_scaled = (vh_scaled+vl_scaled+tau_had0+lep0).M();
            double m_tautau_diff_scaled = std::abs(m_tautau_scaled - Z_MASS);

            if (m_tautau_diff_scaled < 3)
            {
                if(min_m_tautau_diff_scaled > m_tautau_diff_scaled)
                {
                    min_m_tautau_diff_scaled = m_tautau_diff_scaled;
                    pass = 1;
                }
            }
        }
    }
    return pass;
}
    

bool gpu(   float bjet0_pt, float bjet0_eta,    float bjet0_phi,    float bjet0_m,
                            float bjet1_pt, float bjet1_eta,    float bjet1_phi,    float bjet1_m,
                            float lep0_pt,  float lep0_eta,     float lep0_phi,     float lep0_m,
                            float tau0_pt,  float tau0_eta,     float tau0_phi,     float tau0_m,
                            float met_pt,   float met_eta,      float met_phi,      float met_m,
                            float chi0, float chi1)
{
    bool pass=0;

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
    
    //printf("chi0:%f\tchi1:%f\tmbb:%f\n",chi0, chi1, m_bb_scaled);
    if(fabsf(m_bb_scaled - Z_MASS) < 1.0)
    {
        //printf("i:%d\tj:%d\tchi0:%f\tchi1:%f\tPASS1\n",i,j,chi0,chi1);
        float dphi_hl = getdphi(tau0_phi, lep0_phi);
        float met_scaled_phi = getPhiFromPxPy(met_scaled_px, met_scaled_py);
        float met_scaled_pt  = sqrtf(met_scaled_py * met_scaled_py + met_scaled_px * met_scaled_px);
        float dphi_hv_scaled = getdphi(tau0_phi, met_scaled_phi);
        float dphi_lv_scaled = getdphi(lep0_phi, met_scaled_phi);
        bool inside_hl_scaled = (dphi_hl * dphi_hv_scaled > 0) && (fabsf(dphi_hl) > fabsf(dphi_hv_scaled));
        bool close_to_h_scaled = fabsf(dphi_hv_scaled) < 0.17453292f;// 10 degree
        bool close_to_l_scaled = fabsf(dphi_lv_scaled) < 0.17453292f;
        bool v_pos_pass_scaled = inside_hl_scaled || close_to_h_scaled || close_to_l_scaled;
        if(v_pos_pass_scaled)
        {
            //printf("i:%d\tj:%d\tchi0:%f\tchi1:%f\tPASS2\n",i,j,chi0,chi1);
            float vh_scaled_pt =0;
            float vh_scaled_eta=0;
            float vh_scaled_phi=0;
            float vh_scaled_m  =0;

            float vl_scaled_pt =0;
            float vl_scaled_eta=0;
            float vl_scaled_phi=0;
            float vl_scaled_m  =0;
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
            float  m_tautau_scaled = getMass(tautau_e, tautau_px, tautau_py, tautau_pz);
            float m_tautau_diff_scaled = fabsf(m_tautau_scaled - Z_MASS);

            if (m_tautau_diff_scaled < 3)
            {
                pass = 1;
            }
        }
    }
    return pass;
}

int main(){
    ROOT::RDataFrame df("NOMINAL", "debug.root");
    auto ef = df.Filter("cpu_pass");
    auto bjet0 = (ef.Take<ROOT::Math::PtEtaPhiMVector>("bjet_0_p4_n").GetValue())[0];
    auto bjet1 = (ef.Take<ROOT::Math::PtEtaPhiMVector>("bjet_1_p4_n").GetValue())[0];
    auto lep0 =  (ef.Take<ROOT::Math::PtEtaPhiMVector>("lep_0_p4_n").GetValue())[0];
    auto tau0 =  (ef.Take<ROOT::Math::PtEtaPhiMVector>("tau_0_p4_n").GetValue())[0];
    auto met =   (ef.Take<ROOT::Math::PtEtaPhiMVector>("met_reco_p4_n").GetValue())[0];
    double chi0 = (ef.Take<double>("bjet_0_pt_scale_fac").GetValue())[0];
    double chi1 = (ef.Take<double>("bjet_1_pt_scale_fac").GetValue())[0];
    bool cpu_res = cpu(bjet0, bjet1, lep0, tau0, met, chi0, chi1);
    bool gpu_res = gpu(bjet0.Pt(), bjet0.Eta(), bjet0.Phi(), bjet0.M(),
                    bjet1.Pt(), bjet1.Eta(), bjet1.Phi(), bjet1.M(),
                    lep0.Pt(), lep0.Eta(), lep0.Phi(), lep0.M(),
                    tau0.Pt(), tau0.Eta(), tau0.Phi(), tau0.M(),
                    met.Pt(), met.Eta(), met.Phi(), met.M(),
                    chi0, chi1);
}