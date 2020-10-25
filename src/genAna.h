#include <TROOT.h>
#include <TH1.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TMath.h>
#include <ROOT/RDataFrame.hxx>
#include <Math/GenVector/LorentzVector.h>
#include <Math/GenVector/PtEtaPhiM4D.h>
#include <Math/GenVector/PtEtaPhiM4Dfwd.h>
#include <Math/Vector4Dfwd.h>
#include <Math/Vector3Dfwd.h>
#include <Math/GenVector/Boost.h>
#include <TChain.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TLatex.h>

#include "TMatrixD.h"
#include "TVectorD.h"
#include "TMath.h"
#include "Math/Point3D.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/GenVector/Rotation3D.h"
#include "Math/GenVector/EulerAngles.h"
#include "Math/GenVector/AxisAngle.h"
#include "Math/GenVector/Quaternion.h"
#include "Math/GenVector/LorentzRotation.h"
#include "Math/GenVector/Boost.h"
#include "Math/GenVector/Transform3D.h"
#include "Math/GenVector/Plane3D.h"
#include "Math/GenVector/VectorUtil.h"
#include <Math/QuantFuncMathMore.h>
#include <TLorentzVector.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#define Z_MASS 91.1876e3
#define H_MASS 125.18e3
#define GeV 1e3

enum lepType {MUON, ELEC, NO_LEP, MULTI_LEP};

enum vPos {INSIDE_HL_CLOSE_TO_H, INSIDE_HL_CLOSE_TO_L, OUTSIDE_HL_CLOSE_TO_H, OUTSIDE_HL_CLOSE_TO_L, FARAWAY_FROM_HL};


double mTauTauVis(ROOT::Math::PtEtaPhiMVector tau_had, ROOT::Math::PtEtaPhiMVector ele, ROOT::Math::PtEtaPhiMVector muon, uint is_ele, uint is_mu){
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

double getdphi(double v1phi,  double v2phi)
{
    double ret;
    if(std::abs(v1phi - v2phi) > M_PI)
        ret = ((double)(v1phi > 0) -  (double)(v1phi <= 0)) * (2 * M_PI - std::abs(v1phi - v2phi));
    else
        ret = v2phi - v1phi;
    return ret;
}

lepType getLepType(uint is_ele, uint is_mu)
{
    lepType ret;
    if (is_mu && is_ele ) ret = MULTI_LEP;
    else if(is_mu) ret = MUON;
    else if(is_ele) ret = ELEC;
    else ret = NO_LEP;
    return ret;
}

ROOT::Math::PtEtaPhiMVector getLepP4(lepType type, ROOT::Math::PtEtaPhiMVector elec, ROOT::Math::PtEtaPhiMVector muon)
{   
    ROOT::Math::PtEtaPhiMVector lep(0,0,0,0);
    if(type == MUON) lep = muon;
    else if (type == ELEC) lep = elec;
    
    return lep;
}

vPos getVPos(ROOT::Math::PtEtaPhiMVector met, ROOT::Math::PtEtaPhiMVector tau_had, 
                ROOT::Math::PtEtaPhiMVector lep, lepType lep_type)
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

ROOT::Math::PtEtaPhiMVector getVHP4(vPos v_pos, ROOT::Math::PtEtaPhiMVector met, ROOT::Math::PtEtaPhiMVector tau_had, ROOT::Math::PtEtaPhiMVector lep)
{
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

ROOT::Math::PtEtaPhiMVector getVLP4(vPos v_pos, ROOT::Math::PtEtaPhiMVector met, ROOT::Math::PtEtaPhiMVector tau_had, ROOT::Math::PtEtaPhiMVector lep)
{
    ROOT::Math::PtEtaPhiMVector vl(0,0,0,0);
    double dphi_lv = getdphi(lep.phi(),     met.phi());
    double dphi_hl = getdphi(tau_had.phi(), lep.phi());
    double dphi_hv = getdphi(tau_had.phi(), met.phi());
    if (v_pos == OUTSIDE_HL_CLOSE_TO_L) 
    {
        vl = {met.Pt() * std::cos(std::abs(dphi_lv)), lep.Eta(), tau_had.Phi(), 0};
    }
    else if (v_pos == INSIDE_HL_CLOSE_TO_H || v_pos == INSIDE_HL_CLOSE_TO_L)
    {
        double vlpt = met.Pt() * std::sin(std::abs(dphi_hv)) / std::sin(std::abs(dphi_hl));
        vl = {vlpt, lep.Eta(), lep.Phi(), 0};
    }
    return vl;
}

double mTauTauColin(ROOT::Math::PtEtaPhiMVector met, ROOT::Math::PtEtaPhiMVector tau_had, ROOT::Math::PtEtaPhiMVector ele, ROOT::Math::PtEtaPhiMVector muon, uint is_ele, uint is_mu){
    double ret;
    ROOT::Math::PtEtaPhiMVector lep;
    if (is_mu && is_ele) ret = -999;
    else if(is_mu) lep = muon;
    else if(is_ele) lep = ele;

    double dphi_hl = getdphi(tau_had.phi(), lep.phi());
    double dphi_hv = getdphi(tau_had.phi(), met.phi());
    double dphi_lv = getdphi(lep.phi(),     met.phi());
    bool inside_hl = (dphi_hl * dphi_hv > 0) && (std::abs(dphi_hl) > std::abs(dphi_hv));
    bool close_to_h = std::abs(dphi_hv) < std::abs(dphi_lv);
    bool veto = (!inside_hl) && (std::abs(dphi_hv) > M_PI_2 || std::abs(dphi_lv) > M_PI_2);
    if (veto) ret = -999;
    else if(!inside_hl && close_to_h)
    {
        ROOT::Math::PtEtaPhiMVector vh(met.Pt() * std::cos(std::abs(dphi_hv)), tau_had.Eta(), tau_had.Phi(), 0);
        ret = (vh + lep + tau_had).M();
    }
    else if(!inside_hl && !close_to_h)
    {
        ROOT::Math::PtEtaPhiMVector vl(met.Pt() * std::cos(std::abs(dphi_lv)), lep.Eta(), tau_had.Phi(), 0);
        ret = (vl + lep + tau_had).M();
    }
    else
    {
        double vhpt = met.Pt() * std::cos(std::abs(dphi_hv)) - met.Pt() * std::sin(std::abs(dphi_hv)) * (1/std::tan(std::abs(dphi_hl)));
        double vlpt = met.Pt() * std::sin(std::abs(dphi_hv)) / std::sin(std::abs(dphi_hl));
        ROOT::Math::PtEtaPhiMVector vl(vlpt, lep.Eta(), lep.Phi(), 0);
        ROOT::Math::PtEtaPhiMVector vh(vhpt, tau_had.Eta(), tau_had.Phi(), 0);
        ret = (vh + vl + lep + tau_had).M();
    }
    return ret;
}