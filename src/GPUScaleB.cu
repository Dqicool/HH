#include <cstdio>
#include <vector>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <math_constants.h>
#include "GPUScaleB.cuh"
#define Z_MASS 91.1876

#define M_BB_SC     62.5
#define M_TT_SC     62.5
#define MET_SC_1    0.008
#define MET_SC_2    0.0005
#define MET_SC_3    0.000125
#define CHI_SC      0.04

__device__ float getPx(float pt, float phi){
    return pt * cosf(phi);
}

__device__ float getPy(float pt, float phi){
    return pt * sinf(phi);
}

__device__ float getPz(float pt, float eta){
    return pt * sinhf(eta);
}

__device__ float getPScaler(float pt, float eta){
    return pt * coshf(eta);
}

__device__ float getE(float pt, float eta, float m){
    float p = getPScaler(pt, eta);
    return sqrtf(m * m + p * p);
}

__device__ float getMass( float energy, float px, float py, float pz){

    return sqrtf(energy * energy - px * px - py * py - pz * pz);
}

__device__ float getdphi(float v1phi,  float v2phi)
{
    float ret;
    if(fabsf(v1phi - v2phi) > CUDART_PI_F)
        ret = ((float)(v1phi > 0) -  (float)(v1phi <= 0)) * (2 * CUDART_PI_F - fabsf(v1phi - v2phi));
    else
        ret = v2phi - v1phi;
    return ret;
}

__device__ float getPhiFromPxPy(float px, float py)
{
    return  (float)(px>0) * atanf(py / px) + 
            (float)(px<=0 && py > 0) * (atanf(py / px) + CUDART_PI_F) + 
            (float)(px<=0 && py <= 0) * (atanf(py / px) - CUDART_PI_F);
}

__device__ float getPsiMbb(float m_bb, float m_baseline, float scale_const_bb)
{
    return (m_bb - m_baseline) * (m_bb - m_baseline) / scale_const_bb;
}

__device__ float getPsiMtt(float m_tt, float m_baseline, float scale_const_tt)
{
    return (m_tt - m_baseline) * (m_tt - m_baseline) / scale_const_tt;
}

__device__ float getPsiMET(float vl_pt, float vh_pt, float tau_pt, float lep_pt, float omega,float sc_1, float sc_2, float sc_3)
{
    float f_tau_l = vl_pt / (vl_pt + lep_pt);
    float f_tau_h = vh_pt / (vh_pt + tau_pt);
    float A_1 = 0.5 * ((f_tau_l - f_tau_h) / (f_tau_l + f_tau_h) + 1);
    float ret = ((float)(A_1 > 0 && A_1 <= 0.8) * (A_1 - 0.8) * (A_1 - 0.8) / sc_1) +
                ((float)(A_1 > 0.8 && A_1 <= 1) * 0) +
                ((float)(omega > 1)               * (omega - 1) * (omega - 1) / sc_2) +
                ((float)(omega <= 0)              * (80 + (omega * omega / sc_3)));
    return ret;
}

__device__ float getPsiChi(float chi0, float chi1, float sc_chi)
{
    float log2_chi0 = log2f(chi0);
    float log2_chi1 = log2f(chi1);
    float r2 = log2_chi0 * log2_chi0 + log2_chi1 * log2_chi1;
    return r2 / sc_chi;
}

__global__ void mykernel(   float bjet0_pt, float bjet0_eta,    float bjet0_phi,    float bjet0_m,
                            float bjet1_pt, float bjet1_eta,    float bjet1_phi,    float bjet1_m,
                            float lep0_pt,  float lep0_eta,     float lep0_phi,     float lep0_m,
                            float tau0_pt,  float tau0_eta,     float tau0_phi,     float tau0_m,
                            float met_pt,   float met_eta,      float met_phi,      float met_m,
                            bool * pass, float* score, float* score1, float* score2, float* score3, float* score4, int N)
{
    int i = (threadIdx.x + blockIdx.x * blockDim.x);
    int j = (threadIdx.y + blockIdx.y * blockDim.y);
    float chi0 = 0.5 + 0.01 * i;
    float chi1 = 0.5 + 0.01 * j;
    // initialise
    pass[i*N + j] = 0;
    score1[i*N + j]  = -999;
    score2[i*N + j]  = -999;
    score3[i*N + j]  = -999;
    score4[i*N + j]  = -999;
    score[i*N + j]   = -999;
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
    bool close_to_h_scaled = fabsf(dphi_hv_scaled) < CUDART_PIO4_F && fabsf(dphi_hv_scaled) < fabsf(dphi_lv_scaled);
    bool close_to_l_scaled = fabsf(dphi_lv_scaled) < CUDART_PIO4_F && fabsf(dphi_lv_scaled) < fabsf(dphi_hv_scaled);
    bool v_pos_pass_scaled = inside_hl_scaled || close_to_h_scaled || close_to_l_scaled;
    float omega = -999;
    if(!inside_hl_scaled && close_to_l_scaled && dphi_hv_scaled * dphi_hl < 0)
    {
        if(dphi_hl < 0)
            omega = (dphi_hv_scaled - 2 * CUDART_PI_F) / dphi_hl;
        else if (dphi_hv_scaled < 0)
            omega = (dphi_hv_scaled + CUDART_PI_F * 2.0f) / dphi_hl;
    }
    else
    {
        omega = dphi_hv_scaled / dphi_hl;
    }

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
        pass[i*N + j] = 1;
        float score_m_bb = getPsiMbb(m_bb_scaled, Z_MASS, M_BB_SC);
        float score_m_tt = getPsiMtt(m_tautau_scaled, Z_MASS, M_TT_SC);
        float score_met  = getPsiMET(vl_scaled_pt, vh_scaled_pt, tau0_pt, lep0_pt, omega, MET_SC_1, MET_SC_2, MET_SC_3);
        float score_chi  = getPsiChi(chi0, chi1, CHI_SC);
        score1[i*N + j]  = score_m_bb;
        score2[i*N + j]  = score_m_tt;
        score3[i*N + j]  = score_met;
        score4[i*N + j]  = score_chi;
        score[i*N + j] = sqrtf(score_chi * score_chi + score_m_bb * score_m_bb + score_m_tt * score_m_tt + score_met * score_met);
    }
}

std::vector<double> GPUScaleB(float bjet0_pt, float bjet0_eta,    float bjet0_phi,    float bjet0_m,
    float bjet1_pt, float bjet1_eta,    float bjet1_phi,    float bjet1_m,
    float lep0_pt,  float lep0_eta,     float lep0_phi,     float lep0_m,
    float tau0_pt,  float tau0_eta,     float tau0_phi,     float tau0_m,
    float met_pt,   float met_eta,      float met_phi,      float met_m){
    int N = 151;
    size_t size_bool  = N * N * sizeof(bool);
    size_t size_float = N * N * sizeof(float);
    //allocate arrays on CPU
    bool pass[N*N];
    
    float score[N*N];
    float score1[N*N];
    float score2[N*N];
    float score3[N*N];
    float score4[N*N];
    //allocate arrays on GPU
    bool *pass_dev;
    cudaMalloc((void**)&pass_dev, size_bool);

    float *score_dev, *score1_dev, *score2_dev, *score3_dev, *score4_dev; 
    cudaMalloc((void**)&score_dev, size_float);
    cudaMalloc((void**)&score1_dev, size_float);
    cudaMalloc((void**)&score2_dev, size_float);
    cudaMalloc((void**)&score3_dev, size_float);
    cudaMalloc((void**)&score4_dev, size_float);
    //copy stuff to GPU
    // cudaMemcpy(pass_dev,pass,size_bool,cudaMemcpyHostToDevice);
    // cudaMemcpy(score_dev,score,size_float,cudaMemcpyHostToDevice);
    // cudaMemcpy(score1_dev,score1,size_float,cudaMemcpyHostToDevice);
    // cudaMemcpy(score2_dev,score2,size_float,cudaMemcpyHostToDevice);
    // cudaMemcpy(score3_dev,score3,size_float,cudaMemcpyHostToDevice);
    // cudaMemcpy(score4_dev,score4,size_float,cudaMemcpyHostToDevice);
    ////memset version
    // cudaMemset(score_dev, 0, size_float);
    // cudaMemset(score1_dev, 0, size_float);
    // cudaMemset(score2_dev, 0, size_float);
    // cudaMemset(score3_dev, 0, size_float);
    // cudaMemset(score4_dev, 0, size_float);
    // cudaMemset(pass_dev, 0, size_float);

    dim3 grid(N, N);
    mykernel<<<grid, 1>>>(  bjet0_pt,  bjet0_eta,     bjet0_phi,     bjet0_m,
                            bjet1_pt,  bjet1_eta,     bjet1_phi,     bjet1_m,
                            lep0_pt,   lep0_eta,      lep0_phi,      lep0_m,
                            tau0_pt,   tau0_eta,      tau0_phi,      tau0_m,
                            met_pt,    met_eta,       met_phi,       met_m,
                            pass_dev, score_dev, score1_dev, score2_dev, score3_dev, score4_dev, N);
    //cudaDeviceSynchronize();
    cudaMemcpy(pass,pass_dev,size_bool,cudaMemcpyDeviceToHost);
    cudaMemcpy(score,score_dev,size_float,cudaMemcpyDeviceToHost);
    cudaMemcpy(score1,score1_dev,size_float,cudaMemcpyDeviceToHost);
    cudaMemcpy(score2,score2_dev,size_float,cudaMemcpyDeviceToHost);
    cudaMemcpy(score3,score3_dev,size_float,cudaMemcpyDeviceToHost);
    cudaMemcpy(score4,score4_dev,size_float,cudaMemcpyDeviceToHost);
    cudaFree(pass_dev);
    cudaFree(score_dev);
    cudaFree(score1_dev);
    cudaFree(score2_dev);
    cudaFree(score3_dev);
    cudaFree(score4_dev);

    std::vector<double> ret;
    double sf1 = -999;
    double sf2 = -999;
    double min_score = 200;
    double s1 = 100;
    double s2 = 100;
    double s3 = 100;
    double s4 = 100;
    for(int i =0; i < N * N; i++){
        if (pass[i])
        {
            if(score[i] < min_score)
            {
                min_score = score[i];

                s1 = score1[i];
                s2 = score2[i];
                s3 = score3[i];
                s4 = score4[i];

                sf1 = 0.01 * (double)(i/151) + 0.5;
                sf2 = 0.01 * (double)(i%151) + 0.5;
            }
        }
    }
    if (sf1 > 0 && sf2 > 0)
    {
        ret.push_back(sf1);
        ret.push_back(sf2);
        ret.push_back(min_score);
        ret.push_back(s1);
        ret.push_back(s2);
        ret.push_back(s3);
        ret.push_back(s4);
    }

    return ret;
}
