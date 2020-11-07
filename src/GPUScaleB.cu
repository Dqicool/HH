#include <cstdio>
#include <vector>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <math_constants.h>
#include "GPUScaleB.cuh"
#define Z_MASS 91.1876

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

__global__ void mykernel(   float bjet0_pt, float bjet0_eta,    float bjet0_phi,    float bjet0_m,
                            float bjet1_pt, float bjet1_eta,    float bjet1_phi,    float bjet1_m,
                            float lep0_pt,  float lep0_eta,     float lep0_phi,     float lep0_m,
                            float tau0_pt,  float tau0_eta,     float tau0_phi,     float tau0_m,
                            float met_pt,   float met_eta,      float met_phi,      float met_m,
                            bool * pass, float* m_tt, int N)
{
    int i = (threadIdx.x + blockIdx.x * blockDim.x);
    
    int j = (threadIdx.y + blockIdx.y * blockDim.y);
    float chi0 = 0.5 + 0.01 * i;
    
    float chi1 = 0.5 + 0.01 * j;

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

            if (m_tautau_diff_scaled < 3.0)
            {
                //printf("i:%d\tj:%d\tchi0:%f\tchi1:%f\tPASS3\n",i,j,chi0,chi1);
                pass[i*N + j] = 1;
                m_tt[i*N + j] = m_tautau_diff_scaled;
            }
        }
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
    bool *pass;
    pass = (bool*)malloc(size_bool);
    float *m_tt;
    m_tt = (float*)malloc(size_float);

    bool *pass_dev;
    cudaMalloc((void**)&pass_dev, size_bool);

    float *m_tt_dev;
    cudaMalloc((void**)&m_tt_dev, size_float);

    for(int i =0; i < N*N; i++)
    {
        pass[i]=0;
        m_tt[i]=0;
    }
    cudaMemcpy(pass_dev,pass,size_bool,cudaMemcpyHostToDevice);
    cudaMemcpy(m_tt_dev,m_tt,size_float,cudaMemcpyHostToDevice);

    dim3 grid(N, N);
    mykernel<<<grid, 1>>>(  bjet0_pt,  bjet0_eta,     bjet0_phi,     bjet0_m,
                            bjet1_pt,  bjet1_eta,     bjet1_phi,     bjet1_m,
                            lep0_pt,   lep0_eta,      lep0_phi,      lep0_m,
                            tau0_pt,   tau0_eta,      tau0_phi,      tau0_m,
                            met_pt,    met_eta,       met_phi,       met_m,
                            pass_dev, m_tt_dev, N);
    //cudaDeviceSynchronize();
    cudaMemcpy(pass,pass_dev,size_bool,cudaMemcpyDeviceToHost);
    cudaMemcpy(m_tt,m_tt_dev,size_float,cudaMemcpyDeviceToHost);
    cudaFree(pass_dev);
    cudaFree(m_tt_dev);
    double min_mtt_diff = 100000;
    std::vector<double> ret;
    double sf1 = -999;
    double sf2 = -999;
    for(int i =0; i < N * N; i++){
        if (pass[i])
        {
            if(m_tt[i] < min_mtt_diff)
            {
                min_mtt_diff = m_tt[i];
                sf1 = 0.01 * (double)(i/151) + 0.5;
                sf2 = 0.01 * (double)(i%151) + 0.5;
            }
        }
    }
    if (sf1 > 0 && sf2 > 0)
    {
        ret.push_back(sf1);
        ret.push_back(sf2);
    }
    free(pass);
    free(m_tt);

    return ret;
}
