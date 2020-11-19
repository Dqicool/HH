#include "genAna.h"
#define CUDART_PI_F 3.1415926f
#define CUDART_PIO4_F (CUDART_PI_F / 4.0f)
#define Z_MASS 91.1876

#include <cstdio>
#include <vector>


float getPsiChi(float chi0, float chi1, float sc_chi)
{
    // calculate on log plane
    float x = log2f(chi0);
    float y = log2f(chi1);

    float x_0      = 9.62213e-02;
    float y_0      = 1.18202e-01;

    float sigma_x  = 4.81883e-01;
    float sigma_y  = 6.37269e-01;

    // float sigma_x  = 4.81883e-01;
    // float sigma_y  = 4.81883e-01;

    float theta    = 5.40952e-01;


    // // calculate on linear plane
    // float x = chi0;
    // float y = chi1;

    // float x_0      = 1.07186e+00;
    // float y_0      = 1.05827e+00;

    // float sigma_x  = 3.27060e-01 / 2;
    // float sigma_y  = 4.99482e-01 / 2;

    // float theta    = 3.50592e-01;

    float ell1 = (powf((x - x_0) * cosf(theta) + (y - y_0) * sinf(theta), 2)) / powf(sigma_x, 2);
    float ell2 = (powf((x - x_0) * sinf(theta) - (y - y_0) * cosf(theta), 2)) / powf(sigma_y, 2);
    float r2 = ell1+ell2;
    return r2 / sc_chi;
}

int main(){
    TH2D* h = new TH2D("h", "", 1501, 0.5, 2, 1501, 0.5, 2);
    for(double x=0.5;x<=2;x+=0.001){
        for(double y=0.5;y<=2;y+=0.001){
            h->Fill(x, y, getPsiChi(x, y, 1));
        }
    }
    TCanvas c("c","",2000,2000);
    h->SetStats(0);
    h->SetTitle("score of tilted ellipse on log(#chi), teanslated to #chi plane");
    h->SetAxisRange(0,5,"Z");
    h->Draw("COLZ");
    
    c.SaveAs("debug.png");
}