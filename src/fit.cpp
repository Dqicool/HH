#include "genAna.h"
#include <string>
#include <fstream>
#include <sys/stat.h>

#include "TF2.h"
#include "TH2.h"
#include "TMath.h"

bool file_exists_test (const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

TH2D * getDraw2D(const char* file, const char* dist)
{
    TH2D *h = NULL;
    if(file_exists_test(file))
    {
        auto histoname = (std::string)dist;
        TFile *f = TFile::Open(file, "read");
        h = (TH2D*)f->Get(&histoname[0]);
    }
    return h;
}

double Guass2D(double *x, double *par) {
    double xx = x[0];
    double yy = x[1];

    double A        = par[0];
    double x_0      = par[1];
    double y_0      = par[2];
    double sigma_x  = par[3];
    double sigma_y  = par[4];
    double theta    = par[5];


    double ell1 = (TMath::Power((xx - x_0) * TMath::Cos(theta) + (yy - y_0) * TMath::Sin(theta), 2)) / TMath::Power(sigma_x, 2);
    double ell2 = (TMath::Power((xx - x_0) * TMath::Sin(theta) - (yy - y_0) * TMath::Cos(theta), 2)) / TMath::Power(sigma_y, 2);
    double result = A * TMath::Exp(-ell1-ell2);
    return result;
}
 
int main() {
    std::cout<<"\n\n\n\nLINEAR\n";
    const Int_t npar = 6;
    double fparams[npar] =  {0,1.05135,1.05104,1,2,TMath::PiOver4()};
    TF2 *f = new TF2("Gauss2D", Guass2D, 0.5, 2, 0.5, 2, npar);
    // f->SetParLimits(1, 1.05135, 1.05135);
    // f->SetParLimits(2, 1.05104, 1.05104);
    f->SetParLimits(5, -TMath::Pi(), TMath::Pi());
    f->SetParLimits(3, 0, 10);
    f->SetParLimits(4, 0, 10);
    // f->SetParLimits(3, 1, 1);
    // f->SetParLimits(4, 2, 2);
    f->SetParameters(fparams);
    TH2D *h = getDraw2D("output/04_plot_out/363356.txt.root", "h_b_jet_scale_fac_lin");
    TCanvas c1("c1","",2000,2000);
    h->Fit("Gauss2D", "", "COLZ");
    f->Draw("cont1 same");
    c1.SaveAs("debug_lin.png");
    std::cout<<"chi2/ndf = "<<f->GetChisquare()/(h->GetNbinsX() * h->GetNbinsY() - 6)<<std::endl;

    std::cout<<"\n\n\n\nLOG\n";
    double fparams2[npar] = {1,0,0,1,1,TMath::PiOver4()};
    TF2 *ff = new TF2("Gauss2D2", Guass2D, -1, 1, -1, 1, npar);
    ff->SetParameters(fparams2);
    ff->SetParLimits(3, 0, 10);
    ff->SetParLimits(4, 0, 10);
    ff->SetParLimits(5, -TMath::Pi(), TMath::Pi());
    TH2D *hh = getDraw2D("output/04_plot_out/363356.txt.root", "h_b_jet_scale_fac_log");
    TCanvas c2("c2","",2000,2000);
    hh->Fit("Gauss2D2", "", "COLZ");
    ff->Draw("cont1 same");
    c2.SaveAs("debug_log.png");
    std::cout<<"chi2/ndf = "<<ff->GetChisquare()/(h->GetNbinsX() * h->GetNbinsY() - 6)<<std::endl;
    return 0;
}