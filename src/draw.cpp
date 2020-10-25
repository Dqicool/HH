
#include <TCanvas.h>
#include <ROOT/RDataFrame.hxx>
#include <cmath>
#include "genAna.h"

//#define DEBUG

TH1D getHisto(ROOT::RDataFrame df, const char *dist, const char *weight, int nbins, double xmin, double xmax)
{
    auto histoname = "h_" + (std::string)dist;
    auto h = df.Histo1D({&histoname[0], "", nbins, xmin, xmax}, dist, weight);
    return (*h);
}
void draw(const char *in_file, const char *out_histo)
{
    ROOT::EnableImplicitMT();
    //Read files
    ROOT::RDataFrame df("NOMINAL", in_file);

    //book histos
    TH1::SetDefaultSumw2();
    TH1D h_mbb, h_drbb, h_mtauvis, h_mtaucol;
    h_mbb = getHisto(df, "M_bb", "weight", 500, 0, 1000);
    h_drbb = getHisto(df, "Delta_R_bb", "weight", 500, 0,10);
    h_mtauvis = getHisto(df, "M_tau_vis", "weight", 500,0,1000);
    h_mtaucol = getHisto(df, "M_tau_col", "weight", 500,0,2000);

    //save
    TFile *out = TFile::Open(out_histo, "recreate");

    h_mbb.Write();
    h_drbb.Write();
    h_mtauvis.Write();
    h_mtaucol.Write();
    out->Close();
}
#ifndef DEBUG
int main(int argc, char **argv)
{
    const char *in_file = argv[1];
    const char *out_histo = argv[2];

    draw(in_file, out_histo);
}
#else
int main()
{
    const char *in_file = "output/sel_out/361104.txt.root";
    const char *out_histo = "debug.root";
    draw(in_file, out_histo);
}
#endif