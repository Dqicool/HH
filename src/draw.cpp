
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
    TH1D h_mbb, h_drbb, h_mtauvis, h_mtaucol, h_m4body, h_drtt, h_mtop0, h_mtop1, h_top0drcon, h_top1drcon;
    h_mbb = getHisto(df, "m_bb", "weight", 50, 0, 1000);
    h_drbb = getHisto(df, "delta_R_bb", "weight", 50, 0,10);
    h_mtauvis = getHisto(df, "m_tau_vis", "weight", 50,0,1000);
    h_mtaucol = getHisto(df, "m_tau_col", "weight", 50,0,2000);

    h_m4body = getHisto(df, "m_4_body", "weight", 50,0,2000);
    h_drtt = getHisto(df, "delta_R_tautau", "weight", 50, 0,10);
    h_mtop0 = getHisto(df, "m_top_0", "weight", 50,0,500);
    h_mtop1 = getHisto(df, "m_top_1", "weight", 50,0,500);
    h_top0drcon = getHisto(df, "top_0_deltaR_con", "weight", 50, 0,10);
    h_top1drcon = getHisto(df, "top_1_deltaR_con", "weight", 50, 0,10);

    //save
    TFile *out = TFile::Open(out_histo, "recreate");

    h_mbb.Write();
    h_drbb.Write();
    h_mtauvis.Write();
    h_mtaucol.Write();
    h_m4body.Write();
    h_drtt.Write();
    h_mtop0.Write();
    h_mtop1.Write();
    h_top0drcon.Write();
    h_top1drcon.Write();


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