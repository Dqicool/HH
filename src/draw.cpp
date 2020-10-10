
#include <TCanvas.h>
#include <ROOT/RDataFrame.hxx>
#include <cmath>
#include "genAna.h"

#define PRWLUMI (36.21 + 58.45 + 44.31) //fb-1
#define DEBUG

TH1D getHisto(ROOT::RDF::RInterface<ROOT::Detail::RDF::RLoopManager, void> df, const char *dist, const char *weight, int nbins, double xmin, double xmax)
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
    auto ef = df.Define("bjet_0_p4_n_Pt", [](ROOT::Math::PtEtaPhiMVector v){return v.Pt();},{"bjet_0_p4_n"});
    TH1D h_mjj_cut = getHisto(ef, "bjet_0_p4_n_Pt", "weight_total", 1000, 0, 1000);

    //save
    TFile *out = TFile::Open(out_histo, "recreate");

    h_mjj_cut.Write();
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
    const char *in_file = "output/sel_out/410470.txt.root";
    const char *out_histo = "debug.root";
    draw(in_file, out_histo);
}
#endif