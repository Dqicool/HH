
#include <TCanvas.h>
#include <ROOT/RDataFrame.hxx>
#include <cmath>
#include <TH2D.h>
#include "genAna.h"

#define DEBUG

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
    auto ef = df
            .Define("log2_bjet0_pt_sf", [](double chi0){return std::log2(chi0);}, {"bjet_0_pt_scale_fac"})
            .Define("log2_bjet1_pt_sf", [](double chi0){return std::log2(chi0);}, {"bjet_1_pt_scale_fac"})
            ;
    

    //book histos
    TH1::SetDefaultSumw2();
    //TH1D h_mbb;
    TH2D h_sf = (ef.Histo2D({"h_b_jet_scale_fac", "#chi_{0} vs #chi_{1}", 10, -1, 1, 10, -1, 1}, {"log2_bjet0_pt_sf"}, {"log2_bjet1_pt_sf"}, {"weight"})).GetValue();
    //h_mbb = getHisto(df, "m_bb", "weight", 50, 0, 1000);


    //save
    TFile *out = TFile::Open(out_histo, "recreate");

    //h_mbb.Write();
    TCanvas c1("c1","c1",2000,2000);
        h_sf.GetXaxis()->SetTitle("log_{2}(#chi_{0})");
        h_sf.GetYaxis()->SetTitle("log_{2}(#chi_{1})");
        h_sf.SetStats(0);
        h_sf.Draw("COLZ");
        //h_sf.Draw("TEXT SAME");
    c1.SaveAs("debugtt.png");

    h_sf.Write();

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
    //const char *in_file = "output/03_sel_out/363356.txt.root";
    const char *in_file = "output/03_sel_out/410470.txt.root";
    const char *out_histo = "debug.root";
    draw(in_file, out_histo);
}
#endif