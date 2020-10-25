#include "genAna.h"
#define LUMI 139e3
#include<THStack.h>
#define IMAGEX 3000
#define IMAGEY 2000

//#define debug
using namespace std;
TH1D getDraw(const char* file, const char* dist,  Color_t col)
{
    auto histoname = "h_" + (string)dist;
    TFile *f = TFile::Open(file, "read");
    TH1D *h = (TH1D*)f->Get(&histoname[0]);
    h->SetFillColor(col);
    h->SetLineColor(col);
    h->SetMarkerColor(col);
    return (*h);
}


void plotProp(const char* outfile, const char* dist, bool truth){
    ROOT::EnableImplicitMT(24);
    TH1::SetDefaultSumw2();
    //loading data from files
        std::vector<std::vector<std::string>> files;
        std::vector<Color_t> color_vec;
        std::vector<std::string> cata;
        std::vector<std::string> color_name;
        //SINGLE_V
            std::vector<std::string> single_V;
            single_V.push_back("output/plot_out/361100.txt.root");
            //single_V.push_back("output/plot_out/361101.txt.root");
            single_V.push_back("output/plot_out/361102.txt.root");
            single_V.push_back("output/plot_out/361103.txt.root");
            single_V.push_back("output/plot_out/361104.txt.root");
            single_V.push_back("output/plot_out/361105.txt.root");
            single_V.push_back("output/plot_out/361106.txt.root");
            single_V.push_back("output/plot_out/361107.txt.root");
            //single_V.push_back("output/plot_out/361108.txt.root");

        //ZZ
            std::vector<std::string> zz;
            #ifndef debug
            //zz.push_back("output/plot_out/363355.txt.root");
            zz.push_back("output/plot_out/363356.txt.root");
            #endif
        //WZ
            std::vector<std::string> wz;
            #ifndef debug
            //wz.push_back("output/plot_out/363357.txt.root");
            wz.push_back("output/plot_out/363358.txt.root");
            wz.push_back("output/plot_out/363489.txt.root");
            #endif
        //ztt
            std::vector<std::string> ztt;
            #ifndef debug
            //ztt.push_back("output/plot_out/364128.txt.root");
            ztt.push_back("output/plot_out/364129.txt.root");
            ztt.push_back("output/plot_out/364130.txt.root");
            ztt.push_back("output/plot_out/364131.txt.root");
            ztt.push_back("output/plot_out/364132.txt.root");
            ztt.push_back("output/plot_out/364133.txt.root");
            ztt.push_back("output/plot_out/364134.txt.root");
            ztt.push_back("output/plot_out/364135.txt.root");
            ztt.push_back("output/plot_out/364136.txt.root");
            ztt.push_back("output/plot_out/364137.txt.root");
            ztt.push_back("output/plot_out/364138.txt.root");
            ztt.push_back("output/plot_out/364139.txt.root");
            ztt.push_back("output/plot_out/364140.txt.root");
            ztt.push_back("output/plot_out/364141.txt.root");
            #endif
        //qq4l
            std::vector<std::string> qq4l;
            #ifndef debug
            qq4l.push_back("output/plot_out/364250.txt.root");
            qq4l.push_back("output/plot_out/364253.txt.root");
            qq4l.push_back("output/plot_out/364254.txt.root");
            qq4l.push_back("output/plot_out/364255.txt.root");
            #endif
        //ttbar
            std::vector<std::string> ttbar;
            #ifndef debug
            ttbar.push_back("output/plot_out/410470.txt.root");
            #endif
        //single_top
            std::vector<std::string> single_t;
            #ifndef debug
            single_t.push_back("output/plot_out/410644.txt.root");
            single_t.push_back("output/plot_out/410645.txt.root");
            single_t.push_back("output/plot_out/410658.txt.root");
            single_t.push_back("output/plot_out/410659.txt.root");
            #endif
        //Wt
            std::vector<std::string> wt;
            #ifndef debug
            single_t.push_back("output/plot_out/410646.txt.root");
            single_t.push_back("output/plot_out/410647.txt.root");
            #endif
        //DATA
            std::vector<std::string> data;
            data.push_back("output/plot_out/data15_13TeV.txt.root");
            data.push_back("output/plot_out/data16_13TeV.txt.root");
            data.push_back("output/plot_out/data17_13TeV.txt.root");
            data.push_back("output/plot_out/data18_13TeV.txt.root");
            

        files =         { single_V,      zz,     wz,     ztt,        qq4l,       single_t,   wt,        ttbar };
        cata =          {"single_V",    "zz",   "wz",   "ztt",      "qq4l",     "single_t", "wt",       "ttbar" };
        color_vec =     { kOrange,       kRed,   kPink,  kMagenta,   kViolet,    kBlue,      kAzure,     kCyan};
        color_name =    {"kOrange",     "kRed", "kPink","kMagenta", "kViolet",  "kBlue",    "kAzure",   "kCyan"};
        
    //stack MCs
        string st_n = (string)dist + "_stack";
        THStack* stack = new THStack(&st_n[0], "");
        string h_n = (string)dist + "_hist";
        std::vector<TH1D> histo_mc;

        for(uint i = 0; i < files.size(); i++)
        {
            for (uint j = 0; j < (files[i]).size(); j++)
            {
                histo_mc.push_back(getDraw(&(files[i][j])[0], dist, color_vec[i]));
            }
            std::cout<<cata[i] + " is in color " + color_name[i]<<endl;
        }
        
        auto nbins = (histo_mc[0]).GetNbinsX();
        auto xmin  = (histo_mc[0]).GetBinLowEdge(1);
        auto xmax =  (histo_mc[0]).GetBinLowEdge(nbins) + (histo_mc[0]).GetBinWidth(nbins);

        TH1D * h_inc = new TH1D(&h_n[0],"",nbins,xmin,xmax);
        
        for(uint i=0; i<histo_mc.size(); i++)
        {
            stack->Add(&(histo_mc[i]));
            h_inc->Add(&(histo_mc[i]));
        }
    //stack data
        std::vector<TH1D> histo_data;
        string h_n_data = (string)dist + "data_hist";
        for(auto entry:data)
            histo_data.push_back(getDraw(entry.data(), dist, kBlack));
        TH1D * h_inc_data = new TH1D(&h_n_data[0],"",nbins,xmin,xmax);
        for(auto h:histo_data)
            h_inc_data->Add(&h);
        
    //Store
        TFile * out = TFile::Open(outfile,"recreate");
        stack->Write();
        h_inc->Write();
        h_inc_data->Write();
        out->Close();

    // //print stat
    //     cout<<"det\t"<<h_inc_cr->Integral()<<"\t\t"<<(cr_histo_mc[0]).Integral() + (cr_histo_mc[1]).Integral()<<"\t\t"<<h_inc_cr->Integral()/((cr_histo_mc[0]).Integral() + (cr_histo_mc[1]).Integral())<<":1"<<endl;
    //     cout<<"cut\t"<<h_inc_cut->Integral()<<"\t\t"<<(cut_histo_mc[0]).Integral() + (cut_histo_mc[1]).Integral()<<"\t\t"<<h_inc_cut->Integral()/((cut_histo_mc[0]).Integral() + (cut_histo_mc[1]).Integral())<<":1"<<endl;
    //     cout<<"sr\t"<<h_inc_sr->Integral()<<"\t\t"<<(sr_histo_mc[0]).Integral() + (sr_histo_mc[1]).Integral()<<"\t\t"<<h_inc_sr->Integral()/((sr_histo_mc[0]).Integral() + (sr_histo_mc[1]).Integral())<<":1"<<endl;;

    //drawstack
        //sr plot
            TCanvas c1("c1","",IMAGEX,IMAGEY);
            h_inc->SetMarkerColor(kBlack);
            h_inc->SetLineColor(kBlack);
            h_inc->SetFillColor(kBlack);
            h_inc->SetFillStyle(3017);
            //h_inc->SetAxisRange(0,YMAX,"Y");
            stack->Draw("hist");
            h_inc->Draw("E2, same");
            h_inc_data->SetMarkerStyle(kFullCircle);
            h_inc_data->SetMarkerSize(2.0f);
            h_inc_data->SetLineWidth(3.0f);
            h_inc_data->Draw("E1, SAME");
            auto save_name = "plots/" + (string)dist + "_stack.png";
            c1.SaveAs(&save_name[0]);
        // // output
        // double xs, errerr;
        // xs = h_inc->IntegralAndError(1,1000, errerr);
        // cout<<"*************************"<<endl;
        // cout<<"XS = "<<xs<<"\t"<<errerr<<endl;
        // cout<<"*************************"<<endl;

}

int main()
{
    cout<<"mbb:"<<endl;
    plotProp("output/stack_out/mbb.root", "M_bb", 0);
    cout<<"drbb:"<<endl;
    plotProp("output/stack_out/drbb.root", "Delta_R_bb", 0);
    cout<<"mtauvis:"<<endl;
    plotProp("output/stack_out/mtauvis.root", "M_tau_vis", 0);
    cout<<"mtaucol:"<<endl;
    plotProp("output/stack_out/mtaucol.root", "M_tau_col", 0);
    
}
