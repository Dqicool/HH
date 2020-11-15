#include "genAna.h"
#define LUMI 139e3
#include<THStack.h>
#include<TLegend.h>
#define IMAGEX 3000
#define IMAGEY 2000

#include <sys/stat.h>
#include <unistd.h>
#include <string>
#include <fstream>


//#define debug
using namespace std;

bool file_exists_test (const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

TH1D * getDraw(const char* file, const char* dist,  Color_t col)
{
    TH1D *h = NULL;
    if(file_exists_test(file))
    {
        auto histoname = (string)dist;
        TFile *f = TFile::Open(file, "read");
        h = (TH1D*)f->Get(&histoname[0]);
        h->SetFillColor(col);
        h->SetLineColor(col);
        h->SetMarkerColor(col);
    }
    return h;
}


void plotProp(const char* outfile, const char* dist, bool truth){
    ROOT::EnableImplicitMT(24);
    TH1::SetDefaultSumw2();
    //loading data from files
        bool onlyZZ=0;
        bool onlyTT=0;
        std::vector<std::vector<std::string>> files;
        std::vector<Color_t> color_vec;
        std::vector<std::string> cata;
        std::vector<std::string> color_name;
        //SINGLE_V
            std::vector<std::string> single_V;
            if(!onlyZZ && !onlyTT){
                single_V.push_back("output/04_plot_out/361100.txt.root");
                single_V.push_back("output/04_plot_out/361101.txt.root");
                single_V.push_back("output/04_plot_out/361102.txt.root");
                single_V.push_back("output/04_plot_out/361103.txt.root");
                single_V.push_back("output/04_plot_out/361104.txt.root");
                single_V.push_back("output/04_plot_out/361105.txt.root");
                single_V.push_back("output/04_plot_out/361106.txt.root");
                single_V.push_back("output/04_plot_out/361107.txt.root");
                single_V.push_back("output/04_plot_out/361108.txt.root");
            }

        //ZZ
            std::vector<std::string> zz;
            if(!onlyTT){
                zz.push_back("output/04_plot_out/363355.txt.root");
                zz.push_back("output/04_plot_out/363356.txt.root");
            }
            
        //WZ
            std::vector<std::string> wz;
            if(!onlyZZ && !onlyTT){
                wz.push_back("output/04_plot_out/363357.txt.root");
                wz.push_back("output/04_plot_out/363358.txt.root");
                wz.push_back("output/04_plot_out/363489.txt.root");
            }       
        //ztt
            std::vector<std::string> ztt;
            if(!onlyZZ && !onlyTT){
                ztt.push_back("output/04_plot_out/364128.txt.root");
                ztt.push_back("output/04_plot_out/364129.txt.root");
                ztt.push_back("output/04_plot_out/364130.txt.root");
                ztt.push_back("output/04_plot_out/364131.txt.root");
                ztt.push_back("output/04_plot_out/364132.txt.root");
                ztt.push_back("output/04_plot_out/364133.txt.root");
                ztt.push_back("output/04_plot_out/364134.txt.root");
                ztt.push_back("output/04_plot_out/364135.txt.root");
                ztt.push_back("output/04_plot_out/364136.txt.root");
                ztt.push_back("output/04_plot_out/364137.txt.root");
                ztt.push_back("output/04_plot_out/364138.txt.root");
                ztt.push_back("output/04_plot_out/364139.txt.root");
                ztt.push_back("output/04_plot_out/364140.txt.root");
                ztt.push_back("output/04_plot_out/364141.txt.root");
            }
        //qq4l
            std::vector<std::string> qq4l;
            if(!onlyZZ && !onlyTT){
                qq4l.push_back("output/04_plot_out/364250.txt.root");
                qq4l.push_back("output/04_plot_out/364253.txt.root");
                qq4l.push_back("output/04_plot_out/364254.txt.root");
                qq4l.push_back("output/04_plot_out/364255.txt.root");
            }
        //ttbar
            std::vector<std::string> ttbar;
            if(!onlyZZ){
                ttbar.push_back("output/04_plot_out/410470.txt.root");
            }
        //single_top
            std::vector<std::string> single_t;
            if(!onlyZZ && !onlyTT){
                single_t.push_back("output/04_plot_out/410644.txt.root");
                single_t.push_back("output/04_plot_out/410645.txt.root");
                single_t.push_back("output/04_plot_out/410658.txt.root");
                single_t.push_back("output/04_plot_out/410659.txt.root");
            }
        //Wt
            std::vector<std::string> wt;
            if(!onlyZZ && !onlyTT){
                single_t.push_back("output/04_plot_out/410646.txt.root");
                single_t.push_back("output/04_plot_out/410647.txt.root");
            }
        //DATA
            std::vector<std::string> data;
            if(!onlyZZ && !onlyTT){
                data.push_back("output/04_plot_out/data15_13TeV.txt.root");
                data.push_back("output/04_plot_out/data16_13TeV.txt.root");
                data.push_back("output/04_plot_out/data17_13TeV.txt.root");
                data.push_back("output/04_plot_out/data18_13TeV.txt.root");
            }
            

        files =         { zz,       single_V,     wz,     ztt,        qq4l,       single_t,   wt,        ttbar };
        cata =          {"zz",      "single_V",   "wz",   "ztt",      "qq4l",     "single_t", "wt",       "ttbar" };
        color_vec =     { kBlue,    kOrange,       kPink,  kMagenta,   kViolet,    kRed,      kAzure,     kCyan};
        color_name =    {"Blue",    "Orange",      "Pink", "Magenta",  "Violet",   "Red",     "Azure",    "Cyan"};
        
    //stack MCs
        string st_n = (string)dist + "_stack";
        THStack* stack = new THStack(&st_n[0], "");
        string h_n = (string)dist + "_hist";
        std::vector<TH1D> histo_mc;

        for(uint i = 0; i < files.size(); i++)
        {
            for (uint j = 0; j < (files[i]).size(); j++)
            {
                auto h = getDraw(&(files[i][j])[0], dist, color_vec[i]);
                
                if (h != NULL) histo_mc.push_back(*h);
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
        {
            auto h = getDraw(entry.data(), dist, kBlack);
            if (&h != NULL) histo_data.push_back(*h);
        }
            
        TH1D * h_inc_data = new TH1D(&h_n_data[0],"",nbins,xmin,xmax);
        for(auto h:histo_data)
            h_inc_data->Add(&h);
        
    //Store
        TFile * out = TFile::Open(outfile,"recreate");
        stack->Write();
        h_inc->Write();
        h_inc_data->Write();
        out->Close();

        TH1D* h_scaled_zz = (TH1D*)(histo_mc[0].Clone("h_scaled_zz"));
        h_scaled_zz->Scale(150);
        h_scaled_zz->SetFillColorAlpha(kBlue, 0);
        h_scaled_zz->SetStats(0);
        h_scaled_zz->SetLineWidth(3);

    // //print stat
    //     cout<<"det\t"<<h_inc_cr->Integral()<<"\t\t"<<(cr_histo_mc[0]).Integral() + (cr_histo_mc[1]).Integral()<<"\t\t"<<h_inc_cr->Integral()/((cr_histo_mc[0]).Integral() + (cr_histo_mc[1]).Integral())<<":1"<<endl;
    //     cout<<"cut\t"<<h_inc_cut->Integral()<<"\t\t"<<(cut_histo_mc[0]).Integral() + (cut_histo_mc[1]).Integral()<<"\t\t"<<h_inc_cut->Integral()/((cut_histo_mc[0]).Integral() + (cut_histo_mc[1]).Integral())<<":1"<<endl;
    //     cout<<"sr\t"<<h_inc_sr->Integral()<<"\t\t"<<(sr_histo_mc[0]).Integral() + (sr_histo_mc[1]).Integral()<<"\t\t"<<h_inc_sr->Integral()/((sr_histo_mc[0]).Integral() + (sr_histo_mc[1]).Integral())<<":1"<<endl;;

    //drawstack
        //sr plot
            TCanvas c1("c1","",IMAGEX,IMAGEY);
            stack->SetTitle(st_n.data());
            h_inc->SetMarkerColor(kBlack);
            h_inc->SetLineColor(kBlack);
            h_inc->SetFillColor(kBlack);
            h_inc->SetFillStyle(3017);
            //h_inc->SetAxisRange(0,YMAX,"Y");
            h_scaled_zz->SetTitle(dist);
            h_scaled_zz->Draw("hist");
            stack->Draw("hist, same");
            h_inc->Draw("E2, same");
            
            h_inc_data->SetMarkerStyle(kFullCircle);
            h_inc_data->SetMarkerSize(2.0f);
            h_inc_data->SetLineWidth(3.0f);
            h_inc_data->SetLineColor(kBlack);
            h_inc_data->Draw("E1, SAME");
            h_scaled_zz->Draw("hist, same");


            TLegend legend(0.6,0.7,0.9,0.9);
            legend.AddEntry(h_inc_data, "data", "lpfe");
            legend.AddEntry(&histo_mc.back(), "ttbar", "lpf");
            legend.AddEntry(&histo_mc[0], "zz", "lpf");
            legend.AddEntry(h_scaled_zz, "zz x 150", "l");
            legend.Draw();

            auto save_name = "plots/" + (string)dist + "_stack.png";
            if(onlyTT) save_name = "plots/onlyTT/" + (string)dist + "_stack.png";
            else if(onlyZZ) save_name = "plots/onlyZZ/" + (string)dist + "_stack.png";
            c1.SaveAs(&save_name[0]);
        // output
        double xs, errerr;
        xs = h_inc->IntegralAndError(1,50, errerr);
        cout<<"*************************"<<endl;
        cout<<"XS = "<<xs<<"\t"<<errerr<<endl;
        cout<<"*************************"<<endl;

}

int main()
{
    cout<<"mbb:"<<endl;
    //plotProp("output/05_stack_out/score.root", "score", 0);
    plotProp("output/05_stack_out/score_mbb.root", "score_mbb", 0);
    plotProp("output/05_stack_out/score_mtt.root", "score_mtt", 0);
    plotProp("output/05_stack_out/score_met.root", "score_met", 0);
    plotProp("output/05_stack_out/score_chi.root", "score_chi", 0);
}
