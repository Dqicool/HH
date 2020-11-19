#include "selection.h"
#include <chrono>

//#define DEBUG
#define MT

#define LUMI (58.4501+43.5873+36.2369) //fb-1

void selection(const char* infile, const char* outfile)
{
    #ifdef MT 
    ROOT::EnableImplicitMT(4);
    #endif
    ROOT::RDataFrame df("NOMINAL", infile);
    
    auto ef = df
                .Define("evt", getEvent, {"bjet_0_p4_n", "bjet_1_p4_n", "lep_0_p4_n", "tau_0_p4_n", "met_reco_p4_n"})
                .Define("gpu_pass",[](struct event evt){ return !evt.veto;}, {"evt"})
                .Filter([](bool pass){ return pass;}, {"gpu_pass"})
                .Define("bjet_0_pt_scale_fac", [](struct event evt){return evt.bjet0_pt_scale_fac;}, {"evt"})
                .Define("bjet_1_pt_scale_fac", [](struct event evt){return evt.bjet1_pt_scale_fac;}, {"evt"})
                .Define("score", [](struct event evt){return evt.score;}, {"evt"})
                .Define("score_m_bb", [](struct event evt){return evt.score_m_bb;}, {"evt"})
                .Define("score_m_tt", [](struct event evt){return evt.score_m_tt;}, {"evt"})
                .Define("score_met", [](struct event evt){return evt.score_met;}, {"evt"})
                .Define("score_chi", [](struct event evt){return evt.score_chi;}, {"evt"})
                .Filter("")
        ;

    if(df.HasColumn("weight_mc"))
    {
        weight_rela weight_comp = getWeightComps(infile);
        std::cout<<infile<<"\t"<<weight_comp.xspb<<"\t"<<weight_comp.kfac<<"\t"<<weight_comp.fileff<<"\t"<<weight_comp.sumofw<<"\n";
        auto ff = ef.Define("weight", [weight_comp](double weight_mc, double sf){return LUMI*weight_comp.xspb*1000*weight_comp.fileff*weight_comp.kfac/weight_comp.sumofw*weight_mc*sf;}, {"weight_mc", "total_sf"});
        std::vector<std::string> snap_column = {"score","score_m_bb", "score_m_tt", "score_met", "score_chi", "bjet_0_pt_scale_fac", "bjet_1_pt_scale_fac", "weight", "gpu_pass", "bjet_0_p4_n", "bjet_1_p4_n", "lep_0_p4_n", "tau_0_p4_n", "met_reco_p4_n"};
        ff.Snapshot("NOMINAL", outfile, snap_column);
    }
    else
    {
        auto ff = ef.Define("weight", []{return (double)1.0;});
        std::vector<std::string> snap_column = {"score", "score_m_bb", "score_m_tt", "score_met", "score_chi", "bjet_0_pt_scale_fac", "bjet_1_pt_scale_fac", "weight"};
        ff.Snapshot("NOMINAL", outfile, snap_column);
    }
    std::cout<<infile<<"...\t\t\tdone\n";
}

#ifdef DEBUG
int main()
{
    auto start = std::chrono::steady_clock::now();
    //selection("output/02_presel_out/363356.txt.root","debug.root");
    selection("output/02_presel_out/410470.txt.root","debug.root");
    //copyMeta("output/convert_out/361104.txt.root","debug.root");
    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;
    std::cout << (double)std::chrono::duration_cast<std::chrono::microseconds>( end - start ).count()/1000000.<< "s" << std::endl;
    return 0;
}
#else
int main(int argc, char** argv)
{
    selection(argv[1], argv[2]);
    return 0;
}
#endif