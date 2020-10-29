#include "selection.h"

//#define DEBUG

#define LUMI (58.4501+43.5873+36.2369) //fb-1

void selection(const char* infile, const char* outfile)
{
    #ifndef DEBUG
    ROOT::EnableImplicitMT();
    #endif
    ROOT::RDataFrame df("NOMINAL", infile);
    
    auto ef = df
        .Define("m_bb",[](ROOT::Math::PtEtaPhiMVector vb0, ROOT::Math::PtEtaPhiMVector vb1){return (vb0+vb1).M();}, {"bjet_0_p4_n", "bjet_1_p4_n"})
        .Define("delta_R_bb", getDeltaR, {"bjet_0_p4_n", "bjet_1_p4_n"} )
        .Define("v_position", getVPos, {"met_reco_p4_n", "tau_0_p4_n", "lep_0_p4_n"})
        .Define("vh_p4_n", getVHP4, {"v_position", "met_reco_p4_n", "tau_0_p4_n", "lep_0_p4_n"})
        .Define("vl_p4_n", getVLP4, {"v_position", "met_reco_p4_n", "tau_0_p4_n", "lep_0_p4_n"})
        .Define("m_tau_vis", getMTauTauVis, {"tau_0_p4_n","elec_0_p4_n","muon_0_p4_n","elec_0", "muon_0"})
        .Define("m_tau_col", getMTauTauColin, {"tau_0_p4_n", "lep_0_p4_n", "vh_p4_n", "vl_p4_n", "v_position"})
        .Filter([](double mtau){return mtau > 0;}, {"m_tau_col"})
        .Define("delta_R_tautau", getDeltaR, {"lep_0_p4_n", "tau_0_p4_n"})
        .Define("m_4_body", getM4Body, {"tau_0_p4_n", "lep_0_p4_n", "vh_p4_n", "vl_p4_n", "bjet_0_p4_n", "bjet_1_p4_n"})
        .Define("top_type", getTopType, {"tau_0_p4_n", "lep_0_p4_n", "vh_p4_n", "vl_p4_n", "bjet_0_p4_n", "bjet_1_p4_n"})
        .Define("top_0_p4_n", getTop0P4,   {"top_type", "tau_0_p4_n", "lep_0_p4_n", "vh_p4_n", "vl_p4_n", "bjet_0_p4_n", "bjet_1_p4_n"})
        .Define("top_1_p4_n", getTop1P4,   {"top_type", "tau_0_p4_n", "lep_0_p4_n", "vh_p4_n", "vl_p4_n", "bjet_0_p4_n", "bjet_1_p4_n"})
        .Define("m_top_0", [](ROOT::Math::PtEtaPhiMVector top0){return top0.M();}, {"top_0_p4_n"})
        .Define("m_top_1", [](ROOT::Math::PtEtaPhiMVector top1){return top1.M();}, {"top_1_p4_n"})
        .Define("top_0_deltaR_con", getTop0DeltaRBTau, {"top_type", "tau_0_p4_n", "bjet_0_p4_n", "bjet_1_p4_n"})
        .Define("top_1_deltaR_con", getTop1DeltaRBTau, {"top_type", "lep_0_p4_n", "bjet_0_p4_n", "bjet_1_p4_n"})
        
        //.Filter("m_bb > 60 && m_bb < 110")
        //.Filter("m_tau_col > 50 && m_tau_col < 150")
        //.Filter("delta_R_tautau < 2")
        // .Filter("top_0_deltaR_con > 2")
        // .Filter("top_1_deltaR_con > 2")
        ;

    if(df.HasColumn("weight_mc"))
    {
        weight_rela weight_comp = getWeightComps(infile);
        std::cout<<infile<<"\t"<<weight_comp.xspb<<"\t"<<weight_comp.kfac<<"\t"<<weight_comp.fileff<<"\t"<<weight_comp.sumofw<<"\n";
        auto ff = ef.Define("weight", [weight_comp](double weight_mc, double sf){return LUMI*weight_comp.xspb*1000*weight_comp.fileff*weight_comp.kfac/weight_comp.sumofw*weight_mc*sf;}, {"weight_mc", "total_sf"});
        std::vector<std::string> snap_column = {"m_bb", "delta_R_bb", "m_tau_vis", "m_tau_col", "delta_R_tautau", "m_4_body", "m_top_0", "m_top_1", "top_0_deltaR_con", "top_1_deltaR_con", "weight"};
        ff.Snapshot("NOMINAL", outfile, snap_column);
    }
    else
    {
        auto ff = ef.Define("weight", []{return (double)1.0;});
        std::vector<std::string> snap_column = {"m_bb", "delta_R_bb", "m_tau_vis", "m_tau_col", "delta_R_tautau", "m_4_body", "m_top_0", "m_top_1", "top_0_deltaR_con", "top_1_deltaR_con", "weight"};
        ff.Snapshot("NOMINAL", outfile, snap_column);
    }
}

#ifdef DEBUG
int main()
{
    selection("output/convert_out/410470.txt.root","debug.root");
    //copyMeta("output/convert_out/361104.txt.root","debug.root");
    return 0;
}
#else
int main(int argc, char** argv)
{
    selection(argv[1], argv[2]);
    return 0;
}
#endif