#include "genAna.h"

//#define DEBUG

#define LUMI (58.4501+43.5873+36.2369) //fb-1

void copyMeta(const char* infile, const char* outfile)
{
    TFile* in = TFile::Open(infile, "READ");
    TH1D* tmp = (TH1D*)in->Get("h_metadata");
    TFile* out = TFile::Open(outfile, "UPDATE");
    tmp->Write();
    out->Close();
    in->Close();
}

struct weight_rela
{
    double xspb;
    double kfac;
    double fileff;
    double sumofw;
};


struct weight_rela getWeightComps(const char* infile_name){
    std::string sentence(infile_name);
    std::istringstream iss(sentence);
    std::vector<std::string> tokens;
    std::string token;
    while (std::getline(iss, token, '.')) 
    {
        if (!token.empty())
            tokens.push_back(token);
    }
    std::string sentence2(tokens[0]);
    std::istringstream iss2(sentence2);
    std::vector<std::string> tokens2;
    std::string token2;
    while (std::getline(iss2, token2, '/')) 
    {
        if (!token2.empty())
            tokens2.push_back(token2);
    }
    ROOT::RDataFrame df("hehe", "/mnt/NVME/HH/output/weight_related.root");
    auto ef = df.Filter([token2](std::string s){return s.compare(token2) == 0;}, {"SampleID"});
    auto xspb = (ef.Take<double>("xspb").GetValue())[0];
    auto kfac = (ef.Take<double>("kfac").GetValue())[0];
    auto fileff = (ef.Take<double>("fileff").GetValue())[0];

    TFile* in = TFile::Open(infile_name, "READ");
    TH1D* tmp = (TH1D*)in->Get("h_metadata");
    double sumofw = tmp->GetBinContent(8);
    in->Close();
    
    return {xspb, kfac, fileff, sumofw};
}

void selection(const char* infile, const char* outfile)
{
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df("NOMINAL", infile);
    
    auto ef = df
        .Filter([](ROOT::Math::PtEtaPhiMVector v){return v.Pt() > 30;}, {"bjet_0_p4_n"})
        .Filter([](ROOT::Math::PtEtaPhiMVector v){return v.Pt() > 30;}, {"bjet_1_p4_n"})
        .Filter([](ROOT::Math::PtEtaPhiMVector v, uint is_tau){return v.Pt() > 10 && is_tau;}, {"tau_0_p4_n", "tau_0_jet_rnn_tight"})
        .Filter([](ROOT::Math::PtEtaPhiMVector v){return v.Pt() > 10;}, {"met_reco_p4_n"})
        .Define("m_bb",[](ROOT::Math::PtEtaPhiMVector vb0, ROOT::Math::PtEtaPhiMVector vb1){return (vb0+vb1).M();}, {"bjet_0_p4_n", "bjet_1_p4_n"})
        .Define("delta_R_bb", [](ROOT::Math::PtEtaPhiMVector a, ROOT::Math::PtEtaPhiMVector b){return std::sqrt((a.eta()-b.eta())*(a.eta()-b.eta()) + (a.phi()-b.phi())*(a.phi()-b.phi()));}, {"bjet_0_p4_n", "bjet_1_p4_n"} )
        .Define("lep_type", getLepType, {"elec_0", "muon_0"})
        .Define("lep_0_p4_n", getLepP4, {"lep_type", "elec_0_p4_n","muon_0_p4_n"})
        .Define("v_position", getVPos, {"met_reco_p4_n", "tau_0_p4_n", "lep_0_p4_n", "lep_type"})
        .Define("vh_p4_n", getVHP4, {"v_position", "met_reco_p4_n", "tau_0_p4_n", "lep_0_p4_n"})
        .Define("vl_p4_n", getVLP4, {"v_position", "met_reco_p4_n", "tau_0_p4_n", "lep_0_p4_n"})
        .Define("m_tau_vis", mTauTauVis, {"tau_0_p4_n","elec_0_p4_n","muon_0_p4_n","elec_0", "muon_0"})
        .Define("m_tau_col", mTauTauColin, {"met_reco_p4_n", "tau_0_p4_n","elec_0_p4_n","muon_0_p4_n",})
        .Filter([](double mvis){return mvis > 0;}, {"m_tau_col"});

    if(df.HasColumn("weight_mc")){
        weight_rela weight_comp = getWeightComps(infile);
        std::cout<<infile<<"\t"<<weight_comp.xspb<<"\t"<<weight_comp.kfac<<"\t"<<weight_comp.fileff<<"\t"<<weight_comp.sumofw<<"\n";
        auto ff = ef.Define("weight", [weight_comp](double weight_mc){return LUMI*weight_comp.xspb*1000*weight_comp.fileff*weight_comp.kfac/weight_comp.sumofw*weight_mc;}, {"weight_mc"});
        ff.Snapshot("NOMINAL", outfile);

    }
    else
    {
        auto ff = ef.Define("weight", []{return (double)1.0;});
        ff.Snapshot("NOMINAL", outfile);
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