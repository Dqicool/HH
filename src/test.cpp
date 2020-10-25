#include"genAna.h"
#include<ROOT/RCsvDS.hxx>

double myatof(std::string s){
    return std::atof(&s[0]);
}

void getCSV(){
    ROOT::RDataFrame df = ROOT::RDF::MakeCsvDataFrame("/mnt/NVME/HH/data/LepUniv_xsec.csv");
    auto ef = df.Define("xspb", myatof, {"xsectioninpb"}).Define("kfac", myatof, {"kfactor"}).Define("fileff", myatof, {"filterefficiency"});
    //ef.Snapshot("hehe" , "/mnt/NVME/HH/output/weight_related.root", {"SampleID", "xspb", "kfac", "fileff"});
}

std::vector<double> getWeightComps(const char* infile_name){
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
    ef.Snapshot("bababa", "debug.root");
    auto xspb = (ef.Take<double>("xspb").GetValue())[0];
    auto kfac = (ef.Take<double>("kfac").GetValue())[0];
    auto fileff = (ef.Take<double>("fileff").GetValue())[0];
    
    return {xspb, kfac, fileff};
}

int main()
{
    //getCSV();
    //auto b = getWeightComps("/mnt/NVME/HH/output/sel_out/361100.txt.root");
    //for(auto a:b){
        //std::cout<<a<<"\n";
    //}
    double dphi_hl = M_PI * 0.3333333;
    double dphi_hv = M_PI * 0.33333333/2;
    double vhpt = std::cos(std::abs(dphi_hv)) - std::sin(std::abs(dphi_hv)) * (1/std::tan(std::abs(dphi_hl)));
    double vlpt = std::sin(std::abs(dphi_hv)) / std::sin(std::abs(dphi_hl));

    std::cout<<vhpt<<"\t"<<vlpt<<"\n";
}