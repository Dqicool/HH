#include"genAna.h"
#include<ROOT/RCsvDS.hxx>

double myatof(std::string s){
    return std::atof(&s[0]);
}

void getCSV(){
    ROOT::RDataFrame df = ROOT::RDF::MakeCsvDataFrame("/mnt/NVME/HH/data/LepUniv_xsec.csv");
    auto ef = df.Define("xspb", myatof, {"xsectioninpb"}).Define("kfac", myatof, {"kfactor"}).Define("fileff", myatof, {"filterefficiency"});
    ef.Snapshot("hehe" , "/mnt/NVME/HH/output/weight_related.root", {"SampleID", "xspb", "kfac", "fileff"});
}


int main()
{
}