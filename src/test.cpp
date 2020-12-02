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

std::vector<std::string> getFileName(std::string path){
    std::vector<std::string> ret;
    for (const auto & entry : std::filesystem::directory_iterator(path))
    {
        if(entry.path() != path + "/meta.txt")
            ret.push_back(entry.path());
    }
    return ret;
}

int main()
{
    auto a = getFileName("/mnt/NVME/HH/data/diego.v24.mc16_13TeV.410659.PhPy8EG_A14_tchan_BW50_lept_antitop.deriv.DAOD_STDM4.e6671_s3126_r9364_p4097");
    for (auto i:a){
        std::cout<<i<<"\n";
    }
}