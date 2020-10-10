#include <TROOT.h>
#include <TH1.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TMath.h>
#include <ROOT/RDataFrame.hxx>
#include <Math/GenVector/LorentzVector.h>
#include <Math/GenVector/PtEtaPhiM4D.h>
#include <Math/GenVector/PtEtaPhiM4Dfwd.h>
#include <Math/Vector4Dfwd.h>
#include <Math/Vector3Dfwd.h>
#include <Math/GenVector/Boost.h>
#include <TChain.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TLatex.h>

#include "TMatrixD.h"
#include "TVectorD.h"
#include "TMath.h"
#include "Math/Point3D.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/GenVector/Rotation3D.h"
#include "Math/GenVector/EulerAngles.h"
#include "Math/GenVector/AxisAngle.h"
#include "Math/GenVector/Quaternion.h"
#include "Math/GenVector/LorentzRotation.h"
#include "Math/GenVector/Boost.h"
#include "Math/GenVector/Transform3D.h"
#include "Math/GenVector/Plane3D.h"
#include "Math/GenVector/VectorUtil.h"
#include <Math/QuantFuncMathMore.h>
#include <TLorentzVector.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#define Z_MASS 91.1876e3
#define H_MASS 125.18e3
#define GeV 1e3
