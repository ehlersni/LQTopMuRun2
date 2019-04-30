#include "../include/ExtrapolationSRCRRunner.h"
#include <TH1.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TText.h>
#include <TPaveText.h>
#include <TGaxis.h>
#include <TFitResult.h>
#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>
#include <TVectorD.h>
#include <TDecompSVD.h>
#include <sstream>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include <string>
#include <TFile.h>

using namespace std;

//Define function
void ExtrapolationSRCRRunner::SubtractMCfromDATA(){  //bool do_ttbardy, bool do_systematics = false, TString SystType = "PU", TString SystDirection = "up"

bool do_ttbardy = true;
bool do_systematics = false;



TString path     = ExtrapolationSRCRRunner::inpathCRA;
TString path_out = ExtrapolationSRCRRunner::outpathSUB;
// cout << path << endl;
// if(do_systematics){
//   path += "/" + SystType + "_" + SystDirection;
// }
// else path += "/NOMINAL";
path += "uhh2.AnalysisModuleRunner.";
cout << "Input path: " << path+"DATA.DATA.root" << endl;
// if(do_systematics && !SystType.Contains("SCALE")) if (SystDirection != "up" && SystDirection != "down") throw runtime_error("Invalid SystDirection provided for this systematic uncertainty source.");
// if(do_systematics && SystType.Contains("SCALE")) if (SystDirection != "upup" && SystDirection != "downdown" && SystDirection != "upnone" && SystDirection != "noneup" && SystDirection != "downnone" && SystDirection != "nonedown") throw runtime_error("Invalid SystDirection provided for this systematic uncertainty source.");

//path = "/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/35867fb_DibosonNLO/HTLepInvertedSideband/uhh2.AnalysisModuleRunner.";
TFile* file_DATA_Sideband = new TFile(path+"DATA.DATA.root");
TFile* file_TTbar_Sideband = new TFile(path+"MC.TTbar.root");
TFile* file_DY_Sideband = new TFile(path+"MC.DYJets.root");
TFile* file_Diboson_Sideband = new TFile(path+"MC.Diboson.root");
TFile* file_QCD_Sideband = new TFile(path+"MC.QCD.root");
TFile* file_SingleTop_Sideband = new TFile(path+"MC.SingleTop.root");
TFile* file_WJets_Sideband = new TFile(path+"MC.WJets.root");
TFile* file_TTV_Sideband = new TFile(path+"MC.TTV.root");

vector<TH1D*> finalhistos;

// vector<TString> histnames = {"H_T_from350_all_filled_rebin", "N_jets", "Pt_lept1", "pt_jets", "N_pv", "E_Tmiss"};
vector<TString> histnames = {"ST_rebin4", "ST_rebinlimit", "STjets_rebin4", "STlep_rebin4", "M_LQ_rebinlimit"};
//vector<TString> histnames = {"H_T_from350_all_filled_rebin"};

for(unsigned int i=0; i<histnames.size(); i++){

  TH1D *h_HT_DATA, *h_HT_DATA_MCSubtracted, *h_HT_DATA_MCSubtracted_byhand, *h_HT_TTbar, *h_HT_DY, *h_HT_Diboson, *h_HT_QCD, *h_HT_SingleTop, *h_HT_WJets, *h_HT_TTV;

  TString histname = histnames[i];

  TString histpath = "finalselection_alpha/" + histname;
  cout << "histpath: " << histpath << endl;
  /*
  h_HT_DATA = (TH1D*)file_DATA_Sideband->Get("Sideband_weights_applied/H_T_from350_all_filled_rebin");
  h_HT_DATA_MCSubtracted = (TH1D*)h_HT_DATA->Clone("H_T_from350_all_filled_rebin");
  h_HT_DATA_MCSubtracted_byhand = (TH1D*)h_HT_DATA->Clone("h_HT_DATA_MCSubtracted_byhand");

  h_HT_TTbar = (TH1D*)file_TTbar_Sideband->Get("Sideband_weights_applied/H_T_from350_all_filled_rebin");
  h_HT_DY = (TH1D*)file_DY_Sideband->Get("Sideband_weights_applied/H_T_from350_all_filled_rebin");
  h_HT_Diboson = (TH1D*)file_Diboson_Sideband->Get("Sideband_weights_applied/H_T_from350_all_filled_rebin");
  h_HT_QCD = (TH1D*)file_QCD_Sideband->Get("Sideband_weights_applied/H_T_from350_all_filled_rebin");
  h_HT_SingleTop = (TH1D*)file_SingleTop_Sideband->Get("Sideband_weights_applied/H_T_from350_all_filled_rebin");
  h_HT_WJets = (TH1D*)file_WJets_Sideband->Get("Sideband_weights_applied/H_T_from350_all_filled_rebin");
  h_HT_TTV = (TH1D*)file_TTV_Sideband->Get("Sideband_weights_applied/H_T_from350_all_filled_rebin");
  */

  h_HT_DATA = (TH1D*)file_DATA_Sideband->Get(histpath);

  h_HT_DATA_MCSubtracted = (TH1D*)h_HT_DATA->Clone(histname);
  cout << "hier" << endl;

  h_HT_DATA_MCSubtracted_byhand = (TH1D*)h_HT_DATA->Clone("h_HT_DATA_MCSubtracted_byhand");
  cout << "hier" << endl;
  h_HT_TTbar = (TH1D*)file_TTbar_Sideband->Get(histpath);
  cout << "hier" << endl;
  h_HT_DY = (TH1D*)file_DY_Sideband->Get(histpath);
  cout << "hier" << endl;
  h_HT_Diboson = (TH1D*)file_Diboson_Sideband->Get(histpath);
  cout << "hier" << endl;
  h_HT_QCD = (TH1D*)file_QCD_Sideband->Get(histpath);
  cout << "hier" << endl;
  h_HT_SingleTop = (TH1D*)file_SingleTop_Sideband->Get(histpath);
  cout << "hier" << endl;
  h_HT_WJets = (TH1D*)file_WJets_Sideband->Get(histpath);
  cout << "hier" << endl;
  h_HT_TTV = (TH1D*)file_TTV_Sideband->Get(histpath);
  cout << "hier1223" << endl;
  //Add up all MCs and scale them to the integral of data - bin by bin
  TH1D* h_HT_MC_Sum = (TH1D*)h_HT_TTbar->Clone("h_HT_MC_Sum");
  cout << "hier1223" << endl;
  h_HT_MC_Sum->Add(h_HT_DY, 1);
  h_HT_MC_Sum->Add(h_HT_Diboson, 1);
  h_HT_MC_Sum->Add(h_HT_QCD, 1);
  h_HT_MC_Sum->Add(h_HT_SingleTop, 1);
  h_HT_MC_Sum->Add(h_HT_WJets, 1);
  h_HT_MC_Sum->Add(h_HT_TTV, 1);
  /*
  double MC_int = h_HT_MC_Sum->Integral();
  double DATA_int = h_HT_DATA->Integral();
  double ratio = MC_int/DATA_int;
  double SF = 1/ratio;
  */
  cout << "hier" << endl;
  //subtract minor MC backgrounds bin by bin to avoid empty bins
  for(int i=1; i<h_HT_DATA_MCSubtracted->GetNbinsX()+1; i++){
    double mc_bincontent = h_HT_MC_Sum->GetBinContent(i);
    double data_bincontent = h_HT_DATA->GetBinContent(i);
    double data_error = h_HT_DATA->GetBinError(i);
    double sf = data_bincontent/mc_bincontent;

    //scale minor backgrounds by sf and subtract these from data
    double dy_content = h_HT_DY->GetBinContent(i) * sf;
    double diboson_content = h_HT_Diboson->GetBinContent(i) * sf;
    double qcd_content = h_HT_QCD->GetBinContent(i) * sf;
    double singletop_content = h_HT_SingleTop->GetBinContent(i) * sf;
    double wjets_content = h_HT_WJets->GetBinContent(i) * sf;
    double ttv_content = h_HT_TTV->GetBinContent(i) * sf;
    double dy_error = h_HT_DY->GetBinError(i) * sf;
    double diboson_error = h_HT_Diboson->GetBinError(i) * sf;
    double qcd_error = h_HT_QCD->GetBinError(i) * sf;
    double singletop_error = h_HT_SingleTop->GetBinError(i) * sf;
    double wjets_error = h_HT_WJets->GetBinError(i) * sf;
    double ttv_error = h_HT_TTV->GetBinError(i) * sf;

    double new_bincontent;
    if(do_ttbardy) new_bincontent = data_bincontent - diboson_content - qcd_content - singletop_content - wjets_content - ttv_content;
    else new_bincontent = data_bincontent - dy_content - diboson_content - qcd_content - singletop_content - wjets_content - ttv_content;
    double new_error;
    if(do_ttbardy) new_error = sqrt(pow(data_error,2) + pow(diboson_error,2) + pow(qcd_error,2) + pow(singletop_error,2) + pow(wjets_error,2) + pow(ttv_error,2));
    else new_error = sqrt(pow(data_error,2) + pow(dy_error,2) + pow(diboson_error,2) + pow(qcd_error,2) + pow(singletop_error,2) + pow(wjets_error,2) + pow(ttv_error,2));


    if(do_ttbardy){
      cout << "In bin " << i << ", calculating: " << data_bincontent << " - " << diboson_content << " - " << qcd_content << " - " << singletop_content << " - " << wjets_content << " - " << ttv_content << " = " << new_bincontent << endl;
      cout << "SF in this bin was: " << sf << endl << endl;
    }
    else {
      cout << "In bin " << i << ", calculating: " << data_bincontent << " - " << dy_content << " - " << diboson_content << " - " << qcd_content << " - " << singletop_content << " - " << wjets_content << " - " << ttv_content << " = " << new_bincontent << endl;
      cout << "SF in this bin was: " << sf << endl << endl;
    }


    h_HT_DATA_MCSubtracted->SetBinContent(i,new_bincontent);
    h_HT_DATA_MCSubtracted->SetBinError(i,new_error);
  }

  /*
  //subtract MC backgrounds other than TTbar
  if(!do_ttbardy) h_HT_DATA_MCSubtracted->Add(h_HT_DY, -SF);
  h_HT_DATA_MCSubtracted->Add(h_HT_Diboson, -SF);
  h_HT_DATA_MCSubtracted->Add(h_HT_QCD, -SF);
  h_HT_DATA_MCSubtracted->Add(h_HT_SingleTop, -SF);
  h_HT_DATA_MCSubtracted->Add(h_HT_WJets, -SF);
  */

  for(int i=1; i<h_HT_DATA_MCSubtracted->GetNbinsX()+1; i++){
    if (h_HT_DATA->GetBinContent(i) == 0) {
      h_HT_DATA_MCSubtracted->SetBinContent(i,0);
      h_HT_DATA_MCSubtracted->SetBinError(i,0);
    }
    if(h_HT_DATA_MCSubtracted->GetBinContent(i) < 0){
      h_HT_DATA_MCSubtracted->SetBinContent(i,0);
      h_HT_DATA_MCSubtracted->SetBinError(i,1.5*h_HT_DATA_MCSubtracted->GetBinError(i));
    }
  }

  h_HT_DATA_MCSubtracted->SetLineColor(2);
  finalhistos.emplace_back(h_HT_DATA_MCSubtracted);

}

// same values
// if(do_systematics){
//   path_out += SystType + "_" + SystDirection+"/";
// }
// else path_out += "NOMINAL/";
if(do_ttbardy) path_out += "TTbarDYSideband.root";
else           path_out += "TTbarSideband.root";
cout <<"Out_Path: " << path_out << endl;


TFile* f_out = new TFile(path_out,"RECREATE");
f_out->mkdir("finalselection_mlqfalse/");
f_out->cd("finalselection_mlqfalse/");
for(unsigned int i=0; i<finalhistos.size(); i++){
  //path_out = "/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/35867fb_DibosonNLO/HTLepInvertedSideband/TTbarDYSideband.root";
  finalhistos[i]->Write();
}
f_out->Close();







//delete f_out;
/*
delete h_HT_WJets;
delete h_HT_SingleTop;
delete h_HT_QCD;
delete h_HT_Diboson;
delete h_HT_DY;
delete h_HT_TTbar;
delete h_HT_DATA_MCSubtracted_byhand;
delete h_HT_DATA_MCSubtracted;
delete h_HT_DATA;
*/
for(unsigned int i=0; i<finalhistos.size(); i++) delete finalhistos[i];
delete file_WJets_Sideband;
delete file_SingleTop_Sideband;
delete file_QCD_Sideband;
delete file_Diboson_Sideband;
delete file_DY_Sideband;
delete file_TTbar_Sideband;
delete file_DATA_Sideband;








}
