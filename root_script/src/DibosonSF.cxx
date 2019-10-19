#include "../include/DibosonSFRunner.h"
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

// Define function

void DibosonSFRunner::DibosonSF(){

  bool btag = true;


  TString path;
  TString path_out;
  TString path_di = DibosonSFRunner::inpathDi;
  TString path_dibtag = DibosonSFRunner::inpathDibtag;
  TString path_out_di = DibosonSFRunner::outpathDi;
  TString path_out_dibtag = DibosonSFRunner::outpathDibtag;


  if(btag == true){
    path = path_dibtag;
    path_out = path_out_dibtag;
  }
  // else{
  //   path = path_di;
  //   path_out = path_out_dibtag;
  // }



  path += "uhh2.AnalysisModuleRunner.";
  //cout << "Input path: " << path+"DATA.DATA.root" << endl;

  TFile* file_DATA = new TFile(path+"DATA.DATA.root");
  TFile* file_TTbar = new TFile(path+"MC.TTbar.root");
  TFile* file_DY = new TFile(path+"MC.DYJets.root");
  TFile* file_Diboson = new TFile(path+"MC.Diboson.root");

  // TFile* file_Diboson = new TFile(path+"MC.DibosonNLO.root");

  TFile* file_QCD = new TFile(path+"MC.QCD.root");
  TFile* file_SingleTop = new TFile(path+"MC.SingleTop.root");
  // TFile* file_WJets = new TFile(path+"MC.WJets.root");
  TFile* file_TTV = new TFile(path+"MC.TTV.root");

  TH1D *h_data, *h_ttbar, *h_dy, *h_diboson, *h_singletop, *h_ttv, *h_qcd; // *h_WJets,  *h_qcd,

  // TString histpath = "FinalSelection/N_bJets_loose";
  TString histpath = "FinalSelection/sum_event_weights";


  h_data =  (TH1D*)file_DATA->Get(histpath);
  h_ttbar = (TH1D*)file_TTbar->Get(histpath);
  h_dy = (TH1D*)file_DY->Get(histpath);
  h_diboson = (TH1D*)file_Diboson->Get(histpath);
  h_qcd = (TH1D*)file_QCD->Get(histpath);
  h_singletop = (TH1D*)file_SingleTop->Get(histpath);
  // h_wjets = (TH1D*)file_WJets->Get(histpath);
  h_ttv = (TH1D*)file_TTV->Get(histpath);

  TH1D* h_total_before = (TH1D*)h_diboson->Clone("h_total_before");
  h_total_before->Add(h_dy);
  h_total_before->Add(h_singletop);
  h_total_before->Add(h_ttbar);
  h_total_before->Add(h_ttv);
  // h_total_before->Add(h_wjets)
  h_total_before->Add(h_qcd);

  double content_data[3], content_nondb[3], content_remainingdata[3], err_remainingdata[3], content_db[3], err_db[3], sf[3], err_sf[3];
  for(int i=0; i<3; i++){
    content_data[i] = h_data->GetBinContent(i+1);
    content_nondb[i] = h_dy->GetBinContent(i+1) + h_singletop->GetBinContent(i+1) + h_ttbar->GetBinContent(i+1) + h_ttv->GetBinContent(i+1);
    content_remainingdata[i] = content_data[i] - content_nondb[i];
    err_remainingdata[i] = sqrt(pow(h_data->GetBinError(i+1),2) + pow(h_dy->GetBinError(i+1),2) + pow(h_singletop->GetBinError(i+1),2) + pow(h_ttbar->GetBinError(i+1),2) + pow(h_ttv->GetBinError(i+1),2));
    content_db[i] = h_diboson->GetBinContent(i+1);
    // cout << content_db[i] << endl;
    err_db[i] = h_diboson->GetBinError(i+1);
  }

  sf[0] = content_remainingdata[0]/content_db[0];
  err_sf[0] = sqrt(pow(err_remainingdata[0]/content_db[0],2) + pow(content_data[0]/content_db[0]/content_db[0]*err_db[0],2));



  cout << "----- Content Diboson before scaling -----" << endl;
  cout << "0 BTag: " << content_db[0] << endl;
  cout << "1 BTag: " << content_db[1] << endl;
  cout << "2 BTag: " << content_db[2] << endl << endl;

  cout << "----- Content non-Diboson before scaling -----" << endl;
  cout << "0 BTag: " << content_nondb[0] << endl;
  cout << "1 BTag: " << content_nondb[1] << endl;
  cout << "2 BTag: " << content_nondb[2] << endl << endl;

  cout << "----- Content total DATA before scaling -----" << endl;
  cout << "0 BTag: " << content_data[0] << endl;
  cout << "1 BTag: " << content_data[1] << endl;
  cout << "2 BTag: " << content_data[2] << endl << endl;

  cout << "----- Content remaining DATA before scaling -----" << endl;
  cout << "0 BTag: " << content_remainingdata[0] << endl;
  cout << "1 BTag: " << content_remainingdata[1] << endl;
  cout << "2 BTag: " << content_remainingdata[2] << endl << endl;


  h_diboson->Scale(sf[0]);

  for(int i=0; i<3; i++){
    content_db[i] = h_diboson->GetBinContent(i+1);
    err_db[i] = h_diboson->GetBinError(i+1);
    if(i>0){
      sf[i] = content_remainingdata[i]/content_db[i];
      err_sf[i] = sqrt(pow(err_remainingdata[i]/content_db[i],2) + pow(content_data[i]/content_db[i]/content_db[i]*err_db[i],2));
    }
  }


  cout << "sf 0tags: " << sf[0] << endl;
  cout << "sf 1tag : " << sf[1] << endl;
  cout << "sf 2tags: " << sf[2] << endl << endl;


  cout << "----- Content MC after scaling -----" << endl;
  cout << "0 BTag: " << content_db[0]+content_nondb[0] << endl;
  cout << "1 BTag: " << content_db[1]*sf[1]+content_nondb[1] << endl;
  cout << "2 BTag: " << content_db[2]*sf[2]+content_nondb[2] << endl << endl;

  cout << "----- Content DATA after scaling -----" << endl;
  cout << "0 BTag: " << content_data[0] << endl;
  cout << "1 BTag: " << content_data[1] << endl;
  cout << "2 BTag: " << content_data[2] << endl;



  TH1D* Diboson_XSec_SF = new TH1D("Diboson_XSec_SF", "Diboson XSec SF;;SF", 1, -0.5,0.5);
  Diboson_XSec_SF->SetBinContent(1,sf[0]);
  Diboson_XSec_SF->SetBinError(1,err_sf[0]);



  unique_ptr<TFile> out;
  out.reset(new TFile(path_out + "DibosonSF.root", "RECREATE"));
  Diboson_XSec_SF->Write();
  // Diboson_BTag_SF->Write();
  out->Close();


  // delete Diboson_BTag_SF;
  delete Diboson_XSec_SF;
  delete h_total_before;
  delete h_data;
  delete h_ttbar;
  delete h_singletop;
  delete h_dy;
  delete h_diboson;

}
