#include "../include/FakeRateRunner.h"
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




#include <iostream>
#include <TMath.h>
#include <TLine.h>
#include <TLatex.h>

#include <TAxis.h>

#include "../include/cosmetics.h"
#include "../include/FakeRateRunner.h"

#include <TLegend.h>


using namespace std;

void FakeRateRunner::MuFakeRate(){

  cout << "Hello World from MuFakeRate!!" << endl;

  // TString pathele;
  // TString pathele_out;
  TString histpath_h_all_before_data;
  TString histpth_h_matched_after_data;
  TString histpath_fake;
  TString histpath_real;
  TString histpath_match;
  TString histpath_before;

  // TString path_di = DibosonSFRunner::inpathDi;
  // TString path_dibtag = DibosonSFRunner::inpathDibtag;
  // TString path_out_di = DibosonSFRunner::outpathDi;
  // TString path_out_dibtag = DibosonSFRunner::outpathDibtag;
  TString path = pathmu;
  TString path_out = pathmu_out;






  path += "uhh2.AnalysisModuleRunner.";
  //cout << "Input path: " << path+"DATA.DATA.root" << endl;

  TFile* file_DATA = new TFile(path+"DATA.DATA.root");
  TFile* file_TTbar = new TFile(path+"MC.TTbar.root");
  TFile* file_DY = new TFile(path+"MC.DYJets.root");
  TFile* file_Diboson = new TFile(path+"MC.Diboson.root");
  TFile* file_QCD = new TFile(path+"MC.QCD.root");
  TFile* file_SingleTop = new TFile(path+"MC.SingleTop.root");
  // TFile* file_WJets = new TFile(path+"MC.WJets.root");
  TFile* file_TTV = new TFile(path+"MC.TTV.root");

  histpath_h_all_before_data = "FakeRate_DibosonSF/Int_events_denom";
  histpth_h_matched_after_data = "FakeRate_FinalSelection/Int_events_1mu";
  histpath_fake = "FakeRate_FinalSelection/Int_events_1mu_fake";
  histpath_real = "FakeRate_FinalSelection/Int_events_1mu_real";
  histpath_match = "FakeRate_DibosonSF/Int_events_1mu";
  histpath_before = "FakeRate_DibosonSF/Int_events_denom";
  // histpath_real_data = "dummypath";
  // histptah_fake_data = "dummypath";

  TH1D *h_all_before_data, *h_matched_after_data, *h_all_before_dy, *h_1mu_dy, *h_real_after_dy, *h_fake_after_dy;
  TH1D *h_all_before_db, *h_1mu_db, *h_real_after_db, *h_fake_after_db;
  TH1D *h_all_before_ttbar, *h_1mu_ttbar, *h_real_after_ttbar, *h_fake_after_ttbar;
  TH1D *h_all_before_singletop, *h_1mu_singletop, *h_real_after_singletop, *h_fake_after_singletop;
  TH1D *h_all_before_qcd,  *h_1mu_qcd, *h_real_after_qcd, *h_fake_after_qcd;
  TH1D *h_all_before_ttv, *h_1mu_ttv, *h_real_after_ttv, *h_fake_after_ttv;


  h_all_before_data = (TH1D*)file_DATA->Get(histpath_h_all_before_data);
  h_matched_after_data = (TH1D*)file_DATA->Get(histpth_h_matched_after_data);

  h_all_before_dy = (TH1D*)file_DY->Get(histpath_before);
  h_1mu_dy = (TH1D*)file_DY->Get(histpath_match);
  h_real_after_dy = (TH1D*)file_DY->Get(histpath_real);
  h_fake_after_dy = (TH1D*)file_DY->Get(histpath_fake);
  h_all_before_db = (TH1D*)file_Diboson->Get(histpath_before);
  h_1mu_db = (TH1D*)file_Diboson->Get(histpath_match);
  h_real_after_db = (TH1D*)file_Diboson->Get(histpath_real);
  h_fake_after_db = (TH1D*)file_Diboson->Get(histpath_fake);
  h_all_before_ttbar = (TH1D*)file_TTbar->Get(histpath_before);
  h_1mu_ttbar = (TH1D*)file_TTbar->Get(histpath_match);
  h_real_after_ttbar = (TH1D*)file_TTbar->Get(histpath_real);
  h_fake_after_ttbar = (TH1D*)file_TTbar->Get(histpath_fake);
  // h_all_before_wjets = (TH1D*)in_wjets->Get(histpath_before);
  // h_1mu_wjets = (TH1D*)in_wjets->Get(histpath_match);
  // h_real_after_wjets = (TH1D*)in_wjets->Get(histpath_real);
  // h_fake_after_wjets = (TH1D*)in_wjets->Get(histpath_fake);
  h_all_before_singletop = (TH1D*)file_SingleTop->Get(histpath_before);
  h_1mu_singletop = (TH1D*)file_SingleTop->Get(histpath_match);
  h_real_after_singletop = (TH1D*)file_SingleTop->Get(histpath_real);
  h_fake_after_singletop = (TH1D*)file_SingleTop->Get(histpath_fake);
  h_all_before_qcd = (TH1D*)file_QCD->Get(histpath_before);
  h_1mu_qcd = (TH1D*)file_QCD->Get(histpath_match);
  h_real_after_qcd = (TH1D*)file_QCD->Get(histpath_real);
  h_fake_after_qcd = (TH1D*)file_QCD->Get(histpath_fake);
  h_all_before_ttv = (TH1D*)file_TTV->Get(histpath_before);
  h_1mu_ttv = (TH1D*)file_TTV->Get(histpath_match);
  h_real_after_ttv = (TH1D*)file_TTV->Get(histpath_real);
  h_fake_after_ttv = (TH1D*)file_TTV->Get(histpath_fake);


  //1. add all histograms from MC
  TH1D* h_all_before_mc, *h_real_after_mc, *h_fake_after_mc, *h_1mu_mc;
  h_all_before_mc = (TH1D*)h_all_before_dy->Clone("h_all_before_mc");
  h_1mu_mc = (TH1D*)h_1mu_dy->Clone("h_1mu_mc");
  h_real_after_mc = (TH1D*)h_real_after_dy->Clone("h_real_after_mc");
  h_fake_after_mc = (TH1D*)h_fake_after_dy->Clone("h_fake_after_mc");
  h_all_before_mc->Add(h_all_before_db);
  h_1mu_mc->Add(h_1mu_db);
  h_real_after_mc->Add(h_real_after_db);
  h_fake_after_mc->Add(h_fake_after_db);
  h_all_before_mc->Add(h_all_before_ttbar);
  h_1mu_mc->Add(h_1mu_ttbar);
  h_real_after_mc->Add(h_real_after_ttbar);
  h_fake_after_mc->Add(h_fake_after_ttbar);
  // h_all_before_mc->Add(h_all_before_wjets);
  // h_1mu_mc->Add(h_1mu_wjets);
  // h_real_after_mc->Add(h_real_after_wjets);
  // h_fake_after_mc->Add(h_fake_after_wjets);
  h_all_before_mc->Add(h_all_before_singletop);
  h_1mu_mc->Add(h_1mu_singletop);
  h_real_after_mc->Add(h_real_after_singletop);
  h_fake_after_mc->Add(h_fake_after_singletop);
  h_all_before_mc->Add(h_all_before_qcd);
  h_1mu_mc->Add(h_1mu_qcd);
  h_real_after_mc->Add(h_real_after_qcd);
  h_fake_after_mc->Add(h_fake_after_qcd);
  h_all_before_mc->Add(h_all_before_ttv);
  h_1mu_mc->Add(h_1mu_ttv);
  h_real_after_mc->Add(h_real_after_ttv);
  h_fake_after_mc->Add(h_fake_after_ttv);

  //1.1 Calculate DATA/MC SF for 0mu region to get rid of other effects on data/mc ratio
  TH1D* h_0mu_data, *h_0mu_mc;
  h_0mu_mc = (TH1D*)h_all_before_mc->Clone("h_0mu_mc");
  h_0mu_data = (TH1D*)h_all_before_data->Clone("h_0mu_data");
  h_0mu_mc->Add(h_1mu_mc,-1);
  h_0mu_data->Add(h_matched_after_data,-1);

  double SF_0mu = h_0mu_data->GetBinContent(1) / h_0mu_mc->GetBinContent(1);
  double err_0mu = sqrt(pow(h_0mu_data->GetBinError(1)/h_0mu_mc->GetBinContent(1),2) + pow(h_0mu_data->GetBinContent(1)/h_0mu_mc->GetBinContent(1)/h_0mu_mc->GetBinContent(1)*h_0mu_mc->GetBinError(1),2));

  //
  // double new_err_before = sqrt(pow(SF_0mu*h_all_before_mc->GetBinError(1),2) + pow(err_0mu*h_all_before_mc->GetBinContent(1),2));
  // double new_err_real_after = sqrt(pow(SF_0mu*h_real_after_mc->GetBinError(1),2) + pow(err_0mu*h_real_after_mc->GetBinContent(1),2));
  // double new_err_fake_after = sqrt(pow(SF_0mu*h_fake_after_mc->GetBinError(1),2) + pow(err_0mu*h_fake_after_mc->GetBinContent(1),2));
  // h_all_before_mc->SetBinContent(1,h_all_before_mc->GetBinContent(1) * SF_0mu);
  // h_real_after_mc->SetBinContent(1,h_real_after_mc->GetBinContent(1) * SF_0mu);
  // h_fake_after_mc->SetBinContent(1,h_fake_after_mc->GetBinContent(1) * SF_0mu);
  // h_all_before_mc->SetBinError(1,new_err_before);
  // h_real_after_mc->SetBinError(1,new_err_real_after);
  // h_fake_after_mc->SetBinError(1,new_err_fake_after);
  // cout << "scaling MC with " << SF_0mu << " globally" << endl;

  //2.  Calculate (de)nom for data and mc

  TH1D* h_nominator_data, *h_denominator_data, *h_nominator_mc, *h_denominator_mc;
  h_nominator_data  = (TH1D*)h_matched_after_data->Clone("h_nominator_data");
  cout << "total nom, data: " << h_nominator_data->GetBinContent(1) << endl;
  h_nominator_data->Add(h_real_after_mc,-1);
  cout << "real, mc: " << h_real_after_mc->GetBinContent(1) << endl;
  h_denominator_data = (TH1D*)h_all_before_data->Clone("h_denominator_data");
  h_nominator_mc = (TH1D*)h_fake_after_mc->Clone("h_nominator_mc");
  h_denominator_mc = (TH1D*)h_all_before_mc->Clone("h_denominator_mc");

  cout << "nom, data: " << h_nominator_data->GetBinContent(1) << " +- " << h_nominator_data->GetBinError(1) << endl;
  cout << "denom, data: " << h_denominator_data->GetBinContent(1) << " +- " << h_denominator_data->GetBinError(1) << endl;
  cout << "nom, mc: " << h_nominator_mc->GetBinContent(1) << " +- " << h_nominator_mc->GetBinError(1) << endl;
  cout << "denom, mc: " << h_denominator_mc->GetBinContent(1) << " +- " << h_denominator_mc->GetBinError(1) << endl;


  //3. Calculate fake rates for data and mc separately
  TGraphAsymmErrors* gr_fakerate_mc, *gr_fakerate_data;
  gr_fakerate_mc = new TGraphAsymmErrors(h_nominator_mc, h_denominator_mc);
  gr_fakerate_data = new TGraphAsymmErrors(h_nominator_data, h_denominator_data);

  unique_ptr<TCanvas> c1, c2;
  c1.reset(new TCanvas());
  gr_fakerate_mc->Draw();

  c2.reset(new TCanvas());
  gr_fakerate_data->Draw();


  //4. Calculate Data/MC
  double binsx[1] = {0};
  double ex[1] = {0.5};
  double binsy[1], eyl[1], eyh[1];
  for(int i=0; i<1; i++){
    double xmc, ymc, xdata, ydata, ydatalow, ydatahigh, ymclow, ymchigh;
    gr_fakerate_mc->GetPoint(i,xmc,ymc);
    gr_fakerate_data->GetPoint(i,xdata,ydata);
    ymclow = gr_fakerate_mc->GetErrorYlow(i);
    ymchigh = gr_fakerate_mc->GetErrorYhigh(i);
    ydatalow = gr_fakerate_data->GetErrorYlow(i);
    ydatahigh = gr_fakerate_data->GetErrorYhigh(i);

    cout << "DATA: " << ydata << " + " << ydatahigh << " - " << ydatalow << endl;
    cout << "MC: " << ymc << " + " << ymchigh << " - " << ymclow << endl;

    binsy[i] = ydata/ymc;
    eyh[i] = sqrt(pow(ydatahigh/ymc,2) + pow(ydata*ymclow/ymc/ymc,2));
    eyl[i] = sqrt(pow(ydatalow/ymc,2) + pow(ydata*ymchigh/ymc/ymc,2));
    cout << "RATIO: " << binsy[i] << " + " << eyh[i] << " - " << eyl[i] << endl;

  }


  TGraphAsymmErrors* SF;
  SF = new TGraphAsymmErrors(1, binsx, binsy, ex, ex, eyl, eyh);
  unique_ptr<TCanvas> c3;
  c3.reset(new TCanvas());
  SF->Draw();


  TPad* p_top = SetupRatioPadTop();
  TPad* p_bot = SetupRatioPad();
  TCanvas* c_final = new TCanvas("c_final", "Ratio plot", 400, 400);
  p_top->Draw();
  p_bot->Draw();
  p_top->cd();
  HistCosmetics(gr_fakerate_mc,false);
  HistCosmetics(gr_fakerate_data,false);
  HistCosmetics(SF,true);

  gr_fakerate_mc->SetLineColor(kRed);
  gr_fakerate_data->SetLineColor(kBlack);
  gr_fakerate_data->SetMarkerColor(kBlack);
  gr_fakerate_mc->SetMarkerColor(kRed);
  double max1 = TMath::MaxElement(gr_fakerate_data->GetN(), gr_fakerate_data->GetY());
  double max2 = TMath::MaxElement(gr_fakerate_mc->GetN(), gr_fakerate_mc->GetY());
  double maxall  = max(max1,max2);
  gr_fakerate_data->GetHistogram()->SetMaximum(maxall*1.5);
  // gr_fakerate_data->GetHistogram()->SetMinimum(min(0.000002,min(TMath::MinElement(gr_fakerate_data->GetN(), gr_fakerate_data->GetY()),TMath::MinElement(gr_fakerate_mc->GetN(), gr_fakerate_mc->GetY())) - 0.5*min(TMath::MinElement(gr_fakerate_data->GetN(), gr_fakerate_data->GetY()),TMath::MinElement(gr_fakerate_mc->GetN(), gr_fakerate_mc->GetY()))));
  gr_fakerate_data->GetXaxis()->SetLimits(-0.5,0.5);
  gr_fakerate_data->GetHistogram()->SetMinimum(-0.00005);
  // gr_fakerate_data->GetXaxis()->SetLimits(-0.000008, 0.000008);
  gr_fakerate_data->GetYaxis()->SetTitle("Misidentification rate");
  gr_fakerate_data->SetMarkerStyle(20);
  gr_fakerate_data->SetMarkerSize(0.9);

  gr_fakerate_data->Draw("AP");
  gr_fakerate_mc->Draw("P SAME");

  TString cmstext = "CMS";
  TLatex *text2 = new TLatex(3.5, 24, cmstext);
  text2->SetNDC();
  text2->SetTextAlign(13);
  text2->SetX(0.24);
  text2->SetTextFont(62);
  text2->SetTextSize(0.08);
  text2->SetY(0.89);
  text2->Draw();

  TString preltext = "Preliminary";
  TLatex *text3 = new TLatex(3.5, 24, preltext);
  text3->SetNDC();
  text3->SetTextAlign(13);
  text3->SetX(0.24);
  text3->SetTextFont(52);
  text3->SetTextSize(0.06);
  text3->SetY(0.80);
  text3->Draw();

  TLegend* leg_top = new TLegend(0.7,0.7,0.9,0.9);
  leg_top->AddEntry(gr_fakerate_data,"Data","lep");
  leg_top->AddEntry(gr_fakerate_mc,"MC","le");
  leg_top->SetBorderSize(0);
  leg_top->SetFillStyle(0);
  leg_top->Draw();

  p_bot->cd();

  SF->SetLineColor(kBlack);
  SF->GetYaxis()->SetRangeUser(0,5.4);
  SF->GetXaxis()->SetLimits(-0.5,0.5);
  SF->GetXaxis()->SetTitle("");
  SF->GetYaxis()->SetTitle("Data / sim");
  SF->SetMarkerStyle(20);
  SF->SetMarkerSize(0.9);
  SF->Draw("AP");
  TLine* l_unity = new TLine(-0.5, 1, 0.5, 1);
  l_unity->SetLineColor(kRed);
  l_unity->SetLineWidth(2);
  l_unity->SetLineStyle(2);
  l_unity->Draw("SAME");


  c1->SaveAs(path_out + "FakeRateMC_DibosonNLO.eps");
  c2->SaveAs(path_out + "FakeRateDATA_DibosonNLO.eps");
  c3->SaveAs(path_out + "ScaleFactors_DibosonNLO.eps");
  c_final->SaveAs(path_out + "RatioPlot_DibosonNLO.eps");


  SF->SetName("ScaleFactors");
  unique_ptr<TFile> f_out;
  f_out.reset(new TFile(path_out + "MuonFakeRateSF_DibosonNLO.root","RECREATE"));
  SF->Write();
  f_out->Close();


  delete p_top;
  delete p_bot;
  delete c_final;
  delete SF;
  delete gr_fakerate_mc;
  delete gr_fakerate_data;
  delete h_all_before_data;
  delete h_matched_after_data;
  delete h_all_before_dy;
  delete h_real_after_dy;
  delete h_fake_after_dy;
  delete h_all_before_db;
  delete h_real_after_db;
  delete h_fake_after_db;
  delete h_all_before_ttbar;
  delete h_real_after_ttbar;
  delete h_fake_after_ttbar;
  // delete h_all_before_wjets;
  // delete h_real_after_wjets;
  // delete h_fake_after_wjets;
  delete h_all_before_singletop;
  delete h_real_after_singletop;
  delete h_fake_after_singletop;
  delete h_all_before_qcd;
  delete h_real_after_qcd;
  delete h_fake_after_qcd;
  delete h_all_before_mc;
  delete h_real_after_mc;
  delete h_fake_after_mc;
  delete h_nominator_data;
  delete h_denominator_data;
  delete h_nominator_mc;
  delete h_denominator_mc;







}
