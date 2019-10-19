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


void FakeRateRunner::EleFakeRate(){

  cout << "Hello World from EleFakeRate!!" << endl;

  // TString pathele;
  // TString pathele_out;
  TString histpath_fake;
  TString histpath_real;
  TString histpath_match;
  TString histpath_before;

  // TString path_di = DibosonSFRunner::inpathDi;
  // TString path_dibtag = DibosonSFRunner::inpathDibtag;
  // TString path_out_di = DibosonSFRunner::outpathDi;
  // TString path_out_dibtag = DibosonSFRunner::outpathDibtag;
  TString path = pathele;
  TString path_out = pathele_out;






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


  histpath_fake = "FakeRate_FinalSelection/pt_matchingjets_rebin_fakeele";
  histpath_real = "FakeRate_FinalSelection/pt_matchingjets_rebin_realele";
  histpath_match = "FakeRate_FinalSelection/pt_matchingjets_rebin";
  histpath_before = "FakeRate_DibosonSF/pt_alljets_rebin";
  // histpath_real_data = "dummypath";
  // histptah_fake_data = "dummypath";

  TH1D *h_data_match, *h_data_before, *h_data_fake, *h_data_real;
  TH1D *h_ttbar_match, *h_ttbar_before, *h_ttbar_fake, *h_ttbar_real;
  TH1D *h_dy_match, *h_dy_before, *h_dy_fake, *h_dy_real;
  TH1D *h_diboson_match, *h_diboson_before, *h_diboson_fake, *h_diboson_real;
  TH1D *h_qcd_match, *h_qcd_before, *h_qcd_fake, *h_qcd_real;
  TH1D *h_singletop_match, *h_singletop_before, *h_singletop_fake, *h_singletop_real;
  TH1D *h_ttv_match, *h_ttv_before, *h_ttv_fake, *h_ttv_real;


  h_data_match = (TH1D*)file_DATA->Get(histpath_match);
  h_data_before = (TH1D*)file_DATA->Get(histpath_before);
  // h_data_fake = (TH1D*)file_DATA->Get(histpath_fake);
  // h_data_real = (TH1D*)file_DATA->Get(histpath_real);

  h_ttbar_match = (TH1D*)file_TTbar->Get(histpath_match);
  h_ttbar_before = (TH1D*)file_TTbar->Get(histpath_before);
  h_ttbar_fake = (TH1D*)file_TTbar->Get(histpath_fake);
  h_ttbar_real = (TH1D*)file_TTbar->Get(histpath_real);
  h_dy_match = (TH1D*)file_DY->Get(histpath_match);
  h_dy_before = (TH1D*)file_DY->Get(histpath_before);
  h_dy_fake = (TH1D*)file_DY->Get(histpath_fake);
  h_dy_real = (TH1D*)file_DY->Get(histpath_real);
  h_diboson_match = (TH1D*)file_Diboson->Get(histpath_match);
  h_diboson_before = (TH1D*)file_Diboson->Get(histpath_before);
  h_diboson_fake = (TH1D*)file_Diboson->Get(histpath_fake);
  h_diboson_real = (TH1D*)file_Diboson->Get(histpath_real);
  h_qcd_match = (TH1D*)file_QCD->Get(histpath_match);
  h_qcd_before = (TH1D*)file_QCD->Get(histpath_before);
  h_qcd_fake = (TH1D*)file_QCD->Get(histpath_fake);
  h_qcd_real = (TH1D*)file_QCD->Get(histpath_real);
  h_singletop_match = (TH1D*)file_SingleTop->Get(histpath_match);
  h_singletop_before = (TH1D*)file_SingleTop->Get(histpath_before);
  h_singletop_fake = (TH1D*)file_SingleTop->Get(histpath_fake);
  h_singletop_real = (TH1D*)file_SingleTop->Get(histpath_real);
  h_ttv_match = (TH1D*)file_TTV->Get(histpath_match);
  h_ttv_before = (TH1D*)file_TTV->Get(histpath_before);
  h_ttv_fake = (TH1D*)file_TTV->Get(histpath_fake);
  h_ttv_real = (TH1D*)file_TTV->Get(histpath_real);




  TH1D* h_total_match = (TH1D*)h_diboson_match->Clone("h_total_match");
  TH1D* h_total_before = (TH1D*)h_diboson_before->Clone("h_total_before");
  TH1D* h_total_fake = (TH1D*)h_diboson_fake->Clone("h_total_fake");
  TH1D* h_total_real = (TH1D*)h_diboson_real->Clone("h_total_real");
  h_total_match->Add(h_ttbar_match);
  h_total_match->Add(h_dy_match);
  h_total_match->Add(h_qcd_match);
  h_total_match->Add(h_singletop_match);
  h_total_match->Add(h_ttv_match);
  h_total_before->Add(h_ttbar_before);
  h_total_before->Add(h_dy_before);
  h_total_before->Add(h_qcd_before);
  h_total_before->Add(h_singletop_before);
  h_total_before->Add(h_ttv_before);
  h_total_fake->Add(h_ttbar_fake);
  h_total_fake->Add(h_dy_fake);
  h_total_fake->Add(h_qcd_fake);
  h_total_fake->Add(h_singletop_fake);
  h_total_fake->Add(h_ttv_fake);
  h_total_real->Add(h_ttbar_real);
  h_total_real->Add(h_dy_real);
  h_total_real->Add(h_qcd_real);
  h_total_real->Add(h_singletop_real);
  h_total_real->Add(h_ttv_real);



  TH1D* h_nomatch_before_mc = (TH1D*)h_total_before->Clone("h_nomatch_before_mc");
  h_nomatch_before_mc->Add(h_total_match,-1);
  TH1D* h_nomatch_before_data = (TH1D*)h_data_before->Clone("h_nomatch_before_data");
  h_nomatch_before_data->Add(h_data_match,-1);



  double SF_nomatch[3], err_nomatch[3];
  for(unsigned int i=0; i<h_total_before->GetNbinsX(); i++){
    SF_nomatch[i] = h_nomatch_before_data->GetBinContent(i+1) / h_nomatch_before_mc->GetBinContent(i+1);
    err_nomatch[i] = sqrt(pow(h_nomatch_before_data->GetBinError(i+1)/h_nomatch_before_mc->GetBinContent(i+1),2) + pow(h_nomatch_before_data->GetBinContent(i+1)/h_nomatch_before_mc->GetBinContent(i+1)/h_nomatch_before_mc->GetBinContent(i+1)*h_nomatch_before_mc->GetBinError(i+1),2));
    cout << "SF nomatch: " << SF_nomatch[i] << endl;
  }


TH1D* h_nominator_data, *h_denominator_data, *h_nominator_mc, *h_denominator_mc;
h_nominator_data  = (TH1D*)h_data_match->Clone("h_nominator_data");
h_nominator_data->Add(h_total_real,-1);
h_denominator_data = (TH1D*)h_data_before->Clone("h_denominator_data");
h_nominator_mc = (TH1D*)h_total_fake->Clone("h_nominator_mc");
h_denominator_mc = (TH1D*)h_total_before->Clone("h_denominator_mc");

TGraphAsymmErrors* gr_fakerate_mc, *gr_fakerate_data;
gr_fakerate_mc = new TGraphAsymmErrors(h_nominator_mc, h_denominator_mc);
gr_fakerate_data = new TGraphAsymmErrors(h_nominator_data, h_denominator_data);

unique_ptr<TCanvas> c1, c2;
c1.reset(new TCanvas());
gr_fakerate_mc->Draw();

c2.reset(new TCanvas());
gr_fakerate_data->Draw();





//4. Calculate Data/MC
double binsx[3] = {50,150,500};
double ex[3] = {50,50,300};
double binsy[3], eyl[3], eyh[3];
for(int i=0; i<3; i++){
  double xmc, ymc, xdata, ydata, ydatalow, ydatahigh, ymclow, ymchigh;
  gr_fakerate_mc->GetPoint(i,xmc,ymc);
  gr_fakerate_data->GetPoint(i,xdata,ydata);
  ymclow = gr_fakerate_mc->GetErrorYlow(i);
  ymchigh = gr_fakerate_mc->GetErrorYhigh(i);
  ydatalow = gr_fakerate_data->GetErrorYlow(i);
  ydatahigh = gr_fakerate_data->GetErrorYhigh(i);

  binsy[i] = ydata/ymc;
  cout <<"Final SF: " << binsy[i] << endl;
  eyh[i] = sqrt(pow(ydatahigh/ymc,2) + pow(ydata*ymclow/ymc/ymc,2));
  eyl[i] = sqrt(pow(ydatalow/ymc,2) + pow(ydata*ymchigh/ymc/ymc,2));

}


TGraphAsymmErrors* SF;
SF = new TGraphAsymmErrors(3, binsx, binsy, ex, ex, eyl, eyh);
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
gr_fakerate_data->GetHistogram()->SetMaximum(maxall*1.7);
gr_fakerate_data->GetHistogram()->SetMinimum(min(0.000005,min(TMath::MinElement(gr_fakerate_data->GetN(), gr_fakerate_data->GetY()),TMath::MinElement(gr_fakerate_mc->GetN(), gr_fakerate_mc->GetY())) - 0.5*min(TMath::MinElement(gr_fakerate_data->GetN(), gr_fakerate_data->GetY()),TMath::MinElement(gr_fakerate_mc->GetN(), gr_fakerate_mc->GetY()))));
gr_fakerate_data->GetXaxis()->SetLimits(0,800);
gr_fakerate_data->GetYaxis()->SetTitle("Misidentification rate");
gr_fakerate_data->SetMarkerStyle(20);
gr_fakerate_data->SetMarkerSize(0.9);
TGaxis::SetMaxDigits(3);

gr_fakerate_data->Draw("AP");
gr_fakerate_mc->Draw("P SAME");

TLegend* leg_top = new TLegend(0.7,0.7,0.9,0.9);
leg_top->AddEntry(gr_fakerate_data,"Data","lep");
leg_top->AddEntry(gr_fakerate_mc,"Simulation","le");
leg_top->SetBorderSize(0);
leg_top->SetFillStyle(0);
leg_top->Draw();

TString cmstext = "CMS";
TLatex *text2 = new TLatex(3.5, 24, cmstext);
text2->SetNDC();
text2->SetTextAlign(13);
text2->SetX(0.24);
text2->SetTextFont(62);
text2->SetTextSize(0.08);
text2->SetY(0.87);
text2->Draw();

TString preltext = "Preliminary";
TLatex *text3 = new TLatex(3.5, 24, preltext);
text3->SetNDC();
text3->SetTextAlign(13);
text3->SetX(0.24);
text3->SetTextFont(52);
text3->SetTextSize(0.06);
text3->SetY(0.78);
text3->Draw();

p_bot->cd();

SF->SetLineColor(kBlack);
SF->GetYaxis()->SetRangeUser(0,2.9);
SF->GetXaxis()->SetLimits(0,800);
SF->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
SF->GetYaxis()->SetTitle("Data / sim");
SF->SetMarkerStyle(20);
SF->SetMarkerSize(0.9);
SF->Draw("AP");
TLine* l_unity = new TLine(0, 1, 800, 1);
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
f_out.reset(new TFile(path_out + "ElectronFakeRateSF_DibosonNLO.root","RECREATE"));
SF->Write();
f_out->Close();


delete p_top;
delete p_bot;
delete c_final;
delete SF;
delete gr_fakerate_mc;
delete gr_fakerate_data;
// delete h_all_before_data;
// delete h_matched_after_data;

delete h_data_match;
delete h_data_before;




delete h_ttbar_match;
delete h_ttbar_before;
delete h_ttbar_fake;
delete h_ttbar_real;
delete h_dy_match;
delete h_dy_before;
delete h_dy_fake;
delete h_dy_real;
delete h_diboson_match;
delete h_diboson_before;
delete h_diboson_fake;
delete h_diboson_real;
delete h_qcd_match;
delete h_qcd_before;
delete h_qcd_fake;
delete h_qcd_real;
delete h_singletop_match;
delete h_singletop_before;
delete h_singletop_fake;
delete h_singletop_real;
delete h_ttv_match;
delete h_ttv_before;
delete h_ttv_fake;
delete h_ttv_real;





// delete h_all_before_dy;
// delete h_real_after_dy;
// delete h_fake_after_dy;
// delete h_all_before_db;
// delete h_real_after_db;
// delete h_fake_after_db;
// delete h_all_before_ttbar;
// delete h_real_after_ttbar;
// delete h_fake_after_ttbar;
// delete h_all_before_wjets;
// delete h_real_after_wjets;
// delete h_fake_after_wjets;
// delete h_all_before_singletop;
// delete h_real_after_singletop;
// delete h_fake_after_singletop;
// delete h_all_before_qcd;
// delete h_real_after_qcd;
// delete h_fake_after_qcd;
// delete h_all_before_mc;
// delete h_real_after_mc;
// delete h_fake_after_mc;
delete h_nominator_data;
delete h_denominator_data;
delete h_nominator_mc;
delete h_denominator_mc;





}
