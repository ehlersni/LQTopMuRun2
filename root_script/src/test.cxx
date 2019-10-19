// #include "../include/Tools.h"
// #include <TH1D.h>
// #include <TFile.h>
// #include <TGraphAsymmErrors.h>
// #include <TCanvas.h>
// #include <iostream>
//
// using namespace std;
//
// void FakeRateTool::CalculateDibosonSF(){
//
//   unique_ptr<TFile> f_diboson, f_dy, f_singletop, f_ttbar, f_data, f_ttv;
//   f_diboson.reset(new TFile(FakeRateTool::base_path_Diboson + "uhh2.AnalysisModuleRunner.MC.Diboson.root","READ"));
//   f_dy.reset(new TFile(FakeRateTool::base_path_Diboson + "uhh2.AnalysisModuleRunner.MC.DYJets.root","READ"));
//   f_singletop.reset(new TFile(FakeRateTool::base_path_Diboson + "uhh2.AnalysisModuleRunner.MC.SingleTop.root","READ"));
//   f_ttbar.reset(new TFile(FakeRateTool::base_path_Diboson + "uhh2.AnalysisModuleRunner.MC.TTbar.root","READ"));
//   f_ttv.reset(new TFile(FakeRateTool::base_path_Diboson + "uhh2.AnalysisModuleRunner.MC.TTV.root","READ"));
//   f_data.reset(new TFile(FakeRateTool::base_path_Diboson + "uhh2.AnalysisModuleRunner.DATA.DATA.root","READ"));
//
//   TH1D* h_diboson = (TH1D*)f_diboson->Get("FinalSelection/N_bJets_loose");
//   TH1D* h_dy = (TH1D*)f_dy->Get("FinalSelection/N_bJets_loose");
//   TH1D* h_singletop = (TH1D*)f_singletop->Get("FinalSelection/N_bJets_loose");
//   TH1D* h_ttbar = (TH1D*)f_ttbar->Get("FinalSelection/N_bJets_loose");
//   TH1D* h_ttv = (TH1D*)f_ttv->Get("FinalSelection/N_bJets_loose");
//   TH1D* h_data = (TH1D*)f_data->Get("FinalSelection/N_bJets_loose");
//
//   //1. Sum up all mcs
//   TH1D* h_total_before = (TH1D*)h_diboson->Clone("h_total_before");
//   h_total_before->Add(h_dy);
//   h_total_before->Add(h_singletop);
//   h_total_before->Add(h_ttbar);
//   h_total_before->Add(h_ttv);
//
//   //2. calculate x-sec sf from 0-btag bin
//   //subtract non-DB backgrounds from data, then scale DB to data
//   //double content_0tag_mc = h_total_before->GetBinContent(1);
//
//
//   double content_data[3], content_nondb[3], content_remainingdata[3], err_remainingdata[3], content_db[3], err_db[3], sf[3], err_sf[3];
//   for(int i=0; i<3; i++){
//     content_data[i] = h_data->GetBinContent(i+1);
//     content_nondb[i] = h_dy->GetBinContent(i+1) + h_singletop->GetBinContent(i+1) + h_ttbar->GetBinContent(i+1) + h_ttv->GetBinContent(i+1);
//     content_remainingdata[i] = content_data[i] - content_nondb[i];
//     err_remainingdata[i] = sqrt(pow(h_data->GetBinError(i+1),2) + pow(h_dy->GetBinError(i+1),2) + pow(h_singletop->GetBinError(i+1),2) + pow(h_ttbar->GetBinError(i+1),2) + pow(h_ttv->GetBinError(i+1),2));
//     content_db[i] = h_diboson->GetBinContent(i+1);
//     err_db[i] = h_diboson->GetBinError(i+1);
//   }
//
//   sf[0] = content_remainingdata[0]/content_db[0];
//   err_sf[0] = sqrt(pow(err_remainingdata[0]/content_db[0],2) + pow(content_data[0]/content_db[0]/content_db[0]*err_db[0],2));
//
//   //2.2. cross-check
//   cout << "----- Content Diboson before scaling -----" << endl;
//   cout << "0 BTag: " << content_db[0] << endl;
//   cout << "1 BTag: " << content_db[1] << endl;
//   cout << "2 BTag: " << content_db[2] << endl << endl;
//
//   cout << "----- Content non-Diboson before scaling -----" << endl;
//   cout << "0 BTag: " << content_nondb[0] << endl;
//   cout << "1 BTag: " << content_nondb[1] << endl;
//   cout << "2 BTag: " << content_nondb[2] << endl << endl;
//
//   cout << "----- Content total DATA before scaling -----" << endl;
//   cout << "0 BTag: " << content_data[0] << endl;
//   cout << "1 BTag: " << content_data[1] << endl;
//   cout << "2 BTag: " << content_data[2] << endl << endl;
//
//   cout << "----- Content remaining DATA before scaling -----" << endl;
//   cout << "0 BTag: " << content_remainingdata[0] << endl;
//   cout << "1 BTag: " << content_remainingdata[1] << endl;
//   cout << "2 BTag: " << content_remainingdata[2] << endl << endl;
//
//   //3. Scale diboson with this overall SF
//   h_diboson->Scale(sf[0]);
//
//   //4. Obtain sfs and errors for 1-tag and 2-tag bins
//   for(int i=0; i<3; i++){
//     content_db[i] = h_diboson->GetBinContent(i+1);
//     err_db[i] = h_diboson->GetBinError(i+1);
//     if(i>0){
//       sf[i] = content_remainingdata[i]/content_db[i];
//       err_sf[i] = sqrt(pow(err_remainingdata[i]/content_db[i],2) + pow(content_data[i]/content_db[i]/content_db[i]*err_db[i],2));
//     }
//   }
//
//   //5. Cross-check
//   cout << "sf 0tags: " << sf[0] << endl;
//   cout << "sf 1tag : " << sf[1] << endl;
//   cout << "sf 2tags: " << sf[2] << endl << endl;
//
//
//   cout << "----- Content MC after scaling -----" << endl;
//   cout << "0 BTag: " << content_db[0]+content_nondb[0] << endl;
//   cout << "1 BTag: " << content_db[1]*sf[1]+content_nondb[1] << endl;
//   cout << "2 BTag: " << content_db[2]*sf[2]+content_nondb[2] << endl << endl;
//
//   cout << "----- Content DATA after scaling -----" << endl;
//   cout << "0 BTag: " << content_data[0] << endl;
//   cout << "1 BTag: " << content_data[1] << endl;
//   cout << "2 BTag: " << content_data[2] << endl;
//
//   //6. Prepare output
//   //One TGraphAsymmErrors for xsec, one for btag 1&2
//   TH1D* Diboson_XSec_SF = new TH1D("Diboson_XSec_SF", "Diboson XSec SF;;SF", 1, -0.5,0.5);
//   Diboson_XSec_SF->SetBinContent(1,sf[0]);
//   Diboson_XSec_SF->SetBinError(1,err_sf[0]);
//   TH1D* Diboson_BTag_SF = new TH1D("Diboson_BTag_SF","Diboson BTag SF;N_{b-tags};SF",2,0.5,2.5);
//   Diboson_BTag_SF->SetBinContent(1,sf[1]);
//   Diboson_BTag_SF->SetBinContent(2,sf[2]);
//   Diboson_BTag_SF->SetBinError(1,err_sf[1]);
//   Diboson_BTag_SF->SetBinError(2,err_sf[2]);
//
//
//   //7. Write everything
//   unique_ptr<TFile> out;
//   out.reset(new TFile(FakeRateTool::base_path_Diboson + "DibosonSF.root","RECREATE"));
//   Diboson_XSec_SF->Write();
//   Diboson_BTag_SF->Write();
//   out->Close();
//
//
//   delete Diboson_BTag_SF;
//   delete Diboson_XSec_SF;
//   delete h_total_before;
//   delete h_data;
//   delete h_ttbar;
//   delete h_singletop;
//   delete h_dy;
//   delete h_diboson;
//
// }
//
//
//
//
//
//
//
//
//
// #include "../include/Tools.h"
// #include <TH1D.h>
// #include <TFile.h>
// #include <TGraphAsymmErrors.h>
// #include <TCanvas.h>
// #include <iostream>
// #include <TMath.h>
// #include <TLine.h>
// #include <TLatex.h>
// #include <TH1D.h>
// #include <TAxis.h>
// #include <TGaxis.h>
// #include "../include/cosmetics.h"
// #include <TLegend.h>
//
//
// using namespace std;
//
// void FakeRateTool::CalculateMuonFakeRateSF_DibosonNLO(){
//
//   unique_ptr<TFile> in_data, in_dy, in_db, in_ttbar, in_wjets, in_singletop, in_qcd, in_ttv;
//   in_data.reset(new TFile(FakeRateTool::base_path_DY_DibosonSF_Muon + "uhh2.AnalysisModuleRunner.DATA.DATA.root","READ"));
//   in_dy.reset(new TFile(FakeRateTool::base_path_DY_DibosonSF_Muon + "uhh2.AnalysisModuleRunner.MC.DYJets.root","READ"));
//   in_db.reset(new TFile(FakeRateTool::base_path_DY_DibosonSF_Muon + "uhh2.AnalysisModuleRunner.MC.DibosonNLO.root","READ"));
//   in_ttbar.reset(new TFile(FakeRateTool::base_path_DY_DibosonSF_Muon + "uhh2.AnalysisModuleRunner.MC.TTbar.root","READ"));
//   in_wjets.reset(new TFile(FakeRateTool::base_path_DY_DibosonSF_Muon + "uhh2.AnalysisModuleRunner.MC.WJets.root","READ"));
//   in_singletop.reset(new TFile(FakeRateTool::base_path_DY_DibosonSF_Muon + "uhh2.AnalysisModuleRunner.MC.SingleTop.root","READ"));
//   in_qcd.reset(new TFile(FakeRateTool::base_path_DY_DibosonSF_Muon + "uhh2.AnalysisModuleRunner.MC.QCD.root","READ"));
//   in_ttv.reset(new TFile(FakeRateTool::base_path_DY_DibosonSF_Muon + "uhh2.AnalysisModuleRunner.MC.TTV.root","READ"));
//
//   TH1D *h_all_before_data, *h_matched_after_data, *h_all_before_dy, *h_1mu_dy, *h_real_after_dy, *h_fake_after_dy, *h_all_before_db, *h_1mu_db, *h_real_after_db, *h_fake_after_db, *h_all_before_ttbar, *h_1mu_ttbar, *h_real_after_ttbar, *h_fake_after_ttbar, *h_all_before_wjets, *h_1mu_wjets, *h_real_after_wjets, *h_fake_after_wjets, *h_all_before_singletop, *h_1mu_singletop, *h_real_after_singletop, *h_fake_after_singletop, *h_all_before_qcd, *h_1mu_qcd, *h_real_after_qcd, *h_fake_after_qcd, *h_all_before_ttv, *h_1mu_ttv, *h_real_after_ttv, *h_fake_after_ttv;
//
//   h_all_before_data = (TH1D*)in_data->Get("FakeRate_DibosonSF/Int_events_denom");
//   h_matched_after_data = (TH1D*)in_data->Get("FakeRate_FinalSelection/Int_events_1mu");
//
//   h_all_before_dy = (TH1D*)in_dy->Get("FakeRate_DibosonSF/Int_events_denom");
//   h_1mu_dy = (TH1D*)in_dy->Get("FakeRate_DibosonSF/Int_events_1mu");
//   h_real_after_dy = (TH1D*)in_dy->Get("FakeRate_FinalSelection/Int_events_1mu_real");
//   h_fake_after_dy = (TH1D*)in_dy->Get("FakeRate_FinalSelection/Int_events_1mu_fake");
//   h_all_before_db = (TH1D*)in_db->Get("FakeRate_DibosonSF/Int_events_denom");
//   h_1mu_db = (TH1D*)in_db->Get("FakeRate_DibosonSF/Int_events_1mu");
//   h_real_after_db = (TH1D*)in_db->Get("FakeRate_FinalSelection/Int_events_1mu_real");
//   h_fake_after_db = (TH1D*)in_db->Get("FakeRate_FinalSelection/Int_events_1mu_fake");
//   h_all_before_ttbar = (TH1D*)in_ttbar->Get("FakeRate_DibosonSF/Int_events_denom");
//   h_1mu_ttbar = (TH1D*)in_ttbar->Get("FakeRate_DibosonSF/Int_events_1mu");
//   h_real_after_ttbar = (TH1D*)in_ttbar->Get("FakeRate_FinalSelection/Int_events_1mu_real");
//   h_fake_after_ttbar = (TH1D*)in_ttbar->Get("FakeRate_FinalSelection/Int_events_1mu_fake");
//   h_all_before_wjets = (TH1D*)in_wjets->Get("FakeRate_DibosonSF/Int_events_denom");
//   h_1mu_wjets = (TH1D*)in_wjets->Get("FakeRate_DibosonSF/Int_events_1mu");
//   h_real_after_wjets = (TH1D*)in_wjets->Get("FakeRate_FinalSelection/Int_events_1mu_real");
//   h_fake_after_wjets = (TH1D*)in_wjets->Get("FakeRate_FinalSelection/Int_events_1mu_fake");
//   h_all_before_singletop = (TH1D*)in_singletop->Get("FakeRate_DibosonSF/Int_events_denom");
//   h_1mu_singletop = (TH1D*)in_singletop->Get("FakeRate_DibosonSF/Int_events_1mu");
//   h_real_after_singletop = (TH1D*)in_singletop->Get("FakeRate_FinalSelection/Int_events_1mu_real");
//   h_fake_after_singletop = (TH1D*)in_singletop->Get("FakeRate_FinalSelection/Int_events_1mu_fake");
//   h_all_before_qcd = (TH1D*)in_qcd->Get("FakeRate_DibosonSF/Int_events_denom");
//   h_1mu_qcd = (TH1D*)in_qcd->Get("FakeRate_DibosonSF/Int_events_1mu");
//   h_real_after_qcd = (TH1D*)in_qcd->Get("FakeRate_FinalSelection/Int_events_1mu_real");
//   h_fake_after_qcd = (TH1D*)in_qcd->Get("FakeRate_FinalSelection/Int_events_1mu_fake");
//   h_all_before_ttv = (TH1D*)in_ttv->Get("FakeRate_DibosonSF/Int_events_denom");
//   h_1mu_ttv = (TH1D*)in_ttv->Get("FakeRate_DibosonSF/Int_events_1mu");
//   h_real_after_ttv = (TH1D*)in_ttv->Get("FakeRate_FinalSelection/Int_events_1mu_real");
//   h_fake_after_ttv = (TH1D*)in_ttv->Get("FakeRate_FinalSelection/Int_events_1mu_fake");
//
//
//   //1. add all histograms from MC
//   TH1D* h_all_before_mc, *h_real_after_mc, *h_fake_after_mc, *h_1mu_mc;
//   h_all_before_mc = (TH1D*)h_all_before_dy->Clone("h_all_before_mc");
//   h_1mu_mc = (TH1D*)h_1mu_dy->Clone("h_1mu_mc");
//   h_real_after_mc = (TH1D*)h_real_after_dy->Clone("h_real_after_mc");
//   h_fake_after_mc = (TH1D*)h_fake_after_dy->Clone("h_fake_after_mc");
//   h_all_before_mc->Add(h_all_before_db);
//   h_1mu_mc->Add(h_1mu_db);
//   h_real_after_mc->Add(h_real_after_db);
//   h_fake_after_mc->Add(h_fake_after_db);
//   h_all_before_mc->Add(h_all_before_ttbar);
//   h_1mu_mc->Add(h_1mu_ttbar);
//   h_real_after_mc->Add(h_real_after_ttbar);
//   h_fake_after_mc->Add(h_fake_after_ttbar);
//   h_all_before_mc->Add(h_all_before_wjets);
//   h_1mu_mc->Add(h_1mu_wjets);
//   h_real_after_mc->Add(h_real_after_wjets);
//   h_fake_after_mc->Add(h_fake_after_wjets);
//   h_all_before_mc->Add(h_all_before_singletop);
//   h_1mu_mc->Add(h_1mu_singletop);
//   h_real_after_mc->Add(h_real_after_singletop);
//   h_fake_after_mc->Add(h_fake_after_singletop);
//   h_all_before_mc->Add(h_all_before_qcd);
//   h_1mu_mc->Add(h_1mu_qcd);
//   h_real_after_mc->Add(h_real_after_qcd);
//   h_fake_after_mc->Add(h_fake_after_qcd);
//   h_all_before_mc->Add(h_all_before_ttv);
//   h_1mu_mc->Add(h_1mu_ttv);
//   h_real_after_mc->Add(h_real_after_ttv);
//   h_fake_after_mc->Add(h_fake_after_ttv);
//
//   //1.1 Calculate DATA/MC SF for 0mu region to get rid of other effects on data/mc ratio
//   TH1D* h_0mu_data, *h_0mu_mc;
//   h_0mu_mc = (TH1D*)h_all_before_mc->Clone("h_0mu_mc");
//   h_0mu_data = (TH1D*)h_all_before_data->Clone("h_0mu_data");
//   h_0mu_mc->Add(h_1mu_mc,-1);
//   h_0mu_data->Add(h_matched_after_data,-1);
//
//   double SF_0mu = h_0mu_data->GetBinContent(1) / h_0mu_mc->GetBinContent(1);
//   double err_0mu = sqrt(pow(h_0mu_data->GetBinError(1)/h_0mu_mc->GetBinContent(1),2) + pow(h_0mu_data->GetBinContent(1)/h_0mu_mc->GetBinContent(1)/h_0mu_mc->GetBinContent(1)*h_0mu_mc->GetBinError(1),2));
//
//   /*
//   double new_err_before = sqrt(pow(SF_0mu*h_all_before_mc->GetBinError(1),2) + pow(err_0mu*h_all_before_mc->GetBinContent(1),2));
//   double new_err_real_after = sqrt(pow(SF_0mu*h_real_after_mc->GetBinError(1),2) + pow(err_0mu*h_real_after_mc->GetBinContent(1),2));
//   double new_err_fake_after = sqrt(pow(SF_0mu*h_fake_after_mc->GetBinError(1),2) + pow(err_0mu*h_fake_after_mc->GetBinContent(1),2));
//   h_all_before_mc->SetBinContent(1,h_all_before_mc->GetBinContent(1) * SF_0mu);
//   h_real_after_mc->SetBinContent(1,h_real_after_mc->GetBinContent(1) * SF_0mu);
//   h_fake_after_mc->SetBinContent(1,h_fake_after_mc->GetBinContent(1) * SF_0mu);
//   h_all_before_mc->SetBinError(1,new_err_before);
//   h_real_after_mc->SetBinError(1,new_err_real_after);
//   h_fake_after_mc->SetBinError(1,new_err_fake_after);
//   cout << "scaling MC with " << SF_0mu << " globally" << endl;
//   */
//   //2.  Calculate (de)nom for data and mc
//
//   TH1D* h_nominator_data, *h_denominator_data, *h_nominator_mc, *h_denominator_mc;
//   h_nominator_data  = (TH1D*)h_matched_after_data->Clone("h_nominator_data");
//   cout << "total nom, data: " << h_nominator_data->GetBinContent(1) << endl;
//   h_nominator_data->Add(h_real_after_mc,-1);
//   cout << "real, mc: " << h_real_after_mc->GetBinContent(1) << endl;
//   h_denominator_data = (TH1D*)h_all_before_data->Clone("h_denominator_data");
//   h_nominator_mc = (TH1D*)h_fake_after_mc->Clone("h_nominator_mc");
//   h_denominator_mc = (TH1D*)h_all_before_mc->Clone("h_denominator_mc");
//
//   cout << "nom, data: " << h_nominator_data->GetBinContent(1) << " +- " << h_nominator_data->GetBinError(1) << endl;
//   cout << "denom, data: " << h_denominator_data->GetBinContent(1) << " +- " << h_denominator_data->GetBinError(1) << endl;
//   cout << "nom, mc: " << h_nominator_mc->GetBinContent(1) << " +- " << h_nominator_mc->GetBinError(1) << endl;
//   cout << "denom, mc: " << h_denominator_mc->GetBinContent(1) << " +- " << h_denominator_mc->GetBinError(1) << endl;
//
//
//   //3. Calculate fake rates for data and mc separately
//   TGraphAsymmErrors* gr_fakerate_mc, *gr_fakerate_data;
//   gr_fakerate_mc = new TGraphAsymmErrors(h_nominator_mc, h_denominator_mc);
//   gr_fakerate_data = new TGraphAsymmErrors(h_nominator_data, h_denominator_data);
//
//   unique_ptr<TCanvas> c1, c2;
//   c1.reset(new TCanvas());
//   gr_fakerate_mc->Draw();
//
//   c2.reset(new TCanvas());
//   gr_fakerate_data->Draw();
//
//
//   //4. Calculate Data/MC
//   double binsx[1] = {0};
//   double ex[1] = {0.5};
//   double binsy[1], eyl[1], eyh[1];
//   for(int i=0; i<1; i++){
//     double xmc, ymc, xdata, ydata, ydatalow, ydatahigh, ymclow, ymchigh;
//     gr_fakerate_mc->GetPoint(i,xmc,ymc);
//     gr_fakerate_data->GetPoint(i,xdata,ydata);
//     ymclow = gr_fakerate_mc->GetErrorYlow(i);
//     ymchigh = gr_fakerate_mc->GetErrorYhigh(i);
//     ydatalow = gr_fakerate_data->GetErrorYlow(i);
//     ydatahigh = gr_fakerate_data->GetErrorYhigh(i);
//
//     cout << "DATA: " << ydata << " + " << ydatahigh << " - " << ydatalow << endl;
//     cout << "MC: " << ymc << " + " << ymchigh << " - " << ymclow << endl;
//
//     binsy[i] = ydata/ymc;
//     eyh[i] = sqrt(pow(ydatahigh/ymc,2) + pow(ydata*ymclow/ymc/ymc,2));
//     eyl[i] = sqrt(pow(ydatalow/ymc,2) + pow(ydata*ymchigh/ymc/ymc,2));
//     cout << "RATIO: " << binsy[i] << " + " << eyh[i] << " - " << eyl[i] << endl;
//
//   }
//
//
//   TGraphAsymmErrors* SF;
//   SF = new TGraphAsymmErrors(1, binsx, binsy, ex, ex, eyl, eyh);
//   unique_ptr<TCanvas> c3;
//   c3.reset(new TCanvas());
//   SF->Draw();
//
//
//   TPad* p_top = SetupRatioPadTop();
//   TPad* p_bot = SetupRatioPad();
//   TCanvas* c_final = new TCanvas("c_final", "Ratio plot", 400, 400);
//   p_top->Draw();
//   p_bot->Draw();
//   p_top->cd();
//   HistCosmetics(gr_fakerate_mc,false);
//   HistCosmetics(gr_fakerate_data,false);
//   HistCosmetics(SF,true);
//
//   gr_fakerate_mc->SetLineColor(kRed);
//   gr_fakerate_data->SetLineColor(kBlack);
//   gr_fakerate_data->SetMarkerColor(kBlack);
//   gr_fakerate_mc->SetMarkerColor(kRed);
//   double max1 = TMath::MaxElement(gr_fakerate_data->GetN(), gr_fakerate_data->GetY());
//   double max2 = TMath::MaxElement(gr_fakerate_mc->GetN(), gr_fakerate_mc->GetY());
//   double maxall  = max(max1,max2);
//   gr_fakerate_data->GetHistogram()->SetMaximum(maxall*1.5);
//   gr_fakerate_data->GetHistogram()->SetMinimum(min(0.000002,min(TMath::MinElement(gr_fakerate_data->GetN(), gr_fakerate_data->GetY()),TMath::MinElement(gr_fakerate_mc->GetN(), gr_fakerate_mc->GetY())) - 0.5*min(TMath::MinElement(gr_fakerate_data->GetN(), gr_fakerate_data->GetY()),TMath::MinElement(gr_fakerate_mc->GetN(), gr_fakerate_mc->GetY()))));
//   gr_fakerate_data->GetXaxis()->SetLimits(-0.5,0.5);
//   gr_fakerate_data->GetYaxis()->SetTitle("Misidentification rate");
//   gr_fakerate_data->SetMarkerStyle(20);
//   gr_fakerate_data->SetMarkerSize(0.9);
//
//   gr_fakerate_data->Draw("AP");
//   gr_fakerate_mc->Draw("P SAME");
//
//   TString cmstext = "CMS";
//   TLatex *text2 = new TLatex(3.5, 24, cmstext);
//   text2->SetNDC();
//   text2->SetTextAlign(13);
//   text2->SetX(0.24);
//   text2->SetTextFont(62);
//   text2->SetTextSize(0.08);
//   text2->SetY(0.89);
//   text2->Draw();
//
//   TString preltext = "Preliminary";
//   TLatex *text3 = new TLatex(3.5, 24, preltext);
//   text3->SetNDC();
//   text3->SetTextAlign(13);
//   text3->SetX(0.24);
//   text3->SetTextFont(52);
//   text3->SetTextSize(0.06);
//   text3->SetY(0.80);
//   text3->Draw();
//
//   TLegend* leg_top = new TLegend(0.7,0.7,0.9,0.9);
//   leg_top->AddEntry(gr_fakerate_data,"Data","lep");
//   leg_top->AddEntry(gr_fakerate_mc,"MC","le");
//   leg_top->SetBorderSize(0);
//   leg_top->SetFillStyle(0);
//   leg_top->Draw();
//
//   p_bot->cd();
//
//   SF->SetLineColor(kBlack);
//   SF->GetYaxis()->SetRangeUser(0,5.4);
//   SF->GetXaxis()->SetLimits(-0.5,0.5);
//   SF->GetXaxis()->SetTitle("");
//   SF->GetYaxis()->SetTitle("Data / sim");
//   SF->SetMarkerStyle(20);
//   SF->SetMarkerSize(0.9);
//   SF->Draw("AP");
//   TLine* l_unity = new TLine(-0.5, 1, 0.5, 1);
//   l_unity->SetLineColor(kRed);
//   l_unity->SetLineWidth(2);
//   l_unity->SetLineStyle(2);
//   l_unity->Draw("SAME");
//
//
//   c1->SaveAs(FakeRateTool::base_path_DY_DibosonSF_Muon + "FakeRateMC_DibosonNLO.eps");
//   c2->SaveAs(FakeRateTool::base_path_DY_DibosonSF_Muon + "FakeRateDATA_DibosonNLO.eps");
//   c3->SaveAs(FakeRateTool::base_path_DY_DibosonSF_Muon + "ScaleFactors_DibosonNLO.eps");
//   c_final->SaveAs(FakeRateTool::base_path_DY_DibosonSF_Muon + "RatioPlot_DibosonNLO.eps");
//
//
//   SF->SetName("ScaleFactors");
//   unique_ptr<TFile> f_out;
//   f_out.reset(new TFile(FakeRateTool::base_path_DY_DibosonSF_Muon + "MuonFakeRateSF_DibosonNLO.root","RECREATE"));
//   SF->Write();
//   f_out->Close();
//
//
//   delete p_top;
//   delete p_bot;
//   delete c_final;
//   delete SF;
//   delete gr_fakerate_mc;
//   delete gr_fakerate_data;
//   delete h_all_before_data;
//   delete h_matched_after_data;
//   delete h_all_before_dy;
//   delete h_real_after_dy;
//   delete h_fake_after_dy;
//   delete h_all_before_db;
//   delete h_real_after_db;
//   delete h_fake_after_db;
//   delete h_all_before_ttbar;
//   delete h_real_after_ttbar;
//   delete h_fake_after_ttbar;
//   delete h_all_before_wjets;
//   delete h_real_after_wjets;
//   delete h_fake_after_wjets;
//   delete h_all_before_singletop;
//   delete h_real_after_singletop;
//   delete h_fake_after_singletop;
//   delete h_all_before_qcd;
//   delete h_real_after_qcd;
//   delete h_fake_after_qcd;
//   delete h_all_before_mc;
//   delete h_real_after_mc;
//   delete h_fake_after_mc;
//   delete h_nominator_data;
//   delete h_denominator_data;
//   delete h_nominator_mc;
//   delete h_denominator_mc;
//
//
//
//
//
//
//
// }
