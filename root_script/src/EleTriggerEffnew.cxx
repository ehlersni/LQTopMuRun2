// #include "TString.h"
// #include <iostream>
// #include "../include/EleTriggerEffRunner.h"
// #include "TFile.h"
// #include "TH1.h"
// #include "TH1F.h"
// #include "TCanvas.h"
// #include "TGraphAsymmErrors.h"
//
// using namespace std;
//
// void EleTriggerEffRunner::EleTriggerEff(){
//
//
//     TString infilename_ttbar = EleTriggerEffRunner::inpath + "uhh2.AnalysisModuleRunner.MC.TTbar.root";
//     TFile* in_ttbar = new TFile(infilename_ttbar,"READ");
//
//     TString infilename_data = EleTriggerEffRunner::inpath + "uhh2.AnalysisModuleRunner.DATA.DATA.root";
//     TFile* in_data = new TFile(infilename_data,"READ");
//
//
//     TH1F* h_inclusivept_ptfine_ttbar_before = (TH1F*)in_ttbar->Get("EletriggerEff_mutrigger/pt_ele_fine");
//     TH1F* h_inclusivept_ptfine_ttbar_after = (TH1F*)in_ttbar->Get("EletriggerEff_eletriger123/pt_ele_fine");
//
//
//     TH1F* h_inclusivept_ptfine_data_before = (TH1F*)in_data->Get("EletriggerEff_mutrigger/pt_ele_fine");
//     TH1F* h_inclusivept_ptfine_data_after = (TH1F*)in_data->Get("EletriggerEff_eletriger123/pt_ele_fine");
//
//
//
//     TH1F* h_30pt_ptfine_ttbar_before = (TH1F*)in_ttbar->Get("EletriggerEff_trigger30_beforeletrigger/pt_ele_fine");
//     TH1F* h_30pt_ptfine_ttbar_after = (TH1F*)in_ttbar->Get("EletriggerEff_trigger30_aftereletrigger/pt_ele_fine");
//
//     TH1F* h_30pt_ptfine_data_before = (TH1F*)in_data->Get("EletriggerEff_trigger30_beforeletrigger/pt_ele_fine");
//     TH1F* h_30pt_ptfine_data_after = (TH1F*)in_data->Get("EletriggerEff_trigger30_aftereletrigger/pt_ele_fine");
//
//
//     TH1F* h_120pt_ptfine_ttbar_before = (TH1F*)in_ttbar->Get("EletriggerEff_trigger120_beforeletrigger/pt_ele_fine");
//     TH1F* h_120pt_ptfine_ttbar_after = (TH1F*)in_ttbar->Get("EletriggerEff_trigger120_aftereletrigger/pt_ele_fine");
//
//     TH1F* h_120pt_ptfine_data_before = (TH1F*)in_data->Get("EletriggerEff_trigger120_beforeletrigger/pt_ele_fine");
//     TH1F* h_120pt_ptfine_data_after = (TH1F*)in_data->Get("EletriggerEff_trigger120_aftereletrigger/pt_ele_fine");
//
//
//
//     TH1F* h_30ptneg120pt_ptfine_ttbar_before = (TH1F*)in_ttbar->Get("EletriggerEff_trigger30_neg_120_beforeletrigger/pt_ele_fine");
//     TH1F* h_30ptneg120pt_ptfine_ttbar_after = (TH1F*)in_ttbar->Get("EletriggerEff_trigger30_neg_120_aftereletrigger/pt_ele_fine");
//
//     TH1F* h_30ptneg120pt_ptfine_data_before = (TH1F*)in_data->Get("EletriggerEff_trigger30_neg_120_beforeletrigger/pt_ele_fine");
//     TH1F* h_30ptneg120pt_ptfine_data_after = (TH1F*)in_data->Get("EletriggerEff_trigger30_neg_120_aftereletrigger/pt_ele_fine");
//
//
//
//
//
//
//     TGraphAsymmErrors* eff_inclusivept_ptfine_data = new TGraphAsymmErrors(h_inclusivept_ptfine_data_after, h_inclusivept_ptfine_data_before);
//     eff_inclusivept_ptfine_data->SetLineColor(kBlack);
//     eff_inclusivept_ptfine_data->SetMarkerColor(kBlack);
//     eff_inclusivept_ptfine_data->SetLineWidth(2);
//     eff_inclusivept_ptfine_data->GetXaxis()->SetTitle(eff_inclusivept_ptfine_data->GetTitle());
//     eff_inclusivept_ptfine_data->GetYaxis()->SetTitle("Trigger efficiency");
//     eff_inclusivept_ptfine_data->SetTitle("");
//
//     TGraphAsymmErrors* eff_inclusivept_ptfine_ttbar = new TGraphAsymmErrors(h_inclusivept_ptfine_ttbar_after, h_inclusivept_ptfine_ttbar_before);
//     eff_inclusivept_ptfine_ttbar->SetLineColor(kRed+1);
//     eff_inclusivept_ptfine_ttbar->SetMarkerColor(kRed+1);
//     eff_inclusivept_ptfine_ttbar->SetLineWidth(2);
//
//
//
//     TGraphAsymmErrors* eff_30pt_ptfine_data = new TGraphAsymmErrors(h_30pt_ptfine_data_after, h_30pt_ptfine_data_before);
//     eff_30pt_ptfine_data->SetLineColor(kBlack);
//     eff_30pt_ptfine_data->SetMarkerColor(kBlack);
//     eff_30pt_ptfine_data->SetLineWidth(2);
//     eff_30pt_ptfine_data->GetXaxis()->SetTitle(eff_30pt_ptfine_data->GetTitle());
//     eff_30pt_ptfine_data->GetYaxis()->SetTitle("Trigger efficiency");
//     eff_30pt_ptfine_data->SetTitle("");
//
//     TGraphAsymmErrors* eff_30pt_ptfine_ttbar = new TGraphAsymmErrors(h_30pt_ptfine_ttbar_after, h_30pt_ptfine_ttbar_before);
//     eff_30pt_ptfine_ttbar->SetLineColor(kRed+1);
//     eff_30pt_ptfine_ttbar->SetMarkerColor(kRed+1);
//     eff_30pt_ptfine_ttbar->SetLineWidth(2);
//
//
//
//
//
//     TGraphAsymmErrors* eff_120pt_ptfine_data = new TGraphAsymmErrors(h_120pt_ptfine_data_after, h_120pt_ptfine_data_before);
//     eff_120pt_ptfine_data->SetLineColor(kBlack);
//     eff_120pt_ptfine_data->SetMarkerColor(kBlack);
//     eff_120pt_ptfine_data->SetLineWidth(2);
//     eff_120pt_ptfine_data->GetXaxis()->SetTitle(eff_120pt_ptfine_data->GetTitle());
//     eff_120pt_ptfine_data->GetYaxis()->SetTitle("Trigger efficiency");
//     eff_120pt_ptfine_data->SetTitle("");
//
//     TGraphAsymmErrors* eff_120pt_ptfine_ttbar = new TGraphAsymmErrors(h_120pt_ptfine_ttbar_after, h_120pt_ptfine_ttbar_before);
//     eff_120pt_ptfine_ttbar->SetLineColor(kRed+1);
//     eff_120pt_ptfine_ttbar->SetMarkerColor(kRed+1);
//     eff_120pt_ptfine_ttbar->SetLineWidth(2);
//
//
//
//
//
//     TGraphAsymmErrors* eff_30ptneg120pt_ptfine_data = new TGraphAsymmErrors(h_30ptneg120pt_ptfine_data_after, h_30ptneg120pt_ptfine_data_before);
//     eff_30ptneg120pt_ptfine_data->SetLineColor(kBlack);
//     eff_30ptneg120pt_ptfine_data->SetMarkerColor(kBlack);
//     eff_30ptneg120pt_ptfine_data->SetLineWidth(2);
//     eff_30ptneg120pt_ptfine_data->GetXaxis()->SetTitle(eff_30ptneg120pt_ptfine_data->GetTitle());
//     eff_30ptneg120pt_ptfine_data->GetYaxis()->SetTitle("Trigger efficiency");
//     eff_30ptneg120pt_ptfine_data->SetTitle("");
//
//     TGraphAsymmErrors* eff_30ptneg120pt_ptfine_ttbar = new TGraphAsymmErrors(h_30ptneg120pt_ptfine_ttbar_after, h_30ptneg120pt_ptfine_ttbar_before);
//     eff_30ptneg120pt_ptfine_ttbar->SetLineColor(kRed+1);
//     eff_30ptneg120pt_ptfine_ttbar->SetMarkerColor(kRed+1);
//     eff_30ptneg120pt_ptfine_ttbar->SetLineWidth(2);
//
//
//
//
//
//     TCanvas* c1 = new TCanvas("c1", "Nice histogram1", 600, 600);
//     eff_inclusivept_ptfine_data->Draw("AP");
//     eff_inclusivept_ptfine_ttbar->Draw("P SAME");
//
//     TString outname_c1 = EleTriggerEffRunner::outpath + "Plots/eff_inclusivept_ptfine.eps";
//     c1->SaveAs(outname_c1);
//
//
//     TCanvas* c2 = new TCanvas("c2", "Nice histogram2", 600, 600);
//     eff_30pt_ptfine_data->Draw("AP");
//     eff_30pt_ptfine_ttbar->Draw("P SAME");
//
//     TString outname_c2 = EleTriggerEffRunner::outpath + "Plots/eff_30pt_ptfine.eps";
//     c2->SaveAs(outname_c2);
//
//
//
//     TCanvas* c3 = new TCanvas("c3", "Nice histogram3", 600, 600);
//     eff_120pt_ptfine_data->Draw("AP");
//     eff_120pt_ptfine_ttbar->Draw("P SAME");
//
//     TString outname_c3 = EleTriggerEffRunner::outpath + "Plots/eff_120pt_ptfine.eps";
//     c3->SaveAs(outname_c3);
//
//
//
//     TCanvas* c4 = new TCanvas("c4", "Nice histogram4", 600, 600);
//     eff_30ptneg120pt_ptfine_data->Draw("AP");
//     eff_30ptneg120pt_ptfine_ttbar->Draw("P SAME");
//
//     TString outname_c4 = EleTriggerEffRunner::outpath + "Plots/eff_30ptneg120pt_ptfine.eps";
//     c4->SaveAs(outname_c4);
// 
//
//
//
//     // TCanvas* c2 = new TCanvas("c2", "Nice histogram", 600, 600);
//     // eff_inclusivept_ptfine_data->Draw("AP");
//
//
//     // TString outname_c2 = EleTriggerEffRunner::outpath + "Plots/eff_inclusivept_ptfine_data.eps";
//     // c2->SaveAs(outname_c2);
//
//
//
//
//
//
// }
