


#include "TString.h"
#include <iostream>
#include "../include/EleTriggerEffRunner.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include <vector>
#include "../include/cosmetics.h"
#include <TROOT.h>
#include <TH1.h>
#include <TRint.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TObjString.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TLine.h>

using namespace std;

vector<TString> histname;
vector<TString> foldertagbefore;
vector<TString> foldertagafter;
vector<TString> plotname;

void do_stuff(TString inpath, TString outpath, TString histname);
void EleTriggerEffRunner::EleTriggerEff(){


  histname = {"pt_ele_fine", "pt_ele_binned", "eta_ele", "eta_ele_binned", "dR_ele_mu", "dRmin_ele_jet", "dRmin_ele_obj"};
  foldertagbefore = {"EletriggerEff_mutrigger/", "EletriggerEff_trigger30_beforeletrigger/", "EletriggerEff_trigger120_beforeletrigger/", "EletriggerEff_trigger30_neg_120_beforeletrigger/"};
  foldertagafter = {"EletriggerEff_eletriger123/", "EletriggerEff_trigger30_aftereletrigger/", "EletriggerEff_trigger120_aftereletrigger/", "EletriggerEff_trigger30_neg_120_aftereletrigger/"};
  plotname = {"eff_inclusivept_", "eff_30pt_", "eff_120pt_", "eff_30ptneg120pt_"};

  // foldertagbefore = {"EletriggerEff_mutrigger/", "EletriggerEff_trigger30_beforeletrigger/", "EletriggerEff_trigger120_beforeletrigger/", "EletriggerEff_trigger30_neg_120_beforeletrigger/"}
  // foldertagafter = {"EletriggerEff_eletriger123/", "EletriggerEff_trigger30_aftereletrigger/", "EletriggerEff_trigger120_aftereletrigger/", "EletriggerEff_trigger30_neg_120_aftereletrigger/"}
  // plotname = {"eff_inclusivept_", "eff_30pt_", "eff_120pt_", "eff_30ptneg120pt_"}

  for (unsigned int i = 0; i < histname.size(); i++)
  {
    do_stuff(EleTriggerEffRunner::inpath, EleTriggerEffRunner::outpath, histname[i]);
  }
}

void do_stuff(TString inpath, TString outpath, TString histname){



  for (int i = 0; i < 4; i++)
  {
    TString infilename_ttbar = inpath + "uhh2.AnalysisModuleRunner.MC.TTbar.root";
    TFile* in_ttbar = new TFile(infilename_ttbar,"READ");

    TString infilename_data = inpath + "uhh2.AnalysisModuleRunner.DATA.DATA.root";
    TFile* in_data = new TFile(infilename_data,"READ");


    TH1F* h_ttbar_before = (TH1F*)in_ttbar->Get(foldertagbefore[i] + histname);
    TH1F* h_ttbar_after = (TH1F*)in_ttbar->Get(foldertagafter[i] + histname);


    TH1F* h_data_before = (TH1F*)in_data->Get(foldertagbefore[i] + histname);
    TH1F* h_data_after = (TH1F*)in_data->Get(foldertagafter[i] + histname);

    TGraphAsymmErrors* eff_data = new TGraphAsymmErrors(h_data_after, h_data_before);
    TString xaxistitle = eff_data->GetTitle();
    HistCosmetics(eff_data, false);
    eff_data->SetLineColor(kBlack);
    eff_data->SetMarkerColor(kBlack);
    // eff_data->SetLineWidth(2);
    eff_data->GetYaxis()->SetTitle("Trigger efficiency");
    eff_data->SetTitle("");

    TGraphAsymmErrors* eff_ttbar = new TGraphAsymmErrors(h_ttbar_after, h_ttbar_before);
    HistCosmetics(eff_ttbar);
    eff_ttbar->SetLineColor(kRed+1);
    eff_ttbar->SetMarkerColor(kRed+1);
    // eff_ttbar->SetLineWidth(2);

    TCanvas* c1 = new TCanvas("c1", "Nice histogram", 600, 600);
    TPad* pad_top = SetupRatioPadTop();
    TPad* pad_bot = SetupRatioPad();
    pad_top->Draw();
    pad_bot->Draw();
    pad_top->cd();
    eff_data->Draw("APZ");
    // eff_data->SetMarkerStyle(21);
    eff_ttbar->Draw("PZ SAME");
    // eff_ttbar->SetMarkerStyle(21);

    //LEGEND
    //---------------
    //
    // Float_t top = 0.92;
    // Float_t ysize = 0.07;
    // Float_t xleft = 0.65;
    // Float_t xright = 0.92;
    //



    TLegend *leg = new TLegend(.4,.20,.6,.40);
    // TLegend* leg = new TLegend(xleft,top-ysize,xright,top, NULL,"brNDC");
    //leg->SetMarkerStyle(20);
    leg->SetBorderSize(0);
    leg->AddEntry(eff_data,"Data", "lep");
    // leg->AddEntry(eff_data,"Data", "e");
    leg->AddEntry(eff_ttbar,"MC", "lep");
    // leg->AddEntry(eff_ttbar,"MC", "e");

    leg->Draw();
    // oops we forgot the blue line... add it after
    //leg->AddEntry(fun2,
    //"#sqrt{2#pi} P_{T} (#gamma) latex  formula","f");
    // and add a header (or "title") for the legend
    // leg->SetHeader("The Legend Title");
    leg->Draw();




    //----------------

    pad_bot->cd();
    //calculate proper ratio
    const int n_bins = eff_data->GetN();
    int at_in_mc = 0, at_in_data = 0, n_points = 0;
    vector<double> SF_x, SF_y, SF_x_high, SF_x_low, SF_y_high, SF_y_low;
    for(int i=0; i<n_bins; i++){
      double x_MC, y_MC, x_DATA, y_DATA;
      double MC_x_high, MC_x_low, MC_y_high, MC_y_low;
      double DATA_y_high, DATA_y_low;

      eff_ttbar->GetPoint(at_in_mc,x_MC,y_MC);
      eff_data->GetPoint(at_in_data,x_DATA,y_DATA);

      double last_data = -1, last_mc = -1;
      bool skip = false;
      while(x_MC != x_DATA){
        if(x_MC < x_DATA){
          at_in_mc++;
          eff_ttbar->GetPoint(at_in_mc,x_MC,y_MC);
        }
        else if(x_DATA < x_MC){
          at_in_data++;
          eff_data->GetPoint(at_in_data,x_DATA,y_DATA);
        }

        if(x_MC == last_mc && x_DATA == last_data){
          skip = true;
          break;
        }
        last_mc = x_MC;
        last_data = x_DATA;
      }
      if(skip) break;

      MC_x_high = eff_ttbar->GetErrorXhigh(i);
      MC_x_low = eff_ttbar->GetErrorXlow(i);
      MC_y_high = eff_ttbar->GetErrorYhigh(i);
      MC_y_low = eff_ttbar->GetErrorYlow(i);

      DATA_y_high = eff_data->GetErrorYhigh(i);
      DATA_y_low = eff_data->GetErrorYlow(i);

      //gaussian error propagation
      SF_x.push_back(x_MC);
      SF_x_low.push_back(MC_x_low);
      SF_x_high.push_back(MC_x_high);
      SF_y.push_back(y_DATA/y_MC);
      SF_y_low.push_back(sqrt(pow(DATA_y_low/y_MC,2) + pow(y_DATA/y_MC/y_MC*MC_y_high,2)));
      SF_y_high.push_back(sqrt(pow(DATA_y_high/y_MC,2) + pow(y_DATA/y_MC/y_MC*MC_y_low,2)));

      // cout << "x-data: " << x_DATA << ", x-mc: " << x_MC << ", x-sf: " << SF_x.back() << endl;
      at_in_mc++;
      at_in_data++;
      n_points++;
    }


    TGraphAsymmErrors* scale_factor = new TGraphAsymmErrors(n_bins, &SF_x[0], &SF_y[0], &SF_x_low[0], &SF_x_high[0], &SF_y_low[0], &SF_y_high[0]);
    HistCosmetics(scale_factor, true);
    scale_factor->GetXaxis()->SetTitle(xaxistitle);

    // if condition

    if (histname == "eta_ele_binned" && (plotname[i] == "eff_120pt_" || plotname[i] == "eff_30ptneg120pt_" )) {
      /* code */

      TString outname = outpath + plotname[i] + "EleTriggerScaleFactors.root";
      TFile* outf = new TFile(outname,"UPDATE");
      scale_factor->Write();




      outf->Close();
    }
    scale_factor->Draw("APZ");


    TString outname_c1 = outpath + "Plots/" + plotname[i] + histname + ".eps";

    float xmin = scale_factor->GetHistogram()->GetXaxis()->GetXmin();
    float xmax = scale_factor->GetHistogram()->GetXaxis()->GetXmax();
    // float xmin = scale_factor->GetHistogram()->GetXaxis()->GetBinLowEdge(scale_factor->GetHistogram()->GetNbinsX()+1);
    scale_factor->GetYaxis()->SetTitle("Data / MC");
    TLine *l=new TLine(xmin,1.0,xmax,1.0);

    //TLine *l=new TLine(-1,1.0,6,1.0);
    l->SetLineColor(kRed);
    l->Draw();

    c1->SaveAs(outname_c1);
    delete c1;
  }
}
