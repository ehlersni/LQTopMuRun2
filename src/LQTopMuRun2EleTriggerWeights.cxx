#include "UHH2/LQTopMuRun2/include/LQTopMuRun2EleTriggerWeights.h"
#include "UHH2/core/include/Event.h"
#include "TString.h"
#include <iostream>
#include "TFile.h"
#include "TH1.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include <vector>
#include <TROOT.h>
#include <TH1.h>
#include <TRint.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TObjString.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TLine.h>

#include <stdexcept>

using namespace uhh2examples;
using namespace uhh2;
using namespace std;

bool lowpt;





ElectronTriggerWeights::ElectronTriggerWeights(Context & ctx, TString SysDirection_): SysDirection(SysDirection_){

  auto dataset_type = ctx.get("dataset_type");
  bool is_mc = dataset_type == "MC";
  if(!(is_mc == true)){
    cout << "Warning: ElectronTriggerWeights will not have an effect on this non-MC sample (dataset_type = '" + dataset_type + "')" << endl;
    return;
  }

  //unique_ptr<TFile> datalowpt, datahighpt;
  //file.reset(new TFile(path,"READ"));

  // Eff_lowpt.reset((TGraphAsymmErrors*file->Get("gr_lowpt_eta_TTbar_eff"));
  // Eff_highpt.reset((TGraphAsymmErrors*)file->Get("gr_highpt_eta_TTbar_eff"));


  TString infilename_datalowpt = "/nfs/dust/cms/user/ehlersni/LQTopMuRun2/fullselection/eletriggereff/root_script_output/eff_30ptneg120pt_EleTriggerScaleFactors.root";
  TFile* datalowpt = new TFile(infilename_datalowpt,"READ");

  TString infilename_datahighpt = "/nfs/dust/cms/user/ehlersni/LQTopMuRun2/fullselection/eletriggereff/root_script_output/eff_120pt_EleTriggerScaleFactors.root";
  TFile* datahighpt = new TFile(infilename_datahighpt,"READ");

  graph_datalowpt.reset((TGraphAsymmErrors*)datalowpt->Get("Graph"));
  graph_datahighpt.reset((TGraphAsymmErrors*)datahighpt->Get("Graph"));

  //Eff_lowpt_DATA.reset((TGraphAsymmErrors*)file->Get("gr_lowpt_eta_DATA_eff"));
  //Eff_highpt_DATA.reset((TGraphAsymmErrors*)file->Get("gr_highpt_eta_DATA_eff"));

  //if(SysDirection != "nominal" && SysDirection != "up" && SysDirection != "down") throw runtime_error("In LQToTopMuModules.cxx, ElectronTriggerWeights.process(): Invalid SysDirection specified.");

}

bool ElectronTriggerWeights::process(Event & event){

  if(event.isRealData) return true;

  const auto ele = event.electrons->at(0);
  double eta = ele.eta();
  if(fabs(eta) > 2.4) throw runtime_error("In LQToTopMuModules.cxx, ElectronTriggerWeights.process(): Ele-|eta| > 2.4 is not supported at the moment.");


  //find right bin in eta
  int idx = 0;
  //bool lowpt = false;
  if(30 <= ele.pt() && ele.pt() < 120){
    lowpt = true;
    //lowpt trigger
    bool keep_going = true;
    while(keep_going){
      double x,y;
      graph_datalowpt->GetPoint(idx,x,y);
      keep_going = eta > x + graph_datalowpt->GetErrorXhigh(idx);
      if(keep_going) idx++;
    }
  }
  else if(ele.pt() >= 120){
    //highpt trigger
    bool keep_going = true;
    while(keep_going){
      double x,y;
      graph_datahighpt->GetPoint(idx,x,y);
      keep_going = eta > x + graph_datahighpt->GetErrorXhigh(idx);
      if(keep_going) idx++;
    }
  }
  else throw runtime_error("In LQToTopMuModules.cxx, ElectronTriggerWeights.process(): Electron has pt<30. Clean electron collection before applying weights.");

  // access efficiencies for MC and DATA, possibly accout for systematics = statistical + add. 2% up/down
  double sf = -1, dummy_x;
  double stat = -1, tp = 0.02, total_syst = -1; // tp = tag and probe unsicherheit

  if(lowpt){
    graph_datalowpt->GetPoint(idx,dummy_x,sf);

    if(SysDirection == "up"){
      stat = graph_datalowpt->GetErrorYhigh(idx);
      total_syst = sqrt(pow(stat,2) + pow(tp,2));

      sf += total_syst;
    }
    else if(SysDirection == "down"){
      stat = graph_datalowpt->GetErrorYlow(idx);
      total_syst = sqrt(pow(stat,2) + pow(tp,2));

      sf -= graph_datalowpt->GetErrorYlow(idx);
    }
  }
  else{
    graph_datahighpt->GetPoint(idx,dummy_x,sf);

    if(SysDirection == "up"){
      stat = graph_datahighpt->GetErrorYhigh(idx);
      total_syst = sqrt(pow(stat,2) + pow(tp,2));

      sf += graph_datahighpt->GetErrorYhigh(idx);
    }
    else if(SysDirection == "down"){
      stat = graph_datahighpt->GetErrorYlow(idx);
      total_syst = sqrt(pow(stat,2) + pow(tp,2));

      sf -= graph_datahighpt->GetErrorYlow(idx);
    }
  }

  // cout << "Scale Factor:  " << sf << endl;
  event.weight *= sf;

  return true;
}



// ElectronTriggerWeights::ElectronTriggerWeights(Context & ctx, TString path_, TString SysDirection_): path(path_), SysDirection(SysDirection_){
//
//   auto dataset_type = ctx.get("dataset_type");
//   bool is_mc = dataset_type == "MC";
//   if(!is_mc){
//     cout << "Warning: ElectronTriggerWeights will not have an effect on this non-MC sample (dataset_type = '" + dataset_type + "')" << endl;
//     return;
//   }
//
//   unique_ptr<TFile> file;
//   file.reset(new TFile(path,"READ"));
//
//   Eff_lowpt_MC.reset((TGraphAsymmErrors*)file->Get("gr_lowpt_eta_TTbar_eff"));
//   Eff_highpt_MC.reset((TGraphAsymmErrors*)file->Get("gr_highpt_eta_TTbar_eff"));
//   Eff_lowpt_DATA.reset((TGraphAsymmErrors*)file->Get("gr_lowpt_eta_DATA_eff"));
//   Eff_highpt_DATA.reset((TGraphAsymmErrors*)file->Get("gr_highpt_eta_DATA_eff"));
//
//   if(SysDirection != "nominal" && SysDirection != "up" && SysDirection != "down") throw runtime_error("In LQToTopMuModules.cxx, ElectronTriggerWeights.process(): Invalid SysDirection specified.");
//
// }
//
// bool ElectronTriggerWeights::process(Event & event){
//
//   if(event.isRealData) return true;
//
//   const auto ele = event.electrons->at(0);
//   double eta = ele.eta();
//   if(fabs(eta) > 2.4) throw runtime_error("In LQToTopMuModules.cxx, ElectronTriggerWeights.process(): Ele-|eta| > 2.4 is not supported at the moment.");
//
//
//   //find right bin in eta
//   int idx = 0;
//   bool lowpt = false;
//   if(30 <= ele.pt() && ele.pt() < 120){
//     lowpt = true;
//     //lowpt trigger
//     bool keep_going = true;
//     while(keep_going){
//       double x,y;
//       Eff_lowpt_MC->GetPoint(idx,x,y);
//       keep_going = eta > x + Eff_lowpt_MC->GetErrorXhigh(idx);
//       if(keep_going) idx++;
//     }
//   }
//   else if(ele.pt() >= 120){
//     //highpt trigger
//     bool keep_going = true;
//     while(keep_going){
//       double x,y;
//       Eff_highpt_MC->GetPoint(idx,x,y);
//       keep_going = eta > x + Eff_highpt_MC->GetErrorXhigh(idx);
//       if(keep_going) idx++;
//     }
//   }
//   else throw runtime_error("In LQToTopMuModules.cxx, ElectronTriggerWeights.process(): Electron has pt<30. Clean electron collection before applying weights.");
//
//   //access efficiencies for MC and DATA, possibly accout for systematics = statistical + add. 2% up/down
//   double eff_data = -1, eff_mc = -1, dummy_x;
//   double stat_data = -1, stat_mc = -1, tp = 0.02, total_syst_data = -1, total_syst_mc = -1;
//   if(lowpt){
//     Eff_lowpt_MC->GetPoint(idx,dummy_x,eff_mc);
//     Eff_lowpt_DATA->GetPoint(idx,dummy_x,eff_data);
//
//     // if(SysDirection == "up"){
//     //   stat_mc = Eff_lowpt_MC->GetErrorYlow(idx);
//     //   stat_data = Eff_lowpt_DATA->GetErrorYhigh(idx);
//     //   total_syst_mc = sqrt(pow(stat_mc,2) + pow(tp,2));
//     //   total_syst_data = sqrt(pow(stat_data,2) + pow(tp,2));
//     //
//     //   eff_mc -= total_syst_mc;
//     //   eff_data += total_syst_data;
//     // }
//     // else if(SysDirection == "down"){
//     //   stat_mc = Eff_lowpt_MC->GetErrorYhigh(idx);
//     //   stat_data = Eff_lowpt_DATA->GetErrorYlow(idx);
//     //   total_syst_mc = sqrt(pow(stat_mc,2) + pow(tp,2));
//     //   total_syst_data = sqrt(pow(stat_data,2) + pow(tp,2));
//     //
//     //   eff_mc += Eff_lowpt_MC->GetErrorYhigh(idx);
//     //   eff_data -= Eff_lowpt_DATA->GetErrorYlow(idx);
//     // }
//   }
//   else{
//     Eff_highpt_MC->GetPoint(idx,dummy_x,eff_mc);
//     Eff_highpt_DATA->GetPoint(idx,dummy_x,eff_data);
//
//     // if(SysDirection == "up"){
//     //   stat_mc = Eff_highpt_MC->GetErrorYlow(idx);
//     //   stat_data = Eff_highpt_DATA->GetErrorYhigh(idx);
//     //   total_syst_mc = sqrt(pow(stat_mc,2) + pow(tp,2));
//     //   total_syst_data = sqrt(pow(stat_data,2) + pow(tp,2));
//     //
//     //   eff_mc -= Eff_highpt_MC->GetErrorYlow(idx);
//     //   eff_data += Eff_highpt_DATA->GetErrorYhigh(idx);
//     // }
//     // else if(SysDirection == "down"){
//     //   stat_mc = Eff_highpt_MC->GetErrorYhigh(idx);
//     //   stat_data = Eff_highpt_DATA->GetErrorYlow(idx);
//     //   total_syst_mc = sqrt(pow(stat_mc,2) + pow(tp,2));
//     //   total_syst_data = sqrt(pow(stat_data,2) + pow(tp,2));
//     //
//     //   eff_mc += Eff_highpt_MC->GetErrorYhigh(idx);
//     //   eff_data -= Eff_highpt_DATA->GetErrorYlow(idx);
//     // }
//   }
//
//   //Scale weight by (eff_data) / (eff_mc)
//   double SF = eff_data/eff_mc;
//   event.weight *= SF;
//
// return true;
// }
