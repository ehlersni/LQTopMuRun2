#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/LQTopMuRun2/include/LQTopMuRun2FakeRateModule.h"
#include "UHH2/common/include/JetIds.h"
#include "UHH2/common/include/JetCorrections.h"
#include <TH1D.h>


using namespace uhh2;
using namespace std;


JetCorrectorVariable::JetCorrectorVariable(uhh2::Context & ctx, const std::vector<std::string> & JEC_files): JetCorrector(ctx, JEC_files){}

bool JetCorrectorVariable::correct_collection(uhh2::Event & event, std::vector<Jet> & jets){

  //apply jet corrections
  for(auto & jet : jets){
    correct_jet(*corrector, jet, event, jec_uncertainty, direction);
  }
  return true;
};





ElectronFakeRateWeights::ElectronFakeRateWeights(Context & ctx, const std::vector<std::string> & JEC_files, TString path_, TString SysDirection_, const string label_jets, const string label_genjets):  path(path_), SysDirection(SysDirection_){

  cout << "Fake1" << endl;
  auto dataset_type = ctx.get("dataset_type");
  bool is_mc = dataset_type == "MC";
  if(!is_mc){
    cout << "Warning: ElectronFakeRateWeights will not have an effect on this non-MC sample (dataset_type = '" + dataset_type + "')" << endl;
    return;
  }

  unique_ptr<TFile> file;
  file.reset(new TFile(path,"READ"));

  SF.reset((TGraphAsymmErrors*)file->Get("ScaleFactors"));
  //Read the SF file and store all SFs and their pt-range
  n_points = SF->GetN();
  for(int i=0; i<n_points; i++){
    double x,y;
    SF->GetPoint(i,x,y);
    x_low.emplace_back(x-SF->GetErrorXlow(i));
    x_high.emplace_back(x+SF->GetErrorXhigh(i));
    cout << "Fake2" << endl;
  }

  if(SysDirection != "nominal" && SysDirection != "up" && SysDirection != "down") throw runtime_error("In LQToTopMuModules.cxx, ElectronFakeRateWeights.process(): Invalid SysDirection specified.");

  jet_corrector.reset(new JetCorrectorVariable(ctx, JEC_files));

  bool jer_was_applied = ctx.get("meta_jer_applied", "") == "true";
  if(jer_was_applied) ctx.set_metadata("jer_applied", "false", true);
  //jet_smearer.reset(new JetSmearerVariable(ctx, label_jets, label_genjets, JERSmearing::SF_13TeV_2016));
  jet_smearer.reset(new GenericJetResolutionSmearer(ctx, label_jets, label_genjets, true, JERSmearing::SF_13TeV_2016));
  if(!jer_was_applied) ctx.set_metadata("jer_applied", "false", true);

  jet_id = AndId<Jet>(JetPFID(JetPFID::WP_LOOSE), PtEtaCut(30.0, 2.4));

  FakeRateWeightEle = ctx.get_handle<double>("FakeRateWeightEle");
  FakeRateWeightEleUp = ctx.get_handle<double>("FakeRateWeightEleUp");
  FakeRateWeightEleDown = ctx.get_handle<double>("FakeRateWeightEleDown");
  h_jets = ctx.get_handle<std::vector<Jet>>(label_jets);
  cout << "Fake3" << endl;
}

bool ElectronFakeRateWeights::process(Event & event){


  vector<Jet> *jets = &event.get(h_jets);
  jet_corrector->correct_collection(event, *jets);
  jet_smearer->process(event);

  if(event.isRealData || event.electrons->size() < 1){
    event.set(FakeRateWeightEle,1.);
    event.set(FakeRateWeightEleUp,1.);
    event.set(FakeRateWeightEleDown,1.);
    return false;
  }

  //find fake-electrons
  vector<bool> is_fake;

  //if the number of gen and reco-muons are the same, assume muons are real
  unsigned int n_genele = 0;
  for(const auto & gp : *event.genparticles){
    if(fabs(gp.pdgId()) != 11) continue;
    n_genele++;
  }
  //no fakes if number of gen and reco electrons is the same
  if(n_genele == event.electrons->size()){
    event.set(FakeRateWeightEle,1.);
    event.set(FakeRateWeightEleUp,1.);
    event.set(FakeRateWeightEleDown,1.);
    return false;
  }
  else{
    //if ngen and nreco are unequal, try to match muons to muons within 0.1 and to taus within 0.2
    unsigned int n_matched_to_ele = 0, n_matched_to_tau = 0, n_matched_to_b = 0, n_matched_to_c = 0;
    for(const auto & ele : *event.electrons){
      bool is_matched = false;
      for(const auto & gp : *event.genparticles){
        if(fabs(gp.pdgId()) == 11){
          if(deltaR(gp,ele) < 0.1 && !is_matched){
            is_matched = true;
            n_matched_to_ele++;
          }
        }
        else if(fabs(gp.pdgId()) == 15){
          if(deltaR(gp,ele) < 0.2 && !is_matched){
            is_matched = true;
            n_matched_to_tau++;
          }
        }
        else if(fabs(gp.pdgId()) == 5){
          if(deltaR(gp,ele) < 0.2 && !is_matched){
            is_matched = true;
            n_matched_to_b++;
          }
        }
        else if(fabs(gp.pdgId()) == 4){
          if(deltaR(gp,ele) < 0.2 && !is_matched){
            is_matched = true;
            n_matched_to_c++;
          }
        }
      }
      is_fake.push_back(!is_matched);
    }
    if(n_matched_to_tau + n_matched_to_ele + n_matched_to_b + n_matched_to_c == event.electrons->size()){
      event.set(FakeRateWeightEle,1.);
      event.set(FakeRateWeightEleUp,1.);
      event.set(FakeRateWeightEleDown,1.);
      return false;
    }
  }

  //find jets matching fake-electrons
  vector<bool> faking_jet;
  for(auto & jet : *jets){

    bool faking = false;
    double dr_min = 999;
    for(unsigned int i=0; i<event.electrons->size(); i++){

      if(!is_fake[i]) continue;
      //consider only jets to be eligible for faking the electron that fulfill the jetId used when deriving the scale factors
      if(!jet_id(jet,event)) continue;

      if(deltaR(event.electrons->at(i), jet) < 0.4) faking = true;
      if(deltaR(event.electrons->at(i), jet) < dr_min) dr_min = deltaR(event.electrons->at(i), jet);
    }
    faking_jet.push_back(faking);
  }


  //apply SF to the eventweight for each faking jet
  double SF_final = 1, SF_final_up = 1, SF_final_down = 1;
  for(unsigned int i=0; i<jets->size(); i++){
    if(!faking_jet[i]) continue;
    double pt = event.jets->at(i).pt();

    //find right bin
    int bin = -1;
    for(int j=0; j<n_points; j++){
      if(pt >= x_low[j] && pt < x_high[j]) bin = j;
    }

    //possibly accout for systematics = statistical, inflate stat unc. by 50% if >800 GeV
    double x,scale_factor, scale_factor_up, scale_factor_down, mult = 1;
    if(pt > x_high[n_points-1]){
      bin = n_points - 1;
      mult = 1.5;
    }
    else if(bin == -1) throw runtime_error("In ElectronFakeRateWeights::process: Did not find right SF-bin for this jet.");

    SF->GetPoint(bin,x,scale_factor);
    scale_factor_up   = scale_factor + mult * SF->GetErrorYhigh(bin);
    scale_factor_down = scale_factor - mult * SF->GetErrorYlow(bin);

    if(SysDirection == "up") event.weight *= scale_factor_up;
    else if(SysDirection == "down") event.weight *= scale_factor_down;
    else if(SysDirection == "nominal") event.weight *= scale_factor;
    else throw runtime_error("In ElectronFakeRateWeights::process(): SysDirection is not one of the following: ['up', 'down', 'nominal']");

    SF_final      *= scale_factor;
    SF_final_up   *= scale_factor_up;
    SF_final_down *= scale_factor_down;

  }

  event.set(FakeRateWeightEle,SF_final);
  event.set(FakeRateWeightEleUp,SF_final_up);
  event.set(FakeRateWeightEleDown,SF_final_down);


  return true;
}


MuonFakeRateWeights::MuonFakeRateWeights(Context & ctx, TString path_, TString SysDirection_):  path(path_), SysDirection(SysDirection_){

  auto dataset_type = ctx.get("dataset_type");
  bool is_mc = dataset_type == "MC";
  if(!is_mc){
    cout << "Warning: MuonFakeRateWeights will not have an effect on this non-MC sample (dataset_type = '" + dataset_type + "')" << endl;
    return;
  }

  unique_ptr<TFile> file;
  file.reset(new TFile(path,"READ"));
  SF.reset((TGraphAsymmErrors*)file->Get("ScaleFactors"));

  if(SysDirection != "nominal" && SysDirection != "up" && SysDirection != "down") throw runtime_error("In LQToTopMuModules.cxx, MuonFakeRateWeights.process(): Invalid SysDirection specified.");


  FakeRateWeightMu = ctx.get_handle<double>("FakeRateWeightMu");
  FakeRateWeightMuUp = ctx.get_handle<double>("FakeRateWeightMuUp");
  FakeRateWeightMuDown = ctx.get_handle<double>("FakeRateWeightMuDown");

}

bool MuonFakeRateWeights::process(Event & event){

  if(event.isRealData || event.muons->size() < 1){
    event.set(FakeRateWeightMu,1.);
    event.set(FakeRateWeightMuUp,1.);
    event.set(FakeRateWeightMuDown,1.);
    return false;
  }

  //find fake-muons
  vector<bool> is_fake;
  unsigned int n_genmu = 0;
  for(const auto & gp : *event.genparticles){
    if(fabs(gp.pdgId()) != 13) continue;
    n_genmu++;
  }

  if(n_genmu == event.muons->size()){
    event.set(FakeRateWeightMu,1.);
    event.set(FakeRateWeightMuUp,1.);
    event.set(FakeRateWeightMuDown,1.);
    return false;
  }
  else{
    //if ngen and nreco are unequal, try to match muons to muons within 0.1 and to taus within 0.2
    unsigned int n_matched_to_muons = 0, n_matched_to_taus = 0, n_matched_to_b = 0, n_matched_to_c = 0;
    for(const auto & mu : *event.muons){
      bool is_matched = false;
      for(const auto & gp : *event.genparticles){
        if(fabs(gp.pdgId()) == 13){
          if(deltaR(gp,mu) < 0.1 && !is_matched){
            is_matched = true;
            n_matched_to_muons++;
          }
        }
        else if(fabs(gp.pdgId()) == 15){
          if(deltaR(gp,mu) < 0.2 && !is_matched){
            is_matched = true;
            n_matched_to_taus++;
          }
        }
        else if(fabs(gp.pdgId()) == 5){
          if(deltaR(gp,mu) < 0.2 && !is_matched){
            is_matched = true;
            n_matched_to_b++;
          }
        }
        else if(fabs(gp.pdgId()) == 4){
          if(deltaR(gp,mu) < 0.2 && !is_matched){
            is_matched = true;
            n_matched_to_c++;
          }
        }
      }
      is_fake.push_back(!is_matched);
    }
    if(n_matched_to_taus + n_matched_to_muons + n_matched_to_b + n_matched_to_c == event.muons->size()){
      event.set(FakeRateWeightMu,1.);
      event.set(FakeRateWeightMuUp,1.);
      event.set(FakeRateWeightMuDown,1.);
      return false;
    }
  }


  //apply SF to the eventweight for each fake muon
  double SF_final = 1, SF_final_up = 1, SF_final_down = 1;
  for(unsigned int i=0; i<event.muons->size(); i++){
    if(!is_fake[i]) continue;
    int bin = 0;

    double x,scale_factor, scale_factor_up, scale_factor_down;
    SF->GetPoint(bin,x,scale_factor);
    scale_factor_up   = scale_factor + SF->GetErrorYhigh(bin);
    scale_factor_down = scale_factor - SF->GetErrorYlow(bin);

    if(SysDirection == "up") event.weight *= scale_factor_up;
    else if(SysDirection == "down") event.weight *= scale_factor_down;
    else if(SysDirection == "nominal") event.weight *= scale_factor;
    else throw runtime_error("In MuonFakeRateWeights::process(): SysDirection is not one of the following: ['up', 'down', 'nominal']");

    SF_final      *= scale_factor;
    SF_final_up   *= scale_factor_up;
    SF_final_down *= scale_factor_down;

  }

  event.set(FakeRateWeightMu,SF_final);
  event.set(FakeRateWeightMuUp,SF_final_up);
  event.set(FakeRateWeightMuDown,SF_final_down);


  return true;
}
