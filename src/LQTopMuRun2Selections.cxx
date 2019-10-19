#include "UHH2/LQTopMuRun2/include/LQTopMuRun2Selections.h"
#include "UHH2/core/include/Event.h"
#include "TFile.h"
#include <vector>
#include <TROOT.h>
#include <TH1.h>

#include <stdexcept>

using namespace uhh2examples;
using namespace uhh2;
using namespace std;


DijetSelection::DijetSelection(float dphi_min_, float third_frac_max_): dphi_min(dphi_min_), third_frac_max(third_frac_max_){}

bool DijetSelection::passes(const Event & event){
  assert(event.jets); // if this fails, it probably means jets are not read in
  if(event.jets->size() < 2) return false;
  const auto & jet0 = event.jets->at(0);
  const auto & jet1 = event.jets->at(1);
  auto dphi = deltaPhi(jet0, jet1);
  if(dphi < dphi_min) return false;
  if(event.jets->size() == 2) return true;
  const auto & jet2 = event.jets->at(2);
  auto third_jet_frac = jet2.pt() / (0.5 * (jet0.pt() + jet1.pt()));
  return third_jet_frac < third_frac_max;
}

StSelection::StSelection(Context & ctx, float ST_min_, float ST_max_): ST_min(ST_min_), ST_max(ST_max_)
{
  h_ST = ctx.get_handle<float>("ST");
}

bool StSelection::passes(const Event & event)
{
  if(!event.is_valid(h_ST)) throw runtime_error("In LQTopMuRun2Hists.cxx: Trying to fill ST before filling the handle.");
  float ST = event.get(h_ST);

  bool pass = false;
  if(ST >= ST_min && (ST <= ST_max || ST_max < 0)) pass = true;

  return pass;
}

StLeptonSelection::StLeptonSelection(Context & ctx, float STlep_min_, float STlep_max_): STlep_min(STlep_min_), STlep_max(STlep_max_)
{
  h_STlep = ctx.get_handle<float>("STlep");
}

bool StLeptonSelection::passes(const Event & event)
{
  if(!event.is_valid(h_STlep)) throw runtime_error("In LQTopMuRun2Hists.cxx: Trying to fill STlep before filling the handle.");
  float STlep = event.get(h_STlep);

  bool pass = false;
  if(STlep >= STlep_min && (STlep <= STlep_max || STlep_max < 0)) pass = true;

  return pass;
}



InvMassOf2Mu::InvMassOf2Mu(double m_min_, double m_max_):m_min(m_min_), m_max(m_max_){}
bool InvMassOf2Mu::passes(const Event & event)
{
  bool pass = true;
  int Nmuons = event.muons->size();
  double M_mumu;
  LorentzVector muons[Nmuons];
  for(int i = 0; i < Nmuons; i++)
  {
    muons[i] = event.muons->at(i).v4();
  }
  for(int i = 0; i < Nmuons; i++)
  {
    for(int j = 0; j < Nmuons; j++)
    {
      if(j > i)
      {
        M_mumu = (muons[i] + muons[j]).M(); // .M() is a function which calculates the invariant mass of a LorentzVector
        if(!(M_mumu > m_min && (M_mumu < m_max || m_max < 0.)))
        {
          pass = false;
        }
      }
    }
  }
  return pass;
}



InvMassOf2Ele::InvMassOf2Ele(double m_min_, double m_max_):m_min(m_min_), m_max(m_max_){}
bool InvMassOf2Ele::passes(const Event & event)
{
  bool pass = true;
  int Nelectrons = event.electrons->size();
  double M_eleele;
  LorentzVector electrons[Nelectrons];
  for(int i = 0; i < Nelectrons; i++)
  {
    electrons[i] = event.electrons->at(i).v4();
  }
  for(int i = 0; i < Nelectrons; i++)
  {
    for(int j = 0; j < Nelectrons; j++)
    {
      if(j > i)
      {
        M_eleele = (electrons[i] + electrons[j]).M(); // .M() is a function which calculates the invariant mass of a LorentzVector
        if(!(M_eleele > m_min && (M_eleele < m_max || m_max < 0.)))
        {
          pass = false;
        }
      }
    }
  }
  return pass;
}

InvMassEleEleSelection::InvMassEleEleSelection(double m_min_, double m_max_):m_min(m_min_), m_max(m_max_){}
bool InvMassEleEleSelection::passes(const Event & event){

  bool pass = false;
  int Nelectrons = event.electrons->size();
  double M_ee;
  LorentzVector electrons[Nelectrons];
  for(int i=0; i<Nelectrons; i++){
    electrons[i] = event.electrons->at(i).v4();
  }
  for(int i=0; i<Nelectrons; i++){
    for(int j=0; j<Nelectrons; j++){
      if(j>i){
        M_ee = (electrons[i] + electrons[j]).M();
        if(M_ee >= m_min && M_ee <= m_max){
          pass = true;
        }
      }
    }
  }
  return pass;
}

ElectronJetOverlapCleaner::ElectronJetOverlapCleaner(){}

bool ElectronJetOverlapCleaner::process(Event & event, int idx1, int idx2){
  if(event.electrons){
    if(event.electrons->size() < 2){
      cout << "This event does not contain >=2 elctrons, no ElectronJetOverlapCleaning is applied." << endl;
      return false;
    }
  }
  if (idx1 < 0 || idx2 < 0){
    throw runtime_error("ElectronJetOverlapCleaner was not given valid electron indices");
    return false;
  }

  vector<Jet> event_jets_new;
  int Njets = event.jets->size();
  for(int i=0; i < Njets; i++){
    double dR1 = deltaR(event .jets->at(i), event.electrons->at(idx1));
    double dR2 = deltaR(event .jets->at(i), event.electrons->at(idx2));
    if(dR1 > 0.1 && dR2 > 0.1) event_jets_new.push_back(event.jets->at(i));
  }
  if(fabs(Njets - event_jets_new.size()) > 2) throw runtime_error("In LQToTopMuModules.cxx::ElectronJetOverlapCleaner::process(): When deleting jets overlapping with best electrons, more than 2 were deleted. This must not be.");
  swap(event_jets_new, *event.jets);
  return true;
}

ZEEFinder::ZEEFinder(){}

pair<int,int> ZEEFinder::search(Event & event){
  if(event.electrons){
    if(event.electrons->size() < 2){
      cout << "This event does not contain >= 2 electrons, no pair can be found. Going to return (-1,-1)." << endl;
      pair<int,int> dummy(-1,-1);
      return dummy;
    }
  }

  const int Nele = event.electrons->size();
  int idx_best_ele1 = -1;
  int idx_best_ele2 = -1;
  double Mee_best = 99999999;
  for(int i=0; i<Nele; i++){
    for(int j=0; j<Nele; j++){
      if(j>i){
        double Mee = (event.electrons->at(i).v4()+event.electrons->at(j).v4()).M();
        if(fabs(Mee - 91.2) < fabs(Mee_best - 91.2)){
          Mee_best = Mee;
          idx_best_ele1 = i;
          idx_best_ele2 = j;
        }
      }
    }
  }
  pair<int,int> dummy(idx_best_ele1,idx_best_ele2);
  return dummy;
}

pair<int,int> ZEEFinder::search(const Event & event){
  if(event.electrons){
    if(event.electrons->size() < 2){
      cout << "This event does not contain >= 2 electrons, no pair can be found. Going to return (-1,-1)." << endl;
      pair<int,int> dummy(-1,-1);
      return dummy;
    }
  }

  const int Nele = event.electrons->size();
  int idx_best_ele1 = -1;
  int idx_best_ele2 = -1;
  double Mee_best = 99999999;
  for(int i=0; i<Nele; i++){
    for(int j=0; j<Nele; j++){
      if(j>i){
        double Mee = (event.electrons->at(i).v4()+event.electrons->at(j).v4()).M();
        if(fabs(Mee - 91.2) < fabs(Mee_best - 91.2)){
          Mee_best = Mee;
          idx_best_ele1 = i;
          idx_best_ele2 = j;
        }
      }
    }
  }
  pair<int,int> dummy(idx_best_ele1,idx_best_ele2);
  return dummy;
}

DibosonScaleFactors::DibosonScaleFactors(Context & ctx, TString path_, TString SysDirectionXSec_) : path(path_), SysDirectionXSec(SysDirectionXSec_){ //, TString SysDirectionBTag_ , SysDirectionBTag(SysDirectionBTag_){

  auto dataset_type = ctx.get("dataset_type");
  bool is_mc = dataset_type == "MC";
  if(!is_mc){
    cout << 1 << endl;
    cout << "Warning: ElectronTriggerWeights will not have an effect on this non-MC sample (dataset_type = '" + dataset_type + "')" << endl;
    return;
  }
  TString dataset_version = ctx.get("dataset_version");
  if(!dataset_version.Contains("Diboson")){
    is_diboson = false;
    cout << 2 << endl;
    cout << "DibosonScaleFactors will not have an effect on this non-Diboson sample (dataset_version = '" + dataset_version + "')" << endl;
    return;
  }
  else is_diboson = true;
}

bool DibosonScaleFactors::process(Event & event){
  if(event.isRealData) return true;
  if(!is_diboson) return true;

  unique_ptr<TFile> file;
  file.reset(new TFile(path,"READ"));

  unique_ptr<TH1D> XSecSF;//, BTagSF;
  XSecSF.reset((TH1D*)file->Get("Diboson_XSec_SF"));
  // BTagSF.reset((TH1D*)file->Get("Diboson_BTag_SF"));

  double xsec_SF; //, btag1_SF, btag2_SF;
  if(SysDirectionXSec == "nominal"){
    xsec_SF = XSecSF->GetBinContent(1);
  }
  else if(SysDirectionXSec == "up"){
    xsec_SF = XSecSF->GetBinContent(1) + XSecSF->GetBinError(1);
  }
  else if(SysDirectionXSec == "down"){
    xsec_SF = XSecSF->GetBinContent(1) - XSecSF->GetBinError(1);
  }
  else throw runtime_error("In DibosonScaleFactors::process(): Invalid SysDirectionXSec specified.");

  // if(SysDirectionBTag == "nominal"){
  //   btag1_SF = BTagSF->GetBinContent(1);
  //   btag2_SF = BTagSF->GetBinContent(2);
  // }
  // else if(SysDirectionBTag == "up"){
  //   btag1_SF = BTagSF->GetBinContent(1) + BTagSF->GetBinError(1);
  //   btag2_SF = BTagSF->GetBinContent(2) + BTagSF->GetBinError(2);
  // }
  // else if(SysDirectionBTag == "down"){
  //   btag1_SF = BTagSF->GetBinContent(1) - BTagSF->GetBinError(1);
  //   btag2_SF = BTagSF->GetBinContent(2) - BTagSF->GetBinError(2);
  // }
  // else throw runtime_error("In DibosonScaleFactors::process(): Invalid SysDirectionBTag specified.");
  //
  //First apply XSec SF to all Diboson events before applying 1&2-btag SF
  //cout << "Weight before xsec SF: " << event.weight << ", SF: " << xsec_SF << endl;
  event.weight *= xsec_SF;
  //cout << "Weight after xsec SF: " << event.weight << endl;

  //Count number of loose b-jets in the event
  // int n_bjets = 0;
  // CSVBTag Btag_loose = CSVBTag(CSVBTag::WP_LOOSE);
  // for (unsigned int i =0; i<event.jets->size(); ++i) {
  //   if(Btag_loose(event.jets->at(i),event)){
  //     n_bjets++;
  //   }
  // }
  //
  // //cout << "Number of btags in the event: " << n_bjets << endl;
  // //cout << "SF1: " << btag1_SF << ", SF2: " << btag2_SF << endl;
  // //cout << "Weight before applying BTag SF: " << event.weight << endl;
  // //apply 1&2-btag SF
  // if(n_bjets > 2) throw runtime_error("In DibosonScaleFactors::process(): More than 2 b-jets present in the event. Scale factors have only been derived for Nbjets <= 2.");
  // else if(n_bjets == 1) event.weight *= btag1_SF;
  // else if(n_bjets == 2) event.weight *= btag2_SF;
  //cout << "Weight after applying BTag SF: " << event.weight << endl;

  return true;
}

MuonTrkWeights::MuonTrkWeights(Context & ctx, TString path_, TString SysDirection_): path(path_), SysDirection(SysDirection_){

  auto dataset_type = ctx.get("dataset_type");
  bool is_mc = dataset_type == "MC";
  if(!is_mc){
    cout << "Warning: MuonTrkWeights will not have an effect on this non-MC sample (dataset_type = '" + dataset_type + "')" << endl;
    return;
  }

  unique_ptr<TFile> file;
  file.reset(new TFile(path,"READ"));

  Trk_SF.reset((TGraphAsymmErrors*)file->Get("ratio_eff_eta3_dr030e030_corr"));

  if(SysDirection != "nominal" && SysDirection != "up" && SysDirection != "down") throw runtime_error("In LQToTopMuModules.cxx, ElectronTriggerWeights.process(): Invalid SysDirection specified.");

}


bool MuonTrkWeights::process(Event & event){

  if(event.isRealData) return true;

  double SF = 1.0;
  for(const auto & mu : *event.muons){
    double eta = mu.eta();
    //if(fabs(eta) > 2.4) throw runtime_error("In LQToTopMuModules.cxx, MuonTrkWeights::process(): Mu-|eta| > 2.4 is not supported at the moment.");
    if(fabs(eta) > 2.4) eta = 2.39;

    //find right bin in eta
    int idx = 0;
    bool keep_going = true;
    while(keep_going){
      double x,y;
      Trk_SF->GetPoint(idx,x,y);
      keep_going = eta > x + Trk_SF->GetErrorXhigh(idx);
      if(keep_going) idx++;
    }

    double eff_fnl = 1., dummy_x;
    Trk_SF->GetPoint(idx,dummy_x,eff_fnl);

    SF *= eff_fnl;
    if(SysDirection == "up"){
      SF *= 1.005;
    }
    if(SysDirection == "down"){
      SF *= 0.995;
    }

  }

  event.weight *= SF;

  return true;

}

METSelection::METSelection(double met_min_, double met_max_) : met_min(met_min_), met_max(met_max_){}
bool METSelection::passes(const Event & event){

  auto met = event.met->pt();
  bool pass = true;

  pass = met >= met_min && (met <= met_max || met_max < 0);
  return pass;

}
