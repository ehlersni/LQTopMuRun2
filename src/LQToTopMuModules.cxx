#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/LQTopMuRun2/include/LQToTopMuModules.h"
#include "UHH2/common/include/JetIds.h"
#include "UHH2/common/include/JetCorrections.h"
#include <TH1D.h>


using namespace uhh2;
using namespace std;


// WeightsTo14TeV::WeightsTo14TeV(uhh2::Context & ctx, TString pdfname){
//
//   if( ( gSystem->Load("libLHAPDF") )==-1){
//     std::cerr << "libLHAPDF not found, no pdf weights will be applied. To apply pdf re-weighting, add path to libLHAPDF.so to LD_LIBRARY_PATH" << std::endl;
//     m_libvalid=false;
//     return;
//   }
//   m_libvalid=true;
//
//   LHAPDF::initPDFSet(1, (std::string)(pdfname+".LHgrid"));
//
//   pdf = LHAPDF::mkPDF( (std::string) pdfname, 0);
//
//   xmin = LHAPDF::getXmin(0);
//   xmax = LHAPDF::getXmax(0);
//   qmin = sqrt(LHAPDF::getQ2min(0));
//   qmax = sqrt(LHAPDF::getQ2max(0));
//
//   h_x13_1 = ctx.get_handle<double>("x13_1");
//   h_x13_2 = ctx.get_handle<double>("x13_2");
//   h_x14_1 = ctx.get_handle<double>("x14_1");
//   h_x14_2 = ctx.get_handle<double>("x14_2");
//   h_Q = ctx.get_handle<double>("Q");
//   h_xf13_1 = ctx.get_handle<double>("xf13_1");
//   h_xf13_2 = ctx.get_handle<double>("xf13_1");
//   h_xf14_1 = ctx.get_handle<double>("xf14_1");
//   h_xf14_2 = ctx.get_handle<double>("xf14_2");
//   h_f1 = ctx.get_handle<int>("f1");
//   h_f2 = ctx.get_handle<int>("f2");
//   h_weight13 = ctx.get_handle<double>("weight13");
//   h_weight14 = ctx.get_handle<double>("weight14");
//   h_sf = ctx.get_handle<double>("sf");
//
// }
//
// double WeightsTo14TeV::calculateWeight(uhh2::Event & event){
//   if(!m_libvalid) return -1.;
//
//   double x1=event.genInfo->pdf_x1();
//   double x2=event.genInfo->pdf_x2();
//   if(x1 <= 0 || x2 <= 0) throw runtime_error("in LQToTopMuModules.cxx, WeightsTo14TeV::calculateWeight: One of the 13 TeV Bjorken-x is <= 0.");
//
//   int id1 = event.genInfo->pdf_id1();
//   int id2 = event.genInfo->pdf_id2();
//   if(id1 == 21 || id1 == 9) id1 = 0;
//   if(id2 == 21 || id2 == 9) id2 = 0;
//
//   double q = event.genInfo->pdf_scalePDF();
//
//   // Enforce minimum values of q and x
//   if(q < qmin) q = qmin;
//   if(q > qmax) q = qmax;
//   if(x1 < xmin) x1 = xmin;
//   if(x1 > xmax) x1 = xmax;
//   if(x2 < xmin) x2 = xmin;
//   if(x2 > xmax) x2 = xmax;
//
//
//   double xpdf1 = pdf->xfxQ(id1, x1, q);
//   double xpdf2 = pdf->xfxQ(id2, x2, q);
//   //cout <<"x1: " << x1 << ", x2: " << x2 << ", q: " << q << ", xpdf1: " << xpdf1 << ", xpdf2: " << xpdf2 << endl;
//
//   //new values of x can be smaller to achieve the same partonic sqrt(s)
//   double x1_14tev = x1 * (13./14.);
//   double x2_14tev = x2 * (13./14.);
//   if(x1_14tev < xmin) x1_14tev = xmin;
//   if(x2_14tev < xmin) x2_14tev = xmin;
//
//   double xpdf1_14tev = pdf->xfxQ(id1, x1_14tev, q);
//   double xpdf2_14tev = pdf->xfxQ(id2, x2_14tev, q);
//   //cout <<"x1_14tev: " << x1_14tev << ", x2_tev: " << x2_14tev << ", q: " << q << ", xpdf1_tev: " << xpdf1_14tev << ", xpdf2_tev: " << xpdf2_14tev << endl;
//
//   double w0 = xpdf1 * xpdf2;
//   if(w0 == 0) return -1.;
//
//   double weight = xpdf1_14tev * xpdf2_14tev / w0;
//   if(weight == 0) return -1.;
//
//   //cout << "Old weight: " << w0 << ", new weight: " << xpdf1_14tev * xpdf2_14tev << ", SF: " << weight << endl << endl;
//
//   //Set handles
//   event.set(h_x13_1, x1);
//   event.set(h_x13_2, x2);
//   event.set(h_x14_1, x1_14tev);
//   event.set(h_x14_2, x2_14tev);
//   event.set(h_Q, q);
//   event.set(h_xf13_1, xpdf1);
//   event.set(h_xf13_2, xpdf2);
//   event.set(h_xf14_1, xpdf1_14tev);
//   event.set(h_xf14_2, xpdf2_14tev);
//   event.set(h_f1, id1);
//   event.set(h_f2, id2);
//   event.set(h_weight13, w0);
//   event.set(h_weight14, xpdf1_14tev * xpdf2_14tev);
//   event.set(h_sf, weight);
//
//   // This weight has to be applied to the event.weight
//   return weight;
//
// }



JetLeptonOverlapCleaner::JetLeptonOverlapCleaner(double RJet_): RJet(RJet_){}

bool JetLeptonOverlapCleaner::process(Event & event){

   vector<Jet> result_jets;
   vector<Muon> result_muons;
   vector<Electron> result_electrons;

   cout << "############# Start of a new event ##################" << endl << endl;

   cout << "------------- first deal with the jets --------------" << endl << endl;
   //Check, if jets are identical to leptons
   //if so, kick them out
   for(const Jet & thisjet : *event.jets){
     auto thisjet_v4_raw = thisjet.v4() * thisjet.JEC_factor_raw();
     bool keep_thisjet = true;
     cout << "new jet: " << endl;
     for(const Muon & thismu : *event.muons){
       //if true, jet is considered identical to the muon
       cout << "deltaR between this muon and the jet: " << deltaR(thismu,thisjet) << ", mupt/jetptraw: " << thismu.pt()/thisjet_v4_raw.pt() << endl;
       if (deltaR(thismu,thisjet) <= 0.05 && (thismu.pt() <= thisjet_v4_raw.pt()*1.1 && thismu.pt() >= thisjet_v4_raw.pt()*0.8)){
	 keep_thisjet = false;
	 cout << "at least bc/ of this, this jet will not be kept" << endl;
       }
     }//end muons
     for(const Electron & thisele : *event.electrons){
       //if true, jet is considered identical to the electron
      cout << "deltaR between this electron and the jet: " << deltaR(thisele,thisjet) << ", elept/jetptraw: " << thisele.pt()/thisjet_v4_raw.pt() << endl;
       if (deltaR(thisele,thisjet) <= 0.05 && (thisele.pt() <= thisjet_v4_raw.pt()*1.1 && thisele.pt() >= thisjet_v4_raw.pt()*0.8)){
	 keep_thisjet = false;
	 cout << "at least bc/ of this, this jet will not be kept" << endl;
       }
     }//end electrons
     if(keep_thisjet) result_jets.push_back(thisjet);
   }
   //final collection of jets that are not identical to leptons
   swap(*event.jets,result_jets);



   //now check, if leptons are lying inside a jet
   //if so, kick the lepton out

   //first for muons
   cout << "----------------- coming to the muons ------------" << endl << endl;
   for(const Muon & thismu : *event.muons){
     bool keep_thismu = true;
     cout << "new muon: " << endl;
     //already kicked out jets that actually were leptons
     for(const Jet & thisjet : *event.jets){
       cout << "deltaR between this muon and a jet: " <<  deltaR(thismu,thisjet) << endl;
       if(deltaR(thismu,thisjet) <= RJet) {
	 keep_thismu = false;
	 cout << "at least bc/ of this, this muon will not be kept" << endl;
       }
     }
     if(keep_thismu) result_muons.push_back(thismu);
   }

   //the same for electrons
   cout << "----------------- coming to the electrons ------------" << endl << endl;
   for(const Electron & thisele : *event.electrons){
     bool keep_thisele = true;
     cout << "new electron: " << endl;
     //already kicked out jets that actually were leptons
     for(const Jet & thisjet : *event.jets){
       cout << "deltaR between this electron and a jet: " <<  deltaR(thisele,thisjet) << endl;
       if(deltaR(thisele,thisjet) <= RJet){
	 keep_thisele = false;
	 cout << "at least bc/ of this, this muon will not be kept" << endl;
       }
     }
     if(keep_thisele) result_electrons.push_back(thisele);
   }

   swap(*event.muons,result_muons);
   swap(*event.electrons,result_electrons);

   cout << "############### End of the event ###################" << endl << endl;
   return true;
}


ElectronTriggerWeights::ElectronTriggerWeights(Context & ctx, TString path_, TString SysDirection_): path(path_), SysDirection(SysDirection_){

  auto dataset_type = ctx.get("dataset_type");
  bool is_mc = dataset_type == "MC";
  if(!is_mc){
    cout << "Warning: ElectronTriggerWeights will not have an effect on this non-MC sample (dataset_type = '" + dataset_type + "')" << endl;
    return;
  }

  unique_ptr<TFile> file;
  file.reset(new TFile(path,"READ"));

  Eff_lowpt_MC.reset((TGraphAsymmErrors*)file->Get("gr_lowpt_eta_TTbar_eff"));
  Eff_highpt_MC.reset((TGraphAsymmErrors*)file->Get("gr_highpt_eta_TTbar_eff"));
  Eff_lowpt_DATA.reset((TGraphAsymmErrors*)file->Get("gr_lowpt_eta_DATA_eff"));
  Eff_highpt_DATA.reset((TGraphAsymmErrors*)file->Get("gr_highpt_eta_DATA_eff"));

  if(SysDirection != "nominal" && SysDirection != "up" && SysDirection != "down") throw runtime_error("In LQToTopMuModules.cxx, ElectronTriggerWeights.process(): Invalid SysDirection specified.");

}

bool ElectronTriggerWeights::process(Event & event){

  if(event.isRealData) return true;

  const auto ele = event.electrons->at(0);
  double eta = ele.eta();
  if(fabs(eta) > 2.4) throw runtime_error("In LQToTopMuModules.cxx, ElectronTriggerWeights.process(): Ele-|eta| > 2.4 is not supported at the moment.");


  //find right bin in eta
  int idx = 0;
  bool lowpt = false;
  if(30 <= ele.pt() && ele.pt() < 120){
    lowpt = true;
    //lowpt trigger
    bool keep_going = true;
    while(keep_going){
      double x,y;
      Eff_lowpt_MC->GetPoint(idx,x,y);
      keep_going = eta > x + Eff_lowpt_MC->GetErrorXhigh(idx);
      if(keep_going) idx++;
    }
  }
  else if(ele.pt() >= 120){
    //highpt trigger
    bool keep_going = true;
    while(keep_going){
      double x,y;
      Eff_highpt_MC->GetPoint(idx,x,y);
      keep_going = eta > x + Eff_highpt_MC->GetErrorXhigh(idx);
      if(keep_going) idx++;
    }
  }
  else throw runtime_error("In LQToTopMuModules.cxx, ElectronTriggerWeights.process(): Electron has pt<30. Clean electron collection before applying weights.");

  //access efficiencies for MC and DATA, possibly accout for systematics = statistical + add. 2% up/down
  double eff_data = -1, eff_mc = -1, dummy_x;
  double stat_data = -1, stat_mc = -1, tp = 0.02, total_syst_data = -1, total_syst_mc = -1;
  if(lowpt){
    Eff_lowpt_MC->GetPoint(idx,dummy_x,eff_mc);
    Eff_lowpt_DATA->GetPoint(idx,dummy_x,eff_data);

    if(SysDirection == "up"){
      stat_mc = Eff_lowpt_MC->GetErrorYlow(idx);
      stat_data = Eff_lowpt_DATA->GetErrorYhigh(idx);
      total_syst_mc = sqrt(pow(stat_mc,2) + pow(tp,2));
      total_syst_data = sqrt(pow(stat_data,2) + pow(tp,2));

      eff_mc -= total_syst_mc;
      eff_data += total_syst_data;
    }
    else if(SysDirection == "down"){
      stat_mc = Eff_lowpt_MC->GetErrorYhigh(idx);
      stat_data = Eff_lowpt_DATA->GetErrorYlow(idx);
      total_syst_mc = sqrt(pow(stat_mc,2) + pow(tp,2));
      total_syst_data = sqrt(pow(stat_data,2) + pow(tp,2));

      eff_mc += Eff_lowpt_MC->GetErrorYhigh(idx);
      eff_data -= Eff_lowpt_DATA->GetErrorYlow(idx);
    }
  }
  else{
    Eff_highpt_MC->GetPoint(idx,dummy_x,eff_mc);
    Eff_highpt_DATA->GetPoint(idx,dummy_x,eff_data);

    if(SysDirection == "up"){
      stat_mc = Eff_highpt_MC->GetErrorYlow(idx);
      stat_data = Eff_highpt_DATA->GetErrorYhigh(idx);
      total_syst_mc = sqrt(pow(stat_mc,2) + pow(tp,2));
      total_syst_data = sqrt(pow(stat_data,2) + pow(tp,2));

      eff_mc -= Eff_highpt_MC->GetErrorYlow(idx);
      eff_data += Eff_highpt_DATA->GetErrorYhigh(idx);
    }
    else if(SysDirection == "down"){
      stat_mc = Eff_highpt_MC->GetErrorYhigh(idx);
      stat_data = Eff_highpt_DATA->GetErrorYlow(idx);
      total_syst_mc = sqrt(pow(stat_mc,2) + pow(tp,2));
      total_syst_data = sqrt(pow(stat_data,2) + pow(tp,2));

      eff_mc += Eff_highpt_MC->GetErrorYhigh(idx);
      eff_data -= Eff_highpt_DATA->GetErrorYlow(idx);
    }
  }

  //Scale weight by (eff_data) / (eff_mc)
  double SF = eff_data/eff_mc;
  event.weight *= SF;

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


ElectronJetOverlapCleaner::ElectronJetOverlapCleaner(){}

bool ElectronJetOverlapCleaner::process(Event & event, int idx1, int idx2){
  if(event.electrons){
    if(event.electrons->size() < 2){
      cout << "This event does not contain >= 2 electrons, no ElectronJetOverlapCleaning is applied." << endl;
      return false;
    }
  }

  if(idx1 < 0 || idx2 < 0){
    throw runtime_error("ElectronJetOverlapCleaner was not given two valid electron indices.");
    return false;
  }

  vector<Jet> event_jets_new;
  int Njets = event.jets->size();
  for(int i=0; i<Njets; i++){
    double dR1 = deltaR(event.jets->at(i),event.electrons->at(idx1));
    double dR2 = deltaR(event.jets->at(i),event.electrons->at(idx2));
    if(dR1 > 0.1 && dR2 > 0.1) event_jets_new.push_back(event.jets->at(i));
    //cout << "In ElectronJet-Cleaner: This jet has (dR1, dR2): (" << dR1 << ", " << dR2 << "), it has electron-multiplicity: " << event.jets->at(i).electronMultiplicity() << endl;
  }
  if(fabs(Njets - event_jets_new.size()) > 2) throw runtime_error("In LQToTopMuModules.cxx::ElectronJetOverlapCleaner::process(): When deleting jets overlapping with best electrons, more than 2 were deleted. This must not be.");
  swap(event_jets_new, *event.jets);
  return true;
}

DibosonScaleFactors::DibosonScaleFactors(Context & ctx, TString path_, TString SysDirectionXSec_, TString SysDirectionBTag_) : path(path_), SysDirectionXSec(SysDirectionXSec_), SysDirectionBTag(SysDirectionBTag_){

  auto dataset_type = ctx.get("dataset_type");
  bool is_mc = dataset_type == "MC";
  if(!is_mc){
    cout << "Warning: ElectronTriggerWeights will not have an effect on this non-MC sample (dataset_type = '" + dataset_type + "')" << endl;
    return;
  }
  TString dataset_version = ctx.get("dataset_version");
  if(!dataset_version.Contains("Diboson")){
    is_diboson = false;
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

  unique_ptr<TH1D> XSecSF, BTagSF;
  XSecSF.reset((TH1D*)file->Get("Diboson_XSec_SF"));
  BTagSF.reset((TH1D*)file->Get("Diboson_BTag_SF"));

  double xsec_SF, btag1_SF, btag2_SF;
  if(SysDirectionXSec == "nominal")   xsec_SF = XSecSF->GetBinContent(1);
  else if(SysDirectionXSec == "up")   xsec_SF = XSecSF->GetBinContent(1) + XSecSF->GetBinError(1);
  else if(SysDirectionXSec == "down") xsec_SF = XSecSF->GetBinContent(1) - XSecSF->GetBinError(1);
  else throw runtime_error("In DibosonScaleFactors::process(): Invalid SysDirectionXSec specified.");

  if(SysDirectionBTag == "nominal"){
    btag1_SF = BTagSF->GetBinContent(1);
    btag2_SF = BTagSF->GetBinContent(2);
  }
  else if(SysDirectionBTag == "up"){
    btag1_SF = BTagSF->GetBinContent(1) + BTagSF->GetBinError(1);
    btag2_SF = BTagSF->GetBinContent(2) + BTagSF->GetBinError(2);
  }
  else if(SysDirectionBTag == "down"){
    btag1_SF = BTagSF->GetBinContent(1) - BTagSF->GetBinError(1);
    btag2_SF = BTagSF->GetBinContent(2) - BTagSF->GetBinError(2);
  }
  else throw runtime_error("In DibosonScaleFactors::process(): Invalid SysDirectionBTag specified.");

  //First apply XSec SF to all Diboson events before applying 1&2-btag SF
  //cout << "Weight before xsec SF: " << event.weight << ", SF: " << xsec_SF << endl;
  event.weight *= xsec_SF;
  //cout << "Weight after xsec SF: " << event.weight << endl;

  //Count number of loose b-jets in the event
  int n_bjets = 0;
  CSVBTag Btag_loose = CSVBTag(CSVBTag::WP_LOOSE);
  for (unsigned int i =0; i<event.jets->size(); ++i) {
    if(Btag_loose(event.jets->at(i),event)){
      n_bjets++;
    }
  }

  //cout << "Number of btags in the event: " << n_bjets << endl;
  //cout << "SF1: " << btag1_SF << ", SF2: " << btag2_SF << endl;
  //cout << "Weight before applying BTag SF: " << event.weight << endl;
  //apply 1&2-btag SF
  if(n_bjets > 2) throw runtime_error("In DibosonScaleFactors::process(): More than 2 b-jets present in the event. Scale factors have only been derived for Nbjets <= 2.");
  else if(n_bjets == 1) event.weight *= btag1_SF;
  else if(n_bjets == 2) event.weight *= btag2_SF;
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
