#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/CommonModules.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/JetCorrections.h"
#include "UHH2/common/include/JetHists.h"
#include "UHH2/common/include/ElectronHists.h"
#include "UHH2/common/include/ElectronIds.h"
#include "UHH2/common/include/MuonHists.h"
#include "UHH2/common/include/MuonIds.h"
#include "UHH2/common/include/EventHists.h"
#include "UHH2/common/include/LuminosityHists.h"
#include "UHH2/common/include/LumiSelection.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/LQTopMuRun2/include/LQTopMuRun2Selections.h"
#include "UHH2/LQTopMuRun2/include/LQTopMuRun2Hists.h"
#include "UHH2/common/include/PrintingModules.h"
#include "UHH2/common/include/TriggerSelection.h"
#include "UHH2/common/include/MCWeight.h"
#include "UHH2/common/include/AdditionalSelections.h"
//#include "UHH2/LQToTopMu/include/LQToTopMuModules.h"
//#include "UHH2/LQToTopMu/include/HT2dHists.h"
#include "UHH2/LQTopMuRun2/include/LQToTopMuEfficiencyHists.h"

using namespace std;
using namespace uhh2;

namespace uhh2examples {

  class LQTopMuRun2PreselectionLeptonMissidentification: public AnalysisModule {
  public:

    explicit LQTopMuRun2PreselectionLeptonMissidentification(Context & ctx);
    virtual bool process(Event & event) override;

  private:

    unique_ptr<CommonModules> common;
    unique_ptr<JetCleaner> jetcleaner;
    unique_ptr<MuonCleaner> muoncleaner;

    // declare the Selections to use.
    unique_ptr<Selection> muon_selection, nele_selection, lumi_selection, m_ee_selection; //, st_selection; //njet_selection, jet_selection,

    // store the Hists collection as member variables.
    unique_ptr<Hists> h_nocuts, h_jets_nocuts, h_ele_nocuts, h_mu_nocuts, h_event_nocuts, h_topjets_nocuts, h_lumi_nocuts;
    unique_ptr<Hists> h_cleaner, h_jets_cleaner, h_ele_cleaner, h_mu_cleaner, h_event_cleaner, h_topjets_cleaner, h_lumi_cleaner;
    unique_ptr<Hists> h_nmu, h_jets_nmu, h_ele_nmu, h_mu_nmu, h_event_nmu, h_topjets_nmu, h_lumi_nmu;
    //unique_ptr<Hists> h_njets, h_jets_njets, h_ele_njets, h_mu_njets, h_event_njets, h_topjets_njets, h_lumi_njets;

    // unique_ptr<Hists> h_2jets, h_jets_2jets, h_ele_2jets, h_mu_2jets, h_event_2jets, h_topjets_2jets, h_lumi_2jets;

    unique_ptr<Hists> h_2ele, h_jets_2ele, h_ele_2ele, h_mu_2ele, h_event_2ele, h_topjets_2ele, h_lumi_2ele;
    unique_ptr<Hists> h_mee, h_jets_mee, h_ele_mee, h_mu_mee, h_event_mee, h_topjets_mee, h_lumi_mee;
    //unique_ptr<Hists> h_stselec, h_jets_stselec, h_ele_stselec, h_mu_stselec, h_event_stselec, h_topjets_stselec, h_lumi_stselec;



    MuonId MuId;
    ElectronId ElectronID;
    // JetId JetID;

    //uhh2::Event::Handle<float> h_ST;
    uhh2::Event::Handle<bool> h_is_mlq_reconstructed;
    uhh2::Event::Handle<TString> h_mlq_reco_mode;

    bool is_mc, is_ele_channel;


  };

  LQTopMuRun2PreselectionLeptonMissidentification::LQTopMuRun2PreselectionLeptonMissidentification(Context & ctx){

    cout << "Hello from LQTopMuRun2PreselectionLeptonMissidentification" << endl;

    for (auto & kv : ctx.get_all()) {
      cout << " " << kv.first << " = " << kv.second << endl;
    }

    // ElectronID = AndId<Electron>(ElectronID_Spring16_loose, PtEtaCut(30.0, 2.4));
    Year year = extract_year(ctx);
    if (year == Year::is2016v3) {
      ElectronID = AndId<Electron>(ElectronID_Summer16_loose, PtEtaCut(30.0, 2.4));
    }
    else if (year == Year::is2017v1) {
      ElectronID = AndId<Electron>(ElectronID_Fall17_loose, PtEtaCut(30.0, 2.4));
    }
    else if (year == Year::is2018) {
      ElectronID = AndId<Electron>(ElectronID_Fall17_loose, PtEtaCut(30.0, 2.4));
    }




    // MuonID = AndId<Muon>(MuonIDTight(), PtEtaCut(30.0, 2.4), MuonIso(0.15));
    MuId = AndId<Muon>(MuonID(Muon::CutBasedIdTight), PtEtaCut(30., 2.4), MuonIso(0.15));

    //JetID = PtEtaCut(10., 2.4);
    // h_ST = ctx.get_handle<float>("ST");
    is_mc = ctx.get("dataset_type") == "MC";
    is_ele_channel = ctx.get("EleChannel") =="true";
    h_is_mlq_reconstructed = ctx.get_handle<bool>("is_mlq_reconstructed");
    h_mlq_reco_mode = ctx.get_handle<TString>("mlq_reco_mode");
    //-----------------------------------------------------------------------------------

    common.reset(new CommonModules());

    //common->switch_jetPtSorter();

    common->set_electron_id(ElectronID);

    common->set_muon_id(MuId);

    common->switch_metcorrection(true);

    common->init(ctx);

    //common->set_jet_id(JetID);

    //common->init(ctx);

    muoncleaner.reset(new MuonCleaner(MuId));

    jetcleaner.reset(new JetCleaner(ctx, 10.0, 2.5));

    //-----------------------------------------------------------
    // common->switch_metcorrection(true);
    //muoncleaner.reset(new MuonCleaner(MuId));



    if(is_ele_channel) muon_selection.reset(new NMuonSelection(0, 0));
    nele_selection.reset(new NElectronSelection(2, -1));
    // jet_selection.reset(new NJetSelection(2));
    lumi_selection.reset(new LumiSelection(ctx));
    m_ee_selection.reset(new InvMassEleEleSelection(71., 111.));
    // st_selection.reset(new StSelection(ctx, 350.));


    h_nocuts.reset(new LQTopMuRun2Hists(ctx, "NoCuts"));
    h_jets_nocuts.reset(new JetHists(ctx, "Jets_NoCuts"));
    h_ele_nocuts.reset(new ElectronHists(ctx, "Ele_NoCuts"));
    h_mu_nocuts.reset(new MuonHists(ctx, "Mu_NoCuts"));
    h_event_nocuts.reset(new EventHists(ctx, "Event_NoCuts"));
    h_topjets_nocuts.reset(new TopJetHists(ctx, "Topjets_NoCuts"));
    h_lumi_nocuts.reset(new LuminosityHists(ctx, "Lumi_NoCuts"));

    h_nmu.reset(new LQTopMuRun2Hists(ctx, "NMu"));
    h_jets_nmu.reset(new JetHists(ctx, "Jets_NMu"));
    h_ele_nmu.reset(new ElectronHists(ctx, "Ele_NMu"));
    h_mu_nmu.reset(new MuonHists(ctx, "Mu_NMu"));
    h_event_nmu.reset(new EventHists(ctx, "Event_NMu"));
    h_topjets_nmu.reset(new TopJetHists(ctx, "Topjets_NMu"));
    h_lumi_nmu.reset(new LuminosityHists(ctx, "Lumi_NMu"));

    h_cleaner.reset(new LQTopMuRun2Hists(ctx, "Cleaner"));
    h_jets_cleaner.reset(new JetHists(ctx, "Jets_Cleaner"));
    h_ele_cleaner.reset(new ElectronHists(ctx, "Ele_Cleaner"));
    h_mu_cleaner.reset(new MuonHists(ctx, "Mu_Cleaner"));
    h_event_cleaner.reset(new EventHists(ctx, "Event_Cleaner"));
    h_topjets_cleaner.reset(new TopJetHists(ctx, "Topjets_Cleaner"));
    h_lumi_cleaner.reset(new LuminosityHists(ctx, "Lumi_Cleaner"));


    // h_njets.reset(new LQTopMuRun2Hists(ctx, "NJets"));
    // h_jets_njets.reset(new JetHists(ctx, "Jets_NJets"));
    // h_ele_njets.reset(new ElectronHists(ctx, "Ele_NJets"));
    // h_mu_njets.reset(new MuonHists(ctx, "Mu_NJets"));
    // h_event_njets.reset(new EventHists(ctx, "Event_NJets"));
    // h_lumi_njets.reset(new LuminosityHists(ctx, "Lumi_NJets"));

    h_2ele.reset(new LQTopMuRun2Hists(ctx, "2Ele"));
    h_jets_2ele.reset(new JetHists(ctx, "Jets_2Ele"));
    h_ele_2ele.reset(new ElectronHists(ctx, "Ele_2Ele"));
    h_mu_2ele.reset(new MuonHists(ctx, "Mu_2Ele"));
    h_event_2ele.reset(new EventHists(ctx, "Event_2Ele"));
    h_topjets_2ele.reset(new TopJetHists(ctx, "Topjets_2Ele"));
    h_lumi_2ele.reset(new LuminosityHists(ctx, "Lumi_2Ele"));

    h_mee.reset(new LQTopMuRun2Hists(ctx, "Mee"));
    h_jets_mee.reset(new JetHists(ctx, "Jets_Mee"));
    h_ele_mee.reset(new ElectronHists(ctx, "Ele_Mee"));
    h_mu_mee.reset(new MuonHists(ctx, "Mu_Mee"));
    h_event_mee.reset(new EventHists(ctx, "Event_Mee"));
    h_topjets_mee.reset(new TopJetHists(ctx, "Topjets_Mee"));
    h_lumi_mee.reset(new LuminosityHists(ctx, "Lumi_Mee"));

    // h_stselec.reset(new LQTopMuRun2Hists(ctx, "stselec"));
    // h_jets_stselec.reset(new JetHists(ctx, "Jets_stselec"));
    // h_ele_stselec.reset(new ElectronHists(ctx, "Ele_stselec"));
    // h_mu_stselec.reset(new MuonHists(ctx, "Mu_stselec"));
    // h_event_stselec.reset(new EventHists(ctx, "Event_stselec"));
    // h_topjets_stselec.reset(new TopJetHists(ctx, "Topjets_stselec"));
    // h_lumi_stselec.reset(new LuminosityHists(ctx, "Lumi_stselec"));


  }

  bool LQTopMuRun2PreselectionLeptonMissidentification::process(Event & event) {

    event.set(h_is_mlq_reconstructed, false);
    event.set(h_mlq_reco_mode, "none");



    // float ST = 0.;
    // // event.muons->at(i) enstpricht event.muons[i] (an der iten Stelle)
    // for ( unsigned int i = 0; i < event.muons->size() ; i++){
    //   ST += event.muons->at(i).pt();
    // }
    //
    // for ( unsigned int i = 0; i < event.electrons->size() ; i++){
    //   ST += event.electrons->at(i).pt();
    // }
    //
    // for ( unsigned int i = 0; i < event.jets->size() ; i++){
    //   ST += event.jets->at(i).pt();
    // }
    //
    // ST += event.met->pt();
    //
    // // Befuellung des Handle!!!
    // event.set(h_ST, ST);


    if(!is_mc){
      if(!lumi_selection->passes(event)) return false;
    }




    h_nocuts->fill(event);
    h_jets_nocuts->fill(event);
    h_ele_nocuts->fill(event);
    h_mu_nocuts->fill(event);
    h_event_nocuts->fill(event);
    h_topjets_nocuts->fill(event);
    h_lumi_nocuts->fill(event);

    muoncleaner->process(event);
    if(is_ele_channel){
      if(!muon_selection->passes(event)) return false;
    }
    h_nmu->fill(event);
    h_jets_nmu->fill(event);
    h_ele_nmu->fill(event);
    h_mu_nmu->fill(event);
    h_event_nmu->fill(event);
    h_topjets_nmu->fill(event);
    h_lumi_nmu->fill(event);

    bool pass_common = common->process(event);
    if(!pass_common) return false;
    jetcleaner->process(event);


    h_cleaner->fill(event);
    h_jets_cleaner->fill(event);
    h_ele_cleaner->fill(event);
    h_mu_cleaner->fill(event);
    h_event_cleaner->fill(event);
    h_topjets_cleaner->fill(event);
    h_lumi_cleaner->fill(event);

    if(!nele_selection->passes(event)) return false;
    h_2ele->fill(event);
    h_jets_2ele->fill(event);
    h_ele_2ele->fill(event);
    h_mu_2ele->fill(event);
    h_event_2ele->fill(event);
    h_topjets_2ele->fill(event);
    h_lumi_2ele->fill(event);

    if(!m_ee_selection->passes(event)) return false;
    h_mee->fill(event);
    h_jets_mee->fill(event);
    h_ele_mee->fill(event);
    h_mu_mee->fill(event);
    h_event_mee->fill(event);
    h_topjets_mee->fill(event);
    h_lumi_mee->fill(event);



    return true;
  }

  UHH2_REGISTER_ANALYSIS_MODULE(LQTopMuRun2PreselectionLeptonMissidentification)

}
