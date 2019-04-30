#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/CommonModules.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/ElectronHists.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/LQTopMuRun2/include/LQTopMuRun2Selections.h"
#include "UHH2/LQTopMuRun2/include/LQTopMuRun2Hists.h"
#include "UHH2/common/include/JetIds.h"
#include "UHH2/common/include/ElectronIds.h"
#include "UHH2/common/include/MuonIds.h"
#include "UHH2/common/include/MuonHists.h"
#include "UHH2/common/include/JetHists.h"
#include "UHH2/common/include/EventHists.h"
#include "UHH2/common/include/LuminosityHists.h"
#include "UHH2/LQTopMuRun2/include/LQTopMuRun2Hists.h"

using namespace std;
using namespace uhh2;

namespace uhh2examples {

  class LQTopMuRun2PreselectionModule: public AnalysisModule {
  public:

    explicit LQTopMuRun2PreselectionModule(Context & ctx);
    virtual bool process(Event & event) override;
    void book_histograms(Context & ctx, vector<unique_ptr<Hists>> &hists, string tag);

  private:

    unique_ptr<CommonModules> common;

    unique_ptr<Hists> h_nocuts, h_jets_nocuts, h_ele_nocuts, h_mu_nocuts, h_event_nocuts, h_lumi_nocuts;
    unique_ptr<Hists> h_nmuons, h_jets_nmuons, h_ele_nmuons, h_mu_nmuons, h_event_nmuons, h_lumi_nmuons;
    unique_ptr<Hists> h_njets, h_jets_njets, h_ele_njets, h_mu_njets, h_event_njets, h_lumi_njets;

    unique_ptr<Hists> h_stselec, h_jets_stselec, h_ele_stselec, h_mu_stselec, h_event_stselec, h_lumi_stselec;

    unique_ptr<Selection> muon_selection, jet_selection, st_selection;

    MuonId MuonID;
    ElectronId ElectronID;
    JetId JetID;

    uhh2::Event::Handle<float> h_ST;
    uhh2::Event::Handle<bool> h_is_mlq_reconstructed;
    uhh2::Event::Handle<TString> h_mlq_reco_mode;
  };


  LQTopMuRun2PreselectionModule::LQTopMuRun2PreselectionModule(Context & ctx){
    cout << "Hello World from LQTopMuRun2PreselectionModule!" << endl;
    for(auto & kv : ctx.get_all()) cout << " " << kv.first << " = " << kv.second << endl;

    MuonID = AndId<Muon>(MuonIDTight(), PtEtaCut(30., 2.4), MuonIso(0.15));
    // PtEtaCut(float min_pt_, float max_eta_, float max_pt_ =-1, float min_eta_=-1)
    ElectronID = AndId<Electron>(ElectronID_Spring16_loose, PtEtaCut(30., 2.4));
    JetID = PtEtaCut(30., 2.4);
    h_ST = ctx.get_handle<float>("ST");
    h_is_mlq_reconstructed = ctx.get_handle<bool>("is_mlq_reconstructed");
    h_mlq_reco_mode = ctx.get_handle<TString>("mlq_reco_mode");


    // so setzt man unique pointer auf instanz.reset()
    common.reset(new CommonModules());
    common->switch_jetlepcleaner();
    common->switch_jetPtSorter();
    common->set_muon_id(MuonID);
    common->set_electron_id(ElectronID);
    common->set_jet_id(JetID);
    common->init(ctx);

    muon_selection.reset(new NMuonSelection(2));
    jet_selection.reset(new NJetSelection(2));
    st_selection.reset(new StSelection(ctx, 350.));


    h_nocuts.reset(new LQTopMuRun2Hists(ctx, "NoCuts"));
    h_jets_nocuts.reset(new JetHists(ctx, "Jets_NoCuts"));
    h_ele_nocuts.reset(new ElectronHists(ctx, "Ele_NoCuts"));
    h_mu_nocuts.reset(new MuonHists(ctx, "Mu_NoCuts"));
    h_event_nocuts.reset(new EventHists(ctx, "Event_NoCuts"));
    h_lumi_nocuts.reset(new LuminosityHists(ctx, "Lumi_NoCuts"));

    h_nmuons.reset(new LQTopMuRun2Hists(ctx, "NMuons"));
    h_jets_nmuons.reset(new JetHists(ctx, "Jets_NMuons"));
    h_ele_nmuons.reset(new ElectronHists(ctx, "Ele_NMuons"));
    h_mu_nmuons.reset(new MuonHists(ctx, "Mu_NMuons"));
    h_event_nmuons.reset(new EventHists(ctx, "Event_NMuons"));
    h_lumi_nmuons.reset(new LuminosityHists(ctx, "Lumi_NMuons"));

    h_njets.reset(new LQTopMuRun2Hists(ctx, "NJets"));
    h_jets_njets.reset(new JetHists(ctx, "Jets_NJets"));
    h_ele_njets.reset(new ElectronHists(ctx, "Ele_NJets"));
    h_mu_njets.reset(new MuonHists(ctx, "Mu_NJets"));
    h_event_njets.reset(new EventHists(ctx, "Event_NJets"));
    h_lumi_njets.reset(new LuminosityHists(ctx, "Lumi_NJets"));

    h_stselec.reset(new LQTopMuRun2Hists(ctx, "stselec"));
    h_jets_stselec.reset(new JetHists(ctx, "Jets_stselec"));
    h_ele_stselec.reset(new ElectronHists(ctx, "Ele_stselec"));
    h_mu_stselec.reset(new MuonHists(ctx, "Mu_stselec"));
    h_event_stselec.reset(new EventHists(ctx, "Event_stselec"));
    h_lumi_stselec.reset(new LuminosityHists(ctx, "Lumi_stselec"));

  }

  /*
  ██████  ██████   ██████   ██████ ███████ ███████ ███████
  ██   ██ ██   ██ ██    ██ ██      ██      ██      ██
  ██████  ██████  ██    ██ ██      █████   ███████ ███████
  ██      ██   ██ ██    ██ ██      ██           ██      ██
  ██      ██   ██  ██████   ██████ ███████ ███████ ███████
  */

  bool LQTopMuRun2PreselectionModule::process(Event & event) {
    event.set(h_is_mlq_reconstructed, false);
    event.set(h_mlq_reco_mode, "none");


    bool pass_common = common->process(event);
    if(!pass_common) return false;

    float ST = 0.;
    // event.muons->at(i) enstpricht event.muons[i] (an der iten Stelle)
    for ( unsigned int i = 0; i < event.muons->size() ; i++){
      ST += event.muons->at(i).pt();
    }

    for ( unsigned int i = 0; i < event.electrons->size() ; i++){
      ST += event.electrons->at(i).pt();
    }

    for ( unsigned int i = 0; i < event.jets->size() ; i++){
      ST += event.jets->at(i).pt();
    }

    ST += event.met->pt();

    // Befuellung des Handle!!!
    event.set(h_ST, ST);

    h_nocuts->fill(event);
    h_jets_nocuts->fill(event);
    h_ele_nocuts->fill(event);
    h_mu_nocuts->fill(event);
    h_event_nocuts->fill(event);
    h_lumi_nocuts->fill(event);

    if(!muon_selection->passes(event)) return false;
    h_nmuons->fill(event);
    h_jets_nmuons->fill(event);
    h_ele_nmuons->fill(event);
    h_mu_nmuons->fill(event);
    h_event_nmuons->fill(event);
    h_lumi_nmuons->fill(event);

    if(!jet_selection->passes(event)) return false;
    h_njets->fill(event);
    h_jets_njets->fill(event);
    h_ele_njets->fill(event);
    h_mu_njets->fill(event);
    h_event_njets->fill(event);
    h_lumi_njets->fill(event);

    if (!st_selection->passes(event)) return false;
    h_stselec->fill(event);
    h_jets_stselec->fill(event);
    h_ele_stselec->fill(event);
    h_mu_stselec->fill(event);
    h_event_stselec->fill(event);
    h_lumi_stselec->fill(event);




    return true;
  }


  UHH2_REGISTER_ANALYSIS_MODULE(LQTopMuRun2PreselectionModule)

}
