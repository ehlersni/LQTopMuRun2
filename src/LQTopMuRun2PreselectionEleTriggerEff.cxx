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

  class LQTopMuRun2PreselectionEleTriggerEff: public AnalysisModule {
  public:

    explicit LQTopMuRun2PreselectionEleTriggerEff(Context & ctx);
    virtual bool process(Event & event) override;
    void book_histograms(Context & ctx, vector<unique_ptr<Hists>> &hists, string tag);

  private:

    unique_ptr<CommonModules> common;

    unique_ptr<Hists> h_nocuts;
    unique_ptr<Hists> h_nmuons;
    unique_ptr<Hists> h_njets;
    unique_ptr<Hists> h_neles;



    unique_ptr<Selection> muon_selection, jet_selection, ele_selection;

    MuonId MuonID;
    ElectronId ElectronID;
    JetId JetID;

    uhh2::Event::Handle<bool> h_is_mlq_reconstructed;
    uhh2::Event::Handle<TString> h_mlq_reco_mode;


  };

  LQTopMuRun2PreselectionEleTriggerEff::LQTopMuRun2PreselectionEleTriggerEff(Context & ctx)
  {
    cout << "Hello World from LQTopMuRun2PreselectionEleTriggerEff!" << endl;
    for(auto & kv : ctx.get_all()) cout << " " << kv.first << " = " << kv.second << endl;

    ElectronID = AndId<Electron>(ElectronID_Spring16_tight, PtEtaCut(10., 2.4));
    MuonID = AndId<Muon>(MuonIDTight(), PtEtaCut(30., 2.4), MuonIso(0.15));
    JetID = PtEtaCut(30., 2.4);



    common.reset(new CommonModules());
    common->switch_jetlepcleaner();
    common->switch_jetPtSorter();
    common->set_muon_id(MuonID);
    common->set_electron_id(ElectronID);
    common->set_jet_id(JetID);
    common->init(ctx);

    ele_selection.reset(new NElectronSelection(1, 1));
    muon_selection.reset(new NMuonSelection(1, 1));
    jet_selection.reset(new NJetSelection(2));


    h_is_mlq_reconstructed = ctx.get_handle<bool>("is_mlq_reconstructed");
    h_mlq_reco_mode = ctx.get_handle<TString>("mlq_reco_mode");

    //Histograms......









    // so setzt man unique pointer auf instanz.reset()



    h_nocuts.reset(new LQTopMuRun2Hists(ctx, "NoCuts"));


    h_nmuons.reset(new LQTopMuRun2Hists(ctx, "NMuons"));


    h_njets.reset(new LQTopMuRun2Hists(ctx, "NJets"));


    h_neles.reset(new LQTopMuRun2Hists(ctx, "NEles"));


  }

  /*
  ██████  ██████   ██████   ██████ ███████ ███████ ███████
  ██   ██ ██   ██ ██    ██ ██      ██      ██      ██
  ██████  ██████  ██    ██ ██      █████   ███████ ███████
  ██      ██   ██ ██    ██ ██      ██           ██      ██
  ██      ██   ██  ██████   ██████ ███████ ███████ ███████
  */

  bool LQTopMuRun2PreselectionEleTriggerEff::process(Event & event) {
    event.set(h_is_mlq_reconstructed, false);
    event.set(h_mlq_reco_mode, "none");

    bool pass_common = common->process(event);
    if(!pass_common) return false;


    // Befuellung des Handle!!!


    h_nocuts->fill(event);


    if(!muon_selection->passes(event)) return false;
    h_nmuons->fill(event);


    if(!jet_selection->passes(event)) return false;
    h_njets->fill(event);

    if(!ele_selection->passes(event)) return false;
    h_neles->fill(event);



    return true;
  }


  UHH2_REGISTER_ANALYSIS_MODULE(LQTopMuRun2PreselectionEleTriggerEff)

}
