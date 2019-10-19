#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/CommonModules.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/ElectronHists.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/TriggerSelection.h"
#include "UHH2/LQTopMuRun2/include/LQTopMuRun2Selections.h"
#include "UHH2/LQTopMuRun2/include/LQTopMuRun2Hists.h"
#include "UHH2/common/include/JetIds.h"
#include "UHH2/common/include/ElectronIds.h"
#include "UHH2/common/include/MuonIds.h"
#include "UHH2/common/include/MuonHists.h"
#include "UHH2/common/include/JetHists.h"
#include "UHH2/common/include/EventHists.h"
#include "UHH2/common/include/LuminosityHists.h"
#include "UHH2/LQTopMuRun2/include/LQTopMuRun2EleTriggerEffHists.h"
#include "UHH2/common/include/MCWeight.h"
#include "UHH2/LQTopMuRun2/include/LQTopMuRun2EleTriggerWeights.h"
#include "UHH2/common/include/Utils.h"

using namespace std;
using namespace uhh2;

namespace uhh2examples {

  class LQTopMuRun2FullSelectionEleTriggerEff: public AnalysisModule {
  public:

    explicit LQTopMuRun2FullSelectionEleTriggerEff(Context & ctx);
    virtual bool process(Event & event) override;
    void book_histograms(Context & ctx, vector<unique_ptr<Hists>> &hists, string tag);


  private:

    bool is_e_e, is_mc, is_mu_e, apply_EleTriggerSF;

    unique_ptr<CommonModules> common;


    unique_ptr<Selection>  trigger1_mu_selection, trigger2_mu_selection, trigger1_ele_selection, trigger2_ele_selection, trigger3_ele_selection, nele_sel_trigger30, nele_sel_trigger120;

    unique_ptr<Hists> h_nocuts, h_jets_nocuts, h_ele_nocuts, h_mu_nocuts, h_event_nocuts, h_lumi_nocuts, h_eletriggereff_nocuts;
    unique_ptr<Hists> h_mutrigger, h_jets_mutrigger, h_ele_mutrigger, h_mu_mutrigger, h_event_mutrigger, h_lumi_mutrigger, h_eletriggereff_mutrigger;
    unique_ptr<Hists> h_eletrigger123, h_jets_eletrigger123, h_ele_eletrigger123, h_mu_eletrigger123, h_event_eletrigger123, h_lumi_eletrigger123, h_eletriggereff_eletrigger123;
    unique_ptr<Hists> h_trigger30_beforeletrigger, h_jets_trigger30_beforeletrigger, h_ele_trigger30_beforeletrigger, h_mu_trigger30_beforeletrigger, h_event_trigger30_beforeletrigger, h_lumi_trigger30_beforeletrigger, h_eletriggereff_trigger30_beforeletrigger;
    unique_ptr<Hists> h_trigger120_beforeletrigger, h_jets_trigger120_beforeletrigger, h_ele_trigger120_beforeletrigger, h_mu_trigger120_beforeletrigger, h_event_trigger120_beforeletrigger, h_lumi_trigger120_beforeletrigger, h_eletriggereff_trigger120_beforeletrigger;
    unique_ptr<Hists> h_trigger30_neg_120_beforeletrigger, h_jets_trigger30_neg_120_beforeletrigger, h_ele_trigger30_neg_120_beforeletrigger, h_mu_trigger30_neg_120_beforeletrigger, h_event_trigger30_neg_120_beforeletrigger, h_lumi_trigger30_neg_120_beforeletrigger, h_eletriggereff_trigger30_neg_120_beforeletrigger;
    unique_ptr<Hists> h_trigger30_aftereletrigger, h_jets_trigger30_aftereletrigger, h_ele_trigger30_aftereletrigger, h_mu_trigger30_aftereletrigger, h_event_trigger30_aftereletrigger, h_lumi_trigger30_aftereletrigger, h_eletriggereff_trigger30_aftereletrigger;
    unique_ptr<Hists> h_trigger120_aftereletrigger, h_jets_trigger120_aftereletrigger, h_ele_trigger120_aftereletrigger, h_mu_trigger120_aftereletrigger, h_event_trigger120_aftereletrigger, h_lumi_trigger120_aftereletrigger, h_eletriggereff_trigger120_aftereletrigger;
    unique_ptr<Hists> h_trigger30_neg_120_aftereletrigger, h_jets_trigger30_neg_120_aftereletrigger, h_ele_trigger30_neg_120_aftereletrigger, h_mu_trigger30_neg_120_aftereletrigger, h_event_trigger30_neg_120_aftereletrigger, h_lumi_trigger30_neg_120_aftereletrigger, h_eletriggereff_trigger30_neg_120_aftereletrigger;



    unique_ptr<Selection> trigger1_selection, trigger2_selection, trigger3_selection;

    unique_ptr<AnalysisModule> SF_muonID, SF_muonIso, SF_eleReco, SF_eleID;
    unique_ptr<ElectronTriggerWeights> SF_eleTrigger;
    unique_ptr<MCMuonScaleFactor> SF_muonTrigger;




    // now set the unique pointer for the ReconstructionHypothesis !!!
    //unique_ptr<xx> xx;


    JetId BTagID;
    CSVBTag::wp wp_btag_loose;
    BTag::algo btag_algo;



    MuonId MuID;
    ElectronId ElectronID;
    ElectronId ElectronID_trigger30;
    ElectronId ElectronID_trigger120;
    //ElectronId ElectronID_trig50;
    JetId JetID;

    uhh2::Event::Handle<bool> h_is_mlq_reconstructed;
    uhh2::Event::Handle<TString> h_mlq_reco_mode;

    // uhh2::Event::Handle<ElectronId> ElectronID;
    // uhh2::Event::Handle<MuonId> MuonID;
    // uhh2::Event::Handle<JetId> JetID;


  };


  LQTopMuRun2FullSelectionEleTriggerEff::LQTopMuRun2FullSelectionEleTriggerEff(Context & ctx){
    cout << "Hello World from LQTopMuRun2FullSelectionEleTriggerEff!" << endl;
    for(auto & kv : ctx.get_all()) cout << " " << kv.first << " = " << kv.second << endl;



    Year year = extract_year(ctx);
    if (year == Year::is2016v3) {
      ElectronID = AndId<Electron>(ElectronID_Summer16_tight, PtEtaCut(10.0, 2.4));
      ElectronID_trigger30 = AndId<Electron>(ElectronID_Summer16_tight, PtEtaCut(30.0, 2.4));
      ElectronID_trigger120 = AndId<Electron>(ElectronID_Summer16_tight, PtEtaCut(120.0, 2.4));
    }
    else if (year == Year::is2017v1) {
      ElectronID = AndId<Electron>(ElectronID_Fall17_tight, PtEtaCut(10.0, 2.4));
      ElectronID_trigger30 = AndId<Electron>(ElectronID_Fall17_tight, PtEtaCut(30.0, 2.4));
      ElectronID_trigger120 = AndId<Electron>(ElectronID_Fall17_tight, PtEtaCut(120.0, 2.4));
    }
    else if (year == Year::is2018) {
      ElectronID = AndId<Electron>(ElectronID_Fall17_tight, PtEtaCut(10.0, 2.4));
      ElectronID_trigger30 = AndId<Electron>(ElectronID_Fall17_tight, PtEtaCut(30.0, 2.4));
      ElectronID_trigger120 = AndId<Electron>(ElectronID_Fall17_tight, PtEtaCut(120.0, 2.4));
    }

    MuID = AndId<Muon>(MuonID(Muon::CutBasedIdTight), MuonIso(0.15));


    // ElectronID = AndId<Electron>(ElectronID_Spring16_tight, PtEtaCut(10., 2.4));
    // ElectronID_trigger30 = AndId<Electron>(ElectronID_Spring16_tight,PtEtaCut(30.0, 2.4)); //30
    //ElectronID_trigger50 = AndId<Electron>(ElectronID_Spring16_tight,PtEtaCut(50.0, 2.4));
    // ElectronID_trigger120 = AndId<Electron>(ElectronID_Spring16_tight,PtEtaCut(120.0, 2.4));
    // MuonID = AndId<Muon>(MuonIDTight(), PtEtaCut(30., 2.4), MuonIso(0.15));

    JetID = PtEtaCut(30., 2.4);



    common.reset(new CommonModules());
    //common->switch_jetlepcleaner();
    //common->switch_jetPtSorter();
    //common->set_muon_id(MuonID);
    //common->set_electron_id(ElectronID);
    //common->set_jet_id(JetID);
    common->disable_lumisel();
    common->disable_jersmear();
    common->disable_jec();
    common->set_electron_id(ElectronID);
    common->set_muon_id(MuID);
    common->init(ctx);


    h_is_mlq_reconstructed = ctx.get_handle<bool>("is_mlq_reconstructed");
    h_mlq_reco_mode = ctx.get_handle<TString>("mlq_reco_mode");

    is_mc = ctx.get("dataset_type") == "MC";
    is_mu_e = (ctx.get("channel") == "mu_e" || ctx.get("channel") == "e_mu");

    is_e_e = ctx.get("channel") == "e_e";
    if((!is_mu_e && !is_e_e)) throw runtime_error("In SidebandPreselectionModule: Invalid definition of 'channel' in config file, must be 'mu_e', 'e_mu', or 'e_e'.");
    apply_EleTriggerSF = (is_mc && is_e_e);

    SF_muonID.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/ehlersni/RUN2018/CMSSW_10_2_10/src/UHH2/common/data/MuonID_EfficienciesAndSF_average_RunBtoH.root", "MC_NUM_TightID_DEN_genTracks_PAR_pt_eta", 1., "tightID", true, "nominal"));
    SF_muonTrigger.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/ehlersni/RUN2018/CMSSW_10_2_10/src/UHH2/common/data/MuonTrigger_EfficienciesAndSF_average_RunBtoH.root", "IsoMu24_OR_IsoTkMu24_PtEtaBins", 0.5, "trigger", true, "nominal"));
    SF_muonIso.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/ehlersni/RUN2018/CMSSW_10_2_10/src/UHH2/common/data/MuonIso_EfficienciesAndSF_average_RunBtoH.root", "TightISO_TightID_pt_eta", 1., "iso", true, "nominal"));

    SF_eleReco.reset(new MCElecScaleFactor(ctx, "/nfs/dust/cms/user/ehlersni/RUN2018/CMSSW_10_2_10/src/UHH2/common/data/egammaEffi.txt_EGM2D_RecEff_Moriond17.root", 1, "", "nominal"));
    SF_eleID.reset(new MCElecScaleFactor(ctx, "/nfs/dust/cms/user/ehlersni/RUN2018/CMSSW_10_2_10/src/UHH2/common/data/egammaEffi.txt_EGM2D_CutBased_Loose_ID.root", 1, "", "nominal"));






    //Old Trigger Selections
    trigger1_mu_selection.reset(new TriggerSelection("HLT_IsoMu24_v*"));
    trigger2_mu_selection.reset(new TriggerSelection("HLT_IsoTkMu24_v*"));

    // New Trigger Selections
    trigger1_ele_selection.reset(new TriggerSelection("HLT_Ele27_WPTight_Gsf_v*"));
    trigger2_ele_selection.reset(new TriggerSelection("HLT_Ele115_CaloIdVT_GsfTrkIdT_v*"));
    trigger3_ele_selection.reset(new TriggerSelection("HLT_Photon175_v*"));
    nele_sel_trigger30.reset(new NElectronSelection(1, 1, ElectronID_trigger30));
    //nele_sel_trigger50.reset(new NElectronSelection(1, 1, ElectronID_trigger50));
    nele_sel_trigger120.reset(new NElectronSelection(1, 1, ElectronID_trigger120));


    // Hists
    h_nocuts.reset(new LQTopMuRun2Hists(ctx, "NoCuts"));
    h_jets_nocuts.reset(new JetHists(ctx, "Jets_NoCuts"));
    h_ele_nocuts.reset(new ElectronHists(ctx, "Ele_NoCuts"));
    h_mu_nocuts.reset(new MuonHists(ctx, "Mu_NoCuts"));
    h_event_nocuts.reset(new EventHists(ctx, "Event_NoCuts"));
    h_lumi_nocuts.reset(new LuminosityHists(ctx, "Lumi_NoCuts"));
    h_eletriggereff_nocuts.reset(new LQTopMuRun2EleTriggerEffHists(ctx, "EletriggerEff_NoCuts"));

    h_mutrigger.reset(new LQTopMuRun2Hists(ctx, "mutrigger"));
    h_jets_mutrigger.reset(new JetHists(ctx, "Jets_mutrigger"));
    h_ele_mutrigger.reset(new ElectronHists(ctx, "Ele_mutrigger"));
    h_mu_mutrigger.reset(new MuonHists(ctx, "Mu_mutrigger"));
    h_event_mutrigger.reset(new EventHists(ctx, "Event_mutrigger"));
    h_lumi_mutrigger.reset(new LuminosityHists(ctx, "Lumi_mutrigger"));
    h_eletriggereff_mutrigger.reset(new LQTopMuRun2EleTriggerEffHists(ctx, "EletriggerEff_mutrigger"));

    h_eletrigger123.reset(new LQTopMuRun2Hists(ctx, "eletrigger123"));
    h_jets_eletrigger123.reset(new JetHists(ctx, "Jets_eletrigger123"));
    h_ele_eletrigger123.reset(new ElectronHists(ctx, "Ele_eletrigger123"));
    h_mu_eletrigger123.reset(new MuonHists(ctx, "Mu_eletrigger123"));
    h_event_eletrigger123.reset(new EventHists(ctx, "Event_eletrigger123"));
    h_lumi_eletrigger123.reset(new LuminosityHists(ctx, "Lumi_eletrigger123"));
    h_eletriggereff_eletrigger123.reset(new LQTopMuRun2EleTriggerEffHists(ctx, "EletriggerEff_eletriger123"));


    // I dont know just combinations

    h_trigger30_beforeletrigger.reset(new LQTopMuRun2Hists(ctx, "trigger30_beforeletrigger"));
    h_jets_trigger30_beforeletrigger.reset(new JetHists(ctx, "Jets_trigger30_beforeletrigger"));
    h_ele_trigger30_beforeletrigger.reset(new ElectronHists(ctx, "Ele_trigger30_beforeletrigger"));
    h_mu_trigger30_beforeletrigger.reset(new MuonHists(ctx, "Mu_trigger30_beforeletrigger"));
    h_event_trigger30_beforeletrigger.reset(new EventHists(ctx, "Event_trigger30_beforeletrigger"));
    h_lumi_trigger30_beforeletrigger.reset(new LuminosityHists(ctx, "Lumi_trigger30_beforeletrigger"));
    h_eletriggereff_trigger30_beforeletrigger.reset(new LQTopMuRun2EleTriggerEffHists(ctx, "EletriggerEff_trigger30_beforeletrigger"));

    h_trigger120_beforeletrigger.reset(new LQTopMuRun2Hists(ctx, "trigger120_beforeletrigger"));
    h_jets_trigger120_beforeletrigger.reset(new JetHists(ctx, "Jets_trigger120_beforeletrigger"));
    h_ele_trigger120_beforeletrigger.reset(new ElectronHists(ctx, "Ele_trigger120_beforeletrigger"));
    h_mu_trigger120_beforeletrigger.reset(new MuonHists(ctx, "Mu_trigger120_beforeletrigger"));
    h_event_trigger120_beforeletrigger.reset(new EventHists(ctx, "Event_trigger120_beforeletrigger"));
    h_lumi_trigger120_beforeletrigger.reset(new LuminosityHists(ctx, "Lumi_trigger120_beforeletrigger"));
    h_eletriggereff_trigger120_beforeletrigger.reset(new LQTopMuRun2EleTriggerEffHists(ctx, "EletriggerEff_trigger120_beforeletrigger"));

    h_trigger30_neg_120_beforeletrigger.reset(new LQTopMuRun2Hists(ctx, "trigger30_neg_120_beforeletrigger"));
    h_jets_trigger30_neg_120_beforeletrigger.reset(new JetHists(ctx, "Jets_trigger30_neg_120_beforeletrigger"));
    h_ele_trigger30_neg_120_beforeletrigger.reset(new ElectronHists(ctx, "Ele_trigger30_neg_120_beforeletrigger"));
    h_mu_trigger30_neg_120_beforeletrigger.reset(new MuonHists(ctx, "Mu_trigger30_neg_120_beforeletrigger"));
    h_event_trigger30_neg_120_beforeletrigger.reset(new EventHists(ctx, "Event_trigger30_neg_120_beforeletrigger"));
    h_lumi_trigger30_neg_120_beforeletrigger.reset(new LuminosityHists(ctx, "Lumi_trigger30_neg_120_beforeletrigger"));
    h_eletriggereff_trigger30_neg_120_beforeletrigger.reset(new LQTopMuRun2EleTriggerEffHists(ctx, "EletriggerEff_trigger30_neg_120_beforeletrigger"));



    h_trigger30_aftereletrigger.reset(new LQTopMuRun2Hists(ctx, "trigger30_aftereletrigger"));
    h_jets_trigger30_aftereletrigger.reset(new JetHists(ctx, "Jets_trigger30_aftereletrigger"));
    h_ele_trigger30_aftereletrigger.reset(new ElectronHists(ctx, "Ele_trigger30_aftereletrigger"));
    h_mu_trigger30_aftereletrigger.reset(new MuonHists(ctx, "Mu_trigger30_aftereletrigger"));
    h_event_trigger30_aftereletrigger.reset(new EventHists(ctx, "Event_trigger30_aftereletrigger"));
    h_lumi_trigger30_aftereletrigger.reset(new LuminosityHists(ctx, "Lumi_trigger30_aftereletrigger"));
    h_eletriggereff_trigger30_aftereletrigger.reset(new LQTopMuRun2EleTriggerEffHists(ctx, "EletriggerEff_trigger30_aftereletrigger"));

    h_trigger120_aftereletrigger.reset(new LQTopMuRun2Hists(ctx, "trigger120_aftereletrigger"));
    h_jets_trigger120_aftereletrigger.reset(new JetHists(ctx, "Jets_trigger120_aftereletrigger"));
    h_ele_trigger120_aftereletrigger.reset(new ElectronHists(ctx, "Ele_trigger120_aftereletrigger"));
    h_mu_trigger120_aftereletrigger.reset(new MuonHists(ctx, "Mu_trigger120_aftereletrigger"));
    h_event_trigger120_aftereletrigger.reset(new EventHists(ctx, "Event_trigger120_aftereletrigger"));
    h_lumi_trigger120_aftereletrigger.reset(new LuminosityHists(ctx, "Lumi_trigger120_aftereletrigger"));
    h_eletriggereff_trigger120_aftereletrigger.reset(new LQTopMuRun2EleTriggerEffHists(ctx, "EletriggerEff_trigger120_aftereletrigger"));

    h_trigger30_neg_120_aftereletrigger.reset(new LQTopMuRun2Hists(ctx, "trigger30_neg_120_aftereletrigger"));
    h_jets_trigger30_neg_120_aftereletrigger.reset(new JetHists(ctx, "Jets_trigger30_neg_120_aftereletrigger"));
    h_ele_trigger30_neg_120_aftereletrigger.reset(new ElectronHists(ctx, "Ele_trigger30_neg_120_aftereletrigger"));
    h_mu_trigger30_neg_120_aftereletrigger.reset(new MuonHists(ctx, "Mu_trigger30_neg_120_aftereletrigger"));
    h_event_trigger30_neg_120_aftereletrigger.reset(new EventHists(ctx, "Event_trigger30_neg_120_aftereletrigger"));
    h_lumi_trigger30_neg_120_aftereletrigger.reset(new LuminosityHists(ctx, "Lumi_trigger30_neg_120_aftereletrigger"));
    h_eletriggereff_trigger30_neg_120_aftereletrigger.reset(new LQTopMuRun2EleTriggerEffHists(ctx, "EletriggerEff_trigger30_neg_120_aftereletrigger"));




  }

  /*
  ██████  ██████   ██████   ██████ ███████ ███████ ███████
  ██   ██ ██   ██ ██    ██ ██      ██      ██      ██
  ██████  ██████  ██    ██ ██      █████   ███████ ███████
  ██      ██   ██ ██    ██ ██      ██           ██      ██
  ██      ██   ██  ██████   ██████ ███████ ███████ ███████
  */

  bool LQTopMuRun2FullSelectionEleTriggerEff::process(Event & event) {

    event.set(h_is_mlq_reconstructed, false);
    event.set(h_mlq_reco_mode, "none");



    bool pass_common = common->process(event);
    if(!pass_common) return false;

    SF_muonTrigger->process_onemuon(event,0);
    SF_muonID->process(event);
    SF_muonIso->process(event);


    SF_eleReco->process(event);
    SF_eleID->process(event);

    if(apply_EleTriggerSF){
      SF_eleTrigger->process(event);
    }

    //MuonTrigger_EfficienciesAndSF_average_RunBtoH

    h_nocuts->fill(event);
    h_jets_nocuts->fill(event);
    h_ele_nocuts->fill(event);
    h_mu_nocuts->fill(event);
    h_event_nocuts->fill(event);
    h_lumi_nocuts->fill(event);
    h_eletriggereff_nocuts->fill(event);


    if(!(trigger1_mu_selection->passes(event) || trigger2_mu_selection->passes(event))) return false;
    h_mutrigger->fill(event);
    h_jets_mutrigger->fill(event);
    h_ele_mutrigger->fill(event);
    h_mu_mutrigger->fill(event);
    h_event_mutrigger->fill(event);
    h_lumi_mutrigger->fill(event);
    h_eletriggereff_mutrigger->fill(event);


    if(nele_sel_trigger30->passes(event))
    {
      h_trigger30_beforeletrigger->fill(event);
      h_jets_trigger30_beforeletrigger->fill(event);
      h_ele_trigger30_beforeletrigger->fill(event);
      h_mu_trigger30_beforeletrigger->fill(event);
      h_event_trigger30_beforeletrigger->fill(event);
      h_lumi_trigger30_beforeletrigger->fill(event);
      h_eletriggereff_trigger30_beforeletrigger->fill(event);
    }

    if(nele_sel_trigger120->passes(event))
    {
      h_trigger120_beforeletrigger->fill(event);
      h_jets_trigger120_beforeletrigger->fill(event);
      h_ele_trigger120_beforeletrigger->fill(event);
      h_mu_trigger120_beforeletrigger->fill(event);
      h_event_trigger120_beforeletrigger->fill(event);
      h_lumi_trigger120_beforeletrigger->fill(event);
      h_eletriggereff_trigger120_beforeletrigger->fill(event);
    }

    if(nele_sel_trigger30->passes(event) && !(nele_sel_trigger120->passes(event)))
    {
      h_trigger30_neg_120_beforeletrigger->fill(event);
      h_jets_trigger30_neg_120_beforeletrigger->fill(event);
      h_ele_trigger30_neg_120_beforeletrigger->fill(event);
      h_mu_trigger30_neg_120_beforeletrigger->fill(event);
      h_event_trigger30_neg_120_beforeletrigger->fill(event);
      h_lumi_trigger30_neg_120_beforeletrigger->fill(event);
      h_eletriggereff_trigger30_neg_120_beforeletrigger->fill(event);
    }






    if(!(trigger1_ele_selection->passes(event) || trigger2_ele_selection->passes(event) || trigger3_ele_selection->passes(event))) return false;
    h_eletrigger123->fill(event);
    h_jets_eletrigger123->fill(event);
    h_ele_eletrigger123->fill(event);
    h_mu_eletrigger123->fill(event);
    h_event_eletrigger123->fill(event);
    h_lumi_eletrigger123->fill(event);
    h_eletriggereff_eletrigger123->fill(event);

    if(nele_sel_trigger30->passes(event))
    {
      h_trigger30_aftereletrigger->fill(event);
      h_jets_trigger30_aftereletrigger->fill(event);
      h_ele_trigger30_aftereletrigger->fill(event);
      h_mu_trigger30_aftereletrigger->fill(event);
      h_event_trigger30_aftereletrigger->fill(event);
      h_lumi_trigger30_aftereletrigger->fill(event);
      h_eletriggereff_trigger30_aftereletrigger->fill(event);
    }

    if(nele_sel_trigger120->passes(event))
    {
      h_trigger120_aftereletrigger->fill(event);
      h_jets_trigger120_aftereletrigger->fill(event);
      h_ele_trigger120_aftereletrigger->fill(event);
      h_mu_trigger120_aftereletrigger->fill(event);
      h_event_trigger120_aftereletrigger->fill(event);
      h_lumi_trigger120_aftereletrigger->fill(event);
      h_eletriggereff_trigger120_aftereletrigger->fill(event);
    }

    if(nele_sel_trigger30->passes(event) && !(nele_sel_trigger120->passes(event)))
    {
      h_trigger30_neg_120_aftereletrigger->fill(event);
      h_jets_trigger30_neg_120_aftereletrigger->fill(event);
      h_ele_trigger30_neg_120_aftereletrigger->fill(event);
      h_mu_trigger30_neg_120_aftereletrigger->fill(event);
      h_event_trigger30_neg_120_aftereletrigger->fill(event);
      h_lumi_trigger30_neg_120_aftereletrigger->fill(event);
      h_eletriggereff_trigger30_neg_120_aftereletrigger->fill(event);
    }



    return false;
  }


  UHH2_REGISTER_ANALYSIS_MODULE(LQTopMuRun2FullSelectionEleTriggerEff)

}
