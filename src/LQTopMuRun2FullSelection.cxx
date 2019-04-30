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
#include "UHH2/LQTopMuRun2/include/LQTopMuRun2Hists.h"
#include "UHH2/common/include/MCWeight.h"
#include "UHH2/LQTopMuRun2/include/LQReconstruction.h"
#include "UHH2/LQTopMuRun2/include/LQReconstructionHypothesisDiscriminators.h"

using namespace std;
using namespace uhh2;

namespace uhh2examples {

  class LQTopMuRun2FullSelection: public AnalysisModule {
  public:

    explicit LQTopMuRun2FullSelection(Context & ctx);
    virtual bool process(Event & event) override;
    void book_histograms(Context & ctx, vector<unique_ptr<Hists>> &hists, string tag);

  private:

    unique_ptr<CommonModules> common;

    unique_ptr<Hists> h_nocuts, h_jets_nocuts, h_ele_nocuts, h_mu_nocuts, h_event_nocuts, h_lumi_nocuts;
    unique_ptr<Hists> h_trigger, h_jets_trigger, h_ele_trigger, h_mu_trigger, h_event_trigger, h_lumi_trigger;
    unique_ptr<Hists> h_nmuons, h_jets_nmuons, h_ele_nmuons, h_mu_nmuons, h_event_nmuons, h_lumi_nmuons;
    unique_ptr<Hists> h_njets, h_jets_njets, h_ele_njets, h_mu_njets, h_event_njets, h_lumi_njets;


    unique_ptr<Hists> h_stselec, h_jets_stselec, h_ele_stselec, h_mu_stselec, h_event_stselec, h_lumi_stselec;
    unique_ptr<Hists> h_muinvariantmass, h_jets_muinvariantmass, h_ele_muinvariantmass, h_mu_muinvariantmass, h_event_muinvariantmass, h_lumi_muinvariantmass;
    unique_ptr<Hists> h_stlepselec, h_jets_stlepselec, h_ele_stlepselec, h_mu_stlepselec, h_event_stlepselec, h_lumi_stlepselec;
    unique_ptr<Hists> h_btageff;
    unique_ptr<Hists> h_nbtagsel, h_jets_nbtagsel, h_ele_nbtagsel, h_mu_nbtagsel, h_event_nbtagsel, h_lumi_nbtagsel;
    unique_ptr<Hists> h_finalselection, h_jets_finalselection, h_ele_finalselection, h_mu_finalselection, h_event_finalselection, h_lumi_finalselection;
    unique_ptr<Hists> h_finalselection_mlqtrue, h_jets_finalselection_mlqtrue, h_ele_finalselection_mlqtrue, h_mu_finalselection_mlqtrue, h_event_finalselection_mlqtrue, h_lumi_finalselection_mlqtrue;
    unique_ptr<Hists> h_finalselection_mlqfalse, h_jets_finalselection_mlqfalse, h_ele_finalselection_mlqfalse, h_mu_finalselection_mlqfalse, h_event_finalselection_mlqfalse, h_lumi_finalselection_mlqfalse;


    unique_ptr<Selection> trigger1_selection, trigger2_selection, muon_selection, jet_selection, st_selection, mu_invariantmass, stlep_selection, nbtag_selection;

    unique_ptr<AnalysisModule> SF_muonID, SF_muonIso, SF_eleReco, SF_eleID;
    unique_ptr<AnalysisModule> SF_btag;
    unique_ptr<MCMuonScaleFactor> SF_muonTrigger;

    unique_ptr<HighMassInclusiveLQReconstruction> mlq_reco;
    unique_ptr<LQChi2Discriminator> chi2_module;


    // now set the unique pointer for the ReconstructionHypothesis !!!
    //unique_ptr<xx> xx;


    JetId BTagID;
    CSVBTag::wp wp_btag_loose;

    uhh2::Event::Handle<float> h_ST;
    uhh2::Event::Handle<float> h_STlep;
    uhh2::Event::Handle<bool> h_is_mlq_reconstructed;
    uhh2::Event::Handle<TString> h_mlq_reco_mode;


  };


  LQTopMuRun2FullSelection::LQTopMuRun2FullSelection(Context & ctx){
    cout << "Hello World from LQTopMuRun2FullSelection!" << endl;
    for(auto & kv : ctx.get_all()) cout << " " << kv.first << " = " << kv.second << endl;

    h_ST = ctx.get_handle<float>("ST");
    h_STlep = ctx.get_handle<float>("STlep");
    h_is_mlq_reconstructed = ctx.get_handle<bool>("is_mlq_reconstructed");
    h_mlq_reco_mode = ctx.get_handle<TString>("mlq_reco_mode");


    wp_btag_loose = CSVBTag::WP_LOOSE;
    BTagID = CSVBTag(wp_btag_loose);

    // this is how to set unique pointer --> instanz.reset()
    common.reset(new CommonModules());
    // function which disables the jetsmearer
    // this how to set disable functions
    common->disable_jersmear();
    common->disable_jec();
    common->disable_lumisel();
    common->disable_metfilters();
    common->disable_pvfilter();
    common->disable_jetpfidfilter();

    common->init(ctx);

    mlq_reco.reset(new HighMassInclusiveLQReconstruction(ctx, LQNeutrinoReconstruction));
    chi2_module.reset(new LQChi2Discriminator(ctx, "LQHypotheses"));

    SF_muonID.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/ehlersni/CMSSW_8_0_24_patch1/src/UHH2/common/data/MuonID_EfficienciesAndSF_average_RunBtoH.root", "MC_NUM_TightID_DEN_genTracks_PAR_pt_eta", 1., "tightID", true, "nominal"));
    SF_muonTrigger.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/ehlersni/CMSSW_8_0_24_patch1/src/UHH2/common/data/MuonTrigger_EfficienciesAndSF_average_RunBtoH.root", "IsoMu24_OR_IsoTkMu24_PtEtaBins", 0.5, "trigger", true, "nominal"));
    SF_muonIso.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/ehlersni/CMSSW_8_0_24_patch1/src/UHH2/common/data/MuonIso_EfficienciesAndSF_average_RunBtoH.root", "TightISO_TightID_pt_eta", 1., "iso", true, "nominal"));

    SF_eleReco.reset(new MCElecScaleFactor(ctx, "/nfs/dust/cms/user/ehlersni/CMSSW_8_0_24_patch1/src/UHH2/common/data/egammaEffi.txt_EGM2D_RecEff_Moriond17.root", 1, "", "nominal"));
    SF_eleID.reset(new MCElecScaleFactor(ctx, "/nfs/dust/cms/user/ehlersni/CMSSW_8_0_24_patch1/src/UHH2/common/data/egammaEffi.txt_EGM2D_CutBased_Loose_ID.root", 1, "", "nominal"));

    SF_btag.reset(new MCBTagScaleFactor(ctx, wp_btag_loose, "jets", "nominal"));

    trigger1_selection.reset(new TriggerSelection("HLT_IsoMu24_v*"));
    trigger2_selection.reset(new TriggerSelection("HLT_IsoTkMu24_v*"));
    muon_selection.reset(new NMuonSelection(2));
    jet_selection.reset(new NJetSelection(2));
    st_selection.reset(new StSelection(ctx, 350.));
    mu_invariantmass.reset(new InvMassOf2Mu(111.));
    stlep_selection.reset(new StLeptonSelection(ctx, 200.));
    nbtag_selection.reset(new NJetSelection(1, -1, BTagID));


    h_nocuts.reset(new LQTopMuRun2Hists(ctx, "NoCuts"));
    h_jets_nocuts.reset(new JetHists(ctx, "Jets_NoCuts"));
    h_ele_nocuts.reset(new ElectronHists(ctx, "Ele_NoCuts"));
    h_mu_nocuts.reset(new MuonHists(ctx, "Mu_NoCuts"));
    h_event_nocuts.reset(new EventHists(ctx, "Event_NoCuts"));
    h_lumi_nocuts.reset(new LuminosityHists(ctx, "Lumi_NoCuts"));

    h_trigger.reset(new LQTopMuRun2Hists(ctx, "Trigger"));
    h_jets_trigger.reset(new JetHists(ctx, "Jets_Trigger"));
    h_ele_trigger.reset(new ElectronHists(ctx, "Ele_Trigger"));
    h_mu_trigger.reset(new MuonHists(ctx, "Mu_Trigger"));
    h_event_trigger.reset(new EventHists(ctx, "Event_Trigger"));
    h_lumi_trigger.reset(new LuminosityHists(ctx, "Lumi_Trigger"));

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

    // b-tag scale factors
    h_btageff.reset(new BTagMCEfficiencyHists(ctx, "BTagEff", wp_btag_loose));


    h_muinvariantmass.reset(new LQTopMuRun2Hists(ctx, "muinvariantmass"));
    h_jets_muinvariantmass.reset(new JetHists(ctx, "Jets_muinvariantmass"));
    h_ele_muinvariantmass.reset(new ElectronHists(ctx, "Ele_muinvariantmass"));
    h_mu_muinvariantmass.reset(new MuonHists(ctx, "Mu_muinvariantmass"));
    h_event_muinvariantmass.reset(new EventHists(ctx, "Event_muinvariantmass"));
    h_lumi_muinvariantmass.reset(new LuminosityHists(ctx, "Lumi_muinvariantmass"));

    h_stlepselec.reset(new LQTopMuRun2Hists(ctx, "stlepselec"));
    h_jets_stlepselec.reset(new JetHists(ctx, "Jets_stlepselec"));
    h_ele_stlepselec.reset(new ElectronHists(ctx, "Ele_stlepselec"));
    h_mu_stlepselec.reset(new MuonHists(ctx, "Mu_stlepselec"));
    h_event_stlepselec.reset(new EventHists(ctx, "Event_stlepselec"));
    h_lumi_stlepselec.reset(new LuminosityHists(ctx, "Lumi_stlepselec"));

    // for B-Tag


    h_nbtagsel.reset(new LQTopMuRun2Hists(ctx, "nbtagsel"));
    h_jets_nbtagsel.reset(new JetHists(ctx, "Jets_nbtagsel"));
    h_ele_nbtagsel.reset(new ElectronHists(ctx, "Ele_nbtagsel"));
    h_mu_nbtagsel.reset(new MuonHists(ctx, "Mu_nbtagsel"));
    h_event_nbtagsel.reset(new EventHists(ctx, "Event_nbtagsel"));
    h_lumi_nbtagsel.reset(new LuminosityHists(ctx, "Lumi_nbtagsel"));

    // for final StSelection

    h_finalselection.reset(new LQTopMuRun2Hists(ctx, "finalselection"));
    h_jets_finalselection.reset(new JetHists(ctx, "Jets_finalselection"));
    h_ele_finalselection.reset(new ElectronHists(ctx, "Ele_finalselection"));
    h_mu_finalselection.reset(new MuonHists(ctx, "Mu_finalselection"));
    h_event_finalselection.reset(new EventHists(ctx, "Event_finalselection"));
    h_lumi_finalselection.reset(new LuminosityHists(ctx, "Lumi_finalselection"));


    h_finalselection_mlqtrue.reset(new LQTopMuRun2Hists(ctx, "finalselection_mlqtrue"));
    h_jets_finalselection_mlqtrue.reset(new JetHists(ctx, "Jets_finalselection_mlqtrue"));
    h_ele_finalselection_mlqtrue.reset(new ElectronHists(ctx, "Ele_finalselection_mlqtrue"));
    h_mu_finalselection_mlqtrue.reset(new MuonHists(ctx, "Mu_finalselection_mlqtrue"));
    h_event_finalselection_mlqtrue.reset(new EventHists(ctx, "Event_finalselection_mlqtrue"));
    h_lumi_finalselection_mlqtrue.reset(new LuminosityHists(ctx, "Lumi_finalselection_mlqtrue"));

    h_finalselection_mlqfalse.reset(new LQTopMuRun2Hists(ctx, "finalselection_mlqfalse"));
    h_jets_finalselection_mlqfalse.reset(new JetHists(ctx, "Jets_finalselection_mlqfalse"));
    h_ele_finalselection_mlqfalse.reset(new ElectronHists(ctx, "Ele_finalselection_mlqfalse"));
    h_mu_finalselection_mlqfalse.reset(new MuonHists(ctx, "Mu_finalselection_mlqfalse"));
    h_event_finalselection_mlqfalse.reset(new EventHists(ctx, "Event_finalselection_mlqfalse"));
    h_lumi_finalselection_mlqfalse.reset(new LuminosityHists(ctx, "Lumi_finalselection_mlqfalse"));


  }

  /*
  ██████  ██████   ██████   ██████ ███████ ███████ ███████
  ██   ██ ██   ██ ██    ██ ██      ██      ██      ██
  ██████  ██████  ██    ██ ██      █████   ███████ ███████
  ██      ██   ██ ██    ██ ██      ██           ██      ██
  ██      ██   ██  ██████   ██████ ███████ ███████ ███████
  */

  bool LQTopMuRun2FullSelection::process(Event & event) {

    event.set(h_is_mlq_reconstructed, false);
    event.set(h_mlq_reco_mode, "none");



    bool pass_common = common->process(event);
    if(!pass_common) return false;

    SF_muonTrigger->process_onemuon(event,0);
    SF_muonID->process(event);
    SF_muonIso->process(event);
    SF_btag->process(event);


    SF_eleReco->process(event);
    SF_eleID->process(event);


    float ST = 0., STlep = 0.;
    // event.muons->at(i) is eaqual event.muons[i] (the i-th point)
    for ( unsigned int i = 0; i < event.muons->size() ; i++){
      ST += event.muons->at(i).pt();
      STlep += event.muons->at(i).pt();
    }

    for ( unsigned int i = 0; i < event.electrons->size() ; i++){
      ST += event.electrons->at(i).pt();
      STlep += event.electrons->at(i).pt();
    }

    for ( unsigned int i = 0; i < event.jets->size() ; i++){
      ST += event.jets->at(i).pt();
    }

    ST += event.met->pt();

    // fill of the handle!!!
    event.set(h_ST, ST);
    event.set(h_STlep, STlep);

    h_nocuts->fill(event);
    h_jets_nocuts->fill(event);
    h_ele_nocuts->fill(event);
    h_mu_nocuts->fill(event);
    h_event_nocuts->fill(event);
    h_lumi_nocuts->fill(event);

    if(!(trigger1_selection->passes(event) || trigger2_selection->passes(event))) return false;
    h_trigger->fill(event);
    h_jets_trigger->fill(event);
    h_ele_trigger->fill(event);
    h_mu_trigger->fill(event);
    h_event_trigger->fill(event);
    h_lumi_trigger->fill(event);

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

    h_btageff->fill(event);

    if (!mu_invariantmass->passes(event)) return false;
    h_muinvariantmass->fill(event);
    h_jets_muinvariantmass->fill(event);
    h_ele_muinvariantmass->fill(event);
    h_mu_muinvariantmass->fill(event);
    h_event_muinvariantmass->fill(event);
    h_lumi_muinvariantmass->fill(event);

    if (!stlep_selection->passes(event)) return false;
    h_stlepselec->fill(event);
    h_jets_stlepselec->fill(event);
    h_ele_stlepselec->fill(event);
    h_mu_stlepselec->fill(event);
    h_event_stlepselec->fill(event);
    h_lumi_stlepselec->fill(event);

    // for B-Tag

    if (!nbtag_selection->passes(event)) return false;
    h_nbtagsel->fill(event);
    h_jets_nbtagsel->fill(event);
    h_ele_nbtagsel->fill(event);
    h_mu_nbtagsel->fill(event);
    h_event_nbtagsel->fill(event);
    h_lumi_nbtagsel->fill(event);



    //check for at least 1 muon pair with opposite charge
    bool charge_opposite = false;
    for(unsigned int i=0; i<event.muons->size(); i++){
      for(unsigned int j=0; j<event.muons->size(); j++){
        if(j>i){
          if(event.muons->at(i).charge() != event.muons->at(j).charge()) {
            charge_opposite = true;
          }
        }
      }
    }

    if(event.muons->size() >= 2 && event.electrons->size() >= 1 && charge_opposite) event.set(h_mlq_reco_mode, "ele");
    else if(event.electrons->size() == 0 && event.muons->size() >= 3 && charge_opposite) event.set(h_mlq_reco_mode, "muon");
    if(event.get(h_mlq_reco_mode) != "none" && event.get(h_mlq_reco_mode) != "ele" && event.get(h_mlq_reco_mode) != "muon") throw runtime_error("'h_mlq_reco_mode' contains an invalid value!");

    // the funtion LQNeutrinoReconstruction has to reconstruct the Neutrinos for all events (if the event is leptonicaly) that fullfils the LQTopMuRun2FullSelection
    // therefore for the muon decay and the electron decay...

    if(event.get(h_mlq_reco_mode) == "ele" || event.get(h_mlq_reco_mode) == "muon")
    {
      mlq_reco->process(event);
      chi2_module->process(event);
    }


    bool is_mlq_reconstructed = event.get(h_is_mlq_reconstructed);

    h_finalselection->fill(event);
    h_jets_finalselection->fill(event);
    h_ele_finalselection->fill(event);
    h_mu_finalselection->fill(event);
    h_event_finalselection->fill(event);
    h_lumi_finalselection->fill(event);

    if(is_mlq_reconstructed == true)
    {
      h_finalselection_mlqtrue->fill(event);
      h_jets_finalselection_mlqtrue->fill(event);
      h_ele_finalselection_mlqtrue->fill(event);
      h_mu_finalselection_mlqtrue->fill(event);
      h_event_finalselection_mlqtrue->fill(event);
      h_lumi_finalselection_mlqtrue->fill(event);
    }
    else
    {
      h_finalselection_mlqfalse->fill(event);
      h_jets_finalselection_mlqfalse->fill(event);
      h_ele_finalselection_mlqfalse->fill(event);
      h_mu_finalselection_mlqfalse->fill(event);
      h_event_finalselection_mlqfalse->fill(event);
      h_lumi_finalselection_mlqfalse->fill(event);
    }





    return true;
  }


  UHH2_REGISTER_ANALYSIS_MODULE(LQTopMuRun2FullSelection)

}
