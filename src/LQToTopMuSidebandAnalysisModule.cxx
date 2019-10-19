// #include <iostream>
// #include <memory>
//
// #include "UHH2/core/include/AnalysisModule.h"
// #include "UHH2/core/include/Event.h"
// #include "UHH2/common/include/CommonModules.h"
// #include "UHH2/common/include/CleaningModules.h"
// #include "UHH2/common/include/JetCorrections.h"
// #include "UHH2/common/include/EventHists.h"
// #include "UHH2/common/include/JetHists.h"
// #include "UHH2/common/include/JetIds.h"
// #include "UHH2/common/include/TopJetIds.h"
// #include "UHH2/common/include/ElectronHists.h"
// #include "UHH2/common/include/ElectronIds.h"
// #include "UHH2/common/include/TopPtReweight.h"
// #include "UHH2/common/include/MuonHists.h"
// #include "UHH2/common/include/MuonIds.h"
// #include "UHH2/common/include/TauHists.h"
// #include "UHH2/common/include/NSelections.h"
// #include "UHH2/LQToTopMu/include/LQToTopMuSelections.h"
// #include "UHH2/LQToTopMu/include/LQToTopMuHists.h"
// #include "UHH2/LQToTopMu/include/LQToTopMuModules.h"
// #include "UHH2/LQToTopMu/include/LQToTopMuEfficiencyHists.h"
// #include "UHH2/LQToTopMu/include/LQToTopMuPDFHists.h"
// #include "UHH2/LQToTopMu/include/LQToTopMuRecoHists.h"
// #include "UHH2/LQToTopMu/include/LQToTopMuMLQRecoHists.h"
// #include "UHH2/LQToTopMu/include/MET2dHists.h"
// #include "UHH2/LQToTopMu/include/HT2dHists.h"
// #include "UHH2/LQToTopMu/include/HypothesisHistsOwn.h"
// #include "UHH2/common/include/PrintingModules.h"
// #include "UHH2/LQToTopMu/include/LQReconstructionHypothesisDiscriminators.h"
// #include "UHH2/LQToTopMu/include/LQReconstruction.h"
// #include "UHH2/common/include/TTbarGen.h"
// #include "UHH2/LQToTopMu/include/LQGen.h"
// #include "UHH2/common/include/MCWeight.h"
// #include "UHH2/common/include/AdditionalSelections.h"
// #include "UHH2/common/include/TriggerSelection.h"
// #include "UHH2/common/include/LuminosityHists.h"
// #include "TFile.h"
// #include "TGraphAsymmErrors.h"
//
//
// using namespace std;
// using namespace uhh2;
//
// namespace uhh2examples {
//
//   class LQToTopMuSidebandAnalysisModule: public AnalysisModule {
//   public:
//
//     explicit LQToTopMuSidebandAnalysisModule(Context & ctx);
//     virtual bool process(Event & event) override;
//
//   private:
//
//     std::unique_ptr<CommonModules> common;
//     std::unique_ptr<AnalysisModule> Muon_printer, Electron_printer, Jet_printer, GenParticles_printer, syst_module, SF_muonID, SF_muonIso, SF_btag, SF_eleReco, SF_eleID;
//     unique_ptr<AnalysisModule> ttgenprod, LQgenprod;
//     unique_ptr<AnalysisModule> sf_top_pt_reweight;
//     unique_ptr<MCMuonScaleFactor> SF_muonTrigger;
//     std::unique_ptr<ElectronTriggerWeights> SF_eleTrigger;
//     unique_ptr<MuonTrkWeights> SF_muonTrk;
//
//     std::unique_ptr<JetCleaner> jetcleaner;
//
//     // declare the Selections to use.
//     std::unique_ptr<Selection>  nele_sel, nbtag_loose_sel, htlept_sel, mttbar_gen_sel, inv_mass_veto, ele_trigger_sel1, ele_trigger_sel2, n_ele_sel, njet_sel, ht_sel, lq_matchable_sel, ttbar_matchable_sel;
//
//     // store the Hists collection as member variables.
//     std::unique_ptr<Hists> h_nocuts, h_mlqreco_nocuts, h_jets_nocuts, h_ele_nocuts, h_mu_nocuts, h_event_nocuts, h_topjets_nocuts, h_tau_nocuts, h_eff_nocuts, h_lumi_nocuts;
//     std::unique_ptr<Hists> h_hyphists, h_1ele, h_jets_1ele, h_ele_1ele, h_mu_1ele, h_event_1ele, h_topjets_1ele, h_tau_1ele, h_eff_1ele, h_lumi_1ele;
//     std::unique_ptr<Hists> h_1bJetLoose, h_mlqreco_1bJetLoose, h_jets_1bJetLoose, h_mu_1bJetLoose, h_ele_1bJetLoose, h_event_1bJetLoose, h_topjets_1bJetLoose, h_tau_1bJetLoose, h_eff_1bJetLoose, h_lumi_1bJetLoose;
//     std::unique_ptr<Hists> h_InvMassVeto, h_mlqreco_InvMassVeto, h_jets_InvMassVeto, h_mu_InvMassVeto, h_ele_InvMassVeto, h_event_InvMassVeto, h_topjets_InvMassVeto, h_tau_InvMassVeto, h_eff_InvMassVeto, h_lumi_InvMassVeto;
//     std::unique_ptr<Hists> h_htlept200, h_mlqreco_htlept200, h_jets_htlept200, h_ele_htlept200, h_mu_htlept200, h_event_htlept200, h_topjets_htlept200, h_tau_htlept200, h_btageff_htlept200, h_eff_htlept200, h_lumi_htlept200;
//     unique_ptr<Hists> h_ttbarmatchable, h_mlqreco_ttbarmatchable, h_jets_ttbarmatchable, h_ele_ttbarmatchable, h_mu_ttbarmatchable, h_event_ttbarmatchable, h_topjets_ttbarmatchable, h_tau_ttbarmatchable, h_eff_ttbarmatchable, h_lumi_ttbarmatchable;
//     unique_ptr<Hists> h_lqmatchable, h_mlqreco_lqmatchable, h_jets_lqmatchable, h_ele_lqmatchable, h_mu_lqmatchable, h_event_lqmatchable, h_topjets_lqmatchable, h_tau_lqmatchable, h_eff_lqmatchable, h_lumi_lqmatchable;
//     unique_ptr<Hists> h_correctlymatched, h_mlqreco_correctlymatched, h_jets_correctlymatched, h_ele_correctlymatched, h_mu_correctlymatched, h_event_correctlymatched, h_topjets_correctlymatched, h_tau_correctlymatched, h_eff_correctlymatched, h_lumi_correctlymatched;
//     std::unique_ptr<Hists> h_finalSelection, h_mlqreco_finalSelection, h_jets_finalSelection, h_ele_finalSelection, h_mu_finalSelection, h_event_finalSelection, h_topjets_finalSelection, h_tau_finalSelection, h_eff_finalSelection, h_lumi_finalSelection;
//     std::unique_ptr<Hists> h_ht_InvMassVeto, h_ht_finalSelection;
//     std::unique_ptr<Hists> h_Sideband;
//     std::unique_ptr<Hists> h_PDF_variations;
//
//
//     MuonId MuId;
//     ElectronId EleId, EleId_highpt;
//     JetId Btag_loose;
//     CSVBTag::wp wp_btag_loose;
//
//     bool do_scale_variation, is_mc, do_pdf_variations, is_mu_e, is_e_e, apply_EleTriggerSF, apply_alpha;
//
//     Event::Handle<TTbarGen> h_ttbargen;
//     Event::Handle<LQGen> h_LQLQbargen;
//
//
//     std::unique_ptr<TFile> file_alpha;
//     std::unique_ptr<TGraphAsymmErrors> alpha;
//     //std::unique_ptr<TH1D> norm;
//     double EleTriggerSF;
//     string filepath_alpha, Sys_MuonID, Sys_MuonTrk, Sys_BTag, Sys_MuonTrigger, Sys_MuonIso, Sys_PU, Sys_EleID, Sys_EleReco, Sys_EleTrigger, Sys_TTbar, Sys_DY, Sys_ST, Sys_DB, Sys_QCD, Sys_WJ, Sys_TTV;
//     TString dataset_version, scalevariation_process;
//
//     vector<unique_ptr<AnalysisModule>> recomodules, muonic_recomodules;
//     uhh2::Event::Handle<vector<LQReconstructionHypothesis>> h_hyps, h_muonic_hyps;
//
//   };
//
//
//   LQToTopMuSidebandAnalysisModule::LQToTopMuSidebandAnalysisModule(Context & ctx){
//
//     cout << "Hello World from LQToTopMuSidebandAnalysisModule!" << endl;
//
//     for(auto & kv : ctx.get_all()){
//       cout << " " << kv.first << " = " << kv.second << endl;
//     }
//
//     scalevariation_process   = ctx.get("ScaleVariationProcess");
//     do_scale_variation       = (ctx.get("ScaleVariationMuR") == "up" || ctx.get("ScaleVariationMuR") == "down") || (ctx.get("ScaleVariationMuF") == "up" || ctx.get("ScaleVariationMuF") == "down");
//     if(do_scale_variation && (scalevariation_process != "ttbar") && (scalevariation_process != "dy") && (scalevariation_process != "st") && (scalevariation_process != "wj") && (scalevariation_process != "db") && (scalevariation_process != "ttv")) throw runtime_error("In LQToTopMuAnalysisModule.cxx: Invalid process specified for 'ScaleVariationProcess'.");
//     do_pdf_variations = ctx.get("b_PDFUncertainties") == "true";
//     dataset_version = ctx.get("dataset_version");
//
//     is_mc = ctx.get("dataset_type") == "MC";
//     is_mu_e = (ctx.get("channel") == "mu_e" || ctx.get("channel") == "e_mu");
//     is_e_e = ctx.get("channel") == "e_e";
//     if((!is_mu_e && !is_e_e)) throw runtime_error("In SidebandPreselectionModule: Invalid definition of 'channel' in config file, must be 'mu_e', 'e_mu', or 'e_e'.");
//     apply_EleTriggerSF = (is_mc && is_e_e);
//
//     filepath_alpha = ctx.get("filepath_alpha");
//     apply_alpha = ctx.get("Apply_Alpha") == "true";
//     Sys_MuonID = ctx.get("Systematic_MuonID");
//     Sys_MuonTrigger = ctx.get("Systematic_MuonTrigger");
//     Sys_MuonIso = ctx.get("Systematic_MuonIso");
//     Sys_MuonTrk = ctx.get("Systematic_MuonTrk");
//     Sys_EleID = ctx.get("Systematic_EleID");
//     Sys_EleReco = ctx.get("Systematic_EleReco");
//     Sys_EleTrigger = ctx.get("Systematic_EleTrigger");
//     Sys_BTag = ctx.get("Systematic_BTag");
//     Sys_PU = ctx.get("Systematic_PU");
//     Sys_TTbar = ctx.get("Systematic_TTbar");
//     Sys_DY = ctx.get("Systematic_DY");
//     Sys_ST = ctx.get("Systematic_ST");
//     Sys_DB = ctx.get("Systematic_DB");
//     Sys_WJ = ctx.get("Systematic_WJ");
//     Sys_QCD = ctx.get("Systematic_QCD");
//     Sys_TTV = ctx.get("Systematic_TTV");
//     const char* c_filepath_alpha = filepath_alpha.c_str();
//
//     if(apply_alpha){
//       file_alpha.reset(new TFile(c_filepath_alpha,"READ"));
//       alpha.reset((TGraphAsymmErrors*)file_alpha->Get("Graph"));
//       //norm.reset((TH1D*)file_alpha->Get("h_normalization"));
//       //if(!norm) norm.reset((TH1D*)file_alpha->Get("h_normalization_syst_up"));
//       //if(!norm) norm.reset((TH1D*)file_alpha->Get("h_normalization_syst_dn"));
//     }
//
//     // 1. setup other modules. CommonModules and the JetCleaner:
//     Jet_printer.reset(new JetPrinter("Jet-Printer", 0));
//     Electron_printer.reset(new ElectronPrinter("Electron-Printer"));
//     Muon_printer.reset(new MuonPrinter("Muon-Printer"));
//     GenParticles_printer.reset(new GenParticlesPrinter(ctx));
//     MuId = AndId<Muon>(MuonIDTight(),PtEtaCut(30.0, 2.4),MuonIso(0.15));
//     EleId = AndId<Electron>(ElectronID_Spring16_tight,PtEtaCut(30.0, 2.4));
//     EleId_highpt = AndId<Electron>(ElectronID_Spring16_tight,PtEtaCut(120.0, 2.4)); //110
//
//
//     Btag_loose = CSVBTag(CSVBTag::WP_LOOSE);
//     wp_btag_loose = CSVBTag::WP_LOOSE;
//     jetcleaner.reset(new JetCleaner(ctx,30.0, 2.4));
//
//     common.reset(new CommonModules());
//     common->disable_lumisel();
//     common->disable_jersmear();
//     common->disable_jec();
//     common->set_electron_id(EleId);
//     common->set_muon_id(MuId);
//     common->init(ctx,Sys_PU);
//
//     if(is_mu_e){
//       SF_muonID.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/reimersa/CMSSW_8_0_24_patch1/src/UHH2/common/data/MuonID_EfficienciesAndSF_average_RunBtoH.root", "MC_NUM_TightID_DEN_genTracks_PAR_pt_eta", 1., "tightID", true, "nominal"));
//       SF_muonTrigger.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/reimersa/CMSSW_8_0_24_patch1/src/UHH2/common/data/MuonTrigger_EfficienciesAndSF_average_RunBtoH.root", "IsoMu24_OR_IsoTkMu24_PtEtaBins", 0.5, "trigger", true, "nominal"));
//       SF_muonIso.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/reimersa/CMSSW_8_0_24_patch1/src/UHH2/common/data/MuonIso_EfficienciesAndSF_average_RunBtoH.root", "TightISO_TightID_pt_eta", 1., "iso", true, "nominal"));
//       SF_muonTrk.reset(new MuonTrkWeights(ctx, "/nfs/dust/cms/user/reimersa/CMSSW_8_0_24_patch1/src/UHH2/common/data/Tracking_EfficienciesAndSF_BCDEFGH.root", Sys_MuonTrk));
//     }
//
//     SF_eleReco.reset(new MCElecScaleFactor(ctx, "/nfs/dust/cms/user/reimersa/CMSSW_8_0_24_patch1/src/UHH2/common/data/egammaEffi.txt_EGM2D_RecEff_Moriond17.root", 1, "", Sys_EleReco));
//     SF_eleID.reset(new MCElecScaleFactor(ctx, "/nfs/dust/cms/user/reimersa/CMSSW_8_0_24_patch1/src/UHH2/common/data/egammaEffi.txt_EGM2D_CutBased_Tight_ID.root", 1, "", Sys_EleID));
//     if(is_e_e) SF_eleTrigger.reset(new ElectronTriggerWeights(ctx, "/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TagProbe/Optimization/35867fb_Iso27_NonIso115/ElectronEfficiencies.root", Sys_EleTrigger));
//
//     SF_btag.reset(new MCBTagScaleFactor(ctx,wp_btag_loose,"jets",Sys_BTag));
//
//     //sf_top_pt_reweight.reset(new TopPtReweight(ctx, 0.0615, -0.0005, "", "", true, 1.0302)); // https://twiki.cern.ch/twiki/bin/view/CMS/TopPtReweighting
//
//
//     h_hyps = ctx.get_handle<vector<LQReconstructionHypothesis>>("HighMassLQReconstruction");
//     h_muonic_hyps = ctx.get_handle<vector<LQReconstructionHypothesis>>("HighMassMuonicLQReconstruction");
//
//     if(is_mc){
//       ttgenprod.reset(new TTbarGenProducer(ctx, "ttbargen", false));
//       h_ttbargen = ctx.get_handle<TTbarGen>("ttbargen");
//       LQgenprod.reset(new LQGenProducer(ctx, "LQLQbargen", false));
//       h_LQLQbargen = ctx.get_handle<LQGen>("LQLQbargen");
//     }
//
//
//
//
//     // 2. set up selections
//     //Selection
//     njet_sel.reset(new NJetSelection(2, -1));
//     ht_sel.reset(new HtSelection(350., -1));
//     nbtag_loose_sel.reset(new NJetSelection(1, -1, Btag_loose));
//     htlept_sel.reset(new HTLeptSelection(200., -1));
//     if(is_mu_e) inv_mass_veto.reset(new InvMassMuEleVeto(0.,111.));
//     else        inv_mass_veto.reset(new InvMassEleEleVeto(0.,111.));
//     ele_trigger_sel1.reset(new TriggerSelection("HLT_Ele27_WPTight_Gsf_v*"));
//     ele_trigger_sel2.reset(new TriggerSelection("HLT_Ele105_CaloIdVT_GsfTrkIdT_v*"));
//     n_ele_sel.reset(new NElectronSelection(2, -1, EleId_highpt));
//     lq_matchable_sel.reset(new LQSemiLepMatchable());
//     ttbar_matchable_sel.reset(new TTbarSemiLepMatchable());
//
//     //make reconstruction hypotheses
//     recomodules.emplace_back(new LQPrimaryLepton(ctx));
//     recomodules.emplace_back(new HighMassLQReconstruction(ctx,LQNeutrinoReconstruction));
//     recomodules.emplace_back(new LQChi2Discriminator(ctx,"HighMassLQReconstruction"));
//     if(is_mc) recomodules.emplace_back(new LQCorrectMatchDiscriminator(ctx,"HighMassLQReconstruction"));
//     muonic_recomodules.emplace_back(new HighMassMuonicLQReconstruction(ctx,LQNeutrinoReconstruction));
//     muonic_recomodules.emplace_back(new LQChi2Discriminator(ctx,"HighMassMuonicLQReconstruction"));
//     if(is_mc) muonic_recomodules.emplace_back(new LQCorrectMatchDiscriminator(ctx,"HighMassMuonicLQReconstruction"));
//
//     //systematics modules
//     syst_module.reset(new MCScaleVariation(ctx));
//
//     // 3. Set up Hists classes:
//     h_nocuts.reset(new LQToTopMuHists(ctx, "NoCuts"));
//     h_mlqreco_nocuts.reset(new LQToTopMuMLQRecoHists(ctx, "MLQReco_NoCuts"));
//     h_jets_nocuts.reset(new JetHists(ctx, "Jets_NoCuts"));
//     h_ele_nocuts.reset(new ElectronHists(ctx, "Ele_NoCuts"));
//     h_mu_nocuts.reset(new MuonHists(ctx, "Mu_NoCuts"));
//     h_event_nocuts.reset(new EventHists(ctx, "Event_NoCuts"));
//     h_topjets_nocuts.reset(new TopJetHists(ctx, "TopJets_NoCuts"));
//     h_eff_nocuts.reset(new LQToTopMuEfficiencyHists(ctx, "Eff_NoCuts"));
//     h_lumi_nocuts.reset(new LuminosityHists(ctx, "Lumi_NoCuts"));
//
//     h_hyphists.reset(new HypothesisHistsOwn(ctx, "CorrectMatch_Hists", "HighMassLQReconstruction", "CorrectMatch"));
//
//     h_1ele.reset(new LQToTopMuHists(ctx, "1Ele"));
//     h_jets_1ele.reset(new JetHists(ctx, "Jets_1Ele"));
//     h_ele_1ele.reset(new ElectronHists(ctx, "Ele_1Ele"));
//     h_mu_1ele.reset(new MuonHists(ctx, "Mu_1Ele"));
//     h_event_1ele.reset(new EventHists(ctx, "Event_1Ele"));
//     h_topjets_1ele.reset(new TopJetHists(ctx, "TopJets_1Ele"));
//     h_eff_1ele.reset(new LQToTopMuEfficiencyHists(ctx, "Eff_1Ele"));
//     h_lumi_1ele.reset(new LuminosityHists(ctx, "Lumi_1Ele"));
//
//     h_1bJetLoose.reset(new LQToTopMuHists(ctx, "1bJetLoose"));
//     h_mlqreco_1bJetLoose.reset(new LQToTopMuMLQRecoHists(ctx, "MLQReco_1bJetLoose"));
//     h_jets_1bJetLoose.reset(new JetHists(ctx, "Jets_1bJetLoose"));
//     h_ele_1bJetLoose.reset(new ElectronHists(ctx, "Ele_1bJetLoose"));
//     h_mu_1bJetLoose.reset(new MuonHists(ctx, "Mu_1bJetLoose"));
//     h_event_1bJetLoose.reset(new EventHists(ctx, "Event_1bJetLoose"));
//     h_topjets_1bJetLoose.reset(new TopJetHists(ctx, "TopJets_1bJetLoose"));
//     h_eff_1bJetLoose.reset(new LQToTopMuEfficiencyHists(ctx, "Eff_1bJetLoose"));
//     h_lumi_1bJetLoose.reset(new LuminosityHists(ctx, "Lumi_1bJetLoose"));
//
//     h_InvMassVeto.reset(new LQToTopMuHists(ctx, "InvMassVeto"));
//     h_mlqreco_InvMassVeto.reset(new LQToTopMuMLQRecoHists(ctx, "MLQReco_InvMassVeto"));
//     h_jets_InvMassVeto.reset(new JetHists(ctx, "Jets_InvMassVeto"));
//     h_ele_InvMassVeto.reset(new ElectronHists(ctx, "Ele_InvMassVeto"));
//     h_mu_InvMassVeto.reset(new MuonHists(ctx, "Mu_InvMassVeto"));
//     h_event_InvMassVeto.reset(new EventHists(ctx, "Event_InvMassVeto"));
//     h_topjets_InvMassVeto.reset(new TopJetHists(ctx, "TopJets_InvMassVeto"));
//     h_eff_InvMassVeto.reset(new LQToTopMuEfficiencyHists(ctx, "Eff_InvMassVeto"));
//     h_lumi_InvMassVeto.reset(new LuminosityHists(ctx, "Lumi_InvMassVeto"));
//
//
//     h_htlept200.reset(new LQToTopMuHists(ctx, "HTLept200"));
//     h_mlqreco_htlept200.reset(new LQToTopMuMLQRecoHists(ctx, "MLQReco_HTLept200"));
//     h_jets_htlept200.reset(new JetHists(ctx, "Jets_HTLept200"));
//     h_ele_htlept200.reset(new ElectronHists(ctx, "Ele_HTLept200"));
//     h_mu_htlept200.reset(new MuonHists(ctx, "Mu_HTLept200"));
//     h_topjets_htlept200.reset(new TopJetHists(ctx, "TopJets_HTLept200"));
//     h_event_htlept200.reset(new EventHists(ctx, "Event_HTLept200"));
//     h_eff_htlept200.reset(new LQToTopMuEfficiencyHists(ctx, "Eff_HTLept200"));
//     h_lumi_htlept200.reset(new LuminosityHists(ctx, "Lumi_HTLept200"));
//     h_btageff_htlept200.reset(new BTagMCEfficiencyHists(ctx, "BTagEff_HTLept200",wp_btag_loose));
//
//
//     h_ttbarmatchable.reset(new LQToTopMuHists(ctx, "TTbar_Matchable"));
//     h_mlqreco_ttbarmatchable.reset(new LQToTopMuMLQRecoHists(ctx, "MLQReco_TTbar_Matchable"));
//     h_jets_ttbarmatchable.reset(new JetHists(ctx, "Jets_TTbar_Matchable"));
//     h_ele_ttbarmatchable.reset(new ElectronHists(ctx, "Ele_TTbar_Matchable"));
//     h_mu_ttbarmatchable.reset(new MuonHists(ctx, "Mu_TTbar_Matchable"));
//     h_topjets_ttbarmatchable.reset(new TopJetHists(ctx, "TopJets_TTbar_Matchable"));
//     h_event_ttbarmatchable.reset(new EventHists(ctx, "Event_TTbar_Matchable"));
//     h_eff_ttbarmatchable.reset(new LQToTopMuEfficiencyHists(ctx, "Eff_TTbar_Matchable"));
//     h_lumi_ttbarmatchable.reset(new LuminosityHists(ctx, "Lumi_TTbar_Matchable"));
//
//
//     h_lqmatchable.reset(new LQToTopMuHists(ctx, "LQ_Matchable"));
//     h_mlqreco_lqmatchable.reset(new LQToTopMuMLQRecoHists(ctx, "MLQReco_LQ_Matchable"));
//     h_jets_lqmatchable.reset(new JetHists(ctx, "Jets_LQ_Matchable"));
//     h_ele_lqmatchable.reset(new ElectronHists(ctx, "Ele_LQ_Matchable"));
//     h_mu_lqmatchable.reset(new MuonHists(ctx, "Mu_LQ_Matchable"));
//     h_topjets_lqmatchable.reset(new TopJetHists(ctx, "TopJets_LQ_Matchable"));
//     h_event_lqmatchable.reset(new EventHists(ctx, "Event_LQ_Matchable"));
//     h_eff_lqmatchable.reset(new LQToTopMuEfficiencyHists(ctx, "Eff_LQ_Matchable"));
//     h_lumi_lqmatchable.reset(new LuminosityHists(ctx, "Lumi_LQ_Matchable"));
//
//
//     h_correctlymatched.reset(new LQToTopMuHists(ctx, "CorrectlyMatched"));
//     h_mlqreco_correctlymatched.reset(new LQToTopMuMLQRecoHists(ctx, "MLQReco_CorrectlyMatched"));
//     h_jets_correctlymatched.reset(new JetHists(ctx, "Jets_CorrectlyMatched"));
//     h_ele_correctlymatched.reset(new ElectronHists(ctx, "Ele_CorrectlyMatched"));
//     h_mu_correctlymatched.reset(new MuonHists(ctx, "Mu_CorrectlyMatched"));
//     h_topjets_correctlymatched.reset(new TopJetHists(ctx, "TopJets_CorrectlyMatched"));
//     h_event_correctlymatched.reset(new EventHists(ctx, "Event_CorrectlyMatched"));
//     h_eff_correctlymatched.reset(new LQToTopMuEfficiencyHists(ctx, "Eff_CorrectlyMatched"));
//     h_lumi_correctlymatched.reset(new LuminosityHists(ctx, "Lumi_CorrectlyMatched"));
//
//
//     h_finalSelection.reset(new LQToTopMuHists(ctx, "FinalSelection"));
//     h_mlqreco_finalSelection.reset(new LQToTopMuMLQRecoHists(ctx, "MLQReco_FinalSelection"));
//     h_jets_finalSelection.reset(new JetHists(ctx, "Jets_FinalSelection"));
//     h_ele_finalSelection.reset(new ElectronHists(ctx, "Ele_FinalSelection"));
//     h_mu_finalSelection.reset(new MuonHists(ctx, "Mu_FinalSelection"));
//     h_topjets_finalSelection.reset(new TopJetHists(ctx, "TopJets_FinalSelection"));
//     h_event_finalSelection.reset(new EventHists(ctx, "Event_FinalSelection"));
//     h_eff_finalSelection.reset(new LQToTopMuEfficiencyHists(ctx, "Eff_FinalSelection"));
//     h_lumi_finalSelection.reset(new LuminosityHists(ctx, "Lumi_FinalSelection"));
//     h_ht_finalSelection.reset(new HT2dHists(ctx, "HT2d_FinalSelection"));
//     h_PDF_variations.reset(new LQToTopMuPDFHists(ctx, "PDF_variations", true, do_pdf_variations));
//     h_Sideband.reset(new LQToTopMuHists(ctx, "Sideband_weights_applied"));
//
//
//   }
//
//
//   bool LQToTopMuSidebandAnalysisModule::process(Event & event) {
//
//     if(is_mc){
//       double factor_xsec = -1;
//       int control = (dataset_version.Contains("TTbar") && Sys_TTbar != "nominal") + (dataset_version.Contains("DYJets") && Sys_DY != "nominal") + (dataset_version.Contains("SingleTop") && Sys_ST != "nominal") + (dataset_version.Contains("WJets") && Sys_WJ != "nominal") + (dataset_version.Contains("Diboson") && Sys_DB != "nominal") + (dataset_version.Contains("QCD") && Sys_QCD != "nominal")+ (dataset_version.Contains("TTV") && Sys_TTV != "nominal");
//       if(!(control == 0 || control == 1)) throw runtime_error("In LQToTopMuSidebandAnalysisModule.cxx: More than one rate systematic is set to something different than 'nominal'");
//
//       if(control == 0) factor_xsec = 0;
//       else if(control == 1){
// 	if(dataset_version.Contains("TTbar") && Sys_TTbar != "nominal")       factor_xsec = 0.056;
// 	else if(dataset_version.Contains("DYJets") && Sys_DY != "nominal")    factor_xsec = 0.1;
// 	else if(dataset_version.Contains("SingleTop") && Sys_ST != "nominal") factor_xsec = 0.1;
// 	else if(dataset_version.Contains("WJets") && Sys_WJ != "nominal")     factor_xsec = 0.1;
// 	else if(dataset_version.Contains("Diboson") && Sys_DB != "nominal")   factor_xsec = 0.2;
// 	else if(dataset_version.Contains("QCD") && Sys_QCD != "nominal")      factor_xsec = 1;
// 	else if(dataset_version.Contains("TTV") && Sys_TTV != "nominal")      factor_xsec = 0.25;
// 	else if(dataset_version.Contains("LQ"))                               factor_xsec = 0;
//       }
//       double sf_xsec = 1;
//       if(Sys_TTbar == "up" || Sys_DY == "up"|| Sys_ST == "up"|| Sys_DB == "up"|| Sys_WJ == "up"|| Sys_QCD == "up"|| Sys_TTV == "up") sf_xsec += factor_xsec;
//       else if(Sys_TTbar == "down" || Sys_DY == "down"|| Sys_ST == "down"|| Sys_DB == "down"|| Sys_WJ == "down"|| Sys_QCD == "down"|| Sys_TTV == "down") sf_xsec -= factor_xsec;
//       else if(control != 0) throw runtime_error("In LQToTopMuAnalysisModule.cxx: Invalid direction for 'Sys_Rate_YYY' specified.");
//
//       event.weight *= sf_xsec;
//     }
//
//
//
//
//     //check for at least 1 muon pair with opposite charge
//     bool charge_opposite = false;
//     for(unsigned int i=0; i<event.muons->size(); i++){
//       for(unsigned int j=0; j<event.muons->size(); j++){
// 	if(j>i){
// 	  if(event.muons->at(i).charge() != event.muons->at(j).charge()) {
// 	    charge_opposite = true;
// 	  }
// 	}
//       }
//     }
//     double sum_mu_charge = 0;
//     for(const auto & mu : *event.muons) sum_mu_charge += mu.charge();
//     bool reconstruct_mlq_mu = (event.electrons->size() == 0 && event.muons->size() == 3 && fabs(sum_mu_charge) == 1);
//     bool reconstruct_mlq_ele = (event.electrons->size() >= 1 && event.muons->size() >= 2 && charge_opposite);
//     bool reconstruct_mlq = reconstruct_mlq_ele || reconstruct_mlq_mu;
//
//
//     if(is_mc){
//       ttgenprod->process(event);
//       LQgenprod->process(event);
//     }
//
//     double CorrectMatchDiscriminator = 99;
//     // MLQ reco
//     if(reconstruct_mlq_ele){
//       for(auto & m : recomodules){
// 	m->process(event);
//       }
//
//       if(is_mc){
// 	std::vector<LQReconstructionHypothesis> hyps = event.get(h_hyps);
// 	const LQReconstructionHypothesis* hyp = get_best_hypothesis( hyps, "Chi2" );
// 	CorrectMatchDiscriminator = hyp->discriminator("CorrectMatch");
//       }
//
//     }
//
//
//     if(reconstruct_mlq_mu){
//       for(auto & m : muonic_recomodules){
// 	m->process(event);
//       }
//       if(is_mc){
// 	std::vector<LQReconstructionHypothesis> hyps = event.get(h_muonic_hyps);
// 	const LQReconstructionHypothesis* hyp = get_best_hypothesis( hyps, "Chi2" );
// 	CorrectMatchDiscriminator = hyp->discriminator("CorrectMatch");
//       }
//     }
//
//
//
//     //apply muon & electron SFs as in preselection
//     if(is_mu_e){
//       SF_muonTrigger->process_onemuon(event,0);
//       SF_muonID->process(event);
//       SF_muonIso->process(event);
//       SF_muonTrk->process(event);
//     }
//     SF_eleReco->process(event);
//     SF_eleID->process(event);
//     if(apply_EleTriggerSF){
//       SF_eleTrigger->process(event);
//     }
//
//     if(do_scale_variation){
//       if((dataset_version.Contains("TTbar") && scalevariation_process == "ttbar") || (dataset_version.Contains("DYJets") && scalevariation_process == "dy") || (dataset_version.Contains("SingleTop") && scalevariation_process == "st") || (dataset_version.Contains("WJets") && scalevariation_process == "wj") || (dataset_version.Contains("Diboson") && scalevariation_process == "db") || (dataset_version.Contains("TTV") && scalevariation_process == "ttv")){
// 	if(event.genInfo->systweights().size() < 10 && dataset_version.Contains("Diboson")) cout << "SystWeight size: " << event.genInfo->systweights().size() << endl;
// 	else syst_module->process(event);
//       }
//       //Signals not needed for scale variation!
//       else if(dataset_version.Contains("LQtoT")) return false;
//     }
//
//     //sf_top_pt_reweight->process(event);
//
//     bool pass_common = common->process(event);
//     if(!pass_common) return false;
//     jetcleaner->process(event);
//     if(!njet_sel->passes(event)) return false;
//     if(!ht_sel->passes(event))   return false;
//
//     //HT
//     auto met = event.met->pt();
//     double ht = 0.0;
//     double ht_jets = 0.0;
//     double ht_lep = 0.0;
//     for(const auto & jet : *event.jets){
//       ht_jets += jet.pt();
//     }
//     for(const auto & electron : *event.electrons){
//       ht_lep += electron.pt();
//     }
//     for(const auto & muon : *event.muons){
//       ht_lep += muon.pt();
//     }
//     ht = ht_lep + ht_jets + met;
//
//
//     //if(is_mc)GenParticles_printer->process(event);
//
//     if(event.electrons->size() > 0 && event.muons->size() >= 2){
//       for(auto & m : recomodules){
// 	m->process(event);
//       }
//     }
//
//     h_nocuts->fill(event);
//     h_mlqreco_nocuts->fill(event);
//     h_jets_nocuts->fill(event);
//     h_ele_nocuts->fill(event);
//     h_mu_nocuts->fill(event);
//     h_event_nocuts->fill(event);
//     h_topjets_nocuts->fill(event);
//     h_eff_nocuts->fill(event);
//     h_lumi_nocuts->fill(event);
//     if(reconstruct_mlq_ele && is_mc) h_hyphists->fill(event);
//
//     //1 bTag1 loose
//     if(!nbtag_loose_sel->passes(event)) return false;
//     SF_btag->process(event);
//
//     h_jets_1bJetLoose->fill(event);
//     h_1bJetLoose->fill(event);
//     h_ele_1bJetLoose->fill(event);
//     h_mu_1bJetLoose->fill(event);
//     h_event_1bJetLoose->fill(event);
//     h_topjets_1bJetLoose->fill(event);
//     h_eff_1bJetLoose->fill(event);
//     h_lumi_1bJetLoose->fill(event);
//
//     //InvMassVeto
//     if(!inv_mass_veto->passes(event)) return false;
//     h_jets_InvMassVeto->fill(event);
//     h_InvMassVeto->fill(event);
//     h_mlqreco_InvMassVeto->fill(event);
//     h_ele_InvMassVeto->fill(event);
//     h_mu_InvMassVeto->fill(event);
//     h_event_InvMassVeto->fill(event);
//     h_topjets_InvMassVeto->fill(event);
//     h_eff_InvMassVeto->fill(event);
//     h_lumi_InvMassVeto->fill(event);
//
//     if(!htlept_sel->passes(event)) return false;
//     h_jets_htlept200->fill(event);
//     h_htlept200->fill(event);
//     h_mlqreco_htlept200->fill(event);
//     h_ele_htlept200->fill(event);
//     h_mu_htlept200->fill(event);
//     h_event_htlept200->fill(event);
//     h_topjets_htlept200->fill(event);
//     h_eff_htlept200->fill(event);
//     h_lumi_htlept200->fill(event);
//     h_btageff_htlept200->fill(event);
//
//     if(is_mc && reconstruct_mlq && !dataset_version.Contains("LQtoTTau")&& !dataset_version.Contains("LQtoTLep")){
//       if(ttbar_matchable_sel->passes(event)){
// 	h_ttbarmatchable->fill(event);
// 	h_mlqreco_ttbarmatchable->fill(event);
// 	h_jets_ttbarmatchable->fill(event);
// 	h_ele_ttbarmatchable->fill(event);
// 	h_mu_ttbarmatchable->fill(event);
// 	h_event_ttbarmatchable->fill(event);
// 	h_topjets_ttbarmatchable->fill(event);
// 	h_eff_ttbarmatchable->fill(event);
// 	h_lumi_ttbarmatchable->fill(event);
//       }
//
//       if(lq_matchable_sel->passes(event)){
// 	h_lqmatchable->fill(event);
// 	h_mlqreco_lqmatchable->fill(event);
// 	h_jets_lqmatchable->fill(event);
// 	h_ele_lqmatchable->fill(event);
// 	h_mu_lqmatchable->fill(event);
// 	h_event_lqmatchable->fill(event);
// 	h_topjets_lqmatchable->fill(event);
// 	h_eff_lqmatchable->fill(event);
// 	h_lumi_lqmatchable->fill(event);
//
// 	if(CorrectMatchDiscriminator < 99){
// 	  h_correctlymatched->fill(event);
// 	  h_mlqreco_correctlymatched->fill(event);
// 	  h_jets_correctlymatched->fill(event);
// 	  h_ele_correctlymatched->fill(event);
// 	  h_mu_correctlymatched->fill(event);
// 	  h_event_correctlymatched->fill(event);
// 	  h_topjets_correctlymatched->fill(event);
// 	  h_eff_correctlymatched->fill(event);
// 	  h_lumi_correctlymatched->fill(event);
// 	}
//
//       }
//     }
//
//     h_finalSelection->fill(event);
//     h_mlqreco_finalSelection->fill(event);
//     h_jets_finalSelection->fill(event);
//     h_ele_finalSelection->fill(event);
//     h_mu_finalSelection->fill(event);
//     h_event_finalSelection->fill(event);
//     h_topjets_finalSelection->fill(event);
//     h_eff_finalSelection->fill(event);
//     h_lumi_finalSelection->fill(event);
//     h_ht_finalSelection->fill(event);
//
//     h_PDF_variations->fill(event);
//
//
//     double original_weight = event.weight;
//     if(apply_alpha){
//       double d_alpha = alpha->Eval(ht);
//       //double d_norm = norm->GetBinContent(1);
//       //double sideband_weight = d_alpha * d_norm;
//
//       //change weights
//       //event.weight *= sideband_weight;
//       event.weight *= d_alpha;
//     }
//     h_Sideband->fill(event);
//     //restore weights
//     event.weight = original_weight;
//
//     return true;
//   }
//
//
//   UHH2_REGISTER_ANALYSIS_MODULE(LQToTopMuSidebandAnalysisModule)
// }

































// #include "UHH2/core/include/Event.h"
// #include "UHH2/common/include/Utils.h"
// #include "UHH2/LQToTopMu/include/LQToTopMuModules.h"
// #include "UHH2/common/include/JetIds.h"
// #include "UHH2/common/include/JetCorrections.h"
// #include <TH1D.h>
//
//
// using namespace uhh2;
// using namespace std;
//
//
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
//
//
//
// JetLeptonOverlapCleaner::JetLeptonOverlapCleaner(double RJet_): RJet(RJet_){}
//
// bool JetLeptonOverlapCleaner::process(Event & event){
//
//    vector<Jet> result_jets;
//    vector<Muon> result_muons;
//    vector<Electron> result_electrons;
//
//    cout << "############# Start of a new event ##################" << endl << endl;
//
//    cout << "------------- first deal with the jets --------------" << endl << endl;
//    //Check, if jets are identical to leptons
//    //if so, kick them out
//    for(const Jet & thisjet : *event.jets){
//      auto thisjet_v4_raw = thisjet.v4() * thisjet.JEC_factor_raw();
//      bool keep_thisjet = true;
//      cout << "new jet: " << endl;
//      for(const Muon & thismu : *event.muons){
//        //if true, jet is considered identical to the muon
//        cout << "deltaR between this muon and the jet: " << deltaR(thismu,thisjet) << ", mupt/jetptraw: " << thismu.pt()/thisjet_v4_raw.pt() << endl;
//        if (deltaR(thismu,thisjet) <= 0.05 && (thismu.pt() <= thisjet_v4_raw.pt()*1.1 && thismu.pt() >= thisjet_v4_raw.pt()*0.8)){
// 	 keep_thisjet = false;
// 	 cout << "at least bc/ of this, this jet will not be kept" << endl;
//        }
//      }//end muons
//      for(const Electron & thisele : *event.electrons){
//        //if true, jet is considered identical to the electron
//       cout << "deltaR between this electron and the jet: " << deltaR(thisele,thisjet) << ", elept/jetptraw: " << thisele.pt()/thisjet_v4_raw.pt() << endl;
//        if (deltaR(thisele,thisjet) <= 0.05 && (thisele.pt() <= thisjet_v4_raw.pt()*1.1 && thisele.pt() >= thisjet_v4_raw.pt()*0.8)){
// 	 keep_thisjet = false;
// 	 cout << "at least bc/ of this, this jet will not be kept" << endl;
//        }
//      }//end electrons
//      if(keep_thisjet) result_jets.push_back(thisjet);
//    }
//    //final collection of jets that are not identical to leptons
//    swap(*event.jets,result_jets);
//
//
//
//    //now check, if leptons are lying inside a jet
//    //if so, kick the lepton out
//
//    //first for muons
//    cout << "----------------- coming to the muons ------------" << endl << endl;
//    for(const Muon & thismu : *event.muons){
//      bool keep_thismu = true;
//      cout << "new muon: " << endl;
//      //already kicked out jets that actually were leptons
//      for(const Jet & thisjet : *event.jets){
//        cout << "deltaR between this muon and a jet: " <<  deltaR(thismu,thisjet) << endl;
//        if(deltaR(thismu,thisjet) <= RJet) {
// 	 keep_thismu = false;
// 	 cout << "at least bc/ of this, this muon will not be kept" << endl;
//        }
//      }
//      if(keep_thismu) result_muons.push_back(thismu);
//    }
//
//    //the same for electrons
//    cout << "----------------- coming to the electrons ------------" << endl << endl;
//    for(const Electron & thisele : *event.electrons){
//      bool keep_thisele = true;
//      cout << "new electron: " << endl;
//      //already kicked out jets that actually were leptons
//      for(const Jet & thisjet : *event.jets){
//        cout << "deltaR between this electron and a jet: " <<  deltaR(thisele,thisjet) << endl;
//        if(deltaR(thisele,thisjet) <= RJet){
// 	 keep_thisele = false;
// 	 cout << "at least bc/ of this, this muon will not be kept" << endl;
//        }
//      }
//      if(keep_thisele) result_electrons.push_back(thisele);
//    }
//
//    swap(*event.muons,result_muons);
//    swap(*event.electrons,result_electrons);
//
//    cout << "############### End of the event ###################" << endl << endl;
//    return true;
// }
//
//
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
//     if(SysDirection == "up"){
//       stat_mc = Eff_lowpt_MC->GetErrorYlow(idx);
//       stat_data = Eff_lowpt_DATA->GetErrorYhigh(idx);
//       total_syst_mc = sqrt(pow(stat_mc,2) + pow(tp,2));
//       total_syst_data = sqrt(pow(stat_data,2) + pow(tp,2));
//
//       eff_mc -= total_syst_mc;
//       eff_data += total_syst_data;
//     }
//     else if(SysDirection == "down"){
//       stat_mc = Eff_lowpt_MC->GetErrorYhigh(idx);
//       stat_data = Eff_lowpt_DATA->GetErrorYlow(idx);
//       total_syst_mc = sqrt(pow(stat_mc,2) + pow(tp,2));
//       total_syst_data = sqrt(pow(stat_data,2) + pow(tp,2));
//
//       eff_mc += Eff_lowpt_MC->GetErrorYhigh(idx);
//       eff_data -= Eff_lowpt_DATA->GetErrorYlow(idx);
//     }
//   }
//   else{
//     Eff_highpt_MC->GetPoint(idx,dummy_x,eff_mc);
//     Eff_highpt_DATA->GetPoint(idx,dummy_x,eff_data);
//
//     if(SysDirection == "up"){
//       stat_mc = Eff_highpt_MC->GetErrorYlow(idx);
//       stat_data = Eff_highpt_DATA->GetErrorYhigh(idx);
//       total_syst_mc = sqrt(pow(stat_mc,2) + pow(tp,2));
//       total_syst_data = sqrt(pow(stat_data,2) + pow(tp,2));
//
//       eff_mc -= Eff_highpt_MC->GetErrorYlow(idx);
//       eff_data += Eff_highpt_DATA->GetErrorYhigh(idx);
//     }
//     else if(SysDirection == "down"){
//       stat_mc = Eff_highpt_MC->GetErrorYhigh(idx);
//       stat_data = Eff_highpt_DATA->GetErrorYlow(idx);
//       total_syst_mc = sqrt(pow(stat_mc,2) + pow(tp,2));
//       total_syst_data = sqrt(pow(stat_data,2) + pow(tp,2));
//
//       eff_mc += Eff_highpt_MC->GetErrorYhigh(idx);
//       eff_data -= Eff_highpt_DATA->GetErrorYlow(idx);
//     }
//   }
//
//   //Scale weight by (eff_data) / (eff_mc)
//   double SF = eff_data/eff_mc;
//   event.weight *= SF;
//
// return true;
// }
//
//
//
//
//
// JetCorrectorVariable::JetCorrectorVariable(uhh2::Context & ctx, const std::vector<std::string> & JEC_files): JetCorrector(ctx, JEC_files){}
//
// bool JetCorrectorVariable::correct_collection(uhh2::Event & event, std::vector<Jet> & jets){
//
//     //apply jet corrections
//     for(auto & jet : jets){
//       correct_jet(*corrector, jet, event, jec_uncertainty, direction);
//     }
//     return true;
// };
//
//
//
//
// ElectronFakeRateWeights::ElectronFakeRateWeights(Context & ctx, const std::vector<std::string> & JEC_files, TString path_, TString SysDirection_, const string label_jets, const string label_genjets):  path(path_), SysDirection(SysDirection_){
//
//   auto dataset_type = ctx.get("dataset_type");
//   bool is_mc = dataset_type == "MC";
//   if(!is_mc){
//     cout << "Warning: ElectronFakeRateWeights will not have an effect on this non-MC sample (dataset_type = '" + dataset_type + "')" << endl;
//     return;
//   }
//
//   unique_ptr<TFile> file;
//   file.reset(new TFile(path,"READ"));
//
//   SF.reset((TGraphAsymmErrors*)file->Get("ScaleFactors"));
//   //Read the SF file and store all SFs and their pt-range
//   n_points = SF->GetN();
//   for(int i=0; i<n_points; i++){
//     double x,y;
//     SF->GetPoint(i,x,y);
//     x_low.emplace_back(x-SF->GetErrorXlow(i));
//     x_high.emplace_back(x+SF->GetErrorXhigh(i));
//   }
//
//   if(SysDirection != "nominal" && SysDirection != "up" && SysDirection != "down") throw runtime_error("In LQToTopMuModules.cxx, ElectronFakeRateWeights.process(): Invalid SysDirection specified.");
//
//   jet_corrector.reset(new JetCorrectorVariable(ctx, JEC_files));
//
//   bool jer_was_applied = ctx.get("meta_jer_applied", "") == "true";
//   if(jer_was_applied) ctx.set_metadata("jer_applied", "false", true);
//   //jet_smearer.reset(new JetSmearerVariable(ctx, label_jets, label_genjets, JERSmearing::SF_13TeV_2016));
//   jet_smearer.reset(new GenericJetResolutionSmearer(ctx, label_jets, label_genjets, true, JERSmearing::SF_13TeV_2016));
//   if(!jer_was_applied) ctx.set_metadata("jer_applied", "false", true);
//
//   jet_id = AndId<Jet>(JetPFID(JetPFID::WP_LOOSE), PtEtaCut(30.0, 2.4));
//
//   FakeRateWeightEle = ctx.get_handle<double>("FakeRateWeightEle");
//   FakeRateWeightEleUp = ctx.get_handle<double>("FakeRateWeightEleUp");
//   FakeRateWeightEleDown = ctx.get_handle<double>("FakeRateWeightEleDown");
//   h_jets = ctx.get_handle<std::vector<Jet>>(label_jets);
//
// }
//
// bool ElectronFakeRateWeights::process(Event & event){
//
//
//   vector<Jet> *jets = &event.get(h_jets);
//   jet_corrector->correct_collection(event, *jets);
//   jet_smearer->process(event);
//
//   if(event.isRealData || event.electrons->size() < 1){
//     event.set(FakeRateWeightEle,1.);
//     event.set(FakeRateWeightEleUp,1.);
//     event.set(FakeRateWeightEleDown,1.);
//     return false;
//   }
//
//   //find fake-electrons
//   vector<bool> is_fake;
//
//   //if the number of gen and reco-muons are the same, assume muons are real
//   unsigned int n_genele = 0;
//   for(const auto & gp : *event.genparticles){
//     if(fabs(gp.pdgId()) != 11) continue;
//     n_genele++;
//   }
//   //no fakes if number of gen and reco electrons is the same
//   if(n_genele == event.electrons->size()){
//     event.set(FakeRateWeightEle,1.);
//     event.set(FakeRateWeightEleUp,1.);
//     event.set(FakeRateWeightEleDown,1.);
//     return false;
//   }
//   else{
//     //if ngen and nreco are unequal, try to match muons to muons within 0.1 and to taus within 0.2
//     unsigned int n_matched_to_ele = 0, n_matched_to_tau = 0, n_matched_to_b = 0, n_matched_to_c = 0;
//     for(const auto & ele : *event.electrons){
//       bool is_matched = false;
//       for(const auto & gp : *event.genparticles){
// 	if(fabs(gp.pdgId()) == 11){
// 	  if(deltaR(gp,ele) < 0.1 && !is_matched){
// 	    is_matched = true;
// 	    n_matched_to_ele++;
// 	  }
// 	}
// 	else if(fabs(gp.pdgId()) == 15){
// 	  if(deltaR(gp,ele) < 0.2 && !is_matched){
// 	    is_matched = true;
// 	    n_matched_to_tau++;
// 	  }
// 	}
// 	else if(fabs(gp.pdgId()) == 5){
// 	  if(deltaR(gp,ele) < 0.2 && !is_matched){
// 	    is_matched = true;
// 	    n_matched_to_b++;
// 	  }
// 	}
// 	else if(fabs(gp.pdgId()) == 4){
// 	  if(deltaR(gp,ele) < 0.2 && !is_matched){
// 	    is_matched = true;
// 	    n_matched_to_c++;
// 	  }
// 	}
//       }
//       is_fake.push_back(!is_matched);
//     }
//     if(n_matched_to_tau + n_matched_to_ele + n_matched_to_b + n_matched_to_c == event.electrons->size()){
//       event.set(FakeRateWeightEle,1.);
//       event.set(FakeRateWeightEleUp,1.);
//       event.set(FakeRateWeightEleDown,1.);
//       return false;
//     }
//   }
//
//   //find jets matching fake-electrons
//   vector<bool> faking_jet;
//   for(auto & jet : *jets){
//
//     bool faking = false;
//     double dr_min = 999;
//     for(unsigned int i=0; i<event.electrons->size(); i++){
//
//       if(!is_fake[i]) continue;
//       //consider only jets to be eligible for faking the electron that fulfill the jetId used when deriving the scale factors
//       if(!jet_id(jet,event)) continue;
//
//       if(deltaR(event.electrons->at(i), jet) < 0.4) faking = true;
//       if(deltaR(event.electrons->at(i), jet) < dr_min) dr_min = deltaR(event.electrons->at(i), jet);
//     }
//     faking_jet.push_back(faking);
//   }
//
//
//   //apply SF to the eventweight for each faking jet
//   double SF_final = 1, SF_final_up = 1, SF_final_down = 1;
//   for(unsigned int i=0; i<jets->size(); i++){
//     if(!faking_jet[i]) continue;
//     double pt = event.jets->at(i).pt();
//
//     //find right bin
//     int bin = -1;
//     for(int j=0; j<n_points; j++){
//       if(pt >= x_low[j] && pt < x_high[j]) bin = j;
//     }
//
//     //possibly accout for systematics = statistical, inflate stat unc. by 50% if >800 GeV
//     double x,scale_factor, scale_factor_up, scale_factor_down, mult = 1;
//     if(pt > x_high[n_points-1]){
//       bin = n_points - 1;
//       mult = 1.5;
//     }
//     else if(bin == -1) throw runtime_error("In ElectronFakeRateWeights::process: Did not find right SF-bin for this jet.");
//
//     SF->GetPoint(bin,x,scale_factor);
//     scale_factor_up   = scale_factor + mult * SF->GetErrorYhigh(bin);
//     scale_factor_down = scale_factor - mult * SF->GetErrorYlow(bin);
//
//     if(SysDirection == "up") event.weight *= scale_factor_up;
//     else if(SysDirection == "down") event.weight *= scale_factor_down;
//     else if(SysDirection == "nominal") event.weight *= scale_factor;
//     else throw runtime_error("In ElectronFakeRateWeights::process(): SysDirection is not one of the following: ['up', 'down', 'nominal']");
//
//     SF_final      *= scale_factor;
//     SF_final_up   *= scale_factor_up;
//     SF_final_down *= scale_factor_down;
//
//   }
//
//   event.set(FakeRateWeightEle,SF_final);
//   event.set(FakeRateWeightEleUp,SF_final_up);
//   event.set(FakeRateWeightEleDown,SF_final_down);
//
//
//   return true;
// }
//
//
// MuonFakeRateWeights::MuonFakeRateWeights(Context & ctx, TString path_, TString SysDirection_):  path(path_), SysDirection(SysDirection_){
//
//   auto dataset_type = ctx.get("dataset_type");
//   bool is_mc = dataset_type == "MC";
//   if(!is_mc){
//     cout << "Warning: MuonFakeRateWeights will not have an effect on this non-MC sample (dataset_type = '" + dataset_type + "')" << endl;
//     return;
//   }
//
//   unique_ptr<TFile> file;
//   file.reset(new TFile(path,"READ"));
//   SF.reset((TGraphAsymmErrors*)file->Get("ScaleFactors"));
//
//   if(SysDirection != "nominal" && SysDirection != "up" && SysDirection != "down") throw runtime_error("In LQToTopMuModules.cxx, MuonFakeRateWeights.process(): Invalid SysDirection specified.");
//
//
//   FakeRateWeightMu = ctx.get_handle<double>("FakeRateWeightMu");
//   FakeRateWeightMuUp = ctx.get_handle<double>("FakeRateWeightMuUp");
//   FakeRateWeightMuDown = ctx.get_handle<double>("FakeRateWeightMuDown");
//
// }
//
// bool MuonFakeRateWeights::process(Event & event){
//
//   if(event.isRealData || event.muons->size() < 1){
//     event.set(FakeRateWeightMu,1.);
//     event.set(FakeRateWeightMuUp,1.);
//     event.set(FakeRateWeightMuDown,1.);
//     return false;
//   }
//
//   //find fake-muons
//   vector<bool> is_fake;
//   unsigned int n_genmu = 0;
//   for(const auto & gp : *event.genparticles){
//     if(fabs(gp.pdgId()) != 13) continue;
//     n_genmu++;
//   }
//
//   if(n_genmu == event.muons->size()){
//     event.set(FakeRateWeightMu,1.);
//     event.set(FakeRateWeightMuUp,1.);
//     event.set(FakeRateWeightMuDown,1.);
//     return false;
//   }
//   else{
//     //if ngen and nreco are unequal, try to match muons to muons within 0.1 and to taus within 0.2
//     unsigned int n_matched_to_muons = 0, n_matched_to_taus = 0, n_matched_to_b = 0, n_matched_to_c = 0;
//     for(const auto & mu : *event.muons){
//       bool is_matched = false;
//       for(const auto & gp : *event.genparticles){
// 	if(fabs(gp.pdgId()) == 13){
// 	  if(deltaR(gp,mu) < 0.1 && !is_matched){
// 	    is_matched = true;
// 	    n_matched_to_muons++;
// 	  }
// 	}
// 	else if(fabs(gp.pdgId()) == 15){
// 	  if(deltaR(gp,mu) < 0.2 && !is_matched){
// 	    is_matched = true;
// 	    n_matched_to_taus++;
// 	  }
// 	}
// 	else if(fabs(gp.pdgId()) == 5){
// 	  if(deltaR(gp,mu) < 0.2 && !is_matched){
// 	    is_matched = true;
// 	    n_matched_to_b++;
// 	  }
// 	}
// 	else if(fabs(gp.pdgId()) == 4){
// 	  if(deltaR(gp,mu) < 0.2 && !is_matched){
// 	    is_matched = true;
// 	    n_matched_to_c++;
// 	  }
// 	}
//       }
//       is_fake.push_back(!is_matched);
//     }
//     if(n_matched_to_taus + n_matched_to_muons + n_matched_to_b + n_matched_to_c == event.muons->size()){
//     event.set(FakeRateWeightMu,1.);
//     event.set(FakeRateWeightMuUp,1.);
//     event.set(FakeRateWeightMuDown,1.);
//     return false;
//     }
//   }
//
//
//   //apply SF to the eventweight for each fake muon
//   double SF_final = 1, SF_final_up = 1, SF_final_down = 1;
//   for(unsigned int i=0; i<event.muons->size(); i++){
//     if(!is_fake[i]) continue;
//     int bin = 0;
//
//     double x,scale_factor, scale_factor_up, scale_factor_down;
//     SF->GetPoint(bin,x,scale_factor);
//     scale_factor_up   = scale_factor + SF->GetErrorYhigh(bin);
//     scale_factor_down = scale_factor - SF->GetErrorYlow(bin);
//
//     if(SysDirection == "up") event.weight *= scale_factor_up;
//     else if(SysDirection == "down") event.weight *= scale_factor_down;
//     else if(SysDirection == "nominal") event.weight *= scale_factor;
//     else throw runtime_error("In MuonFakeRateWeights::process(): SysDirection is not one of the following: ['up', 'down', 'nominal']");
//
//     SF_final      *= scale_factor;
//     SF_final_up   *= scale_factor_up;
//     SF_final_down *= scale_factor_down;
//
//   }
//
//   event.set(FakeRateWeightMu,SF_final);
//   event.set(FakeRateWeightMuUp,SF_final_up);
//   event.set(FakeRateWeightMuDown,SF_final_down);
//
//
//   return true;
// }
//
//
// ZEEFinder::ZEEFinder(){}
//
// pair<int,int> ZEEFinder::search(Event & event){
//   if(event.electrons){
//     if(event.electrons->size() < 2){
//       cout << "This event does not contain >= 2 electrons, no pair can be found. Going to return (-1,-1)." << endl;
//       pair<int,int> dummy(-1,-1);
//       return dummy;
//     }
//   }
//
//   const int Nele = event.electrons->size();
//   int idx_best_ele1 = -1;
//   int idx_best_ele2 = -1;
//   double Mee_best = 99999999;
//   for(int i=0; i<Nele; i++){
//     for(int j=0; j<Nele; j++){
//       if(j>i){
// 	double Mee = (event.electrons->at(i).v4()+event.electrons->at(j).v4()).M();
// 	if(fabs(Mee - 91.2) < fabs(Mee_best - 91.2)){
// 	  Mee_best = Mee;
// 	  idx_best_ele1 = i;
// 	  idx_best_ele2 = j;
// 	}
//       }
//     }
//   }
//   pair<int,int> dummy(idx_best_ele1,idx_best_ele2);
//   return dummy;
// }
//
// pair<int,int> ZEEFinder::search(const Event & event){
//   if(event.electrons){
//     if(event.electrons->size() < 2){
//       cout << "This event does not contain >= 2 electrons, no pair can be found. Going to return (-1,-1)." << endl;
//       pair<int,int> dummy(-1,-1);
//       return dummy;
//     }
//   }
//
//   const int Nele = event.electrons->size();
//   int idx_best_ele1 = -1;
//   int idx_best_ele2 = -1;
//   double Mee_best = 99999999;
//   for(int i=0; i<Nele; i++){
//     for(int j=0; j<Nele; j++){
//       if(j>i){
// 	double Mee = (event.electrons->at(i).v4()+event.electrons->at(j).v4()).M();
// 	if(fabs(Mee - 91.2) < fabs(Mee_best - 91.2)){
// 	  Mee_best = Mee;
// 	  idx_best_ele1 = i;
// 	  idx_best_ele2 = j;
// 	}
//       }
//     }
//   }
//   pair<int,int> dummy(idx_best_ele1,idx_best_ele2);
//   return dummy;
// }
//
//
// ElectronJetOverlapCleaner::ElectronJetOverlapCleaner(){}
//
// bool ElectronJetOverlapCleaner::process(Event & event, int idx1, int idx2){
//   if(event.electrons){
//     if(event.electrons->size() < 2){
//       cout << "This event does not contain >= 2 electrons, no ElectronJetOverlapCleaning is applied." << endl;
//       return false;
//     }
//   }
//
//   if(idx1 < 0 || idx2 < 0){
//     throw runtime_error("ElectronJetOverlapCleaner was not given two valid electron indices.");
//     return false;
//   }
//
//   vector<Jet> event_jets_new;
//   int Njets = event.jets->size();
//   for(int i=0; i<Njets; i++){
//     double dR1 = deltaR(event.jets->at(i),event.electrons->at(idx1));
//     double dR2 = deltaR(event.jets->at(i),event.electrons->at(idx2));
//     if(dR1 > 0.1 && dR2 > 0.1) event_jets_new.push_back(event.jets->at(i));
//     //cout << "In ElectronJet-Cleaner: This jet has (dR1, dR2): (" << dR1 << ", " << dR2 << "), it has electron-multiplicity: " << event.jets->at(i).electronMultiplicity() << endl;
//   }
//   if(fabs(Njets - event_jets_new.size()) > 2) throw runtime_error("In LQToTopMuModules.cxx::ElectronJetOverlapCleaner::process(): When deleting jets overlapping with best electrons, more than 2 were deleted. This must not be.");
//   swap(event_jets_new, *event.jets);
//   return true;
// }
//
// DibosonScaleFactors::DibosonScaleFactors(Context & ctx, TString path_, TString SysDirectionXSec_, TString SysDirectionBTag_) : path(path_), SysDirectionXSec(SysDirectionXSec_), SysDirectionBTag(SysDirectionBTag_){
//
//   auto dataset_type = ctx.get("dataset_type");
//   bool is_mc = dataset_type == "MC";
//   if(!is_mc){
//     cout << "Warning: ElectronTriggerWeights will not have an effect on this non-MC sample (dataset_type = '" + dataset_type + "')" << endl;
//     return;
//   }
//   TString dataset_version = ctx.get("dataset_version");
//   if(!dataset_version.Contains("Diboson")){
//     is_diboson = false;
//     cout << "DibosonScaleFactors will not have an effect on this non-Diboson sample (dataset_version = '" + dataset_version + "')" << endl;
//     return;
//   }
//   else is_diboson = true;
// }
//
// bool DibosonScaleFactors::process(Event & event){
//   if(event.isRealData) return true;
//   if(!is_diboson) return true;
//
//   unique_ptr<TFile> file;
//   file.reset(new TFile(path,"READ"));
//
//   unique_ptr<TH1D> XSecSF, BTagSF;
//   XSecSF.reset((TH1D*)file->Get("Diboson_XSec_SF"));
//   BTagSF.reset((TH1D*)file->Get("Diboson_BTag_SF"));
//
//   double xsec_SF, btag1_SF, btag2_SF;
//   if(SysDirectionXSec == "nominal")   xsec_SF = XSecSF->GetBinContent(1);
//   else if(SysDirectionXSec == "up")   xsec_SF = XSecSF->GetBinContent(1) + XSecSF->GetBinError(1);
//   else if(SysDirectionXSec == "down") xsec_SF = XSecSF->GetBinContent(1) - XSecSF->GetBinError(1);
//   else throw runtime_error("In DibosonScaleFactors::process(): Invalid SysDirectionXSec specified.");
//
//   if(SysDirectionBTag == "nominal"){
//     btag1_SF = BTagSF->GetBinContent(1);
//     btag2_SF = BTagSF->GetBinContent(2);
//   }
//   else if(SysDirectionBTag == "up"){
//     btag1_SF = BTagSF->GetBinContent(1) + BTagSF->GetBinError(1);
//     btag2_SF = BTagSF->GetBinContent(2) + BTagSF->GetBinError(2);
//   }
//   else if(SysDirectionBTag == "down"){
//     btag1_SF = BTagSF->GetBinContent(1) - BTagSF->GetBinError(1);
//     btag2_SF = BTagSF->GetBinContent(2) - BTagSF->GetBinError(2);
//   }
//   else throw runtime_error("In DibosonScaleFactors::process(): Invalid SysDirectionBTag specified.");
//
//   //First apply XSec SF to all Diboson events before applying 1&2-btag SF
//   //cout << "Weight before xsec SF: " << event.weight << ", SF: " << xsec_SF << endl;
//   event.weight *= xsec_SF;
//   //cout << "Weight after xsec SF: " << event.weight << endl;
//
//   //Count number of loose b-jets in the event
//   int n_bjets = 0;
//   CSVBTag Btag_loose = CSVBTag(CSVBTag::WP_LOOSE);
//   for (unsigned int i =0; i<event.jets->size(); ++i) {
//     if(Btag_loose(event.jets->at(i),event)){
//       n_bjets++;
//     }
//   }
//
//   //cout << "Number of btags in the event: " << n_bjets << endl;
//   //cout << "SF1: " << btag1_SF << ", SF2: " << btag2_SF << endl;
//   //cout << "Weight before applying BTag SF: " << event.weight << endl;
//   //apply 1&2-btag SF
//   if(n_bjets > 2) throw runtime_error("In DibosonScaleFactors::process(): More than 2 b-jets present in the event. Scale factors have only been derived for Nbjets <= 2.");
//   else if(n_bjets == 1) event.weight *= btag1_SF;
//   else if(n_bjets == 2) event.weight *= btag2_SF;
//   //cout << "Weight after applying BTag SF: " << event.weight << endl;
//
//   return true;
// }
//
//
// MuonTrkWeights::MuonTrkWeights(Context & ctx, TString path_, TString SysDirection_): path(path_), SysDirection(SysDirection_){
//
//   auto dataset_type = ctx.get("dataset_type");
//   bool is_mc = dataset_type == "MC";
//   if(!is_mc){
//     cout << "Warning: MuonTrkWeights will not have an effect on this non-MC sample (dataset_type = '" + dataset_type + "')" << endl;
//     return;
//   }
//
//   unique_ptr<TFile> file;
//   file.reset(new TFile(path,"READ"));
//
//   Trk_SF.reset((TGraphAsymmErrors*)file->Get("ratio_eff_eta3_dr030e030_corr"));
//
//   if(SysDirection != "nominal" && SysDirection != "up" && SysDirection != "down") throw runtime_error("In LQToTopMuModules.cxx, ElectronTriggerWeights.process(): Invalid SysDirection specified.");
//
// }
//
//
// bool MuonTrkWeights::process(Event & event){
//
//   if(event.isRealData) return true;
//
//   double SF = 1.0;
//   for(const auto & mu : *event.muons){
//     double eta = mu.eta();
//     //if(fabs(eta) > 2.4) throw runtime_error("In LQToTopMuModules.cxx, MuonTrkWeights::process(): Mu-|eta| > 2.4 is not supported at the moment.");
//     if(fabs(eta) > 2.4) eta = 2.39;
//
//     //find right bin in eta
//     int idx = 0;
//     bool keep_going = true;
//     while(keep_going){
//       double x,y;
//       Trk_SF->GetPoint(idx,x,y);
//       keep_going = eta > x + Trk_SF->GetErrorXhigh(idx);
//       if(keep_going) idx++;
//     }
//
//     double eff_fnl = 1., dummy_x;
//     Trk_SF->GetPoint(idx,dummy_x,eff_fnl);
//
//     SF *= eff_fnl;
//     if(SysDirection == "up"){
//       SF *= 1.005;
//     }
//     if(SysDirection == "down"){
//       SF *= 0.995;
//     }
//
//   }
//
//   event.weight *= SF;
//
//   return true;
//
// }
