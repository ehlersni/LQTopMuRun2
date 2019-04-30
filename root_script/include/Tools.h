// #pragma once
//
// #include <cmath>
// #include <iostream>
// #include <TString.h>
// #include <TFile.h>
//
// class AnalysisTool {
//
//  public:
//
//   // Constructors, destructor
//   AnalysisTool();
//   AnalysisTool(const AnalysisTool &) = default;
//   AnalysisTool & operator = (const AnalysisTool &) = default;
//   ~AnalysisTool() = default;
//
//
//   //main functions
//   void PDFRMS(bool isSR, TString use_case);
//   void PDFRMS_Signal(bool isHT, TString samplename);
//   void PDFRMS_Background(bool isSR, TString use_case);
//   void GenerateToys();
//   void ExtrapolationFunction(bool do_ttbardy, bool do_systematics, TString SystName, TString SystDirection);
//   void NewExtrapolationFunction(bool do_ttbardy, bool do_systematics, TString SystName, TString SystDirection);
//   void NewExtrapolationFunction_lumiprojection(bool do_ttbardy, bool do_systematics, TString SystName, TString SystDirection);
//   void AlternativeExtrapolationFunction(bool do_ttbardy, bool do_systematics, TString SystName, TString SystDirection);
//   void PDFNormalizationFactor(bool is_ttbardy);
//   void ClosureTest();
//   void SubtractMCfromDATA(bool do_ttbardy, bool do_systematics, TString SystType, TString SystDirection);
//   void Prepare_StackPlots_Sideband(bool do_ttbardy);
//   void Prepare_StackPlots_Sideband_ValidationRegion(bool do_ttbardy);
//   void CompareFinalHistograms(bool do_ttbardy);
//   void ExtractSidebandTTbarMC(bool do_systematics, TString SystType, TString SystDirection);
//   void FindBiggestDeviations(bool isMC, bool is_ttbardy, TString process, bool inSR, bool isHT);
//   void FindBiggestDeviations_Split(bool isMC, bool is_ttbardy, TString doforsyst, TString process, bool inSR, bool isHT);
//   void FindBiggestDeviations_Split_AllHists(bool isMC, bool is_ttbardy, TString doforsyst, TString process, bool inSR, bool isHT);
//   void SystematicsSidebandVSMC(bool is_ttbardy, bool isHT);
//   void NormalizeLQ(bool do_systematics, TString samplename, TString SystType, TString SystDirection);
//   void CalculateCorrectMatches();
//
//  private:
//   TString base_path_SR;
//   TString base_path_CR;
//
//
// };
//
//
// class FakeRateTool {
//
//  public:
//
//   // Constructors, destructor
//   FakeRateTool();
//   FakeRateTool(const FakeRateTool &) = default;
//   FakeRateTool & operator = (const FakeRateTool &) = default;
//   ~FakeRateTool() = default;
//
//
//   //main functions
//   void CalculateDibosonSF();
//   void CalculateDibosonSF_NLO();
//   void CalculateElectronFakeRateSF();
//   void CalculateElectronFakeRateSF_DibosonNLO();
//   void CalculateElectronFakeRateSF_SystNonPrompt();
//   void CalculateElectronFakeRateSF_DibosonNLO_SystNonPrompt();
//   void CalculateMuonFakeRateSF();
//   void CalculateMuonFakeRateSF_DibosonNLO();
//   void CalculateMuonFakeRateSF_DibosonNLO_SystNonPrompt();
//
//  private:
//   TString base_path_DY_noSF;
//   TString base_path_DY_DibosonSF;
//   TString base_path_DY_DibosonSF_Muon;
//   TString base_path_Diboson;
//
//
// };
//
//
// class TriggerTool {
//
//  public:
//
//   // Constructors, destructor
//   TriggerTool();
//   TriggerTool(const TriggerTool &) = default;
//   TriggerTool & operator = (const TriggerTool &) = default;
//   ~TriggerTool() = default;
//
//
//   //main functions
//   void CalculateTriggerEfficiencies();
//   void CalculateTriggerEfficienciesData();
//
//  private:
//   TString base_path;
//
//
// };
