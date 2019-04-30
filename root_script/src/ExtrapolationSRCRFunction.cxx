#include "../include/ExtrapolationSRCRRunner.h"
#include <TString.h>
#include <TVirtualFitter.h>
#include <iostream>
#include <TStyle.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TText.h>
#include <TPaveText.h>
#include <TGaxis.h>
#include <TFitResult.h>
#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>
#include <TVectorD.h>
#include <TDecompSVD.h>
#include <sstream>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include <string>
#include <thread>
#include "../include/cosmetics.h"


using namespace std;


void ExtrapolationSRCRRunner::ExtrapolationSRCRFunction(){ //bool do_ttbardy, bool do_systematics, TString SystName, TString SystDirection

  //gStyle->SetOptFit(1100);
  gStyle->SetOptFit(0000);
  gStyle->SetOptStat(000000000);
  TVirtualFitter::SetMaxIterations(10000);

  TString path_sig  = ExtrapolationSRCRRunner::inpathSR;
  TString path_side = ExtrapolationSRCRRunner::inpathCR;
  TString path_out  = ExtrapolationSRCRRunner::outpath;

  path_sig += "/uhh2.AnalysisModuleRunner.MC.TTbarDY.root";
  path_side += "/uhh2.AnalysisModuleRunner.MC.TTbarDY.root";

  //path_sig = "/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/35867fb_DibosonNLO/HTLepInverted/uhh2.AnalysisModuleRunner.MC.TTbarDY.root";
  cout << "path_sig: " << path_sig << endl;
  cout << "path_side: " << path_side << endl;

  unique_ptr<TFile> Sideband( new TFile(path_side,"READ")); // Definitions?
  unique_ptr<TFile> Signal(new TFile(path_sig,"READ"));
  // TFile* Closure = new TFile(inpathCR+"/NOMINAL/uhh2.AnalysisModuleRunner.MC.ToyTTbarDY.root","READ");

  TH1D* h_ht_sideband_tmp, *h_ht_sideband_tmp_rebin, *h_ht_signal_tmp; //*h_ht_sideband_crosscheck_tmp, *h_ht_signal_crosscheck_tmp, *h_ht_sideband_closure_tmp;


  h_ht_signal_tmp = (TH1D*)Signal->Get("finalselection_mlqfalse/ST_rebin4");
  h_ht_sideband_tmp = (TH1D*)Sideband->Get("finalselection/ST_rebin4");
  h_ht_sideband_tmp_rebin = (TH1D*)Sideband->Get("finalselection/ST_rebin4");

  // else{
  //   h_ht_signal_tmp = (TH1D*)Signal->Get("FinalSelection/H_T_comb_NoMLQ_from350_rebin");
  //   h_ht_sideband_tmp = (TH1D*)Sideband->Get("FinalSelection/H_T_from350_rebin");
  //   h_ht_sideband_tmp_rebin = (TH1D*)Sideband->Get("FinalSelection/H_T_from350_rebin");
  //   if(!do_systematics){
  //     h_ht_signal_crosscheck_tmp = (TH1D*)Signal->Get("FinalSelection/H_T_from350_all_filled");
  //     h_ht_sideband_crosscheck_tmp = (TH1D*)Sideband->Get("FinalSelection/H_T_from350_all_filled");
  //     h_ht_sideband_closure_tmp = (TH1D*)Closure->Get("FinalSelection/H_T_from350_rebin");
  //   }
  // }

  const int n_bins = h_ht_signal_tmp->GetNbinsX();
  const double ht_start = h_ht_signal_tmp->GetBinLowEdge(1);
  const double ht_end = h_ht_signal_tmp->GetBinLowEdge(n_bins+1);
  // int n_bins_crosscheck = 1;
  vector<double> bins_crosscheck;
  // if(!do_systematics){
  //   n_bins_crosscheck = h_ht_signal_crosscheck_tmp->GetNbinsX();
  //   for(int i=1; i<n_bins_crosscheck+2; i++){
  //     bins_crosscheck.emplace_back(h_ht_signal_crosscheck_tmp->GetBinLowEdge(i));
  //   }
  // }


  //
  //   //      Define fit functions and set inital parameters
  //   //      ==============================================
  //
  const int n_param = 7;
  const double fit_min = 350.;
  const double fit_max = 3400;
  //const double fit_max = 2000;

  //ALWAYS PUT TWO "  " BEFORE AND AFTER ANY PARAMETER (  [x]  ) FOR LATER REPLACEMENT!


  TString formula = "((  [0]  * TMath::Landau(x,  [1]  ,  [2]  ,1) * (1-(1+TMath::Erf((x-  [3]  )/  [4]  ))/2))  +  ((1+TMath::Erf((x-  [3]  )/  [4]  ))/2 * (exp((  [5]  *x +  [6]  ) ))))";

  TF1 *f_LandauExp_FitSig  = new TF1("f_LandauExp_FitSig", formula, fit_min, fit_max);
  TF1 *f_LandauExp_FitSide = new TF1("f_LandauExp_FitSide", formula, fit_min, fit_max);

  f_LandauExp_FitSig->SetParameter(0,1.2e+06);
  f_LandauExp_FitSig->SetParameter(1,500);
  f_LandauExp_FitSig->SetParameter(2,90);
  f_LandauExp_FitSig->SetParameter(3,1000);
  f_LandauExp_FitSig->SetParameter(4,500);
  f_LandauExp_FitSig->SetParameter(5,-0.003);
  f_LandauExp_FitSig->SetParameter(6,7.5);


  f_LandauExp_FitSide->SetParameter(0,8.15e+05);
  f_LandauExp_FitSide->SetParameter(1,520);
  f_LandauExp_FitSide->SetParameter(2,92);
  f_LandauExp_FitSide->SetParameter(3,857);
  f_LandauExp_FitSide->SetParameter(4,369);
  f_LandauExp_FitSide->SetParameter(5,-0.0034);
  f_LandauExp_FitSide->SetParameter(6,8.5);


  //      Initial fit
  //      ===========


  // Important??
  TH1D* h_ht_sig_initial = new TH1D("h_ht_sig_initial",";S_{T} [GeV];events",n_bins,ht_start,ht_end);
  TH1D* h_ht_sig_initial_ratio = new TH1D("h_ht_sig_initial_ratio",";S_{T} [GeV];events",n_bins,ht_start,ht_end);
  TH1D* h_ht_side_initial = new TH1D("h_ht_side_initial",";S_{T}[GeV];events",n_bins,ht_start,ht_end);
  TH1D* h_ht_side_initial_ratio = new TH1D("h_ht_side_initial_ratio",";S_{T}[GeV];events",n_bins,ht_start,ht_end);
  TH1D* h_ht_sig_intermediate = new TH1D("h_ht_sig_intermediate",";S_{T} [GeV];events",n_bins,ht_start,ht_end);
  TH1D* h_ht_sig_intermediate_ratio = new TH1D("h_ht_sig_intermediate_ratio",";S_{T} [GeV];events",n_bins,ht_start,ht_end);
  TH1D* h_ht_side_intermediate = new TH1D("h_ht_side_intermediate",";S_{T}[GeV];events",n_bins,ht_start,ht_end);
  TH1D* h_ht_side_intermediate_ratio = new TH1D("h_ht_side_intermediate_ratio",";S_{T}[GeV];events",n_bins,ht_start,ht_end);
  TH1D* h_ht_sig_final = new TH1D("h_ht_sig_final",";S_{T} [GeV];events",n_bins,ht_start,ht_end);
  TH1D* h_ht_sig_final_ratio = new TH1D("h_ht_sig_final_ratio",";S_{T} [GeV];events",n_bins,ht_start,ht_end);
  TH1D* h_ht_side_final = new TH1D("h_ht_side_final",";S_{T}[GeV];events",n_bins,ht_start,ht_end);
  TH1D* h_ht_side_final_ratio = new TH1D("h_ht_side_final_ratio",";S_{T}[GeV];events",n_bins,ht_start,ht_end);
  TH1D* h_ht_sig_uncert = new TH1D("h_ht_sig_uncert",";S_{T} [GeV];events",n_bins,ht_start,ht_end);
  TH1D* h_ht_sig_uncert_ratio = new TH1D("h_ht_sig_uncert_ratio",";S_{T} [GeV];events",n_bins,ht_start,ht_end);
  TH1D* h_ht_side_uncert = new TH1D("h_ht_side_uncert",";S_{T}[GeV];events",n_bins,ht_start,ht_end);
  TH1D* h_ht_side_uncert_ratio = new TH1D("h_ht_side_uncert_ratio",";S_{T}[GeV];events",n_bins,ht_start,ht_end);
  // TH1D* h_ht_sig_crosscheck = new TH1D("h_ht_sig_crosscheck",";S_{T} [GeV];events",n_bins_crosscheck,&bins_crosscheck[0]);
  // TH1D* h_ht_sig_crosscheck_ratio = new TH1D("h_ht_sig_crosscheck_ratio",";S_{T} [GeV];events",n_bins_crosscheck,&bins_crosscheck[0]);
  // TH1D* h_ht_side_crosscheck = new TH1D("h_ht_side_crosscheck",";S_{T}[GeV];events",n_bins_crosscheck,&bins_crosscheck[0]);
  // TH1D* h_ht_side_crosscheck_ratio = new TH1D("h_ht_side_crosscheck_ratio",";S_{T}[GeV];events",n_bins_crosscheck,&bins_crosscheck[0]);
  // TH1D* h_ht_sig_closure = new TH1D("h_ht_sig_closure",";S_{T} [GeV];events",n_bins,ht_start,ht_end);
  // TH1D* h_ht_sig_closure_ratio = new TH1D("h_ht_sig_closure_ratio",";S_{T} [GeV];events",n_bins,ht_start,ht_end);
  // TH1D* h_ht_side_closure = new TH1D("h_ht_side_closure",";S_{T}[GeV];events",n_bins,ht_start,ht_end);
  // TH1D* h_ht_side_closure_ratio = new TH1D("h_ht_side_closure_ratio",";S_{T}[GeV];events",n_bins,ht_start,ht_end);
  TH1D* h_ht_sig_side_ratio = new TH1D("h_ht_sig_side_ratio",";S_{T}[GeV];events",n_bins,ht_start,ht_end);
  TH1D* h_ht_sig_side_ratio_rebin = new TH1D("h_ht_sig_side_ratio_rebin",";S_{T}[GeV];events",n_bins,ht_start,ht_end);

  for(int i=1; i<n_bins+1; i++){
    // if(h_ht_signal_tmp->GetBinCenter(i) > 2450 && h_ht_signal_tmp->GetBinCenter(i) < 2550) continue; //!!!! SystName !!! && SystName == "JER"
    h_ht_sig_initial->SetBinContent(i,h_ht_signal_tmp->GetBinContent(i));
    h_ht_sig_initial->SetBinError(i,h_ht_signal_tmp->GetBinError(i));
    h_ht_sig_initial_ratio->SetBinContent(i,h_ht_signal_tmp->GetBinContent(i));
    h_ht_sig_initial_ratio->SetBinError(i,h_ht_signal_tmp->GetBinError(i));
    h_ht_side_initial->SetBinContent(i,h_ht_sideband_tmp->GetBinContent(i));
    h_ht_side_initial->SetBinError(i,h_ht_sideband_tmp->GetBinError(i));
    h_ht_side_initial_ratio->SetBinContent(i,h_ht_sideband_tmp->GetBinContent(i));
    h_ht_side_initial_ratio->SetBinError(i,h_ht_sideband_tmp->GetBinError(i));
    h_ht_sig_intermediate->SetBinContent(i,h_ht_signal_tmp->GetBinContent(i));
    h_ht_sig_intermediate->SetBinError(i,h_ht_signal_tmp->GetBinError(i));
    h_ht_sig_intermediate_ratio->SetBinContent(i,h_ht_signal_tmp->GetBinContent(i));
    h_ht_sig_intermediate_ratio->SetBinError(i,h_ht_signal_tmp->GetBinError(i));
    h_ht_side_intermediate->SetBinContent(i,h_ht_sideband_tmp->GetBinContent(i));
    h_ht_side_intermediate->SetBinError(i,h_ht_sideband_tmp->GetBinError(i));
    h_ht_side_intermediate_ratio->SetBinContent(i,h_ht_sideband_tmp->GetBinContent(i));
    h_ht_side_intermediate_ratio->SetBinError(i,h_ht_sideband_tmp->GetBinError(i));
    h_ht_sig_final->SetBinContent(i,h_ht_signal_tmp->GetBinContent(i));
    h_ht_sig_final->SetBinError(i,h_ht_signal_tmp->GetBinError(i));
    h_ht_sig_final_ratio->SetBinContent(i,h_ht_signal_tmp->GetBinContent(i));
    h_ht_sig_final_ratio->SetBinError(i,h_ht_signal_tmp->GetBinError(i));
    h_ht_side_final->SetBinContent(i,h_ht_sideband_tmp->GetBinContent(i));
    h_ht_side_final->SetBinError(i,h_ht_sideband_tmp->GetBinError(i));
    h_ht_side_final_ratio->SetBinContent(i,h_ht_sideband_tmp->GetBinContent(i));
    h_ht_side_final_ratio->SetBinError(i,h_ht_sideband_tmp->GetBinError(i));
    h_ht_sig_uncert->SetBinContent(i,h_ht_signal_tmp->GetBinContent(i));
    h_ht_sig_uncert->SetBinError(i,h_ht_signal_tmp->GetBinError(i));
    h_ht_sig_uncert_ratio->SetBinContent(i,h_ht_signal_tmp->GetBinContent(i));
    h_ht_sig_uncert_ratio->SetBinError(i,h_ht_signal_tmp->GetBinError(i));
    h_ht_side_uncert->SetBinContent(i,h_ht_sideband_tmp->GetBinContent(i));
    h_ht_side_uncert->SetBinError(i,h_ht_sideband_tmp->GetBinError(i));
    h_ht_side_uncert_ratio->SetBinContent(i,h_ht_sideband_tmp->GetBinContent(i));
    h_ht_side_uncert_ratio->SetBinError(i,h_ht_sideband_tmp->GetBinError(i));
    h_ht_sig_side_ratio->SetBinContent(i,h_ht_signal_tmp->GetBinContent(i));
    h_ht_sig_side_ratio->SetBinError(i,h_ht_signal_tmp->GetBinError(i));
    h_ht_sig_side_ratio_rebin->SetBinContent(i,h_ht_signal_tmp->GetBinContent(i));
    h_ht_sig_side_ratio_rebin->SetBinError(i,h_ht_signal_tmp->GetBinError(i));


    // if(!do_systematics){
    //   h_ht_sig_closure->SetBinContent(i,h_ht_signal_tmp->GetBinContent(i));
    //   h_ht_sig_closure->SetBinError(i,h_ht_signal_tmp->GetBinError(i));
    //   h_ht_sig_closure_ratio->SetBinContent(i,h_ht_signal_tmp->GetBinContent(i));
    //   h_ht_sig_closure_ratio->SetBinError(i,h_ht_signal_tmp->GetBinError(i));
    //   h_ht_side_closure->SetBinContent(i,h_ht_sideband_closure_tmp->GetBinContent(i));
    //   h_ht_side_closure->SetBinError(i,h_ht_sideband_closure_tmp->GetBinError(i));
    //   h_ht_side_closure_ratio->SetBinContent(i,h_ht_sideband_closure_tmp->GetBinContent(i));
    //   h_ht_side_closure_ratio->SetBinError(i,h_ht_sideband_closure_tmp->GetBinError(i));
    // }
  }

  // if(!do_systematics){
  //   for(int i=1; i<n_bins_crosscheck+1; i++){
  //     h_ht_sig_crosscheck->SetBinContent(i,h_ht_signal_crosscheck_tmp->GetBinContent(i));
  //     h_ht_sig_crosscheck->SetBinError(i,h_ht_signal_crosscheck_tmp->GetBinError(i));
  //     h_ht_sig_crosscheck_ratio->SetBinContent(i,h_ht_signal_crosscheck_tmp->GetBinContent(i));
  //     h_ht_sig_crosscheck_ratio->SetBinError(i,h_ht_signal_crosscheck_tmp->GetBinError(i));
  //     h_ht_side_crosscheck->SetBinContent(i,h_ht_sideband_crosscheck_tmp->GetBinContent(i));
  //     h_ht_side_crosscheck->SetBinError(i,h_ht_sideband_crosscheck_tmp->GetBinError(i));
  //     h_ht_side_crosscheck_ratio->SetBinContent(i,h_ht_sideband_crosscheck_tmp->GetBinContent(i));
  //     h_ht_side_crosscheck_ratio->SetBinError(i,h_ht_sideband_crosscheck_tmp->GetBinError(i));
  //   }
  // }
  h_ht_sig_side_ratio->Divide(h_ht_sideband_tmp);

  h_ht_sig_side_ratio_rebin->Rebin(2);
  h_ht_sideband_tmp_rebin->Rebin(2);
  h_ht_sig_side_ratio_rebin->Divide(h_ht_sideband_tmp_rebin);


  TCanvas* c_sig_initial = new TCanvas("c_sig_initial", "Signal region initial fit", 400, 400);
  TPad* p_top = SetupRatioPadTop();
  TPad* p_bot = SetupRatioPad();
  p_top->Draw();
  p_bot->Draw();
  p_top->cd();
  HistCosmetics(h_ht_sig_initial,false);
  HistCosmetics(h_ht_sig_initial_ratio,true);
  TFitResultPtr fit_sig_initial = h_ht_sig_initial->Fit(f_LandauExp_FitSig,"SR");
  fit_sig_initial->Print("V");
  h_ht_sig_initial->SetMinimum(0.003);
  h_ht_sig_initial->GetXaxis()->SetRangeUser(0,fit_max);
  //p_top->SetLogy();


  //calculate ratio
  for(int i=1; i<n_bins+1; i++){
    double nom = h_ht_sig_initial_ratio->GetBinContent(i);
    double err =  h_ht_sig_initial_ratio->GetBinError(i);
    double denom = f_LandauExp_FitSig->Eval(h_ht_sig_initial_ratio->GetBinCenter(i));
    double ratio = nom/denom;
    h_ht_sig_initial_ratio->SetBinContent(i,ratio);
    h_ht_sig_initial_ratio->SetBinError(i,err/denom);
  }
  p_bot->cd();
  h_ht_sig_initial_ratio->GetYaxis()->SetRangeUser(0.,2.4);
  h_ht_sig_initial_ratio->GetXaxis()->SetRangeUser(0.,fit_max);
  h_ht_sig_initial_ratio->GetYaxis()->SetTitle("MC / fit");
  h_ht_sig_initial_ratio->Draw();

  TLine* l_unity = new TLine(fit_min, 1, fit_max, 1);
  l_unity->SetLineColor(kRed);
  l_unity->SetLineWidth(2);
  l_unity->Draw();
  c_sig_initial->SaveAs(path_out + "/Alpha_Output/FitSigOld.eps"); // + SystName + "_" + SystDirection
  delete c_sig_initial;
  //delete p_top;
  //delete p_bot;

  TCanvas* c_side_initial = new TCanvas("c_side_initial", "Control region initial fit", 400, 400);
  p_top = SetupRatioPadTop();
  p_bot = SetupRatioPad();
  p_top->Draw();
  p_bot->Draw();
  p_top->cd();
  HistCosmetics(h_ht_side_initial,false);
  HistCosmetics(h_ht_side_initial_ratio,true);
  TFitResultPtr fit_side_initial = h_ht_side_initial->Fit(f_LandauExp_FitSide,"SR");
  fit_side_initial->Print("V");
  h_ht_side_initial->SetMinimum(0.003);
  h_ht_side_initial->GetXaxis()->SetRangeUser(0,fit_max);
  //p_top->SetLogy();


  //calculate ratio
  for(int i=1; i<n_bins+1; i++){
    double nom = h_ht_side_initial_ratio->GetBinContent(i);
    double err =  h_ht_side_initial_ratio->GetBinError(i);
    double denom = f_LandauExp_FitSide->Eval(h_ht_side_initial_ratio->GetBinCenter(i));
    double ratio = nom/denom;
    h_ht_side_initial_ratio->SetBinContent(i,ratio);
    h_ht_side_initial_ratio->SetBinError(i,err/denom);
  }
  p_bot->cd();
  h_ht_side_initial_ratio->GetYaxis()->SetRangeUser(0.,2.4);
  h_ht_side_initial_ratio->GetXaxis()->SetRangeUser(0.,fit_max);
  h_ht_side_initial_ratio->GetYaxis()->SetTitle("MC / fit");
  h_ht_side_initial_ratio->Draw();
  l_unity->Draw();
  c_side_initial->SaveAs(path_out + "/Alpha_Output/FitSideOld.eps"); // + SystName + "_" + SystDirection
  delete c_side_initial;
  //delete p_top;
  //delete p_bot;




  //
  //   //      Access to the fit result
  //   //      ========================
  //
    TMatrixDSym cov_sig_intermediate = fit_sig_initial->GetCovarianceMatrix();
    TMatrixDSym cov_side_intermediate = fit_side_initial->GetCovarianceMatrix();

    const TMatrixDSymEigen eigen_sig_intermediate(cov_sig_intermediate);
    const TMatrixDSymEigen eigen_side_intermediate(cov_side_intermediate);

    //Eigenvalues on diagonal
    const TVectorD V_Diag_Sig_intermediate = eigen_sig_intermediate.GetEigenValues();
    const TVectorD V_Diag_Side_intermediate = eigen_side_intermediate.GetEigenValues();

    //Normalized eigenvectors --> wanted rotation matrix
    const TMatrixD M_Orth_Sig_intermediate = eigen_sig_intermediate.GetEigenVectors();
    const TMatrixD M_Orth_Side_intermediate = eigen_side_intermediate.GetEigenVectors();

    TMatrixD M_Diag_Sig_intermediate(n_param,n_param);
    TMatrixD M_Diag_Side_intermediate(n_param,n_param);
    M_Diag_Sig_intermediate.UnitMatrix();
    M_Diag_Side_intermediate.UnitMatrix();
    for(int i=0; i<n_param; i++){
      double tmp_sig = V_Diag_Sig_intermediate[i];
      TMatrixDRow(M_Diag_Sig_intermediate,i)(i) = tmp_sig;
      double tmp_side = V_Diag_Side_intermediate[i];
      TMatrixDRow(M_Diag_Side_intermediate,i)(i) = tmp_side;
    }

    TMatrixD RT_Sig_intermediate(M_Orth_Sig_intermediate);
    RT_Sig_intermediate.T();
    TMatrixD RT_Side_intermediate(M_Orth_Side_intermediate);
    RT_Side_intermediate.T();

    vector<double> params_side_initial, params_sig_initial;
    for(int i=0; i<n_param; i++){
      params_side_initial.emplace_back(fit_side_initial->Parameter(i));
      params_sig_initial.emplace_back(fit_sig_initial->Parameter(i));
    }

    const TVectorD V_param_sig_initial(n_param,&params_sig_initial[0]);
    const TVectorD V_param_side_initial(n_param,&params_side_initial[0]);


    //calculate parameters in new basis
    TDecompSVD svd_sig_intermediate(M_Orth_Sig_intermediate);
    TDecompSVD svd_side_intermediate(M_Orth_Side_intermediate);
    Bool_t ok_sig;
    Bool_t ok_side;
    const TVectorD V_param_sig_new_intermediate = svd_sig_intermediate.Solve(V_param_sig_initial,ok_sig);
    const TVectorD V_param_side_new_intermediate = svd_side_intermediate.Solve(V_param_side_initial,ok_side);

    //Transform back new parameters into old ones: express old ones using new ones (for error propagation)
    TDecompSVD svd_sig_intermediate2(RT_Sig_intermediate);
    TDecompSVD svd_side_intermediate2(RT_Side_intermediate);
    Bool_t ok_sig2;
    Bool_t ok_side2;
    const TVectorD V_param_sig_initial_by_new_intermediate = svd_sig_intermediate2.Solve(V_param_sig_new_intermediate,ok_sig2);
    const TVectorD V_param_side_initial_by_new_intermediate = svd_side_intermediate2.Solve(V_param_side_new_intermediate,ok_side2);

    //parametrize old parameters by using new ones with M_Orth_Sig
    TMatrixD M_Dummy_Sig_intermediate(M_Orth_Sig_intermediate);
    TMatrixD M_Dummy_Side_intermediate(M_Orth_Side_intermediate);

    //Fill vector of vectors: 1 vector for each old parameter, containing the values of the new parameters expressing the old one.
    vector<vector<double>> param_side_transformed_intermediate, param_sig_transformed_intermediate;
    for(int i=0; i<n_param; i++){
      vector<double> single_param_sig, single_param_side;
      for(int j=0; j<n_param; j++){
        single_param_sig.emplace_back(TMatrixDRow(M_Dummy_Sig_intermediate,i)[j]);
        single_param_side.emplace_back(TMatrixDRow(M_Dummy_Side_intermediate,i)[j]);
      }
      param_sig_transformed_intermediate.emplace_back(single_param_sig);
      param_side_transformed_intermediate.emplace_back(single_param_side);
      for(unsigned int j=0; j<single_param_sig.size(); j++){
        single_param_sig.pop_back();
        single_param_side.pop_back();
      }
    }


  //
  //   //      Transform the old formula into the new one
  //   //      ==========================================
  //
    TString new_formula_sig_intermediate = formula;
    TString new_formula_side_intermediate = formula;
    for(int i=0; i<n_param; i++){
      TString idx = "  [";
      idx += i;
      idx += "]  ";
      TString new_expr_sig = "( [0] *(";
      new_expr_sig += param_sig_transformed_intermediate[i][0];
      new_expr_sig += ")";
      TString new_expr_side = "( [0] *(";
      new_expr_side += param_side_transformed_intermediate[i][0];
      new_expr_side += ")";
      for(int j=1; j<n_param; j++){
        TString newidx = "";
        newidx += j;
        new_expr_sig += " + ["+newidx+"] *(";
        new_expr_sig += param_sig_transformed_intermediate[i][j];
        new_expr_sig += ")";
        new_expr_side += " + ["+newidx+"] *(";
        new_expr_side += param_side_transformed_intermediate[i][j];
        new_expr_side += ")";
      }
      new_expr_sig += ")";
      new_expr_side += ")";
      new_formula_sig_intermediate.ReplaceAll(idx,new_expr_sig);
      new_formula_side_intermediate.ReplaceAll(idx,new_expr_side);
    }


    TF1* f_LandauExp_NewParam_Sig_intermediate = new TF1("f_LandauExp_NewParam_Sig_intermediate",new_formula_sig_intermediate,fit_min,fit_max);
    TF1* f_LandauExp_NewParam_Side_intermediate = new TF1("f_LandauExp_NewParam_Side_intermediate",new_formula_side_intermediate,fit_min,fit_max);
    for(int i=0; i<n_param; i++){
      f_LandauExp_NewParam_Sig_intermediate->SetParameter(i,V_param_sig_new_intermediate[i]);
      f_LandauExp_NewParam_Side_intermediate->SetParameter(i,V_param_side_new_intermediate[i]);
    }

    TFitResultPtr fit_result_sig_intermediate = h_ht_sig_intermediate->Fit(f_LandauExp_NewParam_Sig_intermediate,"SR");
    TFitResultPtr fit_result_side_intermediate = h_ht_side_intermediate->Fit(f_LandauExp_NewParam_Side_intermediate,"SR");

    fit_result_sig_intermediate->Print("V");
    fit_result_side_intermediate->Print("V");
  //
  //
  //
  //
  //   //      Second diagonalization for numeric reasons
  //   //      ==========================================
  //
  //
    TMatrixDSym cov_sig_final = fit_result_sig_intermediate->GetCovarianceMatrix();
    TMatrixDSym cov_side_final = fit_result_side_intermediate->GetCovarianceMatrix();

    const TMatrixDSymEigen eigen_sig_final(cov_sig_final);
    const TMatrixDSymEigen eigen_side_final(cov_side_final);

    //Eigenvalues on diagonal
    const TVectorD V_Diag_Sig_final = eigen_sig_final.GetEigenValues();
    const TVectorD V_Diag_Side_final = eigen_side_final.GetEigenValues();

    //Normalized eigenvectors --> wanted rotation matrix
    const TMatrixD M_Orth_Sig_final = eigen_sig_final.GetEigenVectors();
    const TMatrixD M_Orth_Side_final = eigen_side_final.GetEigenVectors();

    TMatrixD M_Diag_Sig_final(n_param,n_param);
    TMatrixD M_Diag_Side_final(n_param,n_param);
    M_Diag_Sig_final.UnitMatrix();
    M_Diag_Side_final.UnitMatrix();
    for(int i=0; i<n_param; i++){
      double tmp_sig = V_Diag_Sig_final[i];
      TMatrixDRow(M_Diag_Sig_final,i)(i) = tmp_sig;
      double tmp_side = V_Diag_Side_final[i];
      TMatrixDRow(M_Diag_Side_final,i)(i) = tmp_side;
    }

    TMatrixD RT_Sig_final(M_Orth_Sig_final);
    RT_Sig_final.T();
    TMatrixD RT_Side_final(M_Orth_Side_final);
    RT_Side_final.T();

    vector<double> params_side_initial_final, params_sig_initial_final;
    for(int i=0; i<n_param; i++){
      params_side_initial_final.emplace_back(fit_result_side_intermediate->Parameter(i));
      params_sig_initial_final.emplace_back(fit_result_sig_intermediate->Parameter(i));
    }


    const TVectorD V_param_sig_initial_final(n_param,&params_sig_initial_final[0]);
    const TVectorD V_param_side_initial_final(n_param,&params_side_initial_final[0]);

    //calculate parameters in new basis
    TDecompSVD svd_sig3(M_Orth_Sig_final);
    TDecompSVD svd_side3(M_Orth_Side_final);
    Bool_t ok_sig3;
    Bool_t ok_side3;
    const TVectorD V_param_sig_new_final = svd_sig3.Solve(V_param_sig_initial_final,ok_sig3);
    const TVectorD V_param_side_new_final = svd_side3.Solve(V_param_side_initial_final,ok_side3);

    //Transform back new parameters into old ones: express old ones using new ones (for error propagation)
    TDecompSVD svd_sig4(RT_Sig_intermediate);
    TDecompSVD svd_side4(RT_Side_intermediate);
    Bool_t ok_sig4;
    Bool_t ok_side4;
    const TVectorD V_param_sig_initial_by_new_final = svd_sig4.Solve(V_param_sig_new_final,ok_sig4);
    const TVectorD V_param_side_initial_by_new_final = svd_side4.Solve(V_param_side_new_final,ok_side4);

    //parametrize old parameters by using new ones with M_Orth_Sig
    TMatrixD M_Dummy_Sig_final(M_Orth_Sig_final);
    TMatrixD M_Dummy_Side_final(M_Orth_Side_final);

    //Fill vector of vectors: 1 vector for each old parameter, containing the values of the new parameters expressing the old one.
    vector<vector<double>> param_side_transformed_final, param_sig_transformed_final;
    for(int i=0; i<n_param; i++){
      vector<double> single_param_sig, single_param_side;
      for(int j=0; j<n_param; j++){
        single_param_sig.emplace_back(TMatrixDRow(M_Dummy_Sig_final,i)[j]);
        single_param_side.emplace_back(TMatrixDRow(M_Dummy_Side_final,i)[j]);
      }
      param_sig_transformed_final.emplace_back(single_param_sig);
      param_side_transformed_final.emplace_back(single_param_side);
      for(unsigned int j=0; j<single_param_sig.size(); j++){
        single_param_sig.pop_back();
        single_param_side.pop_back();
      }
    }
  //
  //
  //   //      Transform the old formula into the new one
  //   //      ==========================================
  //
    TString new_formula_sig_final = new_formula_sig_intermediate;
    TString new_formula_side_final = new_formula_side_intermediate;
    for(int i=0; i<n_param; i++){
      TString idx = " [";
      idx += i;
      idx += "] ";
      TString new_expr_sig = "([0]*(";
      new_expr_sig += param_sig_transformed_final[i][0];
      new_expr_sig += ")";
      TString new_expr_side = "([0]*(";
      new_expr_side += param_side_transformed_final[i][0];
      new_expr_side += ")";
      for(int j=1; j<n_param; j++){
        TString newidx = "";
        newidx += j;
        new_expr_sig += " + ["+newidx+"]*(";
        new_expr_sig += param_sig_transformed_final[i][j];
        new_expr_sig += ")";
        new_expr_side += " + ["+newidx+"]*(";
        new_expr_side += param_side_transformed_final[i][j];
        new_expr_side += ")";
      }
      new_expr_sig += ")";
      new_expr_side += ")";
      new_formula_sig_final.ReplaceAll(idx,new_expr_sig);
      new_formula_side_final.ReplaceAll(idx,new_expr_side);
    }


    TF1* f_LandauExp_NewParam_Sig_final = new TF1("f_LandauExp_NewParam_Sig",new_formula_sig_final,fit_min,fit_max);
    TF1* f_LandauExp_NewParam_Side_final = new TF1("f_LandauExp_NewParam_Side",new_formula_side_final,fit_min,fit_max);
    for(int i=0; i<n_param; i++){
      f_LandauExp_NewParam_Sig_final->SetParameter(i,V_param_sig_new_final[i]);
      f_LandauExp_NewParam_Side_final->SetParameter(i,V_param_side_new_final[i]);
    }

    TFitResultPtr fit_result_sig_final = h_ht_sig_final->Fit(f_LandauExp_NewParam_Sig_final,"SR");
    TFitResultPtr fit_result_side_final = h_ht_side_final->Fit(f_LandauExp_NewParam_Side_final,"SR");

    fit_result_sig_final->Print("V");
    fit_result_side_final->Print("V");

  //
  //
  //
  //
  //
  //   //      Worst thing done, now calculate errors of parameters
  //   //      ====================================================
  //
  //
    //Signal first
    TF1* f_Signal_errors = new TF1("f_Signal_errors", new_formula_sig_final, fit_min, fit_max);

    for(int i=0; i<n_param; i++){
      f_Signal_errors->FixParameter(i,fit_result_sig_final->Parameter(i));
    }
    h_ht_sig_final->Fit(f_Signal_errors,"QR");
    double chi2_best_sig = f_Signal_errors->GetChisquare();

    vector<double>  pars_sig, pars_up_sig, pars_dn_sig, err_pars_up_sig, err_pars_dn_sig;
    vector<bool>    pars_up_set_sig, pars_dn_set_sig;
    vector<TGraph*> g_chi2_pars_sig;
    for(int i=0; i<n_param; i++){
      pars_sig.emplace_back(f_Signal_errors->GetParameter(i));
      pars_up_set_sig.emplace_back(false);
      pars_dn_set_sig.emplace_back(false);
      g_chi2_pars_sig.emplace_back(new TGraph(50000));
    }

    int for_percent = 0;
    cout << endl << endl;
    for(int j=0; j<n_param; j++){

      //reset parameters
      for(int k=0; k<n_param; k++) f_Signal_errors->FixParameter(k,fit_result_sig_final->Parameter(k));

      double chi2_diff = 0;
      for(int i=0; i<50000; i++){
        for_percent++;
        double varied_par = -20*pars_sig[j] + pars_sig[j]/1000 * i;//from -20p to +30p in steps of p/1000
        f_Signal_errors->FixParameter(j,varied_par);
        h_ht_sig_final->Fit(f_Signal_errors,"QR");

        double chi2_tmp = f_Signal_errors->GetChisquare();
        chi2_diff = chi2_tmp - chi2_best_sig;

        int divisor = 50000*n_param/100;
        if(i%divisor == 0) cout << "\r"  << for_percent/divisor  << "%, " <<  "Parameter " << j << ", Fit No. " << i  << ", current chi2: " << chi2_tmp << " at current parametervalue: " << varied_par << "           "  << flush;

        if(!pars_dn_set_sig[j] && chi2_diff < 1) {
  	pars_dn_sig.emplace_back(-20*pars_sig[j] + pars_sig[j]/1000 * (i-1));
  	pars_dn_set_sig[j] = true;
        }

        if(pars_dn_set_sig[j] && !pars_up_set_sig[j] && chi2_diff >= 1){
  	pars_up_sig.emplace_back(varied_par);
  	pars_up_set_sig[j] = true;
        }

        //draw chi2 as function of p0
        g_chi2_pars_sig[j]->SetPoint(i,varied_par, chi2_tmp - chi2_best_sig);
      }
      if(pars_dn_sig[j] > pars_up_sig[j]) swap(pars_dn_sig[j],pars_up_sig[j]);

      err_pars_up_sig.emplace_back(fabs(pars_sig[j]-pars_up_sig[j]));
      err_pars_dn_sig.emplace_back(fabs(pars_sig[j]-pars_dn_sig[j]));
      cout << "Parameter " << j << ": " << pars_sig[j] << " + " << pars_up_sig[j] << " - " << pars_dn_sig[j] << endl;
    }
    cout << endl;






    vector<TLine*> l_pars_dn_sig, l_pars_up_sig, l_pars_central_sig;
    for(int i=0; i<n_param; i++){
      l_pars_up_sig.emplace_back(     new TLine(pars_up_sig[i], 10, pars_up_sig[i], -1));
      l_pars_dn_sig.emplace_back(     new TLine(pars_dn_sig[i], 10, pars_dn_sig[i], -1));
      l_pars_central_sig.emplace_back(new TLine(pars_sig[i], 0, pars_sig[i], -1));
      l_pars_up_sig[i]->SetLineWidth(3);
      l_pars_up_sig[i]->SetLineColor(kRed);
      l_pars_dn_sig[i]->SetLineWidth(3);
      l_pars_dn_sig[i]->SetLineColor(kRed);
      l_pars_central_sig[i]->SetLineWidth(3);
      l_pars_central_sig[i]->SetLineColor(kBlue);
    }


    for(int i=0; i<n_param; i++){
      TString name = "p";
      name += i;
      name += "\'";

      g_chi2_pars_sig[i]->SetMinimum(-1);
      g_chi2_pars_sig[i]->SetMaximum(10);
      g_chi2_pars_sig[i]->SetTitle("Signal region");
      g_chi2_pars_sig[i]->GetXaxis()->SetTitleSize(0.05);
      g_chi2_pars_sig[i]->GetXaxis()->SetLabelSize(0.045);
      g_chi2_pars_sig[i]->GetXaxis()->SetTitleOffset(0.9);
      g_chi2_pars_sig[i]->GetXaxis()->SetTitle(name);
      g_chi2_pars_sig[i]->GetYaxis()->SetTitle("#chi^{2} - #chi^{2}_{min}");
      g_chi2_pars_sig[i]->GetYaxis()->SetTitleSize(0.06);
      g_chi2_pars_sig[i]->GetYaxis()->SetLabelSize(0.045);
      g_chi2_pars_sig[i]->GetYaxis()->SetTitleOffset(0.8);
      g_chi2_pars_sig[i]->GetXaxis()->SetLimits(pars_dn_sig[i]-2*err_pars_dn_sig[i],pars_up_sig[i]+2*err_pars_up_sig[i]);
    }

    for(int i=0; i<n_param; i++){
      TCanvas* canv = new TCanvas();
      g_chi2_pars_sig[i]->Draw();
      l_pars_up_sig[i]->Draw("SAME");
      l_pars_dn_sig[i]->Draw("SAME");
      l_pars_central_sig[i]->Draw("SAME");
      TString number = "";
      number += i;
      canv->Print(path_out + "Alpha_Output/Chi2_Sig_p"+number+".eps");
      delete canv;
    }



    //Sideband second
    TF1* f_Sideband_errors = new TF1("f_Sideband_errors", new_formula_side_final, fit_min, fit_max);

    for(int i=0; i<n_param; i++){
      f_Sideband_errors->FixParameter(i,fit_result_side_final->Parameter(i));
    }
    h_ht_side_final->Fit(f_Sideband_errors,"QR");
    double chi2_best_side = f_Sideband_errors->GetChisquare();

    vector<double>  pars_side, pars_up_side, pars_dn_side, err_pars_up_side, err_pars_dn_side;
    vector<bool>    pars_up_set_side, pars_dn_set_side;
    vector<TGraph*> g_chi2_pars_side;
    for(int i=0; i<n_param; i++){
      pars_side.emplace_back(f_Sideband_errors->GetParameter(i));
      pars_up_set_side.emplace_back(false);
      pars_dn_set_side.emplace_back(false);
      g_chi2_pars_side.emplace_back(new TGraph(50000));
    }

    for_percent = 0;
    cout << endl << endl;
    for(int j=0; j<n_param; j++){

      //reset parameters
      for(int k=0; k<n_param; k++) f_Sideband_errors->FixParameter(k,fit_result_side_final->Parameter(k));

      double chi2_diff = 0;
      for(int i=0; i<50000; i++){
        for_percent++;
        double varied_par = -20*pars_side[j] + pars_side[j]/1000 * i;//from -20p to +30p in steps of p/1000
        f_Sideband_errors->FixParameter(j,varied_par);
        h_ht_side_final->Fit(f_Sideband_errors,"QR");

        double chi2_tmp = f_Sideband_errors->GetChisquare();
        chi2_diff = chi2_tmp - chi2_best_side;

        int divisor = 50000*n_param/100;
        if(i%divisor == 0) cout << "\r"  << for_percent/divisor  << "%, " <<  "Parameter " << j << ", Fit No. " << i  << ", current chi2: " << chi2_tmp << " at current parametervalue: " << varied_par << "           "  << flush;

        if(!pars_dn_set_side[j] && chi2_diff < 1) {
  	pars_dn_side.emplace_back(-20*pars_side[j] + pars_side[j]/1000 * (i-1));
  	pars_dn_set_side[j] = true;
        }

        if(pars_dn_set_side[j] && !pars_up_set_side[j] && chi2_diff >= 1){
  	pars_up_side.emplace_back(varied_par);
  	pars_up_set_side[j] = true;
        }

        //draw chi2 as function of p0
        g_chi2_pars_side[j]->SetPoint(i,varied_par, chi2_tmp - chi2_best_side);
      }
      if(pars_dn_side[j] > pars_up_side[j]) swap(pars_dn_side[j],pars_up_side[j]);

      err_pars_up_side.emplace_back(fabs(pars_side[j]-pars_up_side[j]));
      err_pars_dn_side.emplace_back(fabs(pars_side[j]-pars_dn_side[j]));
      cout << "Parameter " << j << ": " << pars_side[j] << " + " << pars_up_side[j] << " - " << pars_dn_side[j] << endl;
    }
    cout << endl;


    vector<TLine*> l_pars_dn_side, l_pars_up_side, l_pars_central_side;
    for(int i=0; i<n_param; i++){
      l_pars_up_side.emplace_back(     new TLine(pars_up_side[i], 10, pars_up_side[i], -1));
      l_pars_dn_side.emplace_back(     new TLine(pars_dn_side[i], 10, pars_dn_side[i], -1));
      l_pars_central_side.emplace_back(new TLine(pars_side[i], 0, pars_side[i], -1));
      l_pars_up_side[i]->SetLineWidth(3);
      l_pars_up_side[i]->SetLineColor(kRed);
      l_pars_dn_side[i]->SetLineWidth(3);
      l_pars_dn_side[i]->SetLineColor(kRed);
      l_pars_central_side[i]->SetLineWidth(3);
      l_pars_central_side[i]->SetLineColor(kBlue);
    }


    for(int i=0; i<n_param; i++){
      TString name = "p";
      name += i;
      name += "\'";

      g_chi2_pars_side[i]->SetMinimum(-1);
      g_chi2_pars_side[i]->SetMaximum(10);
      g_chi2_pars_side[i]->SetTitle("Control region");
      g_chi2_pars_side[i]->GetXaxis()->SetTitleSize(0.05);
      g_chi2_pars_side[i]->GetXaxis()->SetLabelSize(0.045);
      g_chi2_pars_side[i]->GetXaxis()->SetTitleOffset(0.9);
      g_chi2_pars_side[i]->GetXaxis()->SetTitle(name);
      g_chi2_pars_side[i]->GetYaxis()->SetTitle("#chi^{2} - #chi^{2}_{min}");
      g_chi2_pars_side[i]->GetYaxis()->SetTitleSize(0.06);
      g_chi2_pars_side[i]->GetYaxis()->SetLabelSize(0.045);
      g_chi2_pars_side[i]->GetYaxis()->SetTitleOffset(0.8);
      g_chi2_pars_side[i]->GetXaxis()->SetLimits(pars_dn_side[i]-2*err_pars_dn_side[i],pars_up_side[i]+2*err_pars_up_side[i]);
    }

    for(int i=0; i<n_param; i++){
      TCanvas* canv = new TCanvas();
      g_chi2_pars_side[i]->Draw();
      l_pars_up_side[i]->Draw("SAME");
      l_pars_dn_side[i]->Draw("SAME");
      l_pars_central_side[i]->Draw("SAME");
      TString number = "";
      number += i;
      canv->Print(path_out + "Alpha_Output/Chi2_Side_p"+number+".eps");
      delete canv;
    }


  //
  //
  //
  //   //      Add all variations in quadrature
  //   //      ================================
  //
    const int final_range = 4000;

    vector<TF1*> f_up_sig, f_dn_sig, f_up_side, f_dn_side;
    for(int i=0; i<n_param; i++){
      TString number = "";
      number += i;
      TString name_up_sig = "f_up_sig_" + number;
      TString name_dn_sig = "f_dn_sig_" + number;
      TString name_up_side = "f_up_side_" + number;
      TString name_dn_side = "f_dn_side_" + number;
      f_up_sig.emplace_back(new TF1(name_up_sig, new_formula_sig_final, 350, final_range));
      f_dn_sig.emplace_back(new TF1(name_dn_sig, new_formula_sig_final, 350, final_range));
      f_up_side.emplace_back(new TF1(name_up_side, new_formula_side_final, 350, final_range));
      f_dn_side.emplace_back(new TF1(name_dn_side, new_formula_side_final, 350, final_range));
    }

    for(int i=0; i<n_param; i++){
      for(int j=0; j<n_param; j++){
        f_up_sig[i]->FixParameter(j,fit_result_sig_final->Parameter(j));
        f_dn_sig[i]->FixParameter(j,fit_result_sig_final->Parameter(j));
        f_up_side[i]->FixParameter(j,fit_result_side_final->Parameter(j));
        f_dn_side[i]->FixParameter(j,fit_result_side_final->Parameter(j));

        if(i==j){
  	f_up_sig[i]->FixParameter(j,fit_result_sig_final->Parameter(j) + err_pars_up_sig[j]);
  	f_dn_sig[i]->FixParameter(j,fit_result_sig_final->Parameter(j) - err_pars_dn_sig[j]);
  	f_up_side[i]->FixParameter(j,fit_result_side_final->Parameter(j) + err_pars_up_side[j]);
  	f_dn_side[i]->FixParameter(j,fit_result_side_final->Parameter(j) - err_pars_dn_side[j]);
        }
      }
    }

    // Set up vectors to build the fits with uncertainties from
    vector<double> x, y_sig, y_side, quad_err_up_sig, quad_err_dn_sig, quad_err_up_side, quad_err_dn_side, err_x, err_up_sig, err_dn_sig, err_up_side, err_dn_side;

    // Look at variations in steps of 1GeV in ST
    const int n_points = final_range - fit_min;

    //Loop over all n_points points
    for(int i=0; i<n_points; i++){

      x.emplace_back(i+fit_min);
      err_x.emplace_back(0.5);
      y_sig.emplace_back(f_LandauExp_NewParam_Sig_final->Eval(x[i]));
      y_side.emplace_back(f_LandauExp_NewParam_Side_final->Eval(x[i]));
      quad_err_up_sig.emplace_back(0.);
      quad_err_dn_sig.emplace_back(0.);
      quad_err_up_side.emplace_back(0.);
      quad_err_dn_side.emplace_back(0.);

      //Loop over all variations and add them in quadrature
      for(int j=0; j<n_param; j++){
        double tmp_err_up_sig = (f_up_sig[j]->Eval(x[i]) - y_sig[i]);
        double tmp_err_dn_sig = (f_dn_sig[j]->Eval(x[i]) - y_sig[i]);
        double tmp_err_up_side = (f_up_side[j]->Eval(x[i]) - y_side[i]);
        double tmp_err_dn_side = (f_dn_side[j]->Eval(x[i]) - y_side[i]);

        if(tmp_err_up_sig < 0)  quad_err_dn_sig[i] += pow(tmp_err_up_sig,2);
        else                    quad_err_up_sig[i] += pow(tmp_err_up_sig,2);
        if(tmp_err_dn_sig < 0)  quad_err_dn_sig[i] += pow(tmp_err_dn_sig,2);
        else                    quad_err_up_sig[i] += pow(tmp_err_dn_sig,2);
        if(tmp_err_up_side < 0) quad_err_dn_side[i] += pow(tmp_err_up_side,2);
        else                    quad_err_up_side[i] += pow(tmp_err_up_side,2);
        if(tmp_err_dn_side < 0) quad_err_dn_side[i] += pow(tmp_err_dn_side,2);
        else                    quad_err_up_side[i] += pow(tmp_err_dn_side,2);
      }

      //take sqrt to obtain final error
      err_up_sig.emplace_back(sqrt(quad_err_up_sig[i]));
      err_dn_sig.emplace_back(sqrt(quad_err_dn_sig[i]));
      err_up_side.emplace_back(sqrt(quad_err_up_side[i]));
      err_dn_side.emplace_back(sqrt(quad_err_dn_side[i]));

    }


    // Set up vectors containing the ratio-values
    vector<double> y_ratio_sig, y_ratio_side, err_up_ratio_sig, err_dn_ratio_sig, err_up_ratio_side, err_dn_ratio_side;
    for(int i=0; i<n_points; i++){
      y_ratio_sig.emplace_back(1.);
      y_ratio_side.emplace_back(1.);
      err_up_ratio_sig.emplace_back(err_up_sig[i]/y_sig[i]);
      err_dn_ratio_sig.emplace_back(err_dn_sig[i]/y_sig[i]);
      err_up_ratio_side.emplace_back(err_up_side[i]/y_side[i]);
      err_dn_ratio_side.emplace_back(err_dn_side[i]/y_side[i]);
    }

    // Set up TGraphs to hold final fits and ratios with uncertainties
    TGraphAsymmErrors* g_fit_final_sig = new TGraphAsymmErrors(n_points, &x[0], &y_sig[0], &err_x[0], &err_x[0], &err_up_sig[0], &err_dn_sig[0]);
    TGraphAsymmErrors* g_fit_final_side = new TGraphAsymmErrors(n_points, &x[0], &y_side[0], &err_x[0], &err_x[0], &err_up_side[0], &err_dn_side[0]);
    TGraphAsymmErrors* g_fit_final_ratio_sig = new TGraphAsymmErrors(n_points, &x[0], &y_ratio_sig[0], &err_x[0], &err_x[0], &err_up_ratio_sig[0], &err_dn_ratio_sig[0]);
    TGraphAsymmErrors* g_fit_final_ratio_side = new TGraphAsymmErrors(n_points, &x[0], &y_ratio_side[0], &err_x[0], &err_x[0], &err_up_ratio_side[0], &err_dn_ratio_side[0]);

    //Calculate ratios for histo
    for(int i=1; i<n_bins+1; i++){
      double nom_sig = h_ht_sig_uncert_ratio->GetBinContent(i);
      double nom_side = h_ht_side_uncert_ratio->GetBinContent(i);
      double err_sig = h_ht_sig_uncert_ratio->GetBinError(i);
      double err_side = h_ht_side_uncert_ratio->GetBinError(i);
      double denom_sig = g_fit_final_sig->Eval(h_ht_sig_uncert_ratio->GetBinCenter(i));
      double denom_side = g_fit_final_side->Eval(h_ht_side_uncert_ratio->GetBinCenter(i));
      h_ht_sig_uncert_ratio->SetBinContent(i,nom_sig/denom_sig);
      h_ht_sig_uncert_ratio->SetBinError(i,err_sig/denom_sig);
      h_ht_side_uncert_ratio->SetBinContent(i,nom_side/denom_side);
      h_ht_side_uncert_ratio->SetBinError(i,err_side/denom_side);
    }

    TF1* f_const = new TF1("f_const","1", fit_min, fit_max);

    //Do cosmetics
    g_fit_final_sig->SetFillColor(kGray+1);
    g_fit_final_sig->SetFillStyle(1001);
    g_fit_final_side->SetFillColor(kGray+1);
    g_fit_final_side->SetFillStyle(1001);
    g_fit_final_ratio_sig->SetFillColor(kGray+1);
    g_fit_final_ratio_sig->SetFillStyle(1001);
    g_fit_final_ratio_side->SetFillColor(kGray+1);
    g_fit_final_ratio_side->SetFillStyle(1001);

    HistCosmetics(h_ht_sig_uncert,false);
    HistCosmetics(h_ht_side_uncert,false);
    HistCosmetics(h_ht_sig_uncert_ratio,true);
    HistCosmetics(h_ht_side_uncert_ratio,true);

    f_LandauExp_NewParam_Sig_final->SetLineColor(2);
    f_LandauExp_NewParam_Side_final->SetLineColor(2);

    f_const->SetLineColor(2);

    h_ht_sig_uncert->GetXaxis()->SetTitle("S_{T} [GeV]");
    h_ht_sig_uncert->GetXaxis()->SetRangeUser(0.,fit_max);
    h_ht_sig_uncert->GetXaxis()->SetTicks("+-");
    h_ht_sig_uncert->GetYaxis()->SetTitle("events");
    h_ht_sig_uncert_ratio->GetXaxis()->SetTitle("S_{T} [GeV]");
    h_ht_sig_uncert_ratio->GetXaxis()->SetRangeUser(0.,fit_max);
    h_ht_sig_uncert_ratio->GetYaxis()->CenterTitle();
    h_ht_sig_uncert_ratio->GetYaxis()->SetTitle("MC / fit");
    h_ht_sig_uncert_ratio->GetYaxis()->SetRangeUser(0.3,1.7);

    h_ht_side_uncert->GetXaxis()->SetTitle("S_{T} [GeV]");
    h_ht_side_uncert->GetXaxis()->SetRangeUser(0.,fit_max);
    h_ht_side_uncert->GetXaxis()->SetTicks("+-");
    h_ht_side_uncert->GetYaxis()->SetTitle("events");
    h_ht_side_uncert_ratio->GetXaxis()->SetTitle("S_{T} [GeV]");
    h_ht_side_uncert_ratio->GetXaxis()->SetRangeUser(0.,fit_max);
    h_ht_side_uncert_ratio->GetYaxis()->CenterTitle();
    h_ht_side_uncert_ratio->GetYaxis()->SetTitle("MC / fit");
    h_ht_side_uncert_ratio->GetYaxis()->SetRangeUser(0.3,1.7);


    TCanvas* c_sig_uncert = new TCanvas("c_sig_uncert", "Signal region fit with uncertainties", 400, 400);
    p_top = SetupRatioPadTop();
    p_bot = SetupRatioPad();
    p_top->Draw();
    p_bot->Draw();
    p_top->cd();
    h_ht_sig_uncert->Draw("E");
    g_fit_final_sig->Draw("SAME 3");
    f_LandauExp_NewParam_Sig_final->Draw("SAME");
    h_ht_sig_uncert->Draw("E SAME");
    TPaveText* textbox_sig = new TPaveText(0.6,0.77,0.75,0.87,"NDC");
    textbox_sig->SetFillColor(0);
    textbox_sig->SetLineColor(0);
    TText *line_sig = textbox_sig->AddText( "Signal region" );
    line_sig->SetTextColor(1);
    line_sig->SetTextAlign(12);//12
    line_sig->SetTextFont(43);
    line_sig->SetTextSizePixels(20);
    textbox_sig->SetBorderSize(1);
    textbox_sig->Draw();

    p_bot->cd();
    h_ht_sig_uncert_ratio->Draw("E");
    g_fit_final_ratio_sig->Draw("SAME 3");
    f_const->Draw("SAME");
    h_ht_sig_uncert_ratio->Draw("E SAME");
    c_sig_uncert->Print(path_out + "Alpha_Output/HT_Signal_Fit_ParameterUncert.eps");
    delete c_sig_uncert;
    //delete p_top;
    //delete p_bot;
    //delete textbox_sig;

    TCanvas* c_side_uncert = new TCanvas("c_side_uncert", "Control region fit with uncertainties", 400, 400);
    p_top = SetupRatioPadTop();
    p_bot = SetupRatioPad();
    p_top->Draw();
    p_bot->Draw();
    p_top->cd();
    h_ht_side_uncert->Draw("E");
    g_fit_final_side->Draw("SAME 3");
    f_LandauExp_NewParam_Side_final->Draw("SAME");
    h_ht_side_uncert->Draw("E SAME");
    TPaveText* textbox_side = new TPaveText(0.6,0.77,0.75,0.87,"NDC");
    textbox_side->SetFillColor(0);
    textbox_side->SetLineColor(0);
    TText *line_side = textbox_side->AddText( "Control region" );
    line_side->SetTextColor(1);
    line_side->SetTextAlign(12);//12
    line_side->SetTextFont(43);
    line_side->SetTextSizePixels(20);
    textbox_side->SetBorderSize(1);
    textbox_side->Draw();

    p_bot->cd();
    h_ht_side_uncert_ratio->Draw("E");
    g_fit_final_ratio_side->Draw("SAME 3");
    f_const->Draw("SAME");
    h_ht_side_uncert_ratio->Draw("E SAME");
    c_side_uncert->Print(path_out + "Alpha_Output/HT_Sideband_Fit_ParameterUncert.eps");
    delete c_side_uncert;
    //delete p_top;
    //delete p_bot;
    //delete textbox_side;

  //
  //
  //   //      Create all about Alpha
  //   //      ======================
  //
    TString formula_alpha = formula;
    TString formula_denom = formula;
    for(int i=0; i<n_param; i++){
      TString lookup = "[";
      lookup += i;
      lookup += "]";
      TString replace = "[";
      replace += i+n_param;
      replace += "]";
      formula_denom.ReplaceAll(lookup,replace);
    }
    formula_alpha += " / ";
    formula_alpha += formula_denom;
    TF1* f_alpha = new TF1("f_alpha", formula_alpha, fit_min, final_range);
    /*
    TF1* f_alpha_alt = new TF1("f_alpha_alt", "(landaun * (1-(1+TMath::Erf((x-[3])/[4]))/2)+((1+TMath::Erf((x-[3])/[4]))/2 * exp([5]*x+[6]+[7]*x*x))) / (landaun(8) * (1-(1+TMath::Erf((x-[11])/[12]))/2)+((1+TMath::Erf((x-[11])/[12]))/2 * exp([13]*x+[14]+[15]*x*x)))", fit_min, final_range);
    f_alpha_alt->FixParameter(0,1.10629e+06);
    f_alpha_alt->FixParameter(1,508.407);
    f_alpha_alt->FixParameter(2,86.692);
    f_alpha_alt->FixParameter(3,1098.42);
    f_alpha_alt->FixParameter(4,460.719);
    f_alpha_alt->FixParameter(5,-0.00099059);
    f_alpha_alt->FixParameter(6,5.67843);
    f_alpha_alt->FixParameter(7,-4.7281e-07);
    f_alpha_alt->FixParameter(8,780614);
    f_alpha_alt->FixParameter(9,515.981);
    f_alpha_alt->FixParameter(10,90.0114);
    f_alpha_alt->FixParameter(11,886.794);
    f_alpha_alt->FixParameter(12,355.456);
    f_alpha_alt->FixParameter(13,-0.00358004);
    f_alpha_alt->FixParameter(14,8.58543);
    f_alpha_alt->FixParameter(15,6.75657e-08);
    */
    /*
    TF1* f_alpha_alt = new TF1("f_alpha_alt", "(1-TMath::Erf((x-  [0]  )/  [1]  ))*  [2]   +   [3]  ", fit_min, final_range);
    f_alpha_alt->FixParameter(0,1218.33);
    f_alpha_alt->FixParameter(1,631.991);
    f_alpha_alt->FixParameter(2,0.213482);
    f_alpha_alt->FixParameter(3,1.05094);
    */

    // First n_param parameters belong to SR, 2nd half belongs to CR
    for(int i=0; i<n_param; i++)         f_alpha->FixParameter(i,f_LandauExp_FitSig->GetParameter(i));
    for(int i=n_param; i<2*n_param; i++) f_alpha->FixParameter(i,f_LandauExp_FitSide->GetParameter(i-n_param));


    //Set up vectors containing values for TGraphs
    vector<double> val_alpha, val_alpha_up, val_alpha_dn, err_alpha_up, err_alpha_dn;
    for(int i=0; i<n_points; i++){
      val_alpha.emplace_back(f_alpha->Eval(x[i]));
      err_alpha_up.emplace_back(sqrt( pow(err_up_sig[i]/y_side[i],2) + pow(y_sig[i]/y_side[i]/y_side[i] * err_dn_side[i],2) ));
      err_alpha_dn.emplace_back(sqrt( pow(err_dn_sig[i]/y_side[i],2) + pow(y_sig[i]/y_side[i]/y_side[i] * err_up_side[i],2) ));
      //err_alpha_up.emplace_back(sqrt( pow(err_up_sig[i]/y_side[i],2) + pow(y_sig[i]/y_side[i]/y_side[i] * err_dn_side[i],2) + pow(fabs(f_alpha->Eval(x[i]) - f_alpha_alt->Eval(x[i])),2)));
      //err_alpha_dn.emplace_back(sqrt( pow(err_dn_sig[i]/y_side[i],2) + pow(y_sig[i]/y_side[i]/y_side[i] * err_up_side[i],2) + pow(fabs(f_alpha->Eval(x[i]) - f_alpha_alt->Eval(x[i])),2)));
      val_alpha_up.emplace_back(val_alpha[i] + err_alpha_up[i]);
      val_alpha_dn.emplace_back(val_alpha[i] - err_alpha_dn[i]);
    }

    //Set up TGraphs (nom, up, down) of alpha, errors do not matter for up and dn variations, there only nominal values are used

    TGraphAsymmErrors* g_alpha = new TGraphAsymmErrors(n_points, &x[0], &val_alpha[0], &err_x[0], &err_x[0], &err_alpha_dn[0], &err_alpha_up[0]);
    TGraphAsymmErrors* g_alpha_up = new TGraphAsymmErrors(n_points, &x[0], &val_alpha_up[0], &err_x[0], &err_x[0], &err_alpha_dn[0], &err_alpha_up[0]);
    TGraphAsymmErrors* g_alpha_dn = new TGraphAsymmErrors(n_points, &x[0], &val_alpha_dn[0], &err_x[0], &err_x[0], &err_alpha_dn[0], &err_alpha_up[0]);

    //Do cosmetics
    g_alpha->SetFillColor(kGray+1);
    g_alpha->SetFillStyle(1001);
    HistCosmetics(g_alpha);
    HistCosmetics(h_ht_sig_side_ratio);
    HistCosmetics(h_ht_sig_side_ratio_rebin);
    g_alpha->GetXaxis()->SetTitle("S_{T} [GeV]");
    g_alpha->GetYaxis()->SetTitle("SR / CR");
    g_alpha->GetYaxis()->SetRangeUser(0,3);
    g_alpha->GetXaxis()->SetRangeUser(0,fit_max);
    g_alpha->Draw("a3");

    TCanvas* c_alpha = new TCanvas("c_alpha", "Extrapolation function", 400, 400);
    TPad* pad = SetupPad();
    pad->Draw();
    pad->cd();
    g_alpha->Draw("a3");
    f_alpha->Draw("SAME");
    //f_alpha_alt->SetLineColor(kBlue+1);
    //f_alpha_alt->Draw("SAME");
    h_ht_sig_side_ratio->Draw("SAME");

    TLegend* leg = new TLegend(0.25, 0.8, .4, .9);
    leg->SetBorderSize(0);
    leg->SetFillColor(10);
    leg->SetLineColor(1);
    leg->SetTextFont(43);
    leg->SetTextSize(14);
    leg->AddEntry(f_alpha, "Extrapolation function", "l");
    leg->AddEntry(h_ht_sig_side_ratio, "SR / CR (MC)", "le");
    leg->Draw();
    c_alpha->Print(path_out + "Alpha_Output/Alpha.eps");
    delete c_alpha;

    TCanvas* c_alpha_rebin = new TCanvas("c_alpha_rebin", "Extrapolation function", 400, 400);
    pad = SetupPad();
    pad->Draw();
    pad->cd();
    g_alpha->Draw("a3");
    f_alpha->Draw("SAME");
    h_ht_sig_side_ratio_rebin->Draw("SAME");

    leg = new TLegend(0.25, 0.8, .4, .9);
    leg->SetBorderSize(0);
    leg->SetFillColor(10);
    leg->SetLineColor(1);
    leg->SetTextFont(43);
    leg->SetTextSize(14);
    leg->AddEntry(f_alpha, "Extrapolation function", "l");
    leg->AddEntry(h_ht_sig_side_ratio_rebin, "SR / CR (MC)", "le");
    leg->Draw();
    c_alpha_rebin->Print(path_out + "Alpha_Output/Alpha_rebin.eps");
    delete c_alpha_rebin;

  //
  //   //      Apply to histo in final binning (only nominal), and do closure test with toy CR
  //   //      ===============================================================================

    // if(!do_systematics){
    //
    //   //Set up vectors containing information for TGraphs
    //   vector<double> x_crosscheck, x_err_crosscheck, y_sig_crosscheck, y_side_crosscheck, y_sig_stat_up_crosscheck, y_sig_stat_dn_crosscheck, y_side_stat_up_crosscheck, y_side_stat_dn_crosscheck, y_side_total_up_crosscheck, y_side_total_dn_crosscheck;
    //   for(int i=1; i<h_ht_sig_crosscheck->GetNbinsX()+1; i++){
    //     x_crosscheck.emplace_back(h_ht_sig_crosscheck->GetBinCenter(i));
    //     x_err_crosscheck.emplace_back(h_ht_sig_crosscheck->GetBinWidth(i)/2);
    //     y_sig_crosscheck.emplace_back(h_ht_sig_crosscheck->GetBinContent(i));
    //     y_sig_stat_up_crosscheck.emplace_back(h_ht_sig_crosscheck->GetBinError(i));
    //     y_sig_stat_dn_crosscheck.emplace_back(h_ht_sig_crosscheck->GetBinError(i));
    //     y_side_crosscheck.emplace_back(h_ht_side_crosscheck->GetBinContent(i) * f_alpha->Eval(x_crosscheck[i-1]));
    //     y_side_stat_up_crosscheck.emplace_back(h_ht_side_crosscheck->GetBinError(i)* f_alpha->Eval(x_crosscheck[i-1]));
    //     y_side_stat_dn_crosscheck.emplace_back(h_ht_side_crosscheck->GetBinError(i)* f_alpha->Eval(x_crosscheck[i-1]));
    //     double err_from_alpha_up = fabs(h_ht_side_crosscheck->GetBinContent(i) * g_alpha->Eval(x_crosscheck[i-1]) - h_ht_side_crosscheck->GetBinContent(i) * g_alpha_up->Eval(x_crosscheck[i-1]));
    //     double err_from_alpha_dn = fabs(h_ht_side_crosscheck->GetBinContent(i) * g_alpha->Eval(x_crosscheck[i-1]) - h_ht_side_crosscheck->GetBinContent(i) * g_alpha_dn->Eval(x_crosscheck[i-1]));
    //     y_side_total_up_crosscheck.emplace_back(sqrt(pow(y_side_stat_up_crosscheck[i-1],2) + pow(err_from_alpha_up,2)));
    //     y_side_total_dn_crosscheck.emplace_back(sqrt(pow(y_side_stat_dn_crosscheck[i-1],2) + pow(err_from_alpha_dn,2)));
    //   }
    //
    //   vector<double> x_closure, x_err_closure, y_sig_closure, y_side_closure, y_sig_stat_up_closure, y_sig_stat_dn_closure, y_side_stat_up_closure, y_side_stat_dn_closure, y_side_total_up_closure, y_side_total_dn_closure;
    //   for(int i=1; i<n_bins+1; i++){
    //     x_closure.emplace_back(h_ht_sig_closure->GetBinCenter(i));
    //     x_err_closure.emplace_back(h_ht_sig_closure->GetBinWidth(i)/2);
    //     y_sig_closure.emplace_back(h_ht_sig_closure->GetBinContent(i));
    //     y_sig_stat_up_closure.emplace_back(h_ht_sig_closure->GetBinError(i));
    //     y_sig_stat_dn_closure.emplace_back(h_ht_sig_closure->GetBinError(i));
    //     y_side_closure.emplace_back(h_ht_side_closure->GetBinContent(i) * f_alpha->Eval(x_closure[i-1]));
    //     y_side_stat_up_closure.emplace_back(h_ht_side_closure->GetBinError(i)* f_alpha->Eval(x_closure[i-1]));
    //     y_side_stat_dn_closure.emplace_back(h_ht_side_closure->GetBinError(i)* f_alpha->Eval(x_closure[i-1]));
    //     double err_from_alpha_up = fabs(h_ht_side_closure->GetBinContent(i) * g_alpha->Eval(x_closure[i-1]) - h_ht_side_closure->GetBinContent(i) * g_alpha_up->Eval(x_closure[i-1]));
    //     double err_from_alpha_dn = fabs(h_ht_side_closure->GetBinContent(i) * g_alpha->Eval(x_closure[i-1]) - h_ht_side_closure->GetBinContent(i) * g_alpha_dn->Eval(x_closure[i-1]));
    //     y_side_total_up_closure.emplace_back(sqrt(pow(y_side_stat_up_closure[i-1],2) + pow(err_from_alpha_up,2)));
    //     y_side_total_dn_closure.emplace_back(sqrt(pow(y_side_stat_dn_closure[i-1],2) + pow(err_from_alpha_dn,2)));
    //   }
    //
    //   //Calculate ratios
    //   vector<double> y_sig_crosscheck_ratio, y_side_crosscheck_ratio, y_sig_stat_up_crosscheck_ratio, y_sig_stat_dn_crosscheck_ratio, y_side_stat_up_crosscheck_ratio, y_side_stat_dn_crosscheck_ratio, y_side_total_up_crosscheck_ratio, y_side_total_dn_crosscheck_ratio;
    //   for(int i=1; i<h_ht_sig_crosscheck->GetNbinsX()+1; i++){
    //     y_sig_crosscheck_ratio.emplace_back(y_sig_crosscheck[i-1] / y_sig_crosscheck[i-1]);
    //     y_side_crosscheck_ratio.emplace_back(y_side_crosscheck[i-1] / y_sig_crosscheck[i-1]);
    //     y_sig_stat_up_crosscheck_ratio.emplace_back(y_sig_stat_up_crosscheck[i-1] / y_sig_crosscheck[i-1]);
    //     y_sig_stat_dn_crosscheck_ratio.emplace_back(y_sig_stat_dn_crosscheck[i-1] / y_sig_crosscheck[i-1]);
    //     y_side_stat_up_crosscheck_ratio.emplace_back(y_side_stat_up_crosscheck[i-1] / y_sig_crosscheck[i-1]);
    //     y_side_stat_dn_crosscheck_ratio.emplace_back(y_side_stat_dn_crosscheck[i-1] / y_sig_crosscheck[i-1]);
    //     y_side_total_up_crosscheck_ratio.emplace_back(y_side_total_up_crosscheck[i-1] / y_sig_crosscheck[i-1]);
    //     y_side_total_dn_crosscheck_ratio.emplace_back(y_side_total_dn_crosscheck[i-1] / y_sig_crosscheck[i-1]);
    //   }
    //
    //   vector<double> y_sig_closure_ratio, y_side_closure_ratio, y_sig_stat_up_closure_ratio, y_sig_stat_dn_closure_ratio, y_side_stat_up_closure_ratio, y_side_stat_dn_closure_ratio, y_side_total_up_closure_ratio, y_side_total_dn_closure_ratio;
    //   for(int i=1; i<n_bins+1; i++){
    //     y_sig_closure_ratio.emplace_back(y_sig_closure[i-1] / y_sig_closure[i-1]);
    //     y_side_closure_ratio.emplace_back(y_side_closure[i-1] / y_sig_closure[i-1]);
    //     y_sig_stat_up_closure_ratio.emplace_back(y_sig_stat_up_closure[i-1] / y_sig_closure[i-1]);
    //     y_sig_stat_dn_closure_ratio.emplace_back(y_sig_stat_dn_closure[i-1] / y_sig_closure[i-1]);
    //     y_side_stat_up_closure_ratio.emplace_back(y_side_stat_up_closure[i-1] / y_sig_closure[i-1]);
    //     y_side_stat_dn_closure_ratio.emplace_back(y_side_stat_dn_closure[i-1] / y_sig_closure[i-1]);
    //     y_side_total_up_closure_ratio.emplace_back(y_side_total_up_closure[i-1] / y_sig_closure[i-1]);
    //     y_side_total_dn_closure_ratio.emplace_back(y_side_total_dn_closure[i-1] / y_sig_closure[i-1]);
    //   }
    //
    //
    //   TGraphAsymmErrors* g_sig_crosscheck_stat = new TGraphAsymmErrors(h_ht_sig_crosscheck->GetNbinsX(), &x_crosscheck[0], &y_sig_crosscheck[0], &x_err_crosscheck[0], &x_err_crosscheck[0], &y_sig_stat_dn_crosscheck[0], &y_sig_stat_up_crosscheck[0]);
    //   TGraphAsymmErrors* g_side_crosscheck_stat = new TGraphAsymmErrors(h_ht_side_crosscheck->GetNbinsX(), &x_crosscheck[0], &y_side_crosscheck[0], &x_err_crosscheck[0], &x_err_crosscheck[0], &y_side_stat_dn_crosscheck[0], &y_side_stat_up_crosscheck[0]);
    //   TGraphAsymmErrors* g_side_crosscheck_total = new TGraphAsymmErrors(h_ht_side_crosscheck->GetNbinsX(), &x_crosscheck[0], &y_side_crosscheck[0], &x_err_crosscheck[0], &x_err_crosscheck[0], &y_side_total_dn_crosscheck[0], &y_side_total_up_crosscheck[0]);
    //   TGraphAsymmErrors* g_sig_crosscheck_stat_ratio = new TGraphAsymmErrors(h_ht_sig_crosscheck->GetNbinsX(), &x_crosscheck[0], &y_sig_crosscheck_ratio[0], &x_err_crosscheck[0], &x_err_crosscheck[0], &y_sig_stat_dn_crosscheck_ratio[0], &y_sig_stat_up_crosscheck_ratio[0]);
    //   TGraphAsymmErrors* g_side_crosscheck_stat_ratio = new TGraphAsymmErrors(h_ht_side_crosscheck->GetNbinsX(), &x_crosscheck[0], &y_side_crosscheck_ratio[0], &x_err_crosscheck[0], &x_err_crosscheck[0], &y_side_stat_dn_crosscheck_ratio[0], &y_side_stat_up_crosscheck_ratio[0]);
    //   TGraphAsymmErrors* g_side_crosscheck_total_ratio = new TGraphAsymmErrors(h_ht_side_crosscheck->GetNbinsX(), &x_crosscheck[0], &y_side_crosscheck_ratio[0], &x_err_crosscheck[0], &x_err_crosscheck[0], &y_side_total_dn_crosscheck_ratio[0], &y_side_total_up_crosscheck_ratio[0]);
    //
    //
    //
    //   //Do cosmetics
    //   HistCosmetics(g_sig_crosscheck_stat,false);
    //   HistCosmetics(g_side_crosscheck_stat,false);
    //   HistCosmetics(g_side_crosscheck_total,false);
    //   HistCosmetics(g_sig_crosscheck_stat_ratio,true);
    //   HistCosmetics(g_side_crosscheck_stat_ratio,true);
    //   HistCosmetics(g_side_crosscheck_total_ratio,true);
    //
    //   TCanvas* c_crosscheck = new TCanvas("c_crosscheck", "Crosscheck", 400, 400);
    //   p_top = SetupRatioPadTop();
    //   p_bot = SetupRatioPad();
    //   p_top->Draw();
    //   p_bot->Draw();
    //   p_top->cd();
    //   g_sig_crosscheck_stat->GetXaxis()->SetTitle("S_{T} [GeV]");
    //   g_sig_crosscheck_stat->GetYaxis()->SetTitle("events");
    //   g_sig_crosscheck_stat->GetXaxis()->SetRangeUser(0,final_range);
    //   g_sig_crosscheck_stat->SetLineColor(kBlue);
    //   g_side_crosscheck_stat->SetLineColor(kRed);
    //   g_side_crosscheck_total->SetLineColor(kRed);
    //   g_sig_crosscheck_stat_ratio->GetXaxis()->SetTitle("S_{T} [GeV]");
    //   g_sig_crosscheck_stat_ratio->GetYaxis()->CenterTitle();
    //   g_sig_crosscheck_stat_ratio->GetYaxis()->SetTitle("Extr. / SR");
    //   g_sig_crosscheck_stat_ratio->SetMinimum(0.3);
    //   g_sig_crosscheck_stat_ratio->SetMaximum(1.7);
    //   g_sig_crosscheck_stat_ratio->GetXaxis()->SetRangeUser(0,bins_crosscheck[n_bins_crosscheck+1]);
    //   g_sig_crosscheck_stat_ratio->SetLineColor(kBlue);
    //   g_side_crosscheck_stat_ratio->SetLineColor(kRed);
    //   g_side_crosscheck_total_ratio->SetLineColor(kRed);
    //
    //   g_sig_crosscheck_stat->Draw("ap");
    //   g_side_crosscheck_stat->Draw("p same");
    //   g_side_crosscheck_total->Draw("p same");
    //   leg = new TLegend(0.5, 0.73, .85, .9);
    //   leg->SetBorderSize(0);
    //   leg->SetFillColor(10);
    //   leg->SetLineColor(1);
    //   leg->SetTextFont(43);
    //   leg->SetTextSize(12);
    //   leg->AddEntry(g_sig_crosscheck_stat, "SR MC", "l");
    //   leg->AddEntry(g_side_crosscheck_total, "Extrapolation from CR MC", "l");
    //   leg->Draw();
    //   p_bot->cd();
    //   g_sig_crosscheck_stat_ratio->Draw("ap");
    //   g_side_crosscheck_stat_ratio->Draw("p same");
    //   g_side_crosscheck_total_ratio->Draw("p same");
    //
    //   c_crosscheck->Print(path_out + "Alpha_Output/Crosscheck.eps");
    //   delete c_crosscheck;
    //   //delete leg;
    //   //delete p_bot;
    //   //delete p_top;
    //
    //
    //   delete g_side_crosscheck_total_ratio;
    //   delete g_side_crosscheck_stat_ratio;
    //   delete g_sig_crosscheck_stat_ratio;
    //   delete g_side_crosscheck_total;
    //   delete g_side_crosscheck_stat;
    //   delete g_sig_crosscheck_stat;
    //
    //
    //   TGraphAsymmErrors* g_sig_closure_stat = new TGraphAsymmErrors(h_ht_sig_closure->GetNbinsX(), &x_closure[0], &y_sig_closure[0], &x_err_closure[0], &x_err_closure[0], &y_sig_stat_dn_closure[0], &y_sig_stat_up_closure[0]);
    //   TGraphAsymmErrors* g_side_closure_stat = new TGraphAsymmErrors(h_ht_side_closure->GetNbinsX(), &x_closure[0], &y_side_closure[0], &x_err_closure[0], &x_err_closure[0], &y_side_stat_dn_closure[0], &y_side_stat_up_closure[0]);
    //   TGraphAsymmErrors* g_side_closure_total = new TGraphAsymmErrors(h_ht_side_closure->GetNbinsX(), &x_closure[0], &y_side_closure[0], &x_err_closure[0], &x_err_closure[0], &y_side_total_dn_closure[0], &y_side_total_up_closure[0]);
    //   TGraphAsymmErrors* g_sig_closure_stat_ratio = new TGraphAsymmErrors(h_ht_sig_closure->GetNbinsX(), &x_closure[0], &y_sig_closure_ratio[0], &x_err_closure[0], &x_err_closure[0], &y_sig_stat_dn_closure_ratio[0], &y_sig_stat_up_closure_ratio[0]);
    //   TGraphAsymmErrors* g_side_closure_stat_ratio = new TGraphAsymmErrors(h_ht_side_closure->GetNbinsX(), &x_closure[0], &y_side_closure_ratio[0], &x_err_closure[0], &x_err_closure[0], &y_side_stat_dn_closure_ratio[0], &y_side_stat_up_closure_ratio[0]);
    //   TGraphAsymmErrors* g_side_closure_total_ratio = new TGraphAsymmErrors(h_ht_side_closure->GetNbinsX(), &x_closure[0], &y_side_closure_ratio[0], &x_err_closure[0], &x_err_closure[0], &y_side_total_dn_closure_ratio[0], &y_side_total_up_closure_ratio[0]);
    //
    //
    //
    //   //Do cosmetics
    //   HistCosmetics(g_sig_closure_stat,false);
    //   HistCosmetics(g_side_closure_stat,false);
    //   HistCosmetics(g_side_closure_total,false);
    //   HistCosmetics(g_sig_closure_stat_ratio,true);
    //   HistCosmetics(g_side_closure_stat_ratio,true);
    //   HistCosmetics(g_side_closure_total_ratio,true);
    //
    //   TCanvas* c_closure = new TCanvas("c_closure", "Closure", 400, 400);
    //   p_top = SetupRatioPadTop();
    //   p_bot = SetupRatioPad();
    //   p_top->Draw();
    //   p_bot->Draw();
    //   p_top->cd();
    //   g_sig_closure_stat->GetXaxis()->SetTitle("S_{T} [GeV]");
    //   g_sig_closure_stat->GetYaxis()->SetTitle("events");
    //   g_sig_closure_stat->GetXaxis()->SetRangeUser(0,fit_max);
    //   g_sig_closure_stat->SetLineColor(kBlue);
    //   g_side_closure_stat->SetLineColor(kRed);
    //   g_side_closure_total->SetLineColor(kRed);
    //   g_sig_closure_stat_ratio->GetXaxis()->SetTitle("S_{T} [GeV]");
    //   g_sig_closure_stat_ratio->GetYaxis()->CenterTitle();
    //   g_sig_closure_stat_ratio->GetYaxis()->SetTitle("Extr. / SR");
    //   g_sig_closure_stat_ratio->SetMinimum(0.3);
    //   g_sig_closure_stat_ratio->SetMaximum(1.7);
    //   g_sig_closure_stat_ratio->GetXaxis()->SetRangeUser(0,fit_max);
    //   g_sig_closure_stat_ratio->SetLineColor(kBlue);
    //   g_side_closure_stat_ratio->SetLineColor(kRed);
    //   g_side_closure_total_ratio->SetLineColor(kRed);
    //
    //   g_sig_closure_stat->Draw("ap");
    //   g_side_closure_stat->Draw("p same");
    //   g_side_closure_total->Draw("p same");
    //   g_sig_closure_stat->SetMinimum(0.00001);
    //   p_top->SetLogy();
    //   leg = new TLegend(0.5, 0.73, .85, .9);
    //   leg->SetBorderSize(0);
    //   leg->SetFillColor(10);
    //   leg->SetLineColor(1);
    //   leg->SetTextFont(43);
    //   leg->SetTextSize(12);
    //   leg->AddEntry(g_sig_closure_stat, "SR MC", "l");
    //   leg->AddEntry(g_side_closure_total, "Extrapolation from toy MC", "l");
    //   leg->Draw();
    //   p_bot->cd();
    //   g_sig_closure_stat_ratio->Draw("ap");
    //   g_side_closure_stat_ratio->Draw("p same");
    //   g_side_closure_total_ratio->Draw("p same");
    //
    //   c_closure->Print(path_out + "Alpha_Output/ClosureTest.eps");
    //   delete c_closure;
    //   //delete leg;
    //   //delete p_bot;
    //   //delete p_top;
    //
    //
    //   delete g_side_closure_total_ratio;
    //   delete g_side_closure_stat_ratio;
    //   delete g_sig_closure_stat_ratio;
    //   delete g_side_closure_total;
    //   delete g_side_closure_stat;
    //   delete g_sig_closure_stat;
    // }
  //
  //
  //   //      Write output: Alpha function
  //   //      ============================
  //

    //Values corresponding to 'syst'
    // TString outpath;
    // if(do_systematics){
    //   if(is_valid_SystName){
    //     path_out += SystName+"_"+SystDirection;
    //   }
    //   else throw runtime_error("invalid SystName specified");
    //
    //   if(do_ttbardy) outpath = path_out + "/ExtrapolationSRCRFunction_TTbarDY.root";
    //   else           outpath = path_out + "/ExtrapolationSRCRFunction.root";
    //   cout << "outpath = " << outpath << endl;
    // }
    // else{
    //   if(do_ttbardy) outpath = path_out + "/NOMINAL/ExtrapolationSRCRFunction_TTbarDY.root";
    //   else           outpath = path_out + "/NOMINAL/ExtrapolationSRCRFunction.root";
    //   cout << "outpath = " << outpath << endl;
    // }

    TString  outpath = path_out + "ExtrapolationSRCRFunction_TTbarDY.root";
    TFile* output = new TFile(outpath,"RECREATE");
    g_alpha->Write();
    output->Close();
    delete output;
    // bool do_systematics = false;

    //if 'syst' = nominal: Also write syst Alpha variations
    // if(!do_systematics){
    //   // outpath = path_out + "Alpha_Output/ExtrapolationSRCRFunction_TTbarDY.root";
    //
    //
    //   TFile* output_alpha_dn = new TFile(pat,"RECREATE");
    //   g_alpha_dn->Write();
    //   output_alpha_dn->Close();
    //   path_out.ReplaceAll("ALPHA_down", "ALPHA_up");
    //   TFile* output_alpha_up = new TFile(path_out,"RECREATE");
    //   g_alpha_up->Write();
    //   output_alpha_up->Close();
    //
    //   delete output_alpha_dn;
    //   delete output_alpha_up;
    // }




  //

    delete g_alpha_up;
    delete g_alpha_dn;
    delete g_alpha;
    delete f_alpha;
    delete f_const;
    delete g_fit_final_ratio_side;
    delete g_fit_final_ratio_sig;
    delete g_fit_final_side;
    delete g_fit_final_sig;
    for(int i=0; i<n_param; i++){
      delete l_pars_dn_side[i];
      delete l_pars_up_side[i];
      delete l_pars_central_side[i];
    }
    for(unsigned int i=0; i<g_chi2_pars_side.size(); i++) delete g_chi2_pars_side[i];
    delete f_Sideband_errors;
    for(int i=0; i<n_param; i++){
      delete l_pars_dn_sig[i];
      delete l_pars_up_sig[i];
      delete l_pars_central_sig[i];
    }
    for(unsigned int i=0; i<g_chi2_pars_sig.size(); i++) delete g_chi2_pars_sig[i];
    delete f_Signal_errors;
    delete f_LandauExp_NewParam_Sig_final;
    delete f_LandauExp_NewParam_Side_final;
    delete f_LandauExp_NewParam_Side_intermediate;
    delete f_LandauExp_NewParam_Sig_intermediate;
    delete h_ht_sig_side_ratio;
    //delete h_ht_sig_side_ratio_rebin;
    // if(!do_systematics){
    //   delete h_ht_side_closure_ratio;
    //   delete h_ht_sig_closure_ratio;
    //   delete h_ht_side_crosscheck_ratio;
    //   delete h_ht_sig_crosscheck_ratio;
    // }
    delete h_ht_side_uncert_ratio;
    delete h_ht_sig_uncert_ratio;
    delete h_ht_side_final_ratio;
    delete h_ht_sig_final_ratio;
    delete h_ht_side_intermediate_ratio;
    delete h_ht_sig_intermediate_ratio;
    delete h_ht_side_initial_ratio;
    delete h_ht_sig_initial_ratio;
    // if(!do_systematics){
    //   delete h_ht_side_closure;
    //   delete h_ht_sig_closure;
    //   delete h_ht_side_crosscheck;
    //   delete h_ht_sig_crosscheck;
    // }
    delete h_ht_side_uncert;
    delete h_ht_sig_uncert;
    delete h_ht_side_final;
    delete h_ht_sig_final;
    delete h_ht_side_intermediate;
    delete h_ht_sig_intermediate;
    delete h_ht_side_initial;
    delete h_ht_sig_initial;
    delete f_LandauExp_FitSide;
    delete f_LandauExp_FitSig;
    // if(!do_systematics) delete h_ht_sideband_closure_tmp;
    delete h_ht_signal_tmp;
    //delete h_ht_sideband_tmp_rebin;
    delete h_ht_sideband_tmp;
    // delete Closure;
}
