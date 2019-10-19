#pragma once

#include "UHH2/core/include/fwd.h"
#include "UHH2/core/include/Selection.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/common/include/Utils.h"



#include "UHH2/common/include/ObjectIdUtils.h"
#include "UHH2/common/include/JetCorrections.h"
#include "UHH2/common/include/JetIds.h"
#include <TFile.h>
#include <TGraphAsymmErrors.h>
//#include "LHAPDF/LHAPDF.h"
#include "TSystem.h"

#include <vector>

namespace uhh2examples {

  /* Select events with at least two jets in which the leading two jets have deltaphi > 2.7 and the third jet pt is
  * below 20% of the average of the leading two jets, where the minimum deltaphi and
  * maximum third jet pt fraction can be changed in the constructor.
  * The jets are assumed to be sorted in pt.
  */

// my first try :(

  class ElectronTriggerWeights: public uhh2::AnalysisModule {

  public:
    explicit ElectronTriggerWeights(uhh2::Context & ctx, TString SysDirection_);
    virtual bool process(uhh2::Event & event) override;

  private:
    //TString path,
    TString SysDirection;
    std::unique_ptr<TGraphAsymmErrors>  graph_datalowpt, graph_datahighpt;

  };


// class ElectronTriggerWeights: public uhh2::AnalysisModule{
//
//  public:
//   explicit ElectronTriggerWeights(uhh2::Context & ctx, TString path_, TString SysDirection_);
//   virtual bool process(uhh2::Event & event) override;
//
//  private:
//   TString path, SysDirection;
//   std::unique_ptr<TGraphAsymmErrors> Eff_lowpt_MC, Eff_lowpt_DATA, Eff_highpt_MC, Eff_highpt_DATA;
//
// };



}
