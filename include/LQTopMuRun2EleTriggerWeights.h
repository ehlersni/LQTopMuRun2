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
  class ElectronTriggerWeights: public uhh2::AnalysisModule {

  public:
    explicit ElectronTriggerWeights(uhh2::Context & ctx);
    virtual bool process(uhh2::Event & event) override;

  private:
    //TString path,
    std::unique_ptr<TGraphAsymmErrors>  graph_datalowpt, graph_datahighpt;

  };
}
