#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"


namespace uhh2examples {

class LQTopMuRun2EleTriggerEffHists : public uhh2::Hists {

  public:
    LQTopMuRun2EleTriggerEffHists(uhh2::Context & ctx, const std::string & dirname);
    virtual void fill(const uhh2::Event & ev) override;

  protected:
    virtual ~LQTopMuRun2EleTriggerEffHists();

};

}
