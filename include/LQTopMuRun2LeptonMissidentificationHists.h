#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/TTbarGen.h"
#include "UHH2/LQTopMuRun2/include/LQReconstructionHypothesisDiscriminators.h"
#include "UHH2/LQTopMuRun2/include/LQReconstructionHypothesis.h"
#include "UHH2/LQTopMuRun2/include/LQGen.h"
#include "UHH2/LQTopMuRun2/include/LQTopMuRun2Selections.h"


namespace uhh2examples {

class LQTopMuRun2LeptonMissidentificationHists: public uhh2::Hists {
public:
    // use the same constructor arguments as Hists for forwarding:
    LQTopMuRun2LeptonMissidentificationHists(uhh2::Context & ctx, const std::string & dirname, bool is_ele_);

    virtual void fill(const uhh2::Event & ev) override;

  protected:
    bool is_mc, is_ele;
    std::unique_ptr<ZEEFinder> ZEE_finder;


    virtual ~LQTopMuRun2LeptonMissidentificationHists();
};

}
