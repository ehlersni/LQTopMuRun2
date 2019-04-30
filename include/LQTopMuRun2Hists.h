#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/LQTopMuRun2/include/LQReconstructionHypothesis.h"
#include "UHH2/LQTopMuRun2/include/LQReconstructionHypothesisDiscriminators.h"

namespace uhh2examples {

  /**  \brief Example class for booking and filling histograms
  *
  * NOTE: This class uses the 'hist' method to retrieve histograms.
  * This requires a string lookup and is therefore slow if you have
  * many histograms. Therefore, it is recommended to use histogram
  * pointers as member data instead, like in 'common/include/ElectronHists.h'.
  */
  class LQTopMuRun2Hists: public uhh2::Hists {
  public:
    // use the same constructor arguments as Hists for forwarding:
    LQTopMuRun2Hists(uhh2::Context & ctx, const std::string & dirname);

    virtual void fill(const uhh2::Event & ev) override;
    virtual ~LQTopMuRun2Hists();

  private:
    uhh2::Event::Handle<float> h_ST;
    uhh2::Event::Handle<bool> h_is_mlq_reconstructed;
    uhh2::Event::Handle<TString> h_mlq_reco_mode;
    uhh2::Event::Handle<std::vector<LQReconstructionHypothesis>> h_recohyps;
  };

}
