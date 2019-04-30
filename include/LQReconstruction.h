#pragma once

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/LQTopMuRun2/include/LQReconstructionHypothesis.h"
#include "TMinuit.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/ObjectIdUtils.h"
#include "UHH2/common/include/TopJetIds.h"

typedef std::function< std::vector<LorentzVector>  (const LorentzVector & lepton, const LorentzVector & met)> NeutrinoReconstructionMethod;



class LQPrimaryLepton: public uhh2::AnalysisModule {
public:
    explicit LQPrimaryLepton(uhh2::Context & ctx);

    virtual bool process(uhh2::Event & event) override;

    virtual ~LQPrimaryLepton();

private:
    uhh2::Event::Handle<FlavorParticle> h_primlep;
};



class HighMassInclusiveLQReconstruction: public uhh2::AnalysisModule {
public:

    explicit HighMassInclusiveLQReconstruction(uhh2::Context & ctx, const NeutrinoReconstructionMethod & neutrinofunction);

    virtual bool process(uhh2::Event & event) override;

    virtual ~HighMassInclusiveLQReconstruction();

private:
    NeutrinoReconstructionMethod m_neutrinofunction;
    uhh2::Event::Handle<std::vector<LQReconstructionHypothesis>> h_recohyps;
    uhh2::Event::Handle<bool> h_is_mlq_reconstructed;
    uhh2::Event::Handle<TString> h_mlq_reco_mode;
};



std::vector<LorentzVector> LQNeutrinoReconstruction(const LorentzVector & lepton, const LorentzVector & met);
