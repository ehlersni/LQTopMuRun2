#pragma once

#include "UHH2/core/include/fwd.h"
#include "UHH2/core/include/Selection.h"
#include "UHH2/core/include/GenParticle.h"

#include <vector>

namespace uhh2examples {

/* Select events with at least two jets in which the leading two jets have deltaphi > 2.7 and the third jet pt is
 * below 20% of the average of the leading two jets, where the minimum deltaphi and
 * maximum third jet pt fraction can be changed in the constructor.
 * The jets are assumed to be sorted in pt.
 */

  class DijetSelectiona: public uhh2::Selection {
  public:
    DijetSelectiona(float dphi_min = 2.7f, float third_frac_max = 0.2f);
    virtual bool passes(const uhh2::Event & event) override;
  private:
    float dphi_min, third_frac_max;
  };

  class HtSelectiona: public uhh2::Selection {
  public:
    explicit HtSelectiona(double ht_min=0., double ht_max=-1);
    virtual bool passes(const uhh2::Event & event);
  private:
    double ht_min, ht_max;
  };

  class InvMass2MuVetoa: public uhh2::Selection {
  public:
    explicit InvMass2MuVetoa(double m_min = 81., double m_max = 101.);
    virtual bool passes(const uhh2::Event & event);
  private:
    double m_min, m_max;
  };

  class InvMass2MuVetoInverteda: public uhh2::Selection {
  public:
    explicit InvMass2MuVetoInverteda(double m_min = 81., double m_max = 101.);
    virtual bool passes(const uhh2::Event & event);
  private:
    double m_min, m_max;
  };

  class PtLeadingMuonSelectiona: public uhh2::Selection {
  public:
    explicit PtLeadingMuonSelectiona(double pt_min = 30., double m_max = -1);
    virtual bool passes(const uhh2::Event & event);
  private:
    double pt_min, pt_max;
  };

  class Pt2ndMuonSelectiona: public uhh2::Selection {
  public:
    explicit Pt2ndMuonSelectiona(double pt_min = 30., double m_max = -1);
    virtual bool passes(const uhh2::Event & event);
  private:
    double pt_min, pt_max;
  };

  class PtLeadingJetSelectiona: public uhh2::Selection {
  public:
    explicit PtLeadingJetSelectiona(double pt_min = 30., double m_max = -1);
    virtual bool passes(const uhh2::Event & event);
  private:
    double pt_min, pt_max;
  };

  class Pt2ndJetSelectiona: public uhh2::Selection {
  public:
    explicit Pt2ndJetSelectiona(double pt_min = 30., double m_max = -1);
    virtual bool passes(const uhh2::Event & event);
  private:
    double pt_min, pt_max;
  };

  class PtRelMuJetSelectiona : public uhh2::Selection{
  public:
    explicit PtRelMuJetSelectiona(double ptrel_min = 0, double ptrel_max = -1);
    virtual bool passes(const uhh2::Event & event);
  private:
    double ptrel_min, ptrel_max;
  };

  class METSelectiona : public uhh2::Selection{
  public:
    explicit METSelectiona(double met_min = 0, double met_max = -1);
    virtual bool passes(const uhh2::Event & event);
  private:
    double met_min, met_max;
  };

  class GenLvlZMuMuSelectiona : public uhh2::Selection{
  public:
    explicit GenLvlZMuMuSelectiona();
    virtual bool passes(const uhh2::Event & event);
  };

  class GenLvlZEESelectiona : public uhh2::Selection{
  public:
    explicit GenLvlZEESelectiona();
    virtual bool passes(const uhh2::Event & event);
  };

  class GenLvlTopDileptonSelectiona : public uhh2::Selection{
  public:
    explicit GenLvlTopDileptonSelectiona();
    virtual bool passes(const uhh2::Event & event);
  };

  class PtRelMu1JetSelectiona : public uhh2::Selection{
  public:
    explicit PtRelMu1JetSelectiona(double ptrel_min = 0., double ptrel_max = -1);
    virtual bool passes(const uhh2::Event & event);
  private:
    double ptrel_min, ptrel_max;
  };

  class HTJetsSelectiona : public uhh2::Selection{
  public:
    explicit HTJetsSelectiona(double ht_min = 0., double ht_max = -1);
    virtual bool passes(const uhh2::Event & event);
  private:
    double ht_min, ht_max;
  };

  class HTLeptSelectiona : public uhh2::Selection{
  public:
    explicit HTLeptSelectiona(double ht_min = 0., double ht_max = -1);
    virtual bool passes(const uhh2::Event & event);
  private:
    double ht_min, ht_max;
  };

  class EtaLeadingJetSelectiona : public uhh2::Selection{
  public:
    explicit EtaLeadingJetSelectiona(double eta_max = -1);
    virtual bool passes(const uhh2::Event & event);
  private:
    double eta_max;
  };

  class InvMassMuEleVetoa: public uhh2::Selection {
  public:
    explicit InvMassMuEleVetoa(double m_min = 81., double m_max = 101.);
    virtual bool passes(const uhh2::Event & event);
  private:
    double m_min, m_max;
  };

  class InvMassEleEleVetoa: public uhh2::Selection {
  public:
    explicit InvMassEleEleVetoa(double m_min = 81., double m_max = 101.);
    virtual bool passes(const uhh2::Event & event);
  private:
    double m_min, m_max;
  };

  class InvMassEleEleSelectiona: public uhh2::Selection {
  public:
    explicit InvMassEleEleSelectiona(double m_min = 71., double m_max = 111.);
    virtual bool passes(const uhh2::Event & event);
  private:
    double m_min, m_max;
  };

  class dRLeptonJetSelectiona : public uhh2::Selection{
  public:
    explicit dRLeptonJetSelectiona(double dRmin = 0., double dRmax = -1);
    virtual bool passes(const uhh2::Event & event);
  private:
    double dRmin, dRmax;
  };

  class NGenElectronSelectiona : public uhh2::Selection{
  public:
    explicit NGenElectronSelectiona(int n_min_ = 0., int n_max_ = -1);
    virtual bool passes(const uhh2::Event & event);
  private:
    double n_min, n_max;
  };

 class NGenMuonSelectiona : public uhh2::Selection{
  public:
    explicit NGenMuonSelectiona(int n_min_ = 0., int n_max_ = -1);
    virtual bool passes(const uhh2::Event & event);
  private:
    double n_min, n_max;
  };

  class LQSemiLepMatchablea: public uhh2::Selection{
  public:
    LQSemiLepMatchablea();
    ~LQSemiLepMatchablea(){};
    virtual bool passes(const uhh2::Event & event);
  private:
  };

  class TTbarSemiLepMatchablea: public uhh2::Selection{
  public:
    TTbarSemiLepMatchablea();
    ~TTbarSemiLepMatchablea(){};
    virtual bool passes(const uhh2::Event & event);
  private:
  };

}
