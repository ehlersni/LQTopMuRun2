#include "UHH2/LQTopMuRun2/include/LQReconstruction.h"
#include "UHH2/core/include/LorentzVector.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/core/include/Utils.h"

#include <cassert>

using namespace uhh2;
using namespace std;

LQPrimaryLepton::LQPrimaryLepton(Context & ctx) {
    h_primlep = ctx.get_handle<FlavorParticle>("LQPrimaryLepton");
}

bool LQPrimaryLepton::process(uhh2::Event & event) {
  assert(/*event.muons || */event.electrons);
    double ptmax = -infinity;
    FlavorParticle primlep;
    if(event.electrons) {
        for(const auto & ele : *event.electrons) {
            if(ele.pt() > ptmax) {
                ptmax = ele.pt();
                primlep = ele;
            }
        }
    }

    event.set(h_primlep, std::move(primlep));
    return true;
}

LQPrimaryLepton::~LQPrimaryLepton() {}


HighMassInclusiveLQReconstruction::HighMassInclusiveLQReconstruction(Context & ctx, const NeutrinoReconstructionMethod & neutrinofunction): m_neutrinofunction(neutrinofunction) {
  h_recohyps = ctx.get_handle<vector<LQReconstructionHypothesis>>("LQHypotheses");
  h_is_mlq_reconstructed = ctx.get_handle<bool>("is_mlq_reconstructed");
  h_mlq_reco_mode = ctx.get_handle<TString>("mlq_reco_mode");
}

HighMassInclusiveLQReconstruction::~HighMassInclusiveLQReconstruction() {}

bool HighMassInclusiveLQReconstruction::process(uhh2::Event & event) {
  assert(event.jets);
  assert(event.met);
  int n_final_hyps = 0;

  // Requirements:
  // 1) >=2 muons, >=1 electron -->Take hardest electron
  // 2) >=3 muons, ==0 electrons -->Take 3 hardest muons
  TString mode = event.get(h_mlq_reco_mode);
  if(mode != "ele" && mode != "muon") return false;


  std::vector<LQReconstructionHypothesis> recoHyps;
  if(mode == "ele"){

    const Particle & electron = event.electrons->at(0); //always an electron, as only electrons are considered possible primary leptons
    LorentzVector electron_v4 = electron.v4();

    //reconstruct neutrino
    std::vector<LorentzVector> neutrinos = m_neutrinofunction( electron.v4(), event.met->v4());

    unsigned int n_muons = event.muons->size();
    unsigned int n_jets = event.jets->size();
    if(n_jets>7) n_jets=7; //avoid crashes in events with many jets
    // idea: loop over 3^Njet possibilities and write the current loop
    // index j in the 3-base system. The Njets digits represent whether
    // to assign each jet to the hadronic side (0), leptonic side (1),
    // or none of them (2).
    const unsigned int max_j = pow(3, n_jets);

    //loop over neutrino solutions and jet assignments to fill hyotheses
    for(const auto & neutrino_p4 : neutrinos) {
      const LorentzVector wlep_v4 = electron.v4() + neutrino_p4;
      for (unsigned int j=0; j < max_j; j++) {
        LorentzVector tophad_v4;
        LorentzVector toplep_v4 = wlep_v4;
        int hadjets=0;
        int lepjets=0;
        int num = j;
        LQReconstructionHypothesis hyp;
        vector<Particle> had_jets, lep_jets;
        for (unsigned int k=0; k<n_jets; k++) {
          if(num%3==0) {
            tophad_v4 = tophad_v4 + event.jets->at(k).v4();
            had_jets.push_back(event.jets->at(k));
            hadjets++;
          }

          if(num%3==1) {
            toplep_v4 = toplep_v4 + event.jets->at(k).v4();
            lep_jets.push_back(event.jets->at(k));
            lepjets++;
          }
          //in case num%3==2 do not take this jet at all
          //shift the trigits of num to the right:
          num /= 3;
        }

        hyp.set_tophad_jets(had_jets);
        hyp.set_toplep_jets(lep_jets);

        //search jet with highest pt assigned to leptonic top
        int blep_idx(-1);
        float maxpt(-1.);
        for(unsigned int i=0; i<hyp.toplep_jets().size(); ++i){
          if(maxpt< hyp.toplep_jets().at(i).pt()){
            maxpt = hyp.toplep_jets().at(i).pt();
            blep_idx = i;
          }
        }
        if(blep_idx != -1) hyp.set_blep_v4(hyp.toplep_jets().at(blep_idx).v4());

        //fill only hypotheses with at least one jet assigned to each top quark
        if(hadjets>0 && lepjets>0) {
          int max_i = pow(3,n_muons); // analogous to jet combinations
          for(int i=0; i<max_i; i++){ // for each jet comb loop over all possible muon combs
            LorentzVector mu1_v4;
            LorentzVector mu2_v4;
            int hadmu=0;
            int lepmu=0;
            int num = i;
            for(unsigned int k=0; k<n_muons; k++){
              if(num%3==0){
                mu1_v4 = event.muons->at(k).v4();
                hyp.set_mu_had(event.muons->at(k)); // had is only a hypothesis not having considered the charge!!!
                hadmu++;
              }
              if(num%3==1){
                mu2_v4 = event.muons->at(k).v4();
                hyp.set_mu_lep(event.muons->at(k)); // lep is only a hypothesis not having considered the charge!!!
                lepmu++;
              }
              //for num%3==2 do nothing
              num /= 3;
            }

            Particle mu_2 = hyp.mu_lep();
            Particle mu_1 = hyp.mu_had();
            if(hadmu==1 && lepmu==1 && (mu_1.charge()!=mu_2.charge())){ //require exactly 1 muon assigned to each top, opposite charges
              if(mu_2.charge() != electron.charge()){
                hyp.set_mu_had_v4(mu1_v4);
                hyp.set_mu_lep_v4(mu2_v4); // mu2 really is the leptonic one.
                hyp.set_mu_had(hyp.mu_had());
                hyp.set_mu_lep(hyp.mu_lep());
                hyp.set_tophad_v4(tophad_v4);
                hyp.set_toplep_v4(toplep_v4);
                hyp.set_electron(electron);
                hyp.set_electron_v4(electron_v4);
                hyp.set_neutrino_v4(neutrino_p4);
                hyp.set_tophad_jets(had_jets);
                hyp.set_toplep_jets(lep_jets);
              } // charge comparison
              else{
                hyp.set_mu_had_v4(mu2_v4);
                hyp.set_mu_lep_v4(mu1_v4); // mu1 really is the leptonic one, the original hypothesis was wrong.
                hyp.set_mu_had(hyp.mu_lep());
                hyp.set_mu_lep(hyp.mu_had());
                hyp.set_tophad_v4(tophad_v4);
                hyp.set_toplep_v4(toplep_v4);
                hyp.set_electron(electron);
                hyp.set_electron_v4(electron_v4);
                hyp.set_neutrino_v4(neutrino_p4);
                hyp.set_tophad_jets(had_jets);
                hyp.set_toplep_jets(lep_jets);
              } // charge comparison_2

              recoHyps.emplace_back(move(hyp));
              n_final_hyps++;

            } // 1 muon per top
          } // muon combs for-loop
        } // if at least 1 jet is assigned to each top quark
      } // 3^n_jets jet combinations * n_muon muon combinations
    } // neutrinos



  } //ele case
  else if(mode == "muon"){

    //sum should be either +1 or -1, if there is a pair with opposite charges as required by the reco
    vector<Muon> recomuons;
    int idx=0;
    for(const auto & mu : *event.muons){
      if(idx<3){ //Only consider 3 leading muons
        recomuons.emplace_back(mu);
        idx++;
      }
    }

    double sum_charges = 0;
    int inx = 0;
    for(const auto & mu : recomuons){
      sum_charges += mu.charge();
    }
    if(fabs(sum_charges) != 1) throw runtime_error("In MuonicLQReconstruction: All muons have the same charge. Should not be the case.");

    //find muons with same charges, one of these must come from W if our assumption is correct. Must be ==2 candidates
    vector<Particle> muon_candidates;
    vector<bool> used_as_candidate;
    inx=0;
    for(const auto & mu : recomuons){
      bool used = false;
      if(sum_charges > 0){
        if(mu.charge() > 0){
          muon_candidates.push_back(mu);
          used = true;
        }
      }
      else{
        if(mu.charge() < 0){
          muon_candidates.push_back(mu);
          used = true;
        }
      }
      used_as_candidate.push_back(used);
      inx++;
    }
    if(muon_candidates.size() != 2) throw runtime_error("In MuonicLQReconstruction: Not ==2 muon candidates. Should have been prohibited by earlier runtime_error-conditions...");


    //loop over possible muon candidates
    for(int m=0; m<2; m++){

      unsigned int current_muon = 99;
      const Particle muon = muon_candidates.at(m);
      LorentzVector muon_v4 = muon.v4();
      if(m == 0){
        if(used_as_candidate[0]) current_muon = 0;
        else if(used_as_candidate[1]) current_muon = 1;
        else throw runtime_error("In MuonicLQReconstruction: Neither the first nor the second entry of used_as_candidate is 'true'. There are 2 candidates out of 3 muons, two entries should be true, which is not possible anymore.");
      }
      else if(m==1){
        if(used_as_candidate[0] && used_as_candidate[1]) current_muon = 1;
        else if(used_as_candidate[0] && !used_as_candidate[1] && used_as_candidate[2]) current_muon = 2;
        else if(!used_as_candidate[0] && used_as_candidate[1] && used_as_candidate[2]) current_muon = 2;
        else throw runtime_error("In MuonicLQReconstruction: Couldn't find second muon-candidate.");
      }

      //reconstruct neutrino
      std::vector<LorentzVector> neutrinos = m_neutrinofunction( muon.v4(), event.met->v4());

      unsigned int n_muons = recomuons.size();
      unsigned int n_jets = event.jets->size();
      if(n_jets>7) n_jets=7; //avoid crashes in events with many jets
      // idea: loop over 3^Njet possibilities and write the current loop
      // index j in the 3-base system. The Njets digits represent whether
      // to assign each jet to the hadronic side (0), leptonic side (1),
      // or none of them (2).
      const unsigned int max_j = pow(3, n_jets);

      //loop over neutrino solutions and jet assignments to fill hyotheses
      for(const auto & neutrino_p4 : neutrinos) {
        const LorentzVector wlep_v4 = muon.v4() + neutrino_p4;
        for (unsigned int j=0; j < max_j; j++) {
          LorentzVector tophad_v4;
          LorentzVector toplep_v4 = wlep_v4;
          int hadjets=0;
          int lepjets=0;
          int num = j;
          LQReconstructionHypothesis hyp;
          //hyp.set_muon(muon);
          //hyp.set_muon_v4(muon_v4);
          //hyp.set_neutrino_v4(neutrino_p4);
          vector<Particle> had_jets, lep_jets;
          for (unsigned int k=0; k<n_jets; k++) {
            if(num%3==0) {
              tophad_v4 = tophad_v4 + event.jets->at(k).v4();
              //hyp.add_tophad_jet(event.jets->at(k));
              had_jets.push_back(event.jets->at(k));
              hadjets++;
            }

            if(num%3==1) {
              toplep_v4 = toplep_v4 + event.jets->at(k).v4();
              //hyp.add_toplep_jet(event.jets->at(k));
              lep_jets.push_back(event.jets->at(k));
              lepjets++;
            }
            //in case num%3==2 do not take this jet at all
            //shift the trigits of num to the right:
            num /= 3;
          }

          hyp.set_tophad_jets(had_jets);
          hyp.set_toplep_jets(lep_jets);

          //search jet with highest pt assigned to leptonic top
          int blep_idx(-1);
          float maxpt(-1.);
          for(unsigned int i=0; i<hyp.toplep_jets().size(); ++i){
            if(maxpt< hyp.toplep_jets().at(i).pt()){
              maxpt = hyp.toplep_jets().at(i).pt();
              blep_idx = i;
            }
          }
          if(blep_idx != -1) hyp.set_blep_v4(hyp.toplep_jets().at(blep_idx).v4());

          //fill only hypotheses with at least one jet assigned to each top quark
          if(hadjets>0 && lepjets>0) {
            int max_i = pow(3,n_muons-1); // analogous to jet combinations
            for(int i=0; i<max_i; i++){ // for each jet comb loop over all possible muon combs
              LorentzVector mu1_v4;
              LorentzVector mu2_v4;
              int hadmu=0;
              int lepmu=0;
              int num = i;
              for(unsigned int k=0; k<n_muons; k++){
                if(k == current_muon) continue;
                if(num%3==0){
                  mu1_v4 = recomuons.at(k).v4();
                  hyp.set_mu_had(recomuons.at(k)); // had is only a hypothesis not having considered the charge!!!
                  hadmu++;
                }
                if(num%3==1){
                  mu2_v4 = recomuons.at(k).v4();
                  hyp.set_mu_lep(recomuons.at(k)); // lep is only a hypothesis not having considered the charge!!!
                  lepmu++;
                }
                //for num%3==2 do nothing
                num /= 3;
              }

              Particle mu_2 = hyp.mu_lep();
              Particle mu_1 = hyp.mu_had();
              if(hadmu==1 && lepmu==1 && (mu_1.charge()!=mu_2.charge())){ //require exactly 1 muon assigned to each top, opposite charges


                if(mu_2.charge() != muon.charge()){
                  hyp.set_mu_had_v4(mu1_v4);
                  hyp.set_mu_lep_v4(mu2_v4); // mu2 really is the leptonic one.
                  hyp.set_mu_had(hyp.mu_had());
                  hyp.set_mu_lep(hyp.mu_lep());
                  hyp.set_tophad_v4(tophad_v4);
                  hyp.set_toplep_v4(toplep_v4);
                  hyp.set_muon(muon);
                  hyp.set_muon_v4(muon_v4);
                  hyp.set_neutrino_v4(neutrino_p4);
                  hyp.set_tophad_jets(had_jets);
                  hyp.set_toplep_jets(lep_jets);
                } // charge comparison
                else{
                  hyp.set_mu_had_v4(mu2_v4);
                  hyp.set_mu_lep_v4(mu1_v4); // mu1 really is the leptonic one, the original hypothesis was wrong.
                  hyp.set_mu_had(hyp.mu_lep());
                  hyp.set_mu_lep(hyp.mu_had());
                  hyp.set_tophad_v4(tophad_v4);
                  hyp.set_toplep_v4(toplep_v4);
                  hyp.set_muon(muon);
                  hyp.set_muon_v4(muon_v4);
                  hyp.set_neutrino_v4(neutrino_p4);
                  hyp.set_tophad_jets(had_jets);
                  hyp.set_toplep_jets(lep_jets);
                } // charge comparison_2

                recoHyps.emplace_back(move(hyp));
                n_final_hyps++;


              } // 1 muon per top
            } // muon combs for-loop
          } // if at least 1 jet is assigned to each top quark
        } // 3^n_jets jet combinations * n_muon muon combinations
      } // neutrinos
    } //muon candidates



  } //is_muon_case
  else throw runtime_error("In InclusiveLQReconstruction: Neither muon nor electron case. Should have been prohibited by earlier runtime_error-conditions.");
  // cout << "recoHyps.size(): " << recoHyps.size() << endl;

  event.set(h_recohyps, move(recoHyps));
  event.set(h_is_mlq_reconstructed, true);


  return true;
}






std::vector<LorentzVector> LQNeutrinoReconstruction(const LorentzVector & lepton, const LorentzVector & met) {
  TVector3 lepton_pT = toVector(lepton);
  lepton_pT.SetZ(0);
  TVector3 neutrino_pT = toVector(met);
  neutrino_pT.SetZ(0);
  constexpr float mass_w = 80.399f;
  float mu = mass_w * mass_w / 2 + lepton_pT * neutrino_pT;
  float A = - (lepton_pT * lepton_pT);
  float B = mu * lepton.pz();
  float C = mu * mu - lepton.e() * lepton.e() * (neutrino_pT * neutrino_pT);
  float discriminant = B * B - A * C;
  std::vector<LorentzVector> solutions;
  if (0 >= discriminant) {
    // Take only real part of the solution for pz:
    LorentzVectorXYZE solution (met.Px(),met.Py(),-B / A,0);
    solution.SetE(solution.P());
    solutions.emplace_back(toPtEtaPhi(solution));
  }
  else {
    discriminant = sqrt(discriminant);
    LorentzVectorXYZE solution (met.Px(),met.Py(),(-B - discriminant) / A,0);
    solution.SetE(solution.P());
    solutions.emplace_back(toPtEtaPhi(solution));

    LorentzVectorXYZE solution2 (met.Px(),met.Py(),(-B + discriminant) / A,0);
    solution2.SetE(solution2.P());
    solutions.emplace_back(toPtEtaPhi(solution2));
  }
  return solutions;
}
