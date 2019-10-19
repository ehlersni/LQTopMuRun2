#include "UHH2/LQTopMuRun2/include/LQToTopMuSelections.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/Utils.h"

#include <stdexcept>

using namespace uhh2examples;
using namespace uhh2;
using namespace std;

DijetSelectiona::DijetSelectiona(float dphi_min_, float third_frac_max_): dphi_min(dphi_min_), third_frac_max(third_frac_max_){}
bool DijetSelectiona::passes(const Event & event){
    assert(event.jets); // if this fails, it probably means jets are not read in
    if(event.jets->size() < 2) return false;
    const auto & jet0 = event.jets->at(0);
    const auto & jet1 = event.jets->at(1);
    auto dphi = deltaPhi(jet0, jet1);
    if(dphi < dphi_min) return false;
    if(event.jets->size() == 2) return true;
    const auto & jet2 = event.jets->at(2);
    auto third_jet_frac = jet2.pt() / (0.5 * (jet0.pt() + jet1.pt()));
    return third_jet_frac < third_frac_max;
}

HtSelectiona::HtSelectiona(double ht_min_, double ht_max_):ht_min(ht_min_), ht_max(ht_max_){}
bool HtSelectiona::passes(const Event & event){
  auto met = event.met->pt();

  double ht = 0.0;
  double ht_jets = 0.0;
  double ht_lep = 0.0;
  for(const auto & jet : *event.jets){
    ht_jets += jet.pt();
  }
  for(const auto & electron : *event.electrons){
    ht_lep += electron.pt();
  }
  for(const auto & muon : *event.muons){
    ht_lep += muon.pt();
  }

  ht = ht_lep + ht_jets + met;

  bool pass = false;
  pass = ht > ht_min && (ht_max < 0 || ht < ht_max);
  return pass;
}

InvMass2MuVetoa::InvMass2MuVetoa(double m_min_, double m_max_):m_min(m_min_), m_max(m_max_){}
bool InvMass2MuVetoa::passes(const Event & event){

  bool pass = true;
  int Nmuons = event.muons->size();
  double M_mumu;
  LorentzVector muons[Nmuons];
  for(int i=0; i<Nmuons; i++){
    muons[i] = event.muons->at(i).v4();
  }
  for(int i=0; i<Nmuons; i++){
    for(int j=0; j<Nmuons; j++){
      if(j > i){
	M_mumu = (muons[i] + muons[j]).M();
	if(M_mumu > m_min && M_mumu < m_max){
	  pass = false;
	}
      }
    }
  }
  return pass;
}

InvMass2MuVetoInverteda::InvMass2MuVetoInverteda(double m_min_, double m_max_):m_min(m_min_), m_max(m_max_){}
bool InvMass2MuVetoInverteda::passes(const Event & event){

  bool pass = true;
  int Nmuons = event.muons->size();
  double M_mumu;
  LorentzVector muons[Nmuons];
  for(int i=0; i<Nmuons; i++){
    muons[i] = event.muons->at(i).v4();
  }
  for(int i=0; i<Nmuons; i++){
    for(int j=0; j<Nmuons; j++){
      if(j > i){
	M_mumu = (muons[i] + muons[j]).M();
	if(M_mumu < m_min || M_mumu > m_max){
	  pass = false;
	}
      }
    }
  }
  return pass;
}

PtLeadingMuonSelectiona::PtLeadingMuonSelectiona(double pt_min_, double pt_max_):pt_min(pt_min_), pt_max(pt_max_){}
bool PtLeadingMuonSelectiona::passes(const Event & event){

  bool pass = true;
  double pt_leadingmu = event.muons->at(0).pt();
  pass = pt_leadingmu >= pt_min && (pt_leadingmu <= pt_max || pt_max < 0);
  return pass;
}

Pt2ndMuonSelectiona::Pt2ndMuonSelectiona(double pt_min_, double pt_max_):pt_min(pt_min_), pt_max(pt_max_){}
bool Pt2ndMuonSelectiona::passes(const Event & event){

  bool pass = true;
  double pt_2ndmu = event.muons->at(1).pt();
  pass = pt_2ndmu >= pt_min && (pt_2ndmu <= pt_max || pt_max < 0);
  return pass;
}

PtLeadingJetSelectiona::PtLeadingJetSelectiona(double pt_min_, double pt_max_):pt_min(pt_min_), pt_max(pt_max_){}
bool PtLeadingJetSelectiona::passes(const Event & event){

  bool pass = true;
  double pt_leadingjet = event.jets->at(0).pt();
  pass = pt_leadingjet >= pt_min && (pt_leadingjet <= pt_max || pt_max < 0);
  return pass;
}

Pt2ndJetSelectiona::Pt2ndJetSelectiona(double pt_min_, double pt_max_):pt_min(pt_min_), pt_max(pt_max_){}
bool Pt2ndJetSelectiona::passes(const Event & event){

  bool pass = true;
  double pt_2ndjet = event.jets->at(1).pt();
  pass = pt_2ndjet >= pt_min && (pt_2ndjet <= pt_max || pt_max < 0);
  return pass;
}

PtRelMuJetSelectiona::PtRelMuJetSelectiona(double ptrel_min_, double ptrel_max_) : ptrel_min(ptrel_min_), ptrel_max(ptrel_max_){}
bool PtRelMuJetSelectiona::passes(const Event & event){

  bool pass = true;
  double ptrel = 0;
  for(const auto & muon : *event.muons){
    auto nextjet = nextJet(muon, *event.jets);
    ptrel = pTrel(muon, nextjet);
    if(ptrel >= ptrel_min && (ptrel <= ptrel_max || ptrel_max < 0)){
      pass = true;
    }
    else{
      pass = false;
      return pass;
    }
  }
  return pass;
}

METSelectiona::METSelectiona(double met_min_, double met_max_) : met_min(met_min_), met_max(met_max_){}
bool METSelectiona::passes(const Event & event){

  auto met = event.met->pt();
  bool pass = true;

  pass = met >= met_min && (met <= met_max || met_max < 0);
  return pass;

}

GenLvlZMuMuSelectiona::GenLvlZMuMuSelectiona(){}
bool GenLvlZMuMuSelectiona::passes(const Event & event){

  bool pass = true;
  for (const auto & genpart : *event.genparticles){
    if (genpart.pdgId() == 23){
      if(abs(genpart.daughter(event.genparticles, 1)->pdgId()) == 13){
	pass = true;
      }
      else{
	pass = false;
      }
    }
  }
  return pass;
}

GenLvlZEESelectiona::GenLvlZEESelectiona(){}
bool GenLvlZEESelectiona::passes(const Event & event){

  bool pass = true;
  for (const auto & genpart : *event.genparticles){
    if (genpart.pdgId() == 23){
      if(abs(genpart.daughter(event.genparticles, 1)->pdgId()) == 11){
	pass = true;
      }
      else{
	pass = false;
      }
    }
  }
  return pass;
}

GenLvlTopDileptonSelectiona::GenLvlTopDileptonSelectiona(){}
bool GenLvlTopDileptonSelectiona::passes(const Event & event){

  bool pass = false;
  int N_w = 0;
  for(const auto & genpart : *event.genparticles){
    if(abs(genpart.pdgId()) == 24){
      //only W-bosons are left here
      //now check decay mode
      if(abs(genpart.daughter(event.genparticles, 1)->pdgId()) == 11 || abs(genpart.daughter(event.genparticles, 1)->pdgId()) == 13 || abs(genpart.daughter(event.genparticles, 2)->pdgId()) == 11 || abs(genpart.daughter(event.genparticles, 2)->pdgId()) == 13){
	//here, 1 of the daughters of the first W is either a electron or a muon
	//1 leptonic W is found, now allow search for second.
	if(N_w == 0){
	  N_w++;
	}

	//the following is only done, if this if-loop is entered for a second time
	else if(N_w == 1){
	  //the second W is also found to decay in the right way, now set pass = true
	  pass = true;
	  N_w++;
	}

	if(N_w != 2 && pass == true) throw std::runtime_error("in GelLvlTopDibosonSelection: Not exactly 2 W-bosons found despite pass = true");

      }//decay mode
    }//W bosons
  }//genptcls

  return pass;
}//method

PtRelMu1JetSelectiona::PtRelMu1JetSelectiona(double ptrel_min_, double ptrel_max_) : ptrel_min(ptrel_min_), ptrel_max(ptrel_max_){}
bool PtRelMu1JetSelectiona::passes(const Event & event){

  bool pass = true;
  auto muon = event.muons->at(0);
  auto nextjet = nextJet(muon, *event.jets);
  double ptrel = pTrel(muon, nextjet);

  pass = ptrel > ptrel_min && (ptrel < ptrel_max || ptrel_max < 0);
  return pass;
}

HTJetsSelectiona::HTJetsSelectiona(double ht_min_, double ht_max_) : ht_min(ht_min_), ht_max(ht_max_){}
bool HTJetsSelectiona::passes(const Event & event){

  bool pass = true;
  double ht_jets = 0.0;
  for(const auto & jet : *event.jets){
    ht_jets += jet.pt();
    }

  pass = ht_jets > ht_min && (ht_jets < ht_max || ht_max < 0);
  return pass;
}

HTLeptSelectiona::HTLeptSelectiona(double ht_min_, double ht_max_) : ht_min(ht_min_), ht_max(ht_max_){}
bool HTLeptSelectiona::passes(const Event & event){

  bool pass = true;
  double ht_lept = 0.0;
  for(const auto & ele : *event.electrons){
    ht_lept += ele.pt();
    }
  for(const auto & mu : *event.muons){
    ht_lept += mu.pt();
    }

  pass = ht_lept > ht_min && (ht_lept < ht_max || ht_max < 0);
  return pass;
}

EtaLeadingJetSelectiona::EtaLeadingJetSelectiona(double eta_max_) : eta_max(eta_max_){}
bool EtaLeadingJetSelectiona::passes(const Event & event){

  bool pass = true;
  double etajet1 = fabs(event.jets->at(0).eta());

  pass = etajet1 < eta_max || eta_max < 0;
  return pass;
}

InvMassMuEleVetoa::InvMassMuEleVetoa(double m_min_, double m_max_):m_min(m_min_), m_max(m_max_){}
bool InvMassMuEleVetoa::passes(const Event & event){

  bool pass = true;
  int Nmuons = event.muons->size();
  int Nelectrons = event.electrons->size();
  double M_muele;
  LorentzVector muons[Nmuons];
  LorentzVector electrons[Nelectrons];
  for(int i=0; i<Nmuons; i++){
    muons[i] = event.muons->at(i).v4();
  }
  for(int i=0; i<Nelectrons; i++){
    electrons[i] = event.electrons->at(i).v4();
  }
  for(int i=0; i<Nmuons; i++){
    for(int j=0; j<Nelectrons; j++){
      M_muele = (muons[i] + electrons[j]).M();
      if(M_muele > m_min && M_muele < m_max){
	pass = false;
      }
    }
  }
  return pass;
}


InvMassEleEleVetoa::InvMassEleEleVetoa(double m_min_, double m_max_):m_min(m_min_), m_max(m_max_){}
bool InvMassEleEleVetoa::passes(const Event & event){

  bool pass = true;
  int Nelectrons = event.electrons->size();
  double M_ee;
  LorentzVector electrons[Nelectrons];
  for(int i=0; i<Nelectrons; i++){
    electrons[i] = event.electrons->at(i).v4();
  }
  for(int i=0; i<Nelectrons; i++){
    for(int j=0; j<Nelectrons; j++){
      if(j>i){
	M_ee = (electrons[i] + electrons[j]).M();
	if(M_ee > m_min && M_ee < m_max){
	  pass = false;
	}
      }
    }
  }
  return pass;
}

InvMassEleEleSelectiona::InvMassEleEleSelectiona(double m_min_, double m_max_):m_min(m_min_), m_max(m_max_){}
bool InvMassEleEleSelectiona::passes(const Event & event){

  bool pass = false;
  int Nelectrons = event.electrons->size();
  double M_ee;
  LorentzVector electrons[Nelectrons];
  for(int i=0; i<Nelectrons; i++){
    electrons[i] = event.electrons->at(i).v4();
  }
  for(int i=0; i<Nelectrons; i++){
    for(int j=0; j<Nelectrons; j++){
      if(j>i){
	M_ee = (electrons[i] + electrons[j]).M();
	if(M_ee >= m_min && M_ee <= m_max){
	  pass = true;
	}
      }
    }
  }
  return pass;
}

dRLeptonJetSelectiona::dRLeptonJetSelectiona(double dRmin_, double dRmax_) : dRmin(dRmin_), dRmax(dRmax_){}
bool dRLeptonJetSelectiona::passes(const Event & event){

  bool pass = true;

  for(const auto & thismu : *event.muons){
    for(const auto & thisjet : *event.jets){
      if(!(deltaR(thismu,thisjet) >= dRmin && (deltaR(thismu,thisjet) <= dRmax || dRmax < 0))) pass = false;
    }
  }

  for(const auto & thisele : *event.electrons){
    for(const auto & thisjet : *event.jets){
      if(!(deltaR(thisele,thisjet) >= dRmin && (deltaR(thisele,thisjet) <= dRmax || dRmax < 0))) pass = false;
    }
  }

  return pass;
}

NGenElectronSelectiona::NGenElectronSelectiona(int n_min_, int n_max_) : n_min(n_min_), n_max(n_max_){}
bool NGenElectronSelectiona::passes(const Event & event){

  int n_ele = 0;
  for(const auto & gp : *event.genparticles){
    if(abs(gp.pdgId()) == 11) n_ele++;
  }

  return (n_ele >= n_min && (n_ele <= n_max || n_max < 0));
}

NGenMuonSelectiona::NGenMuonSelectiona(int n_min_, int n_max_) : n_min(n_min_), n_max(n_max_){}
bool NGenMuonSelectiona::passes(const Event & event){

  int n_mu = 0;
  for(const auto & gp : *event.genparticles){
    if(abs(gp.pdgId()) == 13) n_mu++;
  }

  return (n_mu >= n_min && (n_mu <= n_max || n_max < 0) );
}


LQSemiLepMatchablea::LQSemiLepMatchablea(){}
bool LQSemiLepMatchablea::passes(const Event & event){
  assert(event.genparticles);
  if(event.isRealData) return false;


  //check, if one LQ decays had and one decays lep
  bool found_had = false, found_lep = false;

  //Loop over genparticles
  for(const auto & gp : *event.genparticles){

    //Find LQs
    if(fabs(gp.pdgId()) == 42){
      bool top_offshell = false;
      auto t = gp.daughter(event.genparticles,1);
      auto muon = gp.daughter(event.genparticles,2);
      if(fabs(t->pdgId()) == 13 && fabs(muon->pdgId()) == 6){
	t = gp.daughter(event.genparticles,2);
	muon = gp.daughter(event.genparticles,1);
      }
      if(!(fabs(t->pdgId()) == 6 && fabs(muon->pdgId()) == 13)) top_offshell = true;

      //check if mu is matched by a reco muon
      bool matched_mu = false;
      for(const auto & mu : *event.muons){
	if(deltaR(mu,*muon) <= 0.1) matched_mu = true;
      }
      if(!matched_mu){
	//cout << "Rejected because LQ muons could not be matched." << endl;
	return false;
      }

      //Get b and W
      auto b = t->daughter(event.genparticles,1);
      auto W = t->daughter(event.genparticles,2);
      if(!top_offshell){
	b = t->daughter(event.genparticles,1);
	W = t->daughter(event.genparticles,2);
	if(fabs(W->pdgId()) == 5 && fabs(b->pdgId()) == 24){
	  b = t->daughter(event.genparticles,2);
	  W = t->daughter(event.genparticles,1);
	}
	if(!(fabs(b->pdgId()) == 5 && fabs(W->pdgId()) == 24)) return false;
      }
      else{
	if(fabs(gp.daughter(event.genparticles,1)->pdgId()) == 13){
	  b = t->daughter(event.genparticles,2);
	  W = t->daughter(event.genparticles,3);
	  if(fabs(W->pdgId()) == 5 && fabs(b->pdgId()) == 24){
	    b = t->daughter(event.genparticles,3);
	    W = t->daughter(event.genparticles,2);
	  }
	  if(!(fabs(b->pdgId()) == 5 && fabs(W->pdgId()) == 24)) return false;
	}
	else if(fabs(gp.daughter(event.genparticles,2)->pdgId()) == 13){
	  b = t->daughter(event.genparticles,1);
	  W = t->daughter(event.genparticles,3);
	  if(fabs(W->pdgId()) == 5 && fabs(b->pdgId()) == 24){
	    b = t->daughter(event.genparticles,3);
	    W = t->daughter(event.genparticles,1);
	  }
	  if(!(fabs(b->pdgId()) == 5 && fabs(W->pdgId()) == 24)) return false;
	}
	else if(fabs(gp.daughter(event.genparticles,3)->pdgId()) == 13){
	  b = t->daughter(event.genparticles,2);
	  W = t->daughter(event.genparticles,3);
	  if(fabs(W->pdgId()) == 5 && fabs(b->pdgId()) == 24){
	    b = t->daughter(event.genparticles,3);
	    W = t->daughter(event.genparticles,2);
	  }
	  if(!(fabs(b->pdgId()) == 5 && fabs(W->pdgId()) == 24)) return false;
	}
	else throw runtime_error("In LQSemiLepMatchablea: LQ does not have muon as a daughter");
      }

      //try to match the b quarks
      bool matched_b = false;
      for(const auto & jet : *event.jets){
	if(deltaR(*b,jet) <= 0.4) matched_b = true;
      }
      if(!matched_b){
	//cout << "Rejected because b could not be matched." << endl;
	return false;
      }

      //Check decaymodes of W
      auto Wd1 = W->daughter(event.genparticles,1);
      auto Wd2 = W->daughter(event.genparticles,2);

      //hadronic
      if(fabs(Wd1->pdgId()) < 7 && fabs(Wd2->pdgId()) < 7){
	if(found_had){
	  //cout << "Rejected because 2nd had top was found." << endl;
	  return false;
	}
	found_had = true;

	//check if both daughters can be matched by jets
	bool matched_d1 = false, matched_d2 = false;
	for(const auto & jet : *event.jets){
	  if(deltaR(*Wd1, jet) <= 0.4) matched_d1 = true;
	  if(deltaR(*Wd2, jet) <= 0.4) matched_d2 = true;
	}
	if(!(matched_d1 && matched_d2)){
	  //cout << "Rejected because W-quarks could not be matched." << endl;
	  return false;
	}
      }

      //leptonic
      else if((fabs(Wd1->pdgId()) == 11 || fabs(Wd1->pdgId()) == 13) || (fabs(Wd2->pdgId()) == 11 || fabs(Wd2->pdgId()) == 13)){
	if(found_lep){
	  //cout << "Rejected because 2nd lep top was found." << endl;
	  return false;
	}
	found_lep = true;

	//Find charged lepton
	auto lep = Wd1;
	auto nu = Wd2;
	if(fabs(Wd2->pdgId()) == 11 || fabs(Wd2->pdgId()) == 13){
	  lep = Wd2;
	  nu = Wd1;
	}
	if(!(fabs(lep->pdgId()) == 11 && fabs(nu->pdgId()) == 12) && !(fabs(lep->pdgId()) == 13 && fabs(nu->pdgId()) == 14)) throw runtime_error("In LQSemiLepMatchable: The leptonic W does not decay into a lepton and its neutrino.");

	//check, if lepton can be matched
	bool matched_lep = false;
	if(fabs(lep->pdgId()) == 11){
	  for(const auto & ele : *event.electrons){
	    if(deltaR(*lep,ele) <= 0.1) matched_lep = true;
	  }
	}
	else if(fabs(lep->pdgId()) == 13){
	  for(const auto & mu : *event.muons){
	    if(deltaR(mu,*lep) <= 0.1) matched_lep = true;
	  }
	}
	else throw runtime_error("In LQSemiLepMatchablea: Lepton from W decay is neither e nor mu.");
	if(!matched_lep){
	  //cout << "Rejected because add lepton could not be matched." << endl;
	  return false;
	}
      }
      //tau-decays
      else{
	//cout << "Rejected because decay mode is neither lep nor had." << endl;
	return false;
      }
    }
  }

  if(!(found_had && found_lep)){
    //cout << "Rejected because not both decay modes were found once." << endl;
    return false;
  }

  return true;
}


TTbarSemiLepMatchablea::TTbarSemiLepMatchablea(){}
bool TTbarSemiLepMatchablea::passes(const Event & event){
  assert(event.genparticles);
  if(event.isRealData) return false;

  //check, if one top decays had and one decays lep
  bool found_had = false, found_lep = false;

  //Loop over genparticles
  for(const auto & gp : *event.genparticles){

    //Get tops
    if(fabs(gp.pdgId()) == 6){

      //Get b and W
      auto b = gp.daughter(event.genparticles,1);
      auto W = gp.daughter(event.genparticles,2);
      if(fabs(W->pdgId()) == 5 && fabs(b->pdgId()) == 24){
	b = gp.daughter(event.genparticles,2);
	W = gp.daughter(event.genparticles,1);
      }
      if(!(fabs(b->pdgId()) == 5 && fabs(W->pdgId()) == 24)) return false;

      //try to match the b quarks
      bool matched_b = false;
      for(const auto & jet : *event.jets){
	if(deltaR(*b,jet) <= 0.4) matched_b = true;
      }
      if(!matched_b) return false;

      //Check decaymodes of W
      auto Wd1 = W->daughter(event.genparticles,1);
      auto Wd2 = W->daughter(event.genparticles,2);

      //hadronic
      if(fabs(Wd1->pdgId()) < 7 && fabs(Wd2->pdgId()) < 7){
	if(found_had) return false;
	found_had = true;

	//check if both daughters can be matched by jets
	bool matched_d1 = false, matched_d2 = false;
	for(const auto & jet : *event.jets){
	  if(deltaR(*Wd1, jet) <= 0.4) matched_d1 = true;
	  if(deltaR(*Wd2, jet) <= 0.4) matched_d2 = true;
	}
	if(!(matched_d1 && matched_d2)) return false;
      }

      //leptonic
      else if((fabs(Wd1->pdgId()) == 11 || fabs(Wd1->pdgId()) == 13) || (fabs(Wd2->pdgId()) == 11 || fabs(Wd2->pdgId()) == 13)){
	if(found_lep) return false;
	found_lep = true;

	//Find charged lepton
	auto lep = Wd1;
	auto nu = Wd2;
	if(fabs(Wd2->pdgId()) == 11 || fabs(Wd2->pdgId()) == 13){
	  lep = Wd2;
	  nu = Wd1;
	}
	if(!(fabs(lep->pdgId()) == 11 && fabs(nu->pdgId()) == 12) && !(fabs(lep->pdgId()) == 13 && fabs(nu->pdgId()) == 14)) throw runtime_error("In TTbarSemiLepMatchable: The leptonic W does not decay into a lepton and its neutrino.");

	//check, if lepton can be matched
	bool matched_lep = false;
	if(fabs(lep->pdgId()) == 11){
	  for(const auto & ele : *event.electrons){
	    if(deltaR(*lep,ele) <= 0.1) matched_lep = true;
	  }
	}
	else if(fabs(lep->pdgId()) == 13){
	  for(const auto & mu : *event.muons){
	    if(deltaR(mu,*lep) <= 0.1) matched_lep = true;
	  }
	}
	else throw runtime_error("In TTbarSemiLepMatchablea: Lepton from W decay is neither e nor mu.");
	if(!matched_lep) return false;
      }
      //tau-decays
      else return false;
    }
  }

  if(!(found_had && found_lep)) return false;

  return true;
}
