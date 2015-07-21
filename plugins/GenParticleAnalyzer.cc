// -*- C++ -*-
//
// Package:    GenProduction/GenParticleAnalyzer
// Class:      GenParticleAnalyzer
// 
/**\class GenParticleAnalyzer GenParticleAnalyzer.cc GenProduction/GenParticleAnalyzer/plugins/GenParticleAnalyzer.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Devdatta Majumder
//         Created:  Thu, 16 Apr 2015 12:47:35 GMT
//
//


// system include files
#include <memory>
#include <map>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h" 

#include <TH1D.h>
#include <TH2D.h>
#include <TLorentzVector.h>

//
// class declaration
//

class GenParticleAnalyzer : public edm::EDAnalyzer {
public:
  explicit GenParticleAnalyzer(const edm::ParameterSet&);
  ~GenParticleAnalyzer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  
  // ----------member data ---------------------------
  edm::Service<TFileService> fs             ; 
  std::map<std::string, TH1D*> h1_          ; 
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
GenParticleAnalyzer::GenParticleAnalyzer(const edm::ParameterSet& iConfig)

{
  //now do what ever initialization is needed

}


GenParticleAnalyzer::~GenParticleAnalyzer()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  
}


//
/// member functions
//

// ------------ method called for each event  ------------
void GenParticleAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace reco;
  using namespace edm;
  using namespace std;
  
  Handle<GenParticleCollection> genParticles;
  iEvent.getByLabel("genParticles", genParticles);
  //cout << "Size = " << genParticles->size() << endl; 
  for(size_t i = 0; i < genParticles->size(); ++ i) { 
    const GenParticle & part = (*genParticles)[i];
    int id = part.pdgId(), st = part.status() ;
    //cout << "PDG Id " << id << " Status " << st << endl;
    int n_stop = 0;
    if ( abs(id) == 1000006  ) {
        n_stop++;
    }
    if ( n_stop > 0 ) cout << " n_stop " << n_stop << endl;
    if ( abs(id) == 1000006  && st == 62 ) {
      TLorentzVector p4l, p4n ; 
      //      const Candidate* mother = gdau->daughter(iggdau) ; 
      //p4l.SetPtEtaPhiM(ggdau->pt(), ggdau->eta(), ggdau->phi(), ggdau->mass()) ; 
      h1_["stop_mass"] -> Fill(part.mass());
      h1_["stop_pt"] -> Fill(part.pt());
      h1_["stop_eta"] -> Fill(part.eta());
     cout << " Stop mass " << part.mass() << " GenPart size " << genParticles->size() << endl;
      unsigned ndau = part.numberOfDaughters() ; 
      cout << " ndau " << ndau << endl;
      if (ndau != 3) continue;
      //if (ndau != 2) continue;
      int dau0 = part.daughter(0)->pdgId();
      int dau1 = part.daughter(1)->pdgId();
      int dau2 = part.daughter(2)->pdgId();	
      int daust0 = part.daughter(0)->status();
      int daust1 = part.daughter(1)->status();
      int daust2 = part.daughter(2)->status();
      
      if ( abs(dau0) == abs(dau1) ) continue; //final state stable particle

      cout << " Id " << id << " Status " << st << " dau0Id " << dau0 << " dau0Status " << daust0 << " dau1Id " << dau1 << " dau1Status " << daust1 << " dau2Id " << dau2 << " dau2Status " << daust2 << endl; 
      //cout << " Id " << id << " Status " << st << " dau0Id " << dau0 << " dau0Status " << daust0 << " dau1Id " << dau1 << " dau1Status " << daust1 << endl; 
      
      if ( abs(dau0) == 1000022 ) {//Neutralino 1
	h1_["N1_mass"] -> Fill(part.daughter(0)->mass());
	h1_["N1_pt"] -> Fill(part.daughter(0)->pt());
	h1_["N1_eta"] -> Fill(part.daughter(0)->eta());
      }

      if ( abs(dau1) == 6 ) {//t-quark
        h1_["t_mass"] -> Fill(part.daughter(1)->mass());
        h1_["t_pt"] -> Fill(part.daughter(1)->pt());
        h1_["t_eta"] -> Fill(part.daughter(1)->eta());
      }

      if ( abs(dau1) == 5 ) {//b-quark
	h1_["b_mass"] -> Fill(part.daughter(1)->mass());
	h1_["b_pt"] -> Fill(part.daughter(1)->pt());
	h1_["b_eta"] -> Fill(part.daughter(1)->eta());
        int b_jet_mult = 0 ;
        if ( part.daughter(1)->pt() > 20 && abs(part.daughter(1)->eta()) < 2.4 ) {
               b_jet_mult++;
        }
        h1_["n_b_jet"] -> Fill(b_jet_mult) ;
      }

      if ( abs(dau2) == 24 ) {//W
	h1_["W_mass"] -> Fill(part.daughter(2)->mass());
	h1_["W_pt"] -> Fill(part.daughter(2)->pt());
	h1_["W_eta"] -> Fill(part.daughter(2)->eta());
      }

      if ( abs(dau2) == 24 && abs(dau1) == 5 ) {
        
        cout << " DeltaR " << reco::deltaR(*part.daughter(1), *part.daughter(2)) << endl ;
        h1_["delta_R_W_b"] -> Fill(reco::deltaR(*part.daughter(1), *part.daughter(2))); 
      } else { cout << " Other t-decays than Wb " << endl;}


    }//stop loop


  }//GenPertticle loop

  //  Handle<reco::GenMETCollection> genMetCalo;
  //  iEvent.getByLabel("genMetCalo", genMetCalo);
  
  Handle<reco::GenJetCollection> h_ak4GenJets;
  iEvent.getByLabel("ak4GenJets", h_ak4GenJets);
  
  // cout << "Size ak4GenJets = " << (h_ak4GenJets.product())->size() << endl; 
  int njet = 0 ;
  std::vector<reco::GenJet>::const_iterator ijet ; 
  for ( ijet = (h_ak4GenJets.product())->begin(); ijet != (h_ak4GenJets.product())->end(); ++ijet ) { 

    std::vector<const GenParticle*> genjetconsts = ijet->getGenConstituents() ; 
    //  std::cout << " genjet constituents " << genjetconsts.size() << std::endl ; 
    
    for(size_t ij = 0; ij < genjetconsts.size(); ++ ij) {
      //const GenParticle & part = (*genParticles)[i];
      //cout << " ij " << ij << endl; 
    }


    if (ijet->pt() <= 40 ) continue;
    if ( njet == 0 ){

      // cout << " Jet Pt0 " << ijet->pt() << " Jet Eta0 " << ijet->eta() << " Jet Phi0 " << ijet->phi() << endl ;
      h1_["jet0_pt"] -> Fill( ijet->pt() ) ;
      h1_["jet0_eta"] -> Fill( ijet->eta() ) ;
      h1_["jet0_phi"] -> Fill( ijet->phi() ) ;
    }
    
    if ( njet == 1 ){ 
      // cout << " Jet Pt1 " << ijet->pt() << " Jet Eta1 " << ijet->eta() << " Jet Phi1 " << ijet->phi() << endl ;
      h1_["jet1_pt"] -> Fill( ijet->pt() ) ;
      h1_["jet1_eta"] -> Fill( ijet->eta() ) ;
      h1_["jet1_phi"] -> Fill( ijet->phi() ) ;
    }
    
    h1_["jet_pt"] -> Fill( ijet->pt() ) ;
    h1_["jet_eta"] -> Fill( ijet->eta() ) ;
    h1_["jet_phi"] -> Fill( ijet->phi() ) ;
    //cout << " Jet Pt " << ijet->pt() << " Jet Eta " << ijet->eta() << " Jet Phi " << ijet->phi() << endl ;
    njet++;
  }
  // cout << " Njets = " << njet << endl ; 
  
  h1_["njets"] -> Fill( njet ) ;
  
  
  
}//end of GenpartAn


// ------------ method called once each job just before starting event loop  ------------
void 
GenParticleAnalyzer::beginJob()
{
  TFileDirectory results = TFileDirectory( fs->mkdir("results") );
  h1_["stop_mass"] = fs->make<TH1D>("stop_mass", "Stop Mass", 1000, 0., 1000.) ;
  h1_["stop_pt"] = fs->make<TH1D>("stop_pt", "Stop Pt", 100, 0., 1600.) ;
  h1_["stop_eta"] = fs->make<TH1D>("stop_eta", "Stop Eta", 50, -5.0, 5.0) ;
  h1_["N1_mass"] = fs->make<TH1D>("N1_mass", "N1 Mass", 100, 0., 1600.) ;
  h1_["N1_pt"] = fs->make<TH1D>("N1_pt", "N1 Pt", 100, 0., 1600.) ;
  h1_["N1_eta"] = fs->make<TH1D>("N1_eta", "N1 Eta", 50, -5.0, 5.0) ;
  h1_["t_mass"] = fs->make<TH1D>("t_mass", "t Mass", 500, 0., 500.) ;
  h1_["t_pt"] = fs->make<TH1D>("t_pt", "t Pt", 100, 0., 400.) ;
  h1_["t_eta"] = fs->make<TH1D>("t_eta", "t Eta", 50, -5.0, 5.0) ;
  h1_["W_mass"] = fs->make<TH1D>("W_mass", "W Mass", 500, 0., 500.) ;
  h1_["W_pt"] = fs->make<TH1D>("W_pt", "W Pt", 50, 0., 800.) ;
  h1_["W_eta"] = fs->make<TH1D>("W_eta", "W Eta", 50, -5.0, 5.0) ;
  h1_["b_mass"] = fs->make<TH1D>("b_mass", "b Mass", 10, 0., 10.) ;
  h1_["b_pt"] = fs->make<TH1D>("b_pt", "b Pt", 100, 0., 400.) ;
  h1_["b_eta"] = fs->make<TH1D>("b_eta", "b Eta", 50, -5.0, 5.0) ;
  h1_["n_b_jet"] = fs->make<TH1D>("n_b_jets", "No. of central b-jets (gen)", 21, -0.5, 20.5) ;
  h1_["delta_R_W_b"] = fs->make<TH1D>("delta_R_W_b", "dR of W,b", 100, 0.0, 10.0) ; 

  h1_["jet_pt"] = fs->make<TH1D>("jet_pt", "All jets Pt > 40", 100, 0., 800.) ;
  h1_["jet_eta"] = fs->make<TH1D>("jet_eta", "Jet Eta", 50, -5.0, 5.0) ;
  h1_["jet_phi"] = fs->make<TH1D>("jet_phi", "Jet Phi", 50, -3.5, 3.5) ;
  h1_["njets"] = fs->make<TH1D>("njets", "No. of Jets", 21, -0.5, 20.5) ;
  
  h1_["jet0_pt"] = fs->make<TH1D>("jet0_pt", "Leading Jet Pt", 100, 0., 800.) ;
  h1_["jet0_eta"] = fs->make<TH1D>("jet0_eta", "Leading Jet Eta", 50, -5.0, 5.0) ;
  h1_["jet0_phi"] = fs->make<TH1D>("jet0_phi", "Leading Jet Phi", 50, -3.5, 3.5) ;
  
  h1_["jet1_pt"] = fs->make<TH1D>("jet1_pt", "Sub-leading Jet Pt", 100, 0., 800.) ;
  h1_["jet1_eta"] = fs->make<TH1D>("jet1_eta", "Sub-leading Jet Eta", 50, -5.0, 5.0) ;
  h1_["jet1_phi"] = fs->make<TH1D>("jet1_phi", "Sub-leading Jet Phi", 50, -3.5, 3.5) ;
  
  //h1_["ptWfromtgen"] = fs->make<TH1D>("ptWfromtgen", ";p_{T}(W);;", 200, 0., 1000.) ;  
  //h1_["MtWfromtgen"] = fs->make<TH1D>("MtWfromtgen", ";M_{T}(W);;", 200, 0., 1000.) ;  
  //h1_["MWfromtgen"] = fs->make<TH1D>("MWfromtgen", ";M(W);;", 100, 0., 100.) ;  
  //h1_["Mlnugen"] = fs->make<TH1D>("Mlnugen", ";M(l#nu);;", 100, 0., 100.) ;  
  //h1_["ptlnugen"] = fs->make<TH1D>("ptlnugen", ";p_{T}(l#nu);;", 200, 0., 1000.) ;  
  //h1_["Mtlnugen"] = fs->make<TH1D>("Mtlnugen", ";M_{T}(l#nu);;", 50, 0., 200.) ;  
  //h1_["ptlgen"] = fs->make<TH1D>("ptlgen", ";p_{T}(l) [GeV]", 200, 0., 1000.) ; 
  //h1_["ptnugen"] = fs->make<TH1D>("ptnugen", ";p_{T}(#nu) [GeV]", 20, 0., 200.) ; 
  //h1_["MWfromHgen"] = fs->make<TH1D>("MWfromHgen", ";M(W);;", 100, 0., 100.) ;  
  //h1_["Mlnust1gen"] = fs->make<TH1D>("Mlnust1gen", ";M(l#nu);;", 100, 0., 100.) ;  
  
}

// ------------ method called once each job just after ending the event loop  ------------
void 
GenParticleAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
  void 
  GenParticleAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when ending the processing of a run  ------------
/*
  void 
  GenParticleAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
  void 
  GenParticleAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
  void 
  GenParticleAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GenParticleAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenParticleAnalyzer);
