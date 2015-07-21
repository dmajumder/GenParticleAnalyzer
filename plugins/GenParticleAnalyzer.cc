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


#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
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
// member functions
//

// ------------ method called for each event  ------------
  void
GenParticleAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace reco;
  using namespace edm;
  using namespace std;

  Handle<GenParticleCollection> genParticles;
  iEvent.getByLabel("genParticles", genParticles);

  for(size_t i = 0; i < genParticles->size(); ++ i) { 
    const GenParticle & part = (*genParticles)[i];
    int id = part.pdgId(), st = part.status() ; 
    if ( abs(id) == 6 && ( abs(st) >=21 && abs(st) <=29) ) {
      TLorentzVector p4l, p4n ; 
      unsigned ndau = part.numberOfDaughters() ; 
      for ( size_t idau = 0; idau < ndau; ++idau) {
        const Candidate* dau = part.daughter(idau) ; 
        int dauid = dau->pdgId() ; 
        if ( abs(dauid) != 24 ) continue ; 
        h1_["MWfromtgen"] -> Fill(dau->mass()) ; 
        TLorentzVector p4w ;
        p4w.SetPtEtaPhiM(dau->pt(), dau->eta(), dau->phi(), dau->mass()) ;
        h1_["MtWfromtgen"] -> Fill(dau->mass()) ; 
        h1_["ptWfromtgen"] -> Fill(p4w.Mt()) ; 
        unsigned ngdau = dau->numberOfDaughters() ; 
        if ( ngdau == 1 && abs(dau->daughter(0)->pdgId()) == 24 ) {
          const Candidate* gdau = dau->daughter(0) ; 
          unsigned nggdau = gdau->numberOfDaughters() ; 
          if ( nggdau != 2 ) continue ; 
          for ( size_t iggdau = 0; iggdau < nggdau; ++iggdau ) {
            const Candidate* ggdau = gdau->daughter(iggdau) ; 
            if ( abs(ggdau->pdgId()) == 11 || abs(ggdau->pdgId()) == 13 || abs(ggdau->pdgId()) == 15 ) 
              p4l.SetPtEtaPhiM(ggdau->pt(), ggdau->eta(), ggdau->phi(), ggdau->mass()) ; 
            if ( abs(ggdau->pdgId()) == 12 || abs(ggdau->pdgId()) == 14 || abs(ggdau->pdgId()) == 16 ) 
              p4n.SetPtEtaPhiM(ggdau->pt(), ggdau->eta(), ggdau->phi(), ggdau->mass()) ; 
          } //// Loop over W daughters 
        } //// W -> W
        else {
          for ( size_t igdau = 0; igdau < ngdau; ++igdau ) {
            const Candidate* gdau = dau->daughter(igdau) ; 
            if ( abs(gdau->pdgId()) == 11 || abs(gdau->pdgId()) == 13 || abs(gdau->pdgId()) == 15 ) 
              p4l.SetPtEtaPhiM(gdau->pt(), gdau->eta(), gdau->phi(), gdau->mass()) ; 
            if ( abs(gdau->pdgId()) == 12 || abs(gdau->pdgId()) == 14 || abs(gdau->pdgId()) == 16 ) 
              p4n.SetPtEtaPhiM(gdau->pt(), gdau->eta(), gdau->phi(), gdau->mass()) ; 
          }
        } //// W -> ff
        if ( p4l.Pt() > 0 && p4n.Pt() > 0  ) {
          h1_["ptlgen"]->Fill(p4l.Pt()) ; 
          h1_["ptnugen"]->Fill(p4n.Pt()) ; 
          h1_["etalgen"]->Fill(p4l.Eta()) ; 
          h1_["etanugen"]->Fill(p4n.Eta()) ; 
          h1_["ptlnugen"]->Fill( (p4l+p4n).Pt() ) ; 
          h1_["Mlnugen"]->Fill( (p4l+p4n).Mag() ) ; 
          h1_["Mtlnugen"]->Fill( sqrt(p4l.Pt()*p4n.Pt()*TMath::Cos( abs(p4l.Theta() - p4n.Theta()) )) ) ;  
        }
      } //// looping over top daughters
    } //// top quark found 
  } //// particle loop


    Handle<reco::GenJetCollection> h_ak4GenJets;
    iEvent.getByLabel("ak4GenJets", h_ak4GenJets);
    std::vector<reco::GenJet>::const_iterator ijet ;
    for ( ijet = (h_ak4GenJets.product())->begin(); ijet != (h_ak4GenJets.product())->end(); ++ijet ) { 

      bool isBjet(false); 
      int nBHadrons(0) ;
      int nDHadrons(0) ; 
      std::vector<const GenParticle*> genjetconsts = ijet->getGenConstituents() ; 
      for (auto part : genjetconsts) {
        int partid = part->pdgId() ; 
        if ( (partid/100)%10 == 5 || (partid/1000)%10 == 5 ) ++nBHadrons ; 
        if (  (partid/100)%10 == 4 || (partid/1000)%10 == 4 ) ++nDHadrons ; 
      }
      if (nBHadrons > 0) isBjet = true ; 
    }


}


// ------------ method called once each job just before starting event loop  ------------
  void 
GenParticleAnalyzer::beginJob()
{
  TFileDirectory results = TFileDirectory( fs->mkdir("results") );
  h1_["ptWfromtgen"] = fs->make<TH1D>("ptWfromtgen", ";p_{T}(W);;", 200, 0., 1000.) ;  
  h1_["MtWfromtgen"] = fs->make<TH1D>("MtWfromtgen", ";M_{T}(W);;", 200, 0., 1000.) ;  
  h1_["MWfromtgen"] = fs->make<TH1D>("MWfromtgen", ";M(W);;", 100, 0., 100.) ;  
  h1_["Mlnugen"] = fs->make<TH1D>("Mlnugen", ";M(l#nu);;", 100, 0., 100.) ;  
  h1_["ptlnugen"] = fs->make<TH1D>("ptlnugen", ";p_{T}(l#nu);;", 200, 0., 1000.) ;  
  h1_["Mtlnugen"] = fs->make<TH1D>("Mtlnugen", ";M_{T}(l#nu);;", 50, 0., 200.) ;  
  h1_["ptlgen"] = fs->make<TH1D>("ptlgen", ";p_{T}(l) [GeV]", 200, 0., 1000.) ; 
  h1_["ptnugen"] = fs->make<TH1D>("ptnugen", ";p_{T}(#nu) [GeV]", 200, 0., 1000.) ; 
  h1_["etalgen"]  = fs->make<TH1D>("etalgen",  ";#eta(l) [GeV]"  , 40, -5., 5.) ; 
  h1_["etanugen"] = fs->make<TH1D>("etanugen", ";#eta(#nu) [GeV]", 40, -5., 5.) ; 


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
