// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"

#include "DataFormats/L1TrackTrigger/interface/L1TkEGTauParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "L1Trigger/Phase2L1ParticleFlow/interface/L1TPFUtils.h"
//#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "CommonTools/Utils/interface/StringObjectFunction.h"

#include <cstdint>
#include <TTree.h>
#include <TLorentzVector.h>

class TauNTuplizer3 : public edm::one::EDAnalyzer<edm::one::SharedResources,edm::one::WatchRuns>  {
    public:
        explicit TauNTuplizer3(const edm::ParameterSet&);
        ~TauNTuplizer3();

    private:
        virtual void beginJob() override;
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
        virtual void beginRun(edm::Run const&, edm::EventSetup const& iSetup) override {}
        virtual void endRun(edm::Run const&, edm::EventSetup const& iSetup) override { } // framework wants this to be implemented

        edm::EDGetTokenT<std::vector<l1t::L1TkEGTauParticle>> L1PFTaus_;
        edm::EDGetTokenT<std::vector<reco::GenParticle>> genparticles_;
        float dr2Max_, minPtRatio_;
        TTree *tree_;
        uint32_t run_, lumi_; uint64_t event_; uint64_t eventcount_;

        struct McVars {
            float pt1, eta1, phi1, dr1;
  	    float pt2, eta2, phi2, dr2;
            int   id;
            void makeBranches(TTree *tree) {
                tree->Branch("gendr1", &dr1, "gendr1/F");
                tree->Branch("genpt1", &pt1, "genpt1/F");
                tree->Branch("geneta1", &eta1, "geneta1/F");
                tree->Branch("genphi1", &phi1, "genphi1/F");

                tree->Branch("gendr2", &dr2, "gendr2/F");
                tree->Branch("genpt2", &pt2, "genpt2/F");
                tree->Branch("geneta2", &eta2, "geneta2/F");
                tree->Branch("genphi2", &phi2, "genphi2/F");
            }
            void clear() {
                pt1 = 0; eta1 = 0; phi1 = 0; dr1 = -999; id = 0;
                pt2 = 0; eta2 = 0; phi2 = 0; dr2 = -999;
            }
	  void fill(const reco::GenParticle &c) { 
            TLorentzVector lVec; 
	    lVec.SetPtEtaPhiM(pt1,eta1,phi1,0);
	    TLorentzVector lNewVec;
	    lNewVec.SetPtEtaPhiM(c.pt(),c.eta(),c.phi(),0);
	    lVec += lNewVec;
	    id  = c.pdgId();
	    pt1  = lVec.Pt(); 
	    eta1 = lVec.Eta(); 
	    phi1 = lVec.Phi();
	  }
	  void fill(TLorentzVector iVec,float iDR) { 
	    if(iVec.Pt() > pt1) {
	      //if(fabs(iDR) < fabs(dr1)) {
	      dr2  = dr1;
	      pt2  = pt1;
	      eta2 = eta1;
	      phi2 = phi1;

	      dr1  = iDR;
	      pt1  = iVec.Pt(); 
	      eta1 = iVec.Eta(); 
	      phi1 = iVec.Phi();
	      //if(iDR < 0) std::cout << "--- " << iDR << " -- " << dr1 << " -- " << fabs(iDR) << " -- " << fabs(dr1) << std::endl;
	    } else if(iVec.Pt() > pt2) { 
	      dr2  = iDR;
	      pt2  = iVec.Pt(); 
	      eta2 = iVec.Eta(); 
	      phi2 = iVec.Phi();
	    }
	  }
        } mc_;

        class RecoVar {
            public:
                RecoVar(const std::string & name, const std::string & expr) : name_(name), expr_(expr,true) {}
                void makeBranch(TTree *tree) {
                    tree->Branch(name_.c_str(), &val_, (name_+"/F").c_str());
                }
                void fill(const reco::Candidate & c) {
                    val_ = expr_(c);
                }
            private:
                std::string name_;
                StringObjectFunction<reco::Candidate> expr_;
                float val_;
        };
        std::vector<RecoVar> reco_;

 
};

TauNTuplizer3::TauNTuplizer3(const edm::ParameterSet& iConfig) :
  L1PFTaus_    (consumes<std::vector<l1t::L1TkEGTauParticle>>(iConfig.getParameter<edm::InputTag>("src"))),
  genparticles_(consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genParticles"))),
  dr2Max_(std::pow(iConfig.getParameter<double>("drMax"), 2)),
  minPtRatio_(float(iConfig.getParameter<double>("minRecoPtOverGenPt"))) {
 
    usesResource("TFileService");
    edm::Service<TFileService> fs;
    tree_ = fs->make<TTree>("tree","tree");
    tree_->Branch("run",  &run_,  "run/i");
    tree_->Branch("lumi", &lumi_, "lumi/i");
    tree_->Branch("event", &event_, "event/l");
    tree_->Branch("event2", &eventcount_, "eventcount/l");

    edm::ParameterSet vars = iConfig.getParameter<edm::ParameterSet>("variables");
    auto reconames = vars.getParameterNamesForType<std::string>();
    for (const std::string & name : reconames) {
        reco_.emplace_back(name, vars.getParameter<std::string>(name));
    }
}

TauNTuplizer3::~TauNTuplizer3() { }
void  TauNTuplizer3::beginJob() {
    mc_.makeBranches(tree_);
    for (auto & v : reco_) v.makeBranch(tree_);
    eventcount_=0;
}
void TauNTuplizer3::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    run_  = iEvent.id().run();
    lumi_ = iEvent.id().luminosityBlock();
    event_ = iEvent.id().event();
    eventcount_++;
    
    edm::Handle<std::vector<reco::GenParticle>> genparticles;
    iEvent.getByToken(genparticles_, genparticles);

    edm::Handle<  l1t::L1TkEGTauParticleCollection > l1PFTaus;
    try { 
      iEvent.getByToken( L1PFTaus_, l1PFTaus);
    } catch(...) { 
      return;
    }
    l1t::L1TkEGTauParticle dummy;
    std::vector<int> matchedGen;
    std::vector<int> matchedTau;
    std::vector<TLorentzVector> taus;
    int pPass  = 1; 
    int lId    = -1;
    while(pPass != 0) { 
      pPass = 0;
      int   pCount = -1;
      float pDR = -1;
      for (const reco::GenParticle &gen : *genparticles) {
	pCount++;
	if(!gen.isDirectPromptTauDecayProductFinalState()) continue;
	if(fabs(gen.pdgId()) > 10 && fabs(gen.pdgId()) < 17) continue;
	bool pMatch = false;
	for(unsigned int i0 = 0; i0 < matchedGen.size(); i0++) if(matchedGen[i0] == pCount) pMatch = true;
	if(pMatch) continue;
	//std::cout << "===> " << gen.pdgId() << " -- " << gen.pt() << " -- " << gen.eta() << " -- " << gen.phi() << " -- " << pDR << "-- " << pPass << " -- " << matchedGen.size() << " -- " << taus.size() << std::endl;
	if(pDR > 0){ 
	  float dr2 = deltaR2(taus[lId].Eta(), taus[lId].Phi(), gen.eta(), gen.phi());
	  if(dr2 > pDR) continue;
	  pPass++;
	  TLorentzVector lVec;
	  lVec.SetPtEtaPhiM(gen.pt(),gen.eta(),gen.phi(),0);
	  taus[lId]+=lVec;
	  matchedGen.push_back(pCount);
	} else { 
	  pDR = 0.25*0.25;
	  TLorentzVector lVec;
	  lVec.SetPtEtaPhiM(gen.pt(),gen.eta(),gen.phi(),0);
	  taus.push_back(lVec);
	  matchedGen.push_back(pCount);
	  lId++;
	}
      }
    }
    if(taus.size() > 2) std::cout << "!!!!===> check " << taus.size() << std::endl;
    for (unsigned int i = 0, n = l1PFTaus->size(); i < n; ++i) {
        const auto & c = (*l1PFTaus)[i];
	//float dr2best = dr2Max_; 
	mc_.clear();
	//int pIndex=-1;
	for(unsigned i0 = 0; i0 < taus.size(); i0++) { 
	  //bool pXMatch = false;
	  //for(unsigned i1 = 0; i1 < matchedTau.size(); i1++) if(int(i0) == matchedTau[i1]) pXMatch = true;
	  //if(pXMatch) continue;
	  float dr2 = deltaR2(taus[i0].Eta(), taus[i0].Phi(), c.eta(), c.phi());
	  //if(dr2 < dr2best) {
	  //dr2best = dr2;
	  mc_.fill(taus[i0],std::sqrt(dr2));
	  //pIndex = i0;
	  //}
	}
	for (auto & v : reco_) v.fill(c);
	//matchedTau.push_back(pIndex);
	tree_->Fill();
	mc_.clear();
    }
    if(l1PFTaus->size() == 0) {
      for(int i0 = 0; i0 < int(taus.size()); i0++) { 
	//bool pMatch = false;
	//for(unsigned i1 = 0; i1 < matchedTau.size(); i1++) if(int(i0) == matchedTau[i1]) pMatch = true;
	//if(pMatch) continue;
	float dr2 = 99;
	mc_.fill(taus[i0],-1.*std::sqrt(dr2));
	//pIndex = i0;
	//pMatch++;
	//}
      }
      for (auto & v : reco_) v.fill(dummy);
      tree_->Fill();
      mc_.clear();
    }
}

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TauNTuplizer3);
