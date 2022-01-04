// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
//#include "TLorentzVector.h"
//#include "TROOT.h"

using HepMC::GenParticle;
using HepMC::GenVertex;

namespace Rivet {

  class HiggsDecayParTemplate_mgc : public Analysis
  {
  public:
    DEFAULT_RIVET_ANALYSIS_CTOR(HiggsDecayParTemplate_mgc);

    void init() {

      book(hist_pT_Higgs,   "pT_Higgs",20,0,400);
      
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
     
      //find Higgs 4l candidates, if true book pT(H) Histo
        //const double weight = 1.0;

	

    double offshell = 999.0; bool findZ1 = false; bool passZ1 = false; bool makeCuts = true; bool pt_pass=false ;
    vector<double> GENlep_id,GENlep_pt, GENlep_eta, GENlep_phi, GENlep_mass,lep_pt_cp;
    const float Zmass = 91.1876;
    int N, L1, L2, L3, L4;
    L1 = 0; L2 = 0;
    
    for ( const GenParticle *prtcl : Rivet::HepMCUtils::particles(event.genEvent()) ){
            if(prtcl->status()==1 && ( prtcl->pdg_id()==13 || prtcl->pdg_id()== - 13 || prtcl->pdg_id()==11 || prtcl->pdg_id()== - 11) ){//this test added because some events had more than four entries in GENlep_[x]
                GENlep_id.push_back(prtcl->pdg_id());
                FourMomentum lep = prtcl->momentum();
                GENlep_mass.push_back(lep.mass());
                GENlep_pt.push_back(lep.pT());
                GENlep_eta.push_back(lep.eta());
                GENlep_phi.push_back(lep.phi());
                
            }
    }

    cout<<"pT"<<" "<<GENlep_pt<<endl;
    cout<<"--------------------"<<endl;
    
    
    /*const Particles FS = apply<FinalState>(event, "FS").particles();
    
    for (const Particle &p : FS){
        //GENlep_id.push_back(p.pdg_id());
        GENlep_pt.push_back(p.pT());
        GENlep_eta.push_back(p.eta());
        GENlep_phi.push_back(p.phi());
        GENlep_mass.push_back(p.mass());
    }*/
    N = GENlep_id.size();
    lep_pt_cp=GENlep_pt;
    sort(lep_pt_cp.begin(),lep_pt_cp.end());
    for(int i=0;i<4;i++){
        cout<<i<<" "<<"pt "<<lep_pt_cp[i]<<endl;
    }
    if(lep_pt_cp[3] > 20 && lep_pt_cp[2] > 10){pt_pass=true;}
   

    if(pt_pass){
    for(int i=0; i<N; i++){
        for(int j=i+1; j<N; j++){
            
            
            if((GENlep_id[i]+GENlep_id[j])!=0) continue;
            
            FourMomentum li, lj;
            li.setEtaPhiMPt(GENlep_eta[i],GENlep_phi[i],GENlep_mass[i],GENlep_pt[i]);
            lj.setEtaPhiMPt(GENlep_eta[j],GENlep_phi[j],GENlep_mass[j],GENlep_pt[j]);
            
            cout<<"gen lep i id: "<<GENlep_id[i]<<" pt: "<<li.pT()<<" lep j id: "<<GENlep_id[j]<<" pt: "<<lj.pT()<<endl;
            
            if (makeCuts) {
                if ( abs(GENlep_id[i]) == 13 && (li.pT() < 5.0 /*|| abs(li.eta()) > 2.4)*/)) continue;
                if ( abs(GENlep_id[i]) == 11 && (li.pT() < 5.0 /*|| abs(li.eta()) > 2.5)*/)) continue;//pt<7 
                
                
                if ( abs(GENlep_id[j]) == 13 && (lj.pT() < 5.0 /*|| abs(lj.eta()) > 2.4)*/)) continue;
                if ( abs(GENlep_id[j]) == 11 && (lj.pT() < 5.0 /*|| abs(lj.eta()) > 2.5)*/)) continue;//pt<7
               
            }
            
            FourMomentum mll = li+lj;
            cout<<"gen mass ij: "<<mll.mass()<<endl;
            
            if(abs(mll.mass()-Zmass)<offshell){
                double mZ1 = mll.mass();
                cout<<"foundZ1"<<endl;
                L1 = i; L2 = j; findZ1 = true; offshell = abs(mZ1-Zmass);
            }
        }
    }
    
    FourMomentum l1, l2;
    l1.setEtaPhiMPt(GENlep_eta[L1],GENlep_phi[L1],GENlep_mass[L1],GENlep_pt[L1]);
    l2.setEtaPhiMPt(GENlep_eta[L2],GENlep_phi[L2],GENlep_mass[L2],GENlep_pt[L2]);
    FourMomentum ml1l2 = l1+l2;
    
    if(ml1l2.mass()>12 && ml1l2.mass()<120 && findZ1) passZ1 = true;//40<m<120
    
    double pTL34 = 0.0; bool findZ2 = false;
    //bool m4lwindow=false; double window_lo=70.0; double window_hi=140.0;
    
    //cout<<"findZ2"<<endl;
    for(int i=0; i<N; i++){
        if(i==L1 || i==L2) continue; // can not be the lep from Z1
        for(int j=i+1; j<N; j++){
            if(j==L1 || j==L2) continue; // can not be the lep from Z1
            if((GENlep_id[i]+GENlep_id[j])!=0) continue;
            
            FourMomentum li, lj;
            li.setEtaPhiMPt(GENlep_eta[i],GENlep_phi[i],GENlep_mass[i],GENlep_pt[i]);
            lj.setEtaPhiMPt(GENlep_eta[j],GENlep_phi[j],GENlep_mass[j],GENlep_pt[j]);
            FourMomentum Z2 = li+lj;

            cout<<"gen lep i id: "<<GENlep_id[i]<<" pt: "<<li.pT()<<" lep j id: "<<GENlep_id[j]<<" pt: "<<lj.pT()<<endl;
            
            if (makeCuts) {
                if ( abs(GENlep_id[i]) == 13 && (li.pT() < 5.0 /*|| abs(li.eta()) > 2.4)*/)) continue;
                if ( abs(GENlep_id[i]) == 11 && (li.pT() < 5.0 /*|| abs(li.eta()) > 2.5)*/)) continue;//pt<7
               
                
                if ( abs(GENlep_id[j]) == 13 && (lj.pT() < 5.0 /*|| abs(lj.eta()) > 2.4)*/)) continue;
                if ( abs(GENlep_id[j]) == 11 && (lj.pT() < 5.0 /*|| abs(lj.eta()) > 2.5)*/)) continue;//pt<7
              
            }
            
            if ( /*(li.pT()+lj.pT())>=pTL34*/ true ) {
                double mZ2 = Z2.mass();
                cout<<"GEN mZ2: "<<mZ2<<endl;
                if( (mZ2>12 && mZ2<120)) {
                    L3 = i; L4 = j; findZ2 = true;
                    pTL34 = li.pT()+lj.pT();
                    cout<<"is the new GEN cand"<<endl;
                    //if (m4l>window_lo && m4l<window_hi) m4lwindow=true;
                } 
            }
            
        } // lj
    } // li

    unsigned int tmp_;
    if(GENlep_pt[L1]<GENlep_pt[L2])    {tmp_=L1;    L1=L2;    L2=tmp_;}
    if(GENlep_pt[L3]<GENlep_pt[L4])    {tmp_=L3;    L3=L4;    L4=tmp_;}
    
    
    if(passZ1 && findZ2){
        for ( const GenParticle *ptcl : Rivet::HepMCUtils::particles(event.genEvent()) ){
            if( PID::isHiggs(ptcl->pdg_id()) ){
                hist_pT_Higgs->fill(ptcl->momentum().perp());
	    }
        }
    }


    
    }//end of pt_pass condition ( 1st > 20 && 2nd > 10 )

    GENlep_id.clear();
    GENlep_eta.clear();
    GENlep_phi.clear();
    GENlep_mass.clear();
    GENlep_pt.clear();

    }//end of analyze event

    /// Normalise histograms etc., after the run
    void finalize() {
      // MSG_DEBUG("crossSection sumOfWeights ");
      MSG_DEBUG("crossSection " << crossSection() << " sumOfWeights " << sumOfWeights());

    }

  private:
      Histo1DPtr hist_pT_Higgs;

  };

  DECLARE_RIVET_PLUGIN(HiggsDecayParTemplate_mgc);
}
