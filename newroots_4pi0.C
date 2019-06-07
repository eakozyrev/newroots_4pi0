
#define newroots_4pi0_cxx
#include "newroots_4pi0.h"

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TProfile.h"
#include "TFile.h"
#include "TMatrixTSym.h"
#include "TMath.h"
#include "TVector.h"
#include "TLorentzVector.h"
#include "TStyle.h"
#include "TBenchmark.h"
#include <iostream>
#include <string>
#define Mpi 139.57
#define mpiz 134.9766

#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TVector.h"
#include "TH1D.h"
#include <fstream>
#include <string>


#include "NewCmdKinFitter/install/inc/CmdKinFitter.h"
#include "NewCmdKinFitter/install/inc/CmdKinFit2Pi4Gamma.h"
#include "NewCmdKinFitter/install/inc/CmdKinFit2Pi2Pi0_6C.h"
#include "NewCmdKinFitter/install/inc/CmdKinFit2Pi2Gamma1Pi0.h"
#include "NewCmdKinFitter/install/inc/CmdKinFit2Pi2GammaPiLoosen.h"
#include "NewCmdKinFitter/install/inc/CmdChargedParticle.h" 
#include "NewCmdKinFitter/install/inc/CmdPhotonParticle.h"

//using namespace cmd3;

void perform(std::string name1){
	ifstream streamm(name1.c_str()); string ss,ss1,ss2;
	TFile *oldfile;
        TTree *oldtree;
	while(streamm.eof()==0){
		streamm>>ss >> ss1 >> ss2;
		if(ss == "stop")break;
		cout << ss +ss1<< endl;
		oldfile =TFile::Open((ss+ss1).c_str());
		oldtree = (TTree*)oldfile->Get("tr_ph");
                newroots_4pi0 a(oldtree);
		bool boolmc = true;
		if(ss2 == "exp")boolmc = false;
		a.Loop(ss1, boolmc); 
	}
}

double tdedx_high(double x){

    return 0.3*(0.5*(4200.*(1+TMath::Exp(-0.0095*(x-190.))+TMath::Exp(-0.024*(x-190.)))+4500)+8000.);

}










void newroots_4pi0::Loop(string namefile, bool boolmc)
{


  
   if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntriesFast();

    TFile newfile(("/store11/eakozyrev/4pi/4pisel/"+namefile).c_str(),"recreate");

    TTree *newtree = fChain->GetTree()->CloneTree(0);
    TLorentzVector Lpip,Lpim,Lph1,Lph2,Lph3,Lph4, Ltr1, Ltr2, Lgamma1, Lgamma2;
    double Chi25C,m2g0,m2g, P_poln, E_poln,mass34ph1tr,theta_miss_pi,mom_miss_pi,m2pirest,Xi2misspi0; 
    double Chi_merge = 300000;
    double Xi2_Ph3=0;
    double m2g_merge;
    double masspi0_with_omega,deltaE1tr;
    double mass34ph1pi0_kin,mass34ph1pi0,misspi0mass;
    int mode1,mode2,mode3,mode4,Err_merge,Err_Chi25C;
    int Ngen = nentries;
    int correspond[20];
    newtree->Branch("Ngen",&Ngen);
    newtree->Branch("Chi25C",&Chi25C);
    newtree->Branch("Err_Chi25C",&Err_Chi25C);
    newtree->Branch("Lpim",&Lpim);
    newtree->Branch("Lpip",&Lpip);
    newtree->Branch("Lph1",&Lph1);
    newtree->Branch("Lph2",&Lph2);
    newtree->Branch("Lph3",&Lph3);
    newtree->Branch("Lph4",&Lph4);
    newtree->Branch("m2g0",&m2g0);
    newtree->Branch("m2g",&m2g);
    newtree->Branch("m2pirest",&m2pirest);
    newtree->Branch("P_poln",&P_poln);
    newtree->Branch("E_poln",&E_poln);
    newtree->Branch("mode1",&mode1);
    newtree->Branch("mode2",&mode2);
    newtree->Branch("mode3",&mode3);
    newtree->Branch("correspond",correspond,"correspond[20]/I");
    newtree->Branch("deltaE1tr",&deltaE1tr);
    newtree->Branch("mass34ph1tr",&mass34ph1tr);
    newtree->Branch("theta_miss_pi",&theta_miss_pi);
    newtree->Branch("mom_miss_pi",&mom_miss_pi);
    newtree->Branch("mass34ph1pi0_kin",&mass34ph1pi0_kin);
    newtree->Branch("mass34ph1pi0",&mass34ph1pi0);
    newtree->Branch("misspi0mass",&misspi0mass);
    newtree->Branch("Xi2misspi0",&Xi2misspi0);
    newtree->Branch("Ltr1",&Ltr1);
    newtree->Branch("Ltr2",&Ltr2);
    newtree->Branch("Lgamma1",&Lgamma1);
    newtree->Branch("Lgamma2",&Lgamma2);
    newtree->Branch("Chi_merge",&Chi_merge);
    newtree->Branch("Xi2_Ph3",&Xi2_Ph3);
    newtree->Branch("m2g_merge",&m2g_merge);
    newtree->Branch("Err_merge",&Err_merge);

    cout << nentries << endl;
    // window for mpi0
    double Cut_mpi0_before = 45.; double Cut_mpi0_above = 255.;
    double th_ph_min = 1.;double th_ph_max = 3.1415 - th_ph_min;
    double th_tr_min = 0.6;double th_tr_max = 3.1415 - th_tr_min;

    TVectorT<double> Pip(5),Pim(5),Ph1(5),Ph2(5),Ph3(5),Ph4(5);
    TMatrixDSym cov_Pip(5),cov_Pim(5),cov_Ph1(5),cov_Ph2(5),cov_Ph3(5),cov_Ph4(5);
   
    int numbers[100000][4],numbers_alone[100000][2];
    int npi0 = 0;int npi0_alone = 0;
    TLorentzVector PPh[1000];TLorentzVector PTr[1000];

   Long64_t nbytes = 0, nb = 0;
   double N_events = 0;
   //bool boolmc = true;
   CmdKinFit2Pi4Gamma event4piw(TLorentzVector(0.,0.,0.,1000));
   CmdKinFit2Pi2Pi0_6C event2Pi2Pi0_6C(TLorentzVector(0.,0.,0.,1000));
   // getchar();
   CmdKinFit2Pi2Gamma1Pi0 event2pi2gamma1Pi0(TLorentzVector(0.,0.,0.,1000));
   // getchar();
   CmdKinFit2Pi2GammaPiLoosen event4pilossenpi0;
   CmdKinFitter Fitter;
   Fitter.SetDebug(0); 
   cout<<"after Fiitter"<<endl;
   // getchar();
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if(jentry > 10000)break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if(finalstate_id == 3)continue;
      if(jentry%500 == 0)  cout << (double)jentry/nentries << "  " << ebeam << "   boolmc = " << boolmc << endl;
      Chi25C = 300.;
      if(fabs(ebeam - emeas)>50)emeas = ebeam;
      TLorentzVector Pebeam(0.,0.,0.,2.*emeas);
      event4piw.SetTotalMoment(Pebeam);
      event2Pi2Pi0_6C.SetTotalMoment(Pebeam);
      event2pi2gamma1Pi0.SetTotalMoment(Pebeam);
      event4pilossenpi0.setMoment(Pebeam);
      m2g = 0;
      m2g0 = 0;
      P_poln = 0;
      E_poln = 0;
      npi0 = 0;
      npi0_alone = 0;
      bool firsttry = true;
      // photons 4-momenta
      for(int k = 0; k < nph; k++){
         PPh[k] = TLorentzVector(phen[k]*TMath::Cos(phphi[k])*TMath::Sin(phth[k]),
				phen[k]*TMath::Sin(phphi[k])*TMath::Sin(phth[k]),
				phen[k]*TMath::Cos(phth[k]),
				phen[k]);
      } 
      // tracks 4-momenta
      for(int k = 0; k < nt; k++){
	PTr[k] = TLorentzVector(tptot[k]*TMath::Cos(tphi[k])*TMath::Sin(tth[k]),
				tptot[k]*TMath::Sin(tphi[k])*TMath::Sin(tth[k]),
				tptot[k]*TMath::Cos(tth[k]),
				sqrt(tptot[k]*tptot[k]+Mpi*Mpi));
      }
    // the list of pi0
    // npi0_alone - the number of pi0 in an event
    for(int i = 0; i < nph-1; i++){
	  for(int j = i + 1; j < nph; j++){
	    if((PPh[i] + PPh[j]). M() < Cut_mpi0_above && (PPh[i] + PPh[j]). M() > Cut_mpi0_before){
	      if(phth[i] < th_ph_min || phth[i] > th_ph_max)continue;
	      if(phth[j] < th_ph_min || phth[j] > th_ph_max)continue;	      
	    numbers_alone[npi0_alone][0] = i; numbers_alone[npi0_alone][1] = j; npi0_alone++;
	  }		
	 }
    }
    correspond[0] = -1;
    correspond[1] = -1;  
    correspond[2] = -1;
    correspond[3] = -1;
    correspond[4] = -1;
    correspond[5] = -1;
    correspond[6] = -1;
    correspond[7] = -1;
    correspond[8] = -1;
    correspond[9] = -1;
    correspond[10] = -1;
    correspond[11] = -1;
    correspond[12] = -1;
    correspond[13] = -1;
    correspond[14] = -1;
    correspond[15] = -1;
    mode1 = 0;    mode2 = 0;    mode3 = 0; mode4 = 0; deltaE1tr = -1;
    double mass = 30;

    if(nt < 1 || npi0_alone < 1){            //-----------------------------REQUIRE-------------------------------------------
     	if(boolmc == true){newtree->Fill();}
	continue;
    }
    int indexes_best_pi01 = 1; 
    for(int i = 0; i < npi0_alone; i++){
        int k = numbers_alone[indexes_best_pi01][0];int l =  numbers_alone[indexes_best_pi01][1];
        int k1 = numbers_alone[i][0];int l1 =  numbers_alone[i][1];
	if(fabs((PPh[k1] + PPh[l1]). M() - mpiz) < fabs((PPh[k] + PPh[l]). M() - mpiz)){indexes_best_pi01 = i;}
    }

    int indexes_best_pi02 = 1; 
    if(indexes_best_pi01 == 1)indexes_best_pi02=0;
    for(int i = 0; i < npi0_alone; i++){
        if(i == indexes_best_pi01)continue;
        int k = numbers_alone[indexes_best_pi02][0];int l =  numbers_alone[indexes_best_pi02][1];
        int k1 = numbers_alone[i][0];int l1 =  numbers_alone[i][1];
	if(fabs((PPh[k1] + PPh[l1]). M() - mpiz) < fabs((PPh[k] + PPh[l]). M() - mpiz)){indexes_best_pi02 = i;}
    }

    int index_best_tr1 = 0;
    for(int one = 0; one < nt; one++){
  	  if(fabs(tz[one]) > 10. || fabs(trho[one]) > 0.2 || tdedx[one] > tdedx_high(tptot[one])){continue;}
	  if(tnhit[one] > tnhit[index_best_tr1])index_best_tr1 = one;
   }    
   if(fabs(tz[index_best_tr1]) > 10. || fabs(trho[index_best_tr1]) > 0.2 || tdedx[index_best_tr1] > tdedx_high(tptot[index_best_tr1])){index_best_tr1=-1;}
   int index_best_tr2 = index_best_tr1;
   for(int one = 0; one < nt; one++){
  	  if(tcharge[one]*tcharge[index_best_tr2] > 0 || fabs(tz[one]) > 10. || fabs(trho[one]) > 0.2 || tdedx[one] > tdedx_high(tptot[one])){continue;}
	  if(tnhit[one] > tnhit[index_best_tr2])index_best_tr2 = one;
   }
   if(index_best_tr2 == index_best_tr1){index_best_tr2=-1;}

    
    for(int i = 0; i < nph-3; i++){
	for(int j = i + 1; j < nph-2; j++){
	  for(int k = j +1; k < nph-1; k++){
	    for(int l = k + 1; l < nph; l++){
	      if(phth[i] < (th_ph_min - 0.2) || phth[i] > (th_ph_max + 0.2))continue;
	      if(phth[j] < (th_ph_min - 0.2) || phth[j] > (th_ph_max + 0.2))continue;	      
	      if(phth[k] < (th_ph_min - 0.2) || phth[k] > (th_ph_max + 0.2))continue;
	      if(phth[l] < (th_ph_min - 0.2) || phth[l] > (th_ph_max + 0.2))continue;	      
	      if((PPh[i] + PPh[j]). M() < Cut_mpi0_above && (PPh[i] + PPh[j]). M() > Cut_mpi0_before
		 && (PPh[k] + PPh[l]). M() < Cut_mpi0_above && (PPh[k] + PPh[l]). M() > Cut_mpi0_before){
		numbers[npi0][0] = i;numbers[npi0][1] = j;numbers[npi0][2] = k;numbers[npi0][3] = l; npi0++;
	      }
	      if((PPh[i] + PPh[k]). M() < Cut_mpi0_above && (PPh[i] + PPh[k]). M() > Cut_mpi0_before
		 && (PPh[j] + PPh[l]). M() < Cut_mpi0_above && (PPh[j] + PPh[l]). M() > Cut_mpi0_before){
		numbers[npi0][0] = i;numbers[npi0][1] = k;numbers[npi0][2] = j;numbers[npi0][3] = l; npi0++;
	      }
	      if((PPh[i] + PPh[l]). M() < Cut_mpi0_above && (PPh[i] + PPh[l]). M() > Cut_mpi0_before
		 && (PPh[k] + PPh[j]). M() < Cut_mpi0_above && (PPh[k] + PPh[j]). M() > Cut_mpi0_before){
		numbers[npi0][0] = i;numbers[npi0][1] = l;numbers[npi0][2] = k;numbers[npi0][3] = j; npi0++;
	      }
	    }
	  }
	}
      }
       
      double sigma2_rho = 0.3;
      double sigma2_z = 0.3;


      //***************************************** 
      //   2pi0 pi+ pi- hypothesis 
      //*****************************************
      Chi25C = 300.;
      Err_Chi25C=10;
      for(int one = 0; one < nt-1; one++){
	  	  for(int two = one + 1; two < nt; two++){
		          if(tcharge[one]*tcharge[two] > 0)continue;
	  	  	  if(fabs(tz[one]) > 13. || fabs(trho[one]) > 0.2 || tdedx[one] > tdedx_high(tptot[one])){continue;}
	  	  	  if(fabs(tz[two]) > 13. || fabs(trho[two]) > 0.2 || tdedx[two] > tdedx_high(tptot[two])){continue;}
	  	  	  for(int j = 0; j < npi0; j++){
	    
			    int k = numbers[j][0];
			    int l = numbers[j][1];
			    int m = numbers[j][2];
			    int n = numbers[j][3];
			    
			    if(fabs((PTr[one] + PTr[two] + PPh[k] + PPh[l]+ PPh[m]+ PPh[n]).E() - Pebeam.E()) > 500){continue;}
			    if(fabs((PTr[one] + PTr[two] + PPh[k] + PPh[l]+ PPh[m]+ PPh[n]).P()) > 500){continue;}
			    
			    Pip(0)= tptot[one]*sin(tth[one]);
			    Pip(1)= cos(tth[one])/sin(tth[one]);
			    Pip(2)= tphi[one];
			    Pip(3) = trho[one];
			    Pip(4) = tz[one];
			    Pim(0)= tptot[two]*sin(tth[two]);
			    Pim(1)= cos(tth[two])/sin(tth[two]);
			    Pim(2)= tphi[two];
			    Pim(3) = trho[two];
			    Pim(4) = tz[two];
			    double Z_photons=0.5*(Pim(4)+Pip(4));
			    //cout<<"Z_Photons = "<<Z_photons<<endl;
			    //getchar();
			    cov_Pip(0,0)= terr0[one][0][0];cov_Pip(1,0)= terr0[one][3][0];cov_Pip(2,0)= terr0[one][1][0];cov_Pip(3,0)=terr0[one][2][0];cov_Pip(4,0)= terr0[one][4][0];
			    cov_Pip(0,1)= terr0[one][0][3];cov_Pip(1,1)= terr0[one][3][3];cov_Pip(2,1)=terr0[one][1][3];cov_Pip(3,1)=terr0[one][2][3];cov_Pip(4,1)= terr0[one][4][3];
			    cov_Pip(0,2)= terr0[one][0][1];cov_Pip(1,2)= terr0[one][3][1];cov_Pip(2,2)=terr0[one][1][1];cov_Pip(3,2)=terr0[one][2][1];cov_Pip(4,2)= terr0[one][4][1];
			    cov_Pip(0,3)= terr0[one][0][2];cov_Pip(1,3)= terr0[one][3][2];cov_Pip(2,3)=terr0[one][1][2];cov_Pip(3,3)=terr0[one][2][2];cov_Pip(4,3)= terr0[one][4][2];
			    cov_Pip(0,4)= terr0[one][0][4];cov_Pip(1,4)= terr0[one][3][4];cov_Pip(2,4)=terr0[one][1][4];cov_Pip(3,4)=terr0[one][2][4];cov_Pip(4,4)= terr0[one][4][4];
			    
			    cov_Pim(0,0)= terr0[two][0][0];cov_Pim(1,0)= terr0[two][3][0];cov_Pim(2,0)= terr0[two][1][0];cov_Pim(3,0)=terr0[two][2][0];cov_Pim(4,0)= terr0[two][4][0];
			    cov_Pim(0,1)= terr0[two][0][3];cov_Pim(1,1)= terr0[two][3][3];cov_Pim(2,1)=terr0[two][1][3];cov_Pim(3,1)=terr0[two][2][3];cov_Pim(4,1)= terr0[two][4][3];
			    cov_Pim(0,2)= terr0[two][0][1];cov_Pim(1,2)= terr0[two][3][1];cov_Pim(2,2)=terr0[two][1][1];cov_Pim(3,2)=terr0[two][2][1];cov_Pim(4,2)= terr0[two][4][1];
			    cov_Pim(0,3)= terr0[two][0][2];cov_Pim(1,3)= terr0[two][3][2];cov_Pim(2,3)=terr0[two][1][2];cov_Pim(3,3)=terr0[two][2][2];cov_Pim(4,3)= terr0[two][4][2];
			    cov_Pim(0,4)= terr0[two][0][4];cov_Pim(1,4)= terr0[two][3][4];cov_Pim(2,4)=terr0[two][1][4];cov_Pim(3,4)=terr0[two][2][4];cov_Pim(4,4)= terr0[two][4][4];
			    
			    double R = phrho[k];
			    Ph1(0)=phen[k];
			    Ph1(1)=R;
			    Ph1(2)=phphi[k];
			    Ph1(3)= z0 + R*cos(phth[k])/sin(phth[k]);
			    Ph1(4)=phth[k];Ph1(4)=atan2(Ph1(1),Ph1(3)-Z_photons);
			    R = phrho[l];
			    Ph2(0)=phen[l];
			    Ph2(1)=R;
			    Ph2(2)=phphi[l];
			    Ph2(3)=z0 + R*cos(phth[l])/sin(phth[l]);
			    Ph2(4)=phth[l];Ph2(4)=atan2(Ph2(1),Ph2(3)-Z_photons);
			    sigma2_rho = 0.3;
			    sigma2_z = 0.3;
			    if(phflag[k] == 3){sigma2_rho = sigma2_z*pow(tan(phth[k]),2) + pow(phrho[k]*pherr[k][1]/sin(phth[k])*cos(phth[k]),2);}
			    else {sigma2_z = sigma2_rho/pow(tan(phth[k]),2) + pow(phrho[k]*pherr[k][1],2)/pow(sin(phth[k]),4);}
			    
			    cov_Ph1(0,0) = 1./(pherr[k][0]*pherr[k][0]);
			    cov_Ph1(1,1) = 1./sigma2_rho;
			    cov_Ph1(2,2) = 1./(pherr[k][2]*pherr[k][2]);
			    cov_Ph1(3,3) = 1./sigma2_z;
			    cov_Ph1(4,4) = 1./(pherr[k][1]*pherr[k][1]);
			    
			    sigma2_rho = 0.3;
			    sigma2_z = 0.3;
			    if(phflag[l] == 3){sigma2_rho = sigma2_z*pow(tan(phth[l]),2) + pow(phrho[l]*pherr[l][1]/sin(phth[l])*cos(phth[l]),2);}
			    else {sigma2_z = sigma2_rho/pow(tan(phth[l]),2) + pow(phrho[l]*pherr[l][1],2)/pow(sin(phth[l]),4);}
			    
			    cov_Ph2(0,0) = 1./(pherr[l][0]*pherr[l][0]);
			    cov_Ph2(1,1) = 1./sigma2_rho;
			    cov_Ph2(2,2) = 1./(pherr[l][2]*pherr[l][2]);
			    cov_Ph2(3,3) = 1./sigma2_z;
			    cov_Ph2(4,4) = 1./(pherr[l][1]*pherr[l][1]);
			    
			    R = phrho[m];
			    Ph3(0)=phen[m];
			    Ph3(1)=R;
			    Ph3(2)=phphi[m];
			    Ph3(3)= z0 + R*cos(phth[m])/sin(phth[m]);//tz[one]+
			    Ph3(4)=phth[m];Ph3(4)=atan2(Ph3(1),Ph3(3)-Z_photons);
			    R = phrho[n];
			    Ph4(0)=phen[n];
			    Ph4(1)=R;
			    Ph4(2)=phphi[n];
			    Ph4(3)=z0 + R*cos(phth[n])/sin(phth[n]);//tz[one]+
			    Ph4(4)=phth[n];Ph4(4)=atan2(Ph4(1),Ph4(3)-Z_photons);
			    
			    sigma2_rho = 0.3;
			    sigma2_z = 0.3;
			    if(phflag[m] == 3){sigma2_rho = sigma2_z*pow(tan(phth[m]),2) + pow(phrho[m]*pherr[m][1]/sin(phth[m])*cos(phth[m]),2);}
			    else {sigma2_z = sigma2_rho/pow(tan(phth[m]),2) + pow(phrho[m]*pherr[m][1],2)/pow(sin(phth[m]),4);}
			    
			    cov_Ph3(0,0) = 1./(pherr[m][0]*pherr[m][0]);
			    cov_Ph3(1,1) = 1./sigma2_rho;
			    cov_Ph3(2,2) = 1./(pherr[m][2]*pherr[m][2]);
			    cov_Ph3(3,3) = 1./sigma2_z;
			    cov_Ph3(4,4) = 1./(pherr[m][1]*pherr[m][1]);
			    
			    sigma2_rho = 0.3;
			    sigma2_z = 0.3;
			    if(phflag[n] == 3){sigma2_rho = sigma2_z*pow(tan(phth[n]),2) + pow(phrho[n]*pherr[n][1]/sin(phth[n])*cos(phth[n]),2);}
			    else {sigma2_z = sigma2_rho/pow(tan(phth[n]),2) + pow(phrho[n]*pherr[n][1],2)/pow(sin(phth[n]),4);}
			    
			    cov_Ph4(0,0) = 1./(pherr[n][0]*pherr[n][0]);
			    cov_Ph4(1,1) = 1./sigma2_rho;
			    cov_Ph4(2,2) = 1./(pherr[n][2]*pherr[n][2]);
			    cov_Ph4(3,3) = 1./sigma2_z;
			    cov_Ph4(4,4) = 1./(pherr[n][1]*pherr[n][1]);
			    
			    // cmd3::CmdChargeParticle ppi("Pi Plus",Mpi,Pip,cov_Pip.Invert());
			    // cmd3::CmdChargeParticle mpi("Pi Minus",Mpi,Pim,cov_Pim.Invert());
			    // cmd3::CmdPhotonParticle ph1("Photon 1",Ph1,cov_Ph1);
			    // cmd3::CmdPhotonParticle ph2("Photon 2",Ph2,cov_Ph2);
			    // cmd3::CmdPhotonParticle ph3("Photon 3",Ph3,cov_Ph3);
			    // cmd3::CmdPhotonParticle ph4("Photon 4",Ph4,cov_Ph4);
			    //double mpi0000000 = (ph3.GetMoment()+ph4.GetMoment()).M();
			    
			    //event4piw.PutParticle(ppi,mpi,ph1,ph2,ph3,ph4);
			    event2Pi2Pi0_6C.GetParticle("PiPl")->setParticleContent(Pip, cov_Pip.Invert());  CmdBaseParticle& ppi= *event2Pi2Pi0_6C.GetParticle("PiPl");
			    event2Pi2Pi0_6C.GetParticle("PiMi")->setParticleContent(Pim, cov_Pim.Invert());  CmdBaseParticle& mpi= *event2Pi2Pi0_6C.GetParticle("PiMi");
			    event2Pi2Pi0_6C.GetParticle("Ph0")->setParticleContent(Ph1, cov_Ph1);            CmdBaseParticle& ph1 = *event2Pi2Pi0_6C.GetParticle("Ph0");
			    event2Pi2Pi0_6C.GetParticle("Ph1")->setParticleContent(Ph2, cov_Ph2);            CmdBaseParticle& ph2 = *event2Pi2Pi0_6C.GetParticle("Ph1");
			    event2Pi2Pi0_6C.GetParticle("Ph2")->setParticleContent(Ph3, cov_Ph3);            CmdBaseParticle& ph3 = *event2Pi2Pi0_6C.GetParticle("Ph2");
			    event2Pi2Pi0_6C.GetParticle("Ph3")->setParticleContent(Ph4, cov_Ph4);            CmdBaseParticle& ph4 = *event2Pi2Pi0_6C.GetParticle("Ph3");
			    double mpi0000000 = (ph3.GetMoment()+ph4.GetMoment()).M();
			    Fitter.Init(&event2Pi2Pi0_6C);
			    //Fitter.Print();
			    //getchar();
			    //Fitter.SetDebug(2);
			    Fitter.davay();
			    int err=Fitter.GetErrCode();
			    if (err!=0) continue;
			    if(Fitter.GetXi2() < Chi25C || firsttry == true){
			      if((ph3.GetMoment()+ph4.GetMoment()).M() < Cut_mpi0_before  || (ph3.GetMoment()+ph4.GetMoment()).M() > Cut_mpi0_above)continue;
			      Err_Chi25C=err;
			      firsttry = false;	
			      Chi25C = Fitter.GetXi2();
			      m2g = (ph3.GetMoment()+ph4.GetMoment()).M();	    	    
			      m2g0 = mpi0000000;
			      Lpip = ppi.GetMoment();
			      Lpim = mpi.GetMoment(); 
			      Lph1 = ph1.GetMoment();
			      Lph2 = ph2.GetMoment();
			      Lph3 = ph3.GetMoment();
			      Lph4 = ph4.GetMoment();
			      m2pirest = (Pebeam - Lpip - Lpim).M();
			      correspond[0] = one;
			      correspond[1] = two;  
			      correspond[2] = k;
			      correspond[3] = l;	
			      correspond[4] = m;
			      correspond[5] = n;
			      E_poln = fabs((PTr[one] + PTr[two] + PPh[k] + PPh[l]+ PPh[m]+ PPh[n]).E() - Pebeam.E());
			      P_poln = fabs((PTr[one] + PTr[two] + PPh[k] + PPh[l]+ PPh[m]+ PPh[n]).P());
			      mode1 = 1;
			      N_events++;
			    }
			  }
		  }
      }
      // end end end end end end end end end end 
      //      2pi0 pi+ pi- hypothesis
      // end end end end end end end end end end

      //****************************************************
      // 2pi0 + 1 track hypothesis
      //****************************************************

      mass = 20;
      for(int one = 0; one < nt; one++){
	  if(fabs(tz[one]) > 12. || fabs(trho[one]) > 0.2 || tdedx[one] > tdedx_high(tptot[one]) || tptot[one] < 50.){continue;} 
	  if(tth[one] < th_tr_min || tth[one] > th_tr_max)continue;
	  for(int j = 0; j < npi0; j++){
			int k = numbers[j][0];
			int l = numbers[j][1];
			int m = numbers[j][2];
			int n = numbers[j][3];
			TLorentzVector miss = Pebeam - PTr[one] - PPh[k] - PPh[l] - PPh[m] - PPh[n];
			if(PPh[k].Theta() > th_ph_max || PPh[l].Theta() > th_ph_max)continue;
			if(PPh[m].Theta() > th_ph_max || PPh[n].Theta() > th_ph_max)continue;	
			if(PPh[k].Theta() < th_ph_min || PPh[l].Theta() < th_ph_min )continue;
			if(PPh[m].Theta() < th_ph_min  || PPh[n].Theta()  < th_ph_min )continue;			
			if(PPh[k].E() < 30. ||PPh[l].E() < 30. || PPh[m].E() < 30. || PPh[n].E() < 30.)continue;
	    		Pip(0)= tptot[one]*sin(tth[one]);Pip(1)= cos(tth[one])/sin(tth[one]);Pip(2)= tphi[one];
	    	        Pip(3) = trho[one];Pip(4) = tz[one];
	    		Pim(0)= miss.P()*sin(miss.Theta());Pim(1)= cos(miss.Theta())/sin(miss.Theta());Pim(2)= miss.Phi();Pim(3) = Pip(3);Pim(4) = Pip(4);
	    		cov_Pip(0,0)= terr0[one][0][0];cov_Pip(1,0)= terr0[one][3][0];cov_Pip(2,0)= terr0[one][1][0];cov_Pip(3,0)=terr0[one][2][0];cov_Pip(4,0)= terr0[one][4][0];
	    		cov_Pip(0,1)= terr0[one][0][3];cov_Pip(1,1)= terr0[one][3][3];cov_Pip(2,1)=terr0[one][1][3];cov_Pip(3,1)=terr0[one][2][3];cov_Pip(4,1)= terr0[one][4][3];
	    		cov_Pip(0,2)= terr0[one][0][1];cov_Pip(1,2)= terr0[one][3][1];cov_Pip(2,2)=terr0[one][1][1];cov_Pip(3,2)=terr0[one][2][1];cov_Pip(4,2)= terr0[one][4][1];
	    		cov_Pip(0,3)= terr0[one][0][2];cov_Pip(1,3)= terr0[one][3][2];cov_Pip(2,3)=terr0[one][1][2];cov_Pip(3,3)=terr0[one][2][2];cov_Pip(4,3)= terr0[one][4][2];
	    		cov_Pip(0,4)= terr0[one][0][4];cov_Pip(1,4)= terr0[one][3][4];cov_Pip(2,4)=terr0[one][1][4];cov_Pip(3,4)=terr0[one][2][4];cov_Pip(4,4)= terr0[one][4][4];
			for(int i = 0; i < 5; i++){for(int zz = 0; zz < 5; zz++){
				cov_Pim(i,zz)=1000.;
				if(i!=zz){cov_Pim(i,zz)=0.;}
			}}
	    		double R = phrho[k];

            		sigma2_rho = 0.3;
            		sigma2_z = 0.3;
	    		if(phflag[k] == 3){sigma2_rho = sigma2_z*pow(tan(phth[k]),2) + pow(phrho[k]*pherr[k][1]/sin(phth[k])*cos(phth[k]),2);}
  	    		else {sigma2_z = sigma2_rho/pow(tan(phth[k]),2) + pow(phrho[k]*pherr[n][1],2)/pow(sin(phth[k]),4);}
	    		Ph1(0)=phen[k];Ph1(1)=R;Ph1(2)=phphi[k]; Ph1(3)=z0+R*cos(phth[k])/sin(phth[k]);Ph1(4)=phth[k];
			cov_Ph1(0,0) = 1./(pherr[k][0]*pherr[k][0]);cov_Ph1(1,1) = 1./sigma2_rho;cov_Ph1(2,2) = 1./(pherr[k][2]*pherr[k][2]);
	    		cov_Ph1(3,3) = 1./sigma2_z;cov_Ph1(4,4) = 1./(pherr[k][1]*pherr[k][1]);

			R = phrho[l];
			Ph2(0)=phen[l];Ph2(1)=R;Ph2(2)=phphi[l];Ph2(3)=z0+R*cos(phth[l])/sin(phth[l]);Ph2(4)=phth[l];
        		sigma2_rho = 0.3;
            		sigma2_z = 0.3;
	    		if(phflag[l] == 3){sigma2_rho = sigma2_z*pow(tan(phth[l]),2) + pow(phrho[l]*pherr[l][1]/sin(phth[l])*cos(phth[l]),2);}
  	    		else {sigma2_z = sigma2_rho/pow(tan(phth[l]),2) + pow(phrho[l]*pherr[l][1],2)/pow(sin(phth[l]),4);}
	    		cov_Ph2(0,0) = 1./(pherr[l][0]*pherr[l][0]);cov_Ph2(1,1) = 1./sigma2_rho;cov_Ph2(2,2) = 1./(pherr[l][2]*pherr[l][2]);
	    		cov_Ph2(3,3) = 1./sigma2_z;cov_Ph2(4,4) = 1./(pherr[l][1]*pherr[l][1]);

	    		R = phrho[m];
	    		Ph3(0)=phen[m];Ph3(1)=R;Ph3(2)=phphi[m];Ph3(3)=z0+R*cos(phth[m])/sin(phth[m]);Ph3(4)=phth[m];
        		sigma2_rho = 0.3;
            		sigma2_z = 0.3;
	    		if(phflag[m] == 3){sigma2_rho = sigma2_z*pow(tan(phth[m]),2) + pow(phrho[m]*pherr[m][1]/sin(phth[m])*cos(phth[m]),2);}
  	    		else {sigma2_z = sigma2_rho/pow(tan(phth[m]),2) + pow(phrho[m]*pherr[m][1],2)/pow(sin(phth[m]),4);}
			cov_Ph3(0,0) = 1./(pherr[m][0]*pherr[m][0]);
	    		cov_Ph3(1,1) = 1./sigma2_rho;cov_Ph3(2,2) = 1./(pherr[m][2]*pherr[m][2]);cov_Ph3(3,3) = 1./sigma2_z;cov_Ph3(4,4) = 1./(pherr[m][1]*pherr[m][1]);

			R = phrho[n];
			Ph4(0)=phen[n];Ph4(1)=R;Ph4(2)=phphi[n];Ph4(3)=z0+R*cos(phth[n])/sin(phth[n]);Ph4(4)=phth[n];   
        		sigma2_rho = 0.3;
            		sigma2_z = 0.3;
	    		if(phflag[n] == 3){sigma2_rho = sigma2_z*pow(tan(phth[n]),2) + pow(phrho[n]*pherr[n][1]/sin(phth[n])*cos(phth[n]),2);}
  	    		else {sigma2_z = sigma2_rho/pow(tan(phth[n]),2) + pow(phrho[n]*pherr[n][1],2)/pow(sin(phth[n]),4);} 		
	    		cov_Ph4(0,0) = 1./(pherr[n][0]*pherr[n][0]);cov_Ph4(1,1) = 1./sigma2_rho;cov_Ph4(2,2) = 1./(pherr[n][2]*pherr[n][2]);
	    		cov_Ph4(3,3) = 1./sigma2_z;cov_Ph4(4,4) = 1./(pherr[n][1]*pherr[n][1]);
	    
	    		// cmd3::CmdChargeParticle ppi("Pi Plus",Mpi,Pip,cov_Pip.Invert());
	    		// cmd3::CmdChargeParticle mpi("Pi Minus",Mpi,Pim,cov_Pim.Invert());
	    		// cmd3::CmdPhotonParticle ph1("Photon 1",Ph1,cov_Ph1);
	    		// cmd3::CmdPhotonParticle ph2("Photon 2",Ph2,cov_Ph2);
	    		// cmd3::CmdPhotonParticle ph3("Photon 3",Ph3,cov_Ph3);
	    		// cmd3::CmdPhotonParticle ph4("Photon 4",Ph4,cov_Ph4);
	    		//event4piw.PutParticle(ppi,mpi,ph1,ph2,ph3,ph4);  
	    		//Fitter.Init();
			
			event4piw.GetParticle("PiPl")->setParticleContent(Pip, cov_Pip.Invert());  CmdBaseParticle& ppi= *event4piw.GetParticle("PiPl"); 
			event4piw.GetParticle("PiMi")->setParticleContent(Pim, cov_Pim.Invert());  CmdBaseParticle& mpi= *event4piw.GetParticle("PiMi");
			event4piw.GetParticle("Ph0")->setParticleContent(Ph1, cov_Ph1);            CmdBaseParticle& ph1 = *event4piw.GetParticle("Ph0");
			event4piw.GetParticle("Ph1")->setParticleContent(Ph2, cov_Ph2);            CmdBaseParticle& ph2 = *event4piw.GetParticle("Ph1");
			event4piw.GetParticle("Ph2")->setParticleContent(Ph3, cov_Ph3);            CmdBaseParticle& ph3 = *event4piw.GetParticle("Ph2");
			event4piw.GetParticle("Ph3")->setParticleContent(Ph4, cov_Ph4);            CmdBaseParticle& ph4 = *event4piw.GetParticle("Ph3");
			
			Fitter.Init(&event4piw); 
 
	    		Fitter.SetDebug(0);
	    		Fitter.davay();

			if(fabs(miss.M()-Mpi) < fabs(mass - Mpi)){
				mode2 = 1;
				mass = (Pebeam - PTr[one] - PPh[k] - PPh[l] - PPh[m] - PPh[n]).M();
 				mass34ph1tr = (ph3.GetMoment()+ph4.GetMoment()).M();
				theta_miss_pi = (Pebeam - ppi.GetMoment() - ph1.GetMoment() - ph2.GetMoment() - ph3.GetMoment() - ph4.GetMoment()).Theta();
				mom_miss_pi = (Pebeam - ppi.GetMoment() - ph1.GetMoment() - ph2.GetMoment() - ph3.GetMoment() - ph4.GetMoment()).P();
				correspond[6] = m;
				correspond[7] = n;
				correspond[8] = k;
				correspond[9] = l;
				correspond[10] = one;
			}
		  }
	  }
 
 
      //**************************************************
      //    1 pi0 pi+ pi- hypothesis
      //**************************************************
      mass = 300.*ebeam/500.;
      Xi2misspi0 = 1000.;
      for(int one = 0; one < nt-1; one++){
	  	  for(int two = one + 1; two < nt; two++){
			  if(tcharge[one]*tcharge[two] > 0){continue;}
			  if(nt > 4)continue;
			  if(nph > 8)continue;
	  	  	  if(fabs(tz[one]) > 12. || fabs(trho[one]) > 0.2 || tdedx[one] > tdedx_high(tptot[one])){continue;}
	  	  	  if(fabs(tz[two]) > 12. || fabs(trho[two]) > 0.2 || tdedx[two] > tdedx_high(tptot[two])){continue;} 
	  	  	  for(int l = 0; l < npi0_alone; l++){
	  	  	  	  int m = numbers_alone[l][0];
	  	  	  	  int n = numbers_alone[l][1];
				  if(fabs((Pebeam - PTr[one] - PTr[two] - PPh[m] - PPh[n]).M() - Mpi) < mass){
					mass = fabs((Pebeam - PTr[one] - PTr[two] - PPh[m] - PPh[n]).M() - Mpi);
					misspi0mass = (Pebeam - PTr[one] - PTr[two] - PPh[m] - PPh[n]).M();
                                        TLorentzVector pi0miss = Pebeam - PTr[one] - PTr[two] - PPh[m] - PPh[n];
					Pip(0)= tptot[one]*sin(tth[one]);
	    				Pip(1)= cos(tth[one])/sin(tth[one]);
	    				Pip(2)= tphi[one];
	    				Pip(3) = trho[one];
	    				Pip(4) = tz[one];
	    				Pim(0)= tptot[two]*sin(tth[two]);
	    				Pim(1)= cos(tth[two])/sin(tth[two]);
	    				Pim(2)= tphi[two];
					Pim(3) = trho[two];
	    				Pim(4) = tz[two];
	    
	    				cov_Pip(0,0)= terr0[one][0][0];cov_Pip(1,0)= terr0[one][3][0];cov_Pip(2,0)= terr0[one][1][0];cov_Pip(3,0)=terr0[one][2][0];cov_Pip(4,0)= terr0[one][4][0];
	    				cov_Pip(0,1)= terr0[one][0][3];cov_Pip(1,1)= terr0[one][3][3];cov_Pip(2,1)=terr0[one][1][3];cov_Pip(3,1)=terr0[one][2][3];cov_Pip(4,1)= terr0[one][4][3];
	    				cov_Pip(0,2)= terr0[one][0][1];cov_Pip(1,2)= terr0[one][3][1];cov_Pip(2,2)=terr0[one][1][1];cov_Pip(3,2)=terr0[one][2][1];cov_Pip(4,2)= terr0[one][4][1];
	    				cov_Pip(0,3)= terr0[one][0][2];cov_Pip(1,3)= terr0[one][3][2];cov_Pip(2,3)=terr0[one][1][2];cov_Pip(3,3)=terr0[one][2][2];cov_Pip(4,3)= terr0[one][4][2];
	    				cov_Pip(0,4)= terr0[one][0][4];cov_Pip(1,4)= terr0[one][3][4];cov_Pip(2,4)=terr0[one][1][4];cov_Pip(3,4)=terr0[one][2][4];cov_Pip(4,4)= terr0[one][4][4];

            				cov_Pim(0,0)= terr0[two][0][0];cov_Pim(1,0)= terr0[two][3][0];cov_Pim(2,0)= terr0[two][1][0];cov_Pim(3,0)=terr0[two][2][0];cov_Pim(4,0)= terr0[two][4][0];
	    				cov_Pim(0,1)= terr0[two][0][3];cov_Pim(1,1)= terr0[two][3][3];cov_Pim(2,1)=terr0[two][1][3];cov_Pim(3,1)=terr0[two][2][3];cov_Pim(4,1)= terr0[two][4][3];
	    				cov_Pim(0,2)= terr0[two][0][1];cov_Pim(1,2)= terr0[two][3][1];cov_Pim(2,2)=terr0[two][1][1];cov_Pim(3,2)=terr0[two][2][1];cov_Pim(4,2)= terr0[two][4][1];
	    				cov_Pim(0,3)= terr0[two][0][2];cov_Pim(1,3)= terr0[two][3][2];cov_Pim(2,3)=terr0[two][1][2];cov_Pim(3,3)=terr0[two][2][2];cov_Pim(4,3)= terr0[two][4][2];
	    				cov_Pim(0,4)= terr0[two][0][4];cov_Pim(1,4)= terr0[two][3][4];cov_Pim(2,4)=terr0[two][1][4];cov_Pim(3,4)=terr0[two][2][4];cov_Pim(4,4)= terr0[two][4][4];
					double R = phrho[m];
	    				Ph1(0)=phen[m];
	    				Ph1(1)=R;
	    				Ph1(2)=phphi[m];
	    				Ph1(3)= z0 + R*cos(phth[m])/sin(phth[m]);
	    				Ph1(4)=phth[m];
	    				R = phrho[n];
	    				Ph2(0)=phen[n];
	    				Ph2(1)=R;
	    				Ph2(2)=phphi[n];
	    				Ph2(3)=z0 + R*cos(phth[n])/sin(phth[n]);
	    				Ph2(4)=phth[n];
            				sigma2_rho = 0.3;
            				sigma2_z = 0.3;
	    				if(phflag[m] == 3){sigma2_rho = sigma2_z*pow(tan(phth[m]),2) + pow(phrho[m]*pherr[m][1]/sin(phth[m])*cos(phth[m]),2);}
  	    				else {sigma2_z = sigma2_rho/pow(tan(phth[m]),2) + pow(phrho[m]*pherr[m][1],2)/pow(sin(phth[m]),4);}

	    				cov_Ph1(0,0) = 1./(pherr[m][0]*pherr[m][0]);
	    				cov_Ph1(1,1) = 1./sigma2_rho;
	    				cov_Ph1(2,2) = 1./(pherr[m][2]*pherr[m][2]);
	    				cov_Ph1(3,3) = 1./sigma2_z;
	    				cov_Ph1(4,4) = 1./(pherr[m][1]*pherr[m][1]);

            				sigma2_rho = 0.3;
            				sigma2_z = 0.3;
	    				if(phflag[n] == 3){sigma2_rho = sigma2_z*pow(tan(phth[n]),2) + pow(phrho[n]*pherr[n][1]/sin(phth[n])*cos(phth[n]),2);}
  	    				else {sigma2_z = sigma2_rho/pow(tan(phth[n]),2) + pow(phrho[n]*pherr[n][1],2)/pow(sin(phth[n]),4);}

	    				cov_Ph2(0,0) = 1./(pherr[n][0]*pherr[n][0]);
	    				cov_Ph2(1,1) = 1./sigma2_rho;
	    				cov_Ph2(2,2) = 1./(pherr[n][2]*pherr[n][2]);
	    				cov_Ph2(3,3) = 1./sigma2_z;
	    				cov_Ph2(4,4) = 1./(pherr[n][1]*pherr[n][1]);

					event4pilossenpi0.GetParticle("PiPl")->setParticleContent(Pip, cov_Pip.Invert());  CmdBaseParticle& ppi= *event4pilossenpi0.GetParticle("PiPl");
					event4pilossenpi0.GetParticle("PiMi")->setParticleContent(Pim, cov_Pim.Invert());  CmdBaseParticle& mpi= *event4pilossenpi0.GetParticle("PiMi");
					event4pilossenpi0.GetParticle("Ph0")->setParticleContent(Ph1, cov_Ph1);            CmdBaseParticle& ph1 = *event4pilossenpi0.GetParticle("Ph0");
					event4pilossenpi0.GetParticle("Ph1")->setParticleContent(Ph2, cov_Ph2);            CmdBaseParticle& ph2 = *event4pilossenpi0.GetParticle("Ph1");
					double massbtv = (ph1.GetMoment()+ph2.GetMoment()).M();
					Fitter.Init(&event4pilossenpi0);
					Fitter.SetDebug(0);
					Fitter.davay();
					if(Fitter.GetErrCode()!=0)continue;
					if(Fitter.GetXi2() < Xi2misspi0){
					  if((ph1.GetMoment()+ph2.GetMoment()).M() < 70 || (ph1.GetMoment()+ph2.GetMoment()).M()  > 230 || ph1.GetMoment().E() < 40 || ph2.GetMoment().E() < 40)continue;
					  if(fabs(1.57-ph1.GetMoment().Theta()) > 0.87)continue;
                                          if(fabs(1.57-ph2.GetMoment().Theta()) > 0.87)continue;
                                          if(fabs(1.57-ppi.GetMoment().Theta()) > 0.77)continue;
                                          if(fabs(1.57-mpi.GetMoment().Theta()) > 0.77)continue;
					  mass34ph1pi0 = massbtv;
					  mass34ph1pi0_kin =  (ph1.GetMoment()+ph2.GetMoment()).M();
					  Xi2misspi0 = Fitter.GetXi2();
					  correspond[11] = m;
					  correspond[12] = n;
					  correspond[13] = one;
					  correspond[14] = two;
					  Ltr1 = ppi.GetMoment();
					  Ltr2 = mpi.GetMoment();
					  Lgamma1 = ph1.GetMoment();
					  Lgamma2 = ph2.GetMoment();
					  mode3 = 1;
					}
				  }
			  }
		  }
      }
      


       // end end end end end end end end end end end end
       //    1 pi0 pi+ pi- hypothesis
       // end end end end end end end end end end end end
 





      if((mode1 + mode2 + mode3) > 0){newtree->Fill();}
      if((mode1 + mode2 + mode3) == 0 && boolmc == true){newtree->Fill();}

   }


    cout << "OK1  " << namefile << " " << (double)newtree->GetEntries()/(double)nentries << endl;
    cout << "N_events = " << N_events << "    " << (double)newtree->GetEntries() << endl;
    newtree->Write();
    delete newtree;
    newfile.Write();
    newfile.Close();


}














      //***************************************** 
      //   pi0 gamma pi+ pi- hypothesis 
      //*****************************************
      /*	
      Chi_merge = 3000000;
      for(int one = 0; one < nt-1; one++){
	  	  for(int two = one + 1; two < nt; two++){
		          if(tcharge[one]*tcharge[two] > 0)continue;
	  	  	  if(fabs(tz[one]) > 13. || fabs(trho[one]) > 0.2 || tdedx[one] > tdedx_high(tptot[one])){continue;}
	  	  	  if(fabs(tz[two]) > 13. || fabs(trho[two]) > 0.2 || tdedx[two] > tdedx_high(tptot[two])){continue;}
			  Err_merge=10;
			  for(int k = 0; k < nph-1; k++){
			    for(int l = k+1; l < nph; l++){
			       for(int m = 0; m < nph; m++){
			         if(m==k || m ==l || phen[m] < 150.)continue;
                                 
				 Pip(0)= tptot[one]*sin(tth[one]);
				 Pip(1)= cos(tth[one])/sin(tth[one]);
				 Pip(2)= tphi[one];
				 Pip(3) = trho[one];
				 Pip(4) = tz[one];
				 Pim(0)= tptot[two]*sin(tth[two]);
				 Pim(1)= cos(tth[two])/sin(tth[two]);
				 Pim(2)= tphi[two];
				 Pim(3) = trho[two];
				 Pim(4) = tz[two];
				 double Z_photons=0.5*(Pim(4)+Pip(4));
				 //cout<<"Z_Photons = "<<Z_photons<<endl;
				 //getchar();
				 cov_Pip(0,0)= terr0[one][0][0];cov_Pip(1,0)= terr0[one][3][0];cov_Pip(2,0)= terr0[one][1][0];cov_Pip(3,0)=terr0[one][2][0];cov_Pip(4,0)= terr0[one][4][0];
				 cov_Pip(0,1)= terr0[one][0][3];cov_Pip(1,1)= terr0[one][3][3];cov_Pip(2,1)=terr0[one][1][3];cov_Pip(3,1)=terr0[one][2][3];cov_Pip(4,1)= terr0[one][4][3];
				 cov_Pip(0,2)= terr0[one][0][1];cov_Pip(1,2)= terr0[one][3][1];cov_Pip(2,2)=terr0[one][1][1];cov_Pip(3,2)=terr0[one][2][1];cov_Pip(4,2)= terr0[one][4][1];
				 cov_Pip(0,3)= terr0[one][0][2];cov_Pip(1,3)= terr0[one][3][2];cov_Pip(2,3)=terr0[one][1][2];cov_Pip(3,3)=terr0[one][2][2];cov_Pip(4,3)= terr0[one][4][2];
				 cov_Pip(0,4)= terr0[one][0][4];cov_Pip(1,4)= terr0[one][3][4];cov_Pip(2,4)=terr0[one][1][4];cov_Pip(3,4)=terr0[one][2][4];cov_Pip(4,4)= terr0[one][4][4];
				 
				 cov_Pim(0,0)= terr0[two][0][0];cov_Pim(1,0)= terr0[two][3][0];cov_Pim(2,0)= terr0[two][1][0];cov_Pim(3,0)=terr0[two][2][0];cov_Pim(4,0)= terr0[two][4][0];
				 cov_Pim(0,1)= terr0[two][0][3];cov_Pim(1,1)= terr0[two][3][3];cov_Pim(2,1)=terr0[two][1][3];cov_Pim(3,1)=terr0[two][2][3];cov_Pim(4,1)= terr0[two][4][3];
				 cov_Pim(0,2)= terr0[two][0][1];cov_Pim(1,2)= terr0[two][3][1];cov_Pim(2,2)=terr0[two][1][1];cov_Pim(3,2)=terr0[two][2][1];cov_Pim(4,2)= terr0[two][4][1];
				 cov_Pim(0,3)= terr0[two][0][2];cov_Pim(1,3)= terr0[two][3][2];cov_Pim(2,3)=terr0[two][1][2];cov_Pim(3,3)=terr0[two][2][2];cov_Pim(4,3)= terr0[two][4][2];
				 cov_Pim(0,4)= terr0[two][0][4];cov_Pim(1,4)= terr0[two][3][4];cov_Pim(2,4)=terr0[two][1][4];cov_Pim(3,4)=terr0[two][2][4];cov_Pim(4,4)= terr0[two][4][4];
				

				 double R = phrho[k];
				 Ph1(0)=phen[k];
				 Ph1(1)=R;
				 Ph1(2)=phphi[k];
				 Ph1(3)= z0 + R*cos(phth[k])/sin(phth[k]);
				 Ph1(4)=phth[k];Ph1(4)=atan2(Ph1(1),Ph1(3)-Z_photons);
				 R = phrho[l];
				 Ph2(0)=phen[l];
				 Ph2(1)=R;
				 Ph2(2)=phphi[l];
				 Ph2(3)=z0 + R*cos(phth[l])/sin(phth[l]);
				 Ph2(4)=phth[l];Ph2(4)=atan2(Ph2(1),Ph2(3)-Z_photons);
				 sigma2_rho = 0.3;
				 sigma2_z = 0.3;
				 if(phflag[k] == 3){sigma2_rho = sigma2_z*pow(tan(phth[k]),2) + pow(phrho[k]*pherr[k][1]/sin(phth[k])*cos(phth[k]),2);}
				 else {sigma2_z = sigma2_rho/pow(tan(phth[k]),2) + pow(phrho[k]*pherr[k][1],2)/pow(sin(phth[k]),4);}
				 
				 cov_Ph1(0,0) = 1./(pherr[k][0]*pherr[k][0]);
				 cov_Ph1(1,1) = 1./sigma2_rho;
				 cov_Ph1(2,2) = 1./(pherr[k][2]*pherr[k][2]);
				 cov_Ph1(3,3) = 1./sigma2_z;
				 cov_Ph1(4,4) = 1./(pherr[k][1]*pherr[k][1]);
				 
				 sigma2_rho = 0.3;
				 sigma2_z = 0.3;
				 if(phflag[l] == 3){sigma2_rho = sigma2_z*pow(tan(phth[l]),2) + pow(phrho[l]*pherr[l][1]/sin(phth[l])*cos(phth[l]),2);}
				 else {sigma2_z = sigma2_rho/pow(tan(phth[l]),2) + pow(phrho[l]*pherr[l][1],2)/pow(sin(phth[l]),4);}
				 
				 cov_Ph2(0,0) = 1./(pherr[l][0]*pherr[l][0]);
				 cov_Ph2(1,1) = 1./sigma2_rho;
				 cov_Ph2(2,2) = 1./(pherr[l][2]*pherr[l][2]);
				 cov_Ph2(3,3) = 1./sigma2_z;
				 cov_Ph2(4,4) = 1./(pherr[l][1]*pherr[l][1]);
				 
				 R = phrho[m];
				 Ph3(0)=phen[m];
				 Ph3(1)=R;
				 Ph3(2)=phphi[m];
				 Ph3(3)= z0 + R*cos(phth[m])/sin(phth[m]);//tz[one]+
				 Ph3(4)=phth[m];Ph3(4)=atan2(Ph3(1),Ph3(3)-Z_photons);

				 sigma2_rho = 0.3;
				 sigma2_z = 0.3;
				 if(phflag[m] == 3){sigma2_rho = sigma2_z*pow(tan(phth[m]),2) + pow(phrho[m]*pherr[m][1]/sin(phth[m])*cos(phth[m]),2);}
				 else {sigma2_z = sigma2_rho/pow(tan(phth[m]),2) + pow(phrho[m]*pherr[m][1],2)/pow(sin(phth[m]),4);}
			    
				 cov_Ph3(0,0) = 1./(pherr[m][0]*pherr[m][0]);
				 cov_Ph3(1,1) = 1./sigma2_rho;
				 cov_Ph3(2,2) = 1./(pherr[m][2]*pherr[m][2]);
				 cov_Ph3(3,3) = 1./sigma2_z;
				 cov_Ph3(4,4) = 1./(pherr[m][1]*pherr[m][1]);
			    
				 event2pi2gamma1Pi0.GetParticle("PiPl")->setParticleContent(Pip, cov_Pip.Invert());  CmdBaseParticle& ppi= *event2pi2gamma1Pi0.GetParticle("PiPl");
				 event2pi2gamma1Pi0.GetParticle("PiMi")->setParticleContent(Pim, cov_Pim.Invert());  CmdBaseParticle& mpi= *event2pi2gamma1Pi0.GetParticle("PiMi");
				 event2pi2gamma1Pi0.GetParticle("Ph0")->setParticleContent(Ph1, cov_Ph1);            CmdBaseParticle& ph1 = *event2pi2gamma1Pi0.GetParticle("Ph0");
				 event2pi2gamma1Pi0.GetParticle("Ph1")->setParticleContent(Ph2, cov_Ph2);            CmdBaseParticle& ph2 = *event2pi2gamma1Pi0.GetParticle("Ph1");
				 event2pi2gamma1Pi0.GetParticle("Pi0")->setParticleContent(Ph3, cov_Ph3);            CmdBaseParticle& pi0 = *event2pi2gamma1Pi0.GetParticle("Pi0");
				 cout << "before " << (pi0.GetMoment()).M() << " " << (pi0.GetMoment()).E() << " " << (pi0.GetMoment()).Theta() << endl;
				 Fitter.Init(&event2pi2gamma1Pi0);
				 //Fitter.Print();
				 //getchar();
				 Fitter.SetDebug(0);
				 Fitter.davay();
				 int err=Fitter.GetErrCode();
				 if (err!=0) continue;
				 //getchar();
				 if(Fitter.GetXi2() < Chi_merge){
				   Xi2_Ph3=pi0.GetXi2();
				   Err_merge = err;
				   Chi_merge = Fitter.GetXi2();
				   m2g_merge = (ph1.GetMoment() + ph2.GetMoment()).M();
				   correspond[15]=one;
				   correspond[16]=two;
				   correspond[17]=k;
				   correspond[18]=l;
				   correspond[19]=m;
                                   cout << "after " << (pi0.GetMoment()).M() << " " << (pi0.GetMoment()).E() << " " << (pi0.GetMoment()).Theta() << endl;
				   
				 }
			       }
			    }
			  }
		  }
      }
      */
      // *************************************************************************
      // ========================================================================
      // ************************************************************************



