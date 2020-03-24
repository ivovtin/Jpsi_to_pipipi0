//ЧЩЮЙУМЕОЙЕ УТЕДОЕЗП ЛПЬЖЖЙГЙЕОФБ УЧЕФПУВПТБ РП РТБЧП/МЕЧП УЮЕФЮЙЛБ ДМС 5-ФЙ ФЙРПЧ УЮЕФЮЙЛПЧ
//СЂРёСЃРѕРІР°РЅРёРµ С‚СЂРµС…РјРµСЂРЅРѕР№ РіРёСЃС‚РѕРіСЂР°РјРјС‹ x,y,z - РѕС‚СЂРёСЃРѕРІРєР° СЃС‡РµС‚С‡РёРєР° РІ С‚СЂРµС…РјРµСЂРЅРѕРј РІРёРґРµ РІ Р»РѕРєР°Р»СЊРЅС‹С… РєРѕРѕСЂРґРёРЅР°С‚Р°С…
#include <Riostream.h>
#include <sstream>
#include <TROOT.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TTree.h>
#include <TF1.h>
#include <TMath.h>
#include "TFile.h"
#include "TVirtualPad.h"
#include <iomanip>
#pragma hdrstop
#include<stdio.h>
#include<stdlib.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TStyle.h>
#include <TLine.h>
#include <TEventList.h>
#include <TProfile.h>
#include <vector>
#include <TChain.h>
#include "TMinuit.h"
#include "TRandom3.h"
#include <math.h>
#include "TVirtualFitter.h"
#include <algorithm>      // std::min_element, std::max_element
#include <assert.h>

using namespace std;
string progname;

int Usage(string status)
{
        cout<<"Usage: "<<progname<<"\t"<<"Range ATC->1,2...  Hit in aer(0,1)  Hit in wls  Lenght  First cnt->0-160  End cnt->0-160 Momentum(350 or 3000)  Data->Exp or Sim(0,1)"<<endl;
        exit(0);
}

int main(int argc, char* argv[])
{  
   progname=argv[0];
    int region;
    int hit_aer;
    int hit_wls;
    int len; 
    int first_cnt;
    int end_cnt;
    int momentum;  
    int fl_sim_exp;
    if( argc>1 )
    {
    region=atoi(argv[1]);
    hit_aer=atoi(argv[2]);
    hit_wls=atoi(argv[3]);
    len=atoi(argv[4]);
    first_cnt=atoi(argv[5]);  
    end_cnt=atoi(argv[6]);
    momentum=atoi(argv[7]);
    fl_sim_exp=atoi(argv[8]);   
    if(region>10){ Usage(progname); return 0;}  
    if(hit_aer>1 || hit_wls>1 || len>4 || fl_sim_exp>1){ Usage(progname); return 0;} 
    if(hit_aer<0 || hit_wls<0 || len<0 || fl_sim_exp<0){ Usage(progname); return 0;} 
    if(first_cnt<0 || first_cnt>160){ Usage(progname); return 0;}  
    if(end_cnt<0 || first_cnt>160){ Usage(progname); return 0;}
//    if(momentum!=350 || momentum!=3000){ Usage(progname); return 0;}
    }
    else
    { Usage(progname);
    }

//if (cnt.phiin!=0 && cnt.phiout!=0 && cnt.aerogel_region0==1 && cnt.wlshit!=1 && cnt.nearwls!=1 ) {
//if ( cnt.phiin!=0 && cnt.phiout!=0 && cnt.single_aerogel_region0!=1 && cnt.wlshit==1 && cnt.nearwls==1 ) {

struct data                                  //СЏв”‚СЏв”ЊСЏв”ЂСЏв”ђРїв•ЁСЏв”ЊСЏв”ЂСЏв”ђРїв•џ СЏв”‚ Рїв•ўРїв•џРїв•«Рїв•«СЏв–ЂРїв•ЄРїв•¦ Рїв–‘Рїв•–Рїв•‘
    {
        int i, t, ntrk, triggered,zero,fitted,estimated,neightrig,wlshit,nearwls,aerogel_region,aerogel_region0,
	aerogel_region5,aerogel_region20,active_region,active_region0,active_region5,active_region20,test,
	single_aerogel_region,single_aerogel_region0,single_aerogel_region5,single_aerogel_region20,single_active_region,
	single_active_region0,single_active_region5,single_active_region20,single_test,
	in_aerogel_region,in_aerogel_region0,in_aerogel_region5,in_aerogel_region20,in_active_region,in_active_region0,in_active_region5,in_active_region20,in_test,
	out_aerogel_region,out_aerogel_region0,out_aerogel_region5,out_aerogel_region20,out_active_region,out_active_region0,out_active_region5,out_active_region20,out_test;
	float amp, rtime,time,chi2,npe,npen,tlen,pathwls,rin,phiin,zin,rout,phiout,zout,rwls,phiwls,zwls,neighnpe,Rin_gl,
	Phiin_gl,Zin_gl,Rout_gl,Phiout_gl,Zout_gl;
    };
data atccr;
data t0atccr0;
data t0atccr1;
data t0atccr2;
data t0atccr3;
data t0atccr4;
data t1atccr0;
data t1atccr1;
data t1atccr2;
data t1atccr3;
data t1atccr4;

struct data2                                 
    { int t,q,ip,nvec,nvecxy,nvecz,nhits,nhitsxy,nhitsz,nhitsvd;
      float p,pt,theta,phi,chi2,rc,xc,yc,zc,za,ph0,ph1,ph2,x0,x1,x2,y0,y1,y2,z0,z1,z2,vx,vy,vz;
      int emc_ncls, atc_ncnt;
    };
//data2 track;
data2 t0;
data2 t1;

struct data3                                 
    { int nhits, dchits, namps, ntimes;
      float time[2], length[2], beta[2], phi[2];
      int type[2];
    };
data3 t0tof;
data3 t1tof;

struct data4                                 
    {
      int ncls,ncls_trk,nlkr,ncsi,nstrcls,nstrtrk;
      float energy,elkr,ecsi;
    };
data4 emc;

struct data5                                 
    {
      int c,lkr,csi,ncells;
      float e,x,y,z,vx,vy,vz,theta,phi,rho,dtheta,dphi,drho,thetastr,phistr;
      int qlty,ncellsbad,str_ncls,str_ntrk,dc_ntrk,emc_ncls;
    };
data5 t0c0;
data5 t0c1;
data5 t0c2;
data5 t1c0;
data5 t1c1;
data5 t1c2;
data5 clgamma0;
data5 clgamma1;
data5 clgamma2;
data5 clgamma3;

struct data7                                 
    {
      int nhits,dcmuhits,octant,layer;
      int status;
    };
data7 mu;

struct data8                                 
    {
      int event,evdaq,run,quality;
      float ebeam;
    };
data8 ev;

struct data9                                 
    {
      	int ntrk, nip, nbeam;
	float x, y, z;
	float sig_x, sig_y, sig_z;
	float theta2t, phi2t;
    };
data9 vrt;

struct data10                                 
    {
      int numHyp;
      float chi2[5],M[5],P1[5],P2[5];
    };
data10 jpsi;


struct data11                               			
	{
	int natc_cr, natc_hits, natc_thr, rawdt;
	float dt;
	};
data11 atcev;

std::vector< float > npe1;
float Kdiff,kx1,kx2,kx3,kx4;

TCanvas *cc1 = new TCanvas();
gStyle->SetOptStat(1111);
gStyle->SetOptFit(1011);
gROOT->SetStyle("Plain");
cc1->cd();

bool sim=fl_sim_exp;
TFile *fout=0;

TString fname;
 if( sim!=1){	
	fname=TString::Format("cnt_%d_%d_%d_%d_%d_exp_jpsito3pi.root",first_cnt,end_cnt,region,len,momentum).Data();
 }else{
	fname=TString::Format("cnt_%d_%d_%d_%d_%d_sim_jpsito3pi.root",first_cnt,end_cnt,region,len,momentum).Data();
 }
cout<<fname<<endl;
fout = new TFile(fname,"RECREATE");

char branchname[1];
char branchname1[161];

TChain *tt=new TChain("et");
if(sim!=1){
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/data_jpsito3pi/jpsi_to_pipipi0_runs_1.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/data_jpsito3pi/jpsi_to_pipipi0_runs_2.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/data_jpsito3pi/jpsi_to_pipipi0_runs_3.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/data_jpsito3pi/jpsi_to_pipipi0_runs_4.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/data_jpsito3pi/jpsi_to_pipipi0_runs_5.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/data_jpsito3pi/jpsi_to_pipipi0_runs_6.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/data_jpsito3pi/jpsi_to_pipipi0_runs_7.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/data_jpsito3pi/jpsi_to_pipipi0_runs_8.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/data_jpsito3pi/jpsi_to_pipipi0_runs_9.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/data_jpsito3pi/jpsi_to_pipipi0_runs_10.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/data_jpsito3pi/jpsi_to_pipipi0_runs_11.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/data_jpsito3pi/jpsi_to_pipipi0_runs_12.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/data_jpsito3pi/jpsi_to_pipipi0_runs_13.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/data_jpsito3pi/jpsi_to_pipipi0_runs_14.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/data_jpsito3pi/jpsi_to_pipipi0_runs_15.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/data_jpsito3pi/jpsi_to_pipipi0_runs_16.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/data_jpsito3pi/jpsi_to_pipipi0_runs_17.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/data_jpsito3pi/jpsi_to_pipipi0_runs_18.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/data_jpsito3pi/jpsi_to_pipipi0_runs_19.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/data_jpsito3pi/jpsi_to_pipipi0_runs_20.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/data_jpsito3pi/jpsi_to_pipipi0_runs_21.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/data_jpsito3pi/jpsi_to_pipipi0_runs_22.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/data_jpsito3pi/jpsi_to_pipipi0_runs_23.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/data_jpsito3pi/jpsi_to_pipipi0_runs_24.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/data_jpsito3pi/jpsi_to_pipipi0_runs_25.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/data_jpsito3pi/jpsi_to_pipipi0_runs_26.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/data_jpsito3pi/jpsi_to_pipipi0_runs_27.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/data_jpsito3pi/jpsi_to_pipipi0_runs_28.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/data_jpsito3pi/jpsi_to_pipipi0_runs_29.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/data_jpsito3pi/jpsi_to_pipipi0_runs_30.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/data_jpsito3pi/jpsi_to_pipipi0_runs_31.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/data_jpsito3pi/jpsi_to_pipipi0_runs_32.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/data_jpsito3pi/jpsi_to_pipipi0_runs_33.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/data_jpsito3pi/jpsi_to_pipipi0_runs_34.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/data_jpsito3pi/jpsi_to_pipipi0_runs_35.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/data_jpsito3pi/jpsi_to_pipipi0_runs_36.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/data_jpsito3pi/jpsi_to_pipipi0_runs_37.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/data_jpsito3pi/jpsi_to_pipipi0_runs_38.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/data_jpsito3pi/jpsi_to_pipipi0_runs_39.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/data_jpsito3pi/jpsi_to_pipipi0_runs_40.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/data_jpsito3pi/jpsi_to_pipipi0_runs_41.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/data_jpsito3pi/jpsi_to_pipipi0_runs_42.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/data_jpsito3pi/jpsi_to_pipipi0_runs_43.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/data_jpsito3pi/jpsi_to_pipipi0_runs_44.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/data_jpsito3pi/jpsi_to_pipipi0_runs_45.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/data_jpsito3pi/jpsi_to_pipipi0_runs_46.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/data_jpsito3pi/jpsi_to_pipipi0_runs_47.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/data_jpsito3pi/jpsi_to_pipipi0_runs_48.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/data_jpsito3pi/jpsi_to_pipipi0_runs_49.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/data_jpsito3pi/jpsi_to_pipipi0_runs_50.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/data_jpsito3pi/jpsi_to_pipipi0_runs_51.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/data_jpsito3pi/jpsi_to_pipipi0_runs_52.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/data_jpsito3pi/jpsi_to_pipipi0_runs_53.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/data_jpsito3pi/jpsi_to_pipipi0_runs_54.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/data_jpsito3pi/jpsi_to_pipipi0_runs_55.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/data_jpsito3pi/jpsi_to_pipipi0_runs_56.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/data_jpsito3pi/jpsi_to_pipipi0_runs_57.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/data_jpsito3pi/jpsi_to_pipipi0_runs_58.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/data_jpsito3pi/jpsi_to_pipipi0_runs_59.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/data_jpsito3pi/jpsi_to_pipipi0_runs_60.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/data_jpsito3pi/jpsi_to_pipipi0_runs_61.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/data_jpsito3pi/jpsi_to_pipipi0_runs_62.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/data_jpsito3pi/jpsi_to_pipipi0_runs_63.root");
}
else{
/*
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/sim_jpsito3pi/jpsito3pi_sim_24.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/sim_jpsito3pi/jpsito3pi_sim_25.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/sim_jpsito3pi/jpsito3pi_sim_26.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/sim_jpsito3pi/jpsito3pi_sim_27.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/sim_jpsito3pi/jpsito3pi_sim_28.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/sim_jpsito3pi/jpsito3pi_sim_29.root");
*/
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/sim_jpsito3pi/jpsito3pi_sim_31.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/sim_jpsito3pi/jpsito3pi_sim_32.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/sim_jpsito3pi/jpsito3pi_sim_33.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/sim_jpsito3pi/jpsito3pi_sim_34.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/sim_jpsito3pi/jpsito3pi_sim_35.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/sim_jpsito3pi/jpsito3pi_sim_36.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/sim_jpsito3pi/jpsito3pi_sim_37.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/sim_jpsito3pi/jpsito3pi_sim_38.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/sim_jpsito3pi/jpsito3pi_sim_39.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/sim_jpsito3pi/jpsito3pi_sim_40.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/sim_jpsito3pi/jpsito3pi_sim_41.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/sim_jpsito3pi/jpsito3pi_sim_42.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/sim_jpsito3pi/jpsito3pi_sim_43.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/sim_jpsito3pi/jpsito3pi_sim_44.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/sim_jpsito3pi/jpsito3pi_sim_45.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/sim_jpsito3pi/jpsito3pi_sim_46.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/sim_jpsito3pi/jpsito3pi_sim_47.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/sim_jpsito3pi/jpsito3pi_sim_48.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/sim_jpsito3pi/jpsito3pi_sim_49.root");
tt->Add("/home/ovtin/development/ATC/Jpsi_to_3pi/sim_jpsito3pi/jpsito3pi_sim_50.root");
}


	Int_t nentr=tt->GetEntries();               //СЏв”¤Рїв•¦СЏв”‚Рїв•©Рїв•¬ СЏв”‚Рїв•¬Рїв• СЏв–ЂСЏв”ЊРїв•¦Рїв•§ Рїв•Ў Рїв•ўРїв•ЈСЏв”ЂРїв•ЈРїв•ЎРїв•Ј
	cout<<"Nentries="<<nentr<<endl;

	tt->SetBranchStatus("*",0);
	tt->SetBranchStatus("ev",1);
	tt->SetBranchStatus("vrt",1);
	tt->SetBranchStatus("emc",1);
	tt->SetBranchStatus("atcev",1);
	tt->SetBranchStatus("t0",1);
	tt->SetBranchStatus("t0atccr0",1);
	tt->SetBranchStatus("t0atccr1",1);
	tt->SetBranchStatus("t0atccr2",1);
	tt->SetBranchStatus("t0atccr3",1);
	tt->SetBranchStatus("t0atccr4",1);
	tt->SetBranchStatus("t0c0",1);
	tt->SetBranchStatus("t0c1",1);
	tt->SetBranchStatus("t0c2",1);
	tt->SetBranchStatus("t0tof",1);
	tt->SetBranchStatus("t1",1);
	tt->SetBranchStatus("t1atccr0",1);
	tt->SetBranchStatus("t1atccr1",1);
	tt->SetBranchStatus("t1atccr2",1);
	tt->SetBranchStatus("t1atccr3",1);
	tt->SetBranchStatus("t1atccr4",1);
	tt->SetBranchStatus("t1c0",1);
	tt->SetBranchStatus("t1c1",1);
	tt->SetBranchStatus("t1c2",1);
	tt->SetBranchStatus("t1tof",1);
	tt->SetBranchStatus("clgamma0",1);
	tt->SetBranchStatus("clgamma1",1);	
	tt->SetBranchStatus("clgamma2",1);
	tt->SetBranchStatus("clgamma3",1);	
	tt->SetBranchStatus("mu",1);
	tt->SetBranchStatus("jpsi",1);
	
	tt->SetBranchAddress("ev",&ev); 
	tt->SetBranchAddress("vrt",&vrt); 
	tt->SetBranchAddress("emc",&emc);  
	tt->SetBranchAddress("atcev",&atcev);  
	tt->SetBranchAddress("t0",&t0);              
	tt->SetBranchAddress("t0atccr0",&t0atccr0);              
	tt->SetBranchAddress("t0atccr1",&t0atccr1);              
	tt->SetBranchAddress("t0atccr2",&t0atccr2);              
	tt->SetBranchAddress("t0atccr3",&t0atccr3);              
	tt->SetBranchAddress("t0atccr4",&t0atccr4); 
	tt->SetBranchAddress("t0c0",&t0c0); 
	tt->SetBranchAddress("t0c1",&t0c1); 
	tt->SetBranchAddress("t0c2",&t0c2); 
	tt->SetBranchAddress("t0tof",&t0tof); 
	tt->SetBranchAddress("t1",&t1);              
	tt->SetBranchAddress("t1atccr0",&t1atccr0);              
	tt->SetBranchAddress("t1atccr1",&t1atccr1);              
	tt->SetBranchAddress("t1atccr2",&t1atccr2);              
	tt->SetBranchAddress("t1atccr3",&t1atccr3);              
	tt->SetBranchAddress("t1atccr4",&t1atccr4); 
	tt->SetBranchAddress("t1c0",&t1c0); 
	tt->SetBranchAddress("t1c1",&t1c1); 
	tt->SetBranchAddress("t1c2",&t1c2); 
	tt->SetBranchAddress("t1tof",&t1tof); 
	tt->SetBranchAddress("clgamma0",&clgamma0); 
	tt->SetBranchAddress("clgamma1",&clgamma1); 
	tt->SetBranchAddress("clgamma2",&clgamma2); 
	tt->SetBranchAddress("clgamma3",&clgamma3); 
	tt->SetBranchAddress("mu",&mu); 
	tt->SetBranchAddress("jpsi",&jpsi); 
	
	char namepr[0], namepr1[0], namepr2[0];
        if(sim==1)sprintf(namepr,"Simulation");
        if(sim!=1)sprintf(namepr,"Experiment");        
        TProfile* prthink=new TProfile(namepr,namepr,50,0,1600,0,500); 
        if(sim==1)sprintf(namepr1,"Simulation - 1 layer");
        if(sim!=1)sprintf(namepr1,"Experiment - 1 layer");
        TProfile* pr1=new TProfile(namepr1,namepr1,50,0,1600,0,500); 
        if(sim==1)sprintf(namepr2,"Simulation - 2 layer");
        if(sim!=1)sprintf(namepr2,"Experiment - 2 layer");
        TProfile* pr2=new TProfile(namepr2,namepr2,50,0,1600,0,500); 
        //TProfile* pr3=new TProfile("N_{ph.e.}&Momentum cut3","N_{ph.e.}&Momentum",100,0,1600,0,50); 
        TH1F* h1=new TH1F("E/p","E/p",100,0.,10.); 
        TH1F* h2=new TH1F("Energy","Energy",1000,0.,5000.); 
        TH1F* h3=new TH1F("E_LKr","Energy LKr",1000,0.,5000.); 
        TH1F* h4=new TH1F("E_CsI","Energy CsI",1000,0.,5000.); 
        TH1F* h5=new TH1F("E_q-","Energy LKr for q-",1000,0.,5000.); 
        TH1F* h6=new TH1F("E_q+","Energy LKr for q+",1000,0.,5000.);                
        TH1F* h7=new TH1F("Momentum","Momentum",1000,0.,5000.); 
        TH1F* h8=new TH1F("theta","t.theta",1000,0.,185.); 
        TH1F* h9=new TH1F("phi","t.phi",1000,0.,380.);                 
        TH1F* h10=new TH1F("theta2t","vrt.theta2t",1000,0.,185.); 
        TH1F* h11=new TH1F("phi2t","vrt.phi2t",1000,0.,185.);                 
        TH1F* h12=new TH1F("cos(theta2t)","cos(theta2t)",100,-1.,1.); 
        TH1F* h13=new TH1F("cos(phi2t)","cos(phi2t)",100,-1.,1.); 
        TH1F* h14=new TH1F("cos(t.theta)","cos(t.theta)",100,-1.,1.); 
        TH1F* h15=new TH1F("cos(t.phi)","cos(t.phi)",100,-1.,1.); 
        TH1F* h16=new TH1F("RecoilMass","RecoilMass",500,0.,2000.); 
        //TH1F* h16=new TH1F("h16","RecoilMass",10000,-1000000.,2000000.); 
        TH1F* h17=new TH1F("InvMass2pi","InvMass(#pi#pi)",200,0.,4000.); 
        TH1F* h18=new TH1F("cospi0","cos(v(#pi^0))",100,-1.,1.); 
        
        TH1F* h19=new TH1F("t0tof.nhits","t0tof.nhits",20,0.,20.); 
        TH1F* h20=new TH1F("t0tof.dchits","t0tof.dchits",20,0.,20.); 
        TH1F* h21=new TH1F("t0tof.namps","t0tof.namps",20,0.,20.); 
        TH1F* h22=new TH1F("t0tof.ntimes","t0tof.ntimes",20,0.,20.); 
        TH1F* h23=new TH1F("t0tof.time","t0tof.time",1000,-150.,150.); 
        TH1F* h24=new TH1F("t0tof.beta","t0tof.beta",1000,-1.,5.); 
        TH1F* h25=new TH1F("t0tof.length","t0tof.length",1000,0.,250.); 

        TH1F* h26=new TH1F("t1tof.nhits","t1tof.nhits",20,0.,20.); 
        TH1F* h27=new TH1F("t1tof.dchits","t1tof.dchits",20,0.,20.); 
        TH1F* h28=new TH1F("t1tof.namps","t1tof.namps",20,0.,20.); 
        TH1F* h29=new TH1F("t1tof.ntimes","t1tof.ntimes",20,0.,20.); 
        TH1F* h30=new TH1F("t1tof.time","t1tof.time",1000,-150.,150.); 
        TH1F* h31=new TH1F("t1tof.beta","t1tof.beta",1000,-1.,5.); 
        TH1F* h32=new TH1F("t1tof.length","t1tof.length",1000,0.,250.); 
        
        TH1F* h33=new TH1F("ratiop1t0p","p1/t0.p",100,0.,10.); 
        TH1F* h34=new TH1F("ratiop2t1p","p2/t1.p",100,0.,10.); 
        TH1F* h35=new TH1F("S","S",150,0.,1.5); 
        TH1F* h36=new TH1F("invMass2pi","#pi#pi Mass",200,0.,4000.); 
        TH1F* h39=new TH1F("pi0Mass","#pi^0 Mass",3000,0.,3000.); 
        TH1F* h40=new TH1F("t0chi2","t0.chi2",1000,0.,1000.); 
        TH1F* h41=new TH1F("t1chi2","t1.chi2",1000,0.,1000.); 
        TH1F* h42=new TH1F("h42","cos_2t",100,-1.,1.); 
        TH1F* h43=new TH1F("ncls","emc.ncls",12,0.,12.); 
        TH1F* h44=new TH1F("ratiop1t0pp2t1p","(p1/t0.p)/(p2/t1.p)",100,0.,10.); 
        TH1F* h45=new TH1F("t0.nhitsxy","t0.nhitsxy",100,0.,100.); 
        TH1F* h46=new TH1F("t1.nhitsxy","t1.nhitsxy",100,0.,100.); 
        TH1F* h47=new TH1F("t0.nvecxy","t0.nvecxy",100,0.,100.); 
        TH1F* h48=new TH1F("t1.nvecxy","t1.nvecxy",100,0.,100.); 
        TH1F* h49=new TH1F("t0.nvec","t0.nvec",100,0.,100.); 
        TH1F* h50=new TH1F("t1.nvec","t1.nvec",100,0.,100.); 
        TH1F* h51=new TH1F("clgamma0theta","clgamma0.theta",1000,0.,185.); 
        TH1F* h52=new TH1F("clgamma1theta","clgamma1.theta",1000,0.,185.); 
        TH1F* h53=new TH1F("cos(cl.theta)","cos(cl.theta)",100,-1.,1.); 
        TH1F* h54=new TH1F("elchi2toPIchi2","elchi2/Pionchi2",1000,0.,5.); 
        TH1F* h55=new TH1F("muchi2toPIchi2","muchi2/Pionchi2",1000,0.,5.); 
        TH1F* h56=new TH1F("Kchi2toPIchi2","Kaonchi2/Pionchi2",1000,0.,5.); 
        TH1F* h57=new TH1F("Pchi2toPIchi2","Protonchi2/Pionchi2",1000,0.,5.); 
        TH1F* h58=new TH1F("munhits","mu.nhits",30,0.,30.);        
        TH1F* h59=new TH1F("ratioP11P22","P11/P22",100,0.,10.);        
        TH1F* hzero=new TH1F("hzero","Jpsi->pi^{+}pi^{-}pi^{0}",1600,0,1600); 
        TH1F* hzero1=new TH1F("hzero1","vrt.theta2t",1000,0.,185.); 
        TH1F* hzero2=new TH1F("hzero2","cos(theta2t)",100,-1.,1.); 
        TH1F* hzero3=new TH1F("hzero3","Kaonchi2/Pionchi2",1000,0.,5.);
        TH1F* hzero4=new TH1F("hzero4","Energy",1000,0.,5000.); 
        TH1F* hzero5=new TH1F("hzero5","emc.ncls",12,0.,12.); 
        TH1F* hzero6=new TH1F("hzero6","pionchi2",4000,0.,4000.); 
           
        
        //TH2F* h21=new TH2F("h21","tof",1600,0,1600,60,0,60); 

        //TProfile* pr4=new TProfile("E/p&Momentum","E/p&Momentum",100,0,3000,0,30); 
        //TProfile* pr5=new TProfile("N_{ph.e.}&#beta#gamma","N_{ph.e.}&#beta#gamma",100,0,15,0,30); 

        vector<float> chi2min; 
        float mom=0;
        
        float eff_teoretic,mu_eff;   
	float eff[15];
	float err_eff[15];
	double num_npenotzero[15];
	double num_npetotal[15];
	float p_all[15];
	float err_p_all[15];
	float num_npezero[15];

	float eff1[15];
	float err_eff1[15];
	double num_npenotzero1[15];
	double num_npetotal1[15];
	float num_npezero1[15];

	float eff2[15];
	float err_eff2[15];
	double num_npenotzero2[15];
	double num_npetotal2[15];
	float num_npezero2[15];

	
	float pi=TMath::Pi();

	for(int ii2=0; ii2<15; ii2++)
	{
	eff[ii2]=0;
	err_eff[ii2]=0;
	num_npenotzero[ii2]=0;
	num_npetotal[ii2]=0;
	p_all[ii2]=0;
	err_p_all[ii2]=0;
	num_npezero[ii2]=0;

	eff1[ii2]=0;
	err_eff1[ii2]=0;
	num_npenotzero1[ii2]=0;
	num_npetotal1[ii2]=0;
	num_npezero1[ii2]=0;

	eff2[ii2]=0;
	err_eff2[ii2]=0;
	num_npenotzero2[ii2]=0;
	num_npetotal2[ii2]=0;
	num_npezero2[ii2]=0;
	}

       int npipi=0;                     
       Double_t m1=139.57, m2=139.57;
       float npetrh=0.7;
       float thicknpetrh=0.7;
       int ncnt1l1=0, ncnt1l2=80;
       int ncnt2l1=80, ncnt2l2=160;
       
       bool verbose=0;
       bool verbose1=0;
       
       int counter=0;
       
       int Natc=0;
       int Nselect=0;
       int Npred=0;
               
        for(int k=0; k<nentr; k++)           //цикл по всем событиям
	{  
   	  tt->GetEntry(k);
   	  
   	  if( (k %100000)==0 )cout<<k<<endl;

           	       Double_t P1, P2;
           	       Double_t P11=0, P22=0;
               	       P1=t0.p;
               	       P2=t1.p;
               	       
               	       //P11=t0.p;
               	       //P22=t1.p;
               	       
                     //cout<<jpsi.P1[0]<<"\t"<<jpsi.P1[1]<<"\t"<<jpsi.P1[2]<<"\t"<<jpsi.P1[3]<<"\t"<<jpsi.P1[4]<<endl;
   	             //cout<<jpsi.P2[0]<<"\t"<<jpsi.P2[1]<<"\t"<<jpsi.P2[2]<<"\t"<<jpsi.P2[3]<<"\t"<<jpsi.P2[4]<<endl;   	             
   	                                                
   	               float en0[3];
                       en0[0]=t0c0.e;
                       en0[1]=t0c1.e;
                       en0[2]=t0c2.e;
                       float en1[3];
                       en1[0]=t1c0.e;
                       en1[1]=t1c1.e;
                       en1[2]=t1c2.e;
                       float engamma[4];
                       engamma[0]=clgamma0.e;
                       engamma[1]=clgamma1.e;
                       engamma[2]=clgamma2.e;
                       engamma[3]=clgamma3.e;
                       float e0=0,e1=0,egamma=0;

    	               for(int i=0; i<t0.emc_ncls; i++){
    	               e0+=en0[i];
   	               }
 		       for(int i=0; i<t1.emc_ncls; i++){
    	               e1+=en1[i];
   	               }
   	               
   	               for(int i=0; i<(emc.ncls-t0.emc_ncls-t1.emc_ncls); i++)
   	               {
   	               egamma+=engamma[i]; 
   	               }
   	               
   	               float Mpi0; 
   	               
   	               float S;
   	               //S=(3/2)*(pow(t0.pt,2)+pow(t1.pt,2))/(pow(t0.p,2)+pow(t1.p,2));
   	                	           
   	                	           
   	          for(int i=0; i<5; i++)
   	  	  {
   	    	    if(jpsi.chi2[i]>0)
   	    		{
              		//cout<<"jpsi.chi2[i]="<<jpsi.chi2[i]<<endl;
              		chi2min.push_back(jpsi.chi2[i]);	    
   	    		}   	  
   	  	  } 
                   
                 int pipi=0; 
                 float kaonchi2=0; 
                 float protonchi2=0; 
                 float electronchi2=0; 
                 float muonchi2=0; 
                 float pionchi2=0; 

                 if(jpsi.numHyp>=1)
                {   
                    for(int i=0; i<5; i++)
   	  	    {		
   	  		//if( *min_element(chi2min.begin(),chi2min.end())==jpsi.chi2[i] && jpsi.M[i]>139. && jpsi.M[i]<140. && jpsi.chi2[i]>0. && jpsi.chi2[i]<4000. )
                        if( jpsi.M[i]>0.5 && jpsi.M[i]<0.6 && jpsi.chi2[i]>0. ){ electronchi2=jpsi.chi2[i]; } 
                        if( jpsi.M[i]>105. && jpsi.M[i]<106. && jpsi.chi2[i]>0. ){ muonchi2=jpsi.chi2[i]; } 
                        if( jpsi.M[i]>492. && jpsi.M[i]<494. && jpsi.chi2[i]>0. ){ kaonchi2=jpsi.chi2[i]; } 
                        if( jpsi.M[i]>938. && jpsi.M[i]<940. && jpsi.chi2[i]>0. ){ protonchi2=jpsi.chi2[i]; } 

   	  		if( jpsi.M[i]>139. && jpsi.M[i]<140. && jpsi.chi2[i]>0. && jpsi.chi2[i]<4000. )
   	  		{
   	  		   npipi++;
   	  		   //cout<<int(k)<<"\t"<<ev.evdaq<<"\t"<<"i="<<i<<"\t"<<"chi2="<<jpsi.chi2[i]<<"\t"<<"M="<<jpsi.M[i]<<"\t"<<"npipi="<<npipi<<endl;
   	  		   /*
   	  		   cout<<jpsi.P1[0]<<"\t"<<jpsi.P1[1]<<"\t"<<jpsi.P1[2]<<"\t"<<jpsi.P1[3]<<"\t"<<jpsi.P1[4]<<endl;  
   	                   cout<<jpsi.P2[0]<<"\t"<<jpsi.P2[1]<<"\t"<<jpsi.P2[2]<<"\t"<<jpsi.P2[3]<<"\t"<<jpsi.P2[4]<<endl;   	  		    	             
   	  		   cout<<jpsi.M[0]<<"\t"<<jpsi.M[1]<<"\t"<<jpsi.M[2]<<"\t"<<jpsi.M[3]<<"\t"<<jpsi.M[4]<<endl; 
   	  		   cout<<jpsi.chi2[0]<<"\t"<<jpsi.chi2[1]<<"\t"<<jpsi.chi2[2]<<"\t"<<jpsi.chi2[3]<<"\t"<<jpsi.chi2[4]<<endl; 
   	  		   cout<<i<<"\t"<<jpsi.M[i]<<"\t"<<jpsi.P1[i]<<"\t"<<jpsi.P2[i]<<"\t"<<jpsi.chi2[i]<<endl;
   	  		   */
   	  		   P11=jpsi.P1[i];
               	           P22=jpsi.P2[i];
               	           
               	           pionchi2=jpsi.chi2[i];
   	  		   pipi=1;
               	           
               	           //P11=t0.p;
               	           //P22=t1.p;                     
 
         	           //cout<<e0/P11<<"\t"<<e1/P22<<endl;   	  		    	  		   
   	  		}
   	            }
   	         }           	           
   	             //S=(3/2)*(pow(t0.pt,2)+pow(t1.pt,2))/(pow(t0.p,2)+pow(t1.p,2)); 	
   	                	   	                	               
   	               //22493   	               
   	               //22973 22969 22968 22949 22948 22944 22945 22946 22937 22936 22934 22930 22923 22912 22893 22890 22889 22887 22885 22883 22882 22881 22880 22879 22878
   	             
//   	       if( cos(pi*(vrt.theta2t)/180)>-0.998 && emc.ncls>2 && emc.ncls<=6 && t0.emc_ncls>=1 && t1.emc_ncls>=1 && e0/P11<0.85 && e1/P22<0.85 && (t0.q+t1.q)==0 && t0tof.nhits<=11 && t1tof.nhits<=11 && (emc.ncls-t0.emc_ncls-t1.emc_ncls)>=1 && mu.nhits<3  ) 
   	       if( cos(pi*(vrt.theta2t)/180)>-0.998 && emc.ncls>2 && emc.ncls<=5 && t0.emc_ncls>=1 && t1.emc_ncls>=1 && e0/t0.p<0.70 && e1/t1.p<0.70 && (t0.q+t1.q)==0 && (emc.ncls-t0.emc_ncls-t1.emc_ncls)>=1 && mu.nhits<3  && en0[0]>20 && en1[0]>20 && engamma[0]>20 ) 
   	       //if( t0.emc_ncls>=1 && t1.emc_ncls>=1 && (t0.q+t1.q)==0 && mu.nhits<3 ) 
   	       {
   	       
   	             Npred++;
                     //if(pipi==1){ cout<<"!!!!!!!--->>>>"<<kaonchi2/pionchi2<<endl; }   	             
   	             //if( electronchi2>0. && (electronchi2/pionchi2<0.80 ) ){ pipi=0; }   	
   	             //if( muonchi2>0. && (muonchi2/pionchi2<0.80 ) ){ pipi=0; }   	
   	      //!       if( kaonchi2>0. && (kaonchi2/pionchi2<1.020 ) ){ pipi=0; }   	
   	             //if( kaonchi2>0. && kaonchi2/pionchi2<0.60 ){ pipi=0; }   	
   	      //!       if( protonchi2>0. && (protonchi2/pionchi2<0.60 ) ){ pipi=0; }  
   	             
   	             h54->Fill(electronchi2/pionchi2);
   	             h55->Fill(muonchi2/pionchi2);
   	             h56->Fill(kaonchi2/pionchi2);
   	             h57->Fill(protonchi2/pionchi2);
   	            	       
   	          //if( pipi==1 && P11/t0.p<1.10 && P22/t1.p<1.10 && P11>0 && P22>0 && P11/P22>0.1 ){
   	          if( pipi==1 && P11>0 && P22>0 && P11/P22>0.1 ){
   	             
   	             Nselect++;
   	       
   	               //float E=1560.869;
           	       // Double_t E=1483.282;
           	       //float E=3096.916/2;            	       
           	       //cout<<ev.ebeam/1000<<endl;           	                  	    
           	       float E;
           	       E=ev.ebeam;
           	       if(sim==1)E=3096.916/2; 
               	       //E=3096.916/2; 
               	       Double_t E1=abs(sqrt(P11*P11+m1*m1));
		       Double_t E2=abs(sqrt(P22*P22+m2*m2));
	               Double_t pprod=abs(P11)*abs(P22)*(t0.vx*t1.vx+t0.vy*t1.vy+t0.vz*t1.vz);
              	       //Double_t pprod=P11*P22*(t0.vx*t1.vx+t0.vy*t1.vy+t0.vz*t1.vz);
              	       
               	       //Double_t rmass2=4*E*(E-E1-E2)+m1*m1+m2*m2+2*(E1*E2-pprod);
               	       Double_t rmass2=4*E*E-4*E*(E1+E2)+m1*m1+m2*m2+2*(E1*E2-pprod);
               	       Double_t InvMass=sqrt(m1*m1+m2*m2+2*(E1*E2-pprod));                                                                
               	       //Double_t InvMass=sqrt(pow((e0+e1),2)-pow(P11,2)-pow(P22,2)-2*P11*P22*(t0.vx*t1.vx+t0.vy*t1.vy+t0.vz*t1.vz));                                                                
                       //Double_t InvMass=sqrt(2*m1*m2+2*e0*e1-2*sqrt((e0*e0-m1*m1)*(e1*e1-m2*m2))*(t0.vx*t1.vx+t0.vy*t1.vy+t0.vz*t1.vz));                                                                                         	           
   	                     	   	                	              	   	           
   	            //h35->Fill(S); 	                  
   	              
   	            Mpi0=sqrt(2*clgamma0.e*clgamma1.e*(1-(clgamma0.vx*clgamma1.vx+clgamma0.vy*clgamma1.vy+clgamma0.vz*clgamma1.vz)));
   	            //Mpi0=sqrt(2*clgamma0.e*clgamma1.e*(1-cos(pi*(clgamma0.theta-clgamma1.theta)/180)));
   	            h39->Fill(Mpi0); 
   	            
   	            h59->Fill(P11/P22);
   	       	                                                    
                     
                    int kk=0;
                    int kk1=0;
                    int kk2=0;
                    int ii1=0;
                    int ii2=0;
                    
                  if(verbose) cout<<ev.run<<"\t"<<ev.evdaq<<"\t"<<"\t"<<t0.p<<"\t"<<t1.p<<"\t"<<P11<<"\t"<<P22<<"\t"<<(P11/t0.p)<<"\t"<<(P22/t1.p)<<"\t"<<(P11/t0.p)/(P22/t1.p)<<"\t"<<(t0.vx*t1.vx+t0.vy*t1.vy+t0.vz*t1.vz)<<"\t"<<clgamma0.vx*clgamma1.vx+clgamma0.vy*clgamma1.vy+clgamma0.vz*clgamma1.vz<<"\t"<<t0c0.e<<"\t"<<t0c1.e<<"\t"<<t1c0.e<<"\t"<<t1c1.e<<"\t"<<clgamma0.e<<"\t"<<clgamma1.e<<"\t"<<emc.ncls<<"\t"<<t0.emc_ncls<<"\t"<<t1.emc_ncls<<"\t"<<emc.ncls-t0.emc_ncls-t1.emc_ncls<<"\t"<<t0tof.nhits<<"\t"<<t1tof.nhits<<endl;
	                
   	            //P11=t0.p;
   	            //P22=t1.p;
                    
                    //if(pipi==1){
                    //cout<<t0atccr0.aerogel_region5<<"\t"<<t0atccr0.npe<<endl;
                    float n00=0, n01=0, n02=0, n03=0, n04=0; 
                    float n10=0, n11=0, n12=0, n13=0, n14=0; 
                   

   		    //if( rmass2>0 && sqrt(rmass2)>100 && sqrt(rmass2)<200 && InvMass<1000 ){                                                                                      
                    //if( rmass2>0 && sqrt(rmass2)>120 && sqrt(rmass2)<160 && InvMass<1000  ) {                   
                    //if( rmass2>0 && sqrt(rmass2)>120 && sqrt(rmass2)<160 ){                                                                                      
                    //if( rmass2>0 && sqrt(rmass2)>100 && sqrt(rmass2)<210 && Mpi0>110 && Mpi0<145 ){                                                                                      
                    //if( Mpi0>115 && Mpi0<145 && pipi==1 && InvMass<1500 ){                                                                                      
                    //if( pipi==1 && rmass2>0 && sqrt(rmass2)>120 && sqrt(rmass2)<160 ){                                                                                      
                    if( pipi==1 ){                                                                                      
                    //if( Mpi0>110 && Mpi0<145 && sqrt(4*E*(E-egamma))<1000 ){                                                                   
                   //if( Mpi0>110 && Mpi0<150 ){                                                                   
                    //if( InvMass>600 && InvMass<850 ){         
                                                                              
   	            if(verbose) cout<<"InvMass="<<InvMass<<"\t"<<"rmass="<<sqrt(rmass2)<<"\t"<<"Mpi0="<<Mpi0<<"\t"<<"e0="<<e0<<"\t"<<"e1="<<e1<<"\t"<<"egamma="<<egamma<<endl; 	                   
                    //=======================  	                                                                                                                                                                             
    if( t0atccr0.aerogel_region0==1 && t0atccr0.wlshit==0 && ( (t0atccr0.i>=ncnt1l1 && t0atccr0.i<ncnt1l2) ) ) { pr1->Fill(P11,t0atccr0.npe); kk=1; counter++;
                               if(verbose1) cout<<"ev.run="<<ev.run<<"\t"<<"cnt.npe="<<t0atccr0.npe<<"\t"<<"t0atccr0.i="<<t0atccr0.i<<endl;
                                             for(int ii=0; ii<=14; ii++)
                        	{
                  			p_all[ii]=100*ii+50;
                  			err_p_all[ii]=50;
                   			if(P11>100*ii&&P11<100+100*ii){
                       				if(t0atccr0.npe<=npetrh){
                                      		num_npezero1[ii]=++num_npezero1[ii];
                                      		//cout<<"cnt.npe="<<cnt.npe<<"\t"<<"t.p="<<t.p<<"\t"<<ii<<"\t"<<"num_npezero[ii]="<<num_npezero[ii]<<endl;
                             			} 
                  				if(t0atccr0.npe>npetrh){num_npenotzero1[ii]=++num_npenotzero1[ii];}
                  				num_npetotal1[ii]=++num_npetotal1[ii];
                  				eff1[ii]=num_npenotzero1[ii]/num_npetotal1[ii];
                  				err_eff1[ii]=sqrt(num_npenotzero1[ii])/num_npetotal1[ii];
                   			}
                 		}       
                    }
                    if( t1atccr0.aerogel_region0==1 && t1atccr0.wlshit==0 && ( (t1atccr0.i>=ncnt1l1 && t1atccr0.i<ncnt1l2) ) ) {pr1->Fill(P22,t1atccr0.npe); kk=1; counter++;
                                 if(verbose1) cout<<"ev.run="<<ev.run<<"\t"<<"cnt.npe="<<t1atccr0.npe<<"\t"<<"t1atccr0.i="<<t1atccr0.i<<endl;
                                                 for(int ii=0; ii<=14; ii++)
                        	{
                  			p_all[ii]=100*ii+50;
                  			err_p_all[ii]=50;
                   			if(P22>100*ii&&P22<100+100*ii){
                       				if(t1atccr0.npe<=npetrh){
                                      		num_npezero1[ii]=++num_npezero1[ii];
                                      		//cout<<"cnt.npe="<<cnt.npe<<"\t"<<"t.p="<<t.p<<"\t"<<ii<<"\t"<<"num_npezero[ii]="<<num_npezero[ii]<<endl;
                             			} 
                  				if(t1atccr0.npe>npetrh){num_npenotzero1[ii]=++num_npenotzero1[ii];}
                  				num_npetotal1[ii]=++num_npetotal1[ii];
                  				eff1[ii]=num_npenotzero1[ii]/num_npetotal1[ii];
                  				err_eff1[ii]=sqrt(num_npenotzero1[ii])/num_npetotal1[ii];
                   			}
                 		}
                    }
                    if( t0atccr1.aerogel_region0==1 && t0atccr1.wlshit==0 && ( (t0atccr1.i>=ncnt1l1 && t0atccr1.i<ncnt1l2) ) ) {pr1->Fill(P11,t0atccr1.npe); kk=1; counter++;
                           if(verbose1) cout<<"ev.run="<<ev.run<<"\t"<<"cnt.npe="<<t0atccr1.npe<<"\t"<<"t0atccr1.i="<<t0atccr1.i<<endl;
                                                     for(int ii=0; ii<=14; ii++)
                        	{
                  			p_all[ii]=100*ii+50;
                  			err_p_all[ii]=50;
                   			if(P11>100*ii&&P11<100+100*ii){
                       				if(t0atccr1.npe<=npetrh){
                                      		num_npezero1[ii]=++num_npezero1[ii];
                                      		//cout<<"cnt.npe="<<cnt.npe<<"\t"<<"t.p="<<t.p<<"\t"<<ii<<"\t"<<"num_npezero[ii]="<<num_npezero[ii]<<endl;
                             			} 
                  				if(t0atccr1.npe>npetrh){num_npenotzero1[ii]=++num_npenotzero1[ii];}
                  				num_npetotal1[ii]=++num_npetotal1[ii];
                  				eff1[ii]=num_npenotzero1[ii]/num_npetotal1[ii];
                  				err_eff1[ii]=sqrt(num_npenotzero1[ii])/num_npetotal1[ii];
                   			}
                 		}
                    }
                    if( t1atccr1.aerogel_region0==1 && t1atccr1.wlshit==0 && ( (t1atccr1.i>=ncnt1l1 && t1atccr1.i<ncnt1l2) ) ) {pr1->Fill(P22,t1atccr1.npe); kk=1; counter++;
                           if(verbose1) cout<<"ev.run="<<ev.run<<"\t"<<"cnt.npe="<<t1atccr1.npe<<"\t"<<"t1atccr1.i="<<t1atccr1.i<<endl;
                                                 for(int ii=0; ii<=14; ii++)
                        	{
                  			p_all[ii]=100*ii+50;
                  			err_p_all[ii]=50;
                   			if(P22>100*ii&&P22<100+100*ii){
                       				if(t1atccr1.npe<=npetrh){
                                      		num_npezero1[ii]=++num_npezero1[ii];
                                      		//cout<<"cnt.npe="<<cnt.npe<<"\t"<<"t.p="<<t.p<<"\t"<<ii<<"\t"<<"num_npezero[ii]="<<num_npezero[ii]<<endl;
                             			} 
                  				if(t1atccr1.npe>npetrh){num_npenotzero1[ii]=++num_npenotzero1[ii];}
                  				num_npetotal1[ii]=++num_npetotal1[ii];
                  				eff1[ii]=num_npenotzero1[ii]/num_npetotal1[ii];
                  				err_eff1[ii]=sqrt(num_npenotzero1[ii])/num_npetotal1[ii];
                   			}
                 		}
                    }
                    if( t0atccr2.aerogel_region0==1 && t0atccr2.wlshit==0 && ( (t0atccr2.i>=ncnt1l1 && t0atccr2.i<ncnt1l2) ) ) {pr1->Fill(P11,t0atccr2.npe); kk=1; counter++;
                      if(verbose1) cout<<"ev.run="<<ev.run<<"\t"<<"cnt.npe="<<t0atccr2.npe<<"\t"<<"t0atccr2.i="<<t0atccr2.i<<endl;
                         
                                                 for(int ii=0; ii<=14; ii++)
                        	{
                  			p_all[ii]=100*ii+50;
                  			err_p_all[ii]=50;
                   			if(P11>100*ii&&P11<100+100*ii){
                       				if(t0atccr2.npe<=npetrh){
                                      		num_npezero1[ii]=++num_npezero1[ii];
                                      		//cout<<"cnt.npe="<<cnt.npe<<"\t"<<"t.p="<<t.p<<"\t"<<ii<<"\t"<<"num_npezero[ii]="<<num_npezero[ii]<<endl;
                             			} 
                  				if(t0atccr2.npe>npetrh){num_npenotzero1[ii]=++num_npenotzero1[ii];}
                  				num_npetotal1[ii]=++num_npetotal1[ii];
                  				eff1[ii]=num_npenotzero1[ii]/num_npetotal1[ii];
                  				err_eff1[ii]=sqrt(num_npenotzero1[ii])/num_npetotal1[ii];
                   			}
                 		} 
                    }
                    if( t1atccr2.aerogel_region0==1 && t1atccr2.wlshit==0 && ( (t1atccr2.i>=ncnt1l1 && t1atccr2.i<ncnt1l2) ) ) {pr1->Fill(P22,t1atccr2.npe);kk=1;  counter++;
                      if(verbose1) cout<<"ev.run="<<ev.run<<"\t"<<"cnt.npe="<<t1atccr2.npe<<"\t"<<"t1atccr2.i="<<t1atccr2.i<<endl;                         
                                                 for(int ii=0; ii<=14; ii++)
                        	{
                  			p_all[ii]=100*ii+50;
                  			err_p_all[ii]=50;
                   			if(P22>100*ii&&P22<100+100*ii){
                       				if(t1atccr2.npe<=npetrh){
                                      		num_npezero1[ii]=++num_npezero1[ii];
                                      		//cout<<"cnt.npe="<<cnt.npe<<"\t"<<"t.p="<<t.p<<"\t"<<ii<<"\t"<<"num_npezero[ii]="<<num_npezero[ii]<<endl;
                             			} 
                  				if(t1atccr2.npe>npetrh){num_npenotzero1[ii]=++num_npenotzero1[ii];}
                  				num_npetotal1[ii]=++num_npetotal1[ii];
                  				eff1[ii]=num_npenotzero1[ii]/num_npetotal1[ii];
                  				err_eff1[ii]=sqrt(num_npenotzero1[ii])/num_npetotal1[ii];
                   			}
                 		}
                    } 
                    if( t0atccr3.aerogel_region0==1 && t0atccr3.wlshit==0 && ( (t0atccr3.i>=ncnt1l1 && t0atccr3.i<ncnt1l2) ) ) {pr1->Fill(P11,t0atccr3.npe); kk=1; counter++;
                      if(verbose1) cout<<"ev.run="<<ev.run<<"\t"<<"cnt.npe="<<t0atccr3.npe<<"\t"<<"t0atccr3.i="<<t0atccr3.i<<endl;                         
                                                 for(int ii=0; ii<=14; ii++)
                        	{
                  			p_all[ii]=100*ii+50;
                  			err_p_all[ii]=50;
                   			if(P11>100*ii&&P11<100+100*ii){
                       				if(t0atccr3.npe<=npetrh){
                                      		num_npezero1[ii]=++num_npezero1[ii];
                                      		//cout<<"cnt.npe="<<cnt.npe<<"\t"<<"t.p="<<t.p<<"\t"<<ii<<"\t"<<"num_npezero[ii]="<<num_npezero[ii]<<endl;
                             			} 
                  				if(t0atccr3.npe>npetrh){num_npenotzero1[ii]=++num_npenotzero1[ii];}
                  				num_npetotal1[ii]=++num_npetotal1[ii];
                  				eff1[ii]=num_npenotzero1[ii]/num_npetotal1[ii];
                  				err_eff1[ii]=sqrt(num_npenotzero1[ii])/num_npetotal1[ii];
                   			}
                 		} 
                    }
                    if( t1atccr3.aerogel_region0==1 && t1atccr3.wlshit==0 && ( (t1atccr3.i>=ncnt1l1 && t1atccr3.i<ncnt1l2) ) ) {pr1->Fill(P22,t1atccr3.npe); kk=1;  counter++;
                      if(verbose1) cout<<"ev.run="<<ev.run<<"\t"<<"cnt.npe="<<t1atccr3.npe<<"\t"<<"t1atccr3.i="<<t1atccr3.i<<endl;                         
                                                 for(int ii=0; ii<=14; ii++)
                        	{
                  			p_all[ii]=100*ii+50;
                  			err_p_all[ii]=50;
                   			if(P22>100*ii&&P22<100+100*ii){
                       				if(t1atccr3.npe<=npetrh){
                                      		num_npezero1[ii]=++num_npezero1[ii];
                                      		//cout<<"cnt.npe="<<cnt.npe<<"\t"<<"t.p="<<t.p<<"\t"<<ii<<"\t"<<"num_npezero[ii]="<<num_npezero[ii]<<endl;
                             			} 
                  				if(t1atccr3.npe>npetrh){num_npenotzero1[ii]=++num_npenotzero1[ii];}
                  				num_npetotal1[ii]=++num_npetotal1[ii];
                  				eff1[ii]=num_npenotzero1[ii]/num_npetotal1[ii];
                  				err_eff1[ii]=sqrt(num_npenotzero1[ii])/num_npetotal1[ii];
                   			}
                 		}
                    }
                   
                    //============================                    
                   
                    if( t0atccr0.aerogel_region0==1 && t0atccr0.wlshit==0 && ( (t0atccr0.i>=ncnt2l1 && t0atccr0.i<ncnt2l2) ) ) { pr2->Fill(P11,t0atccr0.npe); kk=1; counter++;
                               if(verbose1) cout<<"ev.run="<<ev.run<<"\t"<<"cnt.npe="<<t0atccr0.npe<<"\t"<<"t0atccr0.i="<<t0atccr0.i<<endl;
                                             for(int ii=0; ii<=14; ii++)
                                  {   
                  			p_all[ii]=100*ii+50;
                  			err_p_all[ii]=50;
                   			if(P11>100*ii&&P11<100+100*ii){
                       				if(t0atccr0.npe<=npetrh){
                                      		num_npezero2[ii]=++num_npezero2[ii];
                                      		//cout<<"cnt.npe="<<cnt.npe<<"\t"<<"t.p="<<t.p<<"\t"<<ii<<"\t"<<"num_npezero[ii]="<<num_npezero[ii]<<endl;
                             			} 
                  				if(t0atccr0.npe>npetrh){num_npenotzero2[ii]=++num_npenotzero2[ii];}
                  				num_npetotal2[ii]=++num_npetotal2[ii];
                  				eff2[ii]=num_npenotzero2[ii]/num_npetotal2[ii];
                  				err_eff2[ii]=sqrt(num_npenotzero2[ii])/num_npetotal2[ii];
                   			}
                 		}       
                    }
                    if( t1atccr0.aerogel_region0==1 && t1atccr0.wlshit==0 && ( (t1atccr0.i>=ncnt2l1 && t1atccr0.i<ncnt2l2) ) ) {pr2->Fill(P22,t1atccr0.npe); kk=1; counter++;
                                 if(verbose1) cout<<"ev.run="<<ev.run<<"\t"<<"cnt.npe="<<t1atccr0.npe<<"\t"<<"t1atccr0.i="<<t1atccr0.i<<endl;
                                                 for(int ii=0; ii<=14; ii++)
                        	{
                  			p_all[ii]=100*ii+50;
                  			err_p_all[ii]=50;
                   			if(P22>100*ii&&P22<100+100*ii){
                       				if(t1atccr0.npe<=npetrh){
                                      		num_npezero2[ii]=++num_npezero2[ii];
                                      		//cout<<"cnt.npe="<<cnt.npe<<"\t"<<"t.p="<<t.p<<"\t"<<ii<<"\t"<<"num_npezero[ii]="<<num_npezero[ii]<<endl;
                             			} 
                  				if(t1atccr0.npe>npetrh){num_npenotzero2[ii]=++num_npenotzero2[ii];}
                  				num_npetotal2[ii]=++num_npetotal2[ii];
                  				eff2[ii]=num_npenotzero2[ii]/num_npetotal2[ii];
                  				err_eff2[ii]=sqrt(num_npenotzero2[ii])/num_npetotal2[ii];
                   			}
                 		}
                    }
                    if( t0atccr1.aerogel_region0==1 && t0atccr1.wlshit==0 && ( (t0atccr1.i>=ncnt2l1 && t0atccr1.i<ncnt2l2) ) ) {pr2->Fill(P11,t0atccr1.npe); kk=1; counter++;
                           if(verbose1) cout<<"ev.run="<<ev.run<<"\t"<<"cnt.npe="<<t0atccr1.npe<<"\t"<<"t0atccr1.i="<<t0atccr1.i<<endl;
                                                     for(int ii=0; ii<=14; ii++)
                        	{
                   			p_all[ii]=100*ii+50;
                  			err_p_all[ii]=50;
                   			if(P11>100*ii&&P11<100+100*ii){
                      				if(t0atccr1.npe<=npetrh){
                                      		num_npezero2[ii]=++num_npezero2[ii];
                                      		//cout<<"cnt.npe="<<cnt.npe<<"\t"<<"t.p="<<t.p<<"\t"<<ii<<"\t"<<"num_npezero[ii]="<<num_npezero[ii]<<endl;
                             			} 
                  				if(t0atccr1.npe>npetrh){num_npenotzero2[ii]=++num_npenotzero2[ii];}
                  				num_npetotal2[ii]=++num_npetotal2[ii];
                  				eff2[ii]=num_npenotzero2[ii]/num_npetotal2[ii];
                  				err_eff2[ii]=sqrt(num_npenotzero2[ii])/num_npetotal2[ii];
                   			}
                 		}
                    }
                    if( t1atccr1.aerogel_region0==1 && t1atccr1.wlshit==0 && ( (t1atccr1.i>=ncnt2l1 && t1atccr1.i<ncnt2l2) ) ) {pr2->Fill(P22,t1atccr1.npe); kk=1; counter++;
                           if(verbose1) cout<<"ev.run="<<ev.run<<"\t"<<"cnt.npe="<<t1atccr1.npe<<"\t"<<"t1atccr1.i="<<t1atccr1.i<<endl;
                                                 for(int ii=0; ii<=14; ii++)
                        	{
                  			p_all[ii]=100*ii+50;
                  			err_p_all[ii]=50;
                   			if(P22>100*ii&&P22<100+100*ii){
                       				if(t1atccr1.npe<=npetrh){
                                      		num_npezero2[ii]=++num_npezero2[ii];
                                      		//cout<<"cnt.npe="<<cnt.npe<<"\t"<<"t.p="<<t.p<<"\t"<<ii<<"\t"<<"num_npezero[ii]="<<num_npezero[ii]<<endl;
                             			} 
                  				if(t1atccr1.npe>npetrh){num_npenotzero2[ii]=++num_npenotzero2[ii];}
                  				num_npetotal2[ii]=++num_npetotal2[ii];
                  				eff2[ii]=num_npenotzero2[ii]/num_npetotal2[ii];
                  				err_eff2[ii]=sqrt(num_npenotzero2[ii])/num_npetotal2[ii];
                   			}
                 		}
                    }
                    if( t0atccr2.aerogel_region0==1 && t0atccr2.wlshit==0 && ( (t0atccr2.i>=ncnt2l1 && t0atccr2.i<ncnt2l2) ) ) {pr2->Fill(P11,t0atccr2.npe); kk=1; counter++;
                      if(verbose1) cout<<"ev.run="<<ev.run<<"\t"<<"cnt.npe="<<t0atccr2.npe<<"\t"<<"t0atccr2.i="<<t0atccr2.i<<endl;
                         
                                                 for(int ii=0; ii<=14; ii++)
                        	{
                  			p_all[ii]=100*ii+50;
                  			err_p_all[ii]=50;
                   			if(P11>100*ii&&P11<100+100*ii){
                       				if(t0atccr2.npe<=npetrh){
                                      		num_npezero2[ii]=++num_npezero2[ii];
                                      		//cout<<"cnt.npe="<<cnt.npe<<"\t"<<"t.p="<<t.p<<"\t"<<ii<<"\t"<<"num_npezero[ii]="<<num_npezero[ii]<<endl;
                             			} 
                  				if(t0atccr2.npe>npetrh){num_npenotzero2[ii]=++num_npenotzero2[ii];}
                  				num_npetotal2[ii]=++num_npetotal2[ii];
                  				eff2[ii]=num_npenotzero2[ii]/num_npetotal2[ii];
                  				err_eff2[ii]=sqrt(num_npenotzero2[ii])/num_npetotal2[ii];
                   			}
                 		} 
                    }
                    if( t1atccr2.aerogel_region0==1 && t1atccr2.wlshit==0 && ( (t1atccr2.i>=ncnt2l1 && t1atccr2.i<ncnt2l2) ) ) {pr2->Fill(P22,t1atccr2.npe);kk=1;  counter++;
                      if(verbose1) cout<<"ev.run="<<ev.run<<"\t"<<"cnt.npe="<<t1atccr2.npe<<"\t"<<"t1atccr2.i="<<t1atccr2.i<<endl;                         
                                                 for(int ii=0; ii<=14; ii++)
                        	{
                  			p_all[ii]=100*ii+50;
                  			err_p_all[ii]=50;
                   			if(P22>100*ii&&P22<100+100*ii){
                       				if(t1atccr2.npe<=npetrh){
                                      		num_npezero2[ii]=++num_npezero2[ii];
                                      		//cout<<"cnt.npe="<<cnt.npe<<"\t"<<"t.p="<<t.p<<"\t"<<ii<<"\t"<<"num_npezero[ii]="<<num_npezero[ii]<<endl;
                             			} 
                  				if(t1atccr2.npe>npetrh){num_npenotzero2[ii]=++num_npenotzero2[ii];}
                  				num_npetotal2[ii]=++num_npetotal2[ii];
                  				eff2[ii]=num_npenotzero2[ii]/num_npetotal2[ii];
                  				err_eff2[ii]=sqrt(num_npenotzero2[ii])/num_npetotal2[ii];
                   			}
                 		}
                    } 
                    if( t0atccr3.aerogel_region0==1 && t0atccr3.wlshit==0 && ( (t0atccr3.i>=ncnt2l1 && t0atccr3.i<ncnt2l2) ) ) {pr2->Fill(P11,t0atccr3.npe); kk=1; counter++;
                      if(verbose1) cout<<"ev.run="<<ev.run<<"\t"<<"cnt.npe="<<t0atccr3.npe<<"\t"<<"t0atccr3.i="<<t0atccr3.i<<endl;                         
                                                 for(int ii=0; ii<=14; ii++)
                        	{
                  			p_all[ii]=100*ii+50;
                  			err_p_all[ii]=50;
                   			if(P11>100*ii&&P11<100+100*ii){
                       				if(t0atccr3.npe<=npetrh){
                                      		num_npezero2[ii]=++num_npezero2[ii];
                                      		//cout<<"cnt.npe="<<cnt.npe<<"\t"<<"t.p="<<t.p<<"\t"<<ii<<"\t"<<"num_npezero[ii]="<<num_npezero[ii]<<endl;
                             			} 
                  				if(t0atccr3.npe>npetrh){num_npenotzero2[ii]=++num_npenotzero2[ii];}
                  				num_npetotal2[ii]=++num_npetotal2[ii];
                  				eff2[ii]=num_npenotzero2[ii]/num_npetotal2[ii];
                  				err_eff2[ii]=sqrt(num_npenotzero2[ii])/num_npetotal2[ii];
                   			}
                 		} 
                    }
                    if( t1atccr3.aerogel_region0==1 && t1atccr3.wlshit==0 && ( (t1atccr3.i>=ncnt2l1 && t1atccr3.i<ncnt2l2) ) ) {pr2->Fill(P22,t1atccr3.npe); kk=1;  counter++;
                      if(verbose1) cout<<"ev.run="<<ev.run<<"\t"<<"cnt.npe="<<t1atccr3.npe<<"\t"<<"t1atccr3.i="<<t1atccr3.i<<endl;                         
                                                 for(int ii=0; ii<=14; ii++)
                        	{
                  			p_all[ii]=100*ii+50;
                  			err_p_all[ii]=50;
                   			if(P22>100*ii&&P22<100+100*ii){
                       				if(t1atccr3.npe<=npetrh){
                                      		num_npezero2[ii]=++num_npezero2[ii];
                                      		//cout<<"cnt.npe="<<cnt.npe<<"\t"<<"t.p="<<t.p<<"\t"<<ii<<"\t"<<"num_npezero[ii]="<<num_npezero[ii]<<endl;
                             			} 
                  				if(t1atccr3.npe>npetrh){num_npenotzero2[ii]=++num_npenotzero2[ii];}
                  				num_npetotal2[ii]=++num_npetotal2[ii];
                  				eff2[ii]=num_npenotzero2[ii]/num_npetotal2[ii];
                  				err_eff2[ii]=sqrt(num_npenotzero2[ii])/num_npetotal2[ii];
                   			}
                 		}
                    }
                   
                   
                   //================================
                   
                   int cnt11=0, cnt12=0;
                   int cnt21=0, cnt22=0;
                   
                     if( t0atccr0.aerogel_region0==1 && t0atccr0.wlshit==0 && ( (t0atccr0.i>=ncnt1l1 && t0atccr0.i<ncnt1l2) || (t0atccr0.i>=ncnt2l1 && t0atccr0.i<ncnt2l2) ) ) {
                                n00=t0atccr0.npe;   kk1=1;  ii1++; 
                                if(t0atccr0.i>=ncnt1l1 && t0atccr0.i<ncnt1l2)cnt11=1;
                                if(t0atccr0.i>=ncnt2l1 && t0atccr0.i<ncnt2l2)cnt12=1;                                                                                 
                    }
                    
                    if( t0atccr1.aerogel_region0==1 && t0atccr1.wlshit==0 && ( (t0atccr1.i>=ncnt1l1 && t0atccr1.i<ncnt1l2) || (t0atccr1.i>=ncnt2l1 && t0atccr1.i<ncnt2l2) ) ) {
                               n01=t0atccr1.npe;   kk1=1;   ii1++;  
                               if(t0atccr1.i>=ncnt1l1 && t0atccr1.i<ncnt1l2)cnt11=1;
                               if(t0atccr1.i>=ncnt2l1 && t0atccr1.i<ncnt2l2)cnt12=1;                                                                                 
                    }
                    if( t0atccr2.single_aerogel_region0==1 && t0atccr2.wlshit==0 && ( (t0atccr2.i>=ncnt1l1 && t0atccr2.i<ncnt1l2) || (t0atccr2.i>=ncnt2l1 && t0atccr2.i<ncnt2l2) ) ) { 
                               n02=t0atccr2.npe;   kk1=1;   ii1++;
                               if(t0atccr2.i>=ncnt1l1 && t0atccr2.i<ncnt1l2)cnt11=1;
                               if(t0atccr2.i>=ncnt2l1 && t0atccr2.i<ncnt2l2)cnt12=1;                                                                                 
                    }
                    if( t0atccr3.aerogel_region0==1 && t0atccr3.wlshit==0 && ( (t0atccr3.i>=ncnt1l1 && t0atccr3.i<ncnt1l2) || (t0atccr3.i>=ncnt2l1 && t0atccr3.i<ncnt2l2) ) ) {  
                               n03=t0atccr3.npe; kk1=1;    ii1++;
                               if(t0atccr3.i>=ncnt1l1 && t0atccr3.i<ncnt1l2)cnt11=1;
                               if(t0atccr3.i>=ncnt2l1 && t0atccr3.i<ncnt2l2)cnt12=1;                                                                                                                   
                    }
                    if( t0atccr4.aerogel_region0==1 && t0atccr4.wlshit==0 && ( (t0atccr4.i>=ncnt1l1 && t0atccr4.i<ncnt1l2) || (t0atccr4.i>=ncnt2l1 && t0atccr4.i<ncnt2l2) ) ) {  
                               n04=t0atccr4.npe; kk1=1;    ii1++;
                               if(t0atccr4.i>=ncnt1l1 && t0atccr4.i<ncnt1l2)cnt11=1;
                               if(t0atccr4.i>=ncnt2l1 && t0atccr4.i<ncnt2l2)cnt12=1;                                                                                                                   
                    }
                    
                    
                    if( t1atccr0.aerogel_region0==1 && t1atccr0.wlshit==0 && ( (t1atccr0.i>=ncnt1l1 && t1atccr0.i<ncnt1l2) || (t1atccr0.i>=ncnt2l1 && t1atccr0.i<ncnt2l2) ) ) {
                                 n10=t1atccr0.npe;  kk2=1;   ii2++; 
                                if(t1atccr0.i>=ncnt1l1 && t1atccr0.i<ncnt1l2)cnt21=1;
                                if(t1atccr0.i>=ncnt2l1 && t1atccr0.i<ncnt2l2)cnt22=1;                                                                                 
                    }                    
                    if( t1atccr1.aerogel_region0==1 && t1atccr1.wlshit==0 && ( (t1atccr1.i>=ncnt1l1 && t0atccr1.i<ncnt1l2) || (t1atccr1.i>=ncnt2l1 && t1atccr1.i<ncnt2l2) ) ) {
                                 n11=t1atccr1.npe;  kk2=1;   ii2++;  
                                if(t1atccr1.i>=ncnt1l1 && t1atccr1.i<ncnt1l2)cnt21=1;
                                if(t1atccr1.i>=ncnt2l1 && t1atccr1.i<ncnt2l2)cnt22=1;                                                                                 
                     }
                    if( t1atccr2.aerogel_region0==1 && t1atccr2.wlshit==0 && ( (t1atccr2.i>=ncnt1l1 && t1atccr2.i<ncnt1l2) || (t1atccr2.i>=ncnt2l1 && t1atccr2.i<ncnt2l2) ) ) {
                                 n12=t1atccr2.npe;  kk2=1;   ii2++; 
                                if(t1atccr2.i>=ncnt1l1 && t1atccr2.i<ncnt1l2)cnt21=1;
                                if(t1atccr2.i>=ncnt2l1 && t1atccr2.i<ncnt2l2)cnt22=1;                                                                                 
                    }                     
                    if( t1atccr3.aerogel_region0==1 && t1atccr3.wlshit==0 && ( (t1atccr3.i>=ncnt1l1 && t1atccr3.i<ncnt1l2) || (t1atccr3.i>=ncnt2l1 && t1atccr3.i<ncnt2l2) ) ) {
                                 n13=t1atccr3.npe;  kk2=1;   ii2++;                   
                                if(t1atccr3.i>=ncnt1l1 && t1atccr3.i<ncnt1l2)cnt21=1;
                                if(t1atccr3.i>=ncnt2l1 && t1atccr3.i<ncnt2l2)cnt22=1;                                                                                 
                    }
                    if( t1atccr4.aerogel_region0==1 && t1atccr4.wlshit==0 && ( (t1atccr4.i>=ncnt1l1 && t1atccr4.i<ncnt1l2) || (t1atccr4.i>=ncnt2l1 && t1atccr4.i<ncnt2l2) ) ) {
                                 n14=t1atccr4.npe;  kk2=1;   ii2++;                     
                                if(t1atccr4.i>=ncnt1l1 && t1atccr4.i<ncnt1l2)cnt21=1;
                                if(t1atccr4.i>=ncnt2l1 && t1atccr4.i<ncnt2l2)cnt22=1;                                                                                 
                    }
                    
                    
                    //if(kk1==1 && ( (n00+n01+n02+n03)>0 || ( (n00+n01+n02+n03)==0 && ii1==1 ) ) ) 
                    //if( kk1==1 && (n00+n01+n02+n03)>0 ) 
                    if( kk1==1 && ii1>1 && cnt11==1 && cnt12==1 ) 
                   {
                    prthink->Fill(P11,(n00+n01+n02+n03+n04)); 
                    
                    //if( (n00+n01+n02+n03)==0 )przero->Fill(P11,(n00+n01+n02+n03)); 
                    if( (n00+n01+n02+n03+n04)<=0.02 ){
                     hzero->Fill(P11); 
                     hzero1->Fill(vrt.theta2t);               	       
               	     hzero2->Fill(cos(pi*(vrt.theta2t)/180));
               	     hzero3->Fill(kaonchi2/pionchi2); 
               	     hzero4->Fill(emc.energy);  	           
               	     hzero5->Fill(emc.ncls);
               	     //h17->Fill(InvMass); 
                     //h18->Fill(clgamma0.vx*clgamma1.vx+clgamma0.vy*clgamma1.vy+clgamma0.vz*clgamma1.vz);
               	     //h33->Fill(P11/t0.p);
               	     //h34->Fill(P22/t1.p);
	             //S=(3/2)*(pow(t0.pt,2)+pow(t1.pt,2))/(pow(t0.p,2)+pow(t1.p,2)); 	
   	             hzero6->Fill(pionchi2); 
                	        	                          	   
         
                     }

                    Natc++;
                    if(verbose1) cout<<"ev.run="<<ev.run<<"\t"<<"cnt0.npe="<<n00+n01+n02+n03+n04<<"\t"<<endl;                        
                    for(int ii=0; ii<=14; ii++)
                    {
                  			p_all[ii]=100*ii+50;
                  			err_p_all[ii]=50;
                   			if(P11>100*ii&&P11<100+100*ii){
                       				if((n00+n01+n02+n03+n04)<=thicknpetrh){
                                      		num_npezero[ii]=++num_npezero[ii];
                                      		//cout<<"cnt.npe="<<cnt.npe<<"\t"<<"t.p="<<t.p<<"\t"<<ii<<"\t"<<"num_npezero[ii]="<<num_npezero[ii]<<endl;
                             			} 
                  				if((n00+n01+n02+n03+n04)>thicknpetrh){num_npenotzero[ii]=++num_npenotzero[ii];}
                  				num_npetotal[ii]=++num_npetotal[ii];
                  				eff[ii]=num_npenotzero[ii]/num_npetotal[ii];
                  				err_eff[ii]=sqrt(num_npenotzero[ii])/num_npetotal[ii];
                   			}
                    }
                   }
                    //if(kk2==1 && ( (n10+n11+n12+n13)>0 || ( (n10+n11+n12+n13)==0 && ii2==1 ) )  )                    
                    //if( kk2==1 && (n10+n11+n12+n13)>0 )                    
                    if( kk2==1 && ii2>1 && cnt21==1 && cnt22==1 )                    
                   {       
                    prthink->Fill(P22,(n10+n11+n12+n13+n14)); 
                    
                    //if( (n10+n11+n12+n13)==0 )przero->Fill(P22,(n10+n11+n12+n13)); 
                    if( (n10+n11+n12+n13+n14)<=0.02 ){
                     hzero->Fill(P22); 
                     hzero1->Fill(vrt.theta2t);               	       
               	     hzero2->Fill(cos(pi*(vrt.theta2t)/180));
               	     hzero3->Fill(kaonchi2/pionchi2);
               	     hzero4->Fill(emc.energy);  	                          	        	                          	   
                     hzero5->Fill(emc.ncls);  	                          	        	                          	   
                     //h17->Fill(InvMass); 
                     //h18->Fill(clgamma0.vx*clgamma1.vx+clgamma0.vy*clgamma1.vy+clgamma0.vz*clgamma1.vz);
               	     //h33->Fill(P11/t0.p);
               	     //h34->Fill(P22/t1.p);
                     //S=(3/2)*(pow(t0.pt,2)+pow(t1.pt,2))/(pow(t0.p,2)+pow(t1.p,2)); 	
   	             hzero6->Fill(pionchi2); 
           
                    }
                               
                    Natc++;
                    if(verbose1) cout<<"ev.run="<<ev.run<<"\t"<<"cnt1.npe="<<n10+n11+n12+n13<<"\t"<<endl;                        
                    for(int ii=0; ii<=14; ii++)
                    {
                  			p_all[ii]=100*ii+50;
                  			err_p_all[ii]=50;
                   			if(P22>100*ii&&P22<100+100*ii){
                       				if((n10+n11+n12+n13+n14)<=thicknpetrh){
                                      		num_npezero[ii]=++num_npezero[ii];
                                      		//cout<<"cnt.npe="<<cnt.npe<<"\t"<<"t.p="<<t.p<<"\t"<<ii<<"\t"<<"num_npezero[ii]="<<num_npezero[ii]<<endl;
                             			} 
                  				if((n10+n11+n12+n13+n14)>thicknpetrh){num_npenotzero[ii]=++num_npenotzero[ii];}
                  				num_npetotal[ii]=++num_npetotal[ii];
                  				eff[ii]=num_npenotzero[ii]/num_npetotal[ii];
                  				err_eff[ii]=sqrt(num_npenotzero[ii])/num_npetotal[ii];
                   			}
                     }
                    }
                    
                    //===============================
                       
                    }                    
                    
                    //if( ( kk1==1 && (n00+n01+n02+n03)==0 ) || ( kk2==1 && (n10+n11+n12+n13)==0 ) )
                    //if( ( kk1==1 && ii1>1 ) || ( kk2==1 && ii2>1 )
                    //{ 
                       //cout<<int(k)<<"\t"<<ev.evdaq<<"\t"<<t0.p<<"\t"<<t1.p<<endl;
                       //cout<<"\t"<<t0atccr0.i<<"\t"<<t0atccr1.i<<"\t"<<t0atccr2.i<<endl;
                       //cout<<"\t"<<t0atccr0.npe<<"\t"<<t0atccr1.npe<<"\t"<<t0atccr2.npe<<endl;
                       //cout<<"\t"<<t1atccr0.i<<"\t"<<t1atccr1.i<<"\t"<<t1atccr2.i<<endl;
                       //cout<<"\t"<<t1atccr0.npe<<"\t"<<t1atccr1.npe<<"\t"<<t1atccr2.npe<<endl; 
                       //cout<<t0atccr0.wlshit<<endl; 
                       //cout<<t0atccr1.wlshit<<endl; 
                       //cout<<t0atccr2.wlshit<<endl; 
                       //cout<<t1atccr0.wlshit<<endl; 
                       //cout<<t1atccr1.wlshit<<endl; 
                       //cout<<t1atccr2.wlshit<<endl;                        
                                     
                       //h1->Fill(e0/P11); 
                       //h1->Fill(e1/P22); 
                       h1->Fill(e0/t0.p); 
                       h1->Fill(e1/t1.p); 
	               h2->Fill(emc.energy);
		       h3->Fill(emc.elkr);
	               h4->Fill(emc.ecsi);
	     	       if(t0.q<0 && t0c0.lkr>0 )h5->Fill(e0);
	     	       if(t1.q<0 && t1c0.lkr>0 )h5->Fill(e1);
	               if(t0.q>0 && t0c0.lkr>0)h6->Fill(e0);
	               if(t1.q>0 && t1c0.lkr>0)h6->Fill(e1);
	               h7->Fill(P11);	     
               	       h7->Fill(P22);	     
               	       h8->Fill(t0.theta);	            
    	               h8->Fill(t1.theta);	            
    	               h9->Fill(t0.phi);
               	       h9->Fill(t1.phi);
               	       h10->Fill(vrt.theta2t);
               	       h11->Fill(vrt.phi2t);
               	       
               	       h12->Fill(cos(pi*(vrt.theta2t)/180));
               	       h13->Fill(cos(pi*vrt.phi2t/180));
               	       h14->Fill(cos(pi*t0.theta/180));
               	       h15->Fill(cos(pi*t0.phi/180));
              	       h14->Fill(cos(pi*t1.theta/180));
               	       h15->Fill(cos(pi*t1.phi/180));
               	       
               	       if(rmass2>0){h16->Fill(sqrt(rmass2));}
               	       h17->Fill(InvMass);
               	       h18->Fill(clgamma0.vx*clgamma1.vx+clgamma0.vy*clgamma1.vy+clgamma0.vz*clgamma1.vz);
               	       
               	       h19->Fill(t0tof.nhits);
               	       h20->Fill(t0tof.dchits);
               	       h21->Fill(t0tof.namps);
               	       h22->Fill(t0tof.ntimes);
               	       h23->Fill(t0tof.time[0]);
               	       h24->Fill(t0tof.beta[0]);
               	       h25->Fill(t0tof.length[0]);
               	       
               	       h26->Fill(t1tof.nhits);
               	       h27->Fill(t1tof.dchits);
               	       h28->Fill(t1tof.namps);
               	       h29->Fill(t1tof.ntimes);
               	       h30->Fill(t1tof.time[1]);
               	       h31->Fill(t1tof.beta[1]);
               	       h32->Fill(t1tof.length[1]);
               	       
               	       h33->Fill(P11/t0.p);
               	       h34->Fill(P22/t1.p);

   	               //Mpi0=sqrt(pow(egamma,2)-(pow(e0,2)-pow(m1,2))-(pow(e1,2)-pow(m2,2))-2*sqrt((pow(e0,2)-pow(m1,2))*(pow(e1,2)-pow(m2,2)))*(t0.vx*t1.vx+t0.vy*t1.vy+t0.vz*t1.vz));
               	       //h36->Fill(sqrt(pow(egamma,2)-(pow(e0,2)-pow(m1,2))-(pow(e1,2)-pow(m2,2))-2*sqrt((pow(e0,2)-pow(m1,2))*(pow(e1,2)-pow(m2,2)))*(t0.vx*t1.vx+t0.vy*t1.vy+t0.vz*t1.vz)));
               	       //h36->Fill(sqrt(4*E*(E-sqrt(pow(clgamma0.e,2)+pow(clgamma1.e,2)+clgamma0.e*clgamma1.e))));
               	       h36->Fill(sqrt(4*E*(E-egamma)));
               	       
               	       //cout<<sqrt(4*E*(E-sqrt(pow(clgamma0.e,2)+pow(clgamma1.e,2)+clgamma0.e*clgamma1.e)))<<endl;
               	       
             	       S=(3/2)*(pow(t0.pt,2)+pow(t1.pt,2))/(pow(t0.p,2)+pow(t1.p,2));
   	               h35->Fill(S);

                       h40->Fill(t0.chi2); 
                       h41->Fill(t1.chi2); 
                       h42->Fill(t0.vx*t1.vx+t0.vy*t1.vy+t0.vz*t1.vz);
                       h43->Fill(emc.ncls);
                       h44->Fill((P11/t0.p)/(P22/t1.p));
                       h45->Fill(t0.nhitsxy);
                       h46->Fill(t1.nhitsxy);
                       h47->Fill(t0.nvecxy);
                       h48->Fill(t1.nvecxy);
                       h49->Fill(t0.nvec);
                       h50->Fill(t1.nvec);
                       h51->Fill(clgamma0.theta);
                       h52->Fill(clgamma1.theta);
                       h53->Fill(cos(pi*(clgamma0.theta-clgamma1.theta)/180));
                       h58->Fill(mu.nhits);
   	           

                       
                       if(verbose1)cout<<t0c0.theta<<"\t"<<t1c0.theta<<"\t"<<clgamma0.theta<<"\t"<<clgamma1.theta<<"\t"<<(clgamma0.theta+clgamma1.theta)/2<<"\t"<<clgamma0.theta-clgamma1.theta<<"\t"<<clgamma0.vx*clgamma1.vx+clgamma0.vy*clgamma1.vy+clgamma0.vz*clgamma1.vz<<endl;

               	       /*
               	       		for(int ii=0; ii<=14; ii++)
                        	{
                  			p_all[ii]=100*ii+50;
                  			err_p_all[ii]=50;
                   			if(t0.p>100*ii&&t0.p<100+100*ii){
                       				if(cnt.npe<=0.1){
                                      		num_npezero[ii]=++num_npezero[ii];
                                      		//cout<<"cnt.npe="<<cnt.npe<<"\t"<<"t.p="<<t.p<<"\t"<<ii<<"\t"<<"num_npezero[ii]="<<num_npezero[ii]<<endl;
                             			} 
                  				if(cnt.npe>0.1){num_npenotzero[ii]=++num_npenotzero[ii];}
                  				num_npetotal[ii]=++num_npetotal[ii];
                  				eff[ii]=num_npenotzero[ii]/num_npetotal[ii];
                  				err_eff[ii]=sqrt(num_npenotzero[ii])/num_npetotal[ii];
                   			}
                 		}        
               	                */  
               	      //cout<<ev.run<<"\t"<<ev.evdaq<<"\t"<<sqrt(rmass2)<<"\t"<<InvMass<<"\t"<<t0.p<<"\t"<<t1.p<<"\t"<<P1<<"\t"<<P2<<"\t"<<pprod<<"\t"<<(t0.vx*t1.vx+t0.vy*t1.vy+t0.vz*t1.vz)<<"\t"<<t0c0.e<<"\t"<<t0c1.e<<"\t"<<t1c0.e<<"\t"<<t1c1.e<<"\t"<<clgamma0.e<<"\t"<<clgamma1.e<<"\t"<<emc.ncls<<"\t"<<t0.emc_ncls<<"\t"<<t1.emc_ncls<<"\t"<<emc.ncls-t0.emc_ncls-t1.emc_ncls<<endl;

                     //cout<<jpsi.P1[0]<<"\t"<<jpsi.P1[1]<<"\t"<<jpsi.P1[2]<<"\t"<<jpsi.P1[3]<<"\t"<<jpsi.P1[4]<<endl;
   	             //cout<<jpsi.P2[0]<<"\t"<<jpsi.P2[1]<<"\t"<<jpsi.P2[2]<<"\t"<<jpsi.P2[3]<<"\t"<<jpsi.P2[4]<<endl;
   	             //cout<<"======================================================================================"<<endl;
           	   
                    //}
                    //kk1=0;
                    //kk2=0;
                    //}

                
                     //}
                    //}      
                   //}  //pipi
                   chi2min.clear();
                   //memset(&data5,0,sizeof(data5));
               /*    
                   	  clgamma1.reset();                   
             
 	   memset(&t0c0,0,sizeof(data5));
                   memset(&t0c1,0,sizeof(data5));
                   memset(&t1c0,0,sizeof(data5));
                   memset(&t1c1,0,sizeof(data5));
                   memset(&clgamma0,0,sizeof(data5));
                   memset(&clgamma1,0,sizeof(data5));
		*/
                  
	  	
	      } //pipi	
	      }  //if  t emc ...
	   }
	
	cout<<"Natc="<<Natc<<"\t"<<"ATC from selected events, % ->"<<float(Natc*100/(Nselect*2))<<endl;
	cout<<"Npipi="<<npipi<<endl;
	cout<<"Npred="<<Npred<<"\t"<<"Nselect="<<Nselect<<endl;
	
	
	if(verbose1) cout<<counter<<endl;  
	        
        TF1* myfit=new TF1("myfit","[0]+[1]*(x^2-[2]^2)/(x^2)*(TMath::Erf(x-[2])+1)/2",200.,1450.);
	myfit->SetParameter(0,prthink->GetMinimum());
	myfit->SetParameter(1,prthink->GetMaximum());
	myfit->SetParameter(2,430);
	myfit->SetParNames("#mu_{0}","#mu_{max}","P_{thr}");
	if(sim!=1){
        myfit->SetLineColor(kRed);
        prthink->SetLineColor(kRed);
        }
        else{
        myfit->SetLineColor(kBlue);
        prthink->SetLineColor(kBlue);        
        }
	prthink->Fit("myfit","","",200,1450);
        cout<<"pr->GetMinimum()="<<prthink->GetMinimum()<<"\t"<<"pr->GetMaximum()"<<prthink->GetMaximum()<<endl;        
	prthink->SetXTitle("P, MeV/c");
	prthink->SetYTitle("N_{ph.e.}");
        prthink->Draw("prof");


        TF1* myfit1=new TF1("myfit1","[0]+[1]*(x^2-[2]^2)/(x^2)*(TMath::Erf(x-[2])+1)/2",200.,1450.);
	myfit1->SetParameter(0,pr1->GetMinimum());
	myfit1->SetParameter(1,pr1->GetMaximum());
	myfit1->SetParameter(2,90);
	myfit1->SetParNames("#mu_{0}","#mu_{max}","P_{thr}");
	if(sim!=1){
        myfit1->SetLineColor(kRed);
        pr1->SetLineColor(kRed);
        }
        else{
        myfit1->SetLineColor(kBlue);
        pr1->SetLineColor(kBlue);        
        }
	pr1->Fit("myfit1","","",200,1450);
        cout<<"pr1->GetMinimum()="<<pr1->GetMinimum()<<"\t"<<"pr1->GetMaximum()"<<pr1->GetMaximum()<<endl;        
	pr1->SetXTitle("P, MeV/c");
	pr1->SetYTitle("N_{ph.e.}");
        pr1->Draw("prof");


        TF1* myfit2=new TF1("myfit2","[0]+[1]*(x^2-[2]^2)/(x^2)*(TMath::Erf(x-[2])+1)/2",200.,1450.);
	myfit2->SetParameter(0,pr2->GetMinimum());
	myfit2->SetParameter(1,pr2->GetMaximum());
	myfit2->SetParameter(2,90);
	myfit2->SetParNames("#mu_{0}","#mu_{max}","P_{thr}");
        if(sim!=1){
        myfit2->SetLineColor(kRed);
        pr2->SetLineColor(kRed);
        }
        else{
        myfit2->SetLineColor(kBlue);
        pr2->SetLineColor(kBlue);        
        }
	pr2->Fit("myfit2","","",200,1450);
        cout<<"pr2->GetMinimum()="<<pr2->GetMinimum()<<"\t"<<"pr2->GetMaximum()"<<pr2->GetMaximum()<<endl;        
	pr2->SetXTitle("P, MeV/c");
	pr2->SetYTitle("N_{ph.e.}");
        pr2->Draw("prof");

        
       TGraphErrors* gr1=new TGraphErrors(14,p_all,eff,err_p_all,err_eff);
       gr1->SetMarkerStyle(20);
       gr1->SetLineWidth(2);
       if(sim!=1){
       gr1->SetMarkerColor(2);
       gr1->SetLineColor(2);
       }
       else{
       gr1->SetMarkerColor(kBlue);
       gr1->SetLineColor(kBlue);
       }
       gr1->SetTitle("Efficiency&Momentum");
       gr1->GetXaxis()->SetTitle("P, MeV/c");
       gr1->GetYaxis()->SetTitle("Efficiency");
       //gr1->Draw("ap");
       gr1->Write("Think_Efficiency&Momentum"); 


       TGraphErrors* gr2=new TGraphErrors(14,p_all,eff1,err_p_all,err_eff1);
       gr2->SetMarkerStyle(20);
       gr2->SetLineWidth(2);
       if(sim!=1){
       gr2->SetMarkerColor(2);
       gr2->SetLineColor(2);
       }
       else{
       gr2->SetMarkerColor(kBlue);
       gr2->SetLineColor(kBlue);
       }
       gr2->SetTitle("Efficiency&Momentum");
       gr2->GetXaxis()->SetTitle("P, MeV/c");
       gr2->GetYaxis()->SetTitle("Efficiency");
       gr2->Write("1layer_Efficiency&Momentum"); 

       TGraphErrors* gr3=new TGraphErrors(14,p_all,eff2,err_p_all,err_eff2);
       gr3->SetMarkerStyle(20);
       gr3->SetLineWidth(2);
       if(sim!=1){
       gr3->SetMarkerColor(2);
       gr3->SetLineColor(2);
       }
       else{
       gr3->SetMarkerColor(kBlue);
       gr3->SetLineColor(kBlue);
       }
       gr3->SetTitle("Efficiency&Momentum");
       gr3->GetXaxis()->SetTitle("P, MeV/c");
       gr3->GetYaxis()->SetTitle("Efficiency");
       gr3->Write("2layer_Efficiency&Momentum"); 

       
       TH1F *h37 = (TH1F*) h16->Clone();
       h37->SetName("h37");
       Double_t b=h37->GetEntries();
       h37->Scale(1/b);
       h37->GetXaxis()->SetTitle("m(#pi^{0}), MeV/c^{2}");
     
       TH1F *h38 = (TH1F*) h17->Clone();
       h38->SetName("h38");
       Double_t bb=h38->GetEntries();
       h38->Scale(1/bb);
       h38->GetXaxis()->SetTitle("m(#pi#pi), MeV/c^{2}");
              
        if(sim!=1){
        h37->SetLineColor(kRed);
        h38->SetLineColor(kRed);
        }
        else{
        h37->SetLineColor(kBlue);
        h38->SetLineColor(kBlue);
        }
       
 
       TH1F *h60 = (TH1F*) h56->Clone();
       h60->SetName("h60");
       Double_t bbb=h60->GetEntries();
       h60->Scale(1/bbb);
       h60->GetXaxis()->SetTitle("K#chi^{2}/#pi#chi^{2}");       
        if(sim!=1){
        h60->SetLineColor(kRed);
        }
        else{
        h60->SetLineColor(kBlue);
        }
 
       fout->Write();                                                      //ÐÉÛÅÍ ÄÁÎÎÙÅ × ÆÁÊÌ É ÚÁËÒÙ×ÁÅÍ ÅÇÏ
       fout->Close();
} 


