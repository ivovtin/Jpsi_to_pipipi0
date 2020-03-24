#include <unistd.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <list>
#include <algorithm>    // std::min_element, std::max_element

#include "ReadNat/rr_def.h"
#include "ReadNat/re_def.h"
#include "ReadNat/ss_def.h"
#include "VDDCRec/ktracks.h"
#include "VDDCRec/mtofhits.h"
#include "VDDCRec/ToFTrack.hh"
#include "KrToF/tof_system.h"
#include "KEmcRec/emc_struct.h"
#include "KrAtc/atcrec.h"
#include "KrAtc/atc_to_track.h"
#include "KrVDDCMu/dcmu.h"
#include "KrMu/mu_system.h"
#include "KrMu/mu_event.h"

#include "TTree.h"
#include "TFolder.h"
#include "TH1.h"
#include "TBranch.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TBenchmark.h"

#include "KaFramework/kframework.h"
#include "KaFramework/eventbr.h"
#include "KaFramework/vddcbr.h"
#include "KaFramework/emcbr.h"
#include "KaFramework/atcbr.h"
#include "KaFramework/tofbr.h"
#include "KaFramework/mubr.h"
#include "KrAtc/AtcHit.hh"

#include "VDDCRec/kdcvd.h"
#include "VDDCRec/kvd.h"
#include "VDDCRec/kdcpar.h"
#include "VDDCRec/ktracks.h"
#include "VDDCRec/khits.h"
#include "VDDCRec/ktrkhits.h"
#include "VDDCRec/ktrkfpar.h"
#include "VDDCRec/kglobparam.h"
#include "ReadNat/mcdc_def.h"
#include "ReadNat/dcrawhitspar.h"
//#include "KDisplay/kdisplay_event.h"
//#include "KrdEdxPId/KrdEdxPId.hh"
//#include "ReadNat/mc_def.h"
#include "KEmcRec/emc_system.h"
#include "KEmcRec/emc_struct.h"
//#include "ReadNat/mccardindex.h"
//#include "ReadNat/mccard_def.h"
//#include "KrMu/mu_system.h"
#include "KrKRec/KFitOutMultiCPi0Gamma.h"
#include "KrKRec/RecKFitMultiCPi0Gamma.hh"

#include "KcSys/fortran.h"

using namespace std;

const int Ntraks=3;
const int Natccr=6;
const int Nclcr=2;

static const char* progname;

struct ProgramParameters {
	bool read_reco;
	int min_track_number;                 //миниммальное число треков
	int max_track_number;                 //максимальное число треков
	int min_beam_track_number;            //миниммальное число пучковых треков
	int min_ip_track_number;              //миниммальное число треков из точки взаимодействия
	int max_ip_track_number;              //максимальное число треков из точки взаимодействия
        float min_momentum;                   //минимальный импульс                                           1000 МэВ/с
        float max_momentum;                   //максимальный импульс                                          5000 МэВ/с
	float min_cluster_energy; //MeV       //минимальная энергия кластера                                  450 МэВ
	float min_total_energy; //MeV         //минимальная общая энергия                                     2000 МэВ
	int min_cluster_number;               //минимальное кол-во кластеров с энергией больше минимальной    2
	int max_cluster_number;               //минимальное кол-во кластеров с энергией больше минимальной    2
	int min_lkrcl_number;                 //минимальное кол-во кластеров в LKr
	int min_csicl_number;                 //минимальное кол-во кластеров в CsI
	const char* rootfile;
	int MCCalibRunNumber;                 //19862 - номер захода - загрузка из БД калибровок при обработке файла моделирования в заданном интервале заходов
	int NEvents;                          //число обрабатываемых событий
	bool process_only;                    //process one event only
};
//=======================================================================================================================================
//Условия отбора для pi+pi-pi0:
//1)минимальное число треков 2
//2)максимальное число треков 2
//3)миниммальное число пучковых треков 2
//4)минимальное число треков из точки взаимодействия 2
//5)максимальное число треков из точки взаимодействия 2
//6)минимальный импульс 100 МэВ/с
//7)минимальная энергия кластера 50 МэВ
//8)Суммарная энергия по всем кластерам больше 1200 MeV
//9)минимальное число кластеров 3 - для 3pi, 2-для e+e-
//10)максимальное число кластеров 5 - для 3pi
//10)треки из вершины tVertex(t)==1
////11)центральные треки 	fabs(tX0(t)-eVertexX)<1 && fabs(tY0(t)-eVertexY)<1 && fabs(tZ0(t)-eVertexZ)<15 )
//=======================================================================================================================================

//если правильно восстанавливается энергия в калориметре
//отбор по трекам
//static const struct ProgramParameters def_progpar={false,2,2,2,2,2,0,1000000,50,500,3,5,0,0,"/store/users/ovtin/out.root",19862,0,false};
//static const struct ProgramParameters def_progpar={false,2,2,2,2,2,0,1000000,50,200,2,8,0,0,"/store/users/ovtin/out.root",19862,0,false};
//static const struct ProgramParameters def_progpar={false,2,2,2,0,2,0,1000000,50,100,2,8,0,0,"/store/users/ovtin/out.root",19862,0,false};
static const struct ProgramParameters def_progpar={false,2,2,2,0,2,0,1000000,15,45,2,8,0,0,"/store/users/ovtin/out.root",19862,0,false};

static struct ProgramParameters progpar(def_progpar);

enum {                                                                          //перечисление
	CutLongDcRecord=1,
	CutLongVdRecord,
        AtcEventDamageCut,
        AtcIllegalTrackCut,
//      NoAtcOnTrackCut,
        CutTotalEnergyFast,
	CutTrackNumber,
	CutBeamTrackNumber,
	CutIPTrackNumber_min,
	CutIPTrackNumber_max,
        Cut_CPhi2t,            //для отброса событий по расколлинеарности по углам
	Cut_CThe2t,
        EMCthetaCut,
//        Cut_10per_Energy,      //отбрасываю события если выделившаяся энергия в калориметре не приписанная двум рассматриваемым кластерам превышает 10 % от полной
	CutTotalEnergy,
        CutClusterEnergy,
	CutMINClusterNumber,
	CutMAXClusterNumber,
	CutLkrClusterNumber,
	CutCsiClusterNumber,
        Cut_nclntrack,
	OneTrackCut,
	MinMomentumCut,
	MaxMomentumCut,
	NotVertex,
	ChargeCut,
	VectorsZCut,
        HitsVDCut,
	MinXYVectorsCut,
	Chi2Cut,
	tX0Cut,
	tY0Cut,
	tZ0Cut,
        TofCut,
        MUCut,
};

static TTree *eventTree, *trackTree, *emcTowerTree, *lkrStripTree, *lkrStripTrackTree, *ATCTree, *jpsiTree;             //создаем деревья
static struct EventBranch bevent;
static struct VertexBranch bvertex;
//static struct TrackBranch btrack;
//static struct TrackBranch btrack[Ntraks];
//static struct ToFBranch btof[Ntraks];
static struct TrackBranch btrack[3];
static struct ToFBranch btof[3];
static struct MUBranch bmu;
static struct EMCBranch bemc;
//static struct TowerClusterBranch bcluster[Nclcr][Ntraks];  //[2]
//static struct TowerClusterBranch bclgamma[Nclcr];
static struct TowerClusterBranch bcluster[4][3];  //[2]
static struct TowerClusterBranch bclgamma[4];
static struct StripClusterBranch bstrip;
static struct StripTrackBranch bstriptrack;

//static struct ATCCounterBranch bcnt[160];
//static struct ATCCounterBranch bcnt[Natccr][Ntraks];                                                     // создаем АЧС структуру для счетчика
static struct ATCCounterBranch bcnt[6][3];                                                     // создаем АЧС структуру для счетчика
static struct ATCBranch batc;

//TH1F *hpx;
//struct piBranch {
//        Int_t numHyp;
//	Float_t  chi2,M1,P1,M2,P2;
//};
//static struct piBranch bpi;
//static const char* piBranchList="numHyp/I:chi2/F:M1:P1:M2:P2";

typedef struct {Int_t numHyp; Float_t chi2[5],M[5],P1[5],P2[5];} JPSI;
static JPSI jpsi;


int listtr[6],lclpi0[8],lclg[5];
//Rejection prior to reconstruction helps to save CPU time
extern "C" void kemc_energy_(float Ecsi[2],float *Elkr);

int pre_event_rejection()                                                                    //предварительный отброс события
{
    //maximum track number supposed to be 7
	//maximum record size for DC: Ntracks*maxNhits*HitLength + 10% (spare)
	static const int maxDClen=1940; //=1940
	//maximum record size for VD: Ntracks*maxNhits*Length + 2 (for dT) + 10% (spare)
	static const int maxVDlen=156; //=156

	//if( RawLength(SS_DC)>maxDClen ) return CutLongDcRecord; //too long DC event

	//if( RawLength(SS_VD)>maxVDlen ) return CutLongVdRecord; //too long VD event

	float Ecsi[2], Elkr;
	kemc_energy_(Ecsi,&Elkr);

    //Total energy deposition should be more than progpar.min_total_energy
//	if( Ecsi[0]+Ecsi[1]+Elkr<progpar.min_total_energy ) return CutTotalEnergyFast;

	return 0;
}

int vddc_event_rejection()
{
	if( eTracksAll<progpar.min_track_number )       return CutTrackNumber;
	if( eTracksAll>progpar.max_track_number )       return CutTrackNumber;
	if( eTracksBeam<progpar.min_beam_track_number ) return CutBeamTrackNumber;
	if( eTracksIP<progpar.min_ip_track_number )     return CutIPTrackNumber_min;
	if( eTracksIP>progpar.max_ip_track_number )     return CutIPTrackNumber_max;

	if( tP(0)<=progpar.min_momentum )       return MinMomentumCut;
	if( tP(1)<=progpar.min_momentum )       return MinMomentumCut;
	//if( tP(2)<=progpar.min_momentum )       return MinMomentumCut;
	//if( tP(3)<=progpar.min_momentum )       return MinMomentumCut;

	if( tP(0)>=progpar.max_momentum )       return MaxMomentumCut;
	if( tP(1)>=progpar.max_momentum )       return MaxMomentumCut;
	//if( tP(2)>=progpar.max_momentum )       return MaxMomentumCut;
	//if( tP(3)>=progpar.max_momentum )       return MaxMomentumCut;

//	if( tVertex(0)!=1 )      return NotVertex;   //если tVertex(t)!=1, то трек не из вершины        //not in sim  !!!
//	if( tVertex(1)!=1 )      return NotVertex;                                                      //not in sim  !!!
//	if( tVertex(2)!=1 )      return NotVertex;                                                      //not in sim  !!!
//	if( tVertex(3)!=1 )      return NotVertex;                                                      //not in sim  !!!

//	if( tVertex(4)!=1 )      return NotVertex;   //если tVertex(t)!=1, то трек не из вершины
//	if( (tVertex(0)!=1 && tVertex(1)!=1) || (tVertex(0)!=1 && tVertex(2)!=1) || (tVertex(1)!=1 && tVertex(2)!=1) )      return NotVertex;   //если tVertex(t)!=1, то трек не из вершины

	if( (tCharge(0)+tCharge(1))!=0 )         return ChargeCut;   //если сумма зарядов треков не равна 0 то отбрасываем

//        if(tVectorsZ(0)<1 || tVectorsZ(1)<1)    return VectorsZCut;
//	if(tHitsVD(0)<1)                        return HitsVDCut;

	//cout<<"tX0(0)="<<tX0(0)<<"\t"<<"tY0(0)="<<tY0(0)<<"\t"<<"tZ0(0)="<<tZ0(0)<<"\t"<<endl;
	//cout<<"tX0(1)="<<tX0(1)<<"\t"<<"tY0(1)="<<tY0(1)<<"\t"<<"tZ0(1)="<<tZ0(1)<<"\t"<<endl;
        //зажимаем точку взаимодействия
//        if( tX0(0)<-0.4 || tX0(0)>0.1 || tX0(1)<-0.4 || tX0(1)>0.1 )     return tX0Cut;
//        if( tY0(0)<0.6 || tY0(0)>1.0 || tY0(1)<0.6 || tY0(1)>1.0 )       return tY0Cut;
//        if( tZ0(0)<-5 || tZ0(0)>5 || tZ0(1)<-5 || tZ0(1)>5 )             return tZ0Cut;

	/*
	 if( tVectorsXY(0)<3 )                   return MinXYVectorsCut;
	 if( tCh2(0)/tHits(0)>20 )               return Chi2Cut;
	 */
//	if( tCh2(0)>6 || tCh2(1)>6 )              return Chi2Cut;

	//CPhi2t - угол (XY ДК) между треками (начальное направление) для двухтрековых событий
	//CThe2t - угол (в пространстве) между треками (начальное направление) для двухтрековых событий
	//cout<<"CPhi2t="<<acos(CPhi2t)*180./M_PI<<"\t"<<"CThe2t="<<acos(CThe2t)*180./M_PI<<endl;
//	if( (acos(CPhi2t)*180./M_PI)<174. )    return Cut_CPhi2t;       //обрезка по расколлинеарности
//	if( (acos(CThe2t)*180./M_PI)<176. )    return Cut_CThe2t;
//	if( cos(CThe2t)<0.7 )    return Cut_CThe2t;

	return 0;
}

int tof_event_rejection()
{
//    if( kschit_.time_ns[0]<-4 || kschit_.time_ns[0]>4)         return TofCut;
//    if( kschit_.time_ns[1]<-4 || kschit_.time_ns[1]>4)         return TofCut;

    return 0;
}

int emc_event_rejection()
{
	int nemc=0, nlkr=0, ncsi=0;
	float tot_energy=0;
	float two_cluster_energy=0;
        int nclntrack=0;
//        int add_nemc;                                                //дополнительное кол-во сработавших кластеров
//	float add_nemc_energy=0;                                       //энергия доп-х сработавших кластеров

	for(int c=0; c<semc.emc_ncls; c++) {                           //цикл по числу кластеров в калориметре
		tot_energy+=semc.emc_energy[c];                        //подсчет энергии выделившейся в калориметре
		if( semc.emc_energy[c]>progpar.min_cluster_energy ) {  //если энергия кластера больше установленной минимальной энергии
			nemc++;                                        //число кластеров с энергией больше минимальной
                        two_cluster_energy+=semc.emc_energy[c];        //сумма энергий в кластерах с энергией больше минимальной энергии кластера
			if( semc.emc_type[c]==1 )
				nlkr++;
			else
				ncsi++;
		}
		else {
	                return CutClusterEnergy;
                //      add_nemc++;
		//      add_nemc_energy+=semc.emc_energy[c];
		}

		if ( semc.emc_dc_ntrk[c]==1 )            //если кластер привязан к треку то считаем их число      //semc.dc_emc_ncls[t]
		{
		    nclntrack++;
		}

        //if(semc.emc_theta[c]<20 || semc.emc_theta[c]>31 || semc.emc_theta[c]<39 || semc.emc_theta[c]>141 || semc.emc_theta[c]<149 || semc.emc_theta[c]>160) return EMCthetaCut;
	}

        //semc.emc_theta[c];            //угол theta в градусах

//	cout<<"tot_energy="<<tot_energy<<"\t"<<"two_cluster_energy="<<two_cluster_energy<<endl;

	//отбрасываем если какое-либо условие не выполняется
            //        if(tot_energy-two_cluster_energy>0.1*tot_energy)    return  Cut_10per_Energy;
	if( tot_energy<progpar.min_total_energy )           return CutTotalEnergy;
	if( nemc<progpar.min_cluster_number )               return CutMINClusterNumber;
	if( nemc>progpar.max_cluster_number )               return CutMAXClusterNumber;
	//if( nlkr<progpar.min_lkrcl_number )                 return CutLkrClusterNumber;
	//if( ncsi<progpar.min_csicl_number )                 return CutCsiClusterNumber;

	if ( nclntrack<2 || semc.dc_emc_ncls[0]<1 || semc.dc_emc_ncls[1]<1 )  {   //общее число привязанных кластеров должно быть не меньше 2-х и к каждому треку должен быть привязан хотя бы 1 кластер
	    //cout<<"Ev:"<<eNumber<<"\t"<<"nclntrack="<<nclntrack<<"\t"<<"semc.dc_emc_ncls[0]="<<semc.dc_emc_ncls[0]<<"\t"<<"semc.dc_emc_ncls[1]="<<semc.dc_emc_ncls[1]<<endl;
	    return Cut_nclntrack;
	}
	//if(add_nemc<3&&add_nemc_energy>120)               return Cut_add_nemc;
           //if(add_nemc<2&&add_nemc_energy>160)               return Cut_add_nemc;

	return 0;
}

int atc_rejection()          //отброс событий в АЧС
{
    //ATC raw record damaged in any way (including when DeltaT is absent or out of range)     //сырые данные повреждены
    //if( atc_rec.eventdamage ) return AtcEventDamageCut;

    //track determined as illegal by atcrec                                                   //трек не определяется в atcrec
    //if( atc_track.illtrack[0] ) return AtcIllegalTrackCut;

    //no counters on tracks, very strange if occurs                                           //нет счетчика на трек
//    if( atc_track.ncnt_on_track[0]==0 ) return NoAtcOnTrackCut;

    return 0;
}

int mu_event_rejection()
{
    //if( mu_next_event()>0 ) return MUCut;           //отброс космики

    return 0;
}


int analyse_event()         //анализ события
{
        //float EMinPhot=50.;
        float EMinPhot=15.;
	//double  WTotal=3096.916;
	double  WTotal=2*beam_energy;
	if( kedrrun_cb_.Header.RunType == 64 ) { WTotal=3096.916; }
        //cout<<"RunNumber="<<kedrraw_.Header.RunNumber<<"\t"<<"WTotal="<<WTotal<<"\t"<<"Event="<<kdcenum_.EvNum<<"\t"<<"Raw event="<<kedrraw_.Header.Number<<"\t"<<"eTracksAll="<<eTracksAll<<endl;


	//if(eTracksBeam>=2&&eTracksIP>=2&&eTracksAll<=2) {
	if(eTracksBeam>=2&&eTracksIP>=0&&eTracksAll<=2) {
            int ntrfromip=0;

            copy(&bevent);                                //событие
	    bevent.event=kdcenum_.EvNum;                 //take succesive event number from VDDCRec

	    unsigned short nhits=mu_next_event();

	    copy(&bvertex);
	    copy(&bemc);

            for(int t=0; t<eTracksAll; t++)                            //цикл по числу треков  - заполняем однотрековыми
	    {
//		if(sqrt(pow(tX0(t),2)+pow(tY0(t),2)+pow(tZ0(t),2))<30)    {          //cut area in DC       30 - радиус сферы в см
//		if( fabs(tX0(t)-eVertexX)<1 && fabs(tY0(t)-eVertexY)<1 && fabs(tZ0(t)-eVertexZ)<15 )    {          //cut area in DC       30 - радиус сферы в см
		//cout<<"tX0(t)="<<tX0(t)<<"   tY0(t)="<<tY0(t)<<"    tZ0(t)="<<tZ0(t)<<endl;   //coordinates begin track
             	//cout<<"eVertexX="<<eVertexX<<"   eVertexY="<<eVertexY<<"    eVertexZ="<<eVertexZ<<endl;   //координаты восстановленной вершины
                //============================================================
		double xx=tX0IP(t)*tX0IP(t);
		double yy=tY0IP(t)*tY0IP(t);
		double rr=sqrt(xx+yy);
		//cout<<"rr="<<rr<<endl;   //radius begin track

		//if(rr<3) {
		if(rr<5) {
	    	    ntrfromip++;
		    if(ntrfromip<=6){
			listtr[ntrfromip-1]=t+1;                //Список идентификаторов треков ListTr[1:6]
		    }
		}
                //============================================================

		//copy(&btrack,t);
		copy(&btrack[t],t);
		copy(&btof[t],t);
               	//copy(&bemc,t);
		copy(&bmu,t,nhits);
		copy(&batc);
/*                atc_to_track(t);
//		cout<<"Number of counters crossed by the track atctrackinfo.ncnt="<<atctrackinfo.ncnt<<"\n"<<endl;
                for(int k=0; k<atctrackinfo.ncnt; k++)                      //цикл по числу пересеченных счетчиков на трек
		{
                cout<<"Track="<<t<<"\t"<<"atctrackinfo.cnt="<<atctrackinfo.cnt[k]-1<<endl;
//		copy(&bcnt[atctrackinfo.cnt[k]-1],atctrackinfo.cnt[k]-1,t);
		copy(&bcnt[k][t],atctrackinfo.cnt[k]-1,t);
		}
*/
		//for(int k=0; k<Natccr; k++)                    //цикл по числу пересеченных счетчиков на трек
		for(int k=0; k<6; k++)                    //цикл по числу пересеченных счетчиков на трек
		{
//                  cout<<"Event="<<kdcenum_.EvNum<<"\t"<<"Track="<<t<<"\t"<<"i="<<i<<"\t"<<"atc_track.cnt_cross[t][i]="<<atc_track.cnt_cross[t][i]<<"\t"<<"atc_rec.npe="<<atc_rec.npe[atc_track.cnt_cross[t][i]-1]<<endl;
//		    copy(&bcnt[atc_track.cnt_cross[t][i]-1],(atc_track.cnt_cross[t][i]-1),t);
                    //atc_track.cnt_cross[t][i]-1-список пересеченных счетчиков
		    copy(&bcnt[k][t],(atc_track.cnt_cross[t][k]-1),t);
		}

//		ATCTree->Fill();
/*		for(int i=0; i<NATC; i++) {                          //чтобы заполнить ветку для всех счетчиков АЧС
		    copy(&bcnt[i],i,t);
		}
*/
		//trackTree->Fill();

		//eventTree->Fill();
		for(int c=0; c<semc.dc_emc_ncls[t]; c++)                      //dc_emc_ncls[NDCH_TRK] - число кластеров emc, соответствующих данному треку
		{
                    //cout<<"Event="<<kdcenum_.EvNum<<"\t"<<"Track="<<t<<"\t"<<"semc.dc_emc_ncls[t]="<<semc.dc_emc_ncls[t]<<endl;
		    copy(&bcluster[c][t],(semc.dc_emc_cls[t][c]-1)); //?????     //dc_emc_cls[NDCH_TRK][NEMC_CLS] - номера таких кластеров
		    //copy(&bcluster[c],(semc.dc_emc_cls[t][c]-1)); //?????     //dc_emc_cls[NDCH_TRK][NEMC_CLS] - номера таких кластеров
		}
		//emcTowerTree->Fill();


		for(int t=0; t<semc.str_ntracks; t++) {
		    copy(&bstriptrack,t);
		}
		//lkrStripTrackTree->Fill();
//             }        // if( fabs(tX0(t)-eVertexX)<1 && fabs(tY0(t)-eVertexY)<1 && fabs(tZ0(t)-eVertexZ)<15 )    {
	    }

            int nclg=0;
            for(int cl=0; cl<semc.emc_ncls; cl++)       //цикл по числу кластеров
	    {
		if(semc.emc_dc_ntrk[cl]==0)            //если кластер не привязан к треку
		{
		    copy(&bclgamma[nclg],cl);               //	ncls=semc.emc_emc_ncls[c];       ntracks=semc.emc_dc_ntrk[c];
		    nclg++;
		}
                //else
		//{
		//    if(cl=!nclg){copy(&bclgamma[cl],-1);}
		//}
            }

            //cout<<"Event="<<eDaqNumber<<"\t"<<endl;

	    for(int c=0; c<semc.str_ncls; c++) {
		copy(&bstrip,c);
	    }
	    //lkrStripTree->Fill();


            //====================================================================
	    int nclnft=0;
	    for(int cl=0; cl<semc.emc_ncls; cl++)       //цикл по числу кластеров
	    {
		if(semc.emc_dc_ntrk[cl]==0&&            //если кластер не привязан к треку и его энергия больше минимальной (50 МэВ)
		   semc.emc_energy[cl]>EMinPhot){
		    nclnft++;                           //подсчитываем число таких кластеров
		    if(nclnft<=8) {
			lclpi0[nclnft-1]=cl+1;          //Список идентификаторов кластеров для построения  pi0 ListPi0[1:8]

            }
		}
	    }

	    //cout<<"Event="<<eDaqNumber<<"\t"<<"ntrfromip="<<ntrfromip<<"\t"<<"nclnft="<<nclnft<<endl;   //radius begin track
            //cout<<"tP(0)="<<tP(0)<<"\t"<<"tP(1)="<<tP(1)<<endl;

	    jpsi.numHyp=0;
	    for(int i=0; i<5; i++){
	    jpsi.chi2[i]=0;
	    jpsi.M[i]=0;
	    jpsi.P1[i]=0;
	    jpsi.P2[i]=0;
	    }

	    if(ntrfromip==2&&nclnft==0) {             //for selection bhabha events    27/02/2018
//		printf("nclnft==1 event %d: %d tracks, %d from I.P., %d close to beam line\n",eDaqNumber,eTracksAll,eTracksIP,eTracksBeam);
		int qst=0;
		int ntracks=2;
		int npi0=0;
		int ng=0;
		//            lclg[0]=lclpi0[0];
		reckfitmulticpi0gamma(&WTotal,&ntracks,&npi0,&ng,listtr,lclpi0,lclg,&qst);
//		cout<<"qst:"<<qst<<endl;
		if(qst>=0){
//		       cout<<"NumOfHyp:"<<hNumOfHyp<<endl;
		    for(int ih=0; ih<hNumOfHyp; ih++){
//			  cout<<" chi2:"<<hHChi2(ih)<<endl;
                          jpsi.chi2[ih]=hHChi2(ih);
			for(int it=0; it<2; it++){
			    //cout<<"track:"<<it<<"  Mass:"<<hMPc(ih,it)<<"\t"<< "P:"<<hPPc(ih,it)<<endl;
			    if(it==0){
			    jpsi.M[ih]=hMPc(ih,it);
			    jpsi.P1[ih]=hPPc(ih,it);
			    }
			    if(it==1){
			    jpsi.P2[ih]=hPPc(ih,it);
			    }
			}
		    }
		}}

	    if(ntrfromip==2&&nclnft==1) {
//		printf("nclnft==1 event %d: %d tracks, %d from I.P., %d close to beam line\n",eDaqNumber,eTracksAll,eTracksIP,eTracksBeam);
		int qst=0;
		int ntracks=2;
		int npi0=1;
		int ng=0;
		//            lclg[0]=lclpi0[0];
		reckfitmulticpi0gamma(&WTotal,&ntracks,&npi0,&ng,listtr,lclpi0,lclg,&qst);
//		cout<<"qst:"<<qst<<endl;
		if(qst>=0){
//		       cout<<"NumOfHyp:"<<hNumOfHyp<<endl;
		    for(int ih=0; ih<hNumOfHyp; ih++){
//			  cout<<" chi2:"<<hHChi2(ih)<<endl;
                          jpsi.chi2[ih]=hHChi2(ih);
			for(int it=0; it<2; it++){
			    //cout<<"track:"<<it<<"  Mass:"<<hMPc(ih,it)<<"\t"<< "P:"<<hPPc(ih,it)<<endl;
			    if(it==0){
			    jpsi.M[ih]=hMPc(ih,it);
			    jpsi.P1[ih]=hPPc(ih,it);
			    }
			    if(it==1){
			    jpsi.P2[ih]=hPPc(ih,it);
			    }
			}
		    }
		}}

	    if(ntrfromip==2&&nclnft==2) {
//		printf("nclnft==2 event %d: %d tracks, %d from I.P., %d close to beam line\n",eDaqNumber,eTracksAll,eTracksIP,eTracksBeam);
		int qst=0;
		int ntracks=2;
		int npi0=1;
		int ng=0;
		//            lclg[0]=lclpi0[0];
		reckfitmulticpi0gamma(&WTotal,&ntracks,&npi0,&ng,listtr,lclpi0,lclg,&qst);
//		cout<<"qst:"<<qst<<endl;
		if(qst>=0){
//		       cout<<"NumOfHyp:"<<hNumOfHyp<<endl;
		    for(int ih=0; ih<hNumOfHyp; ih++){
//			   cout<<" chi2:"<<hHChi2(ih)<<endl;
			   jpsi.chi2[ih]=hHChi2(ih);
			for(int it=0; it<2; it++){
			    //cout<<"track:"<<it<<"  Mass:"<<hMPc(ih,it)<<"\t"<< "P:"<<hPPc(ih,it)<<endl;
			    if(it==0){
			    jpsi.M[ih]=hMPc(ih,it);
			    jpsi.P1[ih]=hPPc(ih,it);
			    }
			    if(it==1){
			    jpsi.P2[ih]=hPPc(ih,it);
			    }
			}
		    }
		}}

	    if(ntrfromip==2&&nclnft==3) {
//		printf("nclnft==3 event %d: %d tracks, %d from I.P., %d close to beam line\n",eDaqNumber,eTracksAll,eTracksIP,eTracksBeam);
		int qst=0;
		int ntracks=2;
		int npi0=1;
		int ng=0;
		//            lclg[0]=lclpi0[0];
		reckfitmulticpi0gamma(&WTotal,&ntracks,&npi0,&ng,listtr,lclpi0,lclg,&qst);
//		cout<<"qst:"<<qst<<endl;
		if(qst>=0){
//		       cout<<"NumOfHyp:"<<hNumOfHyp<<endl;
		    for(int ih=0; ih<hNumOfHyp; ih++){
//			   cout<<" chi2:"<<hHChi2(ih)<<endl;
			   jpsi.chi2[ih]=hHChi2(ih);
			for(int it=0; it<2; it++){
			    //cout<<"track:"<<it<<"  Mass:"<<hMPc(ih,it)<<"\t"<< "P:"<<hPPc(ih,it)<<endl;
			    if(it==0){
			    jpsi.M[ih]=hMPc(ih,it);
			    jpsi.P1[ih]=hPPc(ih,it);
			    }
			    if(it==1){
			    jpsi.P2[ih]=hPPc(ih,it);
			    }
			}
		    }
		}}

	    if(ntrfromip==2&&nclnft==4) {
//		printf("nclnft==4 event %d: %d tracks, %d from I.P., %d close to beam line\n",eDaqNumber,eTracksAll,eTracksIP,eTracksBeam);
		int qst=0;
		int ntracks=2;
		int npi0=1;
		int ng=0;
		//            lclg[0]=lclpi0[0];
		reckfitmulticpi0gamma(&WTotal,&ntracks,&npi0,&ng,listtr,lclpi0,lclg,&qst);
//		cout<<"qst:"<<qst<<endl;
		if(qst>=0){
//		       cout<<"NumOfHyp:"<<hNumOfHyp<<endl;
		    for(int ih=0; ih<hNumOfHyp; ih++){
//			   cout<<" chi2:"<<hHChi2(ih)<<endl;
			   jpsi.chi2[ih]=hHChi2(ih);
			for(int it=0; it<2; it++){
			    //cout<<"track:"<<it<<"  Mass:"<<hMPc(ih,it)<<"\t"<< "P:"<<hPPc(ih,it)<<endl;
			    if(it==0){
			    jpsi.M[ih]=hMPc(ih,it);
			    jpsi.P1[ih]=hPPc(ih,it);
			    }
			    if(it==1){
			    jpsi.P2[ih]=hPPc(ih,it);
			    }
			}
		    }
		}}

	    if(ntrfromip==2&&nclnft==5) {
//		printf("nclnft==5 event %d: %d tracks, %d from I.P., %d close to beam line\n",eDaqNumber,eTracksAll,eTracksIP,eTracksBeam);
		int qst=0;
		int ntracks=2;
		int npi0=1;
		int ng=0;
		//            lclg[0]=lclpi0[0];
		reckfitmulticpi0gamma(&WTotal,&ntracks,&npi0,&ng,listtr,lclpi0,lclg,&qst);
//		cout<<"qst:"<<qst<<endl;
		if(qst>=0){
//		       cout<<"NumOfHyp:"<<hNumOfHyp<<endl;
		    for(int ih=0; ih<hNumOfHyp; ih++){
//			cout<<" chi2:"<<hHChi2(ih)<<endl;
			jpsi.chi2[ih]=hHChi2(ih);
			for(int it=0; it<2; it++){
			    //cout<<"track:"<<it<<"  Mass:"<<hMPc(ih,it)<<"\t"<< "P:"<<hPPc(ih,it)<<endl;
			    if(it==0){
			    jpsi.M[ih]=hMPc(ih,it);
			    jpsi.P1[ih]=hPPc(ih,it);
			    }
			    if(it==1){
			    jpsi.P2[ih]=hPPc(ih,it);
			    }
			}
		    }
		}}

               if(ntrfromip==2&&nclnft==6) {
//		printf("nclnft==5 event %d: %d tracks, %d from I.P., %d close to beam line\n",eDaqNumber,eTracksAll,eTracksIP,eTracksBeam);
		int qst=0;
		int ntracks=2;
		int npi0=1;
		int ng=0;
		//            lclg[0]=lclpi0[0];
		reckfitmulticpi0gamma(&WTotal,&ntracks,&npi0,&ng,listtr,lclpi0,lclg,&qst);
//		cout<<"qst:"<<qst<<endl;
		if(qst>=0){
//		       cout<<"NumOfHyp:"<<hNumOfHyp<<endl;
		    for(int ih=0; ih<hNumOfHyp; ih++){
//			cout<<" chi2:"<<hHChi2(ih)<<endl;
			jpsi.chi2[ih]=hHChi2(ih);
			for(int it=0; it<2; it++){
			    //cout<<"track:"<<it<<"  Mass:"<<hMPc(ih,it)<<"\t"<< "P:"<<hPPc(ih,it)<<endl;
			    if(it==0){
			    jpsi.M[ih]=hMPc(ih,it);
			    jpsi.P1[ih]=hPPc(ih,it);
			    }
			    if(it==1){
			    jpsi.P2[ih]=hPPc(ih,it);
			    }
			}
		    }
		}}


	    if(ntrfromip==1&&nclnft==1) {
//		printf("nclnft==1 event %d: %d tracks, %d from I.P., %d close to beam line\n",eDaqNumber,eTracksAll,eTracksIP,eTracksBeam);
		int qst=0;
		int ntracks=2;
		int npi0=1;
		int ng=0;
		//            lclg[0]=lclpi0[0];
		reckfitmulticpi0gamma(&WTotal,&ntracks,&npi0,&ng,listtr,lclpi0,lclg,&qst);
//		cout<<"qst:"<<qst<<endl;
		if(qst>=0){
//		       cout<<"NumOfHyp:"<<hNumOfHyp<<endl;
		    for(int ih=0; ih<hNumOfHyp; ih++){
//			  cout<<" chi2:"<<hHChi2(ih)<<endl;
                          jpsi.chi2[ih]=hHChi2(ih);
			for(int it=0; it<2; it++){
			    //cout<<"track:"<<it<<"  Mass:"<<hMPc(ih,it)<<"\t"<< "P:"<<hPPc(ih,it)<<endl;
			    if(it==0){
			    jpsi.M[ih]=hMPc(ih,it);
			    jpsi.P1[ih]=hPPc(ih,it);
			    }
			    if(it==1){
			    jpsi.P2[ih]=hPPc(ih,it);
			    }
			}
		    }
		}}

	    if(ntrfromip==1&&nclnft==2) {
//		printf("nclnft==2 event %d: %d tracks, %d from I.P., %d close to beam line\n",eDaqNumber,eTracksAll,eTracksIP,eTracksBeam);
		int qst=0;
		int ntracks=2;
		int npi0=1;
		int ng=0;
		//            lclg[0]=lclpi0[0];
		reckfitmulticpi0gamma(&WTotal,&ntracks,&npi0,&ng,listtr,lclpi0,lclg,&qst);
//		cout<<"qst:"<<qst<<endl;
		if(qst>=0){
//		       cout<<"NumOfHyp:"<<hNumOfHyp<<endl;
		    for(int ih=0; ih<hNumOfHyp; ih++){
//			   cout<<" chi2:"<<hHChi2(ih)<<endl;
			   jpsi.chi2[ih]=hHChi2(ih);
			for(int it=0; it<2; it++){
			    //cout<<"track:"<<it<<"  Mass:"<<hMPc(ih,it)<<"\t"<< "P:"<<hPPc(ih,it)<<endl;
			    if(it==0){
			    jpsi.M[ih]=hMPc(ih,it);
			    jpsi.P1[ih]=hPPc(ih,it);
			    }
			    if(it==1){
			    jpsi.P2[ih]=hPPc(ih,it);
			    }
			}
		    }
		}}

	    if(ntrfromip==1&&nclnft==3) {
//		printf("nclnft==3 event %d: %d tracks, %d from I.P., %d close to beam line\n",eDaqNumber,eTracksAll,eTracksIP,eTracksBeam);
		int qst=0;
		int ntracks=2;
		int npi0=1;
		int ng=0;
		//            lclg[0]=lclpi0[0];
		reckfitmulticpi0gamma(&WTotal,&ntracks,&npi0,&ng,listtr,lclpi0,lclg,&qst);
//		cout<<"qst:"<<qst<<endl;
		if(qst>=0){
//		       cout<<"NumOfHyp:"<<hNumOfHyp<<endl;
		    for(int ih=0; ih<hNumOfHyp; ih++){
//			   cout<<" chi2:"<<hHChi2(ih)<<endl;
			   jpsi.chi2[ih]=hHChi2(ih);
			for(int it=0; it<2; it++){
			    //cout<<"track:"<<it<<"  Mass:"<<hMPc(ih,it)<<"\t"<< "P:"<<hPPc(ih,it)<<endl;
			    if(it==0){
			    jpsi.M[ih]=hMPc(ih,it);
			    jpsi.P1[ih]=hPPc(ih,it);
			    }
			    if(it==1){
			    jpsi.P2[ih]=hPPc(ih,it);
			    }
			}
		    }
		}}

	    if(ntrfromip==1&&nclnft==4) {
//		printf("nclnft==3 event %d: %d tracks, %d from I.P., %d close to beam line\n",eDaqNumber,eTracksAll,eTracksIP,eTracksBeam);
		int qst=0;
		int ntracks=2;
		int npi0=1;
		int ng=0;
		//            lclg[0]=lclpi0[0];
		reckfitmulticpi0gamma(&WTotal,&ntracks,&npi0,&ng,listtr,lclpi0,lclg,&qst);
//		cout<<"qst:"<<qst<<endl;
		if(qst>=0){
//		       cout<<"NumOfHyp:"<<hNumOfHyp<<endl;
		    for(int ih=0; ih<hNumOfHyp; ih++){
//			   cout<<" chi2:"<<hHChi2(ih)<<endl;
			   jpsi.chi2[ih]=hHChi2(ih);
			for(int it=0; it<2; it++){
			    //cout<<"track:"<<it<<"  Mass:"<<hMPc(ih,it)<<"\t"<< "P:"<<hPPc(ih,it)<<endl;
			    if(it==0){
			    jpsi.M[ih]=hMPc(ih,it);
			    jpsi.P1[ih]=hPPc(ih,it);
			    }
			    if(it==1){
			    jpsi.P2[ih]=hPPc(ih,it);
			    }
			}
		    }
		}}

          	if(ntrfromip==1&&nclnft==5) {
//		printf("nclnft==3 event %d: %d tracks, %d from I.P., %d close to beam line\n",eDaqNumber,eTracksAll,eTracksIP,eTracksBeam);
		int qst=0;
		int ntracks=2;
		int npi0=1;
		int ng=0;
		//            lclg[0]=lclpi0[0];
		reckfitmulticpi0gamma(&WTotal,&ntracks,&npi0,&ng,listtr,lclpi0,lclg,&qst);
//		cout<<"qst:"<<qst<<endl;
		if(qst>=0){
//		       cout<<"NumOfHyp:"<<hNumOfHyp<<endl;
		    for(int ih=0; ih<hNumOfHyp; ih++){
//			   cout<<" chi2:"<<hHChi2(ih)<<endl;
			   jpsi.chi2[ih]=hHChi2(ih);
			for(int it=0; it<2; it++){
			    //cout<<"track:"<<it<<"  Mass:"<<hMPc(ih,it)<<"\t"<< "P:"<<hPPc(ih,it)<<endl;
			    if(it==0){
			    jpsi.M[ih]=hMPc(ih,it);
			    jpsi.P1[ih]=hPPc(ih,it);
			    }
			    if(it==1){
			    jpsi.P2[ih]=hPPc(ih,it);
			    }
			}
		    }
		}}

                if(ntrfromip==1&&nclnft==6) {
//		printf("nclnft==3 event %d: %d tracks, %d from I.P., %d close to beam line\n",eDaqNumber,eTracksAll,eTracksIP,eTracksBeam);
		int qst=0;
		int ntracks=2;
		int npi0=1;
		int ng=0;
		//            lclg[0]=lclpi0[0];
		reckfitmulticpi0gamma(&WTotal,&ntracks,&npi0,&ng,listtr,lclpi0,lclg,&qst);
//		cout<<"qst:"<<qst<<endl;
		if(qst>=0){
//		       cout<<"NumOfHyp:"<<hNumOfHyp<<endl;
		    for(int ih=0; ih<hNumOfHyp; ih++){
//			   cout<<" chi2:"<<hHChi2(ih)<<endl;
			   jpsi.chi2[ih]=hHChi2(ih);
			for(int it=0; it<2; it++){
			    //cout<<"track:"<<it<<"  Mass:"<<hMPc(ih,it)<<"\t"<< "P:"<<hPPc(ih,it)<<endl;
			    if(it==0){
			    jpsi.M[ih]=hMPc(ih,it);
			    jpsi.P1[ih]=hPPc(ih,it);
			    }
			    if(it==1){
			    jpsi.P2[ih]=hPPc(ih,it);
			    }
			}
		    }
		}}

	    //hpx->Fill(hNumOfHyp);
	    jpsi.numHyp=hNumOfHyp;
	    //jpsiTree->Fill();

                if(eNumber%10==0) cout<<"Ev:"<<eNumber<<endl;
	    //==================================================================

            eventTree->Fill();


	    if( progpar.process_only ) return 1;

	    return 0;

	}    //if(eTracksBeam>=2&&eTracksIP>=2&&eTracksAll<=2) {
}


//static const char* optstring="ra:b:p:t:e:c:k:i:o:v:x";
static const char* optstring="ra:d:b:p:h:s:j:t:e:c:l:k:i:o:v:n:x";

void Usage(int status)
{                   //static const struct ProgramParameters def_progpar={false,1,0,0,50,300,1,0,0,"/store/users/ovtin/out.root",false};
	cout.precision(3);
	cout<<"Usage:\n"
		<<progname<<" [OPTIONS] run|nat_file|run_list_file...\n\n"
		<<"Reconstruction of BhaBha events with DC and search of hits in ATC system.\n"
		<<"Options:\n"
		<<"  -r             Read VDDC&EMC reconstructed info from file (reconstruct by default)\n"
		<<"  -a tracks      Minimum total tracks number required (default to "<<def_progpar.min_track_number<<")\n"
		<<"  -d tracks      Maximum total tracks number required (default to "<<def_progpar.max_track_number<<")\n"
		<<"  -b tracks      Minimum beam tracks number required (default to "<<def_progpar.min_beam_track_number<<")\n"
		<<"  -p tracks      Minumum IP tracks number required (default to "<<def_progpar.min_ip_track_number<<")\n"
		<<"  -h tracks      Maximum IP tracks number required (default to "<<def_progpar.max_ip_track_number<<")\n"
            	<<"  -s Momentum    Lower limit of momentum cut, MeV/c (default: "<<def_progpar.min_momentum<<")\n"
            	<<"  -j Momentum    Maximum limit of momentum cut, MeV/c (default: "<<def_progpar.max_momentum<<")\n"
 	        <<"  -t energy      EMC cluster energy threshold (default to "<<def_progpar.min_cluster_energy<<" MeV)\n"
 		<<"  -e energy      Minumum total energy in EMC (default to "<<def_progpar.min_total_energy<<" MeV)\n"
		<<"  -c clusters    Minumum number of clusters in both calorimeters (default to "<<def_progpar.min_cluster_number<<")\n"
		<<"  -l clusters    Maximum number of clusters in both calorimeters (default to "<<def_progpar.max_cluster_number<<")\n"
		<<"  -k clusters    Minumum number of clusters in LKr calorimeter (default to "<<def_progpar.min_lkrcl_number<<")\n"
		<<"  -i clusters    Minumum number of clusters in CsI calorimeter (default to "<<def_progpar.min_csicl_number<<")\n"
	        <<"  -o RootFile    Output ROOT file name (default to "<<def_progpar.rootfile<<")\n"
            	<<"  -v MCCalibRunNumber    MCCalibRunNumber (default to "<<def_progpar.MCCalibRunNumber<<")\n"
            	<<"  -n NEvents     Number events in process "<<def_progpar.NEvents<<")\n"
		<<"  -x             Process the events specified after file exclusively and print debug information"
		<<endl;
	exit(status);
}

int main(int argc, char* argv[])
{
	progname=argv[0];

	//If no arguments print usage help
	if( argc==1 ) Usage(0);

	int opt;
//----------------- Process options -----------------//
	while( (opt=getopt(argc,argv,optstring))>0 ) {
		switch( opt ) {
			case '?': Usage(1); break;
			case 'r': progpar.read_reco=true; break;
			case 'a': progpar.min_track_number=atoi(optarg); break;
			case 'd': progpar.max_track_number=atoi(optarg); break;
			case 'b': progpar.min_beam_track_number=atoi(optarg); break;
			case 'p': progpar.min_ip_track_number=atoi(optarg); break;
			case 'h': progpar.max_ip_track_number=atoi(optarg); break;
			case 's': progpar.min_momentum=atof(optarg); break;
			case 't': progpar.min_cluster_energy=atof(optarg); break;
			case 'e': progpar.min_total_energy=atof(optarg); break;
			case 'c': progpar.min_cluster_number=atoi(optarg); break;
			case 'l': progpar.max_cluster_number=atoi(optarg); break;
			case 'k': progpar.min_lkrcl_number=atoi(optarg); break;
			case 'i': progpar.min_csicl_number=atoi(optarg); break;
		        case 'o': progpar.rootfile=optarg; break;
                        case 'v': progpar.MCCalibRunNumber=atoi(optarg); break;
                        case 'n': progpar.NEvents=atoi(optarg); break;
			case 'x': progpar.process_only=true; break;
			default : Usage(1);
		}
	}

	if( progpar.min_track_number<progpar.min_beam_track_number ||
		progpar.min_beam_track_number<progpar.min_ip_track_number ) {
		cerr<<"Error in parameters specification, should be MinTotalTracks>=MinBeamTracks>=MinIPTracks"<<endl;
		return 1;
	}


	char timestr[161];
	time_t curtime=time(NULL);

	strftime(timestr,161,"%F %T %z",localtime(&curtime));
	cout<<"Current time "<<timestr<<endl;
	cout<<" Writing output trees to "<<progpar.rootfile<<endl;


//----------------- Initialize ROOT file and trees -----------------//
	TFile *fout=0;
	//Create root file if exclusive event processing is not set
	if( !progpar.process_only )
		fout = new TFile(progpar.rootfile,"RECREATE");        //fout = new TFile(progpar.rootfile,"RECREATE","",7);

	eventTree = new TTree("et","Event tree");
	eventTree->SetAutoSave(500000000);  // autosave when 0.5 Gbyte written
	eventTree->Branch("ev",&bevent,eventBranchList);
        eventTree->Branch("vrt",&bvertex,vertexBranchList);
	eventTree->Branch("emc",&bemc,emcBranchList);
	eventTree->Branch("atcev",&batc,atcBranchList);

	//for(int i=0; i<Ntraks; i++) {
	for(int i=0; i<3; i++) {
	//char branchname[Ntraks];
	char branchname[3];
	sprintf(branchname,"t%d",i);
	TBranch* b1=eventTree->Branch(branchname,&btrack[i],trackBranchList);
	  //for(int ii=0; ii<Natccr; ii++) {
	  for(int ii=0; ii<6; ii++) {
	      //char branchname2[Natccr];
	      char branchname2[6];
	      sprintf(branchname2,"t%datccr%d",i,ii);
	      eventTree->Branch(branchname2,&bcnt[ii][i],atcCounterBranchList);
	  }
	  //for (int ii=0; ii<Nclcr; ii++)
	  for (int ii=0; ii<3; ii++)
	  {
	      //char clastername[Nclcr];
	      char clastername[3];
	      sprintf(clastername,"t%dc%d",i,ii);
	      eventTree->Branch(clastername,&bcluster[ii][i],towerClusterBranchList);
	  }
	  //char tofname[Ntraks];
	  char tofname[3];
	  sprintf(tofname,"t%dtof",i);
          eventTree->Branch(tofname,&btof[i],ToFBranchList);
	}

	//for (int ii=0; ii<Nclcr; ii++)
	for (int ii=0; ii<4; ii++)
	{
	    //char clgammaname[Nclcr];
	    char clgammaname[4];
	    sprintf(clgammaname,"clgamma%d",ii);
	    eventTree->Branch(clgammaname,&bclgamma[ii],towerClusterBranchList);
	}

	eventTree->Branch("mu",&bmu,MUBranchList);
	eventTree->Branch("jpsi",&jpsi,"numHyp/I:chi2[5]/F:M[5]:P1[5]:P2[5]");

	eventTree->Branch("strcls",&bstrip,stripClusterBranchList);
	eventTree->Branch("strtrk",&bstriptrack,stripTrackBranchList);

	//trackTree = new TTree("tt","Track tree");
	//trackTree->Branch("t",&btrack,trackBranchList);
	//trackTree->Branch("tof",&btof,ToFBranchList);
	//trackTree->Branch("mu",&bmu,MUBranchList);
	//trackTree->Branch("atc",&batc,atcBranchList);
/*        for(int i=0; i<NATC; i++) {
	char branchname[161];
        sprintf(branchname,"cnt%d",i);
	trackTree->Branch(branchname,&bcnt[i],atcCounterBranchList);
	}
*/
/*
	for(int i=0; i<6; i++) {
	char branchname[6];
//        sprintf(branchname,"crcnt%d",i);
        sprintf(branchname,"cnt%d",i);
	trackTree->Branch(branchname,&bcnt[i],atcCounterBranchList);
	}
*/
/*	ATCTree = new TTree("ATC","ATC tree");
	ATCTree->Branch("atc",&batc,atcBranchList);
        for(int i=0; i<NATC; i++) {                          //чтобы заполнить ветку для всех счетчиков
	char branchname[161];
        sprintf(branchname,"cnt%d",i);
	ATCTree->Branch(branchname,&bcnt[i],atcCounterBranchList);
	}
	*/
        /*
	emcTowerTree = new TTree("emct","EMC tower cluster tree");
	emcTowerTree->Branch("ev",&bevent,eventBranchList);
	emcTowerTree->Branch("emc",&bemc,emcBranchList);
	emcTowerTree->Branch("vrt",&bvertex,vertexBranchList);
	for (int i=0; i<2; i++)
	{
	    char clastername[3];
            sprintf(clastername,"c%d",i);
	    emcTowerTree->Branch(clastername,&bcluster[i],towerClusterBranchList);
	}
        */
        /*
	lkrStripTree = new TTree("strclt","LKr strip cluster tree");
	lkrStripTree->Branch("strcls",&bstrip,stripClusterBranchList);

	lkrStripTrackTree = new TTree("strtrt","LKr strip cluster tree");
	lkrStripTrackTree->Branch("strtrk",&bstriptrack,stripTrackBranchList);
        */
	//jpsiTree = new TTree("jpsi","pi+pi-pi0 tree");
	//jpsiTree->Branch("jpsi",&jpsi,"numHyp/I:chi2[5]/F:M[5]:P1[5]:P2[5]");
	//piTree->Branch("pi",&bpi,piBranchList);
	//hpx=new TH1F("hyp","This is the hyp distribution",6,0,6);

//----------------- Configure kframework -----------------//
	//Set kframework signal handling
	kf_install_signal_handler(1);

	//Set subsystems to be used
//	kf_use(KF_VDDC_SYSTEM|KF_EMC_SYSTEM);
        kf_use(KF_VDDC_SYSTEM|KF_TOF_SYSTEM|KF_ATC_SYSTEM|KF_EMC_SYSTEM|KF_MU_SYSTEM);                             //устанавливаем, какие будем использовать системы

	//Register to kframework used cuts
	char buf[100];
	kf_add_cut(KF_PRE_SEL,CutLongDcRecord,"DC event length >1940 words");
	kf_add_cut(KF_PRE_SEL,CutLongVdRecord,"VD event length >156 words");
	//sprintf(buf,"raw total EMC energy <%d MeV",progpar.min_total_energy);
	//kf_add_cut(KF_PRE_SEL,CutTotalEnergyFast,buf);

	//kf_add_cut(KF_VDDC_SEL,OneTrackCut,"1 cosmic track");
	sprintf(buf,"momentum <= %5.fMeV/c",progpar.min_momentum);
	kf_add_cut(KF_VDDC_SEL,MinMomentumCut,buf);
	sprintf(buf,"momentum >= %5.fMeV/c",progpar.max_momentum);
	kf_add_cut(KF_VDDC_SEL,MaxMomentumCut,buf);
	//kf_add_cut(KF_VDDC_SEL,NotVertex,"Track not from vertex: tVertex(t)!=1");
	kf_add_cut(KF_VDDC_SEL,ChargeCut,"Charge from track: (tCharge(0)+tCharge(1))!=0");
	//kf_add_cut(KF_VDDC_SEL,VectorsZCut,"Number Z vectors: tVectorsZ(0)<1 || tVectorsZ(1)<1");
	//kf_add_cut(KF_VDDC_SEL,HitsVDCut,"Number hits on track VD: tHitsVD(0)<1");
        //kf_add_cut(KF_VDDC_SEL, tX0Cut,"tX0Cut:");
        //kf_add_cut(KF_VDDC_SEL, tY0Cut,"tY0Cut:");
        //kf_add_cut(KF_VDDC_SEL, tZ0Cut,"tX0Cut:");
       	// kf_add_cut(KF_VDDC_SEL,MinXYVectorsCut,"XY vectors <3");
	//kf_add_cut(KF_VDDC_SEL,Chi2Cut,"bad chi-square");
	sprintf(buf,"minimum number of tracks equally %d and maximum number of tracks equally %d<",progpar.min_track_number,progpar.max_track_number);
	kf_add_cut(KF_VDDC_SEL,CutTrackNumber,buf);
	sprintf(buf,"number of beam tracks <%d",progpar.min_beam_track_number);
	kf_add_cut(KF_VDDC_SEL,CutBeamTrackNumber,buf);
	sprintf(buf,"number of IP tracks <%d",progpar.min_ip_track_number);
	kf_add_cut(KF_VDDC_SEL,CutIPTrackNumber_min,"IP tracks cut");
	sprintf(buf,"number of IP tracks >%d",progpar.max_ip_track_number);
	kf_add_cut(KF_VDDC_SEL,CutIPTrackNumber_max,"IP tracks cut");
	//kf_add_cut(KF_VDDC_SEL,EMCthetaCut,"EMCtheta cut");
	//kf_add_cut(KF_VDDC_SEL,Cut_CPhi2t,"CPhi2t cut");
	//kf_add_cut(KF_VDDC_SEL,Cut_CThe2t,"CThe2t cut");

        kf_add_cut(KF_ATC_SEL,AtcEventDamageCut,"damaged record");
        kf_add_cut(KF_ATC_SEL,AtcIllegalTrackCut,"illegal track");

	//    kf_add_cut(KF_ATC_SEL,NoAtcOnTrackCut,"no counters on track");

	sprintf(buf,"total EMC energy <%5.fMeV",progpar.min_total_energy);
	kf_add_cut(KF_EMC_SEL,CutTotalEnergy,buf);
	sprintf(buf,"cluster EMC energy <%5.fMeV",progpar.min_cluster_energy);
	kf_add_cut(KF_EMC_SEL,CutClusterEnergy,buf);
	//kf_add_cut(KF_EMC_SEL,Cut_10per_Energy,"released energy in EMC is not attributed to clusters >10%d of the total EMC energy");
	sprintf(buf,"number of clusters <%d",progpar.min_cluster_number);
	kf_add_cut(KF_EMC_SEL,CutMINClusterNumber,buf);
	sprintf(buf,"number of clusters >%d",progpar.max_cluster_number);
	kf_add_cut(KF_EMC_SEL,CutMAXClusterNumber,buf);
	//sprintf(buf,"number of LKr clusters <%d",progpar.min_lkrcl_number);
	//kf_add_cut(KF_EMC_SEL,CutLkrClusterNumber,buf);
	//sprintf(buf,"number of CsI clusters <%d",progpar.min_csicl_number);
	//kf_add_cut(KF_EMC_SEL,CutCsiClusterNumber,buf);
	kf_add_cut(KF_EMC_SEL,Cut_nclntrack,"nclntrack cut");

        //kf_add_cut(KF_TOF_SEL,TofCut,"TOF system cut");
        //kf_add_cut(KF_MU_SEL,MUCut,"Mu system cut");

	//Register selection routines
	kf_register_selection(KF_PRE_SEL,pre_event_rejection);
	kf_register_selection(KF_VDDC_SEL,vddc_event_rejection);
	kf_register_selection(KF_TOF_SEL,tof_event_rejection);
        kf_register_selection(KF_ATC_SEL,atc_rejection);
	kf_register_selection(KF_EMC_SEL,emc_event_rejection);
	kf_register_selection(KF_MU_SEL,mu_event_rejection);

        kf_MCCalibRunNumber(progpar.MCCalibRunNumber);

	//Set automatic cosmic run determination
//	kf_cosmic(-1);  //auto
	kf_cosmic(0);  //beam

	kf_modify_header(1);     //Modify header flag. REDE energy read from DB will be written to header

	//Register an analysis routine
	kf_register_analysis(analyse_event);
        cout<<"analyse_event="<<analyse_event<<endl;

	//Do not reconstruct, read reconstruction records from file
	kf_reco_from_file(progpar.read_reco);

	//Set exclusive event processing
	kf_process_only(progpar.process_only);

	TBenchmark *benchmark=new TBenchmark;
	cout<<" Starting benchmark test"<<endl;
	benchmark->Start("test");

	//Call analysis job
	cout<<"NEvents="<<progpar.NEvents<<"\t"<<"argv[optind]="<<argv[optind]<<"\t"<<"&argv[optind]="<<&argv[optind]<<"\t"<<"argc-optind="<<argc-optind<<endl;               //argv[optind]:/space/runs/daq021949.nat.bz2

//	kf_process(argc-optind,&argv[optind],0);
	kf_process(argc-optind,&argv[optind],progpar.NEvents);             //установка ограничения на число обрабатываемых событий

	benchmark->Show("test");

	if( fout ) {
		fout->Write();         //пишем данные в файл и закрываем его
		fout->Close();
	}

	return 0;
}
