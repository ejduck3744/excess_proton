#ifndef Shift_hh
#define Shift_hh
//
#include "StRoot/StEpdUtil/StEpdEpFinder.h"
#include "StEpdUtil/StEpdGeom.h"
#include "StEpdUtil/StEpdEpInfo.h"
#include "StMaker.h"
#include <string>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TProfile3D.h>
#include <TRandom3.h>
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"
//
class StPicoDst;
class StPicoDstMaker;
class StPicoEvent;
class StPicoTrack;
class StPileupUtil;
class TFile;
class TTree;
class TH1;
class TH2;
class TProfile;
class TProfile2D;
class TProfile3D;

const Int_t CONST_VZ_BINS = 4; //epd vz bins {-145,-70,0,70,145}
// const Int_t CONST_VZ_BINS = 2; //epd vz bins {positive, negative}
#ifndef ST_NO_NAMESPACES
using std::string;
#endif
//
//  The class declaration. It innherits from StMaker.
class Shift : public StMaker {

public:
  Shift( const Char_t *name, StPicoDstMaker *picoMaker, const Char_t *jobid );   // constructor
  virtual ~Shift(){};                                 // destructor

  virtual void Clear(Option_t *option=""); // called after every event to cleanup
  virtual Int_t  Init();                   // called once at the beginning of your job
  virtual Int_t  Make();                   // invoked for every event
  virtual Int_t  Finish();                 // called once at the end


  // My functions
  Int_t GetRunIndex( const Int_t run );
  Int_t Centrality(Int_t gRefMult);
  Double_t Pileup_rejection(Double_t gRefmult, Int_t gNtofMatch, Double_t gvz,const StPicoEvent *event);
  Double_t proton_efficiency_correction(const StPicoTrack *ptrk,Int_t cent);
  bool zero2twopi(Double_t &tphi);
  bool   isGoodTrack(const StPicoTrack *ptrk);
  bool   isGoodEvent(const StPicoEvent *event);
  bool   isGoodTrigger(const StPicoEvent *event);
  bool	 isProton_test(const StPicoTrack *ptrk, bool req_tof = true, int option = 0);
  bool	 isProton(const StPicoTrack *ptrk, bool req_tof = true, int option = 0);
  bool   isPion(const StPicoTrack *ptrk);
  bool   isKaon(const StPicoTrack *ptrk);
  int    EP_group(int run_number, int cent);
  //double phi_correction(double phi, int which_one);//which one being which sub event plane.
  int    momindex(Double_t mom);

private:
	//because I don't know how to use phys_constants:
	double protonMass;
	double pionMass;
	double kaonMass;
	double y_beam;
	// double y_cm;
	// double y_beam;
	TRandom3 *r1;
	
	
	////////////the new stuff i'm adding
	//primary vertex
	TH2D *hist_VyVx;
	TH1D *hist_Vz;
	TH1D *baryon_yield_histogram;
	//refmult
	TH1D *hist_refmult;
	TH2D *nhits_vs_eta;
	TH2D *refmult_v_tofmult;
	
	//event planes 
	//epd
	TH1D *hist_epd_ab_psi_raw[2];
	TH1D *hist_epd_c_psi_raw[2];
	TH1D *hist_epd_d_psi_raw[2];
	TH1D *hist_epd_abcd_psi_raw[2];
	TH1D *hist_epd_cd_psi_raw[2];
	TH1D *hist_epd_ab_psi_raw_ew;
	TH1D *hist_epd_psi_raw_ew;
	
	TH1D *hist_epd_ab_psi_recentered[2];
	TH1D *hist_epd_c_psi_recentered[2];
	TH1D *hist_epd_d_psi_recentered[2];
	TH1D *hist_epd_abcd_psi_recentered[2];
	TH1D *hist_epd_cd_psi_recentered[2];
	TH1D *hist_epd_ab_psi_recentered_ew;
	TH1D *hist_epd_psi_recentered_ew;
	
	
	TH1D *hist_epd_ab_psi_flattened[2];
	TH1D *hist_epd_c_psi_flattened[2];
	TH1D *hist_epd_d_psi_flattened[2];
	TH1D *hist_epd_abcd_psi_flattened[2];
	TH1D *hist_epd_cd_psi_flattened[2];
	TH1D *hist_epd_ab_psi_flattened_ew;
	TH1D *hist_epd_psi_flattened_ew;
	
	TProfile *v1_vs_eta_east[9][CONST_VZ_BINS];
	TProfile *v1_vs_eta_west[9][CONST_VZ_BINS];
	
	TProfile* epd_eta_weight_east;
	TProfile* epd_eta_weight_east_cent[9][CONST_VZ_BINS];
	TProfile* epd_eta_weight_west;
	TProfile* epd_eta_weight_west_cent[9][CONST_VZ_BINS];
	
	TProfile* nsigma_weight_proton;
	TProfile* nsigma_weight_aproton;
	
	
	// TProfile *profile_epd_qx[4];
	// TProfile *profile_epd_qy[4];
	
	//tpc
	
	TH1D *hist_tpc_A_psi_raw;
	TH1D *hist_tpc_B_psi_raw;
	
	TH1D *hist_tpc_A_psi_recentered;
	TH1D *hist_tpc_B_psi_recentered;
	TH1D *hist_tpc_AB_psi_recentered;
	
	TH1D * hist_tpc_A_psi_flattened;
	TH1D * hist_tpc_B_psi_flattened;
	TH1D * hist_tpc_AB_psi_flattened;
	
	TProfile *v1_vs_eta_tpc[9][CONST_VZ_BINS];
	
	TProfile* tpc_eta_weight_cent[9][CONST_VZ_BINS];
	
	//correlation plots
	TProfile *profile_correlation_abc;
	TProfile *profile_correlation_abd;
	TProfile *profile_correlation_cd;
	
	TProfile *profile_correlation_eABtA;
	TProfile *profile_correlation_eCtA;
	TProfile *profile_correlation_eDtA;
	
	TProfile *profile_correlation_eABtB;
	TProfile *profile_correlation_eCtB;
	TProfile *profile_correlation_eDtB;
	
	TProfile *profile_correlation_tpc_AB;
	
	TProfile *profile_correlation_abew;
	TProfile *profile_correlation_cew;
	TProfile *profile_correlation_dew;
	TProfile *profile_correlation_ew;
	TProfile *profile_correlation_tpc_ew_vz[CONST_VZ_BINS];
	TProfile *profile_correlation_ew_vz[CONST_VZ_BINS];
	TProfile *profile_correlation_eABw_vz[CONST_VZ_BINS];
	TProfile *profile_correlation_eCDw_vz[CONST_VZ_BINS];
	TProfile *profile_correlation_eABeCD_vz[CONST_VZ_BINS];
	TProfile *profile_correlation_ewAB_vz[CONST_VZ_BINS];
	TProfile *profile_correlation_ewCD_vz[CONST_VZ_BINS];
	TProfile *profile_correlation_wABwCD_vz[CONST_VZ_BINS];
	TProfile *profile_correlation_eDw_vz[CONST_VZ_BINS];
	TProfile *profile_correlation_eABeD_vz[CONST_VZ_BINS];
	TProfile *profile_correlation_ewD_vz[CONST_VZ_BINS];
	TProfile *profile_correlation_wABwD_vz[CONST_VZ_BINS];
	TProfile *profile_correlation_TPCEPDe_vz[CONST_VZ_BINS];
	TProfile *profile_correlation_TPCEPDw_vz[CONST_VZ_BINS];
	
	
	
	TH2D *correlation_abc;
	TH2D *correlation_abd;
	TH2D *correlation_cd;
	
	TH2D *correlation_eABtA;
	TH2D *correlation_eCtA;
	TH2D *correlation_eDtA;
	
	TH2D *correlation_eABtB;
	TH2D *correlation_eCtB;
	TH2D *correlation_eDtB;
	
	TH2D *correlation_tpc_AB;
	
	TH2D *correlation_abew;
	TH2D *correlation_cew;
	TH2D *correlation_dew;
	TH2D *correlation_ew;	
	TH2D *correlation_ew_cent[9];	
	
	//flattening coefficients:
	Int_t VZ_BINS = CONST_VZ_BINS;
	
	//for recentering:
	TProfile2D *tpc_recenterx[CONST_VZ_BINS][2];//0=+vz, 1=-vz,A is east, negative, b is west, positive
	TProfile2D *tpc_recentery[CONST_VZ_BINS][2];//0=+vz, 1=-vz,A is east, negative, b is west, positive
	TProfile2D *epd_recenterx[CONST_VZ_BINS][5][2];//vz bins, then ab, c,d,abcd,cd in that order, then 0 is west, positive, 1 is east, negative
	TProfile2D *epd_recentery[CONST_VZ_BINS][5][2];//0=+vz, 1=-vz, then ab, c,d,abcd,cd in that order
	
	TProfile2D * gettpc_recenterx[CONST_VZ_BINS][2];
	TProfile2D * gettpc_recentery[CONST_VZ_BINS][2];
	TProfile2D * getepd_recenterx[CONST_VZ_BINS][5][2];
	TProfile2D * getepd_recentery[CONST_VZ_BINS][5][2];
	
	const int max_flat = 20;
	TProfile2D * epd_flatten_cos[CONST_VZ_BINS][5][20][2];
	TProfile2D * epd_flatten_sin[CONST_VZ_BINS][5][20][2];
	TProfile2D * epd_flatten_ab_ew_cos[CONST_VZ_BINS][20];
	TProfile2D * epd_flatten_ab_ew_sin[CONST_VZ_BINS][20];
	TProfile2D * epd_flatten_ew_cos[CONST_VZ_BINS][20];
	TProfile2D * epd_flatten_ew_sin[CONST_VZ_BINS][20];
	
	TProfile2D * tpc_flatten_cos[CONST_VZ_BINS][3][20]; //1 is A 0 is B,  i know...
	TProfile2D * tpc_flatten_sin[CONST_VZ_BINS][3][20]; //1 is A 0 is B,  i know...	
	
	TProfile2D * getepd_flatten_cos[CONST_VZ_BINS][5][20][2];//then ab, c,d,abcd in that order
	TProfile2D * getepd_flatten_sin[CONST_VZ_BINS][5][20][2];//then ab, c,d,abcd in that order	
	TProfile2D * getepd_flatten_ab_ew_cos[CONST_VZ_BINS][20];
	TProfile2D * getepd_flatten_ab_ew_sin[CONST_VZ_BINS][20];	
	TProfile2D * getepd_flatten_ew_cos[CONST_VZ_BINS][20];
	TProfile2D * getepd_flatten_ew_sin[CONST_VZ_BINS][20];

	TProfile2D * gettpc_flatten_cos[CONST_VZ_BINS][3][20];
	TProfile2D * gettpc_flatten_sin[CONST_VZ_BINS][3][20];	
	//TH1D * events_per_run;
	
	
	// TProfile *epd_flatten_sin_cent[4][20];
	// TProfile *epd_flatten_cos_cent[4][20];
	// TProfile *epd_flatten_sin_none[4][20];
	// TProfile *epd_flatten_cos_none[4][20];	
	// TProfile3D * pp_EPDshiftpar_sin[7];
	// TProfile3D * pp_EPDshiftpar_cos[7];
	
	//shaowei's refine
	// TH2D * h2_pt_y_mc;
	
	//QA tree things
	// TTree *QA_tree;
	Int_t run_number;
	Int_t event_number;
	Double_t refmult;
	Float_t vr;
	Float_t vz;
	Float_t vx;
	Float_t vy;
	// Float_t dedx_protonplus_avg;
	// Float_t dedx_pionplus_avg;
	// Float_t dedx_protonminus_avg;
	// Float_t dedx_pionminus_avg;	
	// Float_t mult_protonplus;
	// Float_t mult_pionplus;
	// Float_t mult_protonminus;
	// Float_t mult_pionminus;	
	
	
	// charged particle v1_eta:
	TProfile *ch_v1_eta_cent[9];
	
	//proton  profile things
	TProfile *p_protonplusv1_y_cent[CONST_VZ_BINS][9][9];
	TProfile *p_protonminusv1_y_cent[CONST_VZ_BINS][9][9];
	TH1D *h_protonminus_y_cent[9][1][10][10][6];
	
	TProfile *p_protonplusv1_yn_cent[9];
	TProfile *p_protonminusv1_yn_cent[9];	
	
	TProfile *p_protonplusv1_eta_cent[9];
	TProfile *p_protonminusv1_eta_cent[9];
	
	TProfile *p_protonplusv1_etan_cent[9];
	TProfile *p_protonminusv1_etan_cent[9];
	
	TProfile *p_protonplusv1_pt_cent[9];
	TProfile *p_protonminusv1_pt_cent[9];
	
	TProfile *p_protonplusv1_y_pip[10];
	TProfile *p_protonminusv1_y_pip[10];
	TProfile *p_protonplusv1_y_mid_pip[10];
	TProfile *p_protonminusv1_y_mid_pip[10];
	
	TProfile *p_protonplusv1_y_mid_vz[4];
	TProfile *p_protonminusv1_y_mid_vz[4];
	
	TH1D* proton_yield[9];
	TH1D* antiproton_yield[9];
	//pion stuff:
	TProfile *p_pionplusv1_y_cent[9];
	TProfile *p_pionminusv1_y_cent[9];
	TH1D* pip_yield[9];
	TH1D* pim_yield[9];
	TH2D* pi_corr_yield;
	TH1D* pip_b_tot_pi_ratio;
	TH2D* pip_b_tot_pi_cent;
	TProfile* pip_b_tot_pi_cent_p;
	TH1D* hist_pion_yield_ratio[9];
	TProfile* profile_pion_yield_ratio;
	
	//kaon stuff:
	TProfile *p_kaonplusv1_y_cent[9];
	TProfile *p_kaonminusv1_y_cent[9];
	//TH1D *h_protonminus_pt_cent[9][10][10][6];
	
	
	// TProfile *p_protonplusv1_y_010;
	// TProfile *p_protonplusv1_y_1040;
	// TProfile *p_protonplusv1_y_4060;
	// TProfile *p_protonplusv1_y_0060;
	
	// TProfile *p_protonpluspx_y_010;
	// TProfile *p_protonpluspx_y_1040;
	// TProfile *p_protonpluspx_y_4060;
	// TProfile *p_protonpluspx_y_0060;
	

	//centrality:
	TH1D *hist_cent;
	// TH1D *hist_cent_weighted;
	//dedx:
	TH2D *hist_dEdx;
	TH2D *hist_dEdx_proton;
	TH2D *hist_dEdx_pion;
	TH2D *hist_dEdx_kaon;
	
	TH2D *hist_dca_vs_p_proton;
	TH2D *hist_dca_vs_p_aproton;
	
	TH2D *hist_pt_centrality;
	TH2D *hist_pt_eta;
	TProfile *hist_Qr_epd[4][2];

	TH2D* hist_tof_pid;
	TH2D* hist_mass2_pid;
	
	TH2D* tof_efficiency_tot;
	TH2D* tof_efficiency_match;
	
	TH2D* hist_tof_pid_proton;
	TH2D *hist_tof_pid_pion;
	TH2D *hist_tof_pid_kaon;
	
	TH2D* hist_mass2_pid_proton;
	TH2D* hist_mass2_pid_pion;
	TH2D* hist_mass2_pid_kaon;
	
	TH2D* hist_nsigma_p_proton;
	TH2D* hist_nsigma_p_aproton;
	TProfile* p_nsigma_proton;
	TProfile* p_nsigma_aproton;
	
	TH2D* hist_nsigma_p_pip;
	TH2D* hist_nsigma_p_pim;
	TProfile* p_nsigma_pip;
	TProfile* p_nsigma_pim;
	
	//purity:
	TH1D* hist_nsigma_proton;
	TH1D* hist_nsigma_aproton;
	TH2D* hist_mass2_proton;
	TH2D* hist_mass2_aproton;
	
	TH2D* hist_proton_accuracy;
	TH2D* hist_aproton_accuracy;
	
	//proton acceptance:
	TH2D *hist_pt_y_proton;
	TH2D *hist_pt_y_antiproton;
	TH2D *hist_pt_y_proton_cent[9];
	TH2D *hist_pt_y_antiproton_cent[9];
	TH2D *hist_pt_y_pip;
	TH2D *hist_pt_y_pim;
	TH2D *hist_pt_y_kp;
	TH2D *hist_pt_y_km;
	
	//proton efficiency:
	TH2D *hist_pt_y_mc_proton_cent[9];
	TH2D *hist_pt_y_match_proton_cent[9];
	TH2D *hist_pt_y_mc_aproton_cent[9];
	TH2D *hist_pt_y_match_aproton_cent[9];
	TH2D *hist_eta_y_tpc_hit;
	TH2D *hist_eta_y_tof_match_hit;
	
	TH2D *hist_pt_y_proton_vz[4];
	TH2D *hist_pt_y_antiproton_vz[4];
	
	TProfile *eta_v_pt_pos_cent[9];
	TProfile *eta_v_pt_neg_cent[9];	
	
	//keeping these because they're definitely useful
	StPicoDstMaker *mPicoDstMaker;
	StPicoDst   *mPicoDst;
	StPicoEvent *picoEvent;
	StEpdGeom *mEpdGeom;
	StEpdEpFinder* mEpFinder; 
	StEpdEpInfo *mEpdEpInfo;
	StPileupUtil* mPileupTool;
	StRefMultCorr *mRefMultCorrUtil;

	TFile *File;
	TFile *re_file;
	TFile *shift_file;
	TFile *eta_weight;
	TFile *efficiency_proton_file;
	TFile *efficiency_aproton_file;
	TFile *efficiency_tof_file;
	
	TString mout_shift;

	ClassDef(Shift,0);
};
#endif
