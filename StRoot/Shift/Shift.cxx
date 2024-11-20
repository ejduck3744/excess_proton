#include <TFile.h>
#include <TTree.h>
#include <StMessMgr.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TROOT.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <StThreeVectorF.hh>
#include <StHelix.hh>
#include <TLorentzVector.h>

#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "StPicoEvent/StPicoETofPidTraits.h"
#include "StPicoEvent/StPicoEpdHit.h"
#include "StBTofUtil/tofPathLength.hh"
#include "StarClassLibrary/StPhysicalHelixD.hh"
#include "StarClassLibrary/StLorentzVectorF.hh"
#include "StEpdUtil/StEpdGeom.h"
#include "StEpdUtil/StEpdEpFinder.h"
#include "StRoot/StPicoEvent/StPicoArrays.h"
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"
#include "StPileupUtil/StPileupUtil.h"
#include "phys_constants.h"

#include "../run/run.h"
#include "../run/badrun.h"
#include "Shift.h"

ClassImp(Shift)

    //__________________________________________________________________________________
    Shift::Shift( const Char_t *name, StPicoDstMaker *picoMaker, const Char_t *jobid ) : StMaker(name) {
        mPicoDstMaker = picoMaker;
        mPicoDst = 0;

        mout_shift=Form("%s.root", jobid);
    }

//__________________________________________________________________________________
Int_t Shift::Init() {
    cout << "Init" << endl;
    const float pi = acos(-1.0);                     
    mEpdGeom = new StEpdGeom();
	mRefMultCorrUtil = new StRefMultCorr("refmult");
    mPileupTool = new StPileupUtil();
    mPileupTool->init();
    cout << "Define Histograms" << endl;
    //===============================                          
    //  Define Histograms                           
    //===============================                                            

    File=new TFile(mout_shift.Data(),"RECREATE");
	cout<<"dataset == "<<dataset<<endl;
	r1 = new TRandom3(); r1->GetSeed();
	y_beam = TMath::ACosH(beam_energy/protonMass);
	TString p_name;
	TString p_title;
	
	//because I don't know how to use phys_constants
	protonMass = 0.93827208816;
	pionMass = 0.1395;
	kaonMass = 0.4936;
	
	//primary vertex
	hist_VyVx = new TH2D("hist_VyVx","V_{Y} vs. V_{X} (pre-cut) ;V_{X} (cm);V_{Y} (cm)",500,-5.0,5.0,500,-5.0,5.0);
	
	if(dataset == "3p85_Gev_2018"){hist_Vz = new TH1D("hist_vz","V_{Z} (post-cut);V_{Z} (cm)",200,198,202);}
	else{hist_Vz = new TH1D("hist_vz","V_{Z} (post-cut);V_{Z} (cm)",1500,-150,150);}

	nhits_vs_eta= new TH2D("nhits_vs_eta","nhits_vs_eta ;#eta;nhits",100,-4,4,100,-0.5,99.5);
	
	if(dataset == "27_Gev_2018")
	{
		if(ep_setting == "recenter"){mEpFinder = new StEpdEpFinder(9*14,"StEpdEpFinderCorrectionHistograms_OUTPUT.root");} // 9 centrality bins (first 0-8 for inner 9-17 for outer)}
		else{mEpFinder = new StEpdEpFinder(9*14,"StEpdEpFinderCorrectionHistograms_OUTPUT.root","/star/u/educkwort/excess_proton/all_energies/flow/EP_finder_27_gev_2024-04-28.root");}
		mEpFinder->SetnMipThreshold(0.3);
		mEpFinder->SetMaxTileWeight(2.0);
		mEpFinder->SetEpdHitFormat(2); 
	}
	//Yrecenter
	if(ep_setting != "recenter" && dataset != "27_Gev_2018")
	{
		re_file = TFile::Open(re_file_name.Data());
		// if(dataset == "9p2_GeV_2020"){re_file = TFile::Open("/star/u/educkwort/excess_proton/all_energies/flatten/recenter_9p2_gev_2024-02-25.root");}
		// else if(dataset == "3p85_Gev_2018"){re_file =TFile::Open("/star/u/educkwort/excess_proton/all_energies/flatten/recenter_3p85_Gev_2023-02-08.root");}
		// else if(dataset == "11p5_gev_2020"){re_file = TFile::Open("/star/u/educkwort/excess_proton/all_energies/flatten/recenter_11p5_gev_2024-04-24.root");}
		// else if(dataset == "14p6_gev_2019"){re_file = TFile::Open("/star/u/educkwort/excess_proton/all_energies/flatten/recenter_14p6_gev_2023-08-17.root");}
		// else if(dataset == "7p7_Gev_2021"){re_file = TFile::Open("/star/u/educkwort/excess_proton/all_energies/flatten/recenter_7p7_2023-08-27.root");}
		// else if(dataset == "17p3_gev_2021"){re_file = TFile::Open("/star/u/educkwort/excess_proton/all_energies/flatten/recenter_17p3_gev_2024-04-24.root");}
		// else if(dataset == "19p6_Gev_2019"){re_file = TFile::Open("/star/u/educkwort/excess_proton/all_energies/flatten/recenter_19p6_gev_vz_145_all_vz_binned_cd_tpc_2024-11-01.root");}
		// else if(dataset == "19p6_Gev_2019"){re_file = TFile::Open("/star/u/educkwort/excess_proton/all_energies/flatten/recenter_19p6_gev_vz_145_4_vz_bin_recenter_2024-09-12.root");}
		// else if(dataset == "19p6_Gev_2019"){re_file = TFile::Open("/star/u/educkwort/excess_proton/all_energies/flatten/recenter_19p6_GeV_vz_145_eta_recenter_vz_binned_2024-10-11.root");}
		// else if(dataset == "19p6_Gev_2019"){re_file = TFile::Open("/star/u/educkwort/excess_proton/all_energies/flatten/recenter_19p6_gev_prelim_match_2024-09-17.root");}
		// else if(dataset == "19p6_Gev_2019"){re_file = TFile::Open("/star/u/educkwort/excess_proton/all_energies/flatten/recenter_19p6_gev_vz_145_4_vz_bin_recenter_eta_wted_2024-09-26.root");}
		// else if(dataset == "19p6_Gev_2019"){re_file = TFile::Open("/star/u/educkwort/excess_proton/all_energies/flatten/recenter_19p6_gev_vz_70_2024-07-24_prelim_pre_eta.root");}
		// else if(dataset == "19p6_Gev_2019"){re_file = TFile::Open("/star/u/educkwort/excess_proton/all_energies/flatten/recenter_19p6_gev_vz_145_full_epd_eta_2024-09-05.root");}
		// else if(dataset == "19p6_Gev_2019"){re_file = TFile::Open("/star/u/educkwort/excess_proton/all_energies/flatten/recenter_19p6_gev_vz_145_2024-08-14.root");}
		// else if(dataset == "19p6_Gev_2019"){re_file = TFile::Open("/star/u/educkwort/excess_proton/all_energies/flatten/recenter_19p6_gev_145_vz_4_vz_bins_all_cd_2024-10-27.root");}
		//else if(dataset == "27_Gev_2018"){re_file = TFile::Open("/star/u/educkwort/excess_proton/all_energies/flatten/recenter_27_gev_2024-03-02.root");}
		for(Int_t i = 0; i<VZ_BINS;i++)
		{
			for(Int_t j = 0; j<5; j++)
			{
				for(Int_t k=0;k<2;k++)
				{
					p_name=Form("epd_recenterx_%d_%d_%d",i,j,k);
					re_file->GetObject(p_name.Data(),getepd_recenterx[i][j][k]);
					p_name=Form("epd_recentery_%d_%d_%d",i,j,k);
					re_file->GetObject(p_name.Data(),getepd_recentery[i][j][k]);
				}
			}
			for(Int_t j = 0; j<2; j++)
			{
				p_name=Form("tpc_recenterx_%d_%d",i,j);
				re_file->GetObject(p_name.Data(),gettpc_recenterx[i][j]);
				p_name=Form("tpc_recentery_%d_%d",i,j);
				re_file->GetObject(p_name.Data(),gettpc_recentery[i][j]);
			}
		}
	}
	else
	{
		for(Int_t i = 0; i<VZ_BINS;i++)
		{
			for(Int_t j = 0; j<2; j++)
			{
				p_name=Form("tpc_recenterx_%d_%d",i,j);
				tpc_recenterx[i][j] = new TProfile2D(p_name.Data(),p_name.Data(),nrun,-0.5,nrun-0.5,9,-0.5,8.5);
				p_name=Form("tpc_recentery_%d_%d",i,j);
				tpc_recentery[i][j] = new TProfile2D(p_name.Data(),p_name.Data(),nrun,-0.5,nrun-0.5,9,-0.5,8.5);
			}
			for(Int_t j = 0; j<5; j++)
			{
				for(Int_t k=0;k<2;k++)
				{
					p_name=Form("epd_recenterx_%d_%d_%d",i,j,k);
					epd_recenterx[i][j][k]= new TProfile2D(p_name.Data(),p_name.Data(),nrun,-0.5,nrun-0.5,9,-0.5,8.5);
					p_name=Form("epd_recentery_%d_%d_%d",i,j,k);
					epd_recentery[i][j][k] = new TProfile2D(p_name.Data(),p_name.Data(),nrun,-0.5,nrun-0.5,9,-0.5,8.5);
				}
			}
		}
	}
	
	if(ep_setting =="flow" && dataset != "27_Gev_2018")
	{
		shift_file = TFile::Open(shift_file_name.Data());
		// if(dataset == "9p2_GeV_2020"){shift_file = TFile::Open("/star/u/educkwort/excess_proton/all_energies/flow/flatten_9p2_gev_2024-02-29.root");}
		// else if(dataset == "7p7_Gev_2021"){shift_file = TFile::Open("/star/u/educkwort/excess_proton/all_energies/flow/flatten_7p7_2023-08-28.root");} 
		// else if(dataset == "3p85_Gev_2018"){shift_file = TFile::Open("/star/u/educkwort/excess_proton/all_energies/flow/flatten_3p85_Gev_2023-02-09.root");}
		// else if(dataset == "11p5_gev_2020"){shift_file = TFile::Open("/star/u/educkwort/excess_proton/all_energies/flow/flatten_11p5_gev_2024-04-25.root");}
		// else if(dataset == "14p6_gev_2019"){shift_file = TFile::Open("/star/u/educkwort/excess_proton/all_energies/flow/flatten_14p6_gev_2023-08-18.root");}
		// else if(dataset == "17p3_gev_2021"){shift_file = TFile::Open("/star/u/educkwort/excess_proton/all_energies/flow/flatten_17p3_gev_2024-04-25.root");}
		// else if(dataset == "19p6_Gev_2019"){shift_file = TFile::Open("/star/u/educkwort/excess_proton/all_energies/flow/flatten_19p6_gev_vz_70_2024-07-25_prelim_pre_eta.root");}
		// else if(dataset == "19p6_Gev_2019"){shift_file = TFile::Open("/star/u/educkwort/excess_proton/all_energies/flow/flatten_19p6_gev_vz_145_4_vz_bin_recenter_2024-09-17.root");}
		// else if(dataset == "19p6_Gev_2019"){shift_file = TFile::Open("/star/u/educkwort/excess_proton/all_energies/flow/flatten_19p6_gev_vz_145_eta_recenter_vz_binned_2024-10-13.root");}
		// else if(dataset == "19p6_Gev_2019"){shift_file = TFile::Open("/star/u/educkwort/excess_proton/all_energies/flow/flatten_19p6_gev_145_vz_4_vz_bins_cd_tpc_2024-11-06.root");}
		// else if(dataset == "19p6_Gev_2019"){shift_file = TFile::Open("/star/u/educkwort/excess_proton/all_energies/flow/flatten_19p6_gev_prelim_match_2024-09-19.root");}
		// else if(dataset == "19p6_Gev_2019"){shift_file = TFile::Open("/star/u/educkwort/excess_proton/all_energies/flow/flatten_19p6_gev_vz_145_4_vz_bin_recenter_eta_wt_2024-09-29.root");}
		// else if(dataset == "19p6_Gev_2019"){shift_file = TFile::Open("/star/u/educkwort/excess_proton/all_energies/flow/flatten_19p6_gev_vz_145_full_epd_eta_wt_2024-09-06.root");}
		// else if(dataset == "19p6_Gev_2019"){shift_file = TFile::Open("/star/u/educkwort/excess_proton/all_energies/flow/flatten_19p6_gev_145_vz_4_vz_bins_all_cd_2024-10-27.root");}
		for(Int_t h=0;h<VZ_BINS;h++)
		{
			for(Int_t j = 0; j<max_flat; j++)
			{
				for(Int_t i = 0; i<5;i++)
				{
					for(Int_t k = 0; k<2; k++)
					{
						p_name=Form("epd_flatten_cos_%d_%d_%d_%d",h,i,j,k);
						shift_file->GetObject(p_name.Data(),getepd_flatten_cos[h][i][j][k]);
						p_name=Form("epd_flatten_sin_%d_%d_%d_%d",h,i,j,k);
						shift_file->GetObject(p_name.Data(),getepd_flatten_sin[h][i][j][k]);
					}
				}
				p_name=Form("epd_flatten_ab_ew_%d_cos_%d",h,j);
				shift_file->GetObject(p_name.Data(),getepd_flatten_ab_ew_cos[h][j]);
				p_name=Form("epd_flatten_ab_ew_%d_sin_%d",h,j);
				shift_file->GetObject(p_name.Data(),getepd_flatten_ab_ew_sin[h][j]);
				
				p_name=Form("epd_flatten_ew_%d_cos_%d",h,j);
				shift_file->GetObject(p_name.Data(),getepd_flatten_ew_cos[h][j]);
				p_name=Form("epd_flatten_ew_%d_sin_%d",h,j);
				shift_file->GetObject(p_name.Data(),getepd_flatten_ew_sin[h][j]);
			}
			for(Int_t i = 0; i<3;i++)
			{
				for(Int_t j = 0; j<max_flat; j++)
				{
					p_name=Form("tpc_flatten_cos_%d_%d_%d",h,i,j);
					shift_file->GetObject(p_name.Data(),gettpc_flatten_cos[h][i][j]);
					p_name=Form("tpc_flatten_sin_%d_%d_%d",h,i,j);
					shift_file->GetObject(p_name.Data(),gettpc_flatten_sin[h][i][j]);
				}
			}
		}
	}
	else
	{
		for(Int_t h=0;h<VZ_BINS;h++)
		{
			for(Int_t j=0;j<max_flat;j++)
			{
				p_name=Form("epd_flatten_ab_ew_%d_cos_%d",h,j);
				epd_flatten_ab_ew_cos[h][j] = new TProfile2D(p_name.Data(),p_name.Data(),nrun,-0.5,nrun-0.5,9,-0.5,8.5);
				p_name=Form("epd_flatten_ab_ew_%d_sin_%d",h,j);
				epd_flatten_ab_ew_sin[h][j] = new TProfile2D(p_name.Data(),p_name.Data(),nrun,-0.5,nrun-0.5,9,-0.5,8.5);
				p_name=Form("epd_flatten_ew_%d_cos_%d",h,j);
				epd_flatten_ew_cos[h][j] = new TProfile2D(p_name.Data(),p_name.Data(),nrun,-0.5,nrun-0.5,9,-0.5,8.5);
				p_name=Form("epd_flatten_ew_%d_sin_%d",h,j);
				epd_flatten_ew_sin[h][j] = new TProfile2D(p_name.Data(),p_name.Data(),nrun,-0.5,nrun-0.5,9,-0.5,8.5);
				for(Int_t i=0;i<5;i++)
				{
					for(Int_t k=0;k<2;k++)
					{
						p_name=Form("epd_flatten_cos_%d_%d_%d_%d",h,i,j,k);
						epd_flatten_cos[h][i][j][k] = new TProfile2D(p_name.Data(),p_name.Data(),nrun,-0.5,nrun-0.5,9,-0.5,8.5);
						p_name=Form("epd_flatten_sin_%d_%d_%d_%d",h,i,j,k);
						epd_flatten_sin[h][i][j][k] = new TProfile2D(p_name.Data(),p_name.Data(),nrun,-0.5,nrun-0.5,9,-0.5,8.5);
					}
				}
				for(Int_t i=0;i<3;i++)
				{
					p_name=Form("tpc_flatten_cos_%d_%d_%d",h,i,j);
					tpc_flatten_cos[h][i][j] = new TProfile2D(p_name.Data(),p_name.Data(),nrun,-0.5,nrun-0.5,9,-0.5,8.5);
					p_name=Form("tpc_flatten_sin_%d_%d_%d",h,i,j);
					tpc_flatten_sin[h][i][j] = new TProfile2D(p_name.Data(),p_name.Data(),nrun,-0.5,nrun-0.5,9,-0.5,8.5);
				}
			}
		}
	}

	if(do_eta_weighting&&(dataset == "9p2_GeV_2020"||dataset == "14p6_gev_2019"||dataset == "17p3_gev_2021" || dataset == "11p5_gev_2020"))
	{
		eta_weight = TFile::Open(eta_weight_name.Data());
		// if(dataset == "9p2_GeV_2020" ){eta_weight = TFile::Open("/star/u/educkwort/excess_proton/all_energies/recenter/eta_9p2_2023-04-03.root" );}
		// if(dataset == "14p6_gev_2019"){eta_weight = TFile::Open("/star/u/educkwort/excess_proton/all_energies/recenter/eta_14p5_2022-07-12.root");}
		// if(dataset == "17p3_gev_2021"){eta_weight = TFile::Open("/star/u/educkwort/excess_proton/all_energies/recenter/eta_17p3_2023-10-10.root");}
		// if(dataset == "11p5_gev_2020"){eta_weight = TFile::Open("/star/u/educkwort/excess_proton/all_energies/recenter/eta_11p5_2023-08-04.root");}
		eta_weight->GetObject("v1_vs_eta_east",epd_eta_weight_east);
		eta_weight->GetObject("v1_vs_eta_west",epd_eta_weight_west);
	}
	if(do_eta_weighting && dataset == "19p6_Gev_2019")
	{
		// if(dataset == "19p6_Gev_2019"){eta_weight = TFile::Open("/star/u/educkwort/excess_proton/all_energies/recenter/eta_19p6_gev_vz_70_2024-07-26.root");}
		// if(dataset == "19p6_Gev_2019"){eta_weight = TFile::Open("/star/u/educkwort/excess_proton/all_energies/recenter/eta_19p6_gev_vz_145_full_epd_eta_wt_2024-09-04.root");}
		// eta_weight = TFile::Open("/star/u/educkwort/excess_proton/all_energies/recenter/eta_19p6_gev_prelim_match_2024-09-11.root");
		// eta_weight = TFile::Open("/star/u/educkwort/excess_proton/all_energies/recenter/eta_19p6_gev_vz_145_4_vz_bin_recenter_2024-09-20.root");
		eta_weight = TFile::Open(eta_weight_name.Data());
		// p_name = Form("v1_vs_eta_east_%d",0);
		// eta_weight->GetObject(p_name.Data(),epd_eta_weight_east);
		// p_name = Form("v1_vs_eta_west_%d",0);
		// eta_weight->GetObject(p_name.Data(),epd_eta_weight_west);
		// for(int i = 1;i<9;i++)
		// {
			// p_name = Form("v1_vs_eta_east_%d",i);
			// eta_weight->GetObject(p_name.Data(),epd_eta_weight_east_cent[i]);
			// p_name = Form("v1_vs_eta_west_%d",i);
			// eta_weight->GetObject(p_name.Data(),epd_eta_weight_west_cent[i]);
			// epd_eta_weight_east->Add(epd_eta_weight_east_cent[i],1.0);
			// epd_eta_weight_west->Add(epd_eta_weight_west_cent[i],1.0);
		// }
		for(int i = 0;i<9;i++)
		{
			for(int j = 0;j<VZ_BINS;j++)
			{
				p_name  = Form("v1_vs_eta_east_%d_%d",i,j);
				eta_weight->GetObject(p_name.Data(),epd_eta_weight_east_cent[i][j]);
				p_name  = Form("epd_eta_weight_east_cent_%d_%d",i,j);
				epd_eta_weight_east_cent[i][j]->SetName(p_name.Data());
				p_name  = Form("v1_vs_eta_west_%d_%d",i,j);
				eta_weight->GetObject(p_name.Data(),epd_eta_weight_west_cent[i][j]);
				p_name  = Form("epd_eta_weight_west_cent_%d_%d",i,j);
				epd_eta_weight_west_cent[i][j]->SetName(p_name.Data());
			}
			
		}
	}
	if(do_eta_weighting && dataset == "27_Gev_2018")
	{
		eta_weight = TFile::Open("/star/u/educkwort/excess_proton/all_energies/recenter/eta_27_gev_2024-04-22.root");
		TH2D wt_n1("Order1etaWeight","Order1etaWeight",60,2.1,5.1,9*14,0,9*14);
		for(int i = 0;i<9;i++)
		{
			for(int j = 0;j<VZ_BINS;j++)
			{
				p_name  = Form("v1_vs_eta_east_%d_%d",i,j);
				eta_weight->GetObject(p_name.Data(),epd_eta_weight_east_cent[i][j]);
				p_name  = Form("v1_vs_eta_west_%d_%d",i,j);
				eta_weight->GetObject(p_name.Data(),epd_eta_weight_west_cent[i][j]);
			}
		}
		for(int i = 0; i<9*14;i++)
		{
			int cent = i%9;
			for(int h = 0;h<VZ_BINS;h++)
			{
				for(int j = 0; j<60;j++)
				{
					double eta = wt_n1.GetXaxis()->GetBinCenter(j);
					double weight = 0.5*(epd_eta_weight_west_cent[cent][h]->GetBinContent(epd_eta_weight_west_cent[cent][h]->FindBin(eta))+
									epd_eta_weight_east_cent[cent][h]->GetBinContent(epd_eta_weight_east_cent[cent][h]->FindBin(-1.0*eta)));
					wt_n1.Fill(eta,i,weight);
				}
			}
		}
		mEpFinder->SetEtaWeights(1,wt_n1);
	}
	
	//EPD stuff:
	hist_epd_ab_psi_raw_ew        = new TH1D("hist_epd_ab_psi_raw_ew"       ,"hist_epd_ab_psi_raw_ew ;#psi^{epd ab_ew}_{1} [Radian];# of events"       ,500,0.0,2.0*pi);
	hist_epd_ab_psi_recentered_ew = new TH1D("hist_epd_ab_psi_recentered_ew","hist_epd_ab_psi_recentered_ew ;#psi^{epd ab_ew}_{1} [Radian];# of events",500,0.0,2.0*pi);
	hist_epd_ab_psi_flattened_ew  = new TH1D("hist_epd_ab_psi_flattened_ew" ,"hist_epd_ab_psi_flattened_ew ;#psi^{epd ab_ew}_{1} [Radian];# of events" ,500,0.0,2.0*pi);
	hist_epd_psi_raw_ew           = new TH1D("hist_epd_psi_raw_ew"          ,"hist_epd_psi_raw_ew ;#psi^{epd ab_ew}_{1} [Radian];# of events"          ,500,0.0,2.0*pi);
	hist_epd_psi_recentered_ew    = new TH1D("hist_epd_psi_recentered_ew"   ,"hist_epd_psi_recentered_ew ;#psi^{epd ab_ew}_{1} [Radian];# of events"   ,500,0.0,2.0*pi);
	hist_epd_psi_flattened_ew     = new TH1D("hist_epd_psi_flattened_ew"    ,"hist_epd_psi_flattened_ew ; #psi^{epd ab_ew}_{1} [Radian]; # of events"  ,500,0.0,2.0*pi);
	for(int i = 0; i<2;i++)
	{
		p_name=Form("hist_epd_ab_psi_raw_%d",i);
		hist_epd_ab_psi_raw[i] = new TH1D(p_name.Data(),p_name.Data(),500,0.0,2.0*pi);
		hist_epd_ab_psi_raw[i]->GetXaxis()->SetTitle("#psi^{epd ab}_{1} [Radian]");
		hist_epd_ab_psi_raw[i]->GetYaxis()->SetTitle("# of events");
		
		p_name=Form("hist_epd_c_psi_raw_%d",i);
		hist_epd_c_psi_raw[i] = new TH1D(p_name.Data(),p_name.Data(),500,0.0,2.0*pi);
		hist_epd_c_psi_raw[i]->GetXaxis()->SetTitle("#psi^{epd c}_{1} [Radian]");
		hist_epd_c_psi_raw[i]->GetYaxis()->SetTitle("# of events");
		
		p_name=Form("hist_epd_d_psi_raw_%d",i);
		hist_epd_d_psi_raw[i] = new TH1D(p_name.Data(),p_name.Data(),500,0.0,2.0*pi);
		hist_epd_d_psi_raw[i]->GetXaxis()->SetTitle("#psi^{epd d}_{1} [Radian]");
		hist_epd_d_psi_raw[i]->GetYaxis()->SetTitle("# of events");
		
		p_name=Form("hist_epd_abcd_psi_raw_%d",i);
		hist_epd_abcd_psi_raw[i] = new TH1D(p_name.Data(),p_name.Data(),500,0.0,2.0*pi);
		hist_epd_abcd_psi_raw[i]->GetXaxis()->SetTitle("#psi^{epd abcd}_{1} [Radian]");
		hist_epd_abcd_psi_raw[i]->GetYaxis()->SetTitle("# of events");
		
		p_name=Form("hist_epd_cd_psi_raw_%d",i);
		hist_epd_cd_psi_raw[i] = new TH1D(p_name.Data(),p_name.Data(),500,0.0,2.0*pi);
		hist_epd_cd_psi_raw[i]->GetXaxis()->SetTitle("#psi^{epd cd}_{1} [Radian]");
		hist_epd_cd_psi_raw[i]->GetYaxis()->SetTitle("# of events");
		
		p_name=Form("hist_epd_ab_psi_recentered_%d",i);
		hist_epd_ab_psi_recentered[i] = new TH1D(p_name.Data(),p_name.Data(),500,0.0,2.0*pi);
		hist_epd_ab_psi_recentered[i]->GetXaxis()->SetTitle("#psi^{epd ab}_{1} [Radian]");
		hist_epd_ab_psi_recentered[i]->GetYaxis()->SetTitle("# of events");
		
		p_name=Form("hist_epd_c_psi_recentered_%d",i);
		hist_epd_c_psi_recentered[i] = new TH1D(p_name.Data(),p_name.Data(),500,0.0,2.0*pi);
		hist_epd_c_psi_recentered[i]->GetXaxis()->SetTitle("#psi^{epd c}_{1} [Radian]");
		hist_epd_c_psi_recentered[i]->GetYaxis()->SetTitle("# of events");
		
		p_name=Form("hist_epd_d_psi_recentered_%d",i);
		hist_epd_d_psi_recentered[i] = new TH1D(p_name.Data(),p_name.Data(),500,0.0,2.0*pi);
		hist_epd_d_psi_recentered[i]->GetXaxis()->SetTitle("#psi^{epd d}_{1} [Radian]");
		hist_epd_d_psi_recentered[i]->GetYaxis()->SetTitle("# of events");
		
		p_name=Form("hist_epd_abcd_psi_recentered_%d",i);
		hist_epd_abcd_psi_recentered[i] = new TH1D(p_name.Data(),p_name.Data(),500,0.0,2.0*pi);
		hist_epd_abcd_psi_recentered[i]->GetXaxis()->SetTitle("#psi^{epd abcd}_{1} [Radian]");
		hist_epd_abcd_psi_recentered[i]->GetYaxis()->SetTitle("# of events");
		
		p_name=Form("hist_epd_cd_psi_recentered_%d",i);
		hist_epd_cd_psi_recentered[i] = new TH1D(p_name.Data(),p_name.Data(),500,0.0,2.0*pi);
		hist_epd_cd_psi_recentered[i]->GetXaxis()->SetTitle("#psi^{epd cd}_{1} [Radian]");
		hist_epd_cd_psi_recentered[i]->GetYaxis()->SetTitle("# of events");
		
		p_name=Form("hist_epd_ab_psi_flattened_%d",i);
		hist_epd_ab_psi_flattened[i] = new TH1D(p_name.Data(),p_name.Data(),500,0.0,2.0*pi);
		hist_epd_ab_psi_flattened[i]->GetXaxis()->SetTitle("#psi^{epd ab}_{1} [Radian]");
		hist_epd_ab_psi_flattened[i]->GetYaxis()->SetTitle("# of events");
		
		p_name=Form("hist_epd_c_psi_flattened_%d",i);
		hist_epd_c_psi_flattened[i] = new TH1D(p_name.Data(),p_name.Data(),500,0.0,2.0*pi);
		hist_epd_c_psi_flattened[i]->GetXaxis()->SetTitle("#psi^{epd c}_{1} [Radian]");
		hist_epd_c_psi_flattened[i]->GetYaxis()->SetTitle("# of events");
		
		p_name=Form("hist_epd_d_psi_flattened_%d",i);
		hist_epd_d_psi_flattened[i] = new TH1D(p_name.Data(),p_name.Data(),500,0.0,2.0*pi);
		hist_epd_d_psi_flattened[i]->GetXaxis()->SetTitle("#psi^{epd d}_{1} [Radian]");
		hist_epd_d_psi_flattened[i]->GetYaxis()->SetTitle("# of events");
		
		p_name=Form("hist_epd_abcd_psi_flattened_%d",i);
		hist_epd_abcd_psi_flattened[i] = new TH1D(p_name.Data(),p_name.Data(),500,0.0,2.0*pi);
		hist_epd_abcd_psi_flattened[i]->GetXaxis()->SetTitle("#psi^{epd abcd}_{1} (Radian)");
		hist_epd_abcd_psi_flattened[i]->GetYaxis()->SetTitle("# of events");
		
		p_name=Form("hist_epd_cd_psi_flattened_%d",i);
		hist_epd_cd_psi_flattened[i] = new TH1D(p_name.Data(),p_name.Data(),500,0.0,2.0*pi);
		hist_epd_cd_psi_flattened[i]->GetXaxis()->SetTitle("#psi^{epd cd}_{1} (Radian)");
		hist_epd_cd_psi_flattened[i]->GetYaxis()->SetTitle("# of events");
		for(Int_t j=0;j<4;j++)
		{
			p_name=Form("hist_Qr_epd_%d_%d",j,i);
			hist_Qr_epd[j][i] = new TProfile(p_name.Data(),p_name.Data(),9,-0.5,8.5);
			hist_Qr_epd[j][i]->GetXaxis()->SetTitle("centrality");
			hist_Qr_epd[j][i]->GetYaxis()->SetTitle("|Q|");
		}
	}
	//eta weighting:
	for(int i = 0; i<9;i++)
	{
		for(int j = 0;j<VZ_BINS;j++)
		{
			//epd eta range is 2.1 to 5.1, but thats only if vz = 0, if vz can be as big as 145 and vr can be as big as 2, then
			// the max eta is ~6.1 and the min eta is ~1.6
			p_name  = Form("v1_vs_eta_%d_%d",i,j);
			p_title = Form("v_{1} vs #eta TPC centrality %d %d;#eta;v_{1}",i,j);
			v1_vs_eta_tpc[i][j] = new TProfile(p_name.Data(),p_title.Data(),40,-1.0,1.0);
			p_name  = Form("v1_vs_eta_east_%d_%d",i,j);
			p_title = Form("v_{1} vs #eta east centrality %d %d;#eta;v_{1}",i,j);
			v1_vs_eta_east[i][j] = new TProfile(p_name.Data(),p_title.Data(),90,-6.1,-1.6);//east is negative
			p_name  = Form("v1_vs_eta_west_%d_%d",i,j);
			p_title = Form("v_{1} vs #eta west centrality %d %d;#eta;v_{1}",i,j);
			v1_vs_eta_west[i][j] = new TProfile(p_name.Data(),p_title.Data(),90, 1.6, 6.1);  //west is positive
		}
	}
		
	//tpc stuff:
	hist_tpc_A_psi_raw         = new TH1D("hist_tpc_A_psi_raw"       ,"TPC A EP raw ; #psi^{tpc A}_{1} [Radian] ; # of events"       ,500,0.0,2.0*pi); //A is negative pseudorapidity
	hist_tpc_B_psi_raw         = new TH1D("hist_tpc_B_psi_raw"       ,"TPC B EP raw ; #psi^{tpc B}_{1} [Radian]; # of events"        ,500,0.0,2.0*pi); //B is positive pseudorapidity
	hist_tpc_A_psi_recentered  = new TH1D("hist_tpc_A_psi_recentered","TPC A EP recentered ; #psi^{tpc A}_{1} [Radian] ; # of events",500,0.0,2.0*pi); //A is negative pseudorapidity
	hist_tpc_B_psi_recentered  = new TH1D("hist_tpc_B_psi_recentered","TPC B EP recentered ;#psi^{tpc B}_{1} [Radian];# of events"   ,500,0.0,2.0*pi); //B is positive pseudorapidity
	hist_tpc_AB_psi_recentered = new TH1D("hist_tpc_AB_psi_recentered","TPC EW EP recentered ;#psi^{tpc}_{1} [Radian];# of events"   ,500,0.0,2.0*pi); 
	hist_tpc_A_psi_flattened   = new TH1D("hist_tpc_A_psi_flattened" ,"TPC A EP flattened ; #psi^{tpc A}_{1} [Radian] ; # of events" ,500,0.0,2.0*pi); //A is negative pseudorapidity
	hist_tpc_B_psi_flattened   = new TH1D("hist_tpc_B_psi_flattened" ,"TPC B EP flattened ;#psi^{tpc B}_{1} [Radian] ;# of events"   ,500,0.0,2.0*pi); //B is positive pseudorapidity
	hist_tpc_AB_psi_flattened  = new TH1D("hist_tpc_AB_psi_flattened" ,"TPC EW EP flattened ;#psi^{tpc}_{1} [Radian] ;# of events"   ,500,0.0,2.0*pi);
	
	//correlation plots:
	profile_correlation_abc    = new TProfile("profile_correlation_abc"   ,"<cos(#psi^{epd_ab}_{1} #minus #psi^{epd_c})> ;Centrality;Correlation"   ,9,-0.5,8.5,0,0,"");
	profile_correlation_abd    = new TProfile("profile_correlation_abd"   ,"<cos(#psi^{epd_ab}_{1} #minus #psi^{epd_d})> ;Centrality;Correlation"   ,9,-0.5,8.5,0,0,"");
	profile_correlation_cd     = new TProfile("profile_correlation_cd"    ,"<cos(#psi^{epd_c}_{1} #minus #psi^{epd_d})> ;Centrality;Correlation"    ,9,-0.5,8.5,0,0,"");
	profile_correlation_eABtB  = new TProfile("profile_correlation_eABtB" ,"<cos(#psi^{epd_AB}_{1} #minus #psi^{tpc_B})> ;Centrality;Correlation"   ,9,-0.5,8.5,0,0,"");
	profile_correlation_eCtB   = new TProfile("profile_correlation_eCtB"  ,"<cos(#psi^{epd_C}_{1} #minus #psi^{tpc_B})> ;Centrality;Correlation"    ,9,-0.5,8.5,0,0,"");
	profile_correlation_eDtB   = new TProfile("profile_correlation_eDtB"  ,"<cos(#psi^{epd_D}_{1} #minus #psi^{tpc_B})> ;Centrality;Correlation"    ,9,-0.5,8.5,0,0,"");
	profile_correlation_eABtA  = new TProfile("profile_correlation_eABtA" ,"<cos(#psi^{epd_AB}_{1} #minus #psi^{tpc_A})> ;Centrality;Correlation"   ,9,-0.5,8.5,0,0,"");
	profile_correlation_eCtA   = new TProfile("profile_correlation_eCtA"  ,"<cos(#psi^{epd_C}_{1} #minus #psi^{tpc_A})> ;Centrality;Correlation"    ,9,-0.5,8.5,0,0,"");
	profile_correlation_eDtA   = new TProfile("profile_correlation_eDtA"  ,"<cos(#psi^{epd_D}_{1} #minus #psi^{tpc_A})> ;Centrality;Correlation"    ,9,-0.5,8.5,0,0,"");
	profile_correlation_tpc_AB = new TProfile("profile_correlation_tpc_AB","<cos(#psi^{tpc_A}_{1} #minus #psi^{tpc_B})> ;Centrality;Correlation"    ,9,-0.5,8.5,0,0,"");
	profile_correlation_abew   = new TProfile("profile_correlation_abew"  ,"<cos(#psi^{epd_abe}_{1} #minus #psi^{epd_abw})> ;Centrality;Correlation",9,-0.5,8.5,0,0,"");
	profile_correlation_cew    = new TProfile("profile_correlation_cew"   ,"<cos(#psi^{epd_ce}_{1} #minus #psi^{epd_cw})> ;Centrality;Correlation"  ,9,-0.5,8.5,0,0,"");
	profile_correlation_dew    = new TProfile("profile_correlation_dew"   ,"<cos(#psi^{epd_de}_{1} #minus #psi^{epd_dw})> ;Centrality;Correlation"  ,9,-0.5,8.5,0,0,"");
	profile_correlation_ew     = new TProfile("profile_correlation_ew"    ,"<cos(#psi^{epd_e}_{1} #minus #psi^{epd_w})> ;Centrality;Correlation"    ,9,-0.5,8.5,0,0,"");

	correlation_abc    = new TH2D("correlation_abc"   ,"correlation_abc ; epd_ab ; epd_c",500,0,2.0*pi,500,0,2.0*pi);
	correlation_abd    = new TH2D("correlation_abd"   ,"correlation_abd ; epd_ab ; epd_d",500,0,2.0*pi,500,0,2.0*pi);
	correlation_cd     = new TH2D("correlation_cd"    ,"correlation_cd ; epd_c ; epd_d",500,0,2.0*pi,500,0,2.0*pi);
	correlation_eABtA  = new TH2D("correlation_eABtA" ,"correlation_eABtA ;epd_ab;tpc_a",500,0,2.0*pi,500,0,2.0*pi);
	correlation_eCtA   = new TH2D("correlation_eCtA"  ,"correlation_eCtA ;epd_c;tpc_a",500,0,2.0*pi,500,0,2.0*pi);
	correlation_eDtA   = new TH2D("correlation_eDtA"  ,"correlation_eDtA ;epd_d;tpc_a",500,0,2.0*pi,500,0,2.0*pi);
	correlation_eABtB  = new TH2D("correlation_eABtB" ,"correlation_eABtB ;epd_ab;tpc_b",500,0,2.0*pi,500,0,2.0*pi);
	correlation_eCtB   = new TH2D("correlation_eCtB"  ,"correlation_eCtB ;epd_c;tpc_b",500,0,2.0*pi,500,0,2.0*pi);
	correlation_eDtB   = new TH2D("correlation_eDtB"  ,"correlation_eDtB ;epd_d;tpc_b",500,0,2.0*pi,500,0,2.0*pi);
	correlation_tpc_AB = new TH2D("correlation_tpc_AB","correlation_tpc_AB ;tpc_a;tpc_b",500,0,2.0*pi,500,0,2.0*pi);	
	correlation_abew   = new TH2D("correlation_abew"  ,"correlation_abew ; epd_abe ; epd_abw",500,0,2.0*pi,500,0,2.0*pi);
	correlation_cew    = new TH2D("correlation_cew"   ,"correlation_cew ; epd_ce ; epd_cw",500,0,2.0*pi,500,0,2.0*pi);
	correlation_dew    = new TH2D("correlation_dew"   ,"correlation_dew ; epd_de ;epd_dw",500,0,2.0*pi,500,0,2.0*pi);
	correlation_ew     = new TH2D("correlation_ew"    ,"correlation east-west EPD;#psi_{epd}^{east} (radian);#psi_{epd}^{west} (radian)",500,0,2.0*pi,500,0,2.0*pi);
	for(int i = 0;i<VZ_BINS;i++)
	{
		p_name=Form("profile_correlation_tpc_ew_vz_%d",i);
		p_title=Form("<cos(#psi^TPC_{e}} #minus #psi^{TPC_{w}})> for V_{z} %d ;Centrality;Correlation",i);
		profile_correlation_tpc_ew_vz[i] = new TProfile(p_name.Data(),p_title.Data(),9,-0.5,8.5,0,0,"");
		p_name=Form("profile_correlation_ew_vz_%d",i);
		p_title=Form("<cos(#psi^{epd_e} #minus #psi^{epd_w})> for V_{z} %d ;Centrality;Correlation",i);
		profile_correlation_ew_vz[i] = new TProfile(p_name.Data(),p_title.Data(),9,-0.5,8.5,0,0,"");
		p_name=Form("profile_correlation_eABw_vz_%d",i);
		p_title=Form("<cos(#psi^{epd_eAB} #minus #psi^{epd_w})> for V_{z} %d ;Centrality;Correlation",i);
		profile_correlation_eABw_vz[i] = new TProfile(p_name.Data(),p_title.Data(),9,-0.5,8.5,0,0,"");
		p_name=Form("profile_correlation_eCDw_vz_%d",i);
		p_title=Form("<cos(#psi^{epd_eCD} #minus #psi^{epd_w})> for V_{z} %d ;Centrality;Correlation",i);
		profile_correlation_eCDw_vz[i] = new TProfile(p_name.Data(),p_title.Data(),9,-0.5,8.5,0,0,"");
		p_name=Form("profile_correlation_eABeCD_vz_%d",i);
		p_title=Form("<cos(#psi^{epd_eCD} #minus #psi^{epd_eAB})> for V_{z} %d ;Centrality;Correlation",i);
		profile_correlation_eABeCD_vz[i] = new TProfile(p_name.Data(),p_title.Data(),9,-0.5,8.5,0,0,"");
		p_name=Form("profile_correlation_ewAB_vz_%d",i);
		p_title=Form("<cos(#psi^{epd_e} #minus #psi^{epd_wAB})> for V_{z} %d ;Centrality;Correlation",i);
		profile_correlation_ewAB_vz[i] = new TProfile(p_name.Data(),p_title.Data(),9,-0.5,8.5,0,0,"");
		p_name=Form("profile_correlation_ewCD_vz_%d",i);
		p_title=Form("<cos(#psi^{epd_e} #minus #psi^{epd_wCD})> for V_{z} %d ;Centrality;Correlation",i);
		profile_correlation_ewCD_vz[i] = new TProfile(p_name.Data(),p_title.Data(),9,-0.5,8.5,0,0,"");
		p_name=Form("profile_correlation_wABwCD_vz_%d",i);
		p_title=Form("<cos(#psi^{epd_wAB} #minus #psi^{epd_wCD})> for V_{z} %d ;Centrality;Correlation",i);
		profile_correlation_wABwCD_vz[i] = new TProfile(p_name.Data(),p_title.Data(),9,-0.5,8.5,0,0,"");
		p_name=Form("profile_correlation_eDw_vz_%d",i);
		p_title=Form("<cos(#psi^{epd_eD} #minus #psi^{epd_w})> for V_{z} %d ;Centrality;Correlation",i);
		profile_correlation_eDw_vz[i] = new TProfile(p_name.Data(),p_title.Data(),9,-0.5,8.5,0,0,"");
		p_name=Form("profile_correlation_eABeD_vz_%d",i);
		p_title=Form("<cos(#psi^{epd_eD} #minus #psi^{epd_w})> for V_{z} %d ;Centrality;Correlation",i);
		profile_correlation_eABeD_vz[i] = new TProfile(p_name.Data(),p_title.Data(),9,-0.5,8.5,0,0,"");
		p_name=Form("profile_correlation_ewD_vz_%d",i);
		p_title=Form("<cos(#psi^{epd_e} #minus #psi^{epd_wD})> for V_{z} %d ;Centrality;Correlation",i);
		profile_correlation_ewD_vz[i] = new TProfile(p_name.Data(),p_title.Data(),9,-0.5,8.5,0,0,"");
		p_name=Form("profile_correlation_wABwD_vz_%d",i);
		p_title=Form("<cos(#psi^{epd_wAB} #minus #psi^{epd_wD})> for V_{z} %d ;Centrality;Correlation",i);
		profile_correlation_wABwD_vz[i] = new TProfile(p_name.Data(),p_title.Data(),9,-0.5,8.5,0,0,"");
		p_name=Form("profile_correlation_TPCEPDe_vz_%d",i);
		p_title=Form("<cos(#psi^{tpc} #minus #psi^{epd_e})> for V_{z} %d ;Centrality;Correlation",i);
		profile_correlation_TPCEPDe_vz[i] = new TProfile(p_name.Data(),p_title.Data(),9,-0.5,8.5,0,0,"");
		p_name=Form("profile_correlation_TPCEPDw_vz_%d",i);
		p_title=Form("<cos(#psi^{tpc} #minus #psi^{epd_w})> for V_{z} %d ;Centrality;Correlation",i);
		profile_correlation_TPCEPDw_vz[i] = new TProfile(p_name.Data(),p_title.Data(),9,-0.5,8.5,0,0,"");
	}
	
	for(int i=0;i<9;i++)
	{
		p_name=Form("correlation_ew_cent_%d",i);
		p_title = Form("correlation east-west EPD for centrality %d ;#psi_{epd}^{east} (radian);#psi_{epd}^{west} (radian)",i);
		correlation_ew_cent[i]     = new TH2D(p_name.Data(),p_title.Data(),500,0,2.0*pi,500,0,2.0*pi);
	}

	//proton profile stuff:
	for(int i=0; i<9; i++)
	{
		p_name=Form("hist_pt_y_proton_cent_%d",i);
		p_title = Form("p_{T} [GeV] vs. y proton for cent %d;y;p_{T} [GeV]",i);
		hist_pt_y_proton_cent[i] = new TH2D(p_name.Data(),p_title.Data(),700,-3.0,3.0,500,0.0,3.5);
		p_name=Form("hist_pt_y_antiproton_cent_%d",i);
		p_title = Form("p_{T} [GeV] vs. y antiproton for cent %d;y;p_{T} [GeV]",i);
	    hist_pt_y_antiproton_cent[i] = new TH2D(p_name.Data(),p_title.Data(),700,-3.0,3.0,500,0.0,3.5);
		
		p_name=Form("proton_yield_%d",i);
		p_title = Form("proton yield for %d centrality;y;dN/dy",i);
		proton_yield[i]  = new TH1D(p_name.Data(),p_title.Data(), 60, -3.0, 3.0);
		p_name=Form("antiproton_yield_%d",i);
		p_title = Form("anti-proton yield for %d centrality;y;dN/dy",i);
		antiproton_yield[i]  = new TH1D(p_name.Data(),p_title.Data(), 60, -3.0, 3.0);
		p_name=Form("pip_yield_%d",i);
		p_title = Form("pi+ yield for %d centrality;y;dN/dy",i);
		pip_yield[i] = new TH1D(p_name.Data(),p_title.Data(), 60, -3.0, 3.0);
		p_name=Form("pim_yield_%d",i);
		p_title = Form("pi- yield for %d centrality;y;dN/dy",i);
		pim_yield[i] = new TH1D(p_name.Data(),p_title.Data(), 60, -3.0, 3.0);
		
		p_name=Form("v1_pt_cent_%d_protonplus",i);
		p_protonplusv1_pt_cent[i]  = new TProfile(p_name.Data(), p_name.Data(), 60, 0.0, 6.0);
		p_protonplusv1_pt_cent[i]->GetXaxis()->SetTitle("P_{T} (GeV)");
		p_protonplusv1_pt_cent[i]->GetYaxis()->SetTitle("v_{1}");
		p_name=Form("v1_pt_cent_%d_protonminus",i);
		p_protonminusv1_pt_cent[i]  = new TProfile(p_name.Data(), p_name.Data(), 60, 0.0, 6.0);
		p_protonminusv1_pt_cent[i]->GetXaxis()->SetTitle("P_{T} (GeV)");
		p_protonminusv1_pt_cent[i]->GetYaxis()->SetTitle("v_{1}");
		
		//pion stuff I put here because there's already a centrality loop:
		p_name=Form("pi_plus_v1_%d",i);
		p_title = Form("#pi^{+} v_{1} for %d centrality;y;v_{1}",i);
		p_pionplusv1_y_cent[i]=  new TProfile(p_name.Data(),p_title.Data(), 60, -3.0, 3.0);
		p_name=Form("pi_min_v1_%d",i);
		p_title = Form("#pi^{-} v_{1} for %d centrality;y;v_{1}",i);
		p_pionminusv1_y_cent[i]= new TProfile(p_name.Data(),p_title.Data(), 60, -3.0, 3.0);
		p_name=Form("hist_pion_yield_ratio_cent_%d",i);
		hist_pion_yield_ratio[i] = new TH1D(p_name.Data(),p_name.Data(),1000,0.0,10);
		hist_pion_yield_ratio[i]->GetXaxis()->SetTitle("#pi^{-}/#pi^{+}");
		hist_pion_yield_ratio[i]->GetYaxis()->SetTitle("Events");
		//kaon
		p_name=Form("p_kaonplusv1_y_%d",i);
		p_title = Form("k^{+} v_{1} for %d centrality;y;v_{1}",i);
		p_kaonplusv1_y_cent[i]=  new TProfile(p_name.Data(),p_title.Data(), 60, -3.0, 3.0);
		p_name=Form("p_kaonminusv1_y_%d",i);
		p_title = Form("k^{-} v_{1} for %d centrality;y;v_{1}",i);
		p_kaonminusv1_y_cent[i]= new TProfile(p_name.Data(),p_title.Data(), 60, -3.0, 3.0);
	}
	for(int h=0; h<9;h++)
	{	
		for(int i=0; i<9; i++)
		{
			for(int j = 0; j<VZ_BINS;j++)
			{
				p_name=Form("v1_y_vz_%d_cent_%d_protonplus_%d",j,i,h);
				p_protonplusv1_y_cent[j][h][i]  = new TProfile(p_name.Data(), p_name.Data(), 60, -3.0, 3.0);
				p_protonplusv1_y_cent[j][h][i]->GetXaxis()->SetTitle("y");
				p_protonplusv1_y_cent[j][h][i]->GetYaxis()->SetTitle("v_{1}");
				p_name=Form("v1_y_vz_%d_cent_%d_protonminus_%d",j,i,h);
				p_protonminusv1_y_cent[j][h][i]  = new TProfile(p_name.Data(), p_name.Data(), 60, -3.0, 3.0);
				p_protonminusv1_y_cent[j][h][i]->GetXaxis()->SetTitle("y");
				p_protonminusv1_y_cent[j][h][i]->GetYaxis()->SetTitle("v_{1}");
			}
			
			for(int j = 0; j<10;j++)//y
			{
				for(int k = 0;k<10;k++)//phi
				{
					for(int m=0;m<6;m++)
					{
						if(i!=0){continue;}
						p_name=Form("protonminus_v1_y_cent_%d_y_%d_phi_%d_p_%d_opt_%d",i,j,k,m,h);
						h_protonminus_y_cent[h][i][j][k][m]=new TH1D(p_name.Data(),p_name.Data(),500,0.0,3.0);
						h_protonminus_y_cent[h][i][j][k][m]->GetXaxis()->SetTitle("mass^{2} (GeV^{2})");
					}
				}
			}
		}
	}

	profile_pion_yield_ratio= new TProfile("profile_pion_yield_ratio","#pi^{-}/#pi^{+} vs. Centrality;Centrality;yield ratio",9,-0.5,8.5);
	
	for(int i = 0;i<10;i++)
	{
		p_name=Form("p_protonplusv1_y_pi_%d",i);
		p_title=Form("proton v1 vs y %d/10<#pi^{+}/(#pi^{+}+#pi^{-})<%d/10",i,i+1);
		p_protonplusv1_y_pip[i]  = new TProfile(p_name.Data(), p_title.Data(), 60, -3.0, 3.0);
		p_name=Form("p_protonminusv1_y_pi_%d",i);
		p_title=Form("anti proton v1 vs y %d/10<#pi^{+}/(#pi^{+}+#pi^{-})<%d/10",i,i+1);
		p_protonminusv1_y_pip[i]  = new TProfile(p_name.Data(), p_title.Data(), 60, -3.0, 3.0);
		
		p_name=Form("p_protonplusv1_y_mid_pip_%d",i);
		p_title=Form("proton v1 vs y mid centrality %d/10<#pi^{+}/(#pi^{+}+#pi^{-})<%d/10",i,i+1);
		p_protonplusv1_y_mid_pip[i]  = new TProfile(p_name.Data(), p_title.Data(), 60, -3.0, 3.0);
		p_name=Form("p_protonminusv1_y_mid_pip_%d",i);
		p_title=Form("anti proton v1 vs y mid centrality %d/10<#pi^{+}/(#pi^{+}+#pi^{-})<%d/10",i,i+1);
		p_protonminusv1_y_mid_pip[i]  = new TProfile(p_name.Data(), p_title.Data(), 60, -3.0, 3.0);
	}
	p_protonplusv1_y_mid_vz[0]  = new TProfile("p_protonplusv1_y_mid_vz0" ,"proton v1 vs y mid cent -145<V_{z}<-70 cm;y;p_{T} [GeV]"    ,60, -3.0, 3.0);
	p_protonminusv1_y_mid_vz[0] = new TProfile("p_protonminusv1_y_mid_vz0","antiproton v1 vs y mid cent -145<V_{z}<-70 cm;y;p_{T} [GeV]",60, -3.0, 3.0);
	p_protonplusv1_y_mid_vz[1]  = new TProfile("p_protonplusv1_y_mid_vz1" ,"proton v1 vs y mid cent -70<V_{z}<0 cm;y;p_{T} [GeV]"       ,60, -3.0, 3.0);
	p_protonminusv1_y_mid_vz[1] = new TProfile("p_protonminusv1_y_mid_vz1","antiproton v1 vs y mid cent -70<V_{z}<0 cm;y;p_{T} [GeV]"   ,60, -3.0, 3.0);
	p_protonplusv1_y_mid_vz[2]  = new TProfile("p_protonplusv1_y_mid_vz2" ,"proton v1 vs y mid cent 0<V_{z}<70 cm;y;p_{T} [GeV]"        ,60, -3.0, 3.0);
	p_protonminusv1_y_mid_vz[2] = new TProfile("p_protonminusv1_y_mid_vz2","antiproton v1 vs y mid cent 0<V_{z}<70 cm;y;p_{T} [GeV]"    ,60, -3.0, 3.0);
	p_protonplusv1_y_mid_vz[3]  = new TProfile("p_protonplusv1_y_mid_vz3" ,"proton v1 vs y mid cent 70<V_{z}<145 cm;y;p_{T} [GeV]"      ,60, -3.0, 3.0);
	p_protonminusv1_y_mid_vz[3] = new TProfile("p_protonminusv1_y_mid_vz3","antiproton v1 vs y mid cent 70<V_{z}<145 cm;y;p_{T} [GeV]"  ,60, -3.0, 3.0);
	
	
	// if(dataset == "3p85_Gev_2018"){cout<<"HEY!!! YOU NEED EFFICIENCY FILE!!!!"<<endl; do_tof_efficiency = 0;}
	// else if(dataset == "7p7_Gev_2021"){efficiency_tof_file = TFile::Open("/star/u/educkwort/excess_proton/all_energies/flow/flow_7p7_2024-05-06_v1_vs_y.root");}
	// else if(dataset == "9p2_GeV_2020"){efficiency_tof_file = TFile::Open("/star/u/educkwort/excess_proton/all_energies/flow/flow_9p2_gev_2024-05-06.root");} 
	// else if(dataset == "11p5_gev_2020"){efficiency_tof_file = TFile::Open("/star/u/educkwort/excess_proton/all_energies/flow/flow_11p5_2024-05-06.root");}
	// else if(dataset == "14p6_gev_2019"){efficiency_tof_file = TFile::Open("/star/u/educkwort/excess_proton/all_energies/flow/flow_14p6_2024-05-07_central_peripheral.root");}
	// else if(dataset == "17p3_gev_2021"){efficiency_tof_file = TFile::Open("/star/u/educkwort/excess_proton/all_energies/flow/flow_17p3_2024-05-07.root");}
	// else if(dataset == "19p6_Gev_2019"){efficiency_tof_file = TFile::Open("/star/u/educkwort/excess_proton/all_energies/flow/flow_19p6_gev_2024-05-31.root");}
	// else if(dataset == "27_Gev_2018"){efficiency_tof_file = TFile::Open("/star/u/educkwort/excess_proton/all_energies/flow/flow_27_gev_2024-05-02.root");}
	if(do_tof_efficiency)
	{
		efficiency_tof_file = TFile::Open(efficiency_tof_file_name.Data());
		efficiency_tof_file->GetObject("tof_efficiency_tot",hist_eta_y_tpc_hit);hist_eta_y_tpc_hit->SetName("hist_eta_y_tpc_hit");
		efficiency_tof_file->GetObject("tof_efficiency_match",hist_eta_y_tof_match_hit);hist_eta_y_tof_match_hit->SetName("hist_eta_y_tof_match_hit");
	}
	// if(dataset == "9p2_GeV_2020")
	// {
		// efficiency_proton_file = TFile::Open("/star/u/educkwort/excess_proton/all_energies/flow/OUT_9p2_gev_proton_picoDst.root");
		// efficiency_aproton_file = TFile::Open("/star/u/educkwort/excess_proton/all_energies/flow/OUT_9p2_gev_antiproton_picoDst.root");
		// for(Int_t i = 0; i<9; i++)
		// {
			// p_name=Form("acceptance_MC_cent_%d",i);
			// efficiency_proton_file->GetObject(p_name.Data(),hist_pt_y_mc_proton_cent[i]);
			// p_title = Form("acceptance_MC_proton_cent_%d",i);
			// hist_pt_y_mc_proton_cent[i]->SetName(p_title.Data());
			// efficiency_aproton_file->GetObject(p_name.Data(),hist_pt_y_mc_aproton_cent[i]);
			// p_title = Form("acceptance_MC_aproton_cent_%d",i);
			// hist_pt_y_mc_aproton_cent[i]->SetName(p_title.Data());
			
			// p_name=Form("acceptance_match_cent_%d",i);
			// efficiency_proton_file->GetObject(p_name.Data(),hist_pt_y_match_proton_cent[i]);
			// p_title = Form("acceptance_match_proton_cent_%d",i);
			// hist_pt_y_match_proton_cent[i]->SetName(p_title.Data());
			// efficiency_aproton_file->GetObject(p_name.Data(),hist_pt_y_match_aproton_cent[i]);
			// p_title = Form("acceptance_match_aproton_cent_%d",i);
			// hist_pt_y_match_aproton_cent[i]->SetName(p_title.Data());
		// }
	// }
	
	//pi ratio stuff
	pi_corr_yield= new TH2D("pi_corr_yield","#pi^{-} vs #pi^{-} multiplicity;#pi^{+};#pi^{-} ",500,0.0,500,500,0.0,500);
	pip_b_tot_pi_ratio = new TH1D("pip_b_tot_pi_ratio","#pi^{+}/(#pi^{+}+#pi^{-}) distribution;#pi^{+}/(#pi^{+}+#pi^{-});entries/1000 bins",1000,0.0,1.0);
	pip_b_tot_pi_cent = new TH2D("pip_b_tot_pi_cent","#pi^{+}/(#pi^{+}+#pi^{-}) distribution vs. centrality;centrality;#pi^{+}/(#pi^{+}+#pi^{-})",9,-0.5,8.5,1000,0.0,1.0);
	pip_b_tot_pi_cent_p = new TProfile("pip_b_tot_pi_cent_p","<#pi^{+}/(#pi^{+}+#pi^{-})> distribution vs. centrality;centrality;#pi^{+}/(#pi^{+}+#pi^{-})",9,-0.5,8.5);
	
	//centrality:
	hist_cent = new TH1D("hist_cent","Centrality",9,-0.5,8.5);
    hist_cent->GetXaxis()->SetTitle("Centrality bin");
    hist_cent->GetYaxis()->SetTitle("# of events");
	
	hist_refmult = new TH1D("hist_refmult","hist_refmult",1000,-0.5,999.5);
	hist_refmult->GetXaxis()->SetTitle("refmult");
    hist_refmult->GetYaxis()->SetTitle("# of events");
	
	refmult_v_tofmult = new TH2D("refmult_v_tofmult","refmult vs. tofmult ;refmult;tofmult",1000,-0.5,999.5,1000,-0.5,999.5);
	
	//dedx:
	hist_dEdx = new TH2D("hist_dEdx","<dE/dx> vs |p|/q ;|p|/q (GeV) ; <dE/dx> ; (keV/cm)",500,-3.0,3.0,500,0.0,10.0);
	hist_dEdx_proton = new TH2D("hist_dEdx_proton","dE/dx vs |p|/q ;|p|/q (GeV/c); <dE/dx> (keV/cm) ",500,-3.0,3.0,500,0.0,10.0);
	hist_dEdx_pion = new TH2D("hist_dEdx_pion","dE/dx vs |p|/q;|p|/q (GeV/c);<dE/dx> (keV/cm)",500,-3.0,3.0,500,0.0,10.0);
	hist_dEdx_kaon = new TH2D("hist_dEdx_kaon","dE/dx vs |p|/q;|p|/q (GeV/c);<dE/dx> (keV/cm)",500,-3.0,3.0,500,0.0,10.0);
	hist_nsigma_p_proton = new TH2D("hist_nsigma_p_proton","hist_n#sigma_p_proton; |p| (GeV); <n#sigma_proton>",500,-3.0,3.0,500,0,10.0);
	hist_nsigma_p_aproton = new TH2D("hist_nsigma_p_aproton","hist_n#sigma_p_aproton; |p| (GeV); <n#sigma_aproton>",500,-3.0,3.0,500,0,10.0);
	p_nsigma_proton = new TProfile("p_nsigma_proton","p_n#sigma_proton ;|p| (GeV); <n#sigma_proton>",30,0.0,2.8);
	p_nsigma_aproton = new TProfile("p_nsigma_aproton","p_n#sigma_aproton ;|p| (GeV); <n#sigma_aproton>",30,0.0,2.8);
	hist_nsigma_p_pip = new TH2D("hist_nsigma_p_pip","hist n#sigma_{#pi+};|p| (GeV);n#sigma_{#pi+}",500,-3.0,3.0,500,0,10.0);
	hist_nsigma_p_pim = new TH2D("hist_nsigma_p_pim","hist n#sigma_{#pi-};|p| (GeV);n#sigma_{#pi-}",500,-3.0,3.0,500,0,10.0);
	p_nsigma_pip = new TProfile("p_nsigma_pip","n#sigma_{#pi+};|p| (GeV);<n#sigma_{#pi+}>",30,0.0,2.8);
	p_nsigma_pim = new TProfile("p_nsigma_pim","<n#sigma_{#pi-}>;|p| (GeV);<n#sigma_{#pi-}>",30,0.0,2.8);
	
	//DCA:
	hist_dca_vs_p_proton = new TH2D("hist_dca_vs_p_proton","DCA of protons;|P| (GeV);DCA (cm)",500,0.0,5.0,300,0,3.0);
	hist_dca_vs_p_aproton = new TH2D("hist_dca_vs_p_aproton","DCA of anti-protons;|P| (GeV);DCA (cm)",500,0.0,5.0,300,0,3.0);
	
	//TOF:
	hist_tof_pid=new TH2D("TOF","TOF ; |p|/q (GeV) ; 1/#beta",500,-4.0,4.0,500,0.0,5.0);
	hist_tof_pid_proton=new TH2D("TOF_proton","TOF proton;|p|/q (GeV/c);1/#beta",500,-4.0,4.0,500,0.0,5.0);
	hist_tof_pid_pion=new TH2D("TOF_pion","TOF #pi;|p|/q (GeV/c);1/#beta",500,-4.0,4.0,500,0.0,5.0);
	hist_tof_pid_kaon=new TH2D("TOF_kaon","TOF k;|p|/q (GeV/c);1/#beta",500,-4.0,4.0,500,0.0,5.0);
	
	hist_mass2_pid =new TH2D("hist_mass2_pid","Mass^{2} ; |p|/q (GeV) ; m^{2} (GeV^{2})",500,-4.0,4.0,500,0.0,5.0);
	hist_mass2_pid_proton =new TH2D("hist_mass2_pid_proton","Mass^{2} proton; |p|/q (GeV) ; m^{2} (GeV^{2})",500,-4.0,4.0,500,0.0,5.0);
	hist_mass2_pid_pion =new TH2D("hist_mass2_pid_pion","Mass^{2} #pi; |p|/q (GeV) ; m^{2} (GeV^{2})",500,-4.0,4.0,500,0.0,5.0);
	hist_mass2_pid_kaon =new TH2D("hist_mass2_pid_kaon","Mass^{2} kaon; |p|/q (GeV) ; m^{2} (GeV^{2})",500,-4.0,4.0,500,0.0,5.0);
	
	hist_mass2_proton = new TH2D("hist_mass2_proton","Mass^{2} Proton ; (m/q)^{2} (GeV^{2}) ; P(GeV)",500,0.0,3.0,500,0.0,5.0);
	hist_mass2_aproton = new TH2D("hist_mass2_aproton","Mass^{2} Anti-Proton ; (m/q)^{2} (GeV^{2}) ; P(GeV)",500,0.0,3.0,500,0.0,5.0);
	
	tof_efficiency_tot = new TH2D("tof_efficiency_tot","TOF Efficiency total;#eta,P_{T} (GeV)",700,-3.0,3.0,500,0.0,3.5);
	tof_efficiency_match = new TH2D("tof_efficiency_match","TOF Efficiency Match;#eta,P_{T} (GeV)",700,-3.0,3.0,500,0.0,3.5);
	
	//purity:
	hist_nsigma_proton = new TH1D("hist_nsigma_proton","n#sigma Proton;n#sigma_proton;Tracks",500,-10.0,10.0);
	hist_nsigma_aproton = new TH1D("hist_nsigma_aproton","n#sigma Anti-Proton;n#sigma_proton;Tracks",500,-10.0,10.0);
	
	//accuracy:
	hist_proton_accuracy = new TH2D("hist_proton_accuracy","Proton Accuracy;n#sigma_proton;(m/q)^{2} (GeV^{2})",500,-3.0,3.0,500,0.0,3.0);
	hist_aproton_accuracy = new TH2D("hist_aproton_accuracy","Anti-Proton Accuracy;n#sigma_proton;(m/q)^{2} (GeV^{2})",500,-3.0,3.0,500,0.0,3.0);
	
	hist_pt_eta = new TH2D("hist_pt_eta","hist_pt_eta;#eta;pt/q",500,-2.0,2,1000,-5.0,5.0);
	hist_pt_centrality = new TH2D("hist_pt_centrality","hist_pt_centrality;p_{T}/q;centrality",500,-5.0,5.0,9,-0.5,8.5);
	
	//ch v1
	for(Int_t i=0; i<9;i++)
	{
		p_name=Form("ch_v1_v_eta_cent_%d",i);
		ch_v1_eta_cent[i] = new TProfile(p_name.Data(),p_name.Data(), 60, -3.0, 3.0);
		ch_v1_eta_cent[i]->GetXaxis()->SetTitle("#eta");
		ch_v1_eta_cent[i]->GetYaxis()->SetTitle("v1");
	}
	
	//proton acceptance:
	hist_pt_y_proton = new TH2D("hist_pt_y_proton","p_{T} [GeV] vs. y proton;y;p_{T} [GeV]",700,-3.0,3.0,500,0.0,3.5);
	hist_pt_y_antiproton = new TH2D("hist_pt_y_antiproton","p_{T} [GeV] vs. y antiproton;y;p_{T} [GeV]",700,-3.0,3.0,500,0.0,3.5);
	hist_pt_y_pip= new TH2D("hist_pt_y_pip","p_{T} [GeV] vs. y #pi+;y;p_{T} [GeV]",700,-3.0,3.0,500,0.0,3.5);
	hist_pt_y_pim = new TH2D("hist_pt_y_pim","p_{T} [GeV] vs. y #pi-;y;p_{T} [GeV]",700,-3.0,3.0,500,0.0,3.5);
	
	
	hist_pt_y_proton_vz[0] = new TH2D("hist_pt_y_proton_vz0","p_{T} [GeV] vs. y -145<V_{z}<-70 cm;y;p_{T} [GeV]",700,-3.0,3.0,500,0.0,3.5);
	hist_pt_y_antiproton_vz[0]  = new TH2D("hist_pt_y_antiproton_vz0","p_{T} [GeV] vs. y -145<V_{z}<-70 cm;y;p_{T} [GeV]",700,-3.0,3.0,500,0.0,3.5);
	hist_pt_y_proton_vz[1] = new TH2D("hist_pt_y_proton_vz1","p_{T} [GeV] vs. y -70<V_{z}<0 cm;y;p_{T} [GeV]",700,-3.0,3.0,500,0.0,3.5);
	hist_pt_y_antiproton_vz[1]  = new TH2D("hist_pt_y_antiproton_vz1","p_{T} [GeV] vs. y -70<V_{z}<0 cm;y;p_{T} [GeV]",700,-3.0,3.0,500,0.0,3.5);
	hist_pt_y_proton_vz[2] = new TH2D("hist_pt_y_proton_vz2","p_{T} [GeV] vs. y 0<V_{z}<70 cm;y;p_{T} [GeV]",700,-3.0,3.0,500,0.0,3.5);
	hist_pt_y_antiproton_vz[2]  = new TH2D("hist_pt_y_antiproton_vz2","p_{T} [GeV] vs. y 0<V_{z}<70 cm;y;p_{T} [GeV]",700,-3.0,3.0,500,0.0,3.5);
	hist_pt_y_proton_vz[3] = new TH2D("hist_pt_y_proton_vz3","p_{T} [GeV] vs. y 70<V_{z}<145 cm;y;p_{T} [GeV]",700,-3.0,3.0,500,0.0,3.5);
	hist_pt_y_antiproton_vz[3]  = new TH2D("hist_pt_y_antiproton_vz3","p_{T} [GeV] vs. y 70<V_{z}<145 cm;y;p_{T} [GeV]",700,-3.0,3.0,500,0.0,3.5);
	
	//kaon acceptance:
	hist_pt_y_kp= new TH2D("hist_pt_y_kp","p_{T} [GeV] vs. y k+;y;p_{T} [GeV]",700,-3.0,3.0,500,0.0,3.5);
	hist_pt_y_km = new TH2D("hist_pt_y_km","p_{T} [GeV] vs. y k-;y;p_{T} [GeV]",700,-3.0,3.0,500,0.0,3.5);

	for(Int_t i=0; i<9;i++)
	{
		p_name=Form("eta_v_pt_pos_cent_%d",i);
		eta_v_pt_pos_cent[i]=new TProfile(p_name.Data(),p_name.Data(),40,-2.0,2.0);
		p_name=Form("eta_v_pt_neg_cent_%d",i);
		eta_v_pt_neg_cent[i]=new TProfile(p_name.Data(),p_name.Data(),40,-2.0,2.0);
	}
	baryon_yield_histogram= new TH1D("baryon_yield_histogram","baryon_yield_histogram;|y|;dN/dy",200,0,6.0);
    cout << "End of Histograms" << endl;
    return kStOK;
}
//__________________________________________________________________________________
void Shift::Clear(Option_t *opt)
{
    StMaker::Clear();
}

//__________________________________________________________________________________
Int_t Shift::Finish() {
    cout << "Shift::Finish()\n";
    //===============================
    //  Write Histograms
    //===============================

    File->cd();
	baryon_yield_histogram->Write();
	hist_VyVx->Write();
	hist_Vz->Write();
	for(Int_t i=0;i<2;i++)
	{
		hist_epd_ab_psi_raw[i]->Write();
		hist_epd_c_psi_raw[i]->Write();
		hist_epd_d_psi_raw[i]->Write();
		hist_epd_abcd_psi_raw[i]->Write();
		hist_epd_cd_psi_raw[i]->Write();
		hist_epd_ab_psi_recentered[i]->Write();
		hist_epd_c_psi_recentered[i]->Write();
		hist_epd_d_psi_recentered[i]->Write();
		hist_epd_abcd_psi_recentered[i]->Write();
		hist_epd_cd_psi_recentered[i]->Write();
		hist_epd_ab_psi_flattened[i]->Write();
		hist_epd_c_psi_flattened[i]->Write();
		hist_epd_d_psi_flattened[i]->Write();
		hist_epd_abcd_psi_flattened[i]->Write();
		hist_epd_cd_psi_flattened[i]->Write();
		for(Int_t j=0;j<4;j++)
		{
			hist_Qr_epd[j][i]->Write();
		}
	}
	hist_epd_ab_psi_raw_ew->Write();
	hist_epd_ab_psi_recentered_ew->Write();
	hist_epd_ab_psi_flattened_ew->Write();

	hist_epd_psi_raw_ew->Write();
	hist_epd_psi_recentered_ew->Write();
	hist_epd_psi_flattened_ew->Write();
	for(int i =0;i<9;i++)
	{
		for(int j =0;j<VZ_BINS;j++)
		{
			v1_vs_eta_east[i][j]->Write();
			v1_vs_eta_west[i][j]->Write();
		}
	}
	hist_tpc_A_psi_raw->Write();
	hist_tpc_B_psi_raw->Write();
	
	hist_tpc_A_psi_recentered->Write();
	hist_tpc_B_psi_recentered->Write();
	hist_tpc_AB_psi_recentered->Write();
	hist_tpc_A_psi_flattened->Write();
	hist_tpc_B_psi_flattened->Write();
	hist_tpc_AB_psi_flattened->Write();
	
	//correlation plots
	profile_correlation_abc->Write();
	profile_correlation_abd->Write();
	profile_correlation_cd->Write();
	
	profile_correlation_eABtB->Write();
	profile_correlation_eCtB->Write();
	profile_correlation_eDtB->Write();
	
	profile_correlation_eABtA->Write();
	profile_correlation_eCtA->Write();
	profile_correlation_eDtA->Write();
	
	profile_correlation_abew->Write();
	profile_correlation_cew->Write();
	profile_correlation_dew->Write();
	
	profile_correlation_ew->Write();
	for(int i = 0;i<VZ_BINS;i++)
	{
		profile_correlation_tpc_ew_vz[i]->Write();
		profile_correlation_ew_vz[i]->Write();
		profile_correlation_eABw_vz[i]->Write();
		profile_correlation_eCDw_vz[i]->Write();
		profile_correlation_eABeCD_vz[i]->Write();
		profile_correlation_ewAB_vz[i]->Write();
		profile_correlation_ewCD_vz[i]->Write();
		profile_correlation_wABwCD_vz[i]->Write();
		profile_correlation_eDw_vz[i]->Write();
		profile_correlation_eABeD_vz[i]->Write();
		profile_correlation_ewD_vz[i]->Write();
		profile_correlation_wABwD_vz[i]->Write();
		profile_correlation_TPCEPDe_vz[i]->Write();
		profile_correlation_TPCEPDw_vz[i]->Write();
	
	}
	
	for(int i = 0;i<9;i++)
	{
		correlation_ew_cent[i]->Write();
	}
	
	profile_correlation_tpc_AB->Write();
	
	correlation_abc->Write();
	correlation_abd->Write();
	correlation_cd->Write();
	
	correlation_eABtB->Write();
	correlation_eCtB->Write();
	correlation_eDtB->Write();
	
	correlation_eABtA->Write();
	correlation_eCtA->Write();
	correlation_eDtA->Write();
	
	correlation_abew->Write();
	correlation_cew->Write();
	correlation_dew->Write();

	correlation_ew->Write();
	
	correlation_tpc_AB->Write();	
	
	// proton_tree->Write();
	
	hist_cent->Write();
	hist_refmult->Write();
	refmult_v_tofmult->Write();

	pip_b_tot_pi_ratio->Write();
	pip_b_tot_pi_cent->Write();
	pip_b_tot_pi_cent_p->Write();
	if(dataset == "27_Gev_2018")
	{
		mEpFinder->Finish();
	}
	if(ep_setting == "recenter")
	{
		for(Int_t i = 0; i<VZ_BINS;i++)
		{
			for(Int_t j = 0; j<5; j++)
			{
				for(Int_t k=0; k<2;k++)
				{
					epd_recenterx[i][j][k]->Write();
					epd_recentery[i][j][k]->Write();
				}
			}
			for(Int_t j = 0; j<2; j++)
			{
				tpc_recenterx[i][j]->Write();
				tpc_recentery[i][j]->Write();
			}
		}
	}
	if(ep_setting == "flatten")
	{
		for(Int_t h=0;h<VZ_BINS;h++)
		{

				
			for(Int_t j=0;j<max_flat;j++)
			{			
				epd_flatten_ab_ew_cos[h][j]->Write();
				epd_flatten_ab_ew_sin[h][j]->Write();
				
				epd_flatten_ew_cos[h][j]->Write();
				epd_flatten_ew_sin[h][j]->Write();
				for(Int_t i=0;i<5;i++)
				{
					for(Int_t k=0;k<2;k++)
					{			
						epd_flatten_cos[h][i][j][k]->Write();
						epd_flatten_sin[h][i][j][k]->Write();
					}
				}
				for(Int_t i=0;i<3;i++)
				{
					tpc_flatten_cos[h][i][j]->Write();
					tpc_flatten_sin[h][i][j]->Write();
				}
			}
		}
	}
	if(ep_setting != "flow"){return kStOK;}
	
	nhits_vs_eta->Write();
	hist_dEdx->Write();
	hist_dEdx_proton->Write();
	hist_dEdx_pion->Write();
	hist_dEdx_kaon->Write();
	hist_dca_vs_p_proton->Write();
	hist_dca_vs_p_aproton->Write();
	hist_pt_eta->Write();
	hist_pt_centrality->Write();
	
	hist_nsigma_p_proton->Write();
	hist_nsigma_p_aproton->Write();
	p_nsigma_proton->Write();
	p_nsigma_aproton->Write();
	
	hist_nsigma_p_pip->Write();
	hist_nsigma_p_pim->Write();
	p_nsigma_pip->Write();
	p_nsigma_pim->Write();
	
	hist_nsigma_proton->Write();
	hist_nsigma_aproton->Write();
	hist_mass2_proton->Write();
	hist_mass2_aproton->Write();
	hist_proton_accuracy->Write();
	hist_aproton_accuracy->Write();
	
	hist_tof_pid->Write();
	hist_tof_pid_proton->Write();
	hist_tof_pid_pion->Write();
	hist_tof_pid_kaon->Write();
	hist_mass2_pid->Write();
	hist_mass2_pid_proton->Write();
	hist_mass2_pid_pion->Write();
	hist_mass2_pid_kaon->Write();
	tof_efficiency_tot->Write();
	tof_efficiency_match->Write();
	
	hist_pt_y_proton->Write();
	hist_pt_y_antiproton->Write();
	hist_pt_y_pip->Write();
	hist_pt_y_pim->Write();
	hist_pt_y_kp->Write();
	hist_pt_y_km->Write();
	hist_pt_y_proton_vz[0]->Write();
	hist_pt_y_antiproton_vz[0]->Write();
	hist_pt_y_proton_vz[1]->Write();
	hist_pt_y_antiproton_vz[1]->Write();
	hist_pt_y_proton_vz[2]->Write();
	hist_pt_y_antiproton_vz[2]->Write();
	hist_pt_y_proton_vz[3]->Write();
	hist_pt_y_antiproton_vz[3]->Write();
	for(Int_t i=0; i<9;i++)
	{
		eta_v_pt_pos_cent[i]->Write();
		eta_v_pt_neg_cent[i]->Write();
	}
	
	for(int i=0; i<9; i++)
    {
		hist_pt_y_proton_cent[i]->Write();
		hist_pt_y_antiproton_cent[i]->Write();
		ch_v1_eta_cent[i]->Write();
		
		proton_yield[i]->Write();
		antiproton_yield[i]->Write();
		pip_yield[i]->Write();
		pim_yield[i]->Write();
		p_pionplusv1_y_cent[i]->Write();
		p_pionminusv1_y_cent[i]->Write();
		p_kaonplusv1_y_cent[i]->Write();
		p_kaonminusv1_y_cent[i]->Write();
		// p_protonplusv1_yn_cent[i]->Write();
		// p_protonminusv1_yn_cent[i]->Write();
		// p_protonplusv1_eta_cent[i]->Write();
		// p_protonminusv1_eta_cent[i]->Write();
		// p_protonplusv1_etan_cent[i]->Write();
		// p_protonminusv1_etan_cent[i]->Write();
		
        p_protonplusv1_pt_cent[i]->Write();
        p_protonminusv1_pt_cent[i]->Write();
		hist_pion_yield_ratio[i]->Write();
	}
	for(int h = 0;h<9;h++)
	{
		for(int i=0; i<9; i++)
		{
			for(int j = 0; j<VZ_BINS;j++)
			{
				p_protonplusv1_y_cent[j][h][i]->Write();
				p_protonminusv1_y_cent[j][h][i]->Write();
			}
			for(int j = 0; j<10;j++)//y
			{
				for(int k = 0;k<10;k++)//phi
				{
					for(int m=0;m<6;m++)
					{
						if(i!=0){continue;}
						h_protonminus_y_cent[h][i][j][k][m]->Write();
					}
				}
			}
		}
	}
	for(int i = 0;i<10;i++)
	{
		p_protonplusv1_y_pip[i]->Write();
		p_protonminusv1_y_pip[i]->Write();
		p_protonplusv1_y_mid_pip[i]->Write();
		p_protonminusv1_y_mid_pip[i]->Write();
	}
	for(int i = 0;i<4;i++)
	{
		p_protonplusv1_y_mid_vz[i]->Write();
		p_protonminusv1_y_mid_vz[i]->Write();
	}
	pi_corr_yield->Write();
	profile_pion_yield_ratio->Write();
	//QA_tree->Write();
    return kStOK;
}
//__________________________________________________________________________________
// Int_t Shift::GetRunIndex( const Int_t run ) {
    // Int_t runindex = -999;
    // for(Int_t i=0; i<nrun; i++){
        // if(run==numbers[i]) runindex = i;
    // }
    // return runindex;
// }
Int_t Shift::GetRunIndex( const Int_t run ) {
    Int_t runindex = -999;
	Int_t run_trunc = run/1000;
    for(Int_t i=0; i<nrun; i++){
        if(run_trunc==numbers[i]) runindex = i;
    }
    return runindex;
}
//---------------------------------------------------------------------------------
Int_t Shift::Centrality(int gRefMult )
{
    int centrality;
    //int centFull declared in run.h
	
    if      (gRefMult>=centFull[8]) centrality=8;// 0-5%
    else if (gRefMult>=centFull[7]) centrality=7;// 5-10%
    else if (gRefMult>=centFull[6]) centrality=6;//10-20%
    else if (gRefMult>=centFull[5]) centrality=5;//20-30%
    else if (gRefMult>=centFull[4]) centrality=4;//30-40%
    else if (gRefMult>=centFull[3]) centrality=3;//40-50%
    else if (gRefMult>=centFull[2]) centrality=2;//50-60%
    else if (gRefMult>=centFull[1]) centrality=1;//60-70%
    else if (gRefMult>=centFull[0]) centrality=0;//70-80%
    else centrality = 9;

    return centrality;
}


//---------------------------------------------------------------------------------
//__________________________________________________________________________________
Int_t Shift::Make() {
    //Begining of Event loop  

    //------------------------------------------------------------------
    if(!mPicoDstMaker) {LOG_WARN << " No PicoDstMaker! Skip! " << endm; return kStWarn;}
    mPicoDst = mPicoDstMaker->picoDst();
    if(!mPicoDst) {LOG_WARN << " No PicoDst! Skip! " << endm;return kStWarn;}
    picoEvent = (StPicoEvent*)mPicoDst->event();
    if( !picoEvent ){LOG_WARN << " No PicoEvent! Skip! " << endm; return kStWarn;}
    //------------------------------------------------------------------
    TVector3 pVertex = picoEvent->primaryVertex();
	vx=pVertex.X(); vy=pVertex.Y(); vz =pVertex.Z();
	vr= sqrt(pVertex.X()*pVertex.X()+pVertex.Y()*pVertex.Y());
	run_number = picoEvent->runId();
    event_number = picoEvent->eventId();
    int runindex = GetRunIndex(run_number);
	const Int_t nTrack = mPicoDst->numberOfTracks();
	Int_t vz_sign = 0; if(vz<0){vz_sign = 1;}	
	refmult= picoEvent->refMult();
	int NtofMatch = picoEvent->nBTOFMatch();
	
	hist_VyVx->Fill(pVertex.X(),pVertex.Y());
	
	if(!isGoodTrigger(picoEvent)){return 0;}
    if(!isGoodEvent(picoEvent)) {return 0;}
	//------ remove bad run---------
	if(dataset == "19p6_Gev_2019" || dataset == "27_Gev_2018" ||dataset == "9p2_GeV_2020"|| dataset == "14p6_gev_2019"||
	   dataset == "17p3_gev_2021" || dataset == "11p5_gev_2020")
	{
		mRefMultCorrUtil->init(run_number);
		if(mRefMultCorrUtil->isBadRun( run_number ) ) {return kStOk;}
	}
    else
	{
		for(int i=0; i<n_badrun; i++)
		{
			if(run_number == badrun[i]) return kStOK;
		}
	}
	
	hist_Vz->Fill(vz);

    //----- Remove pile up events ---------
	if(dataset == "3p85_Gev_2018")
	{
		mPileupTool->initEvent(mPicoDst);
		if(mPileupTool ->isPileupEPD(0)){return 0;}
	}
	else if(dataset == "19p6_Gev_2019" || dataset == "27_Gev_2018" ||dataset == "9p2_GeV_2020"|| dataset == "14p6_gev_2019"||
	   dataset == "17p3_gev_2021" || dataset == "11p5_gev_2020")
	{
		if (mRefMultCorrUtil->isPileUpEvent( refmult, NtofMatch, vz ) ) return kStOK;
		mRefMultCorrUtil->initEvent(refmult, vz, picoEvent->ZDCx());
	}
	else
	{
		refmult = Pileup_rejection(refmult, NtofMatch, vz, picoEvent);if(refmult ==-1){return kStOK;}
	}
	Int_t centrality = 0;
	float gweight = 1.0;
	if(dataset == "19p6_Gev_2019" || dataset == "27_Gev_2018" ||dataset == "9p2_GeV_2020"|| dataset == "14p6_gev_2019"||
	   dataset == "17p3_gev_2021" || dataset == "11p5_gev_2020")
	{
		if (mRefMultCorrUtil->getCentralityBin16() < 0 ||
			mRefMultCorrUtil->getCentralityBin9() < 0) return kStOK;
		centrality = mRefMultCorrUtil->getCentralityBin9();       // 0: 70-80%, 1: 60-70%, ..., 6: 10-20%, 7: 5-10%, 8: 0-5%
	}
	else if(dataset == "3p85_Gev_2018")
	{
		 refmult = mPileupTool->get_refMultPrim();
		 centrality  = mPileupTool->get_centrality9();  if(centrality > 8 || centrality < 0){return kStOK;}
		 gweight     = mPileupTool->get_centralityWeight();
	}	
	else
	{
		centrality = Centrality(refmult);
	}
	if(centrality > 8 || centrality < 0){return kStOK;}
	refmult_v_tofmult->Fill(refmult,NtofMatch);
	hist_cent->Fill(centrality);
	hist_refmult->Fill(refmult);
	
	
	//-----------------Get TPC information------------------------------
	double Qx_rawep_TPC[2]={0}, Qy_rawep_TPC[2]={0}, weight_TPC[2] = {0}, psi1rawep_TPC[2] = {0};//0 is B , 1 is A, I know...
	double Qx_recen_TPC[2]={0}, Qy_recen_TPC[2]={0}, psi1recen_TPC[3] = {0};//0 is B , 1 is A, I know...
    double psi1shift_TPC[3]={0};
	// double mult_pionplus= 0;
	// double mult_pionminus=0;
	Int_t ep_vz_bin = 0;
	if(VZ_BINS==4)
	{
		if(vz<-70){ep_vz_bin = 0;}
		else if(vz<0){ep_vz_bin = 1;}
		else if(vz<70){ep_vz_bin = 2;}
		else{ep_vz_bin = 3;}
	}
	else if(VZ_BINS == 2){ep_vz_bin =vz_sign;}
	if(dataset != "27_Gev_2018")
	{
		for (Int_t itr=0;itr<nTrack;itr++) 
		{
			const StPicoTrack *ptrk = (StPicoTrack*)mPicoDst->track(itr);
			if(!ptrk)  continue;
			if(!ptrk->isPrimary())  continue;  // now selecting primary tracks
			if(!isGoodTrack(ptrk))  continue;
			
			Float_t pt  = ptrk->pMom().Perp();
			Float_t phi = ptrk->pMom().Phi();
			Float_t eta = ptrk->pMom().PseudoRapidity();
			for(int i_group = 0;i_group<2;i_group++)
			{
				if(i_group == 0 && !(eta > 0.01 && eta < 1.0)  ){continue;}// 0: TPC-B
				if(i_group == 1 && !(eta > -1.0 && eta < -0.01)){continue;}// 1: TPC-A
				Qx_rawep_TPC[i_group] += pt*(cos(1.0*phi));
				Qy_rawep_TPC[i_group] += pt*(sin(1.0*phi));
				if(ep_setting !="recenter")
				{
					Qx_recen_TPC[i_group] += pt*(cos(1.0*phi)) - gettpc_recenterx[ep_vz_bin][i_group]->GetBinContent(gettpc_recenterx[ep_vz_bin][i_group]->FindBin(runindex,centrality));
					Qy_recen_TPC[i_group] += pt*(sin(1.0*phi)) - gettpc_recentery[ep_vz_bin][i_group]->GetBinContent(gettpc_recentery[ep_vz_bin][i_group]->FindBin(runindex,centrality));
				}
				else
				{
					tpc_recenterx[ep_vz_bin][i_group]->Fill(runindex,centrality,pt*(cos(1.0*phi)));
					tpc_recentery[ep_vz_bin][i_group]->Fill(runindex,centrality,pt*(sin(1.0*phi)));
				}
				weight_TPC[i_group]++;
			}
		}
		double Qy_recen_TPC_ew = Qy_recen_TPC[1] - Qy_recen_TPC[0];
		double Qx_recen_TPC_ew = Qx_recen_TPC[1] - Qx_recen_TPC[0];
		double weight_TPC_ew   = weight_TPC[1]   +   weight_TPC[0];
		for(int i_group = 0;i_group<2;i_group++)
		{
			//if(dataset == "4p5_gev_FXT_2020" || dataset=="3p85_Gev_2018"){if(weight_TPC[i_group]==0.) {continue;}}
			if(Qy_rawep_TPC[i_group]==0.||Qx_rawep_TPC[i_group]==0.||weight_TPC[i_group]==0.
				||(ep_setting !="recenter" &&(Qx_recen_TPC[i_group]==0.|| Qy_recen_TPC[i_group]==0.))){return 0;}
			
			Qx_rawep_TPC[i_group] /= weight_TPC[i_group]; Qy_rawep_TPC[i_group] /= weight_TPC[i_group];
			Qx_recen_TPC[i_group] /= weight_TPC[i_group]; Qy_recen_TPC[i_group] /= weight_TPC[i_group];

			psi1rawep_TPC[i_group] = (1.0/1.0) * atan2(Qy_rawep_TPC[i_group], Qx_rawep_TPC[i_group]);
			if(zero2twopi(psi1rawep_TPC[i_group])){return 0;}
			if(ep_setting !="recenter")
			{
				psi1recen_TPC[i_group] = (1.0/1.0) * atan2(Qy_recen_TPC[i_group], Qx_recen_TPC[i_group]);
				if(zero2twopi(psi1recen_TPC[i_group])){return 0;}
				psi1shift_TPC[i_group] = psi1recen_TPC[i_group];
			}
		}
		if(ep_setting !="recenter")
		{
			Qy_recen_TPC_ew   /= weight_TPC_ew;
			Qx_recen_TPC_ew   /= weight_TPC_ew;
			psi1recen_TPC[2] = (1.0/1.0) * atan2(Qy_recen_TPC_ew, Qx_recen_TPC_ew);
			if(zero2twopi(psi1recen_TPC[2])){return 0;}
			if(ep_setting == "flatten")
			{
				for(Int_t j=0;j<max_flat;j++)
				{
					for(int k = 0;k<3;k++)
					{
						tpc_flatten_cos[ep_vz_bin][k][j]->Fill(runindex,centrality,cos((1.0*j+1.0)*psi1recen_TPC[k]));
						tpc_flatten_sin[ep_vz_bin][k][j]->Fill(runindex,centrality,sin((1.0*j+1.0)*psi1recen_TPC[k]));
					}
				}
			}
			else
			{
				psi1shift_TPC[2] = psi1recen_TPC[2];
				for(int k = 0;k<3;k++)
				{
					for(Int_t j=0;j<max_flat;j++)
					{
						//see 228.0098-052
						psi1shift_TPC[k]        += (2.0/(1.0*j+1.0))*(-(gettpc_flatten_sin[ep_vz_bin][k][j]   ->GetBinContent(runindex+1,centrality+1))*cos((1.0*j+1.0)*psi1recen_TPC[k]        )+(gettpc_flatten_cos[ep_vz_bin][k][j]   ->GetBinContent(runindex+1,centrality+1))*sin((1.0*j+1.0)*psi1recen_TPC[k]        ));
					}
					if(zero2twopi(psi1shift_TPC[k])){return 0;}
				}
			}
		}
	}
	hist_tpc_A_psi_raw->Fill(psi1rawep_TPC[1]);
	hist_tpc_B_psi_raw->Fill(psi1rawep_TPC[0]);
	hist_tpc_A_psi_recentered->Fill(psi1recen_TPC[1]);
	hist_tpc_B_psi_recentered->Fill(psi1recen_TPC[0]);
	hist_tpc_AB_psi_recentered->Fill(psi1recen_TPC[2]);
	hist_tpc_A_psi_flattened->Fill(psi1shift_TPC[1]);
	hist_tpc_B_psi_flattened->Fill(psi1shift_TPC[0]);	
	hist_tpc_AB_psi_flattened->Fill(psi1shift_TPC[2]);	
    // TPC event plane end
	
	//quick detour to calc pi ratio:
	// int pi_ratio_bin = -1;
	// if( mult_pionminus+mult_pionplus>0)
	// {
		// double pion_ratio =1.0*mult_pionplus/(1.0*mult_pionminus+1.0*mult_pionplus);
		// pip_b_tot_pi_ratio->Fill(pion_ratio);
		// pip_b_tot_pi_cent->Fill(centrality,pion_ratio);
		// pip_b_tot_pi_cent_p->Fill(centrality,pion_ratio);
		// pi_ratio_bin = int(pion_ratio*10);
	// }
	// if(mult_pionplus!=0)
	// {
		// double pi_ratio = 1.0*mult_pionminus/mult_pionplus;
		// hist_pion_yield_ratio[centrality]->Fill(pi_ratio);
		// profile_pion_yield_ratio->Fill(centrality,pi_ratio);
	// }
	// pi_corr_yield->Fill(mult_pionplus,mult_pionminus);
   
   //-----------------Get EPD information------------------------------
    Int_t nepdHits = mPicoDst->numberOfEpdHits();
    StPicoEpdHit *epdHit;
    TVector3 StraightLine_center;
    TVector3 StraightLine_random;
    double phi_epd_center; double phi_epd_random;
    double eta_epd_center; double eta_epd_random;
	double mip; Int_t side = 1;	
    double TileWeight           = {0};
    double Qx_rawep_EPD[5][2]={0.0}, Qy_rawep_EPD[5][2]={0.0};
	double Qx_recenterep_EPD[5][2]={0.0}, Qy_recenterep_EPD[5][2]={0.0};
    double psi1rawep_EPD[5][2]={0.0};
	double psi1recenterep_EPD[5][2]={0.0};
	double psi1flattenep_EPD[5][2]={0.0};
    double weight_EPD[5][2]={0};
	double psi1raw_epd_ab_ew= 0;
	double psi1recentered_epd_ab_ew= 0;
	double psi1flattened_epd_ab_ew = 0;
	double psi1raw_epd_ew = 0;
	double psi1recentered_epd_ew=0;
	double psi1flattened_epd_ew=0;
	
	if(dataset != "27_Gev_2018")
	{
		for(Int_t iHit=0; iHit<nepdHits; iHit++)
		{
			epdHit = mPicoDst->epdHit(iHit);
			mip = epdHit->nMIP();
			int iring = epdHit->row() -1;//(1~16)-1 -> 0-15
			side =1;//east
			if(epdHit->id() > 0 ) {side = 0;}               // unique tile identifier, absolute value is 100*position+tile, sign is +1/-1 for West/East

			int ringgroup = mEpdEpInfo->RingGroup(iring);   //0: 0-7, 1: 8-15    0-> inner most
			if(ringgroup == -1) continue;

			StraightLine_center = mEpdGeom->TileCenter(epdHit->id())        - picoEvent->primaryVertex();
			StraightLine_random = mEpdGeom->RandomPointOnTile(epdHit->id()) - picoEvent->primaryVertex();

			phi_epd_center = StraightLine_center.Phi();
			eta_epd_center = StraightLine_center.Eta();
			phi_epd_random = StraightLine_random.Phi();
			eta_epd_random = StraightLine_random.Eta();

			if(mip < 0.3) continue;
			TileWeight = (mip > 2) ? 2 : mip; //this is in shaowei's code too
			
			//eta weight
			double eta_weight =1.0;
			if(do_eta_weighting && (dataset == "9p2_GeV_2020"|| dataset == "14p6_gev_2019"|| dataset =="17p3_gev_2021" 
			/* ||dataset == "19p6_Gev_2019" */ || dataset == "11p5_gev_2020"))
			{
				if(side){eta_weight=-epd_eta_weight_east->GetBinContent(epd_eta_weight_east->FindBin(eta_epd_random));}//the negative sign is because I forgot to flip sign for calculating weight
				else{eta_weight=-epd_eta_weight_west->GetBinContent(epd_eta_weight_west->FindBin(eta_epd_random));}
			}
			if(do_eta_weighting && dataset == "19p6_Gev_2019")
			{
				if(side){eta_weight=-epd_eta_weight_east_cent[centrality][ep_vz_bin]->GetBinContent(epd_eta_weight_east_cent[centrality][ep_vz_bin]->FindBin(eta_epd_random));}//the negative sign is because I forgot to flip sign for calculating weight
				else{eta_weight=-epd_eta_weight_west_cent[centrality][ep_vz_bin]->GetBinContent(epd_eta_weight_west_cent[centrality][ep_vz_bin]->FindBin(eta_epd_random));}
			}
			TileWeight *= eta_weight;
			for(int i_group = 0;i_group<5;i_group++)
			{
				if(i_group == 0 && !(ringgroup == 0 || ringgroup == 1)){continue;}// 0: EPD-AB
				if(i_group == 1 && !(ringgroup == 2)                  ){continue;}// 1: EPD-C
				if(i_group == 2 && !(ringgroup == 3)                  ){continue;}// 2: EPD-D
				//                                                                // 3: EPD-ABCD
				if(i_group == 4 && !(ringgroup == 2 || ringgroup == 3)){continue;}// 4: EPD-CD
				Qx_rawep_EPD[i_group][side] += TileWeight * cos(1.0*phi_epd_random);
				Qy_rawep_EPD[i_group][side] += TileWeight * sin(1.0*phi_epd_random);
				if(ep_setting !="recenter")
				{
					Qx_recenterep_EPD[i_group][side] += TileWeight * cos(1.0*phi_epd_random) - getepd_recenterx[ep_vz_bin][i_group][side]->GetBinContent(getepd_recenterx[ep_vz_bin][i_group][side]->FindBin(runindex,centrality));
					Qy_recenterep_EPD[i_group][side] += TileWeight * sin(1.0*phi_epd_random) - getepd_recentery[ep_vz_bin][i_group][side]->GetBinContent(getepd_recentery[ep_vz_bin][i_group][side]->FindBin(runindex,centrality));
				}
				else
				{
					epd_recenterx[ep_vz_bin][i_group][side]->Fill(runindex,centrality,TileWeight * cos(1.0*phi_epd_random));
					epd_recentery[ep_vz_bin][i_group][side]->Fill(runindex,centrality,TileWeight * sin(1.0*phi_epd_random));
				}
				weight_EPD[i_group][side] += TileWeight;
			}
		}
		double Qx_raw_epd_ab=Qx_rawep_EPD[0][1]-Qx_rawep_EPD[0][0];
		double Qy_raw_epd_ab=Qy_rawep_EPD[0][1]-Qy_rawep_EPD[0][0];
		double Qx_recentered_epd_ab=Qx_recenterep_EPD[0][1]-Qx_recenterep_EPD[0][0];
		double Qy_recentered_epd_ab=Qy_recenterep_EPD[0][1]-Qy_recenterep_EPD[0][0];
		double weight_epd_ab=weight_EPD[0][0]+weight_EPD[0][1];
		
		//full
		
		double Qx_raw_epd_full=Qx_rawep_EPD[3][1]-Qx_rawep_EPD[3][0];
		double Qy_raw_epd_full=Qy_rawep_EPD[3][1]-Qy_rawep_EPD[3][0];
		double Qx_recentered_epd_full=Qx_recenterep_EPD[3][1]-Qx_recenterep_EPD[3][0];
		double Qy_recentered_epd_full=Qy_recenterep_EPD[3][1]-Qy_recenterep_EPD[3][0];
		double weight_epd_full=weight_EPD[3][0]+weight_EPD[3][1];
		for(int i=0; i<5; i++)
		{
			for(int j=0;j<2;j++)
			{
				if(dataset == "4p5_gev_FXT_2020" || dataset=="3p85_Gev_2018")
				{
					if(weight_EPD[i][0] == 0 && weight_EPD[i][1] == 0) {return 0;}
				}
				else
				{
					if(Qx_rawep_EPD[i][j] ==0. || Qy_rawep_EPD[i][j] == 0. || weight_EPD[i][j] == 0.) {return 0;}
					if(ep_setting !="recenter" &&(Qx_recenterep_EPD[i][j] ==0. || Qy_recenterep_EPD[i][j] ==0.)) {return 0;}
				}
				hist_Qr_epd[i][j]->Fill(centrality,TMath::Sqrt(Qy_recenterep_EPD[i][j]*Qy_recenterep_EPD[i][j]+Qx_recenterep_EPD[i][j]*Qx_recenterep_EPD[i][j]));
				Qx_rawep_EPD[i][j] /= weight_EPD[i][j];
				Qy_rawep_EPD[i][j] /= weight_EPD[i][j];
				Qx_recenterep_EPD[i][j] /= weight_EPD[i][j];
				Qy_recenterep_EPD[i][j] /= weight_EPD[i][j];
				
				psi1rawep_EPD[i][j] = atan2(Qy_rawep_EPD[i][j],Qx_rawep_EPD[i][j])/1.0;
				if(zero2twopi(psi1rawep_EPD[i][j])){return 0;}
				if(ep_setting !="recenter")
				{
					psi1recenterep_EPD[i][j] = atan2(Qy_recenterep_EPD[i][j],Qx_recenterep_EPD[i][j])/1.0;
					if(zero2twopi(psi1recenterep_EPD[i][j])){return 0;}
					psi1flattenep_EPD[i][j]  = psi1recenterep_EPD[i][j];
				}
			}
		}
		if(ep_setting !="recenter")
		{
			Qx_raw_epd_ab /= weight_epd_ab;
			Qy_raw_epd_ab /= weight_epd_ab;
			Qx_recentered_epd_ab /= weight_epd_ab;
			Qy_recentered_epd_ab /= weight_epd_ab;
			psi1raw_epd_ab_ew=atan2(Qy_raw_epd_ab,Qx_raw_epd_ab);
			psi1recentered_epd_ab_ew=atan2(Qy_recentered_epd_ab,Qx_recentered_epd_ab);
			if(zero2twopi(psi1raw_epd_ab_ew)){return 0;}
			if(zero2twopi(psi1recentered_epd_ab_ew)){return 0;}
			psi1flattened_epd_ab_ew = psi1recentered_epd_ab_ew;
			
			Qx_raw_epd_full /= weight_epd_full;
			Qy_raw_epd_full /= weight_epd_full;
			Qx_recentered_epd_full /= weight_epd_full;
			Qy_recentered_epd_full /= weight_epd_full;
			psi1raw_epd_ew=atan2(Qy_raw_epd_full,Qx_raw_epd_full);
			psi1recentered_epd_ew=atan2(Qy_recentered_epd_full,Qx_recentered_epd_full);
			if(zero2twopi(psi1raw_epd_ew)){return 0;}
			if(zero2twopi(psi1recentered_epd_ew)){return 0;}
			psi1flattened_epd_ew = psi1recentered_epd_ew;
			
			if(ep_setting =="flatten")
			{
				for(Int_t j=0;j<max_flat;j++)
				{
					for(Int_t i=0;i<5;i++)
					{
						for(Int_t k=0;k<2;k++)
						{
							epd_flatten_cos[ep_vz_bin][i][j][k]->Fill(runindex,centrality,cos((1.0*j+1.0)*psi1recenterep_EPD[i][k]));
							epd_flatten_sin[ep_vz_bin][i][j][k]->Fill(runindex,centrality,sin((1.0*j+1.0)*psi1recenterep_EPD[i][k]));
						}
					}
					epd_flatten_ab_ew_cos[ep_vz_bin][j]->Fill(runindex,centrality,cos((1.0*j+1.0)*psi1recentered_epd_ab_ew));
					epd_flatten_ab_ew_sin[ep_vz_bin][j]->Fill(runindex,centrality,sin((1.0*j+1.0)*psi1recentered_epd_ab_ew));
					
					epd_flatten_ew_cos[ep_vz_bin][j]->Fill(runindex,centrality,cos((1.0*j+1.0)*psi1recentered_epd_ew));
					epd_flatten_ew_sin[ep_vz_bin][j]->Fill(runindex,centrality,sin((1.0*j+1.0)*psi1recentered_epd_ew));
				}
			}
			else
			{
				for(int j=0; j<max_flat; j++)
				{
					//see 228.0098-052, 228.0099-032
					psi1flattened_epd_ab_ew  += (2.0/(1.0*j+1.0))*(-(getepd_flatten_ab_ew_sin[ep_vz_bin][j]->GetBinContent(runindex+1,centrality+1))*cos((1.0*j+1.0)*psi1recentered_epd_ab_ew)+(getepd_flatten_ab_ew_cos[ep_vz_bin][j]->GetBinContent(runindex+1,centrality+1))*sin((1.0*j+1.0)*psi1recentered_epd_ab_ew));
					psi1flattened_epd_ew  += (2.0/(1.0*j+1.0))*(-(getepd_flatten_ew_sin[ep_vz_bin][j]->GetBinContent(runindex+1,centrality+1))*cos((1.0*j+1.0)*psi1recentered_epd_ew)+(getepd_flatten_ew_cos[ep_vz_bin][j]->GetBinContent(runindex+1,centrality+1))*sin((1.0*j+1.0)*psi1recentered_epd_ew));
					for(int i = 0; i<5;i++)
					{
						for(int k=0;k<2;k++)
						{
							psi1flattenep_EPD[i][k] += (2.0/(1.0*j+1.0))*(-(getepd_flatten_sin[ep_vz_bin][i][j][k]->GetBinContent(runindex+1,centrality+1))*cos((1.0*j+1.0)*psi1recenterep_EPD[i][k])+(getepd_flatten_cos[ep_vz_bin][i][j][k]->GetBinContent(runindex+1,centrality+1))*sin((1.0*j+1.0)*psi1recenterep_EPD[i][k]));
						}
					}
				}
				for(int i=0; i<5; i++){for(int k=0;k<2;k++){if(zero2twopi(psi1flattenep_EPD[i][k])){return 0;}}}
				if(zero2twopi(psi1flattened_epd_ab_ew)){return 0;}
				if(zero2twopi(psi1flattened_epd_ew)){return 0;}
			}
		}
	}
	else if(dataset == "27_Gev_2018")
	{
		int ep_num = EP_group(run_number,centrality);
		StEpdEpInfo EpdResult = mEpFinder->Results(mPicoDst->picoArray(StPicoArrays::EpdHit),pVertex,ep_num);
		psi1rawep_EPD[3][0] = EpdResult.WestRawPsi(1);
		psi1rawep_EPD[3][1]= EpdResult.EastRawPsi(1);
		
		psi1recentered_epd_ew = EpdResult.FullPhiWeightedPsi(1);
		psi1recenterep_EPD[3][0] = EpdResult.WestPhiWeightedPsi(1);
		psi1recenterep_EPD[3][1] = EpdResult.EastPhiWeightedPsi(1);

		psi1flattened_epd_ew = EpdResult.FullPhiWeightedAndShiftedPsi(1); // It returns Psi1(EP angle) from 0 to +2pi;
		psi1flattenep_EPD[3][0] = EpdResult.WestPhiWeightedAndShiftedPsi(1);
		psi1flattenep_EPD[3][1] = EpdResult.EastPhiWeightedAndShiftedPsi(1);
	}
	
	hist_epd_ab_psi_raw_ew->Fill(psi1raw_epd_ab_ew);
	hist_epd_ab_psi_recentered_ew->Fill(psi1recentered_epd_ab_ew);
	hist_epd_ab_psi_flattened_ew->Fill(psi1flattened_epd_ab_ew);
	hist_epd_psi_raw_ew->Fill(psi1raw_epd_ew);
	hist_epd_psi_recentered_ew->Fill(psi1recentered_epd_ew);
	hist_epd_psi_flattened_ew->Fill(psi1flattened_epd_ew);
	
	for(int j=0;j<2;j++)
	{
		hist_epd_ab_psi_raw[j]->Fill(psi1rawep_EPD[0][j]);
		hist_epd_c_psi_raw[j]->Fill(psi1rawep_EPD[1][j]);
		hist_epd_d_psi_raw[j]->Fill(psi1rawep_EPD[2][j]);
		hist_epd_abcd_psi_raw[j]->Fill(psi1rawep_EPD[3][j]);
		hist_epd_cd_psi_raw[j]->Fill(psi1rawep_EPD[4][j]);
		
		hist_epd_ab_psi_recentered[j]->Fill(psi1recenterep_EPD[0][j]);
		hist_epd_c_psi_recentered[j]->Fill(psi1recenterep_EPD[1][j]);
		hist_epd_d_psi_recentered[j]->Fill(psi1recenterep_EPD[2][j]);
		hist_epd_abcd_psi_recentered[j]->Fill(psi1recenterep_EPD[3][j]);
		hist_epd_cd_psi_recentered[j]->Fill(psi1recenterep_EPD[4][j]);
		
		hist_epd_ab_psi_flattened[j]->Fill(psi1flattenep_EPD[0][j]);
		hist_epd_c_psi_flattened[j]->Fill(psi1flattenep_EPD[1][j]);
		hist_epd_d_psi_flattened[j]->Fill(psi1flattenep_EPD[2][j]);
		hist_epd_abcd_psi_flattened[j]->Fill(psi1flattenep_EPD[3][j]);
		hist_epd_cd_psi_flattened[j]->Fill(psi1flattenep_EPD[4][j]);
	}
	
	//////////////////////////////////////eta weighting
	for(Int_t iHit=0; iHit<nepdHits; iHit++)
	{
        epdHit = mPicoDst->epdHit(iHit);
        mip = epdHit->nMIP();
        int iring = epdHit->row() -1;//(1~16)-1 -> 0-15
		side =1;//east
        if(epdHit->id() > 0 ) {side = 0;}               // unique tile identifier, absolute value is 100*position+tile, sign is +1/-1 for West/East

        int ringgroup = mEpdEpInfo->RingGroup(iring);   //0: 0-7, 1: 8-15    0-> inner most
        if(ringgroup == -1) continue;

		StraightLine_random = mEpdGeom->RandomPointOnTile(epdHit->id()) - picoEvent->primaryVertex();
		phi_epd_random = StraightLine_random.Phi();
        eta_epd_random = StraightLine_random.Eta();

        if(mip < 0.3) continue;
        TileWeight = (mip > 2) ? 2 : mip; //this is in shaowei's code too

        // if(ringgroup == 0 || ringgroup == 1) // 0: EPD-AB
		// if(side==0){v1_vs_eta_west[centrality]->Fill(eta_epd_random,TileWeight*cos(phi_epd_random-psi1flattenep_EPD[0][1]));}
		// else{v1_vs_eta_east[centrality]->Fill(eta_epd_random,TileWeight*cos(phi_epd_random-psi1flattenep_EPD[0][0]));}
		// if(side==0){v1_vs_eta_west[centrality]->Fill(eta_epd_random,cos(phi_epd_random-psi1flattenep_EPD[0][1]),TileWeight);}
		// else{v1_vs_eta_east[centrality]->Fill(eta_epd_random,cos(phi_epd_random-psi1flattenep_EPD[0][0]),TileWeight);}
		if(side==0){v1_vs_eta_west[centrality][ep_vz_bin]->Fill(eta_epd_random,cos(phi_epd_random-psi1flattenep_EPD[3][1]),TileWeight);}
		else{v1_vs_eta_east[centrality][ep_vz_bin]->Fill(eta_epd_random,cos(phi_epd_random-psi1flattenep_EPD[3][0]),TileWeight);}
	}
	////////////////////////////////////////end eta weighting
	if(ep_setting !="flow"){return kStOK;}
	//correlation:
	profile_correlation_abc->Fill(centrality,cos(psi1flattenep_EPD[0][1]-psi1flattenep_EPD[1][1]));
	profile_correlation_abd->Fill(centrality,cos(psi1flattenep_EPD[0][1]-psi1flattenep_EPD[2][1]));
	profile_correlation_cd->Fill(centrality,cos(psi1flattenep_EPD[1][1]-psi1flattenep_EPD[2][1]));
	
	profile_correlation_eABtB->Fill(centrality,cos(psi1flattenep_EPD[0][1]-psi1shift_TPC[0]));
	profile_correlation_eCtB->Fill(centrality,cos(psi1flattenep_EPD[1][1]-psi1shift_TPC[0]));
	profile_correlation_eABtB->Fill(centrality,cos(psi1flattenep_EPD[2][1]-psi1shift_TPC[0]));
	
	profile_correlation_eABtA->Fill(centrality,cos(psi1flattenep_EPD[0][1]-psi1shift_TPC[1]));
	profile_correlation_eCtA->Fill(centrality,cos(psi1flattenep_EPD[1][1]-psi1shift_TPC[1]));
	profile_correlation_eDtA->Fill(centrality,cos(psi1flattenep_EPD[2][1]-psi1shift_TPC[1]));
	
	profile_correlation_tpc_AB->Fill(centrality,cos(psi1shift_TPC[1]-psi1shift_TPC[0]));
	
	profile_correlation_abew->Fill(centrality,cos(psi1flattenep_EPD[0][1]-psi1flattenep_EPD[0][0]));
	profile_correlation_cew->Fill(centrality,cos(psi1flattenep_EPD[1][1]-psi1flattenep_EPD[1][0]));
	profile_correlation_dew->Fill(centrality,cos(psi1flattenep_EPD[2][1]-psi1flattenep_EPD[2][0]));
	
	profile_correlation_tpc_ew_vz[ep_vz_bin]->Fill(centrality,cos(psi1shift_TPC[1]-psi1shift_TPC[0]));
	profile_correlation_ew_vz[ep_vz_bin]->Fill(centrality,cos(psi1flattenep_EPD[3][1]-psi1flattenep_EPD[3][0]));
	profile_correlation_eABw_vz[ep_vz_bin]->Fill(centrality,cos(psi1flattenep_EPD[0][1]-psi1flattenep_EPD[3][0]));
	profile_correlation_eCDw_vz[ep_vz_bin]->Fill(centrality,cos(psi1flattenep_EPD[4][1]-psi1flattenep_EPD[3][0]));
	profile_correlation_eABeCD_vz[ep_vz_bin]->Fill(centrality,cos(psi1flattenep_EPD[0][1]-psi1flattenep_EPD[4][1]));
	profile_correlation_ewAB_vz[ep_vz_bin]->Fill(centrality,cos(psi1flattenep_EPD[3][1]-psi1flattenep_EPD[0][0]));
	profile_correlation_ewCD_vz[ep_vz_bin]->Fill(centrality,cos(psi1flattenep_EPD[3][1]-psi1flattenep_EPD[4][0]));
	profile_correlation_wABwCD_vz[ep_vz_bin]->Fill(centrality,cos(psi1flattenep_EPD[0][0]-psi1flattenep_EPD[4][0]));
	
	profile_correlation_eDw_vz[ep_vz_bin]->Fill(centrality,cos(psi1flattenep_EPD[2][1]-psi1flattenep_EPD[3][0]));
	profile_correlation_eABeD_vz[ep_vz_bin]->Fill(centrality,cos(psi1flattenep_EPD[0][1]-psi1flattenep_EPD[2][1]));
	profile_correlation_ewD_vz[ep_vz_bin]->Fill(centrality,cos(psi1flattenep_EPD[3][1]-psi1flattenep_EPD[2][0]));
	profile_correlation_wABwD_vz[ep_vz_bin]->Fill(centrality,cos(psi1flattenep_EPD[0][0]-psi1flattenep_EPD[2][0]));
	profile_correlation_TPCEPDe_vz[ep_vz_bin]->Fill(centrality,cos(psi1shift_TPC[2]-psi1flattenep_EPD[3][1]));
	profile_correlation_TPCEPDw_vz[ep_vz_bin]->Fill(centrality,cos(psi1shift_TPC[2]-psi1flattenep_EPD[3][0]));
	
	profile_correlation_ew->Fill(centrality,cos(psi1flattenep_EPD[3][1]-psi1flattenep_EPD[3][0]));
	
	correlation_abc->Fill(psi1flattenep_EPD[0][1],psi1flattenep_EPD[1][1]);
	correlation_abd->Fill(psi1flattenep_EPD[0][1],psi1flattenep_EPD[2][1]);
	correlation_cd->Fill(psi1flattenep_EPD[1][1],psi1flattenep_EPD[2][1]);
	
	correlation_eABtB->Fill(psi1flattenep_EPD[0][1],psi1shift_TPC[0]);
	correlation_eCtB->Fill(psi1flattenep_EPD[1][1],psi1shift_TPC[0]);
	correlation_eDtB->Fill(psi1flattenep_EPD[2][1],psi1shift_TPC[0]);
	
	correlation_eABtA->Fill(psi1flattenep_EPD[0][1],psi1shift_TPC[1]);
	correlation_eCtA->Fill(psi1flattenep_EPD[1][1],psi1shift_TPC[1]);
	correlation_eDtA->Fill(psi1flattenep_EPD[2][1],psi1shift_TPC[1]);
	
	correlation_tpc_AB->Fill(psi1shift_TPC[1],psi1shift_TPC[0]);
	
	correlation_abew->Fill(psi1flattenep_EPD[0][1],psi1flattenep_EPD[0][0]);
	correlation_cew->Fill(psi1flattenep_EPD[1][1],psi1flattenep_EPD[1][0]);
	correlation_dew->Fill(psi1flattenep_EPD[2][1],psi1flattenep_EPD[2][0]);
	
	correlation_ew->Fill(psi1flattenep_EPD[3][1],psi1flattenep_EPD[3][0]);
	
	correlation_ew_cent[centrality]->Fill(psi1flattenep_EPD[3][1],psi1flattenep_EPD[3][0]);
	
	////////////////////////////////////////////////////PID loop
	for (Int_t itr=0;itr<nTrack;itr++)
	{
        const StPicoTrack *ptrk = (StPicoTrack*)mPicoDst->track(itr);
        if(!ptrk)  continue;
        if(!ptrk->isPrimary())  continue;  // now selecting primary tracks
        if(!isGoodTrack(ptrk))  continue;
		// Get PID parameters
		double pz = ptrk->pMom().Z();
		double pt = ptrk->pPt();
		double trackP = ptrk->pPtot();
		double eta = ptrk->pMom().Eta();
		double phi = ptrk->pMom().Phi();if(zero2twopi(phi)){continue;}
		int charge = ptrk->charge();
		double v1 =cos(1.0*(phi-psi1flattened_epd_ew));
		if(dataset == "3p85_Gev_2018"){v1 =cos(1.0*(phi-psi1flattenep_EPD[0][1]));}// EPD AB east plane
		if(!(v1>-1.0&&v1<1.0)){continue;}
		double max_pt_cut = 2.0;
		double max_ptot_cut = 3.0; if(dataset == "19p6_Gev_2019"){max_ptot_cut = 2.8;}
		// double max_eta_cut = 1.6; 
		double max_eta_cut = 2.0;
		if(dataset == "27_Gev_2018"|| dataset == "19p6_Gev_2019"){max_eta_cut = 1.0;}
		if(dataset == "3p85_Gev_2018")
		{
			max_pt_cut = 2.0;
			max_ptot_cut = 1000000.0;
			max_eta_cut = 10;
		}
		double proton_eff_weight = proton_efficiency_correction(ptrk,centrality);
		//track histograms:
		hist_dEdx->Fill(charge*trackP,ptrk->dEdx());
		ch_v1_eta_cent[centrality]->Fill(eta,v1);
		//purity:
		if(charge>0){hist_nsigma_proton->Fill(ptrk->nSigmaProton());}
		if(charge<0){hist_nsigma_aproton->Fill(ptrk->nSigmaProton());}
		
		//look for tof info
		double Beta=0; double mass2 = 0;
		double btofMatchFlag =0; double btofYLocal = 0;
		double etofMatchFlag =0;
		tof_efficiency_tot->Fill(eta,pt);

		if(ptrk->isTofTrack())
		{
			tof_efficiency_match->Fill(eta,pt);
			StPicoBTofPidTraits *trait = mPicoDst->btofPidTraits(ptrk->bTofPidTraitsIndex());
			if(trait)
			{
				Beta = trait->btofBeta();
				mass2 = trackP*trackP*((1.0/(Beta*Beta))-1.0);
				hist_tof_pid->Fill(charge*trackP,1.0/Beta);
				hist_mass2_pid->Fill(charge*trackP,mass2);
				btofMatchFlag = trait->btofMatchFlag();
				btofYLocal    = trait->btofYLocal();
			}
		}
		if(ptrk->isETofTrack() && include_etof)
		{
			tof_efficiency_match->Fill(eta,pt);
			StPicoETofPidTraits *etrait = mPicoDst->etofPidTraits(ptrk->eTofPidTraitsIndex());
			if(etrait)
			{
				Beta = etrait->beta();
				mass2 = trackP*trackP*((1.0/(Beta*Beta))-1.0);
				hist_tof_pid->Fill(charge*trackP,1.0/Beta);
				hist_mass2_pid->Fill(charge*trackP,mass2);
				etofMatchFlag = etrait->matchFlag();
			}
		}
		//low purity method:
		if(isProton(ptrk,false,0) && ptrk->isTofTrack())
		{
			if(charge>0){hist_mass2_proton->Fill(mass2,trackP);}
			if(charge<0)
			{
				hist_mass2_aproton->Fill(mass2,trackP);
				double energy_P = sqrt(trackP*trackP + protonMass*protonMass);
				double rap_P = -0.5*log( (energy_P + pz) / (energy_P - pz) )-y_cm;
				if(fabs(rap_P)<1.0 && pt>0.4 && pt<3.0 && fabs(eta)<2.0 &&trackP<3.0)
				{
					int y_bin = abs((rap_P+1.0)*5.0);
					double phi_proton = phi-psi1flattened_epd_ew;
					if(phi_proton<-1.0*pi){phi_proton+=2.0*pi;}
					if(phi_proton>1.0*pi){phi_proton-=2.0*pi;}
					int phi_bin = abs(10*phi_proton/pi);
					int p_bin = abs(2*trackP);
					if(3<centrality&&centrality<7){h_protonminus_y_cent[0][0][y_bin][phi_bin][p_bin]->Fill(mass2,proton_eff_weight);}
				}
			}
		}
		// high purity method:
		if(isProton(ptrk,true,0))
		{
			double energy_Proton = sqrt(trackP*trackP + protonMass*protonMass);
			double rap_Proton = -0.5*log( (energy_Proton + pz) / (energy_Proton - pz) )-y_cm;
			// double resolution[9] ={1};
			// if(dataset == "3p85_Gev_2018")
			// {
				// resolution[0]=0.10384;resolution[1]=0.242532;resolution[2]=0.44511;resolution[3]=0.601817;
				// resolution[4]=0.686444;resolution[5]=0.701874;resolution[6]=0.624367;resolution[7]=0.463533;
				// resolution[8]=0.270264;
			// }
			
			//proton generic histograms:
			hist_dEdx_proton->Fill(charge*trackP,ptrk->dEdx());
			hist_pt_centrality->Fill(charge*pt,centrality);
			hist_pt_eta->Fill(eta,charge*pt);
			if(ptrk->isTofTrack()||ptrk->isETofTrack())
			{
				hist_tof_pid_proton->Fill(trackP/charge,1.0/Beta);
				hist_mass2_pid_proton->Fill(trackP/charge,mass2);
			}
			
			if(charge>0)
			{
				baryon_yield_histogram->Fill(TMath::Abs(rap_Proton),1.0);
				proton_yield[centrality]->Fill(rap_Proton);
				hist_dca_vs_p_proton->Fill(trackP,ptrk->gDCA(picoEvent->primaryVertex()).Mag());
				hist_nsigma_p_proton->Fill(ptrk->nSigmaProton(),trackP);
				p_nsigma_proton->Fill(trackP,ptrk->nSigmaProton());
				
				if(pt>0.4 && pt<max_pt_cut && fabs(eta)<max_eta_cut && trackP<max_ptot_cut)
				{
					p_protonplusv1_y_cent[ep_vz_bin][0][centrality]->Fill(rap_Proton,v1,proton_eff_weight);/////////////////HERE!!
					hist_pt_y_proton->Fill(rap_Proton,pt);
					hist_pt_y_proton_cent[centrality]->Fill(rap_Proton,pt);
					eta_v_pt_pos_cent[centrality]->Fill(eta,pt);
					// p_protonplusv1_yn_cent[centrality]->Fill(rap_Proton/y_beam,v1);
					// p_protonplusv1_eta_cent[centrality]->Fill(eta,v1);
					// p_protonplusv1_etan_cent[centrality]->Fill(eta/y_beam,v1);
					
					if(ptrk->isTofTrack()||ptrk->isETofTrack()){hist_proton_accuracy->Fill(ptrk->nSigmaProton(),mass2);}
					// if(pi_ratio_bin!=-1)
					// {
						// p_protonplusv1_y_pip[pi_ratio_bin]->Fill(rap_Proton,v1/resolution[centrality],proton_eff_weight);
						// if(3<centrality&&centrality<7){p_protonplusv1_y_mid_pip[pi_ratio_bin]->Fill(rap_Proton,v1/resolution[centrality],proton_eff_weight);}
					// }
					if(vz<-70)
					{
						hist_pt_y_proton_vz[0]->Fill(rap_Proton,pt);
						if(3<centrality&&centrality<7){p_protonplusv1_y_mid_vz[0]->Fill(rap_Proton,v1,proton_eff_weight);}
					}
					if(vz>-70.0&&vz<0.0)
					{
						hist_pt_y_proton_vz[1]->Fill(rap_Proton,pt);
						if(3<centrality&&centrality<7){p_protonplusv1_y_mid_vz[1]->Fill(rap_Proton,v1,proton_eff_weight);}
					}
					if(vz<70.0&&vz>0.0)
					{
						hist_pt_y_proton_vz[2]->Fill(rap_Proton,pt);
						if(3<centrality&&centrality<7){p_protonplusv1_y_mid_vz[2]->Fill(rap_Proton,v1,proton_eff_weight);}
					}
					if(vz>70.0)
					{
						hist_pt_y_proton_vz[3]->Fill(rap_Proton,pt);
						if(3<centrality&&centrality<7){p_protonplusv1_y_mid_vz[3]->Fill(rap_Proton,v1,proton_eff_weight);}
					}
					if(eta < 0.0){p_protonplusv1_pt_cent[centrality]->Fill(pt,v1,proton_eff_weight);}
					else if(eta > 0.0){p_protonplusv1_pt_cent[centrality]->Fill(pt,-v1,proton_eff_weight);}
				}
			}
			else
			{
				baryon_yield_histogram->Fill(TMath::Abs(rap_Proton),-1.0);
				antiproton_yield[centrality]->Fill(rap_Proton);
				hist_dca_vs_p_aproton->Fill(trackP,ptrk->gDCA(picoEvent->primaryVertex()).Mag());
				hist_nsigma_p_aproton->Fill(ptrk->nSigmaProton(),trackP);
				p_nsigma_aproton->Fill(trackP,ptrk->nSigmaProton());
				
				if(pt>0.4 && pt<max_pt_cut && fabs(eta)<max_eta_cut && trackP<max_ptot_cut)
				{
					p_protonminusv1_y_cent[ep_vz_bin][0][centrality]->Fill(rap_Proton,v1,proton_eff_weight);
					hist_pt_y_antiproton->Fill(rap_Proton,pt);
					hist_pt_y_antiproton_cent[centrality]->Fill(rap_Proton,pt);
					eta_v_pt_neg_cent[centrality]->Fill(eta,pt);
					// p_protonminusv1_yn_cent[centrality]->Fill(rap_Proton/y_beam,v1);
					// p_protonminusv1_eta_cent[centrality]->Fill(eta,v1);
					// p_protonminusv1_etan_cent[centrality]->Fill(eta/y_beam,v1);
					
					if(ptrk->isTofTrack()||ptrk->isETofTrack()){hist_aproton_accuracy->Fill(ptrk->nSigmaProton(),mass2);}
					// if(pi_ratio_bin!=-1)
					// {
						// p_protonminusv1_y_pip[pi_ratio_bin]->Fill(rap_Proton,v1/resolution[centrality],proton_eff_weight);
						// if(3<centrality&&centrality<7){p_protonminusv1_y_mid_pip[pi_ratio_bin]->Fill(rap_Proton,v1/resolution[centrality],proton_eff_weight);}
					// }
					if(vz<-70)
					{
						hist_pt_y_antiproton_vz[0]->Fill(rap_Proton,pt);
						if(3<centrality&&centrality<7){p_protonminusv1_y_mid_vz[0]->Fill(rap_Proton,v1,proton_eff_weight);}
					}
					if(vz>-70.0&&vz<0.0)
					{
						hist_pt_y_antiproton_vz[1]->Fill(rap_Proton,pt);
						if(3<centrality&&centrality<7){p_protonminusv1_y_mid_vz[1]->Fill(rap_Proton,v1,proton_eff_weight);}
					}
					if(vz<70.0&&vz>0.0)
					{
						hist_pt_y_antiproton_vz[2]->Fill(rap_Proton,pt);
						if(3<centrality&&centrality<7){p_protonminusv1_y_mid_vz[2]->Fill(rap_Proton,v1,proton_eff_weight);}
					}
					if(vz>70.0)
					{
						hist_pt_y_antiproton_vz[3]->Fill(rap_Proton,pt);
						if(3<centrality&&centrality<7){p_protonminusv1_y_mid_vz[3]->Fill(rap_Proton,v1,proton_eff_weight);}
					}
					if(eta > -2.0 && eta < 0.0){p_protonminusv1_pt_cent[centrality]->Fill(pt,v1,proton_eff_weight);}
					else if(eta > 0.0 && eta < 2.0){p_protonminusv1_pt_cent[centrality]->Fill(pt,-v1,proton_eff_weight);}
				}
			}
		}
		//////////////////////////////////////////////////////////////////systematics:
		for(int h = 1;h<9;h++)
		{
			if(isProton(ptrk,false,h) && ptrk->isTofTrack()||ptrk->isETofTrack())
			{
				if(charge<0)
				{
					double energy_P = sqrt(trackP*trackP + protonMass*protonMass);
					double rap_P = -0.5*log( (energy_P + pz) / (energy_P - pz) )-y_cm;
					if(fabs(rap_P)<1.0 && pt>0.4 && pt<3.0 && fabs(eta)<2.0 &&trackP<3.0)
					{
						int y_bin = abs((rap_P+1.0)*5.0);
						double phi_proton = phi-psi1flattened_epd_ew;
						if(phi_proton<-1.0*pi){phi_proton+=2.0*pi;}
						if(phi_proton>1.0*pi){phi_proton-=2.0*pi;}
						int phi_bin = abs(10*phi_proton/pi);
						int p_bin = abs(2*trackP);
						if(3<centrality&&centrality<7){h_protonminus_y_cent[h][0][y_bin][phi_bin][p_bin]->Fill(mass2,proton_eff_weight);}
					}
				}
			}
			if(isProton(ptrk,true,h))
			{
				double energy_Proton = sqrt(trackP*trackP + protonMass*protonMass);
				double rap_Proton = -0.5*log( (energy_Proton + pz) / (energy_Proton - pz) )-y_cm;
				if(charge>0)
				{
					if(pt>0.4 && pt<max_pt_cut && fabs(eta)<max_eta_cut && trackP<max_ptot_cut)
					{
						p_protonplusv1_y_cent[ep_vz_bin][h][centrality]->Fill(rap_Proton,v1,proton_eff_weight);
					}
				}
				else
				{
					if(pt>0.4 && pt<max_pt_cut && fabs(eta)<max_eta_cut && trackP<max_ptot_cut)
					{
						p_protonminusv1_y_cent[ep_vz_bin][h][centrality]->Fill(rap_Proton,v1,proton_eff_weight);
					}
				}
			}
		}
		if(isPion(ptrk))
		{
			double energy_pion = sqrt(trackP*trackP + pionMass*pionMass);
			double rap_pion = -0.5*log( (energy_pion + pz) / (energy_pion - pz) )-y_cm;
			hist_dEdx_pion->Fill(charge*trackP,ptrk->dEdx());
			if(ptrk->isTofTrack()||ptrk->isETofTrack())
			{
				hist_tof_pid_pion->Fill(trackP/charge,1.0/Beta);
				hist_mass2_pid_pion->Fill(trackP/charge,mass2);
			}
			if(charge>0)
			{
				p_pionplusv1_y_cent[centrality]->Fill(rap_pion,v1);
				pip_yield[centrality]->Fill(rap_pion);
				p_nsigma_pip->Fill(trackP,ptrk->nSigmaPion());
				hist_nsigma_p_pip->Fill(ptrk->nSigmaPion(),trackP);
				hist_pt_y_pip->Fill(rap_pion,pt);
			}
			if(charge<0)
			{
				p_pionminusv1_y_cent[centrality]->Fill(rap_pion,v1);
				pim_yield[centrality]->Fill(rap_pion);
				p_nsigma_pim->Fill(trackP,ptrk->nSigmaPion());
				hist_nsigma_p_pim->Fill(ptrk->nSigmaPion(),trackP);
				hist_pt_y_pim->Fill(rap_pion,pt);
			}
		}	
		if(isKaon(ptrk))
		{
			double energy_kaon = sqrt(trackP*trackP + kaonMass*kaonMass);
			double rap_kaon = -0.5*log( (energy_kaon + pz) / (energy_kaon - pz) )-y_cm;
			hist_dEdx_kaon->Fill(trackP/charge,ptrk->dEdx());
			if(ptrk->isTofTrack()||ptrk->isETofTrack())
			{
				hist_tof_pid_kaon->Fill(trackP/charge,1.0/Beta);
				hist_mass2_pid_kaon->Fill(trackP/charge,mass2);
			}
			if(charge>0)
			{
				p_kaonplusv1_y_cent[centrality]->Fill(rap_kaon,v1);
				hist_pt_y_kp->Fill(rap_kaon,pt);
			}
			if(charge<0)
			{
				p_kaonminusv1_y_cent[centrality]->Fill(rap_kaon,v1);
				hist_pt_y_km->Fill(rap_kaon,pt);
			}
		}	
	}
	// QA_tree->Fill();
    return kStOK;
}
//end loop

//__________________________________________________________________________________
bool Shift::isGoodEvent(const StPicoEvent *event)
{
    Float_t vx=event->primaryVertex().X();
    Float_t vy=event->primaryVertex().Y();
    Float_t vz=event->primaryVertex().Z();
	if(dataset == "9p2_GeV_2020"||dataset == "7p7_Gev_2021"||dataset == "11p5_gev_2020"|| 
	dataset =="14p6_gev_2019"|| dataset == "17p3_gev_2021" /* || dataset == "19p6_Gev_2019" */)
	{
		if(vz < -145 || vz > 145) return false;
		if(( (vx)*(vx) + (vy)*(vy) ) > 4) return false; //sqrt 4=2 of course
		return true;
	}
	else if(dataset == "27_Gev_2018")
	{
		if(vz < -70 || vz > 70) return false;
		if(( (vx)*(vx) + (vy)*(vy) ) > 1) return false; //sqrt 4=2 of course
		return true;
	}
	else if(dataset == "19p6_Gev_2019")
	{
		if(vz < -145 || vz > 145) return false;
		if(( (vx)*(vx) + (vy)*(vy) ) > 4) return false; //sqrt 4=2 of course
		return true;
	}
	else if(dataset == "3p85_Gev_2018")
	{
		if(vz < 198 || vz > 202)  return false;
		if(( (vx*vx) + (vy+2.0)*(vy+2.0) ) > 4) return false;
		return true;
	}
	return false;
}

bool Shift::isGoodTrack(const StPicoTrack *ptrk) {
    const Float_t pt  = ptrk->pMom().Perp(); 
    const Float_t mom = ptrk->pMom().Mag();
    const Float_t eta = ptrk->pMom().PseudoRapidity();
    const Int_t nHits = ptrk->nHits(); 
    const Float_t dca = ptrk->gDCA( picoEvent->primaryVertex() ).Mag();
    const Int_t nHitsFit = ptrk->nHitsFit();
    const Int_t nHitsPoss = ptrk->nHitsMax();
    const Float_t nHitsDedx = ptrk->nHitsDedx();
    const Float_t quality = (Float_t)nHitsFit/(Float_t)nHitsPoss;

    //if( pt < 0.06 )  return false;
    //if( mom < 0 )   return false;
    //if( eta < -2.5 || 0 < eta) return false;
    if( fabs(dca)>3.0 ) return false;  
    //if( nHits < 10 )  return false;
    if( nHitsFit < 15 )  return false; 
    if( dataset != "19p6_Gev_2019" &&  quality < 0.52 )  return false;
    if( dataset != "19p6_Gev_2019" &&  nHitsDedx < 10 ) return false;

    return true;
}
Double_t Shift::Pileup_rejection(Double_t gRefmult, Int_t gNtofMatch, Double_t gvz,const StPicoEvent *event)
{
	int vz_bin =abs((gvz+145)/10); if(vz_bin<0||vz_bin>28){return -1.0;}
	double b0=10000000,b1=0,b2=0,b3=0,b4=0;
	double d0=0,d1=0,d2=0,d3=0,d4=0;
	double vz_corr_7p7[29]={1.00654,1.00053,0.998852,0.996999,0.996159,0.996542,0.996349,1.00047,1.00468,1.00304,0.99581,0.9983,0.999919,0.999749,
							1,1.00118,0.9992,0.999113,0.995557,1.00442,1.00548,0.999452,0.996127,0.996413,0.996018,0.996163,0.998376,0.998865,1.00318};
	double vz_corr_9p2_trig1[29]={
				0.981917, 0.981996, 0.980511, 0.97507, 0.974313, 0.973339, 0.964442, 0.975903, 0.980922, 0.982125, 0.97381, 0.979041, 0.987953, 0.994754, 1.0,
				1.00207,  1.01171,  1.01019,  1.01176, 1.02353,  1.02679,  1.02377,  1.02388,  1.02206,  1.03417,  1.03343, 1.04365,  1.05574,  1.07789};
	double vz_corr_9p2_trig2[29]={
				1.00894,  1.00431,  1.00199,  1.00086,  1.00139, 1.00079,  1.00023,  1.00587,  1.01,     1.00751,  0.999888, 1.00118,  1.00301, 1.00237, 1.0,
				0.997929, 0.995314, 0.990985, 0.98731, 0.993973, 0.994798, 0.989476, 0.985773, 0.986552, 0.986979, 0.989481, 0.993353, 1.00088, 1.01531};
	gRefmult+=1.0*r1->Uniform()-0.5;
	if(dataset == "7p7_Gev_2021")
	{
		gRefmult = vz_corr_7p7[vz_bin]*gRefmult;
		if(gvz<-87.0)// -145,-87 cm
		{
			b0 = 39.578630496797,b1 = 1.46561577132993,b2 = 0.006515367058115,b3 = -4.06391982010589e-05,b4 = 5.51203917383809e-08;
			d0 = -14.8817460248614,d1 = 0.764539480062978,d2 = 0.00368901349656326,d3 = -1.27602217700865e-05,d4 = 8.02618485000158e-10;
		}
		else if(gvz<-29.0)// -87,-29 cm
		{
			b0 = 26.1841414192908,b1 = 1.73354655107464,b2 = -0.00280668326418846,b3 = 1.22370803379957e-05,b4 = -3.15068617200212e-08;
			d0 = -13.1831127837376,d1 = 0.760227210117286,d2 = 0.00195873375843822,d3 = -2.69378951644624e-06,d4 = -1.05344843941749e-08;
		}
		else if(gvz<29.0)// -29,29 cm
		{
			b0 = 23.3635904884101,b1 = 1.58179764458174,b2 = -0.00100184372825271,b3 = 7.76378744751984e-07,b4 = -6.46469867000365e-09;
			d0 = -11.4340781454132,d1 = 0.72398407747444,d2 = 0.00121092416745035,d3 = 1.17875404059176e-07,d4 = -9.81658682040738e-09;
		}
		else if(gvz<87.0)// 29,87 cm
		{
			b0 = 29.4343991835005,b1 = 1.48353715105631,b2 = 0.00106271734149745,b3 = -9.07835076338586e-06,b4 = 6.7722581625238e-09;
			d0 = -9.97159163811459,d1 = 0.591000613390771,d2 = 0.00449768928484714,d3 = -1.71667412152202e-05,d4 = 1.6467383813372e-08;
		}
		else// 87,145 cm
		{
			b0 = 37.0772875081557,b1 = 1.53484162926915,b2 = 0.00471873506675937,b3 = -2.94958548877277e-05,b4 = 3.60887574265838e-08;
			d0 = -13.3927733032856,d1 = 0.704319390196747,d2 = 0.00485360248820988,d3 = -2.10416804123978e-05,d4 = 1.92342533435503e-08;
		}
	}
	else if(dataset =="9p2_GeV_2020")
	{
		if(event->isTrigger(780010)){gRefmult = vz_corr_9p2_trig1[vz_bin]*gRefmult;}
		else if(event->isTrigger(780020)){gRefmult = vz_corr_9p2_trig2[vz_bin]*gRefmult;}
		if(gvz<-87.0)// -145,-87 cm
		{
			b0=25.6055790979197, b1=2.02528136596901, b2=-0.0058370984051939, b3=2.59602314466234e-05, b4=-5.3014743584261e-08;
			d0=-17.7059596791057, d1=0.614538168662738, d2=0.00534180935164814, d3=-1.79582873880806e-05, d4=1.01623054170579e-08;
		}
		else if(gvz<-29.0)// -87,-29 cm
		{
			b0=23.0160060308621, b1=1.61885832757588, b2=-0.00275873189631398, b3=1.31262550392554e-05, b4=-2.94368020941846e-08;
			d0=-17.3591842617911, d1=0.796170989774258, d2=0.000670722514533827, d3=3.26258075150876e-06, d4=-1.60611460182112e-08;
		}
		else if(gvz<29.0)// -29,29 cm
		{
			b0=16.4277056306649, b1=1.71652229539398, b2=-0.00406847684302521, b3=1.65203560938885e-05, b4=-2.96250329214512e-08;
			d0=-15.7887025834219, d1=0.789786364309292, d2=-0.000637115144252616, d3=1.00019972792727e-05, d4=-2.45208851616324e-08;
		}
		else if(gvz<87.0)// 29,87 cm
		{
			b0=21.2024767158778, b1=1.70521848381614, b2=-0.00352260930859763, b3=1.60905730948817e-05, b4=-3.37443468806432e-08;
			d0=-17.1166088395929, d1=0.814739436616432, d2=0.000227197779215977, d3=6.55397838050604e-06, d4=-2.28812912596058e-08;
		}
		else// 87,145 cm
		{
			b0=26.0970905882739, b1=1.88889714311734, b2=-0.00195374948885512, b3=-6.14244087431038e-06, b4=1.99930095058841e-08;
			d0=-15.6624325989392, d1=0.52385751891358, d2=0.00794996911844969, d3=-4.09239155250494e-05, d4=6.40163739983216e-08;
		}
	}
	else if(dataset =="14p6_gev_2019") //https://drupal.star.bnl.gov/STAR/system/files/Centrality_Study_at_14p6_final_0.pdf
	{
		b0 = 36.4811873257854, b1 = 1.96363692967013, b2 = -0.00491528146300182, b3 = 1.45179464078414e-05, b4 = -1.82634741809226e-08;
		d0 = -16.176117733536, d1 = 0.780745107634961, d2 = -2.03347057620351e-05, d3 = 3.80646723724747e-06,d4 = -9.43403282145648e-09;
	}
	if(gRefmult<d0+
			  d1*gNtofMatch+
			  d2*gNtofMatch*gNtofMatch+
			  d3*gNtofMatch*gNtofMatch*gNtofMatch+
			  d4*gNtofMatch*gNtofMatch*gNtofMatch*gNtofMatch){return -1;}
	if(gRefmult>b0+
			  b1*gNtofMatch+
			  b2*gNtofMatch*gNtofMatch+
			  b3*gNtofMatch*gNtofMatch*gNtofMatch+
			  b4*gNtofMatch*gNtofMatch*gNtofMatch*gNtofMatch){return -1;}
	return gRefmult;
}
bool Shift::isGoodTrigger(const StPicoEvent* event)
{
	if(dataset == "9p2_GeV_2020")
	{
		if(event->isTrigger(780010)|| event->isTrigger(780020)) return true; //9p2GeV minbias only
	}
	else if(dataset == "7p7_Gev_2021")
	{
		if(event->isTrigger(810010)|| event->isTrigger(810020)|| event->isTrigger(810030)|| event->isTrigger(810040)) return true; //7p7GeV minbias only
	}
	else if(dataset == "11p5_gev_2020")
	{
		if(event->isTrigger(710000)|| event->isTrigger(710010)|| event->isTrigger(710020)) return true; //11p5 GeV minbias only
	}
	else if(dataset == "3p85_Gev_2018")
	{
		if(event->isTrigger(620052)) return true;//3p85 GeV minbias
	}
	else if(dataset == "14p6_gev_2019")
	{
		if(event->isTrigger(650000)) return true; //minbias only
	}
	else if(dataset == "17p3_gev_2021")
	{
		if(event->isTrigger(870010)) return true; //minbias only
	}
	else if(dataset == "19p6_Gev_2019")
	{
		if(event->isTrigger(640001) || event->isTrigger(640011) || event->isTrigger(640021) || event->isTrigger(640031) || event->isTrigger(640041) || event->isTrigger(640051)) return true;
	}
	else if(dataset == "27_Gev_2018")
	{
		if(event->isTrigger(1)||event->isTrigger(610001)||event->isTrigger(610011)||event->isTrigger(610021)||event->isTrigger(610031) ||event->isTrigger(610041)||event->isTrigger(610051)) return true;//27 gev
	}
	return false;
}
bool Shift::isProton(const StPicoTrack *ptrk, bool req_tof, int option)
{
	double trackP = ptrk->pPtot();
	double pt = ptrk->pPt();
	Float_t mean_nsigp =0;
	double n_sigma_cut = 3.0; if(dataset == "19p6_Gev_2019"){n_sigma_cut = 2.0;}
	double mass2_min = 0.8;
	double mass2_max = 1.0;
	double dca = ptrk->gDCA( picoEvent->primaryVertex() ).Mag();
	if(dataset == "3p85_Gev_2018")
	{
		double mean_cent5_P[25] = {-1.39109,-0.961178,-0.579124,-0.256328,-0.0156338,0.16313,0.260721,0.42406,0.510594,0.472899,0.524993,0.623829,0.650686,0.601457,0.636739,0.665146,0.680327,0.685059,0.685086,0.68105,0.679174,0.646809,0.644823,0.64806,0.659472,};
		mean_nsigp = mean_cent5_P[momindex(trackP)];
		n_sigma_cut = 2.0;
		mass2_min = 0.6;
		mass2_max = 1.2;
	}
	if(option == 1) //dca<1
	{
		if(fabs(dca)>1.0) return false;
	}
	if(option == 2) //dca<2
	{
		if(fabs(dca)>2.0) return false;
	}
	if(option == 3) //nhitsfit<18
	{
		Int_t nHitsFit = ptrk->nHitsFit();
		if(nHitsFit<18) return false;
	}
	if(option == 4) //nhitsfit<20
	{
		Int_t nHitsFit = ptrk->nHitsFit();
		if(nHitsFit<20) return false;
	}
	if(option == 5) //n_sigma_cut = 2.0;
	{
		n_sigma_cut = 2.0;
		if(dataset == "3p85_Gev_2018"){n_sigma_cut = 2.5;}
	}
	if(option == 6) //n_sigma_cut = 2.5;
	{
		n_sigma_cut = 2.5;
		if(dataset == "3p85_Gev_2018"){n_sigma_cut = 3.0;}
	}
	if(option == 7) //loose mass;
	{
		mass2_min = 0.6;
		mass2_max = 1.2;
	}
	if(option == 8) //eta<1;
	{
		if(TMath::Abs(ptrk->pMom().PseudoRapidity())>1.0) return false;
	}
	if(ptrk->isTofTrack())
	{
		StPicoBTofPidTraits *trait = mPicoDst->btofPidTraits(ptrk->bTofPidTraitsIndex());
		if(trait)
		{
			double btofMatchFlag = trait->btofMatchFlag();
            double btofYLocal    = trait->btofYLocal();
			if(!(fabs(btofYLocal) < 1.8 && btofMatchFlag > 0)) return false;
			
			double Beta = trait->btofBeta();
			double mass2 = trackP*trackP*((1.0/(Beta*Beta))-1.0);
			if((mass2>mass2_min && mass2<mass2_max)||!req_tof)
			{
				if(fabs(ptrk->nSigmaProton() - mean_nsigp) < n_sigma_cut) return true;
			}
			return false;
		}
	}
	else if(ptrk->isETofTrack() && include_etof)
	{
		StPicoETofPidTraits *etrait = mPicoDst->etofPidTraits(ptrk->eTofPidTraitsIndex());
		if(etrait)
		{
			double etofMatchFlag = etrait->matchFlag();
			if(!(etofMatchFlag > 0)) return false;
	
			double Beta = etrait->beta();
			double mass2 = trackP*trackP*((1.0/(Beta*Beta))-1.0);
			if((mass2>mass2_min && mass2<mass2_max)||!req_tof)
			{
				if(fabs(ptrk->nSigmaProton() - mean_nsigp) < n_sigma_cut) return true;
			}
			return false;
		}
	}
	else if(fabs(ptrk->nSigmaProton() - mean_nsigp) < n_sigma_cut &&!req_tof) return true; 
	else if(dataset == "3p85_Gev_2018" && fabs(ptrk->nSigmaProton() - mean_nsigp) < n_sigma_cut && trackP<2.0 ) return true;
	return false;
}
bool Shift::isProton_test(const StPicoTrack *ptrk, bool req_tof, int option)
{
	double trackP = ptrk->pPtot();
	double pt = ptrk->pPt();
	Float_t mean_nsigp =0;
	double n_sigma_cut = 2.5;
	double mass2_min = 0.8;
	double mass2_max = 1.0;
	double dca = ptrk->gDCA( picoEvent->primaryVertex() ).Mag();
	if(dataset == "9p2_GeV_2020" && fabs(dca)>1.0 && option != 1 && option != 2){ return false;  }
	if(dataset == "19p6_Gev_2019")
	{
		Float_t  mean[15]={0.883590 ,0.880394 ,0.879022 ,0.878588 ,0.877887 ,0.876960 ,0.875662 ,0.874262 ,0.873621,0.871457,0.868933,0.866071,0.861743,0.857418,0.848790/*,0.833930,0.835822,0.781637*/};
		Float_t sigma[15]={0.0337246,0.0310994,0.0328318,0.0385098,0.0466531,0.0569312,0.0693882,0.0876769,0.103482 ,0.121335,0.141623,0.164989,0.189501,0.222519,0.254746/*,0.301722,0.343749,0.409971*/};
		if(0.2<pt&&pt<=0.4){mass2_min = mean[0]-2.0*sigma[0];mass2_max = mean[0]+2.0*sigma[0];}
		else if(0.4<pt&&pt<=0.6){mass2_min = mean[0]-2.0*sigma[0];mass2_max = mean[0]+2.0*sigma[0];}
		else if(0.6<pt&&pt<=0.8){mass2_min = mean[1]-2.0*sigma[1];mass2_max = mean[1]+2.0*sigma[1];}
		else if(0.8<pt&&pt<=1.0){mass2_min = mean[2]-2.0*sigma[2];mass2_max = mean[2]+2.0*sigma[2];}
		else if(1.0<pt&&pt<=1.2){mass2_min = 0.6;mass2_max = 1.2;}
		else if(1.2<pt&&pt<=1.4){mass2_min = 0.6;mass2_max = 1.2;}
		else if(1.4<pt&&pt<=1.6){mass2_min = 0.6;mass2_max = 1.2;}
		else if(1.6<pt&&pt<=1.8){mass2_min = 0.6;mass2_max = 1.2;}
		else if(1.8<pt&&pt<=2.0){mass2_min = 0.7;mass2_max = 1.2;}
		else if(2.0<pt&&pt<=2.2){mass2_min = 0.7;mass2_max = 1.2;}
		else if(2.2<pt&&pt<=2.4){mass2_min = 0.7;mass2_max = 1.2;}
		else if(2.4<pt&&pt<=2.6){mass2_min = 0.7;mass2_max = 1.2;}
		else if(2.6<pt&&pt<=2.8){mass2_min = 0.8;mass2_max = 1.4;}
		else if(2.8<pt&&pt<=3.0){mass2_min = 0.8;mass2_max = 1.4;}
		else if(3.0<pt&&pt<=3.2){mass2_min = 0.8;mass2_max = 1.4;}
	}
	if(dataset == "3p85_Gev_2018")
	{
		double mean_cent5_P[25] = {-1.39109,-0.961178,-0.579124,-0.256328,-0.0156338,0.16313,0.260721,0.42406,0.510594,0.472899,0.524993,0.623829,0.650686,0.601457,0.636739,0.665146,0.680327,0.685059,0.685086,0.68105,0.679174,0.646809,0.644823,0.64806,0.659472,};
		mean_nsigp = mean_cent5_P[momindex(trackP)];
		n_sigma_cut = 2.0;
		mass2_min = 0.6;
		mass2_max = 1.2;
	}
	if(dataset == "27_Gev_2018")
	{
		n_sigma_cut = 2.0;
		mass2_min = 0.8;
		mass2_max = 1.0;
	}
	if(option == 1) //dca<1
	{
		if(dataset == "9p2_GeV_2020" && fabs(dca)>3.0){ return false;  }
		else if(dataset != "9p2_GeV_2020" && fabs(dca)>1.0) return false;
	}
	if(option == 2) //dca<2
	{
		if(fabs(dca)>2.0) return false;
	}
	if(option == 3) //nhitsfit<18
	{
		Int_t nHitsFit = ptrk->nHitsFit();
		if(nHitsFit<18) return false;
	}
	if(option == 4) //nhitsfit<20
	{
		Int_t nHitsFit = ptrk->nHitsFit();
		if(nHitsFit<20) return false;
	}
	if(option == 5) //n_sigma_cut = 2.0;
	{
		n_sigma_cut = 2.0;
		if(dataset == "3p85_Gev_2018" ||dataset == "27_Gev_2018"){n_sigma_cut = 2.5;}
	}
	if(option == 6) //n_sigma_cut = 3.0;
	{
		n_sigma_cut = 3.0;
	}
	
	if(option == 7) //tight mass;
	{
		mass2_min = 0.6;
		mass2_max = 1.2;
	}

	if(ptrk->isTofTrack())
	{
		StPicoBTofPidTraits *trait = mPicoDst->btofPidTraits(ptrk->bTofPidTraitsIndex());
		if(trait)
		{
			double btofMatchFlag = trait->btofMatchFlag();
            double btofYLocal    = trait->btofYLocal();
			if(!(fabs(btofYLocal) < 1.8 && btofMatchFlag > 0)) return false;
			
			double Beta = trait->btofBeta();
			double mass2 = trackP*trackP*((1.0/(Beta*Beta))-1.0);
			if((mass2>mass2_min && mass2<mass2_max)||!req_tof)
			{
				if(fabs(ptrk->nSigmaProton() - mean_nsigp) < n_sigma_cut) return true;
			}
			return false;
		}
	}
	else if(ptrk->isETofTrack())
	{
		StPicoETofPidTraits *etrait = mPicoDst->etofPidTraits(ptrk->eTofPidTraitsIndex());
		if(etrait)
		{
			double etofMatchFlag = etrait->matchFlag();
			if(!(etofMatchFlag > 0)) return false;
	
			double Beta = etrait->beta();
			double mass2 = trackP*trackP*((1.0/(Beta*Beta))-1.0);
			if((mass2>mass2_min && mass2<mass2_max)||!req_tof)
			{
				if(fabs(ptrk->nSigmaProton() - mean_nsigp) < n_sigma_cut) return true;
			}
			return false;
		}
	}
	else if(fabs(ptrk->nSigmaProton() - mean_nsigp) < n_sigma_cut && trackP<0.8) return true; 
	else if(fabs(ptrk->nSigmaProton() - mean_nsigp) < n_sigma_cut &&!req_tof) return true; 
	else if(dataset == "3p85_Gev_2018" && fabs(ptrk->nSigmaProton() - mean_nsigp) < n_sigma_cut && trackP<2.0 ) return true;
	return false;
}
Double_t Shift::proton_efficiency_correction(const StPicoTrack *ptrk, Int_t cent)
{
	if(!do_tof_efficiency){return 1.0;}
	double pz = ptrk->pMom().Z();
	double pt = ptrk->pPt();
	double trackP = ptrk->pPtot();
	double eta = ptrk->pMom().Eta();
	int charge = ptrk->charge();
	double energy_P = sqrt(trackP*trackP + protonMass*protonMass);
	double rap_P = -0.5*log( (energy_P + pz) / (energy_P - pz) )-y_cm;
	Double_t numerator = 1.0;
	Double_t denominator = 1.0;
	Double_t efficiency = 1.0;
	denominator = hist_eta_y_tpc_hit->GetBinContent(hist_eta_y_tpc_hit->FindBin(eta,pt));
	numerator = hist_eta_y_tof_match_hit->GetBinContent(hist_eta_y_tof_match_hit->FindBin(eta,pt));
	//TPC efficiency:
	// if(charge>0 && dataset == "9p2_GeV_2020")
	// {
		// denominator = hist_pt_y_mc_proton_cent[cent]->GetBinContent(hist_pt_y_mc_proton_cent[cent]->FindBin(rap_P,pt));
		// numerator = hist_pt_y_match_proton_cent[cent]->GetBinContent(hist_pt_y_match_proton_cent[cent]->FindBin(rap_P,pt));
	// }
	// else if(charge<0 && dataset == "9p2_GeV_2020")
	// {
		// denominator = hist_pt_y_mc_aproton_cent[cent]->GetBinContent(hist_pt_y_mc_aproton_cent[cent]->FindBin(rap_P,pt));
		// numerator = hist_pt_y_match_aproton_cent[cent]->GetBinContent(hist_pt_y_match_aproton_cent[cent]->FindBin(rap_P,pt));
	// }
	if(denominator!=0.0 && numerator!=0.0 && numerator<denominator){efficiency = numerator/denominator;}
	return 1.0/efficiency;
}
bool Shift::isPion(const StPicoTrack *ptrk)
{
	double trackP = ptrk->pPtot();
	Float_t mean_nsigp =0;
	if(ptrk->nHitsDedx()<5){return false;}
	if(fabs(ptrk->nSigmaPion() - mean_nsigp) < 3.0)
	{
		if(ptrk->isTofTrack())
		{
			StPicoBTofPidTraits *trait = mPicoDst->btofPidTraits(ptrk->bTofPidTraitsIndex());
			if(trait)
			{
				double btofMatchFlag = trait->btofMatchFlag();
				double btofYLocal    = trait->btofYLocal();
				if(!(fabs(btofYLocal) < 1.8 && btofMatchFlag > 0)) return false;
		
				double Beta = trait->btofBeta();
				double mass2 = trackP*trackP*((1.0/(Beta*Beta))-1.0);
				if(fabs(mass2-pionMass*pionMass)/(pionMass*pionMass)<0.3){return true;}
				return false;
			}
		}
		else if(ptrk->isETofTrack())
		{
			// Retrieve corresponding trait
			StPicoETofPidTraits *etrait = mPicoDst->etofPidTraits(ptrk->eTofPidTraitsIndex());
			if(etrait)
			{
				double etofMatchFlag = etrait->matchFlag();
				if(!(etofMatchFlag > 0)) return false;
		
				double Beta = etrait->beta();
				double mass2 = trackP*trackP*((1.0/(Beta*Beta))-1.0);
				if(fabs(mass2-pionMass*pionMass)/(pionMass*pionMass)<0.3){return true;}
				return false;
			}
		}
		if(dataset == "3p85_Gev_2018"){if(trackP<1.0){return true;}}
		return false;	
	}
	return false;
}
bool Shift::isKaon(const StPicoTrack *ptrk)
{
	double trackP = ptrk->pPtot();
	Float_t mean_nsigp =0;
	if(ptrk->nHitsDedx()<5){return false;}
	if(fabs(ptrk->nSigmaKaon() - mean_nsigp) < 3.0)
	{
		if(ptrk->isTofTrack())
		{
			StPicoBTofPidTraits *trait = mPicoDst->btofPidTraits(ptrk->bTofPidTraitsIndex());
			if(trait)
			{
				double btofMatchFlag = trait->btofMatchFlag();
				double btofYLocal    = trait->btofYLocal();
				if(!(fabs(btofYLocal) < 1.8 && btofMatchFlag > 0)) return false;
		
				double Beta = trait->btofBeta();
				double mass2 = trackP*trackP*((1.0/(Beta*Beta))-1.0);
				if(fabs(mass2-kaonMass*kaonMass)/(kaonMass*kaonMass)<0.3){return true;}
				return false;
			}
		}
		else if(ptrk->isETofTrack())
		{
			// Retrieve corresponding trait
			StPicoETofPidTraits *etrait = mPicoDst->etofPidTraits(ptrk->eTofPidTraitsIndex());
			if(etrait)
			{
				double etofMatchFlag = etrait->matchFlag();
				if(!( etofMatchFlag > 0)) return false;
		
				double Beta = etrait->beta();
				double mass2 = trackP*trackP*((1.0/(Beta*Beta))-1.0);
				if(fabs(mass2-kaonMass*kaonMass)/(kaonMass*kaonMass)<0.3){return true;}
				return false;
			}
		}
		if(dataset == "3p85_Gev_2018"){if(trackP<1.0){return true;}}
		return false;	
	}
	return false;
}
int Shift::momindex(double mom)
{
	//only needed for 3p85_Gev_2018 proton
    if(mom > 0.0 && mom <= 0.1) return 0;
    else if(mom > 0.1 && mom <= 0.2) return 1;
    else if(mom > 0.2 && mom <= 0.3) return 2;
    else if(mom > 0.3 && mom <= 0.4) return 3;
    else if(mom > 0.4 && mom <= 0.5) return 4;
    else if(mom > 0.5 && mom <= 0.6) return 5;
    else if(mom > 0.6 && mom <= 0.7) return 6;
    else if(mom > 0.7 && mom <= 0.8) return 7;
    else if(mom > 0.8 && mom <= 0.9) return 8;
    else if(mom > 0.9 && mom <= 1.0) return 9;
    else if(mom > 1.0 && mom <= 1.1) return 10;
    else if(mom > 1.1 && mom <= 1.2) return 11;
    else if(mom > 1.2 && mom <= 1.3) return 12;
    else if(mom > 1.3 && mom <= 1.4) return 13;
    else if(mom > 1.4 && mom <= 1.5) return 14;
    else if(mom > 1.5 && mom <= 1.6) return 15;
    else if(mom > 1.6 && mom <= 1.7) return 16;
    else if(mom > 1.7 && mom <= 1.8) return 17;
    else if(mom > 1.8 && mom <= 1.9) return 18;
    else if(mom > 1.9 && mom <= 2.0) return 19;
    else if(mom > 2.0 && mom <= 2.1) return 20;
    else if(mom > 2.1 && mom <= 2.2) return 21;
    else if(mom > 2.2 && mom <= 2.3) return 22;
    else if(mom > 2.3 && mom <= 2.4) return 23;
    else if(mom > 2.4              ) return 24;
}
int Shift::EP_group(int run_number, int cent)
{
	int ep_group_num = 0;
	if(dataset == "27_Gev_2018")
	{
		int event_cuts[14]= {19131037, 19135016, 19137041, 19139063, 19140030, 19141030, 19144012, 19144033,
						   19145034, 19147021, 19147048, 19155057, 19158020, 19268002};
		for(int i = 0;i<14;i++)
		{
			if(run_number<=event_cuts[13-i]){ ep_group_num=(13-i)*9+cent;}
		}
	}
	else
	{
		ep_group_num = GetRunIndex(run_number)*9+cent;
	}
	return ep_group_num;
}
bool Shift::zero2twopi(Double_t &tphi)
{
	if(!(-1e9<tphi&&tphi<1e9)){return 1;}
	if(tphi < 0.0   ){tphi += 2.0*pi;}
	if(tphi > 2.0*pi){tphi -= 2.0*pi;}
	if(0.0 < tphi && tphi < 2.0*pi){return 0;}
	return zero2twopi(tphi);
}
//_________________________________________________________
