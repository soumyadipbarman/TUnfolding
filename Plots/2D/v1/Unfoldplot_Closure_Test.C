#include "TCanvas.h"
#include "TStyle.h"
#include "TH1.h"
#include "TGaxis.h"
#include "TRandom.h"
#include "TAxis.h"
#include "TMath.h"
#include "TGraph.h"
#include "TObject.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string.h>
#include <fstream>
#include <cmath>
#include <TH2D.h>
#include <TTree.h>
#include "TVector.h"
#include <vector>
#include <TF1.h>
#include <TProfile.h>
#include <TStyle.h>
#include "TPostScript.h"
#include "CMS_lumi.C"

#include "TUnfoldBinning.h"
#include "TUnfoldBinningXML.h"
#include "TUnfoldDensity.h"
#include "TUnfold.h"
#include "TUnfoldSys.h"
#include <TVectorD.h>
#include <TDecompSVD.h>

#define PYTHIA
#define MADGRAPH
#define Herwig
#define CLOUSER

void Unfoldplot(){
  static const int unfold_ty =1;  //Unfold method to be plots 
  static const int clos_ty =4;    //Number of Closure Test performed 
  static const int nHLTmx=10;     //HT2 Range
  static const int nusedvar = 5;  //Event Shape variables used
  static const int ntype =2;
  static const int nmc = 4;       //Number of MC sample 
  static const int ndd = 3; 
  static const int umc = 0;       //0 for Py8, 1 for MG , 2 for Herwig : which MC have used in Unfolding

  static const int ndef=3;
  static const int njet=2;
  static const int nkappa=10;

  bool isstat =1;  int irbin=1;
  char histname[100],Title[100], Xaxis[100], Yaxis[100], ratioXaxis[100], ratioYaxis[100],pdfname[100],pdfname1[100],pdfname2[100],LegName[100];
  bool Reco, Gen;
  Int_t color[10] ={2,4,6,5,6,46,3,28,38,42};  //define the color for different histograms
  Int_t var[nusedvar]={3,9,15,18,24};
  Int_t HT2range[nHLTmx+1]={92, 119, 185, 251, 319, 388, 467, 518, 579, 669, 3000};

  const int njetetamn=1;     //eta value used 2.5
  static const int nvar=32;  //Total number of eventshape variables
  
  string  Datadir = "Data2D";
  string folddir = "Folded2D";
  string unfdir = "Unfold2D";
  const char*  dirname[4]={"Pythia8","MG8","HW7","Py8Flat"};
#ifdef CLOUSER
  string mcdir[4]={"Pythia8","Py8Flat","MG8","HW7"}; //Sequence For Closure Test
#else
  string mcdir[4]={"Pythia8","MG8","HW7","Py8Flat"}; //Sequence For Unfold and Refold plot
#endif
  const char* obs_def[3]={"Q_{D}","Q_{L}","Q_{T}"};
  const char* jet_num[2]={"Leading-Jet","Sub-Leading-Jet"};
  const char* njets[2]={"1","2"};
  const char* k_fact[10]={"k=0.1","k=0.2","k=0.3","k=0.4","k=0.5","k=0.6","k=0.7","k=0.8","k=0.9","k=1.0"};
  
  const char* Esvsym[5] = {"#tau_{_{#perp} }", "#rho_{Tot}","Y_{2,3}","B_{ T}","#rho^{T}_{Tot}"};  
  const char* Esvname[5] = {"Complement of transverse thrust", "Total jet mass","Three-jet resolution ","Total Jet broadening","Total transverse jet mass"};
  const char* htrang[10]={"92 < P_{T} < 119", "119 < P_{T} < 185", "185 < P_{T} < 251", "251 < P_{T} < 319", "319 < P_{T} < 388", "388 < P_{T} <467", "467 < P_{T} <518","518 < P_{T} < 579", "579 < P_{T} < 669", "P_{T} > 669"};
  const char* obs_logy[10]={"1/N dN/d","1/N dN/d","1/N dN/d","1/N dN/d","1/N dN/d","1/N dN/d","1/N dN/d","1/N dN/d","1/N dN/d","1/N dN/d"};
  const char* Esvlogx[5] ={"ln(#tau_{ _{#perp } })","ln(#rho_{Tot})","ln(Y_{2,3})", "ln(B_{ T})","ln(#rho^{T}_{Tot})"};
  const char* Esvlogy[5] = {"1/N dN/dln(#tau_{ _{#perp } })","1/N dN/dln(#rho_{Tot})","1/N dN/dln(Y_{2,3})","1/N dN/dln(B_{ T})","1/N dN/dln(#rho^{T}_{Tot})"};

  const  char* Validity_test[4]={"Closure test","Bottom Line test"," Unfolded","Refold"};
  const  char* h2dMat_name[4]={"Covariance matrix","correlation coefficients"," probabilities matrix ","Response matrix"};

#ifdef CLOUSER
  const  char* mcname[4]={"Py8 Gen","Py8Flat Gen","Madgraph Gen","Herwig Gen"};
#else 
  const  char* mcname[4]={"Pythia8","Madgraph","Herwig7","Pythia8 Flat"};
#endif
  const  char* DataEra[3]={"Data","Data","Data"};
  const  char* mcnamerco[4]={"Pythia8 ","Madgraph ","Herwig7","PY8 Flat"};
  const  char* closuretype[4]={"PY8 by PY8","PY8 Flat By PY8 ", "MG By PY8","HW7 By PY8"};
  const  char* itypeN[ntype]={"Jets","Charged Particles"}; 
  const  char* RefoldEra[5]={"Refold Pythia8(No Regularisation) ","Refold Pythia8(L-Curve scan)","Refold Pythia8(Scan Tau)","Refold Pythia8(ScanSURE)","Refold Pythia8(Iterative method)"};
#ifdef CLOUSER
  const  char* Unfoldtype[5]={"Unfold PY8 by Py8 RM","Unfold Py8 Flat by Py8 RM","Unfold Madgraph by PY8 RM","Unfold Herwig by Py8 RM ","TUnfold(Iterative method)"}; 
#else 
  const  char* Unfoldtype[5]={"Unfold by Py8 RM","TUnfold(L-Curve scan)","TUnfold(Scan Tau)","TUnfold(ScanSURE)","TUnfold(Iterative method)"};
#endif
  const  char* Modelnm[3]={"Pythia8","Madgraph","Herwig"};
  const  char* Methodtype[5]={"(No Regularisation) ","(L-Curve scan)","(Scan Tau)","(ScanSURE)","(Iterative method)"};
  const  char* smeared[5]={"TUnfold","Refold","Folded-back","GEN","RECO"};
  static const int iera = 1;
  int iPeriod = 0;  int iPos=10 ;
  
  //TFile *inputbinning = TFile::Open("testUnfold5_binning.root");
 
  TFile *Unf_root[clos_ty];
  TFile *Unf_dataCT[clos_ty];
  
  //******************************************************************************
  //TFile *Unfoldroot = TFile::Open("/home/suman/Paradox/Charged_ESV/Working/Unfolding/Tunfold_2D_HT2Merge/Tunfold_2D/TUnfolding/MC_Unfolded_Result.root");  // Unfolded data 
  //TFile *Unfoldroot = TFile::Open("/home/suman/Paradox/Charged_ESV/Working/Unfolding/TUnfold_check_CT_1Jan21/Test_14Jan21/TUnfolding/Unfolded_Result.root");  // Unfolded data 
  //TFile *Unfoldroot = TFile::Open("Data_Unfolded_Result.root");  // Unfolded data 
  //TFile *Unfoldroot = TFile::Open("Unfolded_Result.root");       // Unfolded data 
  TFile *Unfoldroot = TFile::Open("Unfolded_Result_data.root");
  
  //Unf_root[0] = TFile::Open("/home/suman/Paradox/Charged_ESV/Working/Unfolding/TUnfold_check_CT_1Jan21/Test_14Jan21/TUnfolding/Unfolded_Py8_Py8_14Jan21.root");  // Py- py
  //Unf_root[1] = TFile::Open("/home/suman/Paradox/Charged_ESV/Working/Unfolding/TUnfold_check_CT_1Jan21/Test_14Jan21/TUnfolding/Unfolded_Py8Flat_Py8_14Jan21.root");  // Py -Py flat 
  //Unf_root[2]= TFile::Open("/home/suman/Paradox/Charged_ESV/Working/Unfolding/TUnfold_check_CT_1Jan21/Test_14Jan21/TUnfolding/Unfolded_MG_Py8_14Jan21.root");  // Py MG 
  //Unf_root[3]= TFile::Open("/home/suman/Paradox/Charged_ESV/Working/Unfolding/TUnfold_check_CT_1Jan21/Test_14Jan21/TUnfolding/Unfolded_HW_Py8_14Jan21.root");  // Py MG 
 
  Unf_root[0] = TFile::Open("Unfolded_Result_closure.root");  // Py- py
  Unf_root[1] = TFile::Open("Unfolded_Result_closure.root");  // Py -Py flat
  Unf_root[2]= TFile::Open("Unfolded_Result_closure.root");  // Py MG
  Unf_root[3]= TFile::Open("Unfolded_Result_closure.root");  // Py MG

  //Unf_root[0] = TFile::Open("/home/suman/Paradox/Charged_ESV/Working/Unfolding/TUnfold_check_CT_1Jan21/Test_14Jan21/TUnfolding/Unfolded_PY8_MG_17Jan21.root");  // MG- PY
  //Unf_root[1] = TFile::Open("/home/suman/Paradox/Charged_ESV/Working/Unfolding/TUnfold_check_CT_1Jan21/Test_14Jan21/TUnfolding/Unfolded_PY8_Flat_MG_17Jan21.root");  // MG -Py flat
  //Unf_root[2]= TFile::Open("/home/suman/Paradox/Charged_ESV/Working/Unfolding/TUnfold_check_CT_1Jan21/Test_14Jan21/TUnfolding/Unfolded_MG_MG_17Jan21.root");  // MG MG
  //Unf_root[3]= TFile::Open("/home/suman/Paradox/Charged_ESV/Working/Unfolding/TUnfold_check_CT_1Jan21/Test_14Jan21/TUnfolding/Unfolded_HW_MG_17Jan21.root");  // MG HW

  //Unf_root[0] = TFile::Open("/home/suman/Paradox/Charged_ESV/Working/Unfolding/TUnfold_check_CT_1Jan21/Test_14Jan21/TUnfolding/Unfolded_PY8_PY8_Lcurve_17Jan21.root");  // Py- py
  //Unf_root[1] = TFile::Open("/home/suman/Paradox/Charged_ESV/Working/Unfolding/TUnfold_check_CT_1Jan21/Test_14Jan21/TUnfolding/Unfolded_PY8Flat_PY8_Lcurve_17Jan21.root");//Py -Py flat
  //Unf_root[2]= TFile::Open("/home/suman/Paradox/Charged_ESV/Working/Unfolding/TUnfold_check_CT_1Jan21/Test_14Jan21/TUnfolding/Unfolded_MG_PY8_Lcurve_17Jan21.root");  // Py MG
  //Unf_root[3]= TFile::Open("/home/suman/Paradox/Charged_ESV/Working/Unfolding/TUnfold_check_CT_1Jan21/Test_14Jan21/TUnfolding/Unfolded_HW_PY8_Lcurve_17Jan21.root");  // Py MG
  
  //Unf_dataCT[0]= TFile::Open("/home/suman/Paradox/Charged_ESV/Working/Unfolding/TUnfold_check_CT_1Jan21/Test_14Jan21/TUnfolding/Unfolded_PY8Flat_PY8.root");  // Py MG 
  //Unf_dataCT[1]= TFile::Open("/home/suman/Paradox/Charged_ESV/Working/Unfolding/TUnfold_check_CT_1Jan21/Test_14Jan21/TUnfolding/Unfolded_PY8Flat_MG.root");  // Py MG 
  //Unf_dataCT[2]= TFile::Open("/home/suman/Paradox/Charged_ESV/Working/Unfolding/TUnfold_check_CT_1Jan21/Test_14Jan21/TUnfolding/Unfolded_PY8Flat_HW.root");  // Py MG 
  Unf_dataCT[0]= TFile::Open("Unfolded_Result_closure.root");  // Py MG
  Unf_dataCT[1]= TFile::Open("Unfolded_Result_closure.root");  // Py MG
  Unf_dataCT[2]= TFile::Open("Unfolded_Result_closure.root");  // Py MG

//----------------------------------------------------
//Function declaration
  void Integralhist(TH1D *hist);
  void divBybinWidth(TH1D *hist);
  void Myplotset(TH1D *Myhist,const char* XTitle, const char* YTitle);
  void Set2dHist(TH2D *MyHist, const char* XTitle, const char* YTitle,const char* ZTitle, double titoff[3], double titsize[3]);
  void SetMycanvas(TCanvas *cpt,double bs,double lm, double rm, double tm,double bm);
  void CTLegend(TLegend *legendn, const char* txt1, const char* txt2);
  TCanvas *ratio_can(int Nplot[2],float plegend[7], TH1D* data, TH1D* MC[Nplot[0]], char* lowpadx,const  char* modnam[Nplot[0]],const  char* datanm[1]);
  TCanvas *ratio_canV2(int Nplot[3],float plegend[7], TH1D* data, TH1D* MC[Nplot[0]], char* lowpadx, const char* modnam[Nplot[0]], const  char* datanm[3]);
  void HT2_Normal(TH2D* VarHT, int nht, int htbin[nht+1], TH1D* Var[nht]);
  void HT2_NormalV2(TH2D* VarHT, int nht, int htbin[nht+1], TH1D* Var[nht]);
  TH1D* ReadHist1D(string name,TFile* root, int irbin=1);
  TH2D* ReadHist2D(string name,TFile* root, int irbin=1);
  TLegend* CTLegendV2(float x1, float y1, float x2, float y2, float txtsize, const char* txt1="", const char* txt2="",const char* txt3="",const char* txt4="");
  void Chi2Root(TH1 * data, TH1 * MC, int rebin = 1);
  void MyplotsetV2(TH1D *MyHist, const char* XT, const char* YT, float mx, float min, int ilw, int ilsty, int imsty, int imstysize, int icl);
//----------------------------------------------------
  TH1D *Data_reco[ndef][njet][nkappa];
  TH2D *Data_reco2d[ndef][njet][nkappa];
  TH1D *MC_gen[nmc][ndef][njet][nkappa];      //MC gen
  TH1D *MC_genM[nmc][ndef][njet][nkappa];     //MC gen
  TH2D *MC_gen2d[nmc][ndef][njet][nkappa];    //MC gen
  TH1D *MC_genmiss[nmc][ndef][njet][nkappa];  //MC Gen-miss
  TH1D *MC_reco[nmc][ndef][njet][nkappa];     //MC Reco 
  TH2D *MC_reco2d[nmc][ndef][njet][nkappa];   //MC Reco 
  TH1D *MC_recofake[nmc][ndef][njet][nkappa]; //MC Reco -Fake
  TH2D *MC_Res[nmc][ndef][njet][nkappa];      //RM
  
  TH1D *hist_eff[nmc][ndef][njet][nkappa];
  TH1D *hist_fake[nmc][ndef][njet][nkappa];
  TH1D *hist_purity[nmc][ndef][njet][nkappa];
  TH1D *hist_stbl[nmc][ndef][njet][nkappa];

  TH1D *UnfoldCT[clos_ty][unfold_ty][ndef][njet][nkappa];
  TH1D *UnfoldCT_M[clos_ty][unfold_ty][ndef][njet][nkappa];
  TH1D *RefoldCT[clos_ty][unfold_ty][ndef][njet][nkappa];
  TH1D *UnfoldCT_HT[clos_ty][unfold_ty][ndef][njet][nkappa][nHLTmx];
  TH1D *RefoldCT_HT[clos_ty][unfold_ty][ndef][njet][nkappa][nHLTmx];

  TH1D *UnfoldDataCT[clos_ty][unfold_ty][ndef][njet][nkappa];
  TH1D *UnfoldDataCT_M[clos_ty][unfold_ty][ndef][njet][nkappa];
  TH1D *RefoldDataCT[clos_ty][unfold_ty][ndef][njet][nkappa];
  TH1D *UnfoldDataCT_HT[clos_ty][unfold_ty][ndef][njet][nkappa][nHLTmx];
  TH1D *RefoldDataCT_HT[clos_ty][unfold_ty][ndef][njet][nkappa][nHLTmx];

  TH1D *Unfold[ndef][njet][nkappa];
  TH2D *Unfold2d[ndef][njet][nkappa];
  TH1D *Refold[ndef][njet][nkappa];
  TH1D *Unfold_HT[ndef][njet][nkappa][nHLTmx];
  TH1D *Refold_HT[ndef][njet][nkappa][nHLTmx];
  //HT2 bin
  TH1D *Data_recoHT[ndef][njet][nkappa][nHLTmx];
  TH1D* MC_genHT[nmc][ndef][njet][nkappa][nHLTmx];
  TH1D* MC_recoHT[nmc][ndef][njet][nkappa][nHLTmx];
  TH1D* MC_recofakeHT[nmc][ndef][njet][nkappa][nHLTmx];
  TH1D* MC_genmissHT[nmc][ndef][njet][nkappa][nHLTmx];
  
  TH1D *hist_effHT[nmc][ndef][njet][nkappa][nHLTmx];
  TH1D *hist_fakeHT[nmc][ndef][njet][nkappa][nHLTmx];
  TH1D *hist_purityHT[nmc][ndef][njet][nkappa][nHLTmx];
  TH1D *hist_stblHT[nmc][ndef][njet][nkappa][nHLTmx];
  
  TH1D *BLTDet[ndef][njet][nkappa][nHLTmx];
  TH1D *BLTGen[ndef][njet][nkappa][nHLTmx];
#ifdef CLOUSER  
  TH1D *Psudo_Data_gen[ndef][njet][nkappa];  
  TH1D *Psudo_Data_genHT[ndef][njet][nkappa][nHLTmx];  
  TH1D *Psudo_Data_reco_fake[ndef][njet][nkappa];  
  TH1D *Psudo_Data_gen_miss[ndef][njet][nkappa];
  TH1D *Psudo_Data_reco_fakeHT[ndef][njet][nkappa];
#endif
  TH2D *CorrCT[clos_ty][unfold_ty][ndef][njet][nkappa];
  TH2D *ProbCT[clos_ty][unfold_ty][ndef][njet][nkappa];
  TH2D *EmatrixCT[clos_ty][unfold_ty][ndef][njet][nkappa];

  TH2D *Corr[ndef][njet][nkappa];
  TH2D *Prob[ndef][njet][nkappa];
  TH2D *Ematrix[ndef][njet][nkappa];
//----------------------------------------------------
	for(int id=0; id<ndef; id++){
        	for(int ij=0; ij<njet; ij++){
                	for(int ik =0 ; ik<nkappa ; ik++){
		Data_reco[id][ij][ik] = (TH1D*)ReadHist1D(Datadir+"/dd_reco_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_eta0",Unfoldroot);
          	Data_reco2d[id][ij][ik] = (TH2D*)ReadHist2D(Datadir+"/Edd_reco_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_eta0",Unfoldroot);
          	for(int ipt=0; ipt < nHLTmx; ipt++){
          	Data_recoHT[id][ij][ik][ipt] = (TH1D*)ReadHist1D(Datadir+"/Edd_reco_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_eta0",Unfoldroot);
               	}
#ifdef CLOUSER
          	Psudo_Data_gen[id][ij][ik] =(TH1D*)ReadHist1D(Datadir+"/dd_gen_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_eta0",Unfoldroot); 
	  	for(int ipt=0; ipt < nHLTmx; ipt++){
          	Psudo_Data_genHT[id][ij][ik][ipt] = (TH1D*)ReadHist1D(Datadir+"/Edd_reco_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_eta0",Unfoldroot);
               	}
         	//Psudo_Data_reco_fake[ity][ivar] =(TH1D*)ReadHist1D(Datadir+"/dd_DataRecominusfake"+ to_string(ity)+"_eta0_"+to_string(var[ivar]),Unfoldroot);
         	//for(int ipt=0; ipt < nHLTmx; ipt++){
         	//Psudo_Data_genHT[ity][ivar][ipt] = (TH1D*)ReadHist1D(Datadir+"/Edd_reco_typ_"+ to_string(ity)+"_eta0_"+to_string(var[ivar])+"_pt"+to_string(ipt),Unfoldroot);
//              }
#endif
    			}
		}
	}
//----------------------------------------------------
  	for(int  imc =0; imc < nmc ; imc++){
		for(int id=0; id<ndef; id++){
                	for(int ij=0; ij<njet; ij++){
                        	for(int ik =0 ; ik<nkappa ; ik++){
	  	MC_gen[imc][id][ij][ik] = (TH1D*)ReadHist1D(mcdir[imc]+"/dd_gen_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_eta0",Unfoldroot);  
	  	MC_genM[imc][id][ij][ik] = (TH1D*)ReadHist1D(mcdir[imc]+"/Edd_gen_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_eta0_px",Unfoldroot);  
	  	MC_gen2d[imc][id][ij][ik] = (TH2D*)ReadHist2D(mcdir[imc]+"/Edd_gen_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_eta0",Unfoldroot);  
	  	MC_reco[imc][id][ij][ik] = (TH1D*)ReadHist1D(mcdir[imc]+"/dd_reco_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_eta0",Unfoldroot);
	  	MC_reco2d[imc][id][ij][ik] = (TH2D*)ReadHist2D(mcdir[imc]+"/Edd_reco_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_eta0",Unfoldroot);
          	MC_genmiss[imc][id][ij][ik] = (TH1D*)ReadHist1D(mcdir[imc]+"/dd_Genminusmiss_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_eta0",Unfoldroot);
          	MC_recofake[imc][id][ij][ik] = (TH1D*)ReadHist1D(mcdir[imc]+"/dd_Recominusfake_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_eta0",Unfoldroot);
          
	  	hist_eff[imc][id][ij][ik] = (TH1D*)ReadHist1D(mcdir[imc]+"/dd_miss_rate_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_eta0",Unfoldroot);
          	hist_fake[imc][id][ij][ik] = (TH1D*)ReadHist1D(mcdir[imc]+"/dd_fake_rate_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_eta0",Unfoldroot);
          	hist_purity[imc][id][ij][ik] = (TH1D*)ReadHist1D(mcdir[imc]+"/dd_Purity_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_eta0",Unfoldroot);
          	hist_stbl[imc][id][ij][ik] = (TH1D*)ReadHist1D(mcdir[imc]+"/dd_stability_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_eta0",Unfoldroot);
	  	MC_Res[imc][id][ij][ik] = (TH2D*)ReadHist2D(mcdir[imc]+"/dd_corr_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_eta0",Unfoldroot);

         	for(int ipt=0; ipt < nHLTmx; ipt++){
          	MC_genHT[imc][id][ij][ik][ipt] = (TH1D*)ReadHist1D(mcdir[imc]+"/Edd_gen_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_eta0_pt"+to_string(ipt),Unfoldroot);  
          	MC_recoHT[imc][id][ij][ik][ipt] = (TH1D*)ReadHist1D(mcdir[imc]+"/Edd_reco_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_eta0_pt"+to_string(ipt),Unfoldroot);
          	MC_genmissHT[imc][id][ij][ik][ipt] = (TH1D*)ReadHist1D(mcdir[imc]+"/Edd_Genminusmiss_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_eta0_pt"+to_string(ipt),Unfoldroot);
          	MC_recofakeHT[imc][id][ij][ik][ipt] = (TH1D*)ReadHist1D(mcdir[imc]+"/Edd_Recominusfake_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_eta0_pt"+to_string(ipt),Unfoldroot);
          
          	hist_effHT[imc][id][ij][ik][ipt] = (TH1D*)ReadHist1D(mcdir[imc]+"/Edd_miss_rate_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_eta0_pt"+to_string(ipt),Unfoldroot);
          	hist_fakeHT[imc][id][ij][ik][ipt] = (TH1D*)ReadHist1D(mcdir[imc]+"/Edd_fake_rate_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_eta0_pt"+to_string(ipt),Unfoldroot);
          	hist_purityHT[imc][id][ij][ik][ipt] = (TH1D*)ReadHist1D(mcdir[imc]+"/Edd_Purity_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_eta0_pt"+to_string(ipt),Unfoldroot);
          	hist_stblHT[imc][id][ij][ik][ipt] = (TH1D*)ReadHist1D(mcdir[imc]+"/Edd_stability_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_eta0_pt"+to_string(ipt),Unfoldroot);
	       			}
      			}
    		}
  	}
}//end of one MCinput root file  reading
//----------------------------------------------------
//Read Unfolded MC for CT
	cout << "Read The CT  Result : " << endl;
  	for(int icl=0; icl < clos_ty; icl++){
  		for(int iun=0; iun < unfold_ty; iun++){
                	for(int id=0; id<ndef; id++){
                        	for(int ij=0; ij<njet; ij++){
                                	for(int ik =0 ; ik<nkappa ; ik++){
	  	UnfoldCT[icl][iun][id][ij][ik] = (TH1D*)ReadHist1D(unfdir+"/dd_TUnfold_NoReg_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_eta0",Unf_root[icl]);
	  	UnfoldCT_M[icl][iun][id][ij][ik] = (TH1D*)ReadHist1D(unfdir+"/Edd_TUnfold_NoReg_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_eta0_px",Unf_root[icl]);
	  	RefoldCT[icl][iun][id][ij][ik] = (TH1D*)ReadHist1D(unfdir+"/dd_Refold_NoReg_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_eta0",Unf_root[icl]);
          	CorrCT[icl][iun][id][ij][ik] = (TH2D*)ReadHist2D(unfdir+"/dd_corr_NoReg_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_eta0",Unf_root[icl]);
          	EmatrixCT[icl][iun][id][ij][ik] = (TH2D*)ReadHist2D(unfdir+"/dd_corr_NoReg_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_eta0",Unf_root[icl]);
          	ProbCT[icl][iun][id][ij][ik] = (TH2D*)ReadHist2D(unfdir+"/dd_Prob_NoReg_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_eta0",Unf_root[icl]);
   
          	for(int ipt=0; ipt < nHLTmx; ipt++){
          	UnfoldCT_HT[icl][iun][id][ij][ik][ipt] = (TH1D*)ReadHist1D(unfdir+"/Edd_TUnfold_NoReg_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_eta0_pt"+to_string(ipt),Unf_root[icl]);
          	RefoldCT_HT[icl][iun][id][ij][ik][ipt] = (TH1D*)ReadHist1D(unfdir+"/Edd_Refold_NoReg_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_eta0_pt"+to_string(ipt),Unf_root[icl]);
	     				}
	  			}
         		}
      		}
	}
}
//----------------------------------------------------
//Unfold data with different MC RM
	for(int icl=0; icl < 3; icl++){
  		for(int iun=0; iun < unfold_ty; iun++){
    			for(int id=0; id<ndef; id++){
                                for(int ij=0; ij<njet; ij++){
                                        for(int ik =0 ; ik<nkappa ; ik++){
          UnfoldDataCT[icl][iun][id][ij][ik] = (TH1D*)ReadHist1D(unfdir+"/dd_TUnfold_NoReg_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_eta0",Unf_dataCT[icl]);
          UnfoldDataCT_M[icl][iun][id][ij][ik] = (TH1D*)ReadHist1D(unfdir+"/Edd_TUnfold_NoReg_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_eta0_px",Unf_dataCT[icl]);
          RefoldDataCT[icl][iun][id][ij][ik] = (TH1D*)ReadHist1D(unfdir+"/dd_Refold_NoReg_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_eta0",Unf_dataCT[icl]);

          for(int ipt=0; ipt < nHLTmx; ipt++){
          UnfoldDataCT_HT[icl][iun][id][ij][ik][ipt] = (TH1D*)ReadHist1D(unfdir+"/Edd_TUnfold_NoReg_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_eta0_pt"+to_string(ipt),Unf_dataCT[icl]);
          RefoldDataCT_HT[icl][iun][id][ij][ik][ipt] = (TH1D*)ReadHist1D(unfdir+"/Edd_Refold_NoReg_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_eta0_pt"+to_string(ipt),Unf_dataCT[icl]);
             				}
          			}
         		}
      		}
	}
}
//----------------------------------------------------
//Unfolded data
	for(int id=0; id<ndef; id++){
        	for(int ij=0; ij<njet; ij++){
                	for(int ik =0 ; ik<nkappa ; ik++){
          	Unfold[id][ij][ik] = (TH1D*)ReadHist1D(unfdir+"/dd_TUnfold_NoReg_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_eta0", Unfoldroot);
          	Unfold2d[id][ij][ik] = (TH2D*)ReadHist2D(unfdir+"/Edd_TUnfold_NoReg_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_eta0", Unfoldroot);
          	Refold[id][ij][ik] = (TH1D*)ReadHist1D(unfdir+"/dd_Refold_NoReg_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_eta0", Unfoldroot);
          	Corr[id][ij][ik] = (TH2D*)ReadHist2D(unfdir+"/dd_corr_NoReg_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_eta0", Unfoldroot);
          	Ematrix[id][ij][ik] = (TH2D*)ReadHist2D(unfdir+"/dd_corr_NoReg_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_eta0", Unfoldroot);
          	Prob[id][ij][ik] = (TH2D*)ReadHist2D(unfdir+"/dd_Prob_NoReg_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_eta0", Unfoldroot);

          	for(int ipt=0; ipt < nHLTmx; ipt++){
          	Unfold_HT[id][ij][ik][ipt] = (TH1D*)ReadHist1D(unfdir+"/Edd_TUnfold_NoReg_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_eta0_pt"+to_string(ipt), Unfoldroot);
          	Refold_HT[id][ij][ik][ipt] = (TH1D*)ReadHist1D(unfdir+"/Edd_Refold_NoReg_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_eta0_pt"+to_string(ipt), Unfoldroot);
          
	  	BLTDet[id][ij][ik][ipt] = (TH1D*)ReadHist1D(unfdir+"/EBLTDet_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_eta0_pt"+to_string(ipt), Unfoldroot);
          	BLTGen[id][ij][ik][ipt] = (TH1D*)ReadHist1D(unfdir+"/EBLTgen_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_eta0_pt"+to_string(ipt), Unfoldroot);
            		}
          	}
        }
}
cout <<" Read histogram OK " << endl;
//----------------------------------------------------  
//Reco comparison 
  	for(int id=0; id<ndef; id++){
        	for(int ij=0; ij<njet; ij++){
                	for(int ik =0 ; ik<nkappa ; ik++){
       				for(int ipt=0; ipt < nHLTmx; ipt++){
        	TCanvas *cpt0 = new TCanvas("cpt0", "canvas0", 600,600 );   
    		TH1D *MyHist  = (TH1D*)Data_recoHT[id][ij][ik][ipt]->Clone();  Integralhist(MyHist);	//divBybinWidth(MyHist);
		Myplotset(MyHist,0,0);
	
		MyHist->SetTitle(Form("%s_{%s}^{%s}: %s  ",obs_def[id], jet_num[ij], k_fact[ik], htrang[ipt])); // check ??
		MyHist->GetXaxis()->SetTitle("");
		MyHist->GetYaxis()->SetTitle(Form(" %s %s_{%s}^{%s}" ,obs_logy[ik], obs_def[id], njets[ij], k_fact[ik]));  // check ??
	
		TH1D *MC_input[nmc];
		const char *MCinput_index[nmc+1];
		const char *data_index[1];
		for(int iout = 0 ; iout < nmc ; iout++){
	  		MC_input[iout] = (TH1D*) MC_recoHT[iout][id][ij][ik][ipt]->Clone();  Integralhist(MC_input[iout]);	 //divBybinWidth(MC_input[iout]);
	  		MCinput_index[iout]= mcnamerco[iout]; }
	  	data_index[0]= DataEra[1]; 
	
		char lplot_xtitle[100];
		sprintf(lplot_xtitle," %s_{%s}^{%s}", obs_def[id], njets[ij], k_fact[ik]);
		//float ratio_range1[2]={1.2,0.9};
		int num1[2]={nmc,1} ;
		float lpos1[7] ={.32,0.2,0.55,0.38, .04, 1.5,0.7};
	
		cpt0->cd();
		SetMycanvas(cpt0,0,0,0,0,0);
		cpt0 =(TCanvas*)(ratio_can(num1, lpos1, MyHist, MC_input, lplot_xtitle, MCinput_index,data_index));
		CMS_lumi( cpt0, iPeriod, iPos ); cpt0->Update();
	
		sprintf(pdfname, "%sRecoJC.pdf(","dd_"); sprintf(pdfname1, "%sRecoJC.pdf","dd_"); sprintf(pdfname2, "%sRecoJC.pdf)","dd_"); 
		if(id==0 && ij==0 && ik==0 && ipt ==0){cpt0->Print(pdfname,"pdf");
        	}else if(id==2 && ij==1 && ik==9 && ipt==9) {cpt0->Print(pdfname2,"pdf");
       		}else{cpt0->Print(pdfname1,"pdf");};
		cout <<"Test 1"<<endl;
      			}
		}
	}  
}//end of RECO PLOT  
cout << "End of Reco Plot"<<endl;
//----------------------------------------------------
//Response Matrix
	for(int  imc =0; imc < nmc ; imc++){
		for(int id=0; id<ndef; id++){
                	for(int ij=0; ij<njet; ij++){
                        	for(int ik =0 ; ik<nkappa ; ik++){
        	TCanvas *cpt0 = new TCanvas("cpt0", "canvas0", 600,600 );  
		cpt0->cd(); 
		SetMycanvas(cpt0,0,0.1,0.15,0.05,0.1);
        	char lplot_xtitle[100]; char lplot_ytitle[100];
        	sprintf(lplot_xtitle, "Reco bin ID ( %s_{%s}^{%s}, P_{T})",obs_def[id], njets[ij], k_fact[ik]); 
		sprintf(lplot_ytitle, "Gen  bin ID ( %s_{%s}^{%s}, P_{T}) ",obs_def[id], njets[ij], k_fact[ik]);	
		double titoff1[3]={1.2,1.3,1.0};
        	double titsize1[3] ={0.035,0.035,0.035};
	
        	Set2dHist( MC_Res[imc][id][ij][ik],lplot_xtitle, lplot_ytitle,"",titoff1, titsize1);
        	gPad->SetLogz();
		MC_Res[imc][id][ij][ik]->Draw("colz");
         
        	sprintf(Title,"%s_{%s}^{%s}:   ",obs_def[id], njets[ij], k_fact[ik]);
		TLegend *leg1 = new TLegend(0.05,0.7,0.4,0.8);     
		CTLegend(leg1,Modelnm[imc],Title); leg1->AddEntry((TObject*)0,obs_def[id] , "");leg1->SetTextColor(-8);leg1->Draw();
        	CMS_lumi( cpt0, iPeriod, iPos ); cpt0->Update();
		sprintf(pdfname, "%sResponse_Mat%i.pdf(","dd_",imc); sprintf(pdfname1, "%sResponse_Mat%i.pdf","dd_",imc); sprintf(pdfname2, "%sResponse_Mat%i.pdf)","dd_",imc);
        	if(id==0 && ij==0 && ik==0 ){cpt0->Print(pdfname,"pdf"); }
		else if(id==2 && ij==1 && ik==9) {cpt0->Print(pdfname2,"pdf"); }
		else{  cpt0->Print(pdfname1,"pdf");};
        		}
      		}
	}
}
cout <<" RECO PLOT OK " << endl;
//---------------------------------------------------- 
//Efficiency, Purity,Fake rate, stability
  	for(int id=0; id<ndef; id++){
        	for(int ij=0; ij<njet; ij++){
                	for(int ik =0 ; ik<nkappa ; ik++){
        	TCanvas *cpt1 = new TCanvas("cpt1", "canvas1", 600,600 );  //for 
        	TCanvas *cpt2 = new TCanvas("cpt2", "canvas2", 600,600 );  //for JC
        	TCanvas *cpt3 = new TCanvas("cpt3", "canvas3", 600,600 );  //for JC
        	TCanvas *cpt4 = new TCanvas("cpt4", "canvas4", 800,800 );  //for 

        	SetMycanvas(cpt1,0,0.1,0.02,0.02,0.12);
        	SetMycanvas(cpt2,0,0.1,0.03,0.05,0.12);
        	SetMycanvas(cpt4,0,0.1,0.03,0.05,0.12);
        	SetMycanvas(cpt3,0,0.1,0.03,0.05,0.12);

   		TLegend *leg2 = new TLegend(0.4,0.6,0.7,0.9);
		CTLegend(leg2," ","");
		TLegend *leg1 = new TLegend(0.1,0.6,0.4,0.8);
        	CTLegend(leg1,"", obs_def[id]); 

       		for(int ipt=0; ipt < nHLTmx; ipt++){
        		cpt1->cd();
			if(ipt<=8){sprintf(Title,"%i < H_{T,2} < %i %s", HT2range[ipt] , HT2range[ipt+1] ,"GeV/c" );}
        		else if(ipt==9){ sprintf(Title," H_{T,2} > %i %s", HT2range[ipt] ,"GeV/c" );}

        		//sprintf(Title,"%s:   ",itypeN[ity]);
        		for (int i = 1; i <= hist_effHT[umc][id][ij][ik][ipt]->GetNbinsX(); ++i) {
         			double content = 1 - hist_effHT[umc][id][ij][ik][ipt]->GetBinContent(i);
         			hist_effHT[umc][id][ij][ik][ipt]->SetBinContent(i, content);
       				}
         		hist_effHT[umc][id][ij][ik][ipt]->SetMinimum(-0.01); hist_effHT[umc][id][ij][ik][ipt]->SetMaximum(1.01);
         		hist_fakeHT[umc][id][ij][ik][ipt]->SetMinimum(-0.01); hist_fakeHT[umc][id][ij][ik][ipt]->SetMaximum(1.01);
         		hist_purityHT[umc][id][ij][ik][ipt]->SetMinimum(-0.01); hist_purityHT[umc][id][ij][ik][ipt]->SetMaximum(1.01);
         		hist_stblHT[umc][id][ij][ik][ipt]->SetMinimum(-0.01); hist_stblHT[umc][id][ij][ik][ipt]->SetMaximum(1.01);

        		Myplotset(hist_effHT[umc][id][ij][ik][ipt],obs_def[id],"Efficiency");
        		hist_effHT[umc][id][ij][ik][ipt]->SetLineColor(color[ipt]);
			hist_effHT[umc][id][ij][ik][ipt]->Draw("same hist "); leg1->Draw();
        		leg2->AddEntry(hist_effHT[umc][id][ij][ik][ipt], Title ,"lp");
			if(ipt==9){leg2->Draw();}
        		//leg2->AddEntry(hist_effHT[0][ity][ivar][ipt], Title ,"lp");

        		cpt2->cd();
       			//SetMycanvas(cpt2);
 			Myplotset(hist_purityHT[umc][id][ij][ik][ipt],obs_def[id],"Purity");
 			hist_purityHT[umc][id][ij][ik][ipt]->SetLineColor(color[ipt]);
        		hist_purityHT[umc][id][ij][ik][ipt]->Draw("same hist"); leg1->Draw();
        		if(ipt==9){leg2->Draw();}
			cpt2->Update();

        		cpt3->cd();
        		Myplotset(hist_fakeHT[umc][id][ij][ik][ipt],obs_def[id],"Fake rate");
        		hist_fakeHT[umc][id][ij][ik][ipt]->SetLineColor(color[ipt]);
        		hist_fakeHT[umc][id][ij][ik][ipt]->Draw("same hist");leg1->Draw();
        		if(ipt==9){leg2->Draw();}
       			cpt3->Update();
		
			cpt4->cd();
        		Myplotset(hist_stblHT[umc][id][ij][ik][ipt],obs_def[id],"Stability");
        		hist_stblHT[umc][id][ij][ik][ipt]->SetLineColor(color[ipt]); 
        		hist_stblHT[umc][id][ij][ik][ipt]->Draw("same hist"); leg1->Draw();
        		if(ipt==9){leg2->Draw();}
			cpt4->Update();
        	}
    		sprintf(pdfname, "%seffi_plot.pdf(","dd_"); sprintf(pdfname1, "%seffi_plot.pdf","dd_"); sprintf(pdfname2, "%seffi_plot.pdf)","dd_");
        	if(id==0 && ij==0 && ik==0){cpt1->Print(pdfname,"pdf");
        	}else if(id==2 && ij==1 && ik==9 ) {cpt1->Print(pdfname2,"pdf");
        	}else{ cpt1->Print(pdfname1,"pdf");};
    
		sprintf(pdfname, "%spuri_plot.pdf(","dd_"); sprintf(pdfname1, "%spuri_plot.pdf","dd_"); sprintf(pdfname2, "%spuri_plot.pdf)","dd_"); 
        	if(id==0 && ij==0 && ik==0 ){cpt2->Print(pdfname,"pdf");
        	}else if(id==2 && ij==1 && ik==9 ) {cpt2->Print(pdfname2,"pdf");
        	}else{cpt2->Print(pdfname1,"pdf");};
     
		sprintf(pdfname, "%sfake_plot.pdf(","dd_"); sprintf(pdfname1, "%sfake_plot.pdf","dd_"); sprintf(pdfname2, "%sfake_plot.pdf)","dd_"); 
        	if(id==0 && ij==0 && ik==0 ){cpt3->Print(pdfname,"pdf");
        	}else if(id==2 && ij==1 && ik==9 ) {cpt3->Print(pdfname2,"pdf");
        	}else{cpt3->Print(pdfname1,"pdf");};
     
		sprintf(pdfname, "%sstab_plot.pdf(","dd_"); sprintf(pdfname1, "%sstab_plot.pdf","dd_"); sprintf(pdfname2, "%sstab_plot.pdf)","dd_"); 
        	if(id==0 && ij==0 && ik==0 ){cpt4->Print(pdfname,"pdf");
        	}else if(id==2 && ij==1 && ik==9 ) {cpt4->Print(pdfname2,"pdf");
        	}else{cpt4->Print(pdfname1,"pdf");};
		cpt1->Clear(); cpt2->Clear(); cpt3->Clear(); cpt4->Clear();
     		}
	}
}   
cout << "Effi, Stabilty, fake, purity " <<endl;
//----------------------------------------------------
//Unfold comparisoin
	for(int iun=0; iun < unfold_ty; iun++){
		for(int id=0; id<ndef; id++){
                	for(int ij=0; ij<njet; ij++){
                        	for(int ik =0 ; ik<nkappa ; ik++){
       		char lplot_xtitle[100];
	 	TCanvas *cpt5 = new TCanvas("cpt5", "canvas5", 600,600 );  //for  
	 	TCanvas *cpt6 = new TCanvas("cpt6", "canvas6", 600,600 );  //for  
	 	TCanvas *cpt7 = new TCanvas("cpt7", "canvas7", 600,600 );  //for  
        	for(int ipt=0; ipt < nHLTmx; ipt++){
	 		TCanvas *cpt0 = new TCanvas("cpt0", "canvas0", 600,600 );  //for  
	  		TH1D *MyHist  = (TH1D*)Unfold_HT[id][ij][ik][ipt]->Clone();
	  		Integralhist(MyHist);
	  		divBybinWidth(MyHist);
	  		Myplotset(MyHist,0,0);
	  		sprintf(Yaxis,"%s %s_{%s}^{%s}" ,obs_logy[ik], obs_def[id], njets[ij], k_fact[ik]);
          		sprintf(Title,"%s: %s   ",obs_def[id], htrang[ipt]);
	  		MyHist->SetTitle(Title);
	  		MyHist->GetXaxis()->SetTitle("");
	  		MyHist->GetYaxis()->SetTitle(Yaxis);
	  
	  		TH1D *MC_input[nmc];
	  		const char *MCinput_index[nmc], *data_index[1];
	  		for(int iout = 0 ; iout < nmc ; iout++){
	  			MC_input[iout] = (TH1D*) MC_genHT[iout][id][ij][ik][ipt]->Clone();
	  			Integralhist(MC_input[iout]);
	  			divBybinWidth(MC_input[iout]);
	  			MCinput_index[iout]= mcname[iout]; } 
       				
			data_index[0]= Unfoldtype[iun]; 
       			//data_index[0]= UndoldEra[1]; 
       			sprintf(lplot_xtitle, "%s_{%s}^{%s}",obs_def[id], njets[ij], k_fact[ik]);
       			//float ratio_range1[2]={1.2,0.9};
       			int num1[2]={nmc,1} ;
       			float lpos1[7] ={.32,0.2,0.55,0.38, .04, 1.4,0.7};
	  
       			//cpt0->cd();
       			//SetMycanvas(cpt0,0,0,0,0,0);
       			cpt0 =(TCanvas*)(ratio_can(num1, lpos1, MyHist, MC_input, lplot_xtitle,MCinput_index,data_index));
       			CMS_lumi( cpt0, iPeriod, iPos ); cpt0->Update();
	  
       			sprintf(pdfname, "%sTUnfold_plot_%i.pdf(" ,"dd_",iun); sprintf(pdfname1, "%sTUnfold_plot_%i.pdf" ,"dd_",iun);sprintf(pdfname2, "%sTUnfold_plot_%i.pdf)" ,"dd_",iun);
       			if(id==0 && ij==0 && ik==0 && ipt==0){cpt0->Print(pdfname,"pdf");
	 		}else if(id==2 && ij==1 && ik==9 && ipt ==9) {cpt0->Print(pdfname2,"pdf");
	  		}else{cpt0->Print(pdfname,"pdf");};
		} 
       		cpt5->cd();
       		double titoff1[3]={1.2,1.3,1.0};
       		double titsize1[3] ={0.035,0.035,0.035};
       		SetMycanvas(cpt5,0,0.1,0.15,0.05,0.1);
       		gStyle->SetPaintTextFormat( "4.2f");
       		Set2dHist(Corr[id][ij][ik],lplot_xtitle, lplot_xtitle,"correlation coefficients", titoff1, titsize1);
       		Corr[id][ij][ik]->Draw("colz ");
       		//Corr[iun][ity][ivar]->Draw("colz text");
       		sprintf(Yaxis,"%s_{%s}^{%s}",obs_def[id], njets[ij], k_fact[ik]); // check ??
       		cout <<"OK1" <<endl;
       		TLegend *leg1 = new TLegend(0.05,0.6,0.4,0.8);
       		CTLegend(leg1,"Unfolded with Pythia8",Title); leg1->AddEntry((TObject*)0,obs_def[id] , ""); leg1->AddEntry((TObject*)0,Unfoldtype[iun] , "");leg1->SetTextColor(-8);leg1->Draw();
       
       		CMS_lumi( cpt5, iPeriod, iPos ); cpt5->Update();       
       		sprintf(pdfname, "%sTUnfold_corr_%i.pdf(" ,"dd_",iun); sprintf(pdfname1, "%sTUnfold_corr_%i.pdf" ,"dd_",iun);sprintf(pdfname2, "%sTUnfold_corr_%i.pdf)" ,"dd_",iun); 
          	if(id==0 && ij==0 && ik==0 ){cpt5->Print(pdfname,"pdf");
          	}else if(id==2 && ij==1 && ik==9 ) {cpt5->Print(pdfname2,"pdf");
          	}else{cpt5->Print(pdfname,"pdf");};

       		cpt6->cd();
       		SetMycanvas(cpt6,0,0.1,0.15,0.05,0.1);
       		Set2dHist(Prob[id][ij][ik],lplot_xtitle, lplot_xtitle,"Probability",titoff1, titsize1);
       		gPad->SetLogz();Prob[id][ij][ik]->SetMinimum(1e-6); Prob[id][ij][ik]->SetMaximum(1);
       		Prob[id][ij][ik]->Draw("colz"); leg1->Draw();
       		CMS_lumi( cpt6, iPeriod, iPos ); cpt6->Update();
       		sprintf(pdfname, "%sTUnfold_prob_%i.pdf(" ,"dd_",iun); sprintf(pdfname1, "%sTUnfold_prob_%i.pdf" ,"dd_",iun);sprintf(pdfname2, "%sTUnfold_prob_%i.pdf)" ,"dd_",iun); 
          	if(id==0 && ij==0 && ik==0){cpt6->Print(pdfname,"pdf");
          	}else if(id==2 && ij==1 && ik==9) {cpt6->Print(pdfname2,"pdf");
          	}else{cpt6->Print(pdfname,"pdf");};

       		cpt7->cd();
       		SetMycanvas(cpt7,0,0.1,0.15,0.05,0.1);
       		Set2dHist(Ematrix[id][ij][ik],lplot_xtitle, lplot_xtitle,"",titoff1, titsize1);
       		Ematrix[id][ij][ik]->Draw("colz"); leg1->Draw();
       		CMS_lumi( cpt7, iPeriod, iPos ); cpt7->Update();
       		sprintf(pdfname, "%sTUnfold_COV_%i.pdf(" ,"dd_",iun); sprintf(pdfname1, "%sTUnfold_COV_%i.pdf" ,"dd_",iun);sprintf(pdfname2, "%sTUnfold_COV_%i.pdf)" ,"dd_", iun);
          	if(id==0 && ij==0 && ik==0 ){cpt7->Print(pdfname,"pdf");
          	}else if(id==2 && ij==1 && ik==9) {cpt7->Print(pdfname2,"pdf");
          	}else{cpt7->Print(pdfname,"pdf");};
 			}
		}
	}
}//End of Unfolded plot
//----------------------------------------------------
//Closure Test Ratio Plot
	int markersty[4]={3,26,5,6};
	for(int id=0; id<ndef; id++){
        	for(int ij=0; ij<njet; ij++){
                	for(int ik =0 ; ik<nkappa ; ik++){
        			for(int ipt=0; ipt < nHLTmx; ipt++){
     		TCanvas *cptr = new TCanvas("cptr", "cptr", 500,600 );  SetMycanvas(cptr,0,0.1,0.02,0.04,0.12);
     		TLegend *leg1= CTLegendV2(0.3,0.80,0.5,0.95,0.04, "Closure Test", obs_def[id], htrang[ipt]);
     		TLegend *leg2= CTLegendV2(0.6,0.75,0.9,.95,0.04, "", "");
     		cptr->SetGridy(); cptr->cd();
     		TH1D *histgen[clos_ty]; TH1D *unf[clos_ty]; TH1D *ratio[clos_ty];
     
     		for(int icl=0; icl < clos_ty ; icl++){
     			histgen[icl]= (TH1D*) MC_genHT[icl][id][ij][ik][ipt]->Clone();
     			unf[icl]= (TH1D*) UnfoldCT_HT[icl][0][id][ij][ik][ipt]->Clone();
    			Integralhist(histgen[icl]);     Integralhist(unf[icl]);
     			gStyle->SetPaintTextFormat("2.5f"); 
     			ratio[icl] = (TH1D*)histgen[icl]->Clone();   ratio[icl]->Divide(ratio[icl], unf[icl], 1, 1, "b");
     			//xtilte, ytitle, max, min, linewidth, SetLineStyle, SetMarkerStyle, SetMarkerSize, SetLineColor
     			MyplotsetV2(ratio[icl], obs_def[id], "MC/Unfolded", 1.5, 0.5, 3, 1, 0, 2, color[icl]); // check ??
     			ratio[icl]->Draw("Same p "); //ratio[icl]->Draw("Same p text50");
     			leg2->AddEntry(ratio[icl],closuretype[icl],"lp");
     			}
     		leg2->Draw();leg1->Draw();
     		CMS_lumi( cptr, iPeriod, iPos ); cptr->Update();
     		sprintf(pdfname, "CTratio_plot.pdf("); sprintf(pdfname1, "CTratio_plot.pdf");sprintf(pdfname2, "CTratio_plot.pdf)");
     		if(id==0 && ij==0 && ik==0 && ipt==0){cptr->Print(pdfname,"pdf");}
		else if(id==2 && ij==1 && ik==9 && ipt==9) {cptr->Print(pdfname2,"pdf"); }
		else{cptr->Print(pdfname,"pdf");};
     		cptr->Clear();
				}
        		}
      		}
   	}
//---------------------------------------------------
//The ratio plot for merged HT2
	for(int id=0; id<ndef; id++){
                for(int ij=0; ij<njet; ij++){
                        for(int ik =0 ; ik<nkappa ; ik++){
     		TCanvas *cptr = new TCanvas("cptr", "cptr", 1500,1000 );  
     		TLegend *leg1= CTLegendV2(0.1,0.80,0.4,0.95,0.03, "Closure Test", obs_def[id]);
     		TLegend *leg2= CTLegendV2(0.4,0.78,0.7,1.0,0.03, "", "");
     		SetMycanvas(cptr,0,0.1,0.01,0.05,0.12); cptr->cd(); cptr->SetGridy();
     		TH1D *histgen[clos_ty]; TH1D *unf[clos_ty]; TH1D *ratio[clos_ty];
		//int icl=;
     		for(int icl=0; icl < clos_ty ; icl++){
     			histgen[icl]= (TH1D*) MC_genM[icl][id][ij][ik]->Clone();
     			unf[icl]= (TH1D*) UnfoldCT_M[icl][0][id][ij][ik]->Clone();
    			Integralhist(histgen[icl]);     Integralhist(unf[icl]);

     			ratio[icl] = (TH1D*)histgen[icl]->Clone();   ratio[icl]->Divide(ratio[icl], unf[icl], 1, 1, "b");
     			//xtilte, ytitle, max, min, linewidth, SetLineStyle, SetMarkerStyle, SetMarkerSize, SetLineColor
     			MyplotsetV2(ratio[icl], obs_def[id], "MC/Unfolded", 1.15, 0.85, 2, 1, markersty[icl], 3, color[icl]);
     			ratio[icl]->Draw("Same hist e1");
     			leg2->AddEntry(ratio[icl],closuretype[icl],"lp");
     			}
     		leg2->Draw();leg1->Draw();
     		CMS_lumi( cptr, iPeriod, iPos ); cptr->Update();
     		sprintf(pdfname, "CTratioV3_plot.pdf("); sprintf(pdfname1, "CTratioV3_plot.pdf");sprintf(pdfname2, "CTratioV3_plot.pdf)");
     		if(id==0 && ij==0 && ik==0 ){cptr->Print(pdfname,"pdf"); }
		else if(id==2 && ij==1 && ik==9) {cptr->Print(pdfname2,"pdf");}
		else{cptr->Print(pdfname,"pdf");};
     		cptr->Clear();
			}
      		}
   	}	
//---------------------------------------------------
//closure Test V2
	for(int id=0; id<ndef; id++){
                for(int ij=0; ij<njet; ij++){
                        for(int ik =0 ; ik<nkappa ; ik++){
     		TCanvas *cptr = new TCanvas("cptr", "cptr", 1500,1000 );  //for
     		TLegend *leg1= CTLegendV2(0.1,0.80,0.4,0.95,0.03, "Closure Test", obs_def[id]);
     		TLegend *leg2= CTLegendV2(0.4,0.78,0.7,1.0,0.03, "", "");
     		SetMycanvas(cptr,0,0.1,0.01,0.05,0.12); cptr->cd(); cptr->SetGridy(); 
     		TH1D *histgen[clos_ty]; TH1D *unf[clos_ty]; TH1D *ratio[clos_ty];
		//int icl=;
     		for(int icl=0; icl < clos_ty ; icl++){
     		histgen[icl]= (TH1D*) MC_gen[icl][id][ij][ik]->Clone();
     		unf[icl]= (TH1D*) UnfoldCT[icl][0][id][ij][ik]->Clone();
     		Integralhist(histgen[icl]);     Integralhist(unf[icl]);

     		ratio[icl] = (TH1D*)histgen[icl]->Clone();   ratio[icl]->Divide(ratio[icl], unf[icl], 1, 1, "b");
     		//xtilte, ytitle, max, min, linewidth, SetLineStyle, SetMarkerStyle, SetMarkerSize, SetLineColor
     		MyplotsetV2(ratio[icl], obs_def[id], "MC/Unfolded", 1.15, 0.85, 2, 1, markersty[icl], 3, color[icl]);
     		ratio[icl]->Draw("Same hist e1");
     		leg2->AddEntry(ratio[icl],closuretype[icl],"lp");
     		}
     		leg2->Draw();leg1->Draw();
     		CMS_lumi( cptr, iPeriod, iPos ); cptr->Update();
     		sprintf(pdfname, "CTratioV2_plot.pdf("); sprintf(pdfname1, "CTratioV2_plot.pdf");sprintf(pdfname2, "CTratioV2_plot.pdf)");
     		if(id==0 && ij==0 && ik==0 ){cptr->Print(pdfname,"pdf"); }
		else if(id==2 && ij==1 && ik==9) {cptr->Print(pdfname2,"pdf");}
		else{cptr->Print(pdfname,"pdf");};
     		cptr->Clear();
			}
      		}
   	}
//---------------------------------------------------
//Closure all plot 
	for(int icl=0; icl < clos_ty ; icl++){
  		for(int id=0; id<ndef; id++){
                	for(int ij=0; ij<njet; ij++){
                        	for(int ik =0 ; ik<nkappa ; ik++){
         	char lplot_xtitle[100];
	 	int iun =0;
         	for(int ipt=0; ipt < nHLTmx; ipt++){
         		TCanvas *cpt0 = new TCanvas("cpt0", "canvas0", 600,600 );  //for  
          		TH1D *MyHist  = (TH1D*)UnfoldCT_HT[icl][0][id][ij][ik][ipt]->Clone();; 
          		Integralhist(MyHist);          
			divBybinWidth(MyHist);
          		Myplotset(MyHist,0,0);
          		sprintf(Yaxis,"%s_{%s}^{%s}",obs_def[id], njets[ij], k_fact[ik]);
          		sprintf(Title,"%s: %s   ",obs_def[id], htrang[ipt]);
          		MyHist->SetTitle(Title);          
			MyHist->GetXaxis()->SetTitle("");          
			MyHist->GetYaxis()->SetTitle(Yaxis);

          		TH1D *MC_input[1];
          		const char *MCinput_index[1], *data_index[1];
          		MC_input[0] = (TH1D*) MC_genHT[icl][id][ij][ik][ipt]->Clone();
          		Integralhist(MC_input[0]);          
			divBybinWidth(MC_input[0]);
          		MCinput_index[0] = mcname[icl]; 
        
       			data_index[0]= Unfoldtype[icl];
       			//data_index[0]= UndoldEra[1]; 
       			sprintf(lplot_xtitle, "%s_{%s}^{%s}",obs_def[id], njets[ij], k_fact[ik]);
       			//float ratio_range1[2]={1.2,0.9};
       			int num1[2]={1,1} ;
       			float lpos1[7] ={.32,0.2,0.55,0.38, .04, 1.4,0.7};

       			//cpt0->cd();
       			//SetMycanvas(cpt0,0,0,0,0,0);
       			cpt0 =(TCanvas*)(ratio_can(num1, lpos1, MyHist, MC_input, lplot_xtitle,MCinput_index,data_index));
       			CMS_lumi( cpt0, iPeriod, iPos ); cpt0->Update();

       			sprintf(pdfname, "%i_%sTUnfold_plot_%i.pdf(" ,icl,"dd_",iun); sprintf(pdfname1, "%i_%sTUnfold_plot_%i.pdf" ,icl,"dd_",iun);sprintf(pdfname2, "%i_%sTUnfold_plot_%i.pdf)" ,icl,"dd_",iun);
       			if(id==0 && ij==0 && ik==0 && ipt==0){cpt0->Print(pdfname,"pdf");
         		}else if(id==2 && ij==1 && ik==9 && ipt ==9) {cpt0->Print(pdfname2,"pdf");
          		}else{cpt0->Print(pdfname,"pdf");};
        			}
      			}
    		}
	}
}//End of Unfolded plot
//---------------------------------------------------
//Closure same data with different MC
	for(int id=0; id<ndef; id++){
        	for(int ij=0; ij<njet; ij++){
                	for(int ik =0 ; ik<nkappa ; ik++){
         	char lplot_xtitle[100];
         	int iun =0;
         	for(int ipt=0; ipt < nHLTmx; ipt++){
         		TCanvas *cpt0 = new TCanvas("cpt0", "canvas0", 600,600 );  
          		TH1D *MyHist  = (TH1D*)UnfoldDataCT_HT[0][0][id][ij][ik][ipt]->Clone();;
          		Integralhist(MyHist);          
			divBybinWidth(MyHist);
          		Myplotset(MyHist,0,0);
          		sprintf(Yaxis,"%s_{%s}^{%s}",obs_def[id], njets[ij], k_fact[ik]);
          		sprintf(Title,"%s: %s   ",obs_def[id], htrang[ipt]);
          		MyHist->SetTitle(Title);          
			MyHist->GetXaxis()->SetTitle("");          
			MyHist->GetYaxis()->SetTitle(Yaxis);

          		TH1D *MC_input[2];
          		const char *MCinput_index[2], *data_index[1];
			for(int icl=0; icl < 2 ; icl++){
  	  			MC_input[icl] = (TH1D*) UnfoldDataCT_HT[icl+1][0][id][ij][ik][ipt]->Clone();
          			Integralhist(MC_input[icl]);          
				divBybinWidth(MC_input[icl]);
          			MCinput_index[icl] = mcname[icl];}

       			data_index[0]= Unfoldtype[iun];
       			//data_index[0]= UndoldEra[1];
       			sprintf(lplot_xtitle,"%s_{%s}^{%s}",obs_def[id], njets[ij], k_fact[ik]);
       			//float ratio_range1[2]={1.2,0.9};
       			int num1[2]={2,1} ;
       			float lpos1[7] ={.32,0.2,0.55,0.38, .04, 1.4,0.7};

       			//cpt0->cd();
       			//SetMycanvas(cpt0,0,0,0,0,0);
       			cpt0 =(TCanvas*)(ratio_can(num1, lpos1, MyHist, MC_input, lplot_xtitle,MCinput_index,data_index));
       			CMS_lumi( cpt0, iPeriod, iPos ); cpt0->Update();

       			sprintf(pdfname, "%sDataCT_plot_%i.pdf(" ,"dd_",iun); sprintf(pdfname1, "%sDataCT_plot_%i.pdf" ,"dd_",iun);sprintf(pdfname2, "%sDataCT_plot_%i.pdf)" ,"dd_",iun);
       			if(id==0 && ij==0 && ik==0 && ipt==0){cpt0->Print(pdfname,"pdf");
         		}else if(id==2 && ij==1 && ik==9 && ipt ==9) {cpt0->Print(pdfname2,"pdf");
          		}else{cpt0->Print(pdfname,"pdf");};
        		}
      		}
	}
}//End of Unfolded plot
//---------------------------------------------------
//Closure Test Refold Ratio Plot
	for(int id=0; id<ndef; id++){
                for(int ij=0; ij<njet; ij++){
                        for(int ik =0 ; ik<nkappa ; ik++){
        			for(int ipt=0; ipt < nHLTmx; ipt++){
     		TCanvas *cptr = new TCanvas("cptr", "cptr", 500,500 ); SetMycanvas(cptr,0,0.12,0.02,0.04,0.12); cptr->SetGridy(); //Horizontal grid
     		TLegend *leg1= CTLegendV2(0.3,0.80,0.5,0.95,0.04, "CT Refold", obs_def[id], htrang[ipt]);
     		TLegend *leg2= CTLegendV2(0.6,0.78,0.9,1.0,0.04, "", "");

     		TH1D *histgen[clos_ty];     
		TH1D *unf[clos_ty];     
		TH1D *ratio[clos_ty];
     		cptr->cd();
     		for(int icl=0; icl < clos_ty-1 ; icl++){
     			histgen[icl]= (TH1D*)MC_recofakeHT[icl][id][ij][ik][ipt]->Clone();
     			unf[icl]= (TH1D*) RefoldCT_HT[icl][0][id][ij][ik][ipt]->Clone();
     			Integralhist(histgen[icl]);     
			Integralhist(unf[icl]);

     			ratio[icl] = (TH1D*)histgen[icl]->Clone();   ratio[icl]->Divide(ratio[icl], unf[icl], 1, 1, "b");
     			MyplotsetV2(ratio[icl], obs_def[id], "MC/Refold", 1.15, 0.85, 3, 1, 0, 1, color[icl]); // check ??
     			ratio[icl]->Draw("Same hist e1");
     			leg2->AddEntry(ratio[icl],closuretype[icl],"lp");
     			}
     		leg2->Draw();leg1->Draw();
     		CMS_lumi( cptr, iPeriod, iPos ); cptr->Update();
     		sprintf(pdfname, "CTRefold_plot.pdf("); sprintf(pdfname1, "CTRefold_plot.pdf");sprintf(pdfname2, "CTRefold_plot.pdf)");
     		if(id==0 && ij==0 && ik==0 && ipt==0){cptr->Print(pdfname,"pdf");}
		else if(id==2 && ij==1 && ik==9 && ipt==9) {cptr->Print(pdfname2,"pdf");}
		else{cptr->Print(pdfname,"pdf");};
     		cptr->Clear();
			}
        	}
	}
}
//---------------------------------------------------
//BLT ratio Plot
	for(int id=0; id<ndef; id++){
                for(int ij=0; ij<njet; ij++){
                        for(int ik =0 ; ik<nkappa ; ik++){
         			for(int ipt=0; ipt < nHLTmx; ipt++){
        	TCanvas *cptr = new TCanvas("cptr", "cptr", 1500,1000 );  //for
        	TLegend *leg1= CTLegendV2(0.1,0.80,0.4,0.95,0.03, "BLT", obs_def[id],htrang[ipt]);
        	TLegend *leg2= CTLegendV2(0.4,0.75,0.7,0.98,0.03, "", "");
        	SetMycanvas(cptr,0,0.1,0.15,0.05,0.12);
        	cptr->SetGridy(); //Horizontal grid
     		TH1D* BLT_reco=(TH1D*)MC_recoHT[0][id][ij][ik][ipt]->Clone();
     		TH1D* Data_reco=(TH1D*)Data_recoHT[id][ij][ik][ipt]->Clone();

     		//BLT_reco->Rebin(2);
     		//Data_reco->Rebin(2);
     		Integralhist(BLT_reco);     
		Integralhist(Data_reco);
     		cout << "Det: "; 
     		Chi2Root(BLT_reco,Data_reco);
     		BLT_reco->Divide(Data_reco);
     
     		TH1D* BLT_gen=(TH1D*)MC_genHT[0][id][ij][ik][ipt]->Clone();
     		TH1D* unf=(TH1D*)Unfold_HT[id][ij][ik][ipt]->Clone();
     		Integralhist(BLT_gen);     
		Integralhist(unf);
     		cout << "Gen: "; 
     		Chi2Root(BLT_gen,unf);
     		BLT_gen->Divide(unf);

     		cptr->cd();
     		SetMycanvas(cptr,0,0.1,0.05,0.05,0.12);
     		MyplotsetV2(BLT_reco, obs_def[id], "MC/Data", 1.3, 0.7, 2, 1, 26, 3, 1);
     		BLT_reco->Draw("Same e1");
     
     		MyplotsetV2(BLT_gen, obs_def[id], "MC/Data", 1.3, 0.7, 2, 1, 3, 3, 2);
     		BLT_gen->Draw("Same  e1");
    
     		leg2->AddEntry(BLT_reco,"Detector level","lp"); leg2->AddEntry(BLT_gen,"Particle level","lp"); 
     		leg2->Draw();leg1->Draw();
        	sprintf(pdfname, "%sBLT_plot_%i.pdf(", "dd_",1234); sprintf(pdfname1, "%sBLT_plot_%i.pdf" ,"dd_",1234);sprintf(pdfname2, "%sBLT_plot_%i.pdf)","dd_",1234);
          	if(id==0 && ij==0 && ik==0 && ipt==0){cptr->Print(pdfname,"pdf");}
		else if(id==2 && ij==1 && ik==9 && ipt==9) {cptr->Print(pdfname2,"pdf"); }
		else{cptr->Print(pdfname,"pdf");};
         		} 
        	}
	}
}
//---------------------------------------------------
//BLT V2
/*	
	for(int id=0; id<ndef; id++){
                for(int ij=0; ij<njet; ij++){
                        for(int ik =0 ; ik<nkappa ; ik++){
      				for(int ipt = 0 ; ipt <nHLTmx; ipt++){
     		TCanvas *cptr = new TCanvas("cptr", "cptr", 1500,1000 ); 
     		TLegend *leg1= CTLegendV2(0.1,0.80,0.4,0.95,0.03, "BLT", itypeN[ity], htrang[ipt]);
     		TLegend *leg2= CTLegendV2(0.4,0.78,0.7,1.0,0.03, "", "");
     		SetMycanvas(cptr,0,0.1,0.15,0.05,0.12);
     		cptr->SetGridy();  cptr->cd();
     		SetMycanvas(cptr,0,0.1,0.15,0.05,0.12);

     		TH1D* BLTreco = (TH1D*)BLTDet[ity][ivar][ipt]->Clone();
     		TH1D* BLTgen = (TH1D*)BLTGen[ity][ivar][ipt]->Clone();
    
     		MyplotsetV2(BLTreco, Esvlogx[ivar], "MC/Data", 1.3, 0.7, 2, 1, 26, 3, 1);
     		BLTreco->Draw("Same hist e1");
     		MyplotsetV2(BLTreco, Esvlogx[ivar], "MC/Data", 1.3, 0.7, 2, 1, 3, 3, 2);
     		BLTgen->Draw("Same hist e1");

     		leg2->AddEntry(BLTreco,"Detector level","lp"); leg2->AddEntry(BLTgen,"Particle level","lp");
     		leg1->Draw(); leg2->Draw();
     		cptr->Update();
     		sprintf(pdfname, "BLT2.pdf("); sprintf(pdfname1, "BLT2.pdf");sprintf(pdfname2, "BLT2.pdf)");
     		if(ity==0 && ivar==0 && ipt ==0){cptr->Print(pdfname,"pdf");}else if(ity==1 && ivar==4 && ipt==7) {cptr->Print(pdfname2,"pdf"); }else{  cptr->Print(pdfname1,"pdf");};
     		cptr->Clear();
      			}
		}
	}
}
*/
//---------------------------------------------------
//BLT V3
	for(int id=0; id<ndef; id++){
                for(int ij=0; ij<njet; ij++){
                        for(int ik =0 ; ik<nkappa ; ik++){
			//for(int ipt = 0 ; ipt <nHLTmx; ipt++){
     		TCanvas *cptr = new TCanvas("cptr", "cptr", 1500,1000 );  //for
     		TLegend *leg1= CTLegendV2(0.1,0.80,0.4,0.95,0.03, "BLT", obs_def[id]);
     		TLegend *leg2= CTLegendV2(0.4,0.78,0.7,1.0,0.03, "", "");
     		SetMycanvas(cptr,0,0.1,0.15,0.05,0.12);
     		cptr->SetGridy();  cptr->cd();
     		SetMycanvas(cptr,0,0.1,0.15,0.05,0.12);
     		char Axisname[100];
     		sprintf(Axisname,"var_%s[UO]",obs_def[id]);

     		TH1D* BLTreco = Data_reco2d[id][ij][ik]->ProjectionX();    //BLTreco->Rebin(2);
     		TH1D* MCreco = MC_reco2d[0][id][ij][ik]->ProjectionX();  //MCreco->Rebin(2);

     		TH1D* BLTunf = Unfold2d[id][ij][ik]->ProjectionX(); 
     		TH1D* MCgen = MC_gen2d[0][id][ij][ik]->ProjectionX();

     		Integralhist(BLTreco);     
		Integralhist(MCreco);
     		Integralhist(BLTunf);     
		Integralhist(MCgen);
     		BLTreco->Divide(MCreco);
     		BLTunf->Divide(MCgen);
     		MyplotsetV2(BLTreco, obs_def[id], "MC/Data", 1.3, 0.7, 2, 1, 26, 3, 1);
     		BLTreco->Draw("Same hist e1");
     		MyplotsetV2(BLTunf, obs_def[id], "MC/Data", 1.3, 0.7, 2, 1, 3, 3, 2);
     		BLTunf->Draw("Same hist e1");

     		leg2->AddEntry(BLTreco,"Detector level","lp"); leg2->AddEntry(BLTunf,"Particle level","lp");
     		leg1->Draw(); leg2->Draw();
     		cptr->Update();
     		sprintf(pdfname, "BLT2V3.pdf("); sprintf(pdfname1, "BLT2V3.pdf");sprintf(pdfname2, "BLT2V3.pdf)");
     		if(id==0 && ij==0 && ik==0 ){cptr->Print(pdfname,"pdf");}
		else if(id==2 && ij==1 && ik==9) {cptr->Print(pdfname2,"pdf"); }
		else{  cptr->Print(pdfname1,"pdf");};
     		cptr->Clear();
      				//}	
   			}
		}
	}
//---------------------------------------------------  
}//end of main program
//---------------------------------------------------
//--------------------Functions----------------------
//---------------------------------------------------
//Ratio plot function
TCanvas *ratio_can(int Nplot[2],float plegend[7], TH1D* data, TH1D* MC[Nplot[0]], char* lowpadx, const char* modnam[Nplot[0]], const  char* datanm[1]){
	//Nplot[0] = number of MC enetered
  	//Nplot[1] = place 1 if upper part is log scale needed
  	//plegend[0->3] = x1,y1,x2,y2 of the legend of the upper plot
  	//plegend[4]= text size
  	//plegend[5-6]= ratio plot axis range
  	//data = data histogram
  	//MC = monte carlo histogram array
  	//legendN = name of the legends for mC one by one
  	//lowpadx = x axis title of the ratio plot
  
  	TCanvas *canvas =new TCanvas("cptfun", "canvas_fun", 575,600 );
  	canvas->cd();
  	canvas->SetRightMargin(0.02);
  	canvas->SetTopMargin(0.02);
  
  	char ratioXaxis1[100];
  	char MCindex[100];
  	//float ymax;
  	//canvas->SetBottomMargin(0.1); 
  	data->GetYaxis()->SetLabelSize(0.03);
  	data->GetXaxis()->SetLabelSize(0.03);
  	data->GetYaxis()->SetTitleSize(0.045);
  	data->GetYaxis()->SetTitleOffset(1.0);
  	data->GetXaxis()->SetTitleSize(0.055);
  	data->GetXaxis()->SetTitleOffset(0.12);
  	data->GetYaxis()->CenterTitle();
  	data->GetXaxis()->CenterTitle();
  	data->SetLineWidth(2);
  	data->SetMarkerStyle(9);
  	data->SetMarkerSize(.8);
  	int ifont =42;
  	data->GetXaxis()->SetTitleFont(ifont);
  	data->GetYaxis()->SetTitleFont(ifont);     
  	data->SetTitleFont(ifont);
  
  	//ymax = data->GetMaximum();
  	//Divide the histogram with bin width
  
  	TPad *padfun1 = new TPad("padfun1", "padfun1", 0, 0.30, 1.0, 1.0);
  	padfun1->SetBottomMargin(0.01); //Upper and lower plot are joined
  	padfun1->SetTopMargin(0.05);    //Upper and lowd
  	padfun1->SetRightMargin(.04);
  	padfun1->Draw();                //Draw the upper pad: pad1
  	padfun1->cd();
  	data->SetFillColor(kYellow);
  	data->SetFillStyle(1111);
  
  	//create the legend of ht2 range 
  	char MC_HTindex[100];
  	TLegend *HT_range = new TLegend(.15,.06,.55,.1);
  	HT_range->SetFillStyle(0);
  	HT_range->SetBorderSize(0);
  	HT_range->SetTextSize(0.04);
  	sprintf(MC_HTindex,"%s" , data->GetTitle());   //legend for the HT value range
  	HT_range->AddEntry((TObject*)0, MC_HTindex, "" );
  	HT_range->Draw();
  	data->SetTitle("");
  	//end the legend of ht2 range
  	data->Draw("e2");
  	//data->Draw(" ");
  
  	if(Nplot[1]==1){gPad->SetLogy();} //condition for log scale
  	gStyle->SetOptStat(0);
  	gPad->SetTickx();
  	gPad->SetTicky();
  
  	double chi2rat[Nplot[0]];        //for chi2 plot in legend 
  	double chi2Ndfrat[Nplot[0]];     //for chi2 plot in legend
  
  	int color[40] = {2,62,30,46,8,49,1,41,42,30,46,28,29,38,30,12,37,49,9,32,2,3,6,7,8,49,1,41,42,30,46,28,29,38,30,12,37,49,9,32};
  	int style[40]={1,1,1,1,5,6,7,8,9,9,2,3,4,5,6,7,8,9,9,9,1,2,3,4,5,6,7,8,9,9,2,3,4,5,6,7,8,9,9,9};
  	for(int iup =0; iup < Nplot[0] ; iup++){
    		MC[iup]->SetLineStyle(style[iup]);
    		MC[iup]->SetLineColor(color[iup]);
    		MC[iup]->SetLineWidth(2);
    		//MC[iup]->Draw("same hist "); //gPad->SetLogy();    
    		//Addition fo chi2/Ndf  with legend
    		int nn = data->GetNbinsX();
    		Double_t resl[nn] , chi2l;
    		Int_t ndfl ,igoodl;
    		TH1D *chidatal = (TH1D*)data->Clone("chidatal");   //for chi square test
    		TH1D *chiMCl = (TH1D*)MC[iup]->Clone("chiMCl");    //for chi square test
    		chiMCl->Chi2TestX(chidatal, chi2l, ndfl, igoodl, "WU", resl);
    		//chidata->Chi2TestX(chiMC, chi2, ndf, igood, "WU", res);
    		cout << "Ndf value =" << ndfl << "chi2 value=" << chi2l << endl;
    		chi2rat[iup]=chi2l;
    		chi2Ndfrat[iup]=chi2l/ndfl;
    		//end of chi2/Ndf
    		//Divide the bin contents/error by bin width
    
    		/*int tmpnbn = MC[iup]->GetNbinsX();
      			for(int ix=0; ix<tmpnbn; ix++) {
      				double awidth = MC[iup]->GetBinWidth(ix+1);  //tmpwid;
      				MC[iup]->SetBinContent(ix+1, MC[iup]->GetBinContent(ix+1)/awidth);
      				double error = MC[iup]->GetBinError(ix+1);
      				MC[iup]->SetBinError(ix+1, error/awidth);
      			}*/
    		//end Divide the bin contents/error by bin width
    		MC[iup]->Draw("same hist e1 ");
  	}//end of Montecarlo loop
  
  	/*int tmpnbd = data->GetNbinsX();
    	for (int ix=0; ix<tmpnbd; ix++) {
    		double awidth = data->GetBinWidth(ix+1); //tmpwid;
    		data->SetBinContent(ix+1, data->GetBinContent(ix+1)/awidth);
    		double error = data->GetBinError(ix+1);
    		data->SetBinError(ix+1, error/awidth);
    		}*/
  	data->Draw(" same e1");
  	//Maximum Uncertainty with respect to Monash
  	//make it off if not needed
  /*
    	TH1D *MC_inputerr = (TH1F*)data->Clone("MC_inputerr");
    	int MCbinnum = data->GetNbinsX();
    	double rel_err[MCbinnum];
    	for(int ix =0; ix <MCbinnum ; ix++){rel_err[ix]=0.0; }
    	for(int iup =0; iup < Nplot[0] ; iup++){
    	for (int ix=0; ix<MCbinnum; ix++) {
    		double dbin = MC[0]->GetBinContent(ix+1);
    		//double dbin = data->GetBinContent(ix+1);
    		double ebin = MC[iup]->GetBinContent(ix+1);
    		double rel_error=(fabs(dbin-ebin))/dbin;
    		MC_inputerr->SetBinContent(ix+1,1);
    		if(rel_err[ix+1] < rel_error){
    			rel_err[ix+1]=rel_error;
    			MC_inputerr->SetBinError(ix+1,rel_error);}
    		else{MC_inputerr->SetBinError(ix+1,rel_err[ix+1]);}
    		}//loop over bin
    	}//end of number monte carlo
    
    	MC_inputerr->SetFillColorAlpha(30,0.050);
   	MC_inputerr->SetMarkerStyle(1);
  */
  	//end of Uncertainty
  	//data->Draw("same");
  	//TLegend *legendn = new TLegend(.4,0.20,0.62,0.35);
  	TLegend *legendn = new TLegend(plegend[0],plegend[1],plegend[2],plegend[3]);
  	legendn->SetFillStyle(0);
  	legendn->SetBorderSize(0);
  	legendn->SetTextSize(plegend[4]);
  	legendn->SetTextFont(42);
  
  	sprintf(MCindex,"%s" , datanm[0]);                //use if chi2 is not needed in the legend
  	legendn->AddEntry(data, MCindex,"lp");
  	for(int iup =0 ; iup <  Nplot[0] ; iup++){
    		sprintf(MCindex,"%s" , modnam[iup]);      //use if chi2 is not needed in the legend
    		//sprintf(MCindex,"%s-#chi^{2}/NDF: %.2f" ,legendN[iup],chi2Ndfrat[iup]);   //legend with chi2/Ndf value
    		legendn->AddEntry(MC[iup], MCindex ,"lp");
  		}
  	legendn->Draw();
  	HT_range->Draw();
  	//ratio plot pad
  	canvas->cd();          //Go back to the main canvas before defining pad2
  	TPad *padfun2 = new TPad("padfun2", "padfun2", 0, 0.02,1.0, 0.30);
  	padfun2->SetTopMargin(0);
  	padfun2->SetBottomMargin(.4);
  	padfun2->SetRightMargin(.04);
  	padfun2->SetGridy();  //Horizontal grid
  	padfun2->Draw();
  	padfun2->cd();
  	gPad->SetTickx();
  	gPad->SetTicky();
  	//gStyle->SetErrorX(0.5);
  	for(int ilow =0 ; ilow <  Nplot[0] ; ilow++){  //loop for ratio plot
    		TH1D *rh2;
    		//if(data->GetBinContent(1) > 0 || data->GetBinContent(10) > 0) {
    			//if(data->GetBinContent(10) > 0) {
    				rh2 = (TH1D*)MC[ilow]->Clone("rh2"); 
    				//rh2->Sumw2();
    				rh2->Divide(data);     //MC devide by data
    			/* }else{
       				rh2 = (TH1D*)data->Clone("rh2");
       				rh2->Sumw2();
       				rh2->SetLineColor(MC[ilow]->GetLineColor());
       				rh2->SetLineWidth(MC[ilow]->GetLineWidth());
       				rh2->SetLineStyle(MC[ilow]->GetLineStyle());
       				rh2->Divide(MC[ilow]);
       			}*/
    		rh2->SetTitle("");
    		//rh2->SetLineColor(kBlack);
    		//rh2->SetMaximum(lowYrange[0]);  //.. range
    		//rh2->SetMinimum(lowYrange[1]);  //Define Y ..
    		//cout << "max = " << plegend[5] << "  min = " << plegend[6] << endl;   
    		rh2->SetMinimum(plegend[6]);      //Define Y ..
    		rh2->SetMaximum(plegend[5]);      //.. range
    		rh2->SetStats(0);                 //No statistics on lower plot
    		//rh2->Divide(data);              //MC devide by data
    		//rh2->SetMarkerStyle(21);
    		//rh2->Draw(" same e1");
    		rh2->Draw("same hist e1");
    		//MC_inputerr->Draw("same E2");//Draw the Uncertainty in Ratio plot    
    		rh2->GetXaxis()->SetTitleSize(0.13);
    		//rh2->GetXaxis()->SetTitleFont(43);
    		rh2->GetXaxis()->SetTitleOffset(1.15);
    		rh2->GetXaxis()->SetLabelSize(0.1);
    		rh2->GetXaxis()->CenterTitle();
    		sprintf(ratioXaxis1," %s" ,lowpadx);
    		rh2->GetXaxis()->SetTitle(ratioXaxis1);
    
    		rh2->GetYaxis()->SetTitle("MC/Data");
    		rh2->GetYaxis()->CenterTitle();
    		rh2->GetYaxis()->SetNdivisions(505);
    		rh2->GetYaxis()->SetTitleSize(0.12);
    		rh2->GetXaxis()->SetTitleFont(ifont);
    		rh2->GetYaxis()->SetTitleFont(ifont);
    		rh2->GetYaxis()->SetTitleOffset(0.35);
    		//rh2->GetYaxis()->SetLabelFont(1.0); //Absolute font size in pixel (precision 3)
    		rh2->GetYaxis()->SetLabelSize(0.09);
  		}
  	canvas->Update();
	return canvas;
}//end of  ratio plot function
//---------------------------------------------------
void Integralhist(TH1D *hist){ hist->Scale(1/(hist->Integral()));}
//---------------------------------------------------
void divBybinWidth(TH1D *hist){
	int tmpnbn = hist->GetNbinsX();
  	for (int ix=0; ix<tmpnbn; ix++) {
    		double awidth = hist->GetBinWidth(ix+1); //tmpwid;
    		hist->SetBinContent(ix+1, hist->GetBinContent(ix+1)/awidth);
    		double error = hist->GetBinError(ix+1);
    		hist->SetBinError(ix+1, error/awidth);
  		}
	}
//---------------------------------------------------
void Myplotset(TH1D *MyHist, const char* XTitle, const char* YTitle){
	//int ifornt =102;
  	MyHist->SetTitleOffset(0.4);
  	//MyHist->SetTitleFont(ifornt);
  	MyHist->SetTitleSize(0.02);
  	MyHist->SetStats(0);
  
  	MyHist->GetXaxis()->SetLabelSize(0.03);
  	MyHist->GetXaxis()->SetTitleSize(0.045);
  	MyHist->GetXaxis()->SetTitleOffset(1.0);
  	//MyHist->GetXaxis()->SetTitleFont(ifornt);
  	MyHist->GetXaxis()->CenterTitle();
  	MyHist->GetXaxis()->SetTitle(XTitle);
  
  	MyHist->GetYaxis()->SetLabelSize(0.03);
  	MyHist->GetYaxis()->SetTitleSize(0.040);
  	MyHist->GetYaxis()->SetTitleOffset(1.0);
  	MyHist->GetYaxis()->SetTitle(YTitle);
  	//MyHist->GetYaxis()->SetTitleFont(ifornt);     
  	MyHist->GetYaxis()->CenterTitle();
  	MyHist->SetTitle(""); 
  	//gStyle->SetTitleFontSize(.08);
  	MyHist->SetLineWidth(2);
}
//---------------------------------------------------
void SetMycanvas(TCanvas *cpt,double bs,double lm, double rm, double tm,double bm){
	cpt->SetBorderSize(bs);
        cpt->SetLeftMargin(lm);
        cpt->SetRightMargin(rm);
        cpt->SetTopMargin(tm); 
        cpt->SetBottomMargin(bm); 
}
//---------------------------------------------------
void CTLegend(TLegend *legendn, const char* txt1, const char* txt2){ 
  	legendn->SetFillStyle(0);
  	legendn->SetBorderSize(0);
  	legendn->SetTextSize(0.03);
  	legendn->SetTextFont(42);
  	legendn->AddEntry((TObject*)0, txt1, "");
  	legendn->AddEntry((TObject*)0, txt2, "");  
}
//---------------------------------------------------
TLegend* CTLegendV2(float x1, float y1, float x2, float y2, float txtsize, const char* txt1="", const char* txt2="",const char* txt3="",const char* txt4=""){
  	TLegend *leg = new TLegend(x1,y1,x2,y2);
  	leg->SetFillStyle(0);
  	leg->SetBorderSize(0);
  	leg->SetTextSize(txtsize);
  	leg->SetTextFont(42);
  	leg->AddEntry((TObject*)0, txt1, "");
  	leg->AddEntry((TObject*)0, txt2, "");
  	leg->AddEntry((TObject*)0, txt3, "");
  	leg->AddEntry((TObject*)0, txt4, "");
return leg;
}
//---------------------------------------------------
void Set2dHist(TH2D *MyHist, const char* XTitle, const char* YTitle,const char* ZTitle, double titoff[3], double titsize[3] ){
  	MyHist->SetTitleOffset(0.2);
  	//MyHist->SetTitleFont(102);
  	MyHist->SetTitleSize(0.02);

  	MyHist->GetXaxis()->SetLabelSize(0.03);
  	MyHist->GetXaxis()->SetTitleSize(titsize[0]);
  	MyHist->GetXaxis()->SetTitleOffset(titoff[0]);
  	//MyHist->GetXaxis()->SetTitleFont(102);
  	//MyHist->GetXaxis()->CenterTitle();

  	MyHist->GetYaxis()->SetLabelSize(0.03);
  	MyHist->GetYaxis()->SetTitleSize(titsize[1]);
  	MyHist->GetYaxis()->SetTitleOffset(titoff[1]);

  	MyHist->GetZaxis()->SetLabelSize(0.03);
  	MyHist->GetZaxis()->SetTitleSize(titsize[2]);
  	MyHist->GetZaxis()->SetTitleOffset(titoff[2]);

  	MyHist->SetTitle("");
  	MyHist->GetXaxis()->SetTitle(XTitle);
  	MyHist->GetYaxis()->SetTitle(YTitle);
  	MyHist->GetZaxis()->SetTitle(ZTitle);
}
//---------------------------------------------------
TCanvas *ratio_canV2(int Nplot[3],float plegend[7], TH1D* data, TH1D* MC[Nplot[0]], char* lowpadx, const char* modnam[Nplot[0]], const  char* datanm[3]){
	//Nplot[0] = number of MC enetered
  	//Nplot[1] = place 1 if upper part is log scale needed
  	//Nplot[2] = place 1 if MC/Data or 0 if data/MC
  	//plegend[0->3] = x1,y1,x2,y2 of the legend of the upper plot
  	//plegend[4]= text size
  	//plegend[5-6]= ratio plot axis range
  	//data = data histogram
  	//MC = monte carlo histogram array
  	//legendN = name of the legends for mC one by one
  	//lowpadx = x axis title of the ratio plot
  	//datanm[0] = name of data (or MC in case one MC is used)
  	//datanm[1] = name of MC/data Data term
  	//datanm[2] = name of Mc/Data MC term
  	bool isstat =1;
  	TCanvas *canvas =new TCanvas("cptfun", "canvas_fun", 575,600 );
  	canvas->cd();
  	canvas->SetRightMargin(0.02);
  	canvas->SetTopMargin(0.02);
  	
	char ratioXaxis1[100]; char ratioYaxis[100]; char MCindex[100];
  	//float ymax;
  	//canvas->SetBottomMargin(0.1);
  	data->GetYaxis()->SetLabelSize(0.03);
  	data->GetXaxis()->SetLabelSize(0.03);
  	data->GetYaxis()->SetTitleSize(0.045);
  	data->GetYaxis()->SetTitleOffset(1.0);
  	data->GetXaxis()->SetTitleSize(0.055);
  	data->GetXaxis()->SetTitleOffset(0.12);
  	data->GetYaxis()->CenterTitle();
  	data->GetXaxis()->CenterTitle();
  	data->SetLineWidth(2);
  	data->SetMarkerStyle(9);
  	data->SetMarkerSize(.8);
  	int ifont =42;
  	data->GetXaxis()->SetTitleFont(ifont);
  	data->GetYaxis()->SetTitleFont(ifont);
  	data->SetTitleFont(ifont);
  	//ymax = data->GetMaximum();
  	//Divide the histogram with bin width

  	TPad *padfun1 = new TPad("padfun1", "padfun1", 0, 0.30, 1.0, 1.0);
  	padfun1->SetBottomMargin(0.01); //Upper and lower plot are joined
  	padfun1->SetTopMargin(0.05);    //Upper and lowd
  	padfun1->SetRightMargin(.04);
  	padfun1->Draw();                //Draw the upper pad: pad1
  	padfun1->cd();
  	data->SetFillColor(kYellow);
  	data->SetFillStyle(1111);
  	//create the legend of ht2 range
  	char MC_HTindex[100];
  	TLegend *HT_range = new TLegend(.15,.06,.55,.1);
  	HT_range->SetFillStyle(0);
  	HT_range->SetBorderSize(0);
  	HT_range->SetTextSize(0.04);
  	sprintf(MC_HTindex,"%s" , data->GetTitle());   //legend for the HT value range
  	HT_range->AddEntry((TObject*)0, MC_HTindex, "" );
  	HT_range->Draw();
  	data->SetTitle("");
  	//end the legend of ht2 range
  	data->Draw("");

  	if(Nplot[1]==1){gPad->SetLogy();} //condition for log scale
  	gStyle->SetOptStat(0);
  	gPad->SetTickx();
  	gPad->SetTicky();

  	double chi2rat[Nplot[0]];    //for chi2 plot in legend
  	double chi2Ndfrat[Nplot[0]]; //for chi2 plot in legend

  	int color[40] = {2,62,30,46,8,49,1,41,42,30,46,28,29,38,30,12,37,49,9,32,2,3,6,7,8,49,1,41,42,30,46,28,29,38,30,12,37,49,9,32};
  	int style[40]={1,2,3,4,5,6,7,8,9,9,2,3,4,5,6,7,8,9,9,9,1,2,3,4,5,6,7,8,9,9,2,3,4,5,6,7,8,9,9,9};
  	for(int iup =0; iup < Nplot[0] ; iup++){
    		MC[iup]->SetLineStyle(style[iup]);
    		MC[iup]->SetLineColor(color[iup]);
    		MC[iup]->SetLineWidth(2);
    		//MC[iup]->Draw("same hist "); //gPad->SetLogy();
    		//Addition fo chi2/Ndf  with legend
    		int nn = data->GetNbinsX();
    		Double_t resl[nn] , chi2l;
    		Int_t ndfl ,igoodl;
    		TH1D *chidatal = (TH1D*)data->Clone("chidatal");   //for chi square test

    		TH1D *chiMCl = (TH1D*)MC[iup]->Clone("chiMCl");    //for chi square test
    		chiMCl->Chi2TestX(chidatal, chi2l, ndfl, igoodl, "WU", resl);
		//chidata->Chi2TestX(chiMC, chi2, ndf, igood, "WU", res);
 		//cout << "Ndf value =" << ndfl << "chi2 value=" << chi2l << endl;
    		chi2rat[iup]=chi2l;
    		chi2Ndfrat[iup]=chi2l/ndfl;
    		//end of chi2/Ndf
    		//Divide the bin contents/error by bin width

    		/*int tmpnbn = MC[iup]->GetNbinsX();
      		for (int ix=0; ix<tmpnbn; ix++) {
      			double awidth = MC[iup]->GetBinWidth(ix+1);  //tmpwid;
      			MC[iup]->SetBinContent(ix+1, MC[iup]->GetBinContent(ix+1)/awidth);
      			double error = MC[iup]->GetBinError(ix+1);
      			MC[iup]->SetBinError(ix+1, error/awidth);
      		}*/
    		//end Divide the bin contents/error by bin width
		//MC[iup]->Draw("same hist e1 ");
    		MC[iup]->Draw("same hist");
  	}//end of Montecarlo loop

  	/*int tmpnbd = data->GetNbinsX();
    	for (int ix=0; ix<tmpnbd; ix++) {
    		double awidth = data->GetBinWidth(ix+1);             //tmpwid;
    		data->SetBinContent(ix+1, data->GetBinContent(ix+1)/awidth);
    		double error = data->GetBinError(ix+1);
    		data->SetBinError(ix+1, error/awidth);
    	}*/

   	//data->Draw(" same e1");
   	data->Draw(" same ");

 	//Stat Error(DATA)
   	TH1D *staterr = (TH1D*)data->Clone("staterr");
    	staterr->Reset();

    	for(int ierr=1; ierr < staterr->GetNbinsX()+1; ierr++){
     		double bincont = data->GetBinContent(ierr);
     		double serr = data->GetBinError(ierr);
     		double error1 =sqrt(pow((serr/bincont),2.));
    		staterr->SetBinContent(ierr,1.0);
    		staterr->SetBinError(ierr,error1);
    	}
    	gStyle->SetErrorX(0.5);
    	gStyle->IsTransparent();
    	staterr->SetFillColor(5);
    	staterr->SetMarkerSize(0.1);
    	staterr->SetFillStyle(3002);
    	//staterr->SetFillStyle(1111);
    	//staterr->SetFillColor(21);
   	//staterr->SetMarkerSize(0.1);
  	//Maximum Uncertainty with respect to Monash
 	//make it off if not needed
  /*
    	TH1D *MC_inputerr = (TH1F*)data->Clone("MC_inputerr");
    	int MCbinnum = data->GetNbinsX();
    	double rel_err[MCbinnum];
    	for(int ix =0; ix <MCbinnum ; ix++){rel_err[ix]=0.0; }
    	for(int iup =0; iup < Nplot[0] ; iup++){
    		for (int ix=0; ix<MCbinnum; ix++) {
    			double dbin = MC[0]->GetBinContent(ix+1);
    			//double dbin = data->GetBinContent(ix+1);
    			double ebin = MC[iup]->GetBinContent(ix+1);
    			double rel_error=(fabs(dbin-ebin))/dbin;
    			MC_inputerr->SetBinContent(ix+1,1);

    			if(rel_err[ix+1] < rel_error){
    				rel_err[ix+1]=rel_error;
    				MC_inputerr->SetBinError(ix+1,rel_error);}
    			else{MC_inputerr->SetBinError(ix+1,rel_err[ix+1]);}
    		}//loop over bin
    	}//end of number monte carlo

    	MC_inputerr->SetFillColorAlpha(30,0.050);
    	MC_inputerr->SetMarkerStyle(1);
  */
  	//end of Uncertainty

  	//data->Draw("same");
  	//TLegend *legendn = new TLegend(.4,0.20,0.62,0.35);
  	TLegend *legendn = new TLegend(plegend[0],plegend[1],plegend[2],plegend[3]);
  	legendn->SetFillStyle(0);
  	legendn->SetBorderSize(0);
  	legendn->SetTextSize(plegend[4]);
  	legendn->SetTextFont(42);

  	sprintf(MCindex,"%s" , datanm[0]);                //use if chi2 is not needed in the legend
  	legendn->AddEntry(data, MCindex,"lp");
  	for(int iup =0 ; iup <  Nplot[0] ; iup++){
    		sprintf(MCindex,"%s" , modnam[iup]);      //use if chi2 is not needed in the legend
    		//sprintf(MCindex,"%s-#chi^{2}/NDF: %.2f" ,legendN[iup],chi2Ndfrat[iup]);   //legend with chi2/Ndf value
    		legendn->AddEntry(MC[iup], MCindex ,"lp");
  		}
  	if(isstat)legendn->AddEntry(staterr, "Statistical Uncertainty","f");
  		legendn->Draw();
  		HT_range->Draw();
  	//ratio plot pad
  	canvas->cd();          //Go back to the main canvas before defining pad2
  	TPad *padfun2 = new TPad("padfun2", "padfun2", 0, 0.02,1.0, 0.30);
  	padfun2->SetTopMargin(0);
  	padfun2->SetBottomMargin(.4);
  	padfun2->SetRightMargin(.04);
  	padfun2->SetGridy(); //Horizontal grid
  	padfun2->Draw();
  	padfun2->cd();
  	gPad->SetTickx();
  	gPad->SetTicky();
  	//gStyle->SetErrorX(0.5);
  	for(int ilow =0 ; ilow <  Nplot[0] ; ilow++){  //loop for ratio plot
    		TH1D *rh2;
    		//if(data->GetBinContent(1) > 0 ||  data->GetBinContent(10) > 0) {
    			//if( data->GetBinContent(10) > 0) {
    				rh2 = (TH1D*)MC[ilow]->Clone("rh2");
    				//rh2->Sumw2();
    				rh2->Divide(data);     //MC devide by data
    			/* }else{
       				rh2 = (TH1D*)data->Clone("rh2");
       				rh2->Sumw2();
       				rh2->SetLineColor(MC[ilow]->GetLineColor());
       				rh2->SetLineWidth(MC[ilow]->GetLineWidth());
       				rh2->SetLineStyle(MC[ilow]->GetLineStyle());
       				rh2->Divide(MC[ilow]);
       			}*/
    		rh2->SetTitle("");
    		//rh2->SetLineColor(kBlack);
    		//rh2->SetMaximum(lowYrange[0]);  //.. range
    		//rh2->SetMinimum(lowYrange[1]);  //Define Y ..
    		//cout << "max = " << plegend[5] << "  min = " << plegend[6] << endl;
    		rh2->SetMinimum(plegend[6]);      //Define Y ..
    		rh2->SetMaximum(plegend[5]);      //.. range
    		rh2->SetStats(0);                 //No statistics on lower plot
    		//gStyle->SetStatColor(2);        //No statistics on lower plot
    		//gStyle->SetStatX(0.7);;         //No statistics on lower plot
    		//rh2->Divide(data);              //MC devide by data
    		//rh2->SetMarkerStyle(21);
    		//rh2->Draw(" same e1");
    		rh2->Draw("Same P text ");
  
   		//Stat Error
   		TH1D *mstaterr = (TH1D*)MC[ilow]->Clone("staterr");
    		mstaterr->Reset();
    		for(int ierr=1; ierr < mstaterr->GetNbinsX()+1; ierr++){
     			double bincont = MC[ilow]->GetBinContent(ierr);
     			double serr = MC[ilow]->GetBinError(ierr);
     			double error1 =sqrt(pow((serr/bincont),2.));
    			mstaterr->SetBinContent(ierr,1.0);
    			mstaterr->SetBinError(ierr,error1);
    			}
    		gStyle->SetErrorX(0.5);
    		mstaterr->SetFillColor(5);
    		mstaterr->SetMarkerSize(0.1);
    		mstaterr->SetFillStyle(1111);
    		mstaterr->SetFillStyle(3002);
    		//staterr->SetFillColor(21);
    		//staterr->SetMarkerSize(0.1);
    		//Stat Error
    
    		staterr->Draw("E2 Same");
    		//mstaterr->Draw("E2 Same");
    		//MC_inputerr->Draw("same E2");       //Draw the Uncertainty in Ratio plot

    		rh2->GetXaxis()->SetTitleSize(0.13);
    		//rh2->GetXaxis()->SetTitleFont(43);
    		rh2->GetXaxis()->SetTitleOffset(1.15);
    		rh2->GetXaxis()->SetLabelSize(0.1);
    		rh2->GetXaxis()->CenterTitle();
    		sprintf(ratioXaxis1," %s" ,lowpadx);
    		rh2->GetXaxis()->SetTitle(ratioXaxis1);
    		//rh2->GetYaxis()->SetTitle("MC/Data");
    		if(Nplot[2]==1){ sprintf(ratioYaxis,"%s/%s" ,datanm[1],datanm[2]);}
    		if(Nplot[2]==0){ sprintf(ratioYaxis,"%s/%s" ,datanm[2],datanm[1]);}
    		rh2->GetYaxis()->SetTitle(ratioYaxis);
    		rh2->GetYaxis()->CenterTitle();
    		rh2->GetYaxis()->SetNdivisions(505);
    		rh2->GetYaxis()->SetTitleSize(0.1);
    		rh2->GetXaxis()->SetTitleFont(ifont);
    		rh2->GetYaxis()->SetTitleFont(ifont);
    		rh2->GetYaxis()->SetTitleOffset(0.40);
    		//rh2->GetYaxis()->SetLabelFont(1.0); //Absolute font size in pixel (precision 3)
    		rh2->GetYaxis()->SetLabelSize(0.09);
  	}
  	canvas->Update();
	return canvas;
}//end of  ratio plot function
//---------------------------------------------------
void HT2_Normal(TH2D* VarHT, int nht, int htbin[nht+1], TH1D* Var[nht]){
	TH2D* VarHTin = (TH2D*)VarHT->Clone();
	int nbiny= VarHTin->GetNbinsY();
	int nbinx= VarHTin->GetNbinsX();

	TH1D* h_var = VarHTin->ProjectionX(); h_var->Reset();

	for(int iht =1 ; iht <= nbiny; iht++){
   	h_var->Reset();
   		for(int ivar =1 ; ivar<= nbinx  ; ivar++){
   			h_var->SetBinContent(ivar,VarHTin->GetBinContent(ivar,iht));
   			h_var->SetBinError(ivar,VarHTin->GetBinError(ivar,iht));   
           		}
  		Var[iht-1] =(TH1D*)h_var->Clone();
		cout << " iht :"<< iht << endl;
	}
}
//---------------------------------------------------
TH1D* ReadHist1D(string name, TFile* root, int irbin=1){
	TString histname = name; cout << histname<< endl;
	TH1D* hist=(TH1D*)root->Get(histname);
	//hist->Rebin(irbin);
	//hist->Write();
return hist;
}
//---------------------------------------------------
TH2D* ReadHist2D(string name, TFile* root, int irbin=1){
	TString histname = name; cout << histname<<endl;
	TH2D* hist=(TH2D*)root->Get(histname);
	//hist->RebinY(irbin);
	//hist->Write();
return hist;
}
//---------------------------------------------------
void HT2_NormalV2(TH2D* VarHT, int nht, int htbin[nht+1], TH1D* Var[nht]){
	TH2D* VarHTin = (TH2D*)VarHT->Clone();
	int nbiny= VarHTin->GetNbinsY();
	int nbinx= VarHTin->GetNbinsX();

	for (int iht = 1; iht <= nbiny; ++iht) {
       		const char * name_hbin = Form("%s_pt%i",VarHTin->GetName(), iht);
        	TH1D* h_var  = VarHTin->ProjectionX(name_hbin, iht, iht);
        	Var[iht-1] =(TH1D*)h_var->Clone();        
        	//h_var-Write();
      		}
	}
//---------------------------------------------------
//for Ratio plot
//xtilte, ytitle, max, min, linewidth, SetLineStyle, SetMarkerStyle, SetMarkerSize, SetLineColor 
void MyplotsetV2(TH1D *MyHist, const char* XT, const char* YT, float mx, float min, int ilw, int ilsty, int imsty, int imstysize, int icl){  
	MyHist->GetXaxis()->SetTitle(XT); MyHist->GetYaxis()->SetTitle(YT);
  	MyHist->SetMinimum(min); MyHist->SetMaximum(mx);
  	MyHist->SetLineWidth(ilw); MyHist->SetLineStyle(ilsty);
  	MyHist->SetMarkerStyle(imsty); MyHist->SetMarkerSize(imstysize); MyHist->SetLineColor(icl);
  	MyHist->SetTitle("");
  	MyHist->GetXaxis()->CenterTitle();
  	MyHist->GetYaxis()->CenterTitle();

  	MyHist->GetXaxis()->SetLabelSize(0.03);
  	MyHist->GetXaxis()->SetTitleSize(0.045);
  	MyHist->GetXaxis()->SetTitleOffset(1.0);
  	MyHist->GetYaxis()->SetLabelSize(0.03);
  	MyHist->GetYaxis()->SetTitleSize(0.040);
  	MyHist->GetYaxis()->SetTitleOffset(1.0);
}
//---------------------------------------------------
void Chi2Root(TH1 * data, TH1 * MC, int rebin = 1){
	TH1D *chidata = (TH1D*)data->Clone("chidata");
 	TH1D *chiMC = (TH1D*)MC->Clone("chimc");
	chidata->Rebin(rebin);
	chiMC->Rebin(rebin);
	int n = chidata->GetNbinsX();
	Double_t res[n] , chi2;
	Int_t ndf ,igood;
 	cout << " Numbers of bin = " << n << endl;
 	chiMC->Chi2TestX(chidata, chi2, ndf, igood, "WU", res);
	cout << "Root chi2/ndf = " << chi2 << " / " << ndf << " = " << (chi2/ndf) << endl;
}
//---------------------------------------------------
