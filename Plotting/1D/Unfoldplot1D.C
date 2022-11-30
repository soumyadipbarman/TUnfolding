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

//#define BLT
#define CLOUSER

void Unfoldplot1D(){
  
  static const int unfold_ty =1;   //Unfold method to be plots
  static const int clos_ty = 4;    //Number of Closure Test performed
  static const int nHLTmx=10;      //PT Range
  static const int nmc = 4;        //Number of MC sample
  static const int umc = 0;        //0 for Py8, 1 for Py8flat, 2 for MG , 3 for Herwig : which MC have used in Unfolding
    
  static const int ndef=3;         //3 observable
  static const int njet=2;         //2 jets
  static const int nkappa=10;      //10 kappa values

  bool isstat =1;  int irbin=1;
  
  char histname[100],Title[100], Xaxis[100], Yaxis[100], ratioXaxis[100], ratioYaxis[100], pdfname[100], pdfname1[100], pdfname2[100], LegName[100];
  Int_t color[10] ={2,4,6,5,6,46,3,28,38,42};  //define the color for different histograms
  Int_t HT2range[nHLTmx+1]={92, 119, 185, 251, 319, 388, 467, 518, 579, 669, 3000}; //2017
  const int njetetamn=1;         //eta value used 2.5

  string Datadir = "Data";
  string folddir = "Folded";
  string unfdir = "Unfold";
#ifdef CLOUSER
  string mcdir[4]={"Pythia8","Py8Flat","MG5","HW7"}; //Sequence For Closure Test
#else
  string mcdir[4]={"Pythia8","Py8Flat","MG5","HW7"}; //Sequence For Unfold and Refold plot
#endif
  
  const char* obs_def[3]={"Q_{D}","Q_{L}","Q_{T}"};
  const char* jet_num[2]={"Leading-Jet","Sub-Leading-Jet"};
  const char* njets[2]={"1","2"};
  const char* k_fact[10]={"k=0.1","k=0.2","k=0.3","k=0.4","k=0.5","k=0.6","k=0.7","k=0.8","k=0.9","k=1.0"};

  const char* obs_logy[10]={"1/N dN/d","1/N dN/d","1/N dN/d","1/N dN/d","1/N dN/d","1/N dN/d","1/N dN/d","1/N dN/d","1/N dN/d","1/N dN/d"};
  const char* htrang[10]={"92 < P_{T} < 119", "119 < P_{T} < 185", "185 < P_{T} < 251", "251 < P_{T} < 319", "319 < P_{T} < 388", "388 < P_{T} <467", "467 < P_{T} <518","518 < P_{T} < 579", "579 < P_{T} < 669", "P_{T} > 669"};
  const char* Unit[4]={"GeV","fb^{-1}","Pb","#%"};

  //const char* Validity_test[4]={"Closure test","Bottom Line test"," Unfolded","Refold"};
  //const char* h2dMat_name[4]={"Covariance matrix","correlation coefficients"," probabilities matrix ","Response matrix"};
#ifdef CLOUSER
  const char* mcname[4]={"Py8-gen","Py8 Flat-gen","MadGraph-gen","Herwig7-gen"};
#else
  const char* mcname[4]={"Pythia8","Pythia8 Flat","MadGraph","Herwig7"};
  //const char* mcname[4]={"Py8-gen","Py8 Flat-gen","MadGraph-gen","Herwig7-gen"};
#endif
  const char* DataEra[3]={"Data","Data","Data"};
  const char* RunEra[3]={"2016","2017","2018"};
  static const int iera = 1;  // Run year
  int iPeriod = 0;  int iPos=0 ;

  //const char* mcnamerco[4]={"Pythia8 ","MadGraph ","Herwig7","PY8 Flat"};
  const char* mcnamerco[4]={"Pythia8","PY8 Flat","MadGraph","Herwig7"};

  //const char* closuretype[4]={"PY8 by PY8","PY8 Flat By PY8 ", "MG By PY8","HW7 By PY8"};
  //const char* Refoldtype[5]={"Py8 Refold","Py8 Flat Refold ","MadGraph Refold","Herwig7 Refold","Refold Pythia8(Iterative method)"};
#ifdef CLOUSER
  const char* Unfoldtype[5]={"Py8-unf","Py8 Flat-unf","MadGraph-unf","Herwig7-unf ","Others"};
#else
  const char* Unfoldtype[5]={"Unfold by Py8 RM","Unfold by MG RM","TUnfold(Scan Tau)","TUnfold(ScanSURE)","TUnfold(Iterative method)"};
#endif
  const char* Modelnm[4]={"Pythia8","PY8 Flat","MadGraph","Herwig7"};
  //const char* Methodtype[5]={"(No Regularisation) ","(L-Curve scan)","(Scan Tau)","(ScanSURE)","(Iterative method)"};
  //const char* smeared[5]={"TUnfold","Refold","Folded-back","GEN","RECO"};
  const char* data_unf[5]={"Unf by Py8","Unf by MadGraph","Unf by Herwig7","GEN","RECO"};
  const char* Refoldname[5]={"Py8 Refold","Py8 Flat Refold ","MadGraph Refold","Herwig7 Refold","Refold Pythia8(Iterative method)"};
  const char* Recomcname[4]={"Py8-reco","Py8 Flat-reco","MadGraph-reco","Herwig7-reco"};

  const char* dirname[4]={"Pythia8","Py8Flat","MG5","HW7"};

/*
  const char* regN[4]={"Tunfold_Noreg","Tunfold_lscan_typ_","Tunfold_scantau_typ_","Tunfold_SURE_typ_"};
  const char* reg_refold[4]={"Tunfold_NoReg_Refold","Tunfold_lscan_Refold","Tunfold_scantau_Refold","Tunfold_SURE_Refold"};
  const char* CORR[4]={"Tunfold_Noreg_corr","Tunfold_lscan_corr","Tunfold_scantau_corr","Tunfold_SURE_corr"};
  const char* COVN[4]={"Tunfold_Noreg_Emat","Tunfold_lscan_Emat","Tunfold_scantau_Emat","Tunfold_SURE_Emat"};
  const char* ProbN[4]={"Tunfold_Noreg_probM","Tunfold_lscan_probM","Tunfold_scantau_probM","Tunfold_SURE_probM"};

  const char* Validity_test[4]={"Closure test","Bottom Line test"," Unfolded","Refold"};
  const char* h2dMat_name[4]={"Covariance matrix","correlation coefficients"," probabilities matrix ","Response matrix"};
  const char* mcname[3]={"Pythia8 CP5 Tune","Madgraph","Herwig++"};
  const char* mcnamerco[3]={"Pythia8 RECO","Madgraph RECO","Herwig7 RECO"};
  const char* mcnamegen[3]={"Pythia8 Flat GEN","Madgraph GEN","Herwig7 GEN"};
  const char* DataEra[3]={"Data RECO","Data RECO","Data RECO"};
  const char* UndoldEra[3]={"Unfold 2016","Unfold 2017","Unfold 2018"};
  const char* RefoldEra[5]={"Refold Pythia8(No Regularisation) ","Refold Pythia8(L-Curve scan)","Refold Pythia8(Scan Tau)","Refold Pythia8(ScanSURE)","Refold Pythia8(Iterative method)"};
 const char* RefoldEra[5]={"Refold Madgraph(No Regularisation) ","Refold Madgraph(L-Curve scan)","Refold Madgraph(Scan Tau)","Refold Madgraph(ScanSURE)","Refold Pythia8(Iterative method)"};
  const char* Unfoldtype[5]={"TUnfold(No Regularisation) ","TUnfold(L-Curve scan)","TUnfold(Scan Tau)","TUnfold(ScanSURE)","TUnfold(Iterative method)"};
  const char* closuretype[5]={"Unfolded Flat-Pythia8 (No Regularisation) ","Unfolded Pythia8 (L-curve scan)","Unfolded Pythia8(Scan Tau)","Unfolded Pythia8 (ScanSURE)","Unfolded(Iterative method)"};
  const char* closuretype[5]={"Unfolded Madgraph(No Regularisation) ","Unfolded Madgraph (L-curve scan)","Unfolded Madgraph(Scan Tau)","Unfolded Madgraph (ScanSURE)","Unfolded(Iterative method)"};
  const char* Modelnm[3]={"Pythia8","Madgraph","Herwig"};
  const char* Methodtype[5]={"(No Regularisation) ","(L-Curve scan)","(Scan Tau)","(ScanSURE)","(Iterative method)"};
  const char* smeared[5]={"TUnfold","Refold","Folded-back","GEN","RECO"};
  static const int iera = 1;
  int iPeriod = 0;  
  int iPos=10 ;
*/  

  TFile *Unf_root[clos_ty];
  TFile *Unf_dataCT[clos_ty];
//----------------------------------------------------
  TFile *Unfoldroot = TFile::Open("../../../Unfolded/26Nov2022/1D/Unfolded_HW7_flat_PY8_flat.root");//Unfolded data by Py8

  Unf_root[0] = TFile::Open("../../../Unfolded/26Nov2022/1D/Unfolded_HW7_flat_PY8_flat.root");      //Py- by Py8
  Unf_root[1] = TFile::Open("../../../Unfolded/26Nov2022/1D/Unfolded_HW7_flat_PY8_flat.root");      //PyFlat - By Py8
  Unf_root[2] = TFile::Open("../../../Unfolded/26Nov2022/1D/Unfolded_HW7_flat_PY8_flat.root");      //MG by Py8
  Unf_root[3] = TFile::Open("../../../Unfolded/26Nov2022/1D/Unfolded_HW7_flat_PY8_flat.root");      //HW by Py8

  Unf_dataCT[0] = TFile::Open("../../../Unfolded/26Nov2022/1D/Unfolded_HW7_flat_PY8_flat.root");    //data by Py8
  Unf_dataCT[1] = TFile::Open("../../../Unfolded/26Nov2022/1D/Unfolded_HW7_flat_PY8_flat.root");    //data by Py8flat
  Unf_dataCT[2] = TFile::Open("../../../Unfolded/26Nov2022/1D/Unfolded_HW7_flat_PY8_flat.root");    //data by MG5
  Unf_dataCT[3] = TFile::Open("../../../Unfolded/26Nov2022/1D/Unfolded_HW7_flat_PY8_flat.root");    //data by HW7

//----------------------------------------------------
//Function declaration
  TH1D* ReadHist1D(string name,TFile* root, int irbin=1);
  TH1D* ReadHist1D_v1(string name,TFile* root, int irbin=1);
  TH2D* ReadHist2D(string name,TFile* root, int irbin=1);
  void Integralhist(TH1D *hist);
  void divBybinWidth(TH1D *hist);
  void Myplotset(TH1D *Myhist,const char* XTitle, const char* YTitle);
  void Set2dHist(TH2D *MyHist, const char* XTitle, const char* YTitle,const char* ZTitle, double titoff[3], double titsize[3]);
  void SetMycanvas(TCanvas *cpt,double bs,double lm, double rm, double tm,double bm);
  void CTLegend(TLegend *legendn, const char* txt1, const char* txt2);
  TCanvas *ratio_can(int Nplot[2],float plegend[7], TH1D* data, TH1D* MC[Nplot[0]], char* lowpadx,const  char* modnam[Nplot[0]],const  char* datanm[1]);
  TCanvas *ratio_can1(int Nplot[3],float plegend[7], TH1D* data, TH1D* MC[Nplot[0]], char* lowpadx, const char* modnam[Nplot[0]], const  char* datanm[3]);
  TCanvas *ratio_canV2(int Nplot[3],float plegend[8], TH1D* data, TH1D* MC[Nplot[0]], char* lowpadx, const char* modnam[Nplot[0]], const  char* datanm[3]);
//----------------------------------------------------
  TH1D *Data_reco[ndef][njet][nkappa][nHLTmx];
  TH1D *MC_gen[nmc][ndef][njet][nkappa][nHLTmx];
  TH1D *MC_reco[nmc][ndef][njet][nkappa][nHLTmx];
  TH2D *MC_Res[nmc][ndef][njet][nkappa][nHLTmx];
  //TH1D *Psudo_Data_gen[ndef][njet][nkappa][nHLTmx];

  TH1D *hist_eff[nmc][ndef][njet][nkappa][nHLTmx];
  TH1D *hist_fake[nmc][ndef][njet][nkappa][nHLTmx];
  TH1D *hist_purity[nmc][ndef][njet][nkappa][nHLTmx];
  TH1D *hist_stbl[nmc][ndef][njet][nkappa][nHLTmx];
  
#ifdef CLOUSER
  TH1D *UnfoldCT[clos_ty][nmc][ndef][njet][nkappa][nHLTmx];
  TH1D *RefoldCT[clos_ty][nmc][ndef][njet][nkappa][nHLTmx];
  TH1D *CorrCT[clos_ty][nmc][ndef][njet][nkappa][nHLTmx];
  TH1D *ProbCT[clos_ty][nmc][ndef][njet][nkappa][nHLTmx];
  TH1D *EmatrixCT[clos_ty][nmc][ndef][njet][nkappa][nHLTmx];

  TH1D *UnfoldDataCT[clos_ty][nmc][ndef][njet][nkappa][nHLTmx];
  TH1D *RefoldDataCT[clos_ty][nmc][ndef][njet][nkappa][nHLTmx];
#endif

  TH1D *Unfold[ndef][njet][nkappa][nHLTmx];
  TH1D *Refold[ndef][njet][nkappa][nHLTmx];
  TH2D *Corr[ndef][njet][nkappa][nHLTmx];
  TH2D *Prob[ndef][njet][nkappa][nHLTmx];
  TH2D *Ematrix[ndef][njet][nkappa][nHLTmx];

#ifdef BLT
  TH1D *BLTDet[ndef][njet][nkappa][nHLTmx];
  TH1D *BLTGen[ndef][njet][nkappa][nHLTmx];
#endif

//#ifdef CLOUSER
  //TH1D *UnfoldCT[clos_ty][nmc][ndef][njet][nkappa][nHLTmx];
  //TH1D *RefoldCT[clos_ty][nmc][ndef][njet][nkappa][nHLTmx];
  //TH1D *CorrCT[clos_ty][nmc][ndef][njet][nkappa][nHLTmx];
  //TH1D *ProbCT[clos_ty][nmc][ndef][njet][nkappa][nHLTmx];
  //TH1D *EmatrixCT[clos_ty][nmc][ndef][njet][nkappa][nHLTmx];
  //TH1D *UnfoldDataCT[clos_ty][nmc][ndef][njet][nkappa][nHLTmx];
  //TH1D *RefoldDataCT[clos_ty][nmc][ndef][njet][nkappa][nHLTmx];
//#endif

//----------------------------------------------------
//Read Data
	for(int id=0; id<ndef; id++){
		for(int ij=0; ij<njet; ij++){
      			for(int ik =0 ; ik<nkappa ; ik++){
        			for(int ipt =0; ipt < nHLTmx ; ipt++){    
		Data_reco[id][ij][ik][ipt] = (TH1D*)ReadHist1D(Datadir+"/reco_jc_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_pt"+to_string(ipt)+"_eta0",Unfoldroot);

/*					
#ifdef CLOUSER
					Psudo_Data_gen[id][ij][ik][ipt] = (TH1D*)ReadHist1D(Datadir+"/gen_jc_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_pt"+to_string(ipt)+"_eta0",Unfoldroot);
#else
					Psudo_Data_gen[id][ij][ik][ipt] = (TH1D*)ReadHist1D(dirname[0]+"/gen_jc_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_pt"+to_string(ipt)+"_eta0",Unfoldroot);
#endif
 					//Exclude Underflow overlow
 		        		//TH1D *NewData;
          				//NewData = (TH1D*)Psudo_Data_gen[id][ij][ik][ipt]->Clone();
	  				//int NbinxD = NewData->GetNbinsX();
          				//NewData->Reset();
	  				//for(int ix=1; ix < NbinxD+1 ; ix++){
	  					//NewData->SetBinContent(ix,Psudo_Data_gen[id][ij][ik][ipt]->GetBinContent(ix));
	  					//NewData->SetBinError(ix, sqrt(Psudo_Data_gen[id][ij][ik][ipt]->GetBinError(ix)*Psudo_Data_gen[id][ij][ik][ipt]->GetBinError(ix)));
*/
       					}
				}
      			}
    		}
//	}
//----------------------------------------------------
//Read MC
  	for(int imc =0; imc < nmc ; imc++){
    		for(int id=0; id <ndef; id++){
      			for(int ij =0 ; ij < njet ; ij++){
				for (int ik=0; ik<nkappa; ik++){
					for(int ipt =0; ipt < nHLTmx ; ipt++){    
		MC_gen[imc][id][ij][ik][ipt] = (TH1D*)ReadHist1D(mcdir[imc]+"/gen_jc_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_pt"+to_string(ipt)+"_eta0",Unfoldroot); 
		MC_reco[imc][id][ij][ik][ipt] = (TH1D*)ReadHist1D(mcdir[imc]+"/reco_jc_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_pt"+to_string(ipt)+"_eta0",Unfoldroot);
		MC_Res[imc][id][ij][ik][ipt] = (TH2D*)ReadHist2D(mcdir[imc]+"/RM_jc_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_pt"+to_string(ipt)+"_eta0",Unfoldroot);
		hist_eff[imc][id][ij][ik][ipt] = (TH1D*)ReadHist1D(mcdir[imc]+"/miss_rate_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_pt"+to_string(ipt)+"_eta0",Unfoldroot);
		hist_fake[imc][id][ij][ik][ipt] = (TH1D*)ReadHist1D(mcdir[imc]+"/fake_rate_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_pt"+to_string(ipt)+"_eta0",Unfoldroot);
		hist_purity[imc][id][ij][ik][ipt] = (TH1D*)ReadHist1D(mcdir[imc]+"/Purity_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_pt"+to_string(ipt)+"_eta0",Unfoldroot);
		hist_stbl[imc][id][ij][ik][ipt] = (TH1D*)ReadHist1D(mcdir[imc]+"/stability_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_pt"+to_string(ipt)+"_eta0",Unfoldroot);
		  		}
      			}
    		}
	}
}//end of one MCinput root file  reading
//----------------------------------------------------
//Read Unfolded MC for CT
#ifdef CLOUSER
	cout << "Read The CT Result : "<<endl;
	for (int icl=0; icl < clos_ty; icl++){
		for (int iun=0; iun < unfold_ty; iun++){
			for(int id=0; id <ndef; id++){
	                        for(int ij=0; ij < njet ; ij++){
        	                        for(int ik=0; ik<nkappa; ik++){
                	                        for(int ipt = 0 ; ipt < nHLTmx ; ipt++){
		UnfoldCT[icl][iun][id][ij][ik][ipt] = (TH1D*)ReadHist1D(unfdir+"/Tunfold_Noreg_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_pt"+to_string(ipt)+"_eta0",Unf_root[icl]);
		RefoldCT[icl][iun][id][ij][ik][ipt] = (TH1D*)ReadHist1D(unfdir+"/Tunfold_NoReg_Refold_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_pt"+to_string(ipt)+"_eta0",Unf_root[icl]);
		CorrCT[icl][iun][id][ij][ik][ipt] = (TH1D*)ReadHist1D(unfdir+"/Tunfold_Noreg_corr_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_pt"+to_string(ipt)+"_eta0",Unf_root[icl]);
		EmatrixCT[icl][iun][id][ij][ik][ipt] = (TH1D*)ReadHist1D(unfdir+"/Tunfold_Noreg_Emat_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_pt"+to_string(ipt)+"_eta0",Unf_root[icl]);
		ProbCT[icl][iun][id][ij][ik][ipt] = (TH1D*)ReadHist1D(unfdir+"/Tunfold_Noreg_probM_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_pt"+to_string(ipt)+"_eta0",Unf_root[icl]);
						}
					}
				}
			}
		}
	} 
//----------------------------------------------------
//Unfolded data with different MC RM
	for (int icl=0; icl < clos_ty; icl++){
		for(int iun=0; iun < unfold_ty; iun++){
			for(int id=0; id <ndef; id++){
                        	for(int ij=0; ij < njet ; ij++){
                                	for(int ik=0; ik<nkappa; ik++){
                                        	for(int ipt = 0 ; ipt < nHLTmx ; ipt++){
		UnfoldDataCT[icl][iun][id][ij][ik][ipt] = (TH1D*)ReadHist1D(unfdir+"/Tunfold_Noreg_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_pt"+to_string(ipt)+"_eta0",Unf_dataCT[icl]);
		RefoldDataCT[icl][iun][id][ij][ik][ipt] = (TH1D*)ReadHist1D(unfdir+"/Tunfold_NoReg_Refold_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_pt"+to_string(ipt)+"_eta0",Unf_dataCT[icl]);
						}
					}
				}
			}
		}
	}
#endif
//----------------------------------------------------
//Read Unfolded data wih PY8
	for(int id=0; id <ndef; id++){
        	for(int ij=0; ij < njet ; ij++){
                	for(int ik=0; ik<nkappa; ik++){
                        	for(int ipt = 0 ; ipt < nHLTmx ; ipt++){
		Unfold[id][ij][ik][ipt] = (TH1D*)ReadHist1D(unfdir+"/Tunfold_Noreg_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_pt"+to_string(ipt)+"_eta0",Unfoldroot);
		Refold[id][ij][ik][ipt] = (TH1D*)ReadHist1D(unfdir+"/Tunfold_NoReg_Refold_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_pt"+to_string(ipt)+"_eta0",Unfoldroot);
		Corr[id][ij][ik][ipt] = (TH2D*)ReadHist2D(unfdir+"/Tunfold_Noreg_corr_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_pt"+to_string(ipt)+"_eta0",Unfoldroot);
		Ematrix[id][ij][ik][ipt] = (TH2D*)ReadHist2D(unfdir+"/Tunfold_Noreg_Emat_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_pt"+to_string(ipt)+"_eta0",Unfoldroot);
		Prob[id][ij][ik][ipt] = (TH2D*)ReadHist2D(unfdir+"/Tunfold_Noreg_probM_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_pt"+to_string(ipt)+"_eta0",Unfoldroot);
#ifdef BLT
		BLTDet[id][ij][ik][ipt] = (TH1D*)ReadHist1D(unfdir+"/BLTDet_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_pt"+to_string(ipt)+"_eta0",Unfoldroot);
		BLTGen[id][ij][ik][ipt] = (TH1D*)ReadHist1D(unfdir+"/BLTgen_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_pt"+to_string(ipt)+"_eta0",Unfoldroot);
#endif
				}
			}
		}
	}
cout <<"Read histogram OK" << endl;
//----------------------------------------------------
//PLot canvas declear
/*
  TCanvas *cpt5 = new TCanvas("cpt5", "canvas5", 600,575 );  //for Correlation matrix
  TCanvas *cpt6 = new TCanvas("cpt6", "canvas6", 600,575 );  //for Probabiliy matrix
  TCanvas *cpt7 = new TCanvas("cpt7", "canvas7", 600,575 );  //for Covariance matrix
  TCanvas *cpt8 = new TCanvas("cpt8", "canvas8", 600,575 );  //for Response matrix
  TCanvas *cpt9 = new TCanvas("cpt9", "canvas9", 600,575 );  //for Projection
  */
//----------------------------------------------------
//Reco comparison
  	for(int id=0; id <ndef; id++){
    		for(int ij =0 ; ij < njet ; ij++){
			for(int ik=0; ik<nkappa; ik++){
      				for(int ipt = 0 ; ipt <nHLTmx; ipt++){
			TCanvas *cpt0 = new TCanvas("cpt0", "canvas0", 600,600 );  //for Reco
			//sprintf(histname, "Data/reco_jc_d%i_j%i_k%i_pt%i_eta0", id, ij, ik, ipt); //reco_jc_d0_j0_k0_pt0_eta0
                        //Data_reco[id][ij][ik][ipt] = (TH1D*) Unfoldroot->Get(histname);
                        //cout << histname << endl;
			TH1D *MyHist  = (TH1D*) Data_reco[id][ij][ik][ipt]->Clone();
			Integralhist(MyHist);
			//divBybinWidth(MyHist);
			Myplotset(MyHist,0,0);

			MyHist->SetTitle(Form(" %s %s", htrang[ipt], Unit[0]));
                        MyHist->GetXaxis()->SetTitle("");
                        MyHist->GetYaxis()->SetTitle(Form(" %s %s_{%s}^{%s}", obs_logy[ik], obs_def[id], njets[ij], k_fact[ik]));
	
			TH1D *MC_input[nmc];
			const char *MCinput_index[nmc+1];
			const char *data_index[1];
			for(int iout = 0 ; iout < nmc ; iout++){
	  			MC_input[iout] = (TH1D*) MC_reco[iout][id][ij][ik][ipt]->Clone(); 
	  			Integralhist(MC_input[iout]);
	 			//divBybinWidth(MC_input[iout]);
	  			MCinput_index[iout] = mcnamerco[iout]; }
	  			
			data_index[0]= DataEra[iera]; 
			char lplot_xtitle[100];
			sprintf(lplot_xtitle," %s_{%s}^{%s}", obs_def[id], njets[ij], k_fact[ik]);
			//float ratio_range1[2]={1.2,0.9};
			int num1[2]={nmc,1} ;
			float lpos1[7] ={.32,0.2,0.55,0.38, .04, 1.5,0.7};
	
			cpt0->cd();
			SetMycanvas(cpt0,0,0,0,0,0);

			cpt0 =(TCanvas*)(ratio_can(num1, lpos1, MyHist, MC_input, lplot_xtitle, MCinput_index,data_index));
			SetMycanvas(cpt0,0,0.1,0.02,0.05,0.12);
			CMS_lumi( cpt0, iPeriod, iPos ); cpt0->Update();
	
			sprintf(pdfname, "RecoJC_%s.pdf(",RunEra[iera]); sprintf(pdfname1, "RecoJC_%s.pdf",RunEra[iera]); sprintf(pdfname2, "RecoJC_%s.pdf)",RunEra[iera]);
			if(id==0 && ij==0 && ik==0 && ipt ==0){cpt0->Print(pdfname,"pdf");
			}else if(id==2 && ij==1 && ik==9 && ipt==9) {cpt0->Print(pdfname2,"pdf");
			}else{  cpt0->Print(pdfname1,"pdf");};
			}
      		}
   	 }
}
//end of RECO PLOT    
/*
//----------------------------------------------------
//Reco Projection comparison for MC with Reco-fake
  	for(int id=0; id <ndef; id++){
		for (int ij=0; ij<njet; ij++){
    			for(int ik =0 ; ik < nkappa ; ik++){
      				for(int ipt = 0 ; ipt <nHLTmx; ipt++){
      			sprintf(histname, "Pythia8/%s_d%i_j%i_k%i_pt%i_eta0","ProjectX",id, ij, ik, ipt);  
      			TH1* RMx=(TH1D*) Unfoldroot->Get(histname);
      			sprintf(histname, "Pythia8/%s_d%i_j%i_k%i_pt%i_eta0","Recominusfake",id, ij, ik, ipt);
      			TH1* RecoFake=(TH1D*) Unfoldroot->Get(histname);
      			cpt9->cd();
      			SetMycanvas(cpt9,0,0.1,0.15,0.05,0.12);
      
      			TH1D *MyHist = (TH1D*)RMx->Clone();
      			Integralhist(MyHist);
        		//divBybinWidth(MyHist);
        		Myplotset(MyHist,0,0);

        		if(ipt<=8){sprintf(Title,"%s_{%s}^{%s}:     %i <P_{T}< %i %s", obs_def[id], njets[ij], k_fact[ik], HT2range[ipt] , HT2range[ipt+1] ,"GeV/c"  );}
        		else if(ipt==9){ sprintf(Title,"%s_{%s}^{%s}:      P_{T}> %i %s", obs_def[id], njets[ij], k_fact[ik], HT2range[ipt] ,"GeV/c"  );}
        		//sprintf(Yaxis," %s %s^{%s}" ,obs_logy[ik], obs_def[id], k_fact[ik]);
			sprintf(Yaxis," %s %s_{%s}^{%s}" ,obs_logy[ik], obs_def[id], njets[ij], k_fact[ik]);
        		MyHist->SetTitle(Title);
        		MyHist->GetXaxis()->SetTitle("");
        		MyHist->GetYaxis()->SetTitle(Yaxis);
        		int imc =1;
        			TH1D *MC_input[imc];
        			const char *MCinput_index[imc+1];
        			const char *data_index[1];
        			for(int iout = 0 ; iout < imc ; iout++){
          				MC_input[iout] = (TH1D*) RecoFake->Clone();
          				Integralhist(MC_input[iout]);
         				//divBybinWidth(MC_input[iout]);
          				MCinput_index[iout]= "Reco #minus Fake"; }
          				data_index[0]= "RM ProjectionX";

        			char lplot_xtitle[100];
        			//sprintf(lplot_xtitle," %s %s^{%s} %d" ,jet_num[ij],obs_def[id],k_fact[ik],HT2range[ipt]);
				sprintf(lplot_xtitle," %s_{%s}^{%s}  %s", obs_def[id], njets[ij], k_fact[ik], htrang[ipt]);
        			//float ratio_range1[2]={1.2,0.9};
        			int num1[2]={imc,1} ;
        			float lpos1[7] ={.32,0.2,0.55,0.38, .04, 1.15,0.85};

        			cpt9->Clear();
        			cpt9->cd();
        			SetMycanvas(cpt9,0,0,0,0,0);
	
        			cpt9 =(TCanvas*)(ratio_can(num1, lpos1, MyHist, MC_input, lplot_xtitle, MCinput_index,data_index));
        			CMS_lumi( cpt9, iPeriod, iPos ); cpt9->Update();

      				sprintf(pdfname, "RM_ProjectX.pdf("); sprintf(pdfname1, "RM_ProjectX.pdf"); sprintf(pdfname2, "RM_ProjectX.pdf)");
        			if(id==0 && ij==0 && ik==0 && ipt ==0){cpt9->Print(pdfname,"pdf"); 
        			}else if(id==2 && ij==1 && ik==9 && ipt==9) {cpt9->Print(pdfname2,"pdf");
        			}else{  cpt9->Print(pdfname1,"pdf");};
				}
      			}
    		}
  	}
//----------------------------------------------------
//Gen Projection comparison for MC with Gen-miss
  	for(int id=0; id <ndef; id++){
    		for(int ij =0 ; ij < njet ; ij++){
			for(int ik=0; ik<nkappa; ik++){
      				for(int ipt = 0 ; ipt <nHLTmx; ipt++){
      			sprintf(histname, "Pythia8/%s_d%i_j%i_k%i_pt%i_eta0","ProjectY",id, ij, ik, ipt); 
      			TH1* RMy=(TH1D*) Unfoldroot->Get(histname);
      			sprintf(histname, "Pythia8/%s_d%i_j%i_k%i_pt%i_eta0","Genminusmiss",id, ij, ik, ipt); 
      			TH1* GenMiss=(TH1D*) Unfoldroot->Get(histname);
      			cpt9->cd();
      			SetMycanvas(cpt9,0,0.1,0.15,0.05,0.12);

      			TH1D *MyHist  = (TH1D*)RMy->Clone();
      			Integralhist(MyHist);
        		//divBybinWidth(MyHist);
        		Myplotset(MyHist,0,0);

        		if(ipt<=8){sprintf(Title,"%s_{%s}^{%s}:     %i <P_{T}< %i %s", obs_def[id], njets[ij], k_fact[ik], HT2range[ipt] , HT2range[ipt+1] ,"GeV/c" );}
        		else if(ipt==9){sprintf(Title,"%s_{%s}^{%s}:      P_{T}> %i %s", obs_def[id], njets[ij], k_fact[ik], HT2range[ipt],"GeV/c");}
        		//sprintf(Yaxis," %s %s^{%s}" ,obs_logy[ik], obs_def[id], k_fact[ik]);
			sprintf(Yaxis," %s %s_{%s}^{%s}" ,obs_logy[ik], obs_def[id], njets[ij], k_fact[ik]);
        		MyHist->SetTitle(Title);
        		MyHist->GetXaxis()->SetTitle("");
        		MyHist->GetYaxis()->SetTitle(Yaxis);
        		int imc =1;
        			TH1D *MC_input[imc];
        			const char *MCinput_index[imc+1];
        			const char *data_index[1];
        			for(int iout = 0 ; iout < imc ; iout++){
          				MC_input[iout] = (TH1D*) GenMiss->Clone();
          				Integralhist(MC_input[iout]);
         				//divBybinWidth(MC_input[iout]);
          				MCinput_index[iout] = "Gen #minus Miss"; }
          				data_index[0]= "RM ProjectionY";

        				char lplot_xtitle[100];
        				//sprintf(lplot_xtitle," %s %s^{%s} %d" ,jet_num[ij],obs_def[id],k_fact[ik],HT2range[ipt]);
					sprintf(lplot_xtitle," %s_{%s}^{%s}  %s", obs_def[id], njets[ij], k_fact[ik], htrang[ipt]);
        				//float ratio_range1[2]={1.2,0.9};
        				int num1[2]={imc,1};
        				float lpos1[7] ={.32,0.2,0.55,0.38, .04, 1.15,0.85};

        				cpt9->Clear();
        				cpt9->cd();
        				SetMycanvas(cpt9,0,0,0,0,0);

        				cpt9 =(TCanvas*)(ratio_can(num1, lpos1, MyHist, MC_input, lplot_xtitle, MCinput_index,data_index));
        				CMS_lumi( cpt9, iPeriod, iPos ); cpt9->Update();

	        			sprintf(pdfname, "RM_ProjectY.pdf("); sprintf(pdfname1, "RM_ProjectY.pdf"); sprintf(pdfname2, "RM_ProjectY.pdf)");
        				if(id==0 && ij==0 && ik==0 && ipt==0){cpt9->Print(pdfname,"pdf");
        				}else if(id==2 && ij==1 && ik==9 && ipt==9) {cpt9->Print(pdfname2,"pdf");
        				}else{cpt9->Print(pdfname1,"pdf");};
				}
     			}
    		}
  	}
//----------------------------------------------------
//Fold comparison
  	for(int id=0; id <ndef; id++){
    		for(int ij =0 ; ij < njet ; ij++){
			for (int ik=0; ik<nkappa; ik++){
      				for(int ipt = 0 ; ipt <nHLTmx; ipt++){
      			sprintf(histname, "Folded/%s_d%i_j%i_k%i_pt%i_eta0","Fold",id, ij, ik, ipt); //check from which MC it is folded
      			TH1D* GenFold=(TH1D*)Unfoldroot->Get(histname);
        
        		TH1D *MyHist  = (TH1D*)GenFold->Clone();
        		Integralhist(MyHist);
        		//divBybinWidth(MyHist);
        		Myplotset(MyHist,0,0);

        		if(ipt<=8){sprintf(Title,"%s_{%s}^{%s}:     %i <P_{T}< %i %s", obs_def[id], njets[ij], k_fact[ik], HT2range[ipt] , HT2range[ipt+1] ,"GeV/c" );}
                        else if(ipt==9){sprintf(Title,"%s_{%s}^{%s}:      P_{T}> %i %s", obs_def[id], njets[ij], k_fact[ik], HT2range[ipt] ,"GeV/c");}
                        //sprintf(Yaxis," %s %s^{%s}" ,obs_logy[ik], obs_def[id], k_fact[ik]);
			sprintf(Yaxis," %s %s_{%s}^{%s}" ,obs_logy[ik], obs_def[id], njets[ij], k_fact[ik]);
        		MyHist->SetTitle(Title);
        		MyHist->GetXaxis()->SetTitle("");
        		MyHist->GetYaxis()->SetTitle(Yaxis);
        		int imc =1; // Only Pythia8
        			TH1D *MC_input[imc];
        			const char *MCinput_index[imc+1];
        			const char *data_index[1];
        			for(int iout = 0 ; iout < imc ; iout++){
          				MC_input[iout] = (TH1D*) MC_reco[iout][id][ij][ik][ipt]->Clone();
          				Integralhist(MC_input[iout]);
         				//divBybinWidth(MC_input[iout]);
          				MCinput_index[iout]= "Pythia8 RECO"; }
          				data_index[0]= "Folded";

        				char lplot_xtitle[100];
        				//sprintf(lplot_xtitle," %s %s^{%s} %d" ,jet_num[ij],obs_def[id],k_fact[ik],HT2range[ipt]);
					sprintf(lplot_xtitle," %s_{%s}^{%s}  %s", obs_def[id], njets[ij], k_fact[ik], htrang[ipt]);
        				//float ratio_range1[2]={1.2,0.9};
        				int num1[2]={imc,1} ;
        				float lpos1[7] ={.32,0.2,0.55,0.38, .04, 1.5,0.7};
        
					cpt9->Clear();
        				cpt9->cd();
        				SetMycanvas(cpt9,0,0,0,0,0);

        				cpt9 =(TCanvas*)(ratio_can(num1, lpos1, MyHist, MC_input, lplot_xtitle, MCinput_index,data_index));
        				CMS_lumi( cpt9, iPeriod, iPos ); cpt9->Update();

        				sprintf(pdfname, "genfold.pdf("); sprintf(pdfname1, "genfold.pdf"); sprintf(pdfname2, "genfold.pdf)");
        				if(id==0 && ij==0 && ik==0 && ipt ==0){cpt9->Print(pdfname,"pdf");
        				}else if(id==2 && ij==1 && ik==9 && ipt==9) {cpt9->Print(pdfname2,"pdf");
        				}else{cpt9->Print(pdfname1,"pdf");};
				}
      			}
    		}
  	}	 
*/
//----------------------------------------------------
//Response Matrix
	for(int imc =0; imc < nmc ; imc++){
  		for(int id=0; id <ndef; id++){
			for(int ij=0; ij<njet; ij++){
    				for(int ik =0 ; ik < nkappa ; ik++){
      					for(int ipt = 0 ; ipt <nHLTmx; ipt++){
			TCanvas *cpt0 = new TCanvas("cpt0", "canvas0", 600,600 );
			SetMycanvas(cpt0,0,0.1,0.15,0.05,0.1);
			//MC_Res[imc][ity][ivar][ipt]->RebinY(2);
        		char lplot_xtitle[800]; char lplot_ytitle[800];
			sprintf(lplot_xtitle, "RECO     %s_{%s}^{%s}", obs_def[id], njets[ij], k_fact[ik]); 
			sprintf(lplot_ytitle, "GEN      %s_{%s}^{%s}", obs_def[id], njets[ij], k_fact[ik]);
			double titoff1[3]={1.2,1.3,1.0};
        		double titsize1[3] ={0.035,0.035,0.035};
	
        		Set2dHist( MC_Res[imc][id][ij][ik][ipt],lplot_xtitle, lplot_ytitle,"",titoff1, titsize1);
			gPad->SetLogz();
			MC_Res[imc][id][ij][ik][ipt]->Draw("colz");
         
			//if(ipt<=8){sprintf(Title,"%s_{%s}^{%s}:     %i <P_{T}< %i %s", obs_def[id], jet_num[ij], k_fact[ik], HT2range[ipt] , HT2range[ipt+1] ,"GeV/c" );}
                        //else if(ipt==9){sprintf(Title,"%s_{%s}^{%s}:      P_{T}> %i %s", obs_def[id], jet_num[ij], k_fact[ik], HT2range[ipt] ,"GeV/c");}

			sprintf(Title, "%s:   ",obs_def[id]);
			TLegend *leg1 = new TLegend(0.05,0.7,0.4,0.8);     
			CTLegend(leg1,Modelnm[imc],Title); leg1->AddEntry((TObject*)0,obs_def[id] , "");leg1->SetTextColor(-8);//leg1->Draw();
        		CMS_lumi( cpt0, iPeriod, iPos ); cpt0->Update();
			sprintf(pdfname, "Response_Mat_%s_%s.pdf(",dirname[imc],RunEra[iera]); sprintf(pdfname1, "Response_Mat_%s_%s.pdf",dirname[imc],RunEra[iera]); sprintf(pdfname2, "Response_Mat_%s_%s.pdf)",dirname[imc],RunEra[iera]);
        		if(id==0 && ij==0 && ik==0 && ipt ==0){cpt0->Print(pdfname,"pdf");
        		}else if(id==2 && ij==1 && ik==9 && ipt==9) {cpt0->Print(pdfname2,"pdf");
        		}else{cpt0->Print(pdfname1,"pdf");};
        			}
			}
    		}
  	}
}
cout <<" RECO PLOT OK " << endl;
//----------------------------------------------------
//Efficiency, Purity,Fake rate, stability
  	for(int id=0; id <ndef; id++){
    		for(int ij =0 ; ij < njet ; ij++){
			for(int ik=0; ik<nkappa; ik++){
  		TCanvas *cpt1 = new TCanvas("cpt1", "canvas1", 600,700 );  //for 
  		TCanvas *cpt2 = new TCanvas("cpt2", "canvas2", 600,700 );  //for JCs
  		TCanvas *cpt3 = new TCanvas("cpt3", "canvas3", 600,700 );  //for JCs
  		TCanvas *cpt4 = new TCanvas("cpt4", "canvas4", 600,700 );  //for 

		SetMycanvas(cpt1,0,0.11,0.02,0.05,0.12);
		SetMycanvas(cpt2,0,0.11,0.03,0.05,0.12);
		SetMycanvas(cpt3,0,0.11,0.03,0.05,0.12);
		SetMycanvas(cpt4,0,0.11,0.03,0.05,0.12);

   		TLegend *leg2 = new TLegend(0.4,0.15,0.7,0.55);
		CTLegend(leg2," ","");
		leg2->SetTextFont(132);
		leg2->SetTextSize(0.035);

		TLegend *leg1 = new TLegend(0.1,0.6,0.4,0.8);
        	CTLegend(leg1,"", obs_def[id]); 
	
		for(int ipt = 0 ; ipt <nHLTmx; ipt++){
        		cpt1->cd();
			if(ipt<=8){sprintf(Title,"%i < P_{T} < %i %s", HT2range[ipt] , HT2range[ipt+1] ,"GeV/c" );}
                        else if(ipt==9){sprintf(Title," P_{T} > %i %s", HT2range[ipt],"GeV/c");}

        		for (int i = 1; i <= hist_eff[umc][id][ij][ik][ipt]->GetNbinsX(); ++i) { 
				double content = 1 - hist_eff[umc][id][ij][ik][ipt]->GetBinContent(i);
         			hist_eff[umc][id][ij][ik][ipt]->SetBinContent(i, content);
       						}

        			hist_eff[umc][id][ij][ik][ipt]->SetMinimum(-0.01); hist_eff[umc][id][ij][ik][ipt]->SetMaximum(1.1);
				hist_fake[umc][id][ij][ik][ipt]->SetMinimum(-0.01); hist_fake[umc][id][ij][ik][ipt]->SetMaximum(1.01);
				hist_purity[umc][id][ij][ik][ipt]->SetMinimum(-0.01); hist_purity[umc][id][ij][ik][ipt]->SetMaximum(1.01);
				hist_stbl[umc][id][ij][ik][ipt]->SetMinimum(-0.01); hist_stbl[umc][id][ij][ik][ipt]->SetMaximum(1.01);
   
				Myplotset(hist_eff[umc][id][ij][ik][ipt],obs_def[id],"Efficiency");
				hist_eff[umc][id][ij][ik][ipt]->SetLineColor(color[ipt]);
				hist_eff[umc][id][ij][ik][ipt]->Draw("same hist"); //leg1->Draw();
				leg2->AddEntry(hist_eff[umc][id][ij][ik][ipt], Title, "lp");
				if(ipt==9){leg2->Draw();}

				cpt2->cd();
 				Myplotset(hist_purity[umc][id][ij][ik][ipt],obs_def[id],"Purity");
 				hist_purity[umc][id][ij][ik][ipt]->SetLineColor(color[ipt]);
        			hist_purity[umc][id][ij][ik][ipt]->Draw("same hist"); //leg1->Draw();
	 			if(ipt==9){leg2->Draw();}
				cpt2->Update();

        			cpt3->cd();
        			Myplotset(hist_fake[umc][id][ij][ik][ipt],obs_def[id],"Fake rate");
        			hist_fake[umc][id][ij][ik][ipt]->SetLineColor(color[ipt]);
        			hist_fake[umc][id][ij][ik][ipt]->Draw("same hist");//leg1->Draw();
        			if(ipt==9){leg2->Draw();}
         			cpt3->Update();

				cpt4->cd();
        			Myplotset(hist_stbl[umc][id][ij][ik][ipt],obs_def[id],"Stability");
        			hist_stbl[umc][id][ij][ik][ipt]->SetLineColor(color[ipt]); 
        			hist_stbl[umc][id][ij][ik][ipt]->Draw("same hist"); //leg1->Draw();
        			if(ipt==9){leg2->Draw();}
   	 			cpt4->Update();
    				}

				CMS_lumi( cpt1, iPeriod, iPos ); cpt1->Update();
				CMS_lumi( cpt2, iPeriod, iPos ); cpt2->Update();
				CMS_lumi( cpt3, iPeriod, iPos ); cpt3->Update();
				CMS_lumi( cpt4, iPeriod, iPos ); cpt4->Update();

    				sprintf(pdfname, "efficiency_%s.pdf(",RunEra[iera]); sprintf(pdfname1, "efficiency_%s.pdf",RunEra[iera]); sprintf(pdfname2, "efficiency_%s.pdf)",RunEra[iera]);
        			if(id==0 && ij==0 && ik==0){cpt1->Print(pdfname,"pdf");
        			}else if(id==2 && ij==1 && ik==9 ) {cpt1->Print(pdfname2,"pdf");
        			}else{  cpt1->Print(pdfname1,"pdf");};

    				sprintf(pdfname, "purity_%s.pdf(",RunEra[iera]); sprintf(pdfname1, "purity_%s.pdf",RunEra[iera]); sprintf(pdfname2, "purity_%s.pdf)",RunEra[iera]); 
        			if(id==0 && ij==0 && ik==0 ){cpt2->Print(pdfname,"pdf");
        			}else if(id==2 && ij==1 && ik==9 ) {cpt2->Print(pdfname2,"pdf");
        			}else{  cpt2->Print(pdfname1,"pdf");};

     				sprintf(pdfname, "fakerate_%s.pdf(",RunEra[iera]); sprintf(pdfname1, "fakerate_%s.pdf",RunEra[iera]); sprintf(pdfname2, "fakerate_%s.pdf)",RunEra[iera]); 
        			if(id==0 && ij==0 && ik==0){cpt3->Print(pdfname,"pdf");
                                }else if(id==2 && ij==1 && ik==9 ) {cpt3->Print(pdfname2,"pdf");
                                }else{  cpt3->Print(pdfname1,"pdf");};

     				sprintf(pdfname, "stability_%s.pdf(",RunEra[iera]); sprintf(pdfname1, "stability_%s.pdf",RunEra[iera]); sprintf(pdfname2, "stability_%s.pdf)",RunEra[iera]); 
        			if(id==0 && ij==0 && ik==0){cpt4->Print(pdfname,"pdf");
        			}else if(id==2 && ij==1 && ik==9 ) {cpt4->Print(pdfname2,"pdf");
        			}else{  cpt4->Print(pdfname1,"pdf");};

				cpt1->Clear(); cpt2->Clear(); cpt3->Clear(); cpt4->Clear();
			}
     		}
  	}
cout << "Efficiency, Stability, Fakerate, Purity "<<endl;
//----------------------------------------------------
//Unfold comparisoin by Py8
  	for(int iun=0; iun < unfold_ty; iun++){
    		for(int id=0; id <ndef; id++){
			for (int ij=0; ij<njet; ij++){
      				for(int ik =0 ; ik < nkappa ; ik++){
					for(int ipt = 0 ; ipt <nHLTmx; ipt++){
		char lplot_xtitle[100]; char lplot_ytitle[100];
		TCanvas *cpt5 = new TCanvas("cpt5", "canvas5", 600,600);
		TCanvas *cpt6 = new TCanvas("cpt5", "canvas5", 600,600);
		TCanvas *cpt7 = new TCanvas("cpt5", "canvas5", 600,600);

			TCanvas *cpt0 = new TCanvas("cpt0", "canvas0", 600,700);
	  		//TH1D *MyHist  = (TH1D*)Unfold[iun][id][ij][ik][ipt]->Clone(); // unfolded histogram of variables : Tunfold_Noreg_d2_j1_k9_pt7_eta0
			TH1D *MyHist  = (TH1D*)Unfold[id][ij][ik][ipt]->Clone();
	  		Integralhist(MyHist);
	  		divBybinWidth(MyHist);
	  		Myplotset(MyHist,0,0);

	  		//if(ipt<=8){sprintf(Title,"%s_{%s}^{%s}:     %i <P_{T}< %i %s", obs_def[id], jet_num[ij], k_fact[ik], HT2range[ipt] , HT2range[ipt+1] ,"GeV/c" );}
                        //else if(ipt==9){sprintf(Title,"%s_{%s}^{%s}:      P_{T}> %i %s", obs_def[id], jet_num[ij], k_fact[ik], HT2range[ipt] ,"GeV/c");}
			
			if(ipt<=8){MyHist->SetTitle(Form(" %i < P_{T} < %i %s", HT2range[ipt], HT2range[ipt+1], Unit[0]));}
			else if(ipt==9){MyHist->SetTitle(Form(" P_{T} > %i %s", HT2range[ipt],Unit[0]));}
			MyHist->GetXaxis()->SetTitle("");
			MyHist->GetYaxis()->SetTitle(Form(" %s %s_{%s}^{%s}", obs_logy[ik], obs_def[id], njets[ij], k_fact[ik]));
			
	  		TH1D *MC_input[nmc];
	  		const char *MCinput_index[nmc], *data_index[1];
	  		for(int iout = 0 ; iout < nmc ; iout++){
	    			MC_input[iout] = (TH1D*) MC_gen[iout][id][ij][ik][ipt]->Clone(); // Gen MCs
	    			Integralhist(MC_input[iout]);
	    			divBybinWidth(MC_input[iout]);
	    			MCinput_index[iout]= mcname[iout]; }
	 
	 	 	data_index[0]= DataEra[iera];
	  		//data_index[0]= UndoldEra[1]; 
	  		//sprintf(lplot_xtitle," %s %s^{%s} %d" ,jet_num[ij],obs_def[id],k_fact[ik],HT2range[ipt]);
			//sprintf(lplot_xtitle," %s_{%s}^{%s}  %s", obs_def[id], njets[ij], k_fact[ik], htrang[ipt]);
			sprintf(lplot_xtitle," %s_{%s}^{%s}", obs_def[id], njets[ij], k_fact[ik]);
	  		//float ratio_range1[2]={1.2,0.9};
	  		int num1[2]={nmc,1} ;
	  		float lpos1[7] ={.32,0.2,0.55,0.38, .04, 1.5,0.7};
	  
	  		//cpt0->cd();
          		//SetMycanvas(cpt0,0,0,0,0,0);
	  		cpt0 =(TCanvas*)(ratio_can(num1, lpos1, MyHist, MC_input, lplot_xtitle,MCinput_index,data_index));
			SetMycanvas(cpt0,0,0.1,0.02,0.05,0.12);
	  		CMS_lumi( cpt0, iPeriod, iPos ); cpt0->Update();
	  
	  		sprintf(pdfname, "TUnfold_%s_%s.pdf(" ,RunEra[iera],dirname[iun]); sprintf(pdfname1, "TUnfold_%s_%s.pdf" ,RunEra[iera],dirname[iun]);sprintf(pdfname2, "TUnfold_%s_%s.pdf)" ,RunEra[iera],dirname[iun]);
	  		if(id==0 && ij==0 && ik==0 && ipt ==0){cpt0->Print(pdfname,"pdf");
                        }else if(id==2 && ij==1 && ik==9 && ipt==9) {cpt0->Print(pdfname2,"pdf");
                        }else{cpt0->Print(pdfname,"pdf");};
	
       			cpt5->cd();
			sprintf(lplot_xtitle, "Reco (%s_{%s}^{%s})",obs_def[id], njets[ij], k_fact[ik]);
                	sprintf(lplot_ytitle, "Gen (%s_{%s}^{%s})",obs_def[id], njets[ij], k_fact[ik]);
       			double titoff1[3]={1.2,1.3,1.0};
       			double titsize1[3] ={0.035,0.035,0.035};
       			SetMycanvas(cpt5,0,0.1,0.15,0.05,0.1);
       			gStyle->SetPaintTextFormat( "4.2f");
       			//Set2dHist(Corr[iun][id][ij][ik][ipt],lplot_xtitle, lplot_ytitle,"correlation coefficients", titoff1, titsize1);
			//Corr[iun][id][ij][ik][ipt]->Draw("colz");
			Set2dHist(Corr[id][ij][ik][ipt],lplot_xtitle, lplot_ytitle,"correlation coefficients", titoff1, titsize1);
       			Corr[id][ij][ik][ipt]->Draw("colz");
       			//Corr[iun][ity][ivar][ipt]->Draw("colz text");
			sprintf(Yaxis, " %s %s_{%s}^{%s}", obs_logy[ik], obs_def[id], njets[ij], k_fact[ik]);
			cout <<"OK1"<<endl;

       			//if(ipt<=8){sprintf(Title,"%s_{%s}^{%s}:     %i <P_{T}< %i %s", obs_def[id], jet_num[ij], k_fact[ik], HT2range[ipt] , HT2range[ipt+1] ,"GeV/c" );}
                        //else if(ipt==9){sprintf(Title,"%s_{%s}^{%s}:      P_{T}> %i %s", obs_def[id], jet_num[ij], k_fact[ik], HT2range[ipt] ,"GeV/c");}
                        //sprintf(Yaxis," %s %s^{%s}" ,obs_logy[ik], obs_def[id], k_fact[ik]);
			//sprintf(Yaxis," %s %s_{%s}^{%s}" ,obs_logy[ik], obs_def[id], njets[ij], k_fact[ik]);

        		TLegend *leg1 = new TLegend(0.05,0.6,0.4,0.8);
        		CTLegend(leg1,"Unfolded with Pythia8",Title); leg1->AddEntry((TObject*)0,obs_def[id] , ""); leg1->AddEntry((TObject*)0,Unfoldtype[iun] , "");leg1->SetTextColor(-8);leg1->Draw();
       
       			CMS_lumi( cpt5, iPeriod, iPos ); cpt5->Update();       
       			sprintf(pdfname, "TUnfold_corr_%s_%s.pdf(",RunEra[iera] ,dirname[iun]); sprintf(pdfname1, "TUnfold_corr_%s_%s.pdf" ,RunEra[iera],dirname[iun]);sprintf(pdfname2, "TUnfold_corr_%s_%s.pdf)" ,RunEra[iera],dirname[iun]);
          		if(id==0 && ij==0 && ik==0 && ipt ==0){cpt5->Print(pdfname,"pdf");
                        }else if(id==2 && ij==1 && ik==9 && ipt==9) {cpt5->Print(pdfname2,"pdf");
                        }else{cpt5->Print(pdfname,"pdf");};

       			cpt6->cd();
       			SetMycanvas(cpt6,0,0.1,0.15,0.05,0.1);
       			//Set2dHist(Prob[iun][id][ij][ik][ipt],lplot_xtitle, lplot_ytitle,"Probability",titoff1, titsize1);
       			//Prob[iun][id][ij][ik][ipt]->Draw("colz"); leg1->Draw();
			Set2dHist(Prob[id][ij][ik][ipt],lplot_xtitle, lplot_ytitle,"Probability",titoff1, titsize1);
			gPad->SetLogz();Prob[id][ij][ik][ipt]->SetMinimum(1e-6); Prob[id][ij][ik][ipt]->SetMaximum(1);// from 2D setup
			Prob[id][ij][ik][ipt]->Draw("colz"); leg1->Draw();
       			CMS_lumi( cpt6, iPeriod, iPos ); cpt6->Update();
       			sprintf(pdfname, "TUnfold_prob_%s_%s.pdf(" ,RunEra[iera],dirname[iun]); sprintf(pdfname1, "TUnfold_prob_%s_%s.pdf" ,RunEra[iera],dirname[iun]);sprintf(pdfname2, "TUnfold_prob_%s_%s.pdf)" ,RunEra[iera],dirname[iun]);
          		if(id==0 && ij==0 && ik==0 && ipt ==0){cpt6->Print(pdfname,"pdf");
          		}else if(id==2 && ij==1 && ik==9 && ipt==9) {cpt6->Print(pdfname2,"pdf");
          		}else{cpt6->Print(pdfname,"pdf");};

       			cpt7->cd();
       			SetMycanvas(cpt7,0,0.1,0.15,0.05,0.1);
       			//Set2dHist(Ematrix[iun][id][ij][ik][ipt],lplot_xtitle, lplot_ytitle,"Covariance",titoff1, titsize1);
       			//Ematrix[iun][id][ij][ik][ipt]->Draw("colz"); leg1->Draw();
			Set2dHist(Ematrix[id][ij][ik][ipt],lplot_xtitle, lplot_ytitle,"Covariance",titoff1, titsize1);
                        Ematrix[id][ij][ik][ipt]->Draw("colz"); leg1->Draw();
       			CMS_lumi( cpt7, iPeriod, iPos ); cpt7->Update();
       			sprintf(pdfname, "TUnfold_cov_%s_%s.pdf(" ,RunEra[iera],dirname[iun]); sprintf(pdfname1, "TUnfold_cov_%s_%s.pdf" ,RunEra[iera],dirname[iun]);sprintf(pdfname2, "TUnfold_cov_%s_%s.pdf)" ,RunEra[iera],dirname[iun]);
          		if(id==0 && ij==0 && ik==0 && ipt ==0){cpt7->Print(pdfname,"pdf");
          		}else if(id==2 && ij==1 && ik==9 && ipt==9) {cpt7->Print(pdfname2,"pdf");
          		}else{cpt7->Print(pdfname,"pdf");};
				}
			}
      		}
	}
}//End of Unfolded plot
//----------------------------------------------------
//Closure Test Plot with histogram
#ifdef CLOUSER
        int markersty[4]={3,26,5,6};
        for(int icl=0; icl < clos_ty ; icl++){
                for(int id=0; id<ndef; id++){
                        for(int ij=0; ij<njet; ij++){
                                for(int ik =0 ; ik<nkappa ; ik++){
                char lplot_xtitle[100];
                int iun =0;
                for(int ipt=0; ipt < nHLTmx; ipt++){
                        TCanvas *cpt0 = new TCanvas("cpt0", "canvas0", 600,600 );
                        TH1D *MyHist  = (TH1D*)UnfoldCT[icl][0][id][ij][ik][ipt]->Clone();
                        
			Integralhist(MyHist);
                        divBybinWidth(MyHist);
                        Myplotset(MyHist,0,0);
                        MyHist->GetXaxis()->SetTitle("");
                        if(ipt<=8){MyHist->SetTitle(Form(" %i < P_{T} < %i %s", HT2range[ipt] , HT2range[ipt+1],Unit[0]));}
                        else if(ipt==9){MyHist->SetTitle(Form(" P_{T} > %i %s", HT2range[ipt],Unit[0]));}
                        MyHist->GetYaxis()->SetTitle(Form("%s %s_{%s}^{%s}" ,obs_logy[ik], obs_def[id], njets[ij], k_fact[ik]));

                        TH1D *MC_input[1];
                        const char *MCinput_index[1], *data_index[3];
                        MC_input[0] = (TH1D*) MC_gen[icl][id][ij][ik][ipt]->Clone();
                        Integralhist(MC_input[0]);
                        divBybinWidth(MC_input[0]);
                        MCinput_index[0] = mcname[icl];

                        data_index[0]= Unfoldtype[icl];
                        data_index[1]= Unfoldtype[icl];
                        data_index[2]= mcname[icl];
                        sprintf(lplot_xtitle,"%s_{%s}^{%s}",obs_def[id], njets[ij], k_fact[ik]);
                        //float ratio_range1[2]={1.2,0.9};
                        int num1[3]={1,1,0} ;
                        float lpos1[8] ={.32,0.2,0.55,0.45, .055, 1.6,0.2,0.1};

                        //cpt0 =(TCanvas*)(ratio_can(num1, lpos1, MyHist, MC_input, lplot_xtitle,MCinput_index,data_index));
                        cpt0 =(TCanvas*)(ratio_canV2(num1, lpos1, MyHist, MC_input, lplot_xtitle,MCinput_index,data_index));
                        SetMycanvas(cpt0,0,0.1,0.02,0.05,0.12);
                        CMS_lumi( cpt0, iPeriod, iPos ); cpt0->Update();

                        sprintf(pdfname, "%i_ClosureHist_%s_%i.pdf(" ,icl,RunEra[iera],iun); sprintf(pdfname1, "%i_ClosureHist_%s_%i.pdf" ,icl,RunEra[iera],iun);sprintf(pdfname2, "%i_ClosureHist_%s_%i.pdf)" ,icl,RunEra[iera],iun); // check ??
                        if(id==0 && ij==0 && ik==0 && ipt==0){cpt0->Print(pdfname,"pdf");
                        }else if(id==2 && ij==1 && ik==9 && ipt==9) {cpt0->Print(pdfname2,"pdf");
                        }else{cpt0->Print(pdfname,"pdf");};
                                }
                        }
                }
        }
}//End of Unfolded plot
//----------------------------------------------------
//Closure same data with different MC
        for(int id=0; id<ndef; id++){
                for(int ij=0; ij<njet; ij++){
                        for(int ik =0 ; ik<nkappa ; ik++){
                char lplot_xtitle[100];
                int iun =0;
                for(int ipt=0; ipt < nHLTmx; ipt++){
                        TCanvas *cpt0 = new TCanvas("cpt0", "canvas0", 600,600 );  //for
                        TH1D *MyHist  = (TH1D*)UnfoldDataCT[0][0][id][ij][ik][ipt]->Clone();;
                        
			Integralhist(MyHist);
                        divBybinWidth(MyHist);
                        Myplotset(MyHist,0,0);
                        MyHist->GetXaxis()->SetTitle("");
                        if(ipt<=8){MyHist->SetTitle(Form(" %i < P_{T} < %i %s", HT2range[ipt] , HT2range[ipt+1],Unit[0]));}
                        else if(ipt==9){MyHist->SetTitle(Form(" P_{T} > %i %s", HT2range[ipt],Unit[0]));}
                        MyHist->GetYaxis()->SetTitle(Form(" %s %s_{%s}^{%s}" ,obs_logy[ik], obs_def[id], njets[ij], k_fact[ik]));

                        TH1D *MC_input[2];
                        const char *MCinput_index[2], *data_index[3];

                        for(int icl=0; icl < 3 ; icl++){
                                MC_input[icl] = (TH1D*) UnfoldDataCT[icl+1][0][id][ij][ik][ipt]->Clone();
                                Integralhist(MC_input[icl]);
                                divBybinWidth(MC_input[icl]);
                                MCinput_index[icl] = data_unf[icl+1];}

                        data_index[0]= data_unf[iun];
                        data_index[1]= "Unf by other";
                        data_index[2]= data_unf[iun];
                        sprintf(lplot_xtitle, "%s_{%s}^{%s}",obs_def[id], njets[ij], k_fact[ik]);
                        //float ratio_range1[2]={1.2,0.9};
                        int num1[3]={2,1,0} ;
                        float lpos1[8] ={.32,0.2,0.55,0.45, .05, 1.45,0.7,0.1};

                        //cpt0 =(TCanvas*)(ratio_can(num1, lpos1, MyHist, MC_input, lplot_xtitle,MCinput_index,data_index));
                        cpt0 =(TCanvas*)(ratio_canV2(num1, lpos1, MyHist, MC_input, lplot_xtitle,MCinput_index,data_index));
                        SetMycanvas(cpt0,0,0.1,0.02,0.05,0.12);
                        CMS_lumi( cpt0, iPeriod, iPos ); cpt0->Update();

                        sprintf(pdfname, "DataCT_%s_%s.pdf(" ,RunEra[iera],dirname[iun]); sprintf(pdfname1, "DataCT_%s_%s.pdf" ,RunEra[iera],dirname[iun]);sprintf(pdfname2, "DataCT_%s_%s.pdf)" ,RunEra[iera],dirname[iun]);
                        if(id==0 && ij==0 && ik==0 && ipt==0){cpt0->Print(pdfname,"pdf");
                        }else if(id==2 && ij==1 && ik==9 && ipt ==9) {cpt0->Print(pdfname2,"pdf");
                        }else{cpt0->Print(pdfname,"pdf");};
                        }
                }
        }
}//End of Unfolded plot
//----------------------------------------------------
//Refold Plot with histogram
        for(int icl=0; icl < clos_ty ; icl++){
                for(int id=0; id<ndef; id++){
                        for(int ij=0; ij<njet; ij++){
                                for(int ik =0 ; ik<nkappa ; ik++){
                char lplot_xtitle[100];
                int iun =0;
                for(int ipt=0; ipt < nHLTmx; ipt++){
                        TCanvas *cpt0 = new TCanvas("cpt0", "canvas0", 600,600 );
                        TH1D *MyHist  = (TH1D*)RefoldCT[icl][0][id][ij][ik][ipt]->Clone();;
                        
			Integralhist(MyHist);
                        divBybinWidth(MyHist);
                        Myplotset(MyHist,0,0);
                        MyHist->GetXaxis()->SetTitle("");
                        if(ipt<=8){MyHist->SetTitle(Form(" %i < P_{T} < %i %s", HT2range[ipt] , HT2range[ipt+1],Unit[0]));}
                        else if(ipt==9){MyHist->SetTitle(Form(" P_{T} > %i %s", HT2range[ipt],Unit[0]));}
                        MyHist->GetYaxis()->SetTitle(Form("%s %s_{%s}^{%s}" ,obs_logy[ik], obs_def[id], njets[ij], k_fact[ik]));

                        TH1D *MC_input[1];
                        const char *MCinput_index[1], *data_index[3];
                        MC_input[0] = (TH1D*) MC_reco[icl][id][ij][ik][ipt]->Clone();
                        Integralhist(MC_input[0]);
                        divBybinWidth(MC_input[0]);
                        MCinput_index[0] = Recomcname[icl];

                        data_index[0]= Refoldname[icl];
                        data_index[1]= Refoldname[icl];
                        data_index[2]= Recomcname[icl];
                        sprintf(lplot_xtitle,"%s_{%s}^{%s}",obs_def[id], njets[ij], k_fact[ik]);
                        //float ratio_range1[2]={1.2,0.9};
                        int num1[3]={1,1,0} ;
                        float lpos1[8] ={.32,0.2,0.55,0.45, .055, 1.6,0.2,0.1};

                        //cpt0 =(TCanvas*)(ratio_can(num1, lpos1, MyHist, MC_input, lplot_xtitle,MCinput_index,data_index));
                        cpt0 =(TCanvas*)(ratio_canV2(num1, lpos1, MyHist, MC_input, lplot_xtitle,MCinput_index,data_index));
                        SetMycanvas(cpt0,0,0.1,0.02,0.05,0.12);
                        CMS_lumi( cpt0, iPeriod, iPos ); cpt0->Update();

                        sprintf(pdfname, "%i_RefoldHist_%s_%s.pdf(" ,icl,RunEra[iera],dirname[iun]); sprintf(pdfname1, "%i_RefoldHist_%s_%s.pdf" ,icl,RunEra[iera],dirname[iun]);sprintf(pdfname2, "%i_RefoldHist_%s_%s.pdf)" ,icl,RunEra[iera],dirname[iun]);
                        if(id==0 && ij==0 && ik==0 && ipt==0){cpt0->Print(pdfname,"pdf");
                        }else if(id==2 && ij==1 && ik==9 && ipt ==9) {cpt0->Print(pdfname2,"pdf");
                        }else{cpt0->Print(pdfname,"pdf");};
                                }
                        }
                }
        }
}//End of Unfolded plot
#endif
/*
//----------------------------------------------------
//All Unfold in one plot
int tmc = 0; //0 for Py8, 1 for MG , 2 for Herwig
    	for(int id=0; id <ndef; id++){
      		for(int ij =0 ; ij < njet ; ij++){
			for(int ik=0; ik<nkappa; ik++){
        			for(int ipt = 0 ; ipt <nHLTmx; ipt++){
#ifdef CLOUSER
          		TH1D *MyHist  = (TH1D*) Psudo_Data_gen[id][ij][ik][ipt]->Clone();
#else
	  		TH1D *MyHist  = (TH1D*) MC_gen[tmc][id][ij][ik][ipt]->Clone();
#endif
    	  		Integralhist(MyHist);
         		//divBybinWidth(MyHist);
          		Myplotset(MyHist,0,0);

          		if(ipt<=8){sprintf(Title,"%s_{%s}^{%s}:     %i <P_{T}< %i %s", obs_def[id], jet_num[ij], k_fact[ik], HT2range[ipt] , HT2range[ipt+1] ,"GeV/c" );}
                        else if(ipt==9){sprintf(Title,"%s_{%s}^{%s}:      P_{T}> %i %s", obs_def[id], jet_num[ij], k_fact[ik], HT2range[ipt] ,"GeV/c");}
                        //sprintf(Yaxis," %s %s^{%s}" ,obs_logy[ik], obs_def[id], k_fact[ik]);
			sprintf(Yaxis," %s %s_{%s}^{%s}" ,obs_logy[ik], obs_def[id], njets[ij], k_fact[ik]);
          		MyHist->SetTitle(Title);
          		MyHist->GetXaxis()->SetTitle("");
          		MyHist->GetYaxis()->SetTitle(Yaxis);

          		TH1D *unfold_input[unfold_ty];
          		const char *MCinput_index[unfold_ty], *data_index[3];
          		for(int iout = 0 ; iout < unfold_ty ; iout++){
            			unfold_input[iout] = (TH1D*) Unfold[id][ij][ik][ipt]->Clone(); // unfolded histogram of variables : Tunfold_Noreg_d2_j1_k9_pt7_eta0
           			//MC_input[iout]->Rebin(2);
            			Integralhist(unfold_input[iout]);
           			//divBybinWidth(unfold_input[iout]);
         			//MCinput_index[iout]= Unfoldtype[iout]; }
            			//MCinput_index[iout]= closuretype[iout]; }
#ifdef CLOUSER
           		MCinput_index[iout]= closuretype[iout];
#else
           		MCinput_index[iout]= Unfoldtype[iout];
#endif
      			}

      			for(int iout = 0 ; iout < unfold_ty ; iout++){
            			cout << unfold_input[iout]->GetTitle()<< endl;
            			cout <<  MCinput_index[iout] << endl;
          		}

	        data_index[0]= mcnamegen[0];
        	data_index[1]= "MC";
         	data_index[2]= "Unfolded";
          	char lplot_xtitle[100];
          	//sprintf(lplot_xtitle," %s %s^{%s} %d" ,jet_num[ij],obs_def[id],k_fact[ik],HT2range[ipt]);
		sprintf(lplot_xtitle," %s_{%s}^{%s}  %s", obs_def[id], njets[ij], k_fact[ik], htrang[ipt]);
          	//float ratio_range1[2]={1.2,0.9};
          	int num1[3]={unfold_ty,1,0};
          	float lpos1[7] ={.32,0.2,0.55,0.38, .04, 1.3,0.75};
          	//if(id==1 && ivar==2 && ipt ==0){lpos1[0] =.45; lpos1[1]=0.65; lpos1[2]=0.8; lpos1[3]=0.9;}

          	cpt0->cd();
          	cpt0->SetBorderSize(0);
          	cpt0->SetRightMargin(0.0);
          	cpt0->SetTopMargin(0.0);
          	cpt0 =(TCanvas*)(ratio_can1(num1, lpos1, MyHist, unfold_input, lplot_xtitle,MCinput_index,data_index));
          	CMS_lumi( cpt0, iPeriod, iPos ); cpt0->Update();

          	sprintf(pdfname, "unfold_%i.pdf(", 1234); sprintf(pdfname1, "unfold_%i.pdf" ,1234);sprintf(pdfname2, "unfold_%i.pdf)",1234);
          	if(id==0 && ij==0 && ik==0 && ipt ==0){cpt0->Print(pdfname,"pdf");
          	}else if(id==2 && ij==1 && ik==9 && ipt==9) {cpt0->Print(pdfname2,"pdf");
          	}else{cpt0->Print(pdfname,"pdf");};
      			}
    		}
	}
}  
//----------------------------------------------------
//Refold comparisoin
  	for(int iun=0; iun < unfold_ty; iun++){
    		for(int id=0; id <ndef; id++){
			for(int ij=0; ij<njet; ij++){
      				for(int ik =0 ; ik < nkappa ; ik++){
					for(int ipt = 0 ; ipt <nHLTmx; ipt++){ 
	  		TH1D *MyHist  = (TH1D*) Refold[id][ij][ik][ipt]->Clone();
	  		Integralhist(MyHist);
	  		//divBybinWidth(MyHist);
	  		Myplotset(MyHist,0,0);
	  
	  		if(ipt<=8){sprintf(Title,"%s_{%s}^{%s}:     %i <P_{T}< %i %s", obs_def[id], jet_num[ij], k_fact[ik], HT2range[ipt] , HT2range[ipt+1] ,"GeV/c" );}
                        else if(ipt==9){sprintf(Title,"%s_{%s}^{%s}:      P_{T}> %i %s", obs_def[id], jet_num[ij], k_fact[ik], HT2range[ipt] ,"GeV/c");}
                        //sprintf(Yaxis," %s %s^{%s}" ,obs_logy[ik], obs_def[id], k_fact[ik]);
			sprintf(Yaxis," %s %s_{%s}^{%s}" ,obs_logy[ik], obs_def[id], njets[ij], k_fact[ik]);
	  		MyHist->SetTitle(Title);
	  		MyHist->GetXaxis()->SetTitle("");
	  		MyHist->GetYaxis()->SetTitle(Yaxis);
	 
	  		TH1D *MC_input[nmc];
	  		const char *MCinput_index[nmc], *data_index[1];
	  		for(int iout = 0 ; iout < nmc ; iout++){
	    			MC_input[iout] = (TH1D*) MC_reco[iout][id][ij][ik][ipt]->Clone(); // All MC reco histograms
	   			//MC_input[iout]->Rebin(2);
	    			Integralhist(MC_input[iout]);
	   			//divBybinWidth(MC_input[iout]);
	    			MCinput_index[iout]= mcnamerco[iout]; }
	  
	  			data_index[0]= RefoldEra[iun];
	  			char lplot_xtitle[100];
	  			//sprintf(lplot_xtitle," %s %s^{%s} %d" ,jet_num[ij],obs_def[id],k_fact[ik],HT2range[ipt]);
				sprintf(lplot_xtitle," %s_{%s}^{%s}  %s", obs_def[id], njets[ij], k_fact[ik], htrang[ipt]);
	  			//float ratio_range1[2]={1.2,0.9};
	  			int num1[2]={nmc,1};
	  			float lpos1[7] ={.32,0.2,0.55,0.38, .04, 1.5,0.7};
	  
	  			cpt0->cd();
	  			cpt0->SetBorderSize(0);
	  			cpt0->SetRightMargin(0.0);
	  			cpt0->SetTopMargin(0.0);
	  			cpt0 =(TCanvas*)(ratio_can(num1, lpos1, MyHist, MC_input, lplot_xtitle,MCinput_index,data_index));
	  			CMS_lumi( cpt0, iPeriod, iPos ); cpt0->Update();
	  
	  			sprintf(pdfname, "refold_%i.pdf(" ,iun); sprintf(pdfname1, "refold_%i.pdf" ,iun);sprintf(pdfname2, "refold_%i.pdf)" ,iun); 
	  			if(id==0 && ij==0 && ik==0 && ipt ==0){cpt0->Print(pdfname,"pdf");
	  			}else if(id==2 && ij==1 && ik==9 && ipt==9) {cpt0->Print(pdfname2,"pdf");
	  			}else{cpt0->Print(pdfname,"pdf");};
	  			}
			}
      		}
    	}
}//End of Refolded plot  
//----------------------------------------------------
//Refold in one plot
    	for(int id=0; id <ndef; id++){
		for(int ij=0; ij<njet; ij++ ){
    	  		for(int ik =0 ; ik < nkappa ; ik++){
        			for(int ipt = 0 ; ipt <nHLTmx; ipt++){
          		TH1D *MyHist  = (TH1D*) MC_reco[tmc][id][ij][ik][ipt]->Clone(); // MC Reco histograms
          		Integralhist(MyHist);
         		//divBybinWidth(MyHist);
          		Myplotset(MyHist,0,0);

          		if(ipt<=8){sprintf(Title,"%s_{%s}^{%s}:     %i <P_{T}< %i %s", obs_def[id], jet_num[ij], k_fact[ik], HT2range[ipt] , HT2range[ipt+1] ,"GeV/c" );}
                        else if(ipt==9){sprintf(Title,"%s_{%s}^{%s}:      P_{T}> %i %s", obs_def[id], jet_num[ij], k_fact[ik], HT2range[ipt] ,"GeV/c");}
                        //sprintf(Yaxis," %s %s^{%s}" ,obs_logy[ik], obs_def[id], k_fact[ik]);
			sprintf(Yaxis," %s %s_{%s}^{%s}" ,obs_logy[ik], obs_def[id], njets[ij], k_fact[ik]);
          		MyHist->SetTitle(Title);
          		MyHist->GetXaxis()->SetTitle("");
          		MyHist->GetYaxis()->SetTitle(Yaxis);

          		TH1D *Refold_input[unfold_ty];
          		const char *MCinput_index[unfold_ty], *data_index[3];
          		for(int iout = 0 ; iout < unfold_ty ; iout++){
            			Refold_input[iout] = (TH1D*) Refold[id][ij][ik][ipt]->Clone();
           			//MC_input[iout]->Rebin(2);
            			Integralhist(Refold_input[iout]);
           			//divBybinWidth(Refold_input[iout]);
            			MCinput_index[iout]= RefoldEra[iout]; }

          			data_index[0]= mcnamerco[tmc];
          			data_index[1]= "MC";
          			data_index[2]= "Refolded";
          			char lplot_xtitle[100];
          			//sprintf(lplot_xtitle," %s %s^{%s} %d" ,jet_num[ij],obs_def[id],k_fact[ik],HT2range[ipt]);
				sprintf(lplot_xtitle," %s_{%s}^{%s}  %s", obs_def[id], njets[ij], k_fact[ik], htrang[ipt]);
          			//float ratio_range1[2]={1.2,0.9};
          			int num1[3]={unfold_ty,1,0};
          			float lpos1[7] ={.32,0.2,0.55,0.38, .04, 1.5,0.7};

          			cpt0->cd();
          			cpt0->SetBorderSize(0);
          			cpt0->SetRightMargin(0.0);
          			cpt0->SetTopMargin(0.0);
          			cpt0 =(TCanvas*)(ratio_can1(num1, lpos1, MyHist, Refold_input, lplot_xtitle,MCinput_index,data_index));
          			CMS_lumi( cpt0, iPeriod, iPos ); cpt0->Update();

          			sprintf(pdfname, "refold_%i.pdf(", 1234); sprintf(pdfname1, "refold_%i.pdf" ,1234);sprintf(pdfname2, "refold_%i.pdf)",1234);
          			if(id==0 && ij==0 && ik==0 && ipt ==0){cpt0->Print(pdfname,"pdf");
          			}else if(id==2 && ij==1 && ik==9 && ipt==9) {cpt0->Print(pdfname2,"pdf");
          			}else{cpt0->Print(pdfname,"pdf");};
			}
      		}
    	}
} 
*/
}//end of main program
//-------------------------------------------
//-----------------Functions-----------------
//-------------------------------------------
//Ratio plot function
TCanvas *ratio_can(int Nplot[2],float plegend[7], TH1D* data, TH1D* MC[Nplot[0]], char* lowpadx, const char* modnam[Nplot[0]], const  char* datanm[1]){
  //Nplot[0] = number of MC enetered
  //Nplot[1] = place 1 if upper part is log scale needed
  //plegend[0->3] = x1,y1,x2,y2 of the legend of the upper plot
  //plegend[4]= text size
  //plegend[5-6]= ratio plot axis range
  //data = data histogram
  // MC = monte carlo histogram array
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
  //ymax = data->GetMaximum(); //Divide the histogram with bin width
  
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
  //TLegend *HT_range = new TLegend(.15,.06,.55,.1);
  TLegend *HT_range = new TLegend(plegend[0]-0.08,.06,plegend[2],.1);
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
  
  //if(Nplot[1]==1){gPad->SetLogy();} //condition for log scale
  gStyle->SetOptStat(0);
  gPad->SetTickx();
  gPad->SetTicky();
  
  double chi2rat[Nplot[0]];        //for chi2 plot in legend 
  double chi2Ndfrat[Nplot[0]];     //for chi2 plot in legend
  
  //int color[40] = {2,4,6,7,8,49,1,41,42,30,46,28,29,38,30,12,37,49,9,32,2,3,6,7,8,49,1,41,42,30,46,28,29,38,30,12,37,49,9,32};
  //int style[40]={1,2,3,4,5,6,7,8,9,9,2,3,4,5,6,7,8,9,9,9,1,2,3,4,5,6,7,8,9,9,2,3,4,5,6,7,8,9,9,9};
  int color[40] = {2,62,30,46,8,49,1,41,42,30,46,28,29,38,30,12,37,49,9,32,2,3,6,7,8,49,1,41,42,30,46,28,29,38,30,12,37,49,9,32};
        int style[40]={1,1,1,1,5,6,7,8,9,9,2,3,4,5,6,7,8,9,9,9,1,2,3,4,5,6,7,8,9,9,2,3,4,5,6,7,8,9,9,9};
  for(int iup =0; iup < Nplot[0] ; iup++){
  	MC[iup]->SetLineStyle(style[iup]);
    	MC[iup]->SetLineColor(color[iup]);
    	MC[iup]->SetLineWidth(2);
    	//MC[iup]->Draw("same hist ");// gPad->SetLogy();        
  //Addition fo chi2/Ndf  with legend
  int nn = data->GetNbinsX();
  Double_t resl[nn] , chi2l;
  Int_t ndfl ,igoodl;
  TH1D *chidatal = (TH1D*)data->Clone("chidatal");   //for chi square test   
  TH1D *chiMCl = (TH1D*)MC[iup]->Clone("chiMCl");    //for chi square test
  chiMCl->Chi2TestX(chidatal, chi2l, ndfl, igoodl, "WU", resl);
    
  //chidata->Chi2TestX(chiMC, chi2, ndf, igood, "WU", res);
  //cout << "Ndf value =" << ndfl << endl;
  //cout << "chi2 value=" << chi2l << endl;
  chi2rat[iup]=chi2l;
  chi2Ndfrat[iup]=chi2l/ndfl;
  //end of chi2/Ndf
  //Divide the bin contents/error by bin width
  /*  
  int tmpnbn = MC[iup]->GetNbinsX();
  for (int ix=0; ix<tmpnbn; ix++) {
  	double awidth = MC[iup]->GetBinWidth(ix+1); // tmpwid;
      	MC[iup]->SetBinContent(ix+1, MC[iup]->GetBinContent(ix+1)/awidth);
      	double error = MC[iup]->GetBinError(ix+1);
      	MC[iup]->SetBinError(ix+1, error/awidth);
      }
  */
  //end Divide the bin contents/error by bin width    
  MC[iup]->Draw("same hist e1 ");
  }//end of Montecarlo loop
  /*
  int tmpnbd = data->GetNbinsX();
  for (int ix=0; ix<tmpnbd; ix++) {
  	double awidth = data->GetBinWidth(ix+1); // /tmpwid;
    	data->SetBinContent(ix+1, data->GetBinContent(ix+1)/awidth);
    	double error = data->GetBinError(ix+1);
    	data->SetBinError(ix+1, error/awidth);
   	}
   */
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
    				else{MC_inputerr->SetBinError(ix+1,rel_err[ix+1]);
    				}
    			} //loop over bin
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
  sprintf(MCindex,"%s" , datanm[0]);              //use if chi2 is not needed in the legend
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
  padfun2->SetGridy();   //Horizontal grid
  padfun2->Draw();
  padfun2->cd();
  gPad->SetTickx();
  gPad->SetTicky();
  //gStyle->SetErrorX(0.5);
  for(int ilow =0 ; ilow <  Nplot[0] ; ilow++){  //loop for ratio plot
  	TH1D *rh2;
    	//if(data->GetBinContent(1) > 0 ||  data->GetBinContent(10) > 0) {
    	//if(  data->GetBinContent(10) > 0) {
    	rh2 = (TH1D*)MC[ilow]->Clone("rh2"); 
    	//rh2->Sumw2();
    	rh2->Divide(data);     //MC divide by data
    /* 	}else{
       	rh2 = (TH1D*)data->Clone("rh2");
       	rh2->Sumw2();
       	rh2->SetLineColor(MC[ilow]->GetLineColor());
       	rh2->SetLineWidth(MC[ilow]->GetLineWidth());
       	rh2->SetLineStyle(MC[ilow]->GetLineStyle());
       	rh2->Divide(MC[ilow]);
       }
       */
    	rh2->SetTitle("");
    	//rh2->SetLineColor(kBlack);
    	//rh2->SetMaximum(lowYrange[0]);       //.. range
    	//rh2->SetMinimum(lowYrange[1]);       //Define Y ..
    	//cout << "max = " << plegend[5] << "  min = " << plegend[6] << endl;   
    	rh2->SetMinimum(plegend[6]);           //Define Y ..
    	rh2->SetMaximum(plegend[5]);           //.. range
    	rh2->SetStats(0);                      //No statistics on lower plot
    	//rh2->Divide(data);                   //MC divide by data
    	//rh2->SetMarkerStyle(21);
    	//rh2->Draw(" same e1");
    	rh2->Draw("same e1");
    	//MC_inputerr->Draw("same E2");        //Draw the Uncertainty in Ratio plot
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
    	//rh2->GetYaxis()->SetLabelFont(1.0);   //Absolute font size in pixel (precision 3)
    	rh2->GetYaxis()->SetLabelSize(0.09);
  	}
  	canvas->Update();
  return canvas;
}//end of ratio plot function
//----------------------------------------------------
TCanvas *ratio_canV2(int Nplot[3],float plegend[8], TH1D* data, TH1D* MC[Nplot[0]], char* lowpadx, const char* modnam[Nplot[0]], const  char* datanm[3]){
        //Nplot[0] = number of MC enetered
        //Nplot[1] = place 1 if upper part is log scale needed
        //Nplot[2] = place 1 if MC/Data or 0 if data/MC
        //plegend[0->3] = x1,y1,x2,y2 of the legend of the upper plot
        //plegend[4]= text size
        //plegend[5-6]= ratio plot axis range
        //plegend[7]= ratio plot Yaxis Title size
        //data = data histogram
        //MC = monte carlo histogram array
        //legendN = name of the legends for mC one by one
        //lowpadx = x axis title of the ratio plot
        //datanm[0] = name of data (or MC in case one MC is used)
        //datanm[1] = name of MC/data Data term
        //datanm[2] = name of Mc/Data MC term
        bool isstat =0;
        TCanvas *canvas =new TCanvas("cptfun", "canvas_fun", 680,750 );
        canvas->cd();
        canvas->SetRightMargin(0.02);
        canvas->SetTopMargin(0.04);

        char ratioXaxis1[100]; char ratioYaxis[100]; char MCindex[100];
        //float ymax;
        //canvas->SetBottomMargin(0.1);
        data->GetYaxis()->SetLabelSize(0.055);
        data->GetXaxis()->SetLabelSize(0.0);
        data->GetYaxis()->SetTitleSize(0.065);
        data->GetYaxis()->SetTitleOffset(1.33);
        data->GetXaxis()->SetTitleSize(0.053);
        data->GetXaxis()->SetTitleOffset(0.12);
        data->GetYaxis()->CenterTitle();
        data->GetXaxis()->CenterTitle();
        data->GetXaxis()->SetTickLength(0.05);
        data->GetYaxis()->SetTickLength(0.05);

 	data->SetLineWidth(2);
        data->SetMarkerStyle(9);
        data->SetMarkerSize(.9);
        int ifont =42;
        data->GetXaxis()->SetTitleFont(ifont);
        data->GetYaxis()->SetTitleFont(ifont);
        data->SetTitleFont(ifont);
        //ymax = data->GetMaximum();
        //Divide the histogram with bin width

        TPad *padfun1 = new TPad("padfun1", "padfun1", 0.05, 0.35, 1.0, 1.0);
        padfun1->SetBottomMargin(0.0); // Upper and lower plot are joined
        //padfun1->SetTopMargin(0.05); // Upper and lowd
        padfun1->SetLeftMargin(0.17); // Left
        padfun1->SetRightMargin(.02);
        padfun1->Draw();             // Draw the upper pad: pad1
        padfun1->cd();
        //data->SetFillColor(kYellow);
        //data->SetFillStyle(1111);

        //create the legend of ht2 range
        char MC_HTindex[100];
        //TLegend *HT_range = new TLegend(.15,.06,.55,.1);
        TLegend *HT_range = new TLegend(plegend[0]-0.08,.06,plegend[2],.1);
        HT_range->SetFillStyle(0);
        HT_range->SetBorderSize(0);
        HT_range->SetTextSize(0.06);
        sprintf(MC_HTindex,"%s" , data->GetTitle());   //legend for the HT value range
        HT_range->AddEntry((TObject*)0, MC_HTindex, "" );
        HT_range->Draw();
        data->SetTitle("");
        //end the legend of ht2 range
        data->Draw("e2");

	//if(Nplot[1]==1){gPad->SetLogy();} //condition for log scale
        gStyle->SetOptStat(0);
        gPad->SetTickx();
        gPad->SetTicky();

        double chi2rat[Nplot[0]];    // for chi2 plot in legend
        double chi2Ndfrat[Nplot[0]]; // for chi2 plot in legend

        int color[40] = {2,62,30,46,8,49,1,41,42,30,46,28,29,38,30,12,37,49,9,32,2,3,6,7,8,49,1,41,42,30,46,28,29,38,30,12,37,49,9,32};
        int style[40]={1,2,3,4,5,6,7,8,9,9,2,3,4,5,6,7,8,9,9,9,1,2,3,4,5,6,7,8,9,9,2,3,4,5,6,7,8,9,9,9};
        for(int iup =0; iup < Nplot[0] ; iup++){
                MC[iup]->SetLineStyle(style[iup]);
                MC[iup]->SetLineColor(color[iup]);
                MC[iup]->SetLineWidth(3);
                //MC[iup]->Draw("same hist ");// gPad->SetLogy();

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

        MC[iup]->Draw("same hist e1 ");
        //MC[iup]->Draw("same hist");
        }//end of Montecarlo loop

        if(isstat){data->Draw(" same e1");}
        else{data->Draw("same ");}

	//Stat Error(DATA)
        TH1D *datastaterr = (TH1D*)data->Clone("staterr"); datastaterr->Reset();
        for(int ierr=1; ierr < datastaterr->GetNbinsX()+1; ierr++){
                double bincont = data->GetBinContent(ierr);
                double serr = data->GetBinError(ierr);
                double relerr =sqrt(pow((serr/bincont),2.));
                datastaterr->SetBinContent(ierr,1.0);
                datastaterr->SetBinError(ierr,relerr);
                }
        gStyle->SetErrorX(0.5);
        gStyle->IsTransparent();
        datastaterr->SetFillColor(5);
        datastaterr->SetMarkerSize(0.1);
        datastaterr->SetFillStyle(3002);
        //datastaterr->SetFillStyle(1111);
        //datastaterr->SetFillColor(21);
        //datastaterr->SetMarkerSize(0.1);

        TLegend *legendn = new TLegend(plegend[0],plegend[1],plegend[2],plegend[3]);
        legendn->SetFillStyle(0);
        legendn->SetBorderSize(0);
        legendn->SetTextSize(plegend[4]);
        legendn->SetTextFont(42);
        sprintf(MCindex,"%s" , datanm[0]);      // use if chi2 is not needed in the legend
        legendn->AddEntry(data, MCindex,"lp");
        for(int iup =0 ; iup <  Nplot[0] ; iup++){
                sprintf(MCindex,"%s" , modnam[iup]);      // use if chi2 is not needed in the legend
                //sprintf(MCindex,"%s-#chi^{2}/NDF: %.2f" ,legendN[iup],chi2Ndfrat[iup]);   //legend with chi2/Ndf value
                legendn->AddEntry(MC[iup], MCindex ,"lp");
                }
        if(isstat)legendn->AddEntry(datastaterr, "Statistical Uncertainty","f"); //Stat unc
        legendn->Draw();
        HT_range->Draw();

	//ratio plot pad
        canvas->cd();          // Go back to the main canvas before defining pad2
        TPad *padfun2 = new TPad("padfun2", "padfun2", 0.05, 0.02,1.0, 0.35);
        padfun2->SetTopMargin(0);
        padfun2->SetBottomMargin(.4);
        padfun2->SetLeftMargin(.17);
        padfun2->SetRightMargin(.02);
        padfun2->SetGridy(); // Horizontal grid
        padfun2->Draw();
        padfun2->cd();
        gPad->SetTickx();
        gPad->SetTicky();
        //gStyle->SetErrorX(0.5);

        for(int ilow =0 ; ilow <  Nplot[0] ; ilow++){  //loop over MC for ratio plot
                TH1D *rh2;
                //if(data->GetBinContent(1) > 0 ||  data->GetBinContent(10) > 0) {
                //if(data->GetBinContent(10) > 0) {
                rh2 = (TH1D*)MC[ilow]->Clone("rh2");
                //rh2->Sumw2();
                rh2->Divide(data);     //MC divide by data
    /*          }else{
                        rh2 = (TH1D*)data->Clone("rh2");
                        rh2->Sumw2();
                        rh2->SetLineColor(MC[ilow]->GetLineColor());
                        rh2->SetLineWidth(MC[ilow]->GetLineWidth());
                        rh2->SetLineStyle(MC[ilow]->GetLineStyle());
                        rh2->Divide(MC[ilow]);
       }*/
                rh2->SetTitle("");
                //rh2->SetLineColor(kBlack);
                rh2->SetMinimum(plegend[6]);  // Define Y ..
                rh2->SetMaximum(plegend[5]);  // .. range
                rh2->SetStats(0);      // No statistics on lower plot
                //gStyle->SetStatColor(2);      // No statistics on lower plot
                //gStyle->SetStatX(0.7);;      // No statistics on lower plot
                //rh2->SetMarkerStyle(21);
                rh2->SetMarkerStyle(kFullCircle); //For Closure Test
                rh2->Draw("same p e1");
                //rh2->Draw("Same P text ");
                //rh2->Draw("Same P ");

		//Stat Error Monte Carlo
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
                if(isstat) datastaterr->Draw("E2 Same");
                if(isstat) mstaterr->Draw("E2 Same");
                //MC_inputerr->Draw("same E2");//Draw the Uncertainty in Ratio plot
                rh2->GetXaxis()->SetTitleSize(0.12);
                //rh2->GetXaxis()->SetTitleFont(43);
                rh2->GetXaxis()->SetTitleOffset(1.33);
                rh2->GetXaxis()->SetLabelSize(0.13);
                rh2->GetXaxis()->CenterTitle();
                sprintf(ratioXaxis1," %s" ,lowpadx);
                rh2->GetXaxis()->SetTitle(ratioXaxis1);
                //rh2->GetYaxis()->SetTitle("#frac{Pythia8-unf}{Pythia8-gen}");
                //if(Nplot[2]==1){ sprintf(ratioYaxis,"%s/%s" ,datanm[1],datanm[2]);}
                //if(Nplot[2]==0){ sprintf(ratioYaxis,"%s/%s" ,datanm[2],datanm[1]);}
                if(Nplot[2]==1){ sprintf(ratioYaxis,"#frac{%s}{%s}" ,datanm[1],datanm[2]);}
                if(Nplot[2]==0){ sprintf(ratioYaxis,"#frac{%s}{%s}" ,datanm[2],datanm[1]);}

		rh2->GetYaxis()->SetTitle(ratioYaxis);
                rh2->GetYaxis()->CenterTitle();
                rh2->GetYaxis()->SetNdivisions(505);
                rh2->GetYaxis()->SetTitleSize(plegend[7]);
                rh2->GetXaxis()->SetTitleFont(ifont);
                rh2->GetYaxis()->SetTitleFont(ifont);
                rh2->GetYaxis()->SetTitleOffset(0.6);
                rh2->GetXaxis()->SetTickLength(0.06);
                rh2->GetYaxis()->SetTickLength(0.06);
                //rh2->GetYaxis()->SetLabelFont(1.0); // Absolute font size in pixel (precision 3)
                rh2->GetYaxis()->SetLabelSize(0.09);
        }
        canvas->Update();
return canvas;
}//end of ratio plot function
//----------------------------------------------------
void Integralhist(TH1D *hist){ hist->Scale(1/(hist->Integral()));}
//----------------------------------------------------
void divBybinWidth(TH1D *hist){
	int tmpnbn = hist->GetNbinsX();
  	for (int ix=0; ix<tmpnbn; ix++) {
    	double awidth = hist->GetBinWidth(ix+1); //tmpwid;
    		hist->SetBinContent(ix+1, hist->GetBinContent(ix+1)/awidth);
    		double error = hist->GetBinError(ix+1);
    			hist->SetBinError(ix+1, error/awidth);
  			}
		}
//----------------------------------------------------
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
//----------------------------------------------------
void SetMycanvas(TCanvas *cpt,double bs,double lm, double rm, double tm,double bm){
	cpt->SetBorderSize(bs);
        cpt->SetLeftMargin(lm);
        cpt->SetRightMargin(rm);
        cpt->SetTopMargin(tm); 
       	cpt->SetBottomMargin(bm); 
}
//----------------------------------------------------
void CTLegend(TLegend *legendn, const char* txt1, const char* txt2){ 
  	legendn->SetFillStyle(0);
  	legendn->SetBorderSize(0);
  	legendn->SetTextSize(0.03);
  	legendn->SetTextFont(42);
  	legendn->AddEntry((TObject*)0, txt1, "");
  	legendn->AddEntry((TObject*)0, txt2, "");  
}
//----------------------------------------------------
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
//----------------------------------------------------
TCanvas *ratio_can1(int Nplot[3],float plegend[7], TH1D* data, TH1D* MC[Nplot[0]], char* lowpadx, const char* modnam[Nplot[0]], const  char* datanm[3]){
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
  //ymax = data->GetMaximum(); //Divide the histogram with bin width

  TPad *padfun1 = new TPad("padfun1", "padfun1", 0, 0.30, 1.0, 1.0);
  padfun1->SetBottomMargin(0.01); //Upper and lower plot are joined
  padfun1->SetTopMargin(0.05);    //Upper and lowd
  padfun1->SetRightMargin(.04);
  padfun1->Draw("");                //Draw the upper pad: pad1
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
  //data->Draw("e2");
  data->Draw(" ");

  //if(Nplot[1]==1){gPad->SetLogy();}  //condition for log scale
  gStyle->SetOptStat(0);
  gPad->SetTickx();
  gPad->SetTicky();

  double chi2rat[Nplot[0]];          //for chi2 plot in legend
  double chi2Ndfrat[Nplot[0]];       //for chi2 plot in legend
  int color[40] = {2,4,6,7,8,49,1,41,42,30,46,28,29,38,30,12,37,49,9,32,2,3,6,7,8,49,1,41,42,30,46,28,29,38,30,12,37,49,9,32};
  int style[40]={1,2,3,4,5,6,7,8,9,9,2,3,4,5,6,7,8,9,9,9,1,2,3,4,5,6,7,8,9,9,2,3,4,5,6,7,8,9,9,9};
  for(int iup =0; iup < Nplot[0] ; iup++){
  	MC[iup]->SetLineStyle(style[iup]);
    	MC[iup]->SetLineColor(color[iup]);
    	MC[iup]->SetLineWidth(2);
    	//MC[iup]->Draw("same hist ");// gPad->SetLogy();
	//Addition fo chi2/Ndf  with legend
    	int nn = data->GetNbinsX();
    	Double_t resl[nn] , chi2l;
    	Int_t ndfl ,igoodl;
    	TH1D *chidatal = (TH1D*)data->Clone("chidatal");   //for chi square test

    	TH1D *chiMCl = (TH1D*)MC[iup]->Clone("chiMCl");    //for chi square test
    	chiMCl->Chi2TestX(chidatal, chi2l, ndfl, igoodl, "WU", resl);

    	//chidata->Chi2TestX(chiMC, chi2, ndf, igood, "WU", res);
    	//cout << "Ndf value =" << ndfl << endl;
    	//cout << "chi2 value=" << chi2l << endl;
    	chi2rat[iup]=chi2l;
    	chi2Ndfrat[iup]=chi2l/ndfl;
	//end of chi2/Ndf
	//Divide the bin contents/error by bin width
	/*
    	int tmpnbn = MC[iup]->GetNbinsX();
      	for (int ix=0; ix<tmpnbn; ix++) {
      	double awidth = MC[iup]->GetBinWidth(ix+1); // tmpwid;
      		MC[iup]->SetBinContent(ix+1, MC[iup]->GetBinContent(ix+1)/awidth);
      		double error = MC[iup]->GetBinError(ix+1);
      			MC[iup]->SetBinError(ix+1, error/awidth);
      			}
      	*/
	//end Divide the bin contents/error by bin width

    	MC[iup]->Draw("same hist e1 ");
	//MC[iup]->Draw("same hist ");
  	}//end of Montecarlo loop

  /*
  int tmpnbd = data->GetNbinsX();
  for (int ix=0; ix<tmpnbd; ix++) {
  	double awidth = data->GetBinWidth(ix+1); //tmpwid;
    		data->SetBinContent(ix+1, data->GetBinContent(ix+1)/awidth);
    	double error = data->GetBinError(ix+1);
    		data->SetBinError(ix+1, error/awidth);
    		}
  */
  data->Draw(" same e1");
  //data->Draw(" same ");
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
    		else{MC_inputerr->SetBinError(ix+1,rel_err[ix+1]);
    		}
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
  	sprintf(MCindex,"%s" , modnam[iup]);      // use if chi2 is not needed in the legend
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
  padfun2->SetGridy();   //Horizontal grid
  padfun2->Draw();
  padfun2->cd();
  gPad->SetTickx();
  gPad->SetTicky();
  //gStyle->SetErrorX(0.5);
  for(int ilow =0 ; ilow <  Nplot[0] ; ilow++){  //loop for ratio plot
  	TH1D *rh2;
    	//if(data->GetBinContent(1) > 0 ||  data->GetBinContent(10) > 0) {
    		//if(data->GetBinContent(10) > 0) {
    			rh2 = (TH1D*)MC[ilow]->Clone("rh2");
    			//rh2->Sumw2();
    			rh2->Divide(data);     //MC divide by data
    /* 		}else{
       			rh2 = (TH1D*)data->Clone("rh2");
       			rh2->Sumw2();
       			rh2->SetLineColor(MC[ilow]->GetLineColor());
       			rh2->SetLineWidth(MC[ilow]->GetLineWidth());
       			rh2->SetLineStyle(MC[ilow]->GetLineStyle());
       			rh2->Divide(MC[ilow]);
       			}*/
    		rh2->SetTitle("");
    		//rh2->SetLineColor(kBlack);
    	//rh2->SetMaximum(lowYrange[0]);  // .. range
    	//rh2->SetMinimum(lowYrange[1]);  // Define Y ..
    	//cout << "max = " << plegend[5] << "  min = " << plegend[6] << endl;
    	rh2->SetMinimum(plegend[6]);  // Define Y ..
    	rh2->SetMaximum(plegend[5]);      // .. range
    	rh2->SetStats(0);                 //No statistics on lower plot
    	//rh2->Divide(data);              //MC divide by data
    	//rh2->SetMarkerStyle(21);
    	//rh2->Draw(" same e1");
    	rh2->Draw("same P e1");
    	//MC_inputerr->Draw("same E2");   //Draw the Uncertainty in Ratio plot

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
}
//end of ratio plot function
//----------------------------------------------------
TH1D* ReadHist1D(string name, TFile* root, int irbin=1){
        TString histname = name; cout << histname<< endl;
        TH1D* hist=(TH1D*)root->Get(histname);
        //hist->Rebin(irbin);
        //hist->Write();
return hist;
}
//----------------------------------------------------
TH1D* ReadHist1D_v1(string name, TFile* root, int irbin=1){
        TString histname = name; //cout << histname<< endl;
        TH1D* hist=(TH1D*)root->Get(histname);
        //hist->Rebin(irbin);
        //hist->Write();
return hist;
}
//----------------------------------------------------
TH2D* ReadHist2D(string name, TFile* root, int irbin=1){
        TString histname = name; cout << histname<<endl;
        TH2D* hist=(TH2D*)root->Get(histname);
        //hist->RebinY(irbin);
        //hist->Write();
return hist;
}
//----------------------------------------------------
