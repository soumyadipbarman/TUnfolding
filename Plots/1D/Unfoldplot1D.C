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

#define CLOUSER

void Unfoldplot1D(){
  
  static const int unfold_ty =1;   //Unfold method to be plots 
  static const int nHLTmx=10;       //HT2 Range
  static const int nmc = 3;
    
  static const int ndef=3;
  static const int njet=2;
  static const int nkappa=10;

  TH1D *MC_gen[nmc][ndef][njet][nkappa][nHLTmx];
  TH1D *MC_reco[nmc][ndef][njet][nkappa][nHLTmx];
  TH2D *MC_Res[nmc][ndef][njet][nkappa][nHLTmx];
  TH1D *Data_reco[ndef][njet][nkappa][nHLTmx];
  TH1D *Psudo_Data_gen[ndef][njet][nkappa][nHLTmx];
  
  TH1D *hist_eff[nmc][ndef][njet][nkappa][nHLTmx];
  TH1D *hist_fake[nmc][ndef][njet][nkappa][nHLTmx];
  TH1D *hist_purity[nmc][ndef][njet][nkappa][nHLTmx];
  TH1D *hist_stbl[nmc][ndef][njet][nkappa][nHLTmx];
  
  TH1D *Unfold[unfold_ty][ndef][njet][nkappa][nHLTmx];
  TH1D *Refold[unfold_ty][ndef][njet][nkappa][nHLTmx];
  TH2D *Corr[unfold_ty][ndef][njet][nkappa][nHLTmx];
  TH2D *Prob[unfold_ty][ndef][njet][nkappa][nHLTmx];
  TH2D *Ematrix[unfold_ty][ndef][njet][nkappa][nHLTmx];
  
  char histname[100],Title[100], Xaxis[100], Yaxis[100], ratioXaxis[100], ratioYaxis[100],pdfname[100],pdfname1[100],pdfname2[100],LegName[100];
  Int_t color[10] ={1,2,4,5,6,46,3,28,38,42};  // define the color for different histograms
  Int_t HT2range[nHLTmx+1]={92, 119, 185, 251, 319, 388, 467, 518, 579, 669, 3000};
  const int njetetamn=1;  //eta value used 2.5

  const char* obs_def[3]={"Q","Q_{L}","Q_{T}"};
  const char* jet_num[2]={"Leading-Jet","Sub-Leading-Jet"};
  const char* k_fact[10]={"k=0.1","k=0.2","k=0.3","k=0.4","k=0.5","k=0.6","k=0.7","k=0.8","k=0.9","k=1.0"};
  
  const char* htrang[10]={"92 < P_{T} < 119", "119 < P_{T} < 185", "185 < P_{T} < 251", "251 < P_{T} < 319", "319 < P_{T} < 388", "388 < P_{T} <467", "467 < P_{T} <518","518 < P_{T} < 579", "579 < P_{T} < 669", "P_{T} > 669"};
  const char* obs_logy[10]={"1/N dN/d","1/N dN/d","1/N dN/d","1/N dN/d","1/N dN/d","1/N dN/d","1/N dN/d","1/N dN/d","1/N dN/d","1/N dN/d"};

  const char* regN[4]={"Tunfold_Noreg","Tunfold_lscan_typ_","Tunfold_scantau_typ_","Tunfold_SURE_typ_"};
  const char* reg_refold[4]={"Tunfold_NoReg_Refold","Tunfold_lscan_Refold","Tunfold_scantau_Refold","Tunfold_SURE_Refold"};
  const char* CORR[4]={"Tunfold_Noreg_corr","Tunfold_lscan_corr","Tunfold_scantau_corr","Tunfold_SURE_corr"};
  const char* COVN[4]={"Tunfold_Noreg_Emat","Tunfold_lscan_Emat","Tunfold_scantau_Emat","Tunfold_SURE_Emat"};
  const char* ProbN[4]={"Tunfold_Noreg_probM","Tunfold_lscan_probM","Tunfold_scantau_probM","Tunfold_SURE_probM"};
  const char* dirname[3]={"Pythia8","MG8","HW7"};

  const char* Validity_test[4]={"Closure test","Bottom Line test"," Unfolded","Refold"};
  const char* h2dMat_name[4]={"Covariance matrix","correlation coefficients"," probabilities matrix ","Response matrix"};
  const char* mcname[3]={"Pythia8 CP5 Tune","Madgraph","Herwig++"};
  const char* mcnamerco[3]={"Pythia8 RECO","Madgraph RECO","Herwig++ RECO"};
  const char* mcnamegen[3]={"Pythia8 GEN","Madgraph GEN","Herwig++ GEN"};
  const char* DataEra[3]={"Data RECO","Data RECO","Data RECO"};
  const char* UndoldEra[3]={"Unfold 2016","Unfold 2017","Unfold 2018"};
  const char* RefoldEra[5]={"Refold Pythia8(No Regularisation) ","Refold Pythia8(L-Curve scan)","Refold Pythia8(Scan Tau)","Refold Pythia8(ScanSURE)","Refold Pythia8(Iterative method)"};
 //const char* RefoldEra[5]={"Refold Madgraph(No Regularisation) ","Refold Madgraph(L-Curve scan)","Refold Madgraph(Scan Tau)","Refold Madgraph(ScanSURE)","Refold Pythia8(Iterative method)"};
  const char* Unfoldtype[5]={"TUnfold(No Regularisation) ","TUnfold(L-Curve scan)","TUnfold(Scan Tau)","TUnfold(ScanSURE)","TUnfold(Iterative method)"};
  const char* closuretype[5]={"Unfolded Pythia8(No Regularisation) ","Unfolded Pythia8 (L-curve scan)","Unfolded Pythia8(Scan Tau)","Unfolded Pythia8 (ScanSURE)","Unfolded(Iterative method)"};
  //const char* closuretype[5]={"Unfolded Madgraph(No Regularisation) ","Unfolded Madgraph (L-curve scan)","Unfolded Madgraph(Scan Tau)","Unfolded Madgraph (ScanSURE)","Unfolded(Iterative method)"};
  const char* Modelnm[3]={"Pythia8","Madgraph","Herwig"};
  const char* Methodtype[5]={"(No Regularisation) ","(L-Curve scan)","(Scan Tau)","(ScanSURE)","(Iterative method)"};
  const char* smeared[5]={"TUnfold","Refold","Folded-back","GEN","RECO"};
  static const int iera = 1;
  int iPeriod = 0;  
  int iPos=10 ;
  
  //input root file
  TFile *Unfoldroot = TFile::Open("/home/soumyadip/Package/TUnfold/JetCharge/TUnfold1D/Unfolded_Result.root");  // Unfolded data 

//----------------------------------------------------
//Function declear
  void Integralhist(TH1D *hist);
  void divBybinWidth(TH1D *hist);
  void Myplotset(TH1D *Myhist,const char* XTitle, const char* YTitle);
  void Set2dHist(TH2D *MyHist, const char* XTitle, const char* YTitle,const char* ZTitle, double titoff[3], double titsize[3]);
  void SetMycanvas(TCanvas *cpt,double bs,double lm, double rm, double tm,double bm);
  void CTLegend(TLegend *legendn, const char* txt1, const char* txt2);
  TCanvas *ratio_can(int Nplot[2],float plegend[7], TH1D* data, TH1D* MC[Nplot[0]], char* lowpadx,const  char* modnam[Nplot[0]],const  char* datanm[1]);
  TCanvas *ratio_can1(int Nplot[3],float plegend[7], TH1D* data, TH1D* MC[Nplot[0]], char* lowpadx, const char* modnam[Nplot[0]], const  char* datanm[3]);
//----------------------------------------------------
	for(int id=0; id<ndef; id++){
		for(int ij=0; ij<njet; ij++){
      			for(int ik =0 ; ik<nkappa ; ik++){
        			for(Int_t ipt =0; ipt < nHLTmx ; ipt++){     
#ifdef CLOUSER
          		sprintf(histname, "Data/gen_jc_d%i_j%i_k%i_pt%i_eta0",id, ij, ik, ipt); //gen_jc_d0_j0_k0_pt0_eta0
#else
          		sprintf(histname, "%s/gen_jc_d%i_j%i_k%i_pt%i_eta0",dirname[0], id, ij, ik, ipt); //gen_jc_d0_j0_k0_pt0_eta0
#endif
          		Psudo_Data_gen[id][ij][ik][ipt] = (TH1D*) Unfoldroot->Get(histname);
          		//Psudo_Data_gen[ity][ivar][ipt]->Rebin(2);  // check rebin part
          		cout << histname << endl;
 			//Exclude Underflow overlow
 		        TH1D *NewData;
          		NewData = (TH1D*)Psudo_Data_gen[id][ij][ik][ipt]->Clone();
	  		int NbinxD = NewData->GetNbinsX();
          		NewData->Reset();
	  		for(int ix=1; ix < NbinxD+1 ; ix++){
	  			NewData->SetBinContent(ix,Psudo_Data_gen[id][ij][ik][ipt]->GetBinContent(ix));
	  			NewData->SetBinError(ix, sqrt(Psudo_Data_gen[id][ij][ik][ipt]->GetBinError(ix)*Psudo_Data_gen[id][ij][ik][ipt]->GetBinError(ix)));
       					}
				}
      			}
    		}
	}
//----------------------------------------------------
  	for(int  imc =0; imc < nmc ; imc++){
    		for(int id=0; id <ndef; id++){
      			for(int ij =0 ; ij < njet ; ij++){
				for (int ik=0; ik<nkappa; ik++){
					for(Int_t ipt =0; ipt < nHLTmx ; ipt++){     
	  		sprintf(histname, "%s/gen_jc_d%i_j%i_k%i_pt%i_eta0",dirname[imc], id, ij, ik, ipt); //reco_typ_1_pt6_eta0_15
	  		MC_gen[imc][id][ij][ik][ipt] = (TH1D*) Unfoldroot->Get(histname);
	  		//MC_gen[imc][id][ij][ik][ipt]->Rebin(2);   // check rebin 

	  		sprintf(histname, "%s/reco_jc_d%i_j%i_k%i_pt%i_eta0", dirname[imc], id, ij, ik,  ipt); //reco_typ_1_pt6_eta0_15
	  		MC_reco[imc][id][ij][ik][ipt] = (TH1D*) Unfoldroot->Get(histname);
	  		cout << histname << endl;
			//Read Reso
			sprintf(histname, "%s/RM_jc_d%i_j%i_k%i_pt%i_eta0",dirname[imc], id, ij, ik,  ipt); //reco_typ_1_pt6_eta0_15
	  		MC_Res[imc][id][ij][ik][ipt] = (TH2D*) Unfoldroot->Get(histname);
	  		cout << histname << endl;

	  		sprintf(histname,"Unfold/genmiss_jc_d%i_j%i_k%i_pt%i_eta0", id, ij, ik, ipt);
	  		hist_eff[imc][id][ij][ik][ipt] = (TH1D*) Unfoldroot->Get(histname);
	  		sprintf(histname,"Unfold/recofake_jc_d%i_j%i_k%i_pt%i_eta0", id, ij, ik, ipt);
	  		hist_fake[imc][id][ij][ik][ipt] = (TH1D*) Unfoldroot->Get(histname);
	  		sprintf(histname,"Unfold/Purity1_d%i_j%i_k%i_pt%i_eta0", id, ij, ik, ipt);
	  		hist_purity[imc][id][ij][ik][ipt] = (TH1D*) Unfoldroot->Get(histname);
	  		sprintf(histname,"Unfold/stability_d%i_j%i_k%i_pt%i_eta0", id, ij, ik, ipt);
	  		hist_stbl[imc][id][ij][ik][ipt] = (TH1D*) Unfoldroot->Get(histname);
		  		}
      			}
    		}
	}
}//end of one MCinput root file  reading
//----------------------------------------------------
//Read Data
	for(int id=0; id <ndef; id++){
    		for(int ij =0 ; ij < njet ; ij++){
			for (int ik=0; ik<nkappa; ik++){
      				for(Int_t ipt =0; ipt < nHLTmx ; ipt++){
			sprintf(histname, "Data/reco_jc_d%i_j%i_k%i_pt%i_eta0", id, ij, ik, ipt); //reco_typ_1_pt6_eta0_15
			Data_reco[id][ij][ik][ipt] = (TH1D*) Unfoldroot->Get(histname);
			cout << histname << endl;     
      				}
			}
    		}	 
	}
//----------------------------------------------------
//Read Unfolded
  	for(int iun=0; iun < unfold_ty; iun++){
    		for(int id=0; id <ndef; id++){
      			for(int ij=0; ij < njet ; ij++){
				for(int ik=0; ik<nkappa; ik++){
					for(int ipt = 0 ; ipt < nHLTmx ; ipt++){
	  		sprintf(histname, "Unfold/%s_d%i_j%i_k%i_pt%i_eta0",regN[iun], id, ij, ik, ipt); //Tunfolded_typ_0_pt0_eta0_3
	  		Unfold[iun][id][ij][ik][ipt] = (TH1D*) Unfoldroot->Get(histname);
	  		cout << histname << endl;
	  		//Exclude Underflow and overflow bins
	  		TH1D *NewData;
          		NewData = (TH1D*) Unfold[iun][id][ij][ik][ipt]->Clone();
          		int NbinxD = NewData->GetNbinsX();
          		NewData->Reset();
          		for(int ix=1; ix < NbinxD+1 ; ix++){
          		NewData->SetBinContent(ix, Unfold[iun][id][ij][ik][ipt]->GetBinContent(ix));
          		NewData->SetBinError(ix, sqrt( Unfold[iun][id][ij][ik][ipt]->GetBinError(ix)* Unfold[iun][id][ij][ik][ipt]->GetBinError(ix)));
       			}
	 
	  		sprintf(histname, "Unfold/%s_d%i_j%i_k%i_pt%i_eta0",CORR[iun], id, ij, ik, ipt); //Tunfolded_typ_0_pt0_eta0_3
			//sprintf(histname, "Unfold/Tunfold_Noreg_corr_d%i_j%i_k%i_pt%i_eta0",id, ij, ik, ipt);
          		Corr[iun][id][ij][ik][ipt]=(TH2D*) Unfoldroot->Get(histname);
	  		sprintf(histname, "Unfold/%s_d%i_j%i_k%i_pt%i_eta0",COVN[iun], id, ij, ik, ipt); //Tunfolded_typ_0_pt0_eta0_3
			//sprintf(histname, "Unfold/Tunfold_Noreg_Emat_d%i_j%i_k%i_pt%i_eta0",COVN[iun], id, ij, ik, ipt); //Tunfolded_typ_0_pt0_eta0_3
           		Ematrix[iun][id][ij][ik][ipt]= (TH2D*) Unfoldroot->Get(histname);
	  		sprintf(histname, "Unfold/%s_d%i_j%i_k%i_pt%i_eta0",ProbN[iun], id, ij, ik, ipt); //Tunfolded_typ_0_pt0_eta0_3
	   		Prob[iun][id][ij][ik][ipt] = (TH2D*) Unfoldroot->Get(histname);
				}
      			}
    		} //unfolded histUnfold
	}
}
//----------------------------------------------------
//Read Folded Back
  	for(int iun=0; iun < unfold_ty; iun++){
    		for(int id=0; id <ndef; id++){
      			for(int ij=0; ij < njet ; ij++){
				for(int ik=0; ik<nkappa; ik++){
					for(int ipt = 0 ; ipt < nHLTmx ; ipt++){
	  		sprintf(histname, "Unfold/%s_d%i_j%i_k%i_pt%i_eta0",reg_refold[iun], id, ij, ik, ipt); //Tunfolded_typ_0_pt0_eta0_3
	  		Refold[iun][id][ij][ik][ipt] = (TH1D*) Unfoldroot->Get(histname);
	  		cout << histname << endl;
				}
      			}
		}
	}//unfolded histUnfold
}
cout <<"Read histogram OK" << endl;
//----------------------------------------------------
//PLot canvas declear
  TCanvas *cpt0 = new TCanvas("cpt0", "canvas0", 700,600 );  //for Reco
  TCanvas *cpt5 = new TCanvas("cpt5", "canvas5", 600,575 );  //for Corr
  TCanvas *cpt6 = new TCanvas("cpt6", "canvas6", 600,575 );  //for Prob
  TCanvas *cpt7 = new TCanvas("cpt7", "canvas7", 600,575 );  //for COV
  TCanvas *cpt8 = new TCanvas("cpt8", "canvas8", 600,575 );  //for Response
  TCanvas *cpt9 = new TCanvas("cpt9", "canvas9", 600,575 );  //for Projection
//----------------------------------------------------
//Reco comparison 
  	for(int id=0; id <ndef; id++){
    		for(int ij =0 ; ij < njet ; ij++){
			for(int ik=0; ik<nkappa; ik++){
      				for(int ipt = 0 ; ipt <nHLTmx; ipt++){
			sprintf(histname, "Data/reco_jc_d%i_j%i_k%i_pt%i_eta0", id, ij, ik, ipt); //reco_typ_1_pt6_eta0_15
                        Data_reco[id][ij][ik][ipt] = (TH1D*) Unfoldroot->Get(histname);
                        cout << histname << endl;
			TH1D *MyHist  = (TH1D*) Data_reco[id][ij][ik][ipt]->Clone();
			Integralhist(MyHist);
			//divBybinWidth(MyHist);
			Myplotset(MyHist,0,0);
	
			if(ipt<=8){sprintf(Title,"%s_{%s}^{%s}:     %i <P_{T}< %i %s", obs_def[id], jet_num[ij], k_fact[ik], HT2range[ipt] , HT2range[ipt+1] ,"GeV/c" );}
			else if(ipt==9){ sprintf(Title,"%s_{%s}^{%s}:     <P_{T}> %i %s", obs_def[id], jet_num[ij], k_fact[ik], HT2range[ipt] ,"GeV/c"  );}
				sprintf(Yaxis," %s %s^{%s}" ,obs_logy[ik], obs_def[id], k_fact[ik]);
				MyHist->SetTitle(Title);
				MyHist->GetXaxis()->SetTitle("");
				MyHist->GetYaxis()->SetTitle(Yaxis);
	
			TH1D *MC_input[nmc];
			const char *MCinput_index[nmc+1];
			const char *data_index[1];
			for(int iout = 0 ; iout < nmc ; iout++){
	  			MC_input[iout] = (TH1D*) MC_reco[iout][id][ij][ik][ipt]->Clone(); 
	  			Integralhist(MC_input[iout]);
	 			//divBybinWidth(MC_input[iout]);
	  			MCinput_index[iout]= mcnamerco[iout]; }
	  			data_index[0]= DataEra[1]; 
	
			char lplot_xtitle[100];
			//sprintf(lplot_xtitle, "%s",Esvlogx[ivar]);
			sprintf(lplot_xtitle," %s %s^{%s} %d" ,jet_num[ij],obs_def[id],k_fact[ik],HT2range[ipt]);
			//float ratio_range1[2]={1.2,0.9};
			int num1[2]={nmc,1} ;
			float lpos1[7] ={.32,0.2,0.55,0.38, .04, 1.5,0.7};
	
			cpt0->cd();
			SetMycanvas(cpt0,0,0,0,0,0);

			cpt0 =(TCanvas*)(ratio_can(num1, lpos1, MyHist, MC_input, lplot_xtitle, MCinput_index,data_index));
			CMS_lumi( cpt0, iPeriod, iPos ); cpt0->Update();
	
			sprintf(pdfname, "RecoEVS_Plot.pdf("); sprintf(pdfname1, "RecoEVS_Plot.pdf"); sprintf(pdfname2, "RecoEVS_Plot.pdf)"); //TData_recoed_typ_0_pt0_eta0_3
			if(id==0 && ij==0 && ik==0 && ipt ==0){cpt0->Print(pdfname,"pdf");
			}else if(id==2 && ij==1 && ik==9 && ipt==9) {cpt0->Print(pdfname2,"pdf");
			}else{  cpt0->Print(pdfname1,"pdf");};
			}
      		}
   	 }//end of phase space cut and variable loop
}
//end of RECO PLOT    
//----------------------------------------------------
//Reco Projection comparison
  	for(int id=0; id <ndef; id++){
		for (int ij=0; ij<njet; ij++){
    			for(int ik =0 ; ik < nkappa ; ik++){
      				for(int ipt = 0 ; ipt <nHLTmx; ipt++){
      			sprintf(histname, "Pythia8/%s_d%i_j%i_k%i_pt%i_eta0","ProjectX",id, ij, ik, ipt); //Tunfolded_typ_0_pt0_eta0_3
      			TH1* RMx=(TH1D*) Unfoldroot->Get(histname);
      			sprintf(histname, "Pythia8/%s_d%i_j%i_k%i_pt%i_eta0","Recominusfake",id, ij, ik, ipt); //Tunfolded_typ_0_pt0_eta0_3
      			TH1* RecoFake=(TH1D*) Unfoldroot->Get(histname);
      			cpt9->cd();
      			SetMycanvas(cpt9,0,0.1,0.15,0.05,0.12);
      
      			TH1D *MyHist  = (TH1D*)RMx->Clone();
      			Integralhist(MyHist);
        		//divBybinWidth(MyHist);
        		Myplotset(MyHist,0,0);

        		if(ipt<=8){sprintf(Title,"%s_{%s}^{%s}:     %i <P_{T}< %i %s", obs_def[id], jet_num[ij], k_fact[ik], HT2range[ipt] , HT2range[ipt+1] ,"GeV/c"  );}
        		else if(ipt==9){ sprintf(Title,"%s_{%s}^{%s}:      <P_{T}> %i %s", obs_def[id], jet_num[ij], k_fact[ik], HT2range[ipt] ,"GeV/c"  );}
        		sprintf(Yaxis," %s %s^{%s}" ,obs_logy[ik], obs_def[id], k_fact[ik]);
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
        			sprintf(lplot_xtitle," %s %s^{%s} %d" ,jet_num[ij],obs_def[id],k_fact[ik],HT2range[ipt]);
        			//float ratio_range1[2]={1.2,0.9};
        			int num1[2]={imc,1} ;
        			float lpos1[7] ={.32,0.2,0.55,0.38, .04, 1.15,0.85};

        			cpt9->Clear();
        			cpt9->cd();
        			SetMycanvas(cpt9,0,0,0,0,0);
	
        			cpt9 =(TCanvas*)(ratio_can(num1, lpos1, MyHist, MC_input, lplot_xtitle, MCinput_index,data_index));
        			CMS_lumi( cpt9, iPeriod, iPos ); cpt9->Update();

      				sprintf(pdfname, "RM_ProjectX_Plot.pdf("); sprintf(pdfname1, "RM_ProjectX_Plot.pdf"); sprintf(pdfname2, "RM_ProjectX_Plot.pdf)"); //TData_recoed_typ_0_pt0_eta0_3
        			if(id==0 && ij==0 && ik==0 && ipt ==0){cpt9->Print(pdfname,"pdf");  // check ??
        			}else if(id==2 && ij==1 && ik==9 && ipt==9) {cpt9->Print(pdfname2,"pdf");
        			}else{  cpt9->Print(pdfname1,"pdf");};
				}
      			}
    		}
  	}
//----------------------------------------------------
//Gen Projection comparison
  	for(int id=0; id <ndef; id++){
    		for(int ij =0 ; ij < njet ; ij++){
			for(int ik=0; ik<nkappa; ik++){
      				for(int ipt = 0 ; ipt <nHLTmx; ipt++){
      			sprintf(histname, "Pythia8/%s_d%i_j%i_k%i_pt%i_eta0","ProjectY",id, ij, ik, ipt);     //Tunfolded_typ_0_pt0_eta0_3
      			TH1* RMx=(TH1D*) Unfoldroot->Get(histname);
      			sprintf(histname, "Pythia8/%s_d%i_j%i_k%i_pt%i_eta0","Genminusmiss",id, ij, ik, ipt); //Tunfolded_typ_0_pt0_eta0_3
      			TH1* GenMiss=(TH1D*) Unfoldroot->Get(histname);
      			cpt9->cd();
      			SetMycanvas(cpt9,0,0.1,0.15,0.05,0.12);

      			TH1D *MyHist  = (TH1D*)RMx->Clone();
      			Integralhist(MyHist);
        		//divBybinWidth(MyHist);
        		Myplotset(MyHist,0,0);

        		if(ipt<=8){sprintf(Title,"%s_{%s}^{%s}:     %i <P_{T}< %i %s", obs_def[id], jet_num[ij], k_fact[ik], HT2range[ipt] , HT2range[ipt+1] ,"GeV/c" );}
        		else if(ipt==9){sprintf(Title,"%s_{%s}^{%s}:      <P_{T}> %i %s", obs_def[id], jet_num[ij], k_fact[ik], HT2range[ipt],"GeV/c");}
        		sprintf(Yaxis," %s %s^{%s}" ,obs_logy[ik], obs_def[id], k_fact[ik]);
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
        				sprintf(lplot_xtitle," %s %s^{%s} %d" ,jet_num[ij],obs_def[id],k_fact[ik],HT2range[ipt]);
        				//float ratio_range1[2]={1.2,0.9};
        				int num1[2]={imc,1};
        				float lpos1[7] ={.32,0.2,0.55,0.38, .04, 1.15,0.85};

        				cpt9->Clear();
        				cpt9->cd();
        				SetMycanvas(cpt9,0,0,0,0,0);

        				cpt9 =(TCanvas*)(ratio_can(num1, lpos1, MyHist, MC_input, lplot_xtitle, MCinput_index,data_index));
        				CMS_lumi( cpt9, iPeriod, iPos ); cpt9->Update();

	        			sprintf(pdfname, "RM_ProjectY_Plot.pdf("); sprintf(pdfname1, "RM_ProjectY_Plot.pdf"); sprintf(pdfname2, "RM_ProjectY_Plot.pdf)"); //TData_recoed_typ_0_pt0_eta0_3
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
      			sprintf(histname, "fold/%s_d%i_j%i_k%i_pt%i_eta0","Fold",id, ij, ik, ipt); //Tunfolded_typ_0_pt0_eta0_3
      			TH1D* GenFold=(TH1D*)Unfoldroot->Get(histname);
        
        		TH1D *MyHist  = (TH1D*)GenFold->Clone();
        		Integralhist(MyHist);
        		//divBybinWidth(MyHist);
        		Myplotset(MyHist,0,0);

        		if(ipt<=8){sprintf(Title,"%s_{%s}^{%s}:     %i <P_{T}< %i %s", obs_def[id], jet_num[ij], k_fact[ik], HT2range[ipt] , HT2range[ipt+1] ,"GeV/c" );}
                        else if(ipt==9){sprintf(Title,"%s_{%s}^{%s}:      <P_{T}> %i %s", obs_def[id], jet_num[ij], k_fact[ik], HT2range[ipt] ,"GeV/c");}
                        sprintf(Yaxis," %s %s^{%s}" ,obs_logy[ik], obs_def[id], k_fact[ik]);
        		MyHist->SetTitle(Title);
        		MyHist->GetXaxis()->SetTitle("");
        		MyHist->GetYaxis()->SetTitle(Yaxis);
        		int imc =1;
        			TH1D *MC_input[imc];
        			const char *MCinput_index[imc+1];
        			const char *data_index[1];
        			for(int iout = 0 ; iout < imc ; iout++){
          				MC_input[iout] = (TH1D*) MC_reco[iout][id][ij][ik][ipt]->Clone();
          				Integralhist(MC_input[iout]);
         				//divBybinWidth(MC_input[iout]);
          				MCinput_index[iout]= "Pythia RECO"; }
          				data_index[0]= "Folded";

        				char lplot_xtitle[100];
        				sprintf(lplot_xtitle," %s %s^{%s} %d" ,jet_num[ij],obs_def[id],k_fact[ik],HT2range[ipt]);
        				//float ratio_range1[2]={1.2,0.9};
        				int num1[2]={imc,1} ;
        				float lpos1[7] ={.32,0.2,0.55,0.38, .04, 1.5,0.7};
        
					cpt9->Clear();
        				cpt9->cd();
        				SetMycanvas(cpt9,0,0,0,0,0);

        				cpt9 =(TCanvas*)(ratio_can(num1, lpos1, MyHist, MC_input, lplot_xtitle, MCinput_index,data_index));
        				CMS_lumi( cpt9, iPeriod, iPos ); cpt9->Update();

        				sprintf(pdfname, "Genfold_Plot.pdf("); sprintf(pdfname1, "Genfold_Plot.pdf"); sprintf(pdfname2, "Genfold_Plot.pdf)"); //TData_recoed_typ_0_pt0_eta0_3
        				if(id==0 && ij==0 && ik==0 && ipt ==0){cpt9->Print(pdfname,"pdf");
        				}else if(id==2 && ij==1 && ik==9 && ipt==9) {cpt9->Print(pdfname2,"pdf");
        				}else{cpt9->Print(pdfname1,"pdf");};
				}
      			}
    		}//end of phase space cut and variable loop
  	}	 
//----------------------------------------------------
//Response Matrix
	for(int imc =0; imc < nmc ; imc++){
  		for(int id=0; id <ndef; id++){
			for(int ij=0; ij<njet; ij++){
    				for(int ik =0 ; ik < nkappa ; ik++){
      					for(int ipt = 0 ; ipt <nHLTmx; ipt++){
			cpt8->cd();

			//MC_Res[imc][ity][ivar][ipt]->RebinY(2);
        		char lplot_xtitle[800]; char lplot_ytitle[800];
			sprintf(lplot_xtitle, "RECO     %s %s^{%s}" ,jet_num[ij],obs_def[id],k_fact[ik]); 
			sprintf(lplot_ytitle, "GEN      %s %s^{%s}" ,jet_num[ij],obs_def[id],k_fact[ik]);	
			double titoff1[3]={1.2,1.3,1.0};
        		double titsize1[3] ={0.035,0.035,0.035};
	
			SetMycanvas(cpt8,0,0.1,0.15,0.05,0.1);
        		Set2dHist( MC_Res[imc][id][ij][ik][ipt],lplot_xtitle, lplot_ytitle,"",titoff1, titsize1);
			MC_Res[imc][id][ij][ik][ipt]->Draw("colz");
         
			if(ipt<=8){sprintf(Title,"%s_{%s}^{%s}:     %i <P_{T}< %i %s", obs_def[id], jet_num[ij], k_fact[ik], HT2range[ipt] , HT2range[ipt+1] ,"GeV/c" );}
                        else if(ipt==9){sprintf(Title,"%s_{%s}^{%s}:      <P_{T}> %i %s", obs_def[id], jet_num[ij], k_fact[ik], HT2range[ipt] ,"GeV/c");}

			TLegend *leg1 = new TLegend(0.05,0.7,0.4,0.8);     
			CTLegend(leg1,Modelnm[imc],Title); leg1->AddEntry((TObject*)0,obs_def[id] , "");leg1->SetTextColor(-8);leg1->Draw();
        		CMS_lumi( cpt8, iPeriod, iPos ); cpt8->Update();
			sprintf(pdfname, "Response_Mat%i.pdf(",imc); sprintf(pdfname1, "Response_Mat%i.pdf",imc); sprintf(pdfname2, "Response_Mat%i.pdf)",imc); //TData_recoed_typ_0_pt0_eta0_3
        		if(id==0 && ij==0 && ik==0 && ipt ==0){cpt8->Print(pdfname,"pdf");
        		}else if(id==2 && ij==1 && ik==9 && ipt==9) {cpt8->Print(pdfname2,"pdf");
        		}else{cpt8->Print(pdfname1,"pdf");};
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
  				TCanvas *cpt1 = new TCanvas("cpt1", "canvas1", 600,600 );  //for 
  				TCanvas *cpt2 = new TCanvas("cpt2", "canvas2", 600,600 );  //for ESVs
  				TCanvas *cpt3 = new TCanvas("cpt3", "canvas3", 600,600 );  //for ESVs
  				TCanvas *cpt4 = new TCanvas("cpt4", "canvas4", 800,800 );  //for 

   				TLegend *leg2 = new TLegend(0.4,0.5,0.7,0.8);
				CTLegend(leg2," ","");
				TLegend *leg1 = new TLegend(0.1,0.6,0.4,0.8);
        			CTLegend(leg1,"", obs_def[id]); 
	
				for(int ipt = 0 ; ipt <nHLTmx; ipt++){
        				cpt1->cd();

					if(ipt<=8){sprintf(Title,"%s_{%s}^{%s}:     %i <P_{T}< %i %s", obs_def[id], jet_num[ij], k_fact[ik], HT2range[ipt] , HT2range[ipt+1] ,"GeV/c" );}
                        		else if(ipt==9){sprintf(Title,"%s_{%s}^{%s}:      <P_{T}> %i %s", obs_def[id], jet_num[ij], k_fact[ik], HT2range[ipt],"GeV/c");}

        				for (int i = 1; i <= hist_eff[0][id][ij][ik][ipt]->GetNbinsX(); ++i) {
         					double content = 1 - hist_eff[0][id][ij][ik][ipt]->GetBinContent(i);
         					hist_eff[0][id][ij][ik][ipt]->SetBinContent(i, content);
       						}

        			SetMycanvas(cpt1,0,0.1,0.15,0.05,0.12);
        			Myplotset(hist_eff[0][id][ij][ik][ipt],obs_def[id],"Efficiency");
        			hist_eff[0][id][ij][ik][ipt]->SetLineColor(color[ipt]);
				hist_eff[0][id][ij][ik][ipt]->Draw("same"); leg1->Draw();
   
        			leg2->AddEntry(hist_eff[0][id][ij][ik][ipt], Title ,"lp");
				if(ipt==7){leg2->Draw();}

        			cpt2->cd();
				//SetMycanvas(cpt2);
        			SetMycanvas(cpt2,0,0.1,0.15,0.05,0.12);
 				Myplotset(hist_purity[0][id][ij][ik][ipt],obs_def[id],"Purity");
 				hist_purity[0][id][ij][ik][ipt]->SetLineColor(color[ipt]);
        			hist_purity[0][id][ij][ik][ipt]->Draw("same"); leg1->Draw();
	 			if(ipt==7){leg2->Draw();}
				cpt2->Update();

        			cpt3->cd();
        			SetMycanvas(cpt3,0,0.1,0.1,0.05,0.12);
        			Myplotset(hist_fake[0][id][ij][ik][ipt],obs_def[id],"Fake rate");
        			hist_fake[0][id][ij][ik][ipt]->SetLineColor(color[ipt]);
        			hist_fake[0][id][ij][ik][ipt]->Draw("same");leg1->Draw();
        			if(ipt==7){leg2->Draw();}
         			cpt3->Update();
				cpt4->cd();
        			SetMycanvas(cpt4,0,0.1,0.1,0.05,0.12);
        			Myplotset(hist_stbl[0][id][ij][ik][ipt],obs_def[id],"Stability");
        			hist_stbl[0][id][ij][ik][ipt]->SetLineColor(color[ipt]); 
        			hist_stbl[0][id][ij][ik][ipt]->Draw("same"); leg1->Draw();
        			if(ipt==7){leg2->Draw();}
   	 			cpt4->Update();
    				}
    				sprintf(pdfname, "effi_plot.pdf("); sprintf(pdfname1, "effi_plot.pdf"); sprintf(pdfname2, "effi_plot.pdf)");
        			if(id==0 && ij==0 && ik==0){cpt1->Print(pdfname,"pdf");
        			}else if(id==2 && ij==1 && ik==9 ) {cpt1->Print(pdfname2,"pdf");
        			}else{  cpt1->Print(pdfname1,"pdf");};
				cpt1->Clear();
    				sprintf(pdfname, "puri_plot.pdf("); sprintf(pdfname1, "puri_plot.pdf"); sprintf(pdfname2, "puri_plot.pdf)"); 
        			if(id==0 && ij==0 && ik==0 ){cpt2->Print(pdfname,"pdf");
        			}else if(id==2 && ij==1 && ik==9 ) {cpt2->Print(pdfname2,"pdf");
        			}else{  cpt2->Print(pdfname1,"pdf");};
				cpt2->Clear();
     				sprintf(pdfname, "fake_plot.pdf("); sprintf(pdfname1, "fake_plot.pdf"); sprintf(pdfname2, "fake_plot.pdf)"); 
        			if(id==0 && ij==0 && ik==0){cpt4->Print(pdfname,"pdf");
                                }else if(id==2 && ij==1 && ik==9 ) {cpt4->Print(pdfname2,"pdf");
                                }else{  cpt3->Print(pdfname1,"pdf");};
				cpt3->Clear();
     				sprintf(pdfname, "stab_plot.pdf("); sprintf(pdfname1, "stab_plot.pdf"); sprintf(pdfname2, "stab_plot.pdf)"); 
        			if(id==0 && ij==0 && ik==0){cpt4->Print(pdfname,"pdf");
        			}else if(id==2 && ij==1 && ik==9 ) {cpt4->Print(pdfname2,"pdf");
        			}else{  cpt4->Print(pdfname1,"pdf");};
				cpt4->Clear();
			}
     		}
  	}
//----------------------------------------------------
//Unfold comparisoin
  	for(int iun=0; iun < unfold_ty; iun++){
    		for(int id=0; id <ndef; id++){
			for (int ij=0; ij<njet; ij++){
      				for(int ik =0 ; ik < nkappa ; ik++){
					for(int ipt = 0 ; ipt <nHLTmx; ipt++){
	  		TH1D *MyHist  = (TH1D*)Unfold[iun][id][ij][ik][ipt]->Clone();
	  		Integralhist(MyHist);
	  		//divBybinWidth(MyHist);
	  		Myplotset(MyHist,0,0);
	  		if(ipt<=8){sprintf(Title,"%s_{%s}^{%s}:     %i <P_{T}< %i %s", obs_def[id], jet_num[ij], k_fact[ik], HT2range[ipt] , HT2range[ipt+1] ,"GeV/c" );}
                        else if(ipt==9){sprintf(Title,"%s_{%s}^{%s}:      <P_{T}> %i %s", obs_def[id], jet_num[ij], k_fact[ik], HT2range[ipt] ,"GeV/c");}
                        sprintf(Yaxis," %s %s^{%s}" ,obs_logy[ik], obs_def[id], k_fact[ik]);
	  		MyHist->SetTitle(Title);
	  		MyHist->GetXaxis()->SetTitle("");
	  		MyHist->GetYaxis()->SetTitle(Yaxis);
	  
	  		TH1D *MC_input[nmc];
	  		const char *MCinput_index[nmc], *data_index[1];
	  			for(int iout = 0 ; iout < nmc ; iout++){
	    				MC_input[iout] = (TH1D*) MC_gen[iout][id][ij][ik][ipt]->Clone();
	    				Integralhist(MC_input[iout]);
	    				//divBybinWidth(MC_input[iout]);
	    				MCinput_index[iout]= mcnamegen[iout]; }
	 
#ifdef CLOUSER
	 	 	data_index[0]= closuretype[iun]; 
#else
	  		data_index[0]= Unfoldtype[iun]; 
#endif	 
	  		//data_index[0]= UndoldEra[1]; 
	  		char lplot_xtitle[100];
	  		sprintf(lplot_xtitle," %s %s^{%s} %d" ,jet_num[ij],obs_def[id],k_fact[ik],HT2range[ipt]);
	  		//float ratio_range1[2]={1.2,0.9};
	  		int num1[2]={nmc,1} ;
	  		float lpos1[7] ={.32,0.2,0.55,0.38, .04, 1.5,0.7};
	  
	  		cpt0->cd();
          		SetMycanvas(cpt0,0,0,0,0,0);
	  		cpt0 =(TCanvas*)(ratio_can(num1, lpos1, MyHist, MC_input, lplot_xtitle,MCinput_index,data_index));
	  		CMS_lumi( cpt0, iPeriod, iPos ); cpt0->Update();
	  
	  		sprintf(pdfname, "TUnfold_plot_%i.pdf(" ,iun); sprintf(pdfname1, "TUnfold_plot_%i.pdf" ,iun);
			sprintf(pdfname2, "TUnfold_plot_%i.pdf)" ,iun); //Tunfolded_typ_0_pt0_eta0_3
	  		if(id==0 && ij==0 && ik==0 && ipt ==0){cpt6->Print(pdfname,"pdf");
                        }else if(id==2 && ij==1 && ik==9 && ipt==9) {cpt6->Print(pdfname2,"pdf");
                        }else{cpt6->Print(pdfname,"pdf");};
	 
       			cpt5->cd();
       			double titoff1[3]={1.2,1.3,1.0};
       			double titsize1[3] ={0.035,0.035,0.035};
       			SetMycanvas(cpt5,0,0.1,0.15,0.05,0.1);
       			gStyle->SetPaintTextFormat( "4.2f");
       			Set2dHist(Corr[iun][id][ij][ik][ipt],lplot_xtitle, lplot_xtitle,"correlation coefficients", titoff1, titsize1);
       			Corr[iun][id][ij][ik][ipt]->Draw("colz ");
       			//Corr[iun][ity][ivar][ipt]->Draw("colz text");
       			if(ipt<=8){sprintf(Title,"%s_{%s}^{%s}:     %i <P_{T}< %i %s", obs_def[id], jet_num[ij], k_fact[ik], HT2range[ipt] , HT2range[ipt+1] ,"GeV/c" );}
                        else if(ipt==9){sprintf(Title,"%s_{%s}^{%s}:      <P_{T}> %i %s", obs_def[id], jet_num[ij], k_fact[ik], HT2range[ipt] ,"GeV/c");}
                        sprintf(Yaxis," %s %s^{%s}" ,obs_logy[ik], obs_def[id], k_fact[ik]);

        		TLegend *leg1 = new TLegend(0.05,0.6,0.4,0.8);
        		CTLegend(leg1,"Unfolded with Pythia8",Title); leg1->AddEntry((TObject*)0,obs_def[id] , ""); leg1->AddEntry((TObject*)0,Unfoldtype[iun] , "");leg1->SetTextColor(-8);leg1->Draw();
       
       			CMS_lumi( cpt5, iPeriod, iPos ); cpt5->Update();       
       			sprintf(pdfname, "TUnfold_corr_%i.pdf(" ,iun); sprintf(pdfname1, "TUnfold_corr_%i.pdf" ,iun);sprintf(pdfname2, "TUnfold_corr_%i.pdf)" ,iun); //Tunfolded_typ_0_pt0_eta0_3
          		if(id==0 && ij==0 && ik==0 && ipt ==0){cpt6->Print(pdfname,"pdf");
                        }else if(id==2 && ij==1 && ik==9 && ipt==9) {cpt6->Print(pdfname2,"pdf");
                        }else{cpt6->Print(pdfname,"pdf");};

       			cpt6->cd();
       			SetMycanvas(cpt6,0,0.1,0.15,0.05,0.1);
       			Set2dHist(Prob[iun][id][ij][ik][ipt],lplot_xtitle, lplot_xtitle,"Probability",titoff1, titsize1);
       			Prob[iun][id][ij][ik][ipt]->Draw("colz"); leg1->Draw();
       			CMS_lumi( cpt6, iPeriod, iPos ); cpt6->Update();
       			sprintf(pdfname, "TUnfold_prob_%i.pdf(" ,iun); sprintf(pdfname1, "TUnfold_prob_%i.pdf" ,iun);sprintf(pdfname2, "TUnfold_prob_%i.pdf)" ,iun); //Tunfolded_typ_0_pt0_eta0_3
          		if(id==0 && ij==0 && ik==0 && ipt ==0){cpt6->Print(pdfname,"pdf");
          		}else if(id==2 && ij==1 && ik==9 && ipt==9) {cpt6->Print(pdfname2,"pdf");
          		}else{cpt6->Print(pdfname,"pdf");};

       			cpt7->cd();
       			SetMycanvas(cpt7,0,0.1,0.15,0.05,0.1);
       			Set2dHist(Ematrix[iun][id][ij][ik][ipt],lplot_xtitle, lplot_xtitle,"",titoff1, titsize1);
       			Ematrix[iun][id][ij][ik][ipt]->Draw("colz"); leg1->Draw();
       			CMS_lumi( cpt7, iPeriod, iPos ); cpt7->Update();
       			sprintf(pdfname, "TUnfold_COV_%i.pdf(" ,iun); sprintf(pdfname1, "TUnfold_COV_%i.pdf" ,iun);sprintf(pdfname2, "TUnfold_COV_%i.pdf)" ,iun); //Tunfolded_typ_0_pt0_eta0_3
          		if(id==0 && ij==0 && ik==0 && ipt ==0){cpt7->Print(pdfname,"pdf");
          		}else if(id==2 && ij==1 && ik==9 && ipt==9) {cpt7->Print(pdfname2,"pdf");
          		}else{cpt7->Print(pdfname,"pdf");};
				}
			}
      		}//end of phase space cut and variable loop
	}
}//End of Unfolded plot
//----------------------------------------------------
int tmc = 0; //0 for Py8, 1 for MG , 2 for Herwig
//All Unfold in one plot
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
                        else if(ipt==9){sprintf(Title,"%s_{%s}^{%s}:      <P_{T}> %i %s", obs_def[id], jet_num[ij], k_fact[ik], HT2range[ipt] ,"GeV/c");}
                        sprintf(Yaxis," %s %s^{%s}" ,obs_logy[ik], obs_def[id], k_fact[ik]);
          		MyHist->SetTitle(Title);
          		MyHist->GetXaxis()->SetTitle("");
          		MyHist->GetYaxis()->SetTitle(Yaxis);

          		TH1D *unfold_input[unfold_ty];
          		const char *MCinput_index[unfold_ty], *data_index[3];
          		for(int iout = 0 ; iout < unfold_ty ; iout++){
            			unfold_input[iout] = (TH1D*) Unfold[iout][id][ij][ik][ipt]->Clone();
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
          	sprintf(lplot_xtitle," %s %s^{%s} %d" ,jet_num[ij],obs_def[id],k_fact[ik],HT2range[ipt]);
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

          	sprintf(pdfname, "unfold_plot_%i.pdf(", 1234); sprintf(pdfname1, "unfold_plot_%i.pdf" ,1234);sprintf(pdfname2, "unfold_plot_%i.pdf)",1234); //TRefolded_typ_0_pt0_eta0_3
          	if(id==0 && ij==0 && ik==0 && ipt ==0){cpt0->Print(pdfname,"pdf");
          	}else if(id==2 && ij==1 && ik==9 && ipt==9) {cpt0->Print(pdfname2,"pdf");
          	}else{cpt0->Print(pdfname,"pdf");};
      			}//end of phase space cut and variable loop
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
	  		TH1D *MyHist  = (TH1D*) Refold[iun][id][ij][ik][ipt]->Clone();
	  		Integralhist(MyHist);
	  		//divBybinWidth(MyHist);
	  		Myplotset(MyHist,0,0);
	  
	  		if(ipt<=8){sprintf(Title,"%s_{%s}^{%s}:     %i <P_{T}< %i %s", obs_def[id], jet_num[ij], k_fact[ik], HT2range[ipt] , HT2range[ipt+1] ,"GeV/c" );}
                        else if(ipt==9){sprintf(Title,"%s_{%s}^{%s}:      <P_{T}> %i %s", obs_def[id], jet_num[ij], k_fact[ik], HT2range[ipt] ,"GeV/c");}
                        sprintf(Yaxis," %s %s^{%s}" ,obs_logy[ik], obs_def[id], k_fact[ik]);
	  		MyHist->SetTitle(Title);
	  		MyHist->GetXaxis()->SetTitle("");
	  		MyHist->GetYaxis()->SetTitle(Yaxis);
	 
	  		TH1D *MC_input[nmc];
	  		const char *MCinput_index[nmc], *data_index[1];
	  		for(int iout = 0 ; iout < nmc ; iout++){
	    			MC_input[iout] = (TH1D*) MC_reco[iout][id][ij][ik][ipt]->Clone();
	   			//MC_input[iout]->Rebin(2);
	    			Integralhist(MC_input[iout]);
	   			//divBybinWidth(MC_input[iout]);
	    			MCinput_index[iout]= mcnamerco[iout]; }
	  
	  			data_index[0]= RefoldEra[iun];
	  			char lplot_xtitle[100];
	  			sprintf(lplot_xtitle," %s %s^{%s} %d" ,jet_num[ij],obs_def[id],k_fact[ik],HT2range[ipt]);
	  			//float ratio_range1[2]={1.2,0.9};
	  			int num1[2]={nmc,1};
	  			float lpos1[7] ={.32,0.2,0.55,0.38, .04, 1.5,0.7};
	  
	  			cpt0->cd();
	  			cpt0->SetBorderSize(0);
	  			cpt0->SetRightMargin(0.0);
	  			cpt0->SetTopMargin(0.0);
	  			cpt0 =(TCanvas*)(ratio_can(num1, lpos1, MyHist, MC_input, lplot_xtitle,MCinput_index,data_index));
	  			CMS_lumi( cpt0, iPeriod, iPos ); cpt0->Update();
	  
	  			sprintf(pdfname, "Refold_plot_%i.pdf(" ,iun); sprintf(pdfname1, "Refold_plot_%i.pdf" ,iun);sprintf(pdfname2, "Refold_plot_%i.pdf)" ,iun); //TRefolded_typ_0_pt0_eta0_3
	  			if(id==0 && ij==0 && ik==0 && ipt ==0){cpt0->Print(pdfname,"pdf");
	  			}else if(id==2 && ij==1 && ik==9 && ipt==9) {cpt0->Print(pdfname2,"pdf");
	  			}else{cpt0->Print(pdfname,"pdf");};
	  			}
			}
      		}  //end of phase space cut and variable loop
    	}
}//End of Refolded plot  
//----------------------------------------------------
//Refold in one plot
    	for(int id=0; id <ndef; id++){
		for(int ij=0; ij<njet; ij++ ){
    	  		for(int ik =0 ; ik < nkappa ; ik++){
        			for(int ipt = 0 ; ipt <nHLTmx; ipt++){
          		TH1D *MyHist  = (TH1D*) MC_reco[tmc][id][ij][ik][ipt]->Clone();
          		Integralhist(MyHist);
         		//divBybinWidth(MyHist);
          		Myplotset(MyHist,0,0);

          		if(ipt<=8){sprintf(Title,"%s_{%s}^{%s}:     %i <P_{T}< %i %s", obs_def[id], jet_num[ij], k_fact[ik], HT2range[ipt] , HT2range[ipt+1] ,"GeV/c" );}
                        else if(ipt==9){sprintf(Title,"%s_{%s}^{%s}:      <P_{T}> %i %s", obs_def[id], jet_num[ij], k_fact[ik], HT2range[ipt] ,"GeV/c");}
                        sprintf(Yaxis," %s %s^{%s}" ,obs_logy[ik], obs_def[id], k_fact[ik]);
          		MyHist->SetTitle(Title);
          		MyHist->GetXaxis()->SetTitle("");
          		MyHist->GetYaxis()->SetTitle(Yaxis);

          		TH1D *Refold_input[unfold_ty];
          		const char *MCinput_index[unfold_ty], *data_index[3];
          		for(int iout = 0 ; iout < unfold_ty ; iout++){
            			Refold_input[iout] = (TH1D*) Refold[iout][id][ij][ik][ipt]->Clone();
           			//MC_input[iout]->Rebin(2);
            			Integralhist(Refold_input[iout]);
           			//divBybinWidth(Refold_input[iout]);
            			MCinput_index[iout]= RefoldEra[iout]; }

          			data_index[0]= mcnamerco[tmc];
          			data_index[1]= "MC";
          			data_index[2]= "Refolded";
          			char lplot_xtitle[100];
          			sprintf(lplot_xtitle," %s %s^{%s} %d" ,jet_num[ij],obs_def[id],k_fact[ik],HT2range[ipt]);
          			//float ratio_range1[2]={1.2,0.9};
          			int num1[3]={unfold_ty,1,0};
          			float lpos1[7] ={.32,0.2,0.55,0.38, .04, 1.5,0.7};

          			cpt0->cd();
          			cpt0->SetBorderSize(0);
          			cpt0->SetRightMargin(0.0);
          			cpt0->SetTopMargin(0.0);
          			cpt0 =(TCanvas*)(ratio_can1(num1, lpos1, MyHist, Refold_input, lplot_xtitle,MCinput_index,data_index));
          			CMS_lumi( cpt0, iPeriod, iPos ); cpt0->Update();

          			sprintf(pdfname, "Refold_plot_%i.pdf(", 1234); sprintf(pdfname1, "Refold_plot_%i.pdf" ,1234);sprintf(pdfname2, "Refold_plot_%i.pdf)",1234); //TRefolded_typ_0_pt0_eta0_3
          			if(id==0 && ij==0 && ik==0 && ipt ==0){cpt0->Print(pdfname,"pdf");
          			}else if(id==2 && ij==1 && ik==9 && ipt==9) {cpt0->Print(pdfname2,"pdf");
          			}else{cpt0->Print(pdfname,"pdf");};
			}
      		}//end of phase space cut and variable loop
    	}
} 
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

  	if(Nplot[1]==1){gPad->SetLogy();}  //condition for log scale
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
    /* 			}else{
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
