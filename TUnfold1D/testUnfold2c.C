#include <iostream>
#include <fstream>
#include <cmath>
#include <array>
#include <vector>
#include <TArrayC.h>
#include <string>
#include "TH1.h"
#include "TH2D.h"
#include <TH1D.h>
#include "TH1F.h"
#include "TH2F.h"
#include "TVirtualPad.h"
#include "TApplication.h"
#include "TBranch.h"
#include "TBranchElement.h"
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include <TMath.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TFitter.h>
#include <TF1.h>
#include <TStyle.h>
#include <TVector.h>
#include <TGraph.h>
#include <map>
#include <cmath>
#include <TFile.h>
#include <TDirectoryFile.h>
#include <TTree.h>
#include <TProfile.h>
#include <TH2D.h>
#include <THStack.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLine.h>
#include <TLorentzVector.h>
#include "TUnfoldDensity.h"
#include "TUnfoldIterativeEM.h"

#include "TUnfoldBinning.h"
#include "TUnfoldBinningXML.h"
#include "TUnfoldDensity.h"
#include "TUnfold.h"
#include "TUnfoldIterativeEM.h"
#include "TUnfoldSys.h"
#include <TVectorD.h>
#include <TDecompSVD.h>

#define CLOUSER

using namespace std;

static const auto feps = numeric_limits<float>::epsilon();

void testUnfold2c()
{
  //switch on histogram errors
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  
  //Input Data and MC histogram
  TFile *inputData=new TFile("Data_UL2017.root");
  TFile *inputMC=new TFile("PY8_bin.root");
  
  TFile *inputMC1=new TFile("MG5_PY8_bin.root");
  TFile *inputMC2=new TFile("HW_CH3_QCD_Pt-15to7000_TuneCH3_Flat_13TeV_herwig7.root");
  
  //Unfolded Data and Covarince matrix, efficincy,fake rate, purity, stability
  TFile *outputFile=new TFile("testunfold2c_unfolded.root","recreate");
  
  //int irbin = 2;
  int irbin =1;
  
  char histname[100], name[100];
  const int nHLTmx=10; //HT2 Range
  const int njetptmn = nHLTmx;
  double leadingPtThreshold[njetptmn+1] = {92, 119, 185, 251, 319, 388, 467, 518, 579, 669, 3000}; //Fit Value dijet trigger

  static const int ndef = 3;
  static const int njet = 2;
  static const int nkappa = 10;
  const char* obs_def[ndef]={"Q","Q_{L}","Q_{T}"};
  const char* jet_num[njet]={"Leading-Jet","Sub-Leading-Jet"};
  const char* k_fact[nkappa]={"k=0.1","k=0.2","k=0.3","k=0.4","k=0.5","k=0.6","k=0.7","k=0.8","k=0.9","k=1.0"};

  //----------------------for reco and gen bin number
  int rnbinsx[ndef][njet][nkappa][nHLTmx];
  int gnbinsx[ndef][njet][nkappa][nHLTmx];
  
  TH1D *MC_Reco[ndef][njet][nkappa][njetptmn];  //Reconstructed MC
  TH1D *MC_fake[ndef][njet][nkappa][njetptmn];  //Fake :  Reco but No Gen
  TH1D *MC_Gen[ndef][njet][nkappa][njetptmn];   //Generator MC
  TH1D *MC_miss[ndef][njet][nkappa][njetptmn];   //Miss:  No Reco but in Gen
  TH1D *Data_Reco[ndef][njet][nkappa][njetptmn];    //Reconstructed Data
  TH1D *PsudoData_Gen[ndef][njet][nkappa][njetptmn];    //Gen Level Psudo Data i.e MC is treated as Data
  TH2D *h2dGenDetMC[ndef][njet][nkappa][njetptmn];   // MC Gen Vs Reco
 
  TH1D *MC_Reco1[ndef][njet][nkappa][njetptmn];  //Reconstructed MC copy
  TH1D *MC_fake1[ndef][njet][nkappa][njetptmn];  //Fake :  Reco but No Gen Copy
  TH1D *MC_Gen1[ndef][njet][nkappa][njetptmn];   //Generator level MC copy 
  TH1D *MC_miss1[ndef][njet][nkappa][njetptmn];   ///Miss:  No Reco but in Gen Copy
  TH1D *Data_Reco1[ndef][njet][nkappa][njetptmn];    //Reconstructed Data copy
  TH2D *h2dGenDetMC1[ndef][njet][nkappa][njetptmn];   // MC Gen Vs Reco copy
 
  TH1D *MC_Reco2[ndef][njet][nkappa][njetptmn];  //Reconstructed MC copy
  TH1D *MC_fake2[ndef][njet][nkappa][njetptmn];  //Fake :  Reco but No Gen Copy
  TH1D *MC_Gen2[ndef][njet][nkappa][njetptmn];   //Generator level MC copy
  TH1D *MC_miss2[ndef][njet][nkappa][njetptmn];   //Miss:  No Reco but in Gen Copy
  TH1D *Data_Reco2[ndef][njet][nkappa][njetptmn];    //Reconstructed Data copy
  TH2D *h2dGenDetMC2[ndef][njet][nkappa][njetptmn];  // MC Gen Vs Reco copy
  
  //---------------------MG5
  TH1D *MG5_MC_Reco[ndef][njet][nkappa][njetptmn];  //Reconstructed MC
  TH1D *MG5_MC_Gen[ndef][njet][nkappa][njetptmn];   //Generator level MC
  TH2D *MG5_h2dGenDetMC[ndef][njet][nkappa][njetptmn];   // MC Gen Vs Reco
 
  TH1D *MG5_MC_Reco1[ndef][njet][nkappa][njetptmn];  //Reconstructed MC copy
  TH1D *MG5_MC_Gen1[ndef][njet][nkappa][njetptmn];   //Generator level MC copy
  TH2D *MG5_h2dGenDetMC1[ndef][njet][nkappa][njetptmn];   // MC Gen Vs Reco copy
 
  TH1D *MG5_MC_Reco2[ndef][njet][nkappa][njetptmn];  //Reconstructed MC copy
  TH1D *MG5_MC_Gen2[ndef][njet][nkappa][njetptmn];   //Generator level MC copy
  TH2D *MG5_h2dGenDetMC2[ndef][njet][nkappa][njetptmn]; // MC Gen Vs Reco copy
 
  //----------------------HW
  TH1D *HW7_MC_Reco[ndef][njet][nkappa][njetptmn];  //Reconstructed MC
  TH1D *HW7_MC_Gen[ndef][njet][nkappa][njetptmn];   //Generator level MC
  TH2D *HW7_h2dGenDetMC[ndef][njet][nkappa][njetptmn];   // MC Gen Vs Reco
 
  TH1D *HW7_MC_Reco1[ndef][njet][nkappa][njetptmn];  //Reconstructed MC copy
  TH1D *HW7_MC_Gen1[ndef][njet][nkappa][njetptmn];   //Generator level MC copy
  TH2D *HW7_h2dGenDetMC1[ndef][njet][nkappa][njetptmn];   // MC Gen Vs Reco copy
 
  TH1D *HW7_MC_Reco2[ndef][njet][nkappa][njetptmn];  //Reconstructed MC copy
  TH1D *HW7_MC_Gen2[ndef][njet][nkappa][njetptmn];   //Generator level MC copy
  TH2D *HW7_h2dGenDetMC2[ndef][njet][nkappa][njetptmn]; // MC Gen Vs Reco copy
 
  TH1D *unfold_NoReg[ndef][njet][nkappa][njetptmn];
  TH1D *unfold_Tau[ndef][njet][nkappa][njetptmn];
  TH1D *unfold_L[ndef][njet][nkappa][njetptmn];
 
  TH1D* hist_eff[ndef][njet][nkappa][njetptmn];
  TH1D* hist_fake[ndef][njet][nkappa][njetptmn];
  TH1D* hist_purity[ndef][njet][nkappa][njetptmn];
  TH1D* hist_stbl[ndef][njet][nkappa][njetptmn];
 
  TH1D* hist_eff1[ndef][njet][nkappa][njetptmn];
  TH1D* hist_fake1[ndef][njet][nkappa][njetptmn];
  TH1D* hist_purity1[ndef][njet][nkappa][njetptmn];
  TH1D* hist_stbl1[ndef][njet][nkappa][njetptmn];
 
  TH2D* COV_Mat_NoReg[ndef][njet][nkappa][njetptmn];
  TH2D* COV_Mat_Tau[ndef][njet][nkappa][njetptmn];
  TH2D* COV_Mat_L[ndef][njet][nkappa][njetptmn];
 
  TH2D* corr_mat_NoReg[ndef][njet][nkappa][njetptmn];
 
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetPadRightMargin(0.06);
  gStyle->SetPadBottomMargin(0.18);
  gStyle->SetPadLeftMargin(0.19);
  gStyle->SetPadTopMargin(0.07);
  int font=42;
  gStyle->SetLabelFont(font,"XYZ");
  gStyle->SetLegendFont(font);
  gStyle->SetStatFont(font);
  gStyle->SetTextFont(font);
  gStyle->SetTitleFont(font,"XYZ");
  gStyle->SetTitleOffset(1.5,"y");
  gStyle->SetTitleOffset(1.2,"x");
  gStyle->SetTitleSize(0.8,"P");
  gStyle->SetTitleSize(0.06,"xy");
  gStyle->SetLabelSize(0.05,"xy");
  gStyle->SetLabelOffset(0.012,"xy");
 
  TDirectoryFile *inputDir0=new TDirectoryFile("Data","Inputs Data");
  TDirectoryFile *inputDir1=new TDirectoryFile("Pythia8"," Pythia8 , MC and Probability Matrix");
  TDirectoryFile *inputDir2=new TDirectoryFile("MG8","Madgraph, MC and Probability Matrix");
  TDirectoryFile *inputDir3=new TDirectoryFile("HW7","Herwig7 MC and Probability Matrix");
  TDirectoryFile *foldpy8=new TDirectoryFile("fold","folded with Probablility matrix");
  TDirectoryFile *Unfold=new TDirectoryFile("Unfold","Unfolded, Refold, correlation");
 
  int subtract_background(TH2D* h2d_correl, TH1D* reco, TH1D* gen, TH1D* data, double* fakerate, double* effi, double* purity, double* stbl);
  void Fold(TH2D* HistoMatrix, TH1D* HistReco, TH1D* HistoGen, TH1D* miss, TH1D* fake, TH1D* HistoCorrect);
  //void Condition (TH2* RM, TH1* miss);
  double Condition (TH2 * RM);
  double ConditionV2(TH2* RM); 

//Read Input Data MC and Response matrix
for(int id=0; id<ndef; id++){
 	for (int ij=0; ij<njet; ij++){
		for (int ik=0; ik<nkappa; ik++){
			for(int ipt=0; ipt<njetptmn; ipt++){
       		//Reco Data
       	 	inputDir0->cd();
       		sprintf(histname, "analyzeBasicPat/reco_jc_d%i_j%i_k%i_pt%i_eta0", id, ij, ik, ipt); //reco_jc_d0_j0_k0_pt0_eta0
      		Data_Reco1[id][ij][ik][ipt] = (TH1D*) inputData->Get(histname);

       		for (int ibin =1 ; ibin <  Data_Reco1[id][ij][ik][ipt]->GetNbinsX()+1; ibin++ ){ 
       			if(Data_Reco1[id][ij][ik][ipt]->GetBinContent(ibin) == 0) {cout << "Data Reco Bin is Zero for bin number : ***** "<<  ibin  << endl; }
       			}      
#ifdef CLOUSER
       		sprintf(histname, "analyzeBasicPat/gen_jc_d%i_j%i_k%i_pt%i_eta0", id, ij, ik, ipt); //gen_jc_d0_j0_k0_pt0_eta0  ?? CHECK
       		PsudoData_Gen[id][ij][ik][ipt] = (TH1D*) inputData->Get(histname);
		cout << "Closure input histogram: "<<histname<<endl;
		cout <<"Entries:"<<PsudoData_Gen[id][ij][ik][ipt]->GetNbinsX()<<endl;
#else
       		PsudoData_Gen[id][ij][ik][ipt] = (TH1D*) inputMC->Get(histname);
		//cout << "Closure extra histogram: "<<histname<<endl;
#endif
       		inputDir1->cd();
       		sprintf(histname, "analyzeBasicPat/reco_jc_d%i_j%i_k%i_pt%i_eta0", id, ij, ik, ipt); //reco_jc_d0_j0_k0_pt0_eta0
       		MC_Reco1[id][ij][ik][ipt] = (TH1D*) inputMC->Get(histname);
		//cout << histname <<  "  Reco =" << MC_Reco1[id][ij][ik][ipt]->GetEntries() <<endl;
       
       		int recobins = MC_Reco1[id][ij][ik][ipt]->GetNbinsX();
       		rnbinsx[id][ij][ik][ipt]=MC_Reco1[id][ij][ik][ipt]->GetNbinsX();
       		if(recobins !=Data_Reco1[id][ij][ik][ipt]->GetNbinsX()) {cout << "reco Bin miss Match, Check bins"<<endl;}
       
       		//MC Fake
       		sprintf(histname, "analyzeBasicPat/recofake_jc_d%i_j%i_k%i_pt%i_eta0", id, ij, ik, ipt); //recofake_jc_d0_j0_k0_pt0_eta0
       		MC_fake1[id][ij][ik][ipt] = (TH1D*) inputMC->Get(histname);
       
       		cout << "Fake = " <<MC_fake1[id][ij][ik][ipt]->GetEntries() <<" Reco-fake: " <<(MC_Reco1[id][ij][ik][ipt]->GetEntries() - MC_fake1[id][ij][ik][ipt]->GetEntries())<<endl;
       
       		//Gen MC
      		sprintf(histname, "analyzeBasicPat/gen_jc_d%i_j%i_k%i_pt%i_eta0", id, ij, ik, ipt); 
       		MC_Gen1[id][ij][ik][ipt] = (TH1D*) inputMC->Get(histname);
   
       		int genbins = MC_Gen1[id][ij][ik][ipt]->GetNbinsX();
      	 	gnbinsx[id][ij][ik][ipt]=MC_Gen1[id][ij][ik][ipt]->GetNbinsX();
       		for(int ibin = 1; ibin < MC_Gen1[id][ij][ik][ipt]->GetNbinsX()+1; ibin++ ){
       			if (MC_Gen1[id][ij][ik][ipt]->GetBinContent(ibin) == 0) {cout << "MC gen Bin is Zero for bin number :***** "<<  ibin  << endl; }
       			}
       		//cout <<" Gen; " << MC_Gen1[id][ij][ik][ipt]->GetEntries()<<endl;

       		//MC miss 
       		sprintf(histname, "analyzeBasicPat/genmiss_jc_d%i_j%i_k%i_pt%i_eta0", id, ij, ik, ipt); //genmiss_jc_d0_j0_k0_pt0_eta0
       		MC_miss1[id][ij][ik][ipt] = (TH1D*) inputMC->Get(histname);
      
       		//int GenMiss= ((MC_Gen1[id][ij][ik][ipt]->GetEntries()) - (MC_miss1[id][ij][ik][ipt]->GetEntries()));
       		//cout << " Miss= " << MC_miss1[id][ij][ik][ipt]->GetEntries() <<" Gen-Miss: " << GenMiss<<endl;
       		cout << "Miss = " << MC_miss1[id][ij][ik][ipt]->GetEntries() <<" Gen-Miss: " << (MC_Gen1[id][ij][ik][ipt]->GetEntries() - MC_miss1[id][ij][ik][ipt]->GetEntries())<<endl;
       
       		//Response Matrix
       		sprintf(histname, "analyzeBasicPat/RM_jc_d%i_j%i_k%i_pt%i_eta0", id, ij, ik, ipt); //RM_jc_d0_j0_k0_pt0_eta0    
       		h2dGenDetMC1[id][ij][ik][ipt] = (TH2D*) inputMC->Get(histname);   //Xgen(coarse) , Yreco(fine) IS IT FINE ??
       		cout << "Corr = "  << h2dGenDetMC1[id][ij][ik][ipt]->GetEntries() <<endl;
       		cout << "Corr NbinsX = "<<h2dGenDetMC1[id][ij][ik][ipt]->GetNbinsX() <<endl;
       		cout << "Corr NbinsY = "<<h2dGenDetMC1[id][ij][ik][ipt]->GetNbinsY() <<endl;
     
	        // Clone histos	
       		Data_Reco2[id][ij][ik][ipt] = (TH1D*)Data_Reco1[id][ij][ik][ipt]->Clone();
       		MC_Reco2[id][ij][ik][ipt] = (TH1D*)MC_Reco1[id][ij][ik][ipt]->Clone();
       		MC_Gen2[id][ij][ik][ipt] = (TH1D*)MC_Gen1[id][ij][ik][ipt]->Clone();
       		h2dGenDetMC2[id][ij][ik][ipt] = (TH2D*)h2dGenDetMC1[id][ij][ik][ipt]->Clone();
       
       		MC_miss2[id][ij][ik][ipt] = (TH1D*)MC_miss1[id][ij][ik][ipt]->Clone();
       		MC_fake2[id][ij][ik][ipt] = (TH1D*)MC_fake1[id][ij][ik][ipt]->Clone();
      
//----------------------Check RM Projection with Reco(gen)-Fake(miss)  : Patrick 1 Sep20		
		//Saved in Pythia8 directory
       		TH1* RMx = h2dGenDetMC1[id][ij][ik][ipt]->ProjectionX();
       		TH1* RMy = h2dGenDetMC1[id][ij][ik][ipt]->ProjectionY();

       		sprintf(name,"ProjectX_d%i_j%i_k%i_pt%i_eta0", id, ij, ik, ipt);
       		RMx->SetNameTitle(name,name);
       		RMx->Write();

       		sprintf(name,"Recominusfake_d%i_j%i_k%i_pt%i_eta0", id, ij, ik, ipt);
       		TH1* RecoFakeCorrect = (TH1D*)MC_Reco1[id][ij][ik][ipt]->Clone();  // IS IT CORRECTED ??
       		RecoFakeCorrect->Reset();
       		RecoFakeCorrect->SetNameTitle(name,name);
       
       		sprintf(name,"ProjectY_d%i_j%i_k%i_pt%i_eta0", id, ij, ik, ipt);
       		RMy->SetNameTitle(name,name);
       		RMy->Write();
       
       		sprintf(name,"Genminusmiss_d%i_j%i_k%i_pt%i_eta0", id, ij, ik, ipt);  // IS IT CORRECTED ??
       		TH1* GenMissCorrect = (TH1D*)MC_Gen1[id][ij][ik][ipt]->Clone();
       		GenMissCorrect->Reset();
       		GenMissCorrect->SetNameTitle(name,name);

        	for (int i = 1; i <= RecoFakeCorrect->GetNbinsX(); ++i) {
         		double content = MC_Reco1[id][ij][ik][ipt]->GetBinContent(i);
         		double factor = MC_fake1[id][ij][ik][ipt]->GetBinContent(i);
         		content -= factor;
			//cout << " fake Factor " << factor <<endl;
         		RecoFakeCorrect->SetBinContent(i, content);
			//cout <<"RecoFakeCorrect: "<<RecoFakeCorrect->GetBinContent(i)<<endl;
       			}
         
		//RecoFakeCorrect->Divide(RMx);
        	//RecoFakeCorrect->SetMinimum(0.5); MC_fake1[ity][ivar][ipt]->SetMaximum(1.6);
		//cout <<"RecoFakeCorrect: "<<RecoFakeCorrect->GetEntries()<<endl;
		//cout <<"RecoFakeCorrect: "<<RecoFakeCorrect->GetBinContent()<<endl;
        	RecoFakeCorrect->Write();

        	for (int i = 1; i <= GenMissCorrect->GetNbinsX(); ++i) {
         		double content = MC_Gen1[id][ij][ik][ipt]->GetBinContent(i);
         		double factor = MC_miss1[id][ij][ik][ipt]->GetBinContent(i);
         		content -= factor;
			//cout << " Miss Factor " << factor << endl;
         		GenMissCorrect->SetBinContent(i, content);
			//cout <<"GenMissCorrect: "<<GenMissCorrect->GetBinContent(i)<<endl;
       			}
         
		//GenMissCorrect->Divide(RMy);
        	//GenMissCorrect->SetMinimum(0.85); MC_fake1[ity][ivar][ipt]->SetMaximum(1.15);
		//cout <<"GenMissCorrect: "<<GenMissCorrect->GetEntries()<<endl;
		//cout <<"GenMissCorrect: "<<GenMissCorrect->GetBinContent()<<endl;
	 	GenMissCorrect->Write();

//--------------------------------Fake and Miss rate : Thanks to Patrick  14Aug20
//
       		MC_fake1[id][ij][ik][ipt]->Divide(MC_fake1[id][ij][ik][ipt], MC_Reco1[id][ij][ik][ipt], 1, 1, "b");
       		//MC_fake1[ity][ivar][ipt]->Divide(MC_Reco1[ity][ivar][ipt]);
      		MC_miss1[id][ij][ik][ipt]->Divide(MC_miss1[id][ij][ik][ipt], MC_Gen1[id][ij][ik][ipt], 1, 1, "b");
       		//MC_miss1[ity][ivar][ipt]->Divide(MC_Gen1[ity][ivar][ipt]);
       
       		MC_fake1[id][ij][ik][ipt]->SetMinimum(-0.05); MC_fake1[id][ij][ik][ipt]->SetMaximum(1.01);
       		MC_miss1[id][ij][ik][ipt]->SetMinimum(-0.05); MC_miss1[id][ij][ik][ipt]->SetMaximum(1.01);

       		Data_Reco[id][ij][ik][ipt] =(TH1D*)Data_Reco1[id][ij][ik][ipt]->Clone();
       		MC_Reco[id][ij][ik][ipt] = (TH1D*)MC_Reco1[id][ij][ik][ipt]->Clone();
       		MC_Gen[id][ij][ik][ipt] = (TH1D*)MC_Gen1[id][ij][ik][ipt]->Clone();
       		MC_fake[id][ij][ik][ipt] = (TH1D*)MC_fake1[id][ij][ik][ipt]->Clone();
       		MC_miss[id][ij][ik][ipt] = (TH1D*)MC_miss1[id][ij][ik][ipt]->Clone();
       		h2dGenDetMC[id][ij][ik][ipt] = (TH2D*)h2dGenDetMC1[id][ij][ik][ipt]->Clone();
       
       		inputDir0->cd();
       		Data_Reco[id][ij][ik][ipt]->Write();
       		PsudoData_Gen[id][ij][ik][ipt]->Write();
       		inputDir1->cd();
       		MC_Reco[id][ij][ik][ipt]->Write();
       		MC_Gen[id][ij][ik][ipt]->Write();
       		MC_fake[id][ij][ik][ipt]->Write();
       		MC_miss[id][ij][ik][ipt]->Write();
       		h2dGenDetMC[id][ij][ik][ipt]->Write();      
			} 
     		}
  	}
}
cout << "PY8 Histos OK" <<endl;

//------------------Read Madgraph
cout << "Entering MG+PY8 Loop ..."<<endl;
for(int id=0; id <ndef; id++){
	for (int ij=0; ij<njet; ij++){
        	for (int ik=0; ik<nkappa; ik++){
     			for(int ipt = 0 ; ipt < njetptmn ; ipt++){
     
       		inputDir2->cd();
       		sprintf(histname, "analyzeBasicPat/reco_jc_d%i_j%i_k%i_pt%i_eta0", id, ij, ik, ipt); //reco_typ_1_pt4_eta0_24 reco_jc_D3_j1_k5_pt1_eta0
       		MG5_MC_Reco1[id][ij][ik][ipt] = (TH1D*) inputMC1->Get(histname);
       		//cout << histname << endl;
       		MG5_MC_Reco1[id][ij][ik][ipt]->Write();
       
       		//Gen MC
       		sprintf(histname, "analyzeBasicPat/gen_jc_d%i_j%i_k%i_pt%i_eta0", id, ij, ik, ipt); //reco_typ_1_pt4_eta0_24
       		MG5_MC_Gen1[id][ij][ik][ipt] = (TH1D*) inputMC1->Get(histname);
       		//MC_Gen1[ity][ivar][ipt]->Rebin(2);
       		MG5_MC_Gen1[id][ij][ik][ipt]->Write();

       		//Response Matrix
       		sprintf(histname, "analyzeBasicPat/RM_jc_d%i_j%i_k%i_pt%i_eta0", id, ij, ik, ipt); //corr_typ_0_pt2_eta0_3
       		MG5_h2dGenDetMC1[id][ij][ik][ipt] = (TH2D*) inputMC1->Get(histname);   //Xgen(coarse) , Yreco(fine)
       		//MG5_h2dGenDetMC1[ity][ivar][ipt]->RebinY(2);
       		MG5_h2dGenDetMC1[id][ij][ik][ipt]->Write();   //Xgen(coarse) , Yreco(fine)
	 
      	 	MG5_MC_Reco2[id][ij][ik][ipt] = (TH1D*)MG5_MC_Reco1[id][ij][ik][ipt]->Clone();
       		MG5_MC_Gen2[id][ij][ik][ipt] = (TH1D*)MG5_MC_Gen1[id][ij][ik][ipt]->Clone();
       		MG5_h2dGenDetMC2[id][ij][ik][ipt] = (TH2D*)MG5_h2dGenDetMC1[id][ij][ik][ipt]->Clone();
			}
     		}
   	}
 }
cout << "MG+PY8 Histos OK" <<endl;
 
//----------------HW7
cout<<"Entering HW7 Loop ..."<<endl;
for(int id=0; id <ndef; id++){
	for (int ij=0; ij<njet; ij++){
       		for (int ik=0; ik<nkappa; ik++){
     			for(int ipt = 0 ; ipt < njetptmn ; ipt++){
       
         	inputDir3->cd();
         	sprintf(histname, "analyzeBasicPat/reco_jc_d%i_j%i_k%i_pt%i_eta0", id, ij, ik, ipt); //reco_typ_1_pt4_eta0_24
         	HW7_MC_Reco1[id][ij][ik][ipt] = (TH1D*) inputMC2->Get(histname);
	 	//cout << histname << endl;
         	HW7_MC_Reco1[id][ij][ik][ipt]->Write();
       
		//Gen MC
         	sprintf(histname, "analyzeBasicPat/gen_jc_d%i_j%i_k%i_pt%i_eta0", id, ij, ik, ipt); //reco_typ_1_pt4_eta0_24
         	HW7_MC_Gen1[id][ij][ik][ipt] = (TH1D*) inputMC2->Get(histname);
	 	//HW7_MC_Gen1[ity][ivar][ipt]->Rebin(2);
         	HW7_MC_Gen1[id][ij][ik][ipt]->Write();
	 
         	//Response Matrix
         	sprintf(histname, "analyzeBasicPat/RM_jc_d%i_j%i_k%i_pt%i_eta0", id, ij, ik, ipt); //corr_typ_0_pt2_eta0_3
         	HW7_h2dGenDetMC1[id][ij][ik][ipt] = (TH2D*) inputMC2->Get(histname);   //Xgen(coarse) , Yreco(fine)
	 	//HW7_h2dGenDetMC1[ity][ivar][ipt]->RebinY(2);
	 	HW7_h2dGenDetMC1[id][ij][ik][ipt]->Write();   //Xgen(coarse) , Yreco(fine)

         	HW7_MC_Reco2[id][ij][ik][ipt] = (TH1D*)HW7_MC_Reco1[id][ij][ik][ipt]->Clone();
         	HW7_MC_Gen2[id][ij][ik][ipt] = (TH1D*)HW7_MC_Gen1[id][ij][ik][ipt]->Clone();
         	HW7_h2dGenDetMC2[id][ij][ik][ipt] = (TH2D*)HW7_h2dGenDetMC1[id][ij][ik][ipt]->Clone();	 
			}
     		}
   	}
}
cout << "HW7 Histos OK" <<endl;

//------------------Fold check : Patrick 1 Sep 20
//Get Probability Matrix  
//(Multiply the gen level by the probability matrix. Of course, don't forget to account for miss and fake entries (if applicable).
foldpy8->cd();
for(int id=0; id <ndef; id++){
     	for (int ij=0; ij<njet; ij++){
                for (int ik=0; ik<nkappa; ik++){
     			for(int ipt = 0; ipt < njetptmn ; ipt++){

      		TH2D* RM  = (TH2D*)h2dGenDetMC1[id][ij][ik][ipt]->Clone();
      		TH1D* Reco = (TH1D*)MC_Reco1[id][ij][ik][ipt]->Clone();
      		TH1D* Gen  = (TH1D*)MC_Gen1[id][ij][ik][ipt]->Clone();
      		TH1D* fake = (TH1D*)MC_fake2[id][ij][ik][ipt]->Clone();
      		TH1D* miss = (TH1D*)MC_miss2[id][ij][ik][ipt]->Clone();
      
      		//RM->RebinY(2);Gen->Rebin(2);miss->Rebin(2);  // Have to check rebin part
		//RM->RebinX(2);Gen->Rebin(2);miss->Rebin(2);
      		TH1D* Folded = (TH1D*)MC_Reco1[id][ij][ik][ipt]->Clone(); Folded->Reset();
       
      		Fold(RM, Reco, Gen, miss, fake, Folded);
      		sprintf(name,"Fold_d%i_j%i_k%i_pt%i_eta0", id, ij, ik, ipt);
      		Folded->SetNameTitle(name,name);

      		Folded->Write();
			}
     		}
   	}
}

//Condition of Probability Matrix
for(int id=0; id <ndef; id++){
        for (int ij=0; ij<njet; ij++){
                for (int ik=0; ik<nkappa; ik++){
     			for(int ipt = 0 ; ipt < njetptmn ; ipt++){
      
		TH2D* RM  = (TH2D*)h2dGenDetMC1[id][ij][ik][ipt]->Clone();
      		//TH1D* miss = (TH1D*)MC_miss2[id][ij][ik][ipt]->Clone();
      		RM->RebinX(2); //miss->Rebin(2);
		//cout << "Checking Condition Number: "<<endl;
      		cout <<setw(3) << id <<setw(3) << ij << setw(3) << ik << setw(3) << ipt <<'\n';
      		//Condition(RM, miss);
		Condition(RM);
		//ConditionV2(RM);
			}
     		}
   	}
}

//--------------------------------------------------------------
Unfold->cd();
for(int id=0; id <ndef; id++){
	for (int ij=0; ij<njet; ij++){
        	for (int ik=0; ik<nkappa; ik++){
     			for(int ipt = 0 ; ipt < njetptmn ; ipt++){
       
       		//Print the reco bins
       		//cout <<"Type = " << ity <<" var " << ivar << " HT2 =" << ipt << endl;
       		//cout <<" Reco Bin : {"; 
       		double rxbins[MC_Reco[id][ij][ik][ipt]->GetNbinsX()+1]={0};
       		for (int ix=0; ix<MC_Reco[id][ij][ik][ipt]->GetNbinsX()+1; ix++) {
	 		rxbins[ix] = MC_Reco[id][ij][ik][ipt]->GetXaxis()->GetBinLowEdge(ix+1);
			//cout <<  rxbins[ix] <<", " ; 
       			} 
		//cout << endl; 
      		//cout <<" Gen Bin : {"; 
       		//Get Gen bins
       		double gxbins[MC_Gen[id][ij][ik][ipt]->GetNbinsX()+1]={0};
       		for (int ix=0; ix<MC_Gen[id][ij][ik][ipt]->GetNbinsX()+1; ix++) {
	 		gxbins[ix] = MC_Gen[id][ij][ik][ipt]->GetXaxis()->GetBinLowEdge(ix+1);
			//cout <<  gxbins[ix] <<", "; 
       			}
      		//cout << endl; 
      		//cout <<" 2d Reco  Bin : {"; 
       		//Get Gen bins
       		for (int ix=0; ix<h2dGenDetMC[id][ij][ik][ipt]->GetNbinsX()+1; ix++) {
			//cout <<  h2dGenDetMC[ity][ivar][ipt]->GetXaxis()->GetBinLowEdge(ix+1) <<" , ";
       			}
     		//cout << endl; 
      		//cout <<" 2D Gen  Bin : {"; 
       		for (int ix=0; ix<h2dGenDetMC[id][ij][ik][ipt]->GetNbinsY()+1; ix++) {
			//cout <<  h2dGenDetMC[ity][ivar][ipt]->GetYaxis()->GetBinLowEdge(ix+1) <<" , ";
       			}
       
       		cout << endl; 
       		sprintf(name,"Effi_d%i_j%i_k%i_pt%i_eta0", id, ij, ik, ipt);
       		hist_eff[id][ij][ik][ipt] = new TH1D(name,name,MC_Gen[id][ij][ik][ipt]->GetNbinsX(),gxbins);
       		hist_eff[id][ij][ik][ipt]->Sumw2();
       		sprintf(name,"Fake_d%i_j%i_k%i_pt%i_eta0", id, ij, ik, ipt);
       		hist_fake[id][ij][ik][ipt] = new TH1D(name,name, MC_Reco[id][ij][ik][ipt]->GetNbinsX(),rxbins);
       		hist_fake[id][ij][ik][ipt]->Sumw2();
       		sprintf(name,"Purity_d%i_j%i_k%i_pt%i_eta0", id, ij, ik, ipt);
       		hist_purity[id][ij][ik][ipt] =  new TH1D(name,name, MC_Reco[id][ij][ik][ipt]->GetNbinsX(),rxbins);
       		hist_purity[id][ij][ik][ipt]->Sumw2();
       		sprintf(name,"stability_d%i_j%i_k%i_pt%i_eta0", id, ij, ik, ipt);
       		hist_stbl[id][ij][ik][ipt] =  new TH1D(name,name, MC_Reco[id][ij][ik][ipt]->GetNbinsX(),rxbins);
       		hist_stbl[id][ij][ik][ipt]->Sumw2();
       
       		//-----------------------------------------------for cross check
       
		sprintf(name,"Effi1_d%i_j%i_k%i_pt%i_eta0", id, ij, ik, ipt);
       		hist_eff1[id][ij][ik][ipt] = new TH1D(name,name,MC_Gen[id][ij][ik][ipt]->GetNbinsX(),gxbins);
       		hist_eff1[id][ij][ik][ipt]->Sumw2();
       		sprintf(name,"Fake1_d%i_j%i_k%i_pt%i_eta0", id, ij, ik, ipt);
       		hist_fake1[id][ij][ik][ipt] = new TH1D(name,name, MC_Reco[id][ij][ik][ipt]->GetNbinsX(),rxbins);
       		hist_fake1[id][ij][ik][ipt]->Sumw2();
       		sprintf(name,"Purity1_d%i_j%i_k%i_pt%i_eta0", id, ij, ik, ipt);
       		hist_purity1[id][ij][ik][ipt] =  new TH1D(name,name, MC_Reco[id][ij][ik][ipt]->GetNbinsX(),rxbins);
       		hist_purity1[id][ij][ik][ipt]->Sumw2();
       		sprintf(name,"stability1_d%i_j%i_k%i_pt%i_eta0", id, ij, ik, ipt);
       		hist_stbl1[id][ij][ik][ipt] =  new TH1D(name,name, MC_Reco[id][ij][ik][ipt]->GetNbinsX(),rxbins);
       		hist_stbl1[id][ij][ik][ipt]->Sumw2();
       
		//-----------------------------------------------for cross check
       
       		sprintf(name,"No_Reg_d%i_j%i_k%i_pt%i_eta0", id, ij, ik, ipt);
       		unfold_NoReg[id][ij][ik][ipt] = new TH1D(name,name, MC_Gen[id][ij][ik][ipt]->GetNbinsX(),gxbins);
       		unfold_NoReg[id][ij][ik][ipt]->Sumw2();
       		//unfold_NoReg[ity][ivar][ipt]->Rebin(2);
       
       		sprintf(name,"TauScan_d%i_j%i_k%i_pt%i_eta0", id, ij, ik, ipt);
       		unfold_Tau[id][ij][ik][ipt] = new TH1D(name,name, MC_Gen[id][ij][ik][ipt]->GetNbinsX(),gxbins);
       		unfold_Tau[id][ij][ik][ipt]->Sumw2();
       		//unfold_Tau[ity][ivar][ipt]->Rebin(2);
       
       		sprintf(name,"Lcurve_d%i_j%i_k%i_pt%i_eta0", id, ij, ik, ipt);
       		unfold_L[id][ij][ik][ipt] = new TH1D(name,name, MC_Gen[id][ij][ik][ipt]->GetNbinsX(),gxbins);
       		unfold_L[id][ij][ik][ipt]->Sumw2();
       		//unfold_L[ity][ivar][ipt]->Rebin(2);
       
       
       		sprintf(name,"Cov_Matrix_NoReg_d%i_j%i_k%i_pt%i_eta0", id, ij, ik, ipt);
       		COV_Mat_NoReg[id][ij][ik][ipt] = new TH2D(name,name, MC_Reco[id][ij][ik][ipt]->GetNbinsX(), rxbins, MC_Gen[id][ij][ik][ipt]->GetNbinsX(), gxbins);
       		COV_Mat_NoReg[id][ij][ik][ipt]->Sumw2();
       		//COV_Mat_NoReg[ity][ivar][ipt]->RebinY(2);
       
       		sprintf(name,"Corr_Matrix_NoReg_d%i_j%i_k%i_pt%i_eta0", id, ij, ik, ipt);
       		corr_mat_NoReg[id][ij][ik][ipt] = new TH2D(name,name, MC_Gen[id][ij][ik][ipt]->GetNbinsX(), gxbins, MC_Gen[id][ij][ik][ipt]->GetNbinsX(), gxbins);
       		corr_mat_NoReg[id][ij][ik][ipt]->Sumw2();
       		//corr_mat_NoReg[ity][ivar][ipt]->RebinY(2);
       		//corr_mat_NoReg[ity][ivar][ipt]->RebinX(2);
			}
     		}
   	}
}
 
//cross check efficincy and purity calculation (Tunfold example code)
for(int id=0; id <ndef; id++){
	for (int ij=0; ij<njet; ij++){
        	for (int ik=0; ik<nkappa; ik++){
     			for(int ipt = 0 ; ipt < njetptmn ; ipt++){
       //----------------------------------------
       		for(int binGen=0;binGen<= h2dGenDetMC[id][ij][ik][ipt]->GetNbinsY()+1;binGen++) {
         		double sum0=0.;
         		double sum1=0.;
         		for(int binRec=0;binRec<= h2dGenDetMC[id][ij][ik][ipt]->GetNbinsX()+1; binRec++) {
	   			//double c=  h2dGenDetMC[ity][ivar][ipt]->GetBinContent(binGen,binRec);
	   			double c=  h2dGenDetMC[id][ij][ik][ipt]->GetBinContent(binRec,binGen);
	   			sum0+=c;
	   			if((binRec>0)&&(binRec<=h2dGenDetMC[id][ij][ik][ipt]->GetNbinsX())) {
	     				sum1+=c;
	   				}
         			}
         			if(sum0>0.0) {
	   				hist_eff1[id][ij][ik][ipt]->SetBinContent(binGen,sum1/sum0);
         				}
       			}
       //---------------------------
       		hist_eff1[id][ij][ik][ipt]->SetMinimum(0.9); hist_eff1[id][ij][ik][ipt]->SetMaximum(1.1);
       		hist_eff1[id][ij][ik][ipt]->Write();
       
       		for(int binRec=0; binRec<= hist_purity1[id][ij][ik][ipt]->GetNbinsX()+1; binRec++) {
       		double sum=0.;
         	for(int binGen=0; binGen<=hist_purity1[id][ij][ik][ipt]->GetNbinsX()+1; binGen++) {
	   		//sum += h2dGenDetMC[ity][ivar][ipt]->GetBinContent(binGen,binRec);
	   		sum += h2dGenDetMC[id][ij][ik][ipt]->GetBinContent(binRec,binGen);
         		}
         		double p=0.;
         		if(sum>0.0) {
	   			p = h2dGenDetMC[id][ij][ik][ipt]->GetBinContent(binRec,binRec)/sum;
         			}
         		hist_purity1[id][ij][ik][ipt]->SetBinContent(binRec,p);
       			}

       		hist_purity1[id][ij][ik][ipt]->SetMinimum(-0.05); hist_purity1[id][ij][ik][ipt]->SetMaximum(1.01);
       		hist_purity1[id][ij][ik][ipt]->Write();
     			}
   		}
  	}
}
//cross check efficincy and purity calculation
for(int id=0; id <ndef; id++){
	for (int ij=0; ij<njet; ij++){
        	for (int ik=0; ik<nkappa; ik++){
     			for(int ipt = 0 ; ipt < njetptmn ; ipt++){

       		double fakerate[100] = {0.};
       		double eff[100] ;
       		for(int jj=0; jj<100; jj++){
	 		eff[jj] = 1.;
       			}
       		double purity[100] = {0.};
       		double stbl[100] = {0.};
       
       		//lastbin[ij][jk] = LastBin_Counter(hist_data[ij][jk]);
       		//For Purity and Stability 
       
		subtract_background(h2dGenDetMC[id][ij][ik][ipt],MC_Reco[id][ij][ik][ipt],MC_Gen[id][ij][ik][ipt],Data_Reco[id][ij][ik][ipt],fakerate,eff,purity,stbl); //No correction to Reco or Gen
       
       		for(int bn=0; bn<(MC_Gen[id][ij][ik][ipt]->GetNbinsX()); bn++){
         	//hist_eff[ity][ivar][ipt]->SetBinContent(bn+1,eff[bn+1]);
        	//hist_fake[ity][ivar][ipt]->SetBinContent(bn+1,fakerate[bn+1]) ;
         		hist_purity[id][ij][ik][ipt]->SetBinContent(bn+1,purity[bn+1]);
         		hist_stbl[id][ij][ik][ipt]->SetBinContent(bn+1,stbl[bn+1]);
       			}

       		hist_fake[id][ij][ik][ipt] = (TH1D*)MC_fake[id][ij][ik][ipt]->Clone();
      	 	hist_eff[id][ij][ik][ipt]  = (TH1D*)MC_miss[id][ij][ik][ipt]->Clone();
       
       		hist_eff[id][ij][ik][ipt]->SetMinimum(-0.01); hist_eff[id][ij][ik][ipt]->SetMaximum(1.01);
       		hist_fake[id][ij][ik][ipt]->SetMinimum(-0.1); hist_fake[id][ij][ik][ipt]->SetMaximum(1.01);
       		hist_purity[id][ij][ik][ipt]->SetMinimum(-0.01); hist_purity[id][ij][ik][ipt]->SetMaximum(1.01);
       		hist_stbl[id][ij][ik][ipt]->SetMinimum(-0.01); hist_stbl[id][ij][ik][ipt]->SetMaximum(1.01);
       
       		hist_eff[id][ij][ik][ipt]->Write();
       		hist_fake[id][ij][ik][ipt]->Write();
       		hist_purity[id][ij][ik][ipt]->Write();
       		hist_stbl[id][ij][ik][ipt]->Write();
     			}
   		}
 	}
}

for(int id=0; id <ndef; id++){
	for (int ij=0; ij<njet; ij++){
        	for (int ik=0; ik<nkappa; ik++){
     			for(int ipt = 0 ; ipt < njetptmn ; ipt++){

       		//if (ivar==2) continue; 
       		//cout <<"type "<< ity << " : Variables " << ivar << " HT2 Bin : " << ipt << endl;
       		//file <<"["<< id << "," << ij << "," << ik << ","<< ipt <<"] --->" << endl;
     
     	  	//Get reco bins
       		double rxbins[MC_Reco[id][ij][ik][ipt]->GetNbinsX()+1]={0};
       		for (int ix=0; ix<MC_Reco[id][ij][ik][ipt]->GetNbinsX()+1; ix++) {
	 		rxbins[ix] = MC_Reco[id][ij][ik][ipt]->GetXaxis()->GetBinLowEdge(ix+1); 
			}
       
       		//Get Gen bins
       		double gxbins[MC_Gen[id][ij][ik][ipt]->GetNbinsX()+1]={0};
       		for (int ix=0; ix<MC_Gen[id][ij][ik][ipt]->GetNbinsX()+1; ix++) {
	 		gxbins[ix] = MC_Gen[id][ij][ik][ipt]->GetXaxis()->GetBinLowEdge(ix+1); }
       
       		//Rebin for match the condition of reco vs gen bin 
          	h2dGenDetMC[id][ij][ik][ipt]->RebinY(irbin);  // check rebin part
          	MC_Gen[id][ij][ik][ipt]->Rebin(irbin);
       
		//h2dGenDetMC[ity][ivar][ipt]->Rebin(1,2);
         
       		TH2* hist_migrationCoarseFine_MC = (TH2D*)h2dGenDetMC[id][ij][ik][ipt]->Clone();
       		TH1* input = (TH1D*)Data_Reco[id][ij][ik][ipt]->Clone();
       		TH1* mcgen = (TH1D*)MC_Gen[id][ij][ik][ipt]->Clone();
       
       		//correction for Fake as background subtraction : Patrick
       		TH1* mcbackground = (TH1D*)MC_fake[id][ij][ik][ipt]->Clone();
       		mcbackground->Multiply(input);  // multiply or subtraction

       		//double biasScale =bias[ivar] ;
       		double biasScale = 0.0;
       		const char *REGULARISATION_DISTRIBUTION=0;
       		const char *REGULARISATION_AXISSTEERING="*[UOB]";
       		char Rhoname[100], Rhotitle[100], Lcursure[100], Lcurtitle[100], TauLsure[100], TuaLtitle[100], probMat[100], probmat_title[100], Rhoname2d[100], Rhotitle2d[100], Ematrix[100], Ematrixtitle[100], foldback[100], foldback_title[100];
       
       		//preserve the area
       		//TUnfold::EConstraint constraintMode= TUnfold::kEConstraintArea;
       		//TUnfoldDensity::EConstraint constraintMode= TUnfoldDensity::kEConstraintArea;
       		//TUnfold::EConstraint constraintMode= TUnfold::kEConstraintNone;
       		//TUnfold::ERegMode regMode = TUnfold::kRegModeCurvature;
       		//TUnfold::ERegMode regMode = TUnfold::kRegModeSize;
       		//TUnfoldDensity::ERegMode regMode = TUnfoldDensity::kRegModeCurvature;
       		//TUnfoldDensity::EDensityMode densityFlags = TUnfoldDensity::kDensityModeNone;
       		//TUnfoldDensity::EDensityMode densityFlags = TUnfoldDensity::kDensityModeBinWidth;
       
       		//https://root.cern.ch/doc/master/testUnfold5d_8C.html      : this get input covariance matrix : Data covariance matrix
       		sprintf(name,"Data_covariance_d%i_j%i_k%i_pt%i_eta0", id, ij, ik, ipt);
       		TH2D* covarianceM = new  TH2D(name,name, input->GetNbinsX(), rxbins, input->GetNbinsX(), rxbins);
       		covarianceM->Sumw2();
       		for (int ix=1; ix<input->GetNbinsX()+1; ix++) {
          	double err = input->GetBinError(ix);
	  		covarianceM->SetBinContent(ix,ix,err*err);
       			}
       
       		covarianceM->Write();
       
       		double taumx = 0.; double taumi = 0.;
       		double tau = 1e-4;
       
       		//TUnfoldDensity class
       		//No regularisation --------------------------------
       		TUnfoldDensity tunfoldNoRegularisation(hist_migrationCoarseFine_MC, // Response matrix
					      TUnfold::kHistMapOutputVert,          // truth level on y-axis of response matrix
					      TUnfoldDensity::kRegModeNone,         // without regularisation
					      TUnfoldDensity::kEConstraintNone,     // no constrain on area
					      TUnfoldDensity::kDensityModeNone);//,0,0,REGULARISATION_DISTRIBUTION,REGULARISATION_AXISSTEERING);//,binningCoarseGen, binningFineReco);
       
       		// TUnfoldDensity::kDensityModeBinWidthAndUser);//,0,0,REGULARISATION_DISTRIBUTION,REGULARISATION_AXISSTEERING);//,binningCoarseGen, binningFineReco);
       
       		tunfoldNoRegularisation.SubtractBackground(mcbackground, "Background", 1.0, 0.03); // hist,name,scale, scale error 
       		int status = tunfoldNoRegularisation.SetInput(input,biasScale,0,covarianceM);
       
       		int nBadErrors = status%10000, nUnconstrOutBins = status/10000;
       		cout << nBadErrors << " bad errors and " << nUnconstrOutBins << " unconstrained output bins" << endl;
       		//tunfoldNoRegularisation.SubtractBackground(mcbackground,"Background", 1.0,0.04); // hist,name,scale, scale error 
       
       		//if(>=10000) { std::cout<<"Unfolding result may be wrong\n";  }
       
       		//the initial bias vector is determined from the response matrix
      		//but may be changed by using this method https://root.cern.ch/doc/master/classTUnfold.html#a58a869050370480d020ece2df3eb2688
       		//tunfoldNoRegularisation.SetBias(mcgen);   //not much affect on result
       
       		//if(ity==0 && ivar == 0) {tunfoldNoRegularisation.RegularizeBins(7,1,4,regMode);}  //Test for Regularization in few  unmatched bins
       		//if(ity==0 && ivar == 0) {tunfoldNoRegularisation.RegularizeCurvature(10,9,11,1.0,0.6);}
 
       		//Choose a value of tau to unfold with tau,0 means no regularization
       		tunfoldNoRegularisation.DoUnfold(0.0);//,input,biasScale);
       
       		/* //Binmaps : Thanks to Suman(Tifr)
       		int binnos = rnbinsx[ity][ivar][ipt];
       		Int_t *binMap=new Int_t[binnos+2];
       		for(Int_t i=1;i<=binnos;i++) binMap[i]=i;
       		binMap[0]=-1;   binMap[binnos+1]=-1;
       		*/

       		char unfoldhist[100], title[100], NoReg_RhoIJ[100], NoReg_RhoIJ_tit[100], NoReg_Ematrix[100], NoReg_Ematrix_tit[100];// NoReg_prob[100], NoReg_probtitle[100];
       		sprintf(unfoldhist, "Tunfold_Noreg_D%i_j%i_k%i_pt%i_eta0", id, ij, ik, ipt); //unfolded_typ_0_pt2_eta0_3
       		sprintf(title, "Tunfolded Noreg %i 2.5 %s %s %s", int(leadingPtThreshold[ipt]), obs_def[id], jet_num[ij], k_fact[ik] );
       		sprintf(NoReg_RhoIJ, "Tunfold_Noreg_corr_D%i_j%i_k%i_pt%i_eta0", id, ij, ik, ipt); //unfolded_typ_0_pt2_eta0_3
       		sprintf(NoReg_RhoIJ_tit, "2D correlation coefficients No Regularisation %i 2.5 %s %s %s", int(leadingPtThreshold[ipt]), obs_def[id], jet_num[ij], k_fact[ik]);
       		sprintf(NoReg_Ematrix, "Tunfold_Noreg_Emat_typ_D%i_j%i_k%i_pt%i_eta0", id, ij, ik, ipt); //unfolded_typ_0_pt2_eta0_3
       		sprintf(NoReg_Ematrix_tit, "EMatrix No Regularisation %i 2.5 %s %s %s", int(leadingPtThreshold[ipt]), obs_def[id], jet_num[ij], k_fact[ik]);
       		sprintf(probMat, "Tunfold_Noreg_probM_D%i_j%i_k%i_pt%i_eta0", id, ij, ik, ipt);
       		sprintf(probmat_title, "Probability matrix Noreg %i 2.5 %s %s %s", int(leadingPtThreshold[ipt]), obs_def[id], jet_num[ij], k_fact[ik]);
       		sprintf(foldback, "Tunfold_NoReg_Refold_D%i_j%i_k%i_pt%i_eta0", id, ij, ik, ipt); //unfolded_typ_0_pt2_eta0_3
       		sprintf(foldback_title, "TunfoldFolded back  NoReg %i 2.5 %s %s %s", int(leadingPtThreshold[ipt]), obs_def[id], jet_num[ij], k_fact[ik]);

	       	//unfolding result, signal only
       		TH1 *hist_PTunfolded_noRegularisation = tunfoldNoRegularisation.GetOutput(unfoldhist,title);//,0,"*[UO]" ,true);//,"","*[UO]");//,"signal");
       		//TH1 *hist_PTunfolded_noRegularisation = tunfoldNoRegularisation.GetOutput(unfoldhist,title);//,"","*[UO]");//,"signal");
       		//TH1 *hist_PTunfolded_noRegularisation = tunfoldNoRegularisation.GetOutput("hist_PTunfolded_noRegularisation", "P_{T,unfolded} [GeV]","signal");

       		TH2 *hist_RhoIJ_noRegularisation = tunfoldNoRegularisation.GetRhoIJtotal(NoReg_RhoIJ, NoReg_RhoIJ_tit);//,"signal");
       		TH2 *hist_Rho2D_noRegularisation = tunfoldNoRegularisation.GetEmatrixTotal(NoReg_Ematrix, NoReg_Ematrix_tit);//,"signal");
       		//tunfoldNoRegularisation.GetEmatrix(COV_Mat_NoReg[ity][ivar][ipt],binMap);//,"signal");
       		TH1 *hist_foldedback_NoReg = tunfoldNoRegularisation.GetFoldedOutput(foldback, foldback_title);//,"signal");
       		TH2 *hist_prob_Noreg = tunfoldNoRegularisation.GetProbabilityMatrix(probMat, probmat_title,TUnfold::kHistMapOutputVert);//,"signal");

       		//correction for Miss entries  : Partick 
       		for (int i = 1; i <= hist_PTunfolded_noRegularisation->GetNbinsX(); ++i) {
	 		double content = hist_PTunfolded_noRegularisation->GetBinContent(i);
	 		double factor = 1;
	 		factor += MC_miss1[id][ij][ik][ipt]->GetBinContent(i);
         		content *= factor;
	 		hist_PTunfolded_noRegularisation->SetBinContent(i, content);
       			}
      
       		/*  //Quick check for closure Ratio Plots
       		hist_PTunfolded_noRegularisation->Scale(1/(hist_PTunfolded_noRegularisation->Integral()));
       		mcgen->Scale(1/mcgen->Integral());
       		hist_PTunfolded_noRegularisation->Divide(mcgen);
       		hist_PTunfolded_noRegularisation->SetMinimum(0.85); hist_PTunfolded_noRegularisation->SetMaximum(1.15);
       		*/

       		hist_PTunfolded_noRegularisation->Write();
       		hist_Rho2D_noRegularisation->Write();
       		hist_RhoIJ_noRegularisation->Write();
       		hist_foldedback_NoReg->Write();
       		hist_prob_Noreg->Write();
       
       		cout << "Without Regularization finish " << endl;
       
	     		}
   		}
	}
}
 
delete outputFile;
//file.close();
}

//--------------FUNCTIONS-------------------

//Copy of RooUnfold Subtract Background method with no correction to Reco or Gen --->Can be used  for purity and stability
int subtract_background(TH2D* h2d_correl, TH1D* reco, TH1D* gen, TH1D* data, double* fakerate, double* effi, double* purity, double* stbl) {
	int nbinx = h2d_correl->GetNbinsX();
  	int nbiny = h2d_correl->GetNbinsY(); 
  	const int nbinmx = 100 ;
  	double totalgen[nbinmx]={0.};
  	double totalreco[nbinmx]={0.};
  	for (int ix=0; ix<nbinx+1; ix++) {
    		for (int iy=0; iy<nbiny+1; iy++) {      
      			if(ix==0&&iy==0) continue ;
      				totalreco[ix] +=h2d_correl->GetBinContent(ix, iy);
      			if (iy==0) fakerate[ix-1] = h2d_correl->GetBinContent(ix,iy);
      				totalgen[iy] +=h2d_correl->GetBinContent(ix, iy);
      			if (ix==0) effi[iy-1] =h2d_correl->GetBinContent(ix, iy);  
    			}//iy
  		}//ix
  
  	for (int iy=0; iy<nbiny; iy++) {
    		effi[iy] = (totalgen[iy+1] - effi[iy])/max(1.e-10, totalgen[iy+1]);
   		//if(iy>10 && effi[iy]>0.1 && effi[iy]<0.9){
   			//effi[iy] = effi[iy-1]-0.001;
   			//}
  		//if (effi[iy]<1.e-6) effi[iy]=1.e-6;
   		} //iy
  
  	for (int ix=0; ix<nbinx; ix++){
    		if(totalreco[ix+1]>1.e-6) { fakerate[ix] /= totalreco[ix+1]; }
    		else { fakerate[ix] = 0.; }
    
    		//if((ix>20&&fakerate[ix]>0.05) || (ix>10&&fakerate[ix]>0.2)) { fakerate[ix] = 1.e-3;  }
  		}//ix
  
 	//for(int ix=0; ix <((reco->GetNbinsX())+1); ix++){
    		//reco->SetBinContent(ix+1,(1-fakerate[ix])*(reco->GetBinContent(ix+1)));
    		//reco->SetBinError(ix+1,sqrt(1-fakerate[ix])*(reco->GetBinError(ix+1)));
 		//}
  
 	//for(int ix=0; ix <((data->GetNbinsX())+1); ix++){
    		//data->SetBinContent(ix+1,(1-fakerate[ix])*(data->GetBinContent(ix+1)));
    		//data->SetBinError(ix+1,sqrt(1-fakerate[ix])*(data->GetBinError(ix+1)));
 		//}
  
	//for(int ix=0; ix<(gen->GetNbinsX()+1); ix++){
    		//gen->SetBinContent(ix+1,(gen->GetBinContent(ix+1))*effi[ix]) ;
    		//gen->SetBinError(ix+1,sqrt(gen->GetBinError(ix+1))*effi[ix]) ;
		//}
  
  	double tot_reco_in[nbinmx] = {0.};
  	double tot_gen_in[nbinmx] = {0.};
  
  	for(int bn=0; bn<nbinx; bn++){
    		for(int am=0; am<nbiny; am++){
      			tot_gen_in[am]+=h2d_correl->GetBinContent(bn+1,am+1);
      			tot_reco_in[bn]+=h2d_correl->GetBinContent(bn+1,am+1);
    			}
  		}
  
  	for(int bn=0; bn<nbinx; bn++){
    		purity[bn] = (tot_reco_in[bn]>1.e-7)?h2d_correl->GetBinContent(bn+1,bn+1)*1./tot_reco_in[bn]:0.;
    		stbl[bn] = (tot_gen_in[bn]>1.e-7)?h2d_correl->GetBinContent(bn+1,bn+1)*1./tot_gen_in[bn]:0.;
  		}
  	
	return 0;
} //end substract func

void Fold(TH2D* HistoMatrix, TH1D* HistReco, TH1D* HistGen, TH1D* miss, TH1D* fake, TH1D* HistoCorrect){
	TH2D* Histprob = (TH2D*) HistoMatrix->Clone(); Histprob->Reset();
	//calculate Probability Matrix
	for(int ij=0; ij<(HistoMatrix->GetNbinsY()+2); ij++){
	double row_sum = 0.;
		for(int jk=0; jk<(HistoMatrix->GetNbinsX()+2); jk++){
  			row_sum+=HistoMatrix->GetBinContent(jk,ij);
			}//jk
		if(row_sum>1.e-10){
	 		for(int jk=0; jk<(HistoMatrix->GetNbinsX()+2); jk++){
   				Histprob->SetBinContent(jk,ij,(HistoMatrix->GetBinContent(jk,ij)*1./row_sum)) ; //Probability
  				}//jk
			}
		}//ij

//folding gen level to Reco
	for(int i=0;i<Histprob->GetNbinsX()+2;i++){
     	double sum=0.; double Err =0.;
       		for(int j=0;j<HistGen->GetNbinsX()+2;j++){
       		double misscorr = (HistGen->GetBinContent(j))-(miss->GetBinContent(j)); //Miss correction
       			sum += Histprob->GetBinContent(i,j)*misscorr;
       			Err += (Histprob->GetBinContent(i,j)*HistGen->GetBinError(j))*(Histprob->GetBinContent(i,j)*(HistGen->GetBinError(j)));
           		}
   			//sum = sum +(fake->GetBinContent(i)); //fake correction
      			HistoCorrect->SetBinContent(i,sum);
     	 		HistoCorrect->SetBinError(i,sqrt(Err));
    			}
}//end Fold
/*
//Condition number calculation
//void Condition (TH2* RM, TH1* miss){
void Condition (TH2* RM){
	const int Nx = RM->GetNbinsX(),
       		  Ny = RM->GetNbinsY();
    
    	cout << Nx <<"  " << Ny << endl;
    	cout << "Testing 0"<<endl;
    	//if (Ny*2 != Nx) { cout << Nx << ' ' << Ny << endl;  continue; }
		cout << "Testing 1" <<endl;
    	//TH1D* RMy = RM->ProjectionY("RMy", 0, -1); //Gen Projection

    	//normalisation & condition
 
    	TMatrixD m(Ny,Nx);
    		for (int i = 1; i <=Nx; ++i) {
			for (int j = 1; j <=Ny; j++ ){
				m(j-1,i-1) = RM->GetBinContent(i,j);
				}
			}
        	//double normalisation = RMy->GetBinContent(i);
               	//normalisation += miss->GetBinContent(i);
        	//	if (normalisation > 0)
		//		cout <<"Testing 2"<<endl;
        	//	for (int j = 1; j <= Nx; ++j) {
		//		cout <<"Testing 3"<<endl;
            	//		double content = RM->GetBinContent(j,i);
            	//			content /= normalisation;
            	//			m((j-1)/2,i-1) += content;
        	//		}
    		//	}
 
	//TMatrixD  matd(Ny, Ny);
	//for( int i =0; i <Ny ; i++){
    	//	for( int j =0; j <Ny ; j++){
	//		matd[i][j] = RM->GetBinContent( j+1, i+1);
    	//		}
	//	}
    	cout << "Testing 4"<<endl;
//TDecompSVD svd(m);
//bool ok = svd.Decompose();
//TMatrixD b;
//if (ok)
//b = svd.Invert();
//else {
//cout << "SVD failed, condition: " << svd.Condition() <<endl;
//m.Print();
//}
    	TDecompSVD svd(m);
    	//TDecompSVD svd(Ny,Nx);
    	//svd.SetMatrix(matd);
    //	TVectorD v = svd.GetSig();
    	cout << "Condition Loop"<<endl;
    	cout << "Condition : " << svd.Condition() << endl;
    //	double Max = v[0];
    //	for(int k = 0; k < Ny; ++k) {
//		if (abs(v[k]) < feps) break;
        		//cout << setw(5) << k << setw(15) << v[k] << setw(15) << Max/v[k] << '\n';
  //  		}
}// End condition number
*/
/*
void Condition (TH2* RM, TH1* miss){
    const int Nx = RM->GetNbinsX(),
              Ny = RM->GetNbinsY();
    cout << Nx <<"  " << Ny << endl;
    //if (Ny*2 != Nx) { cout << Nx << ' ' << Ny << endl;  return; }

    TH1D* RMy = RM->ProjectionY("RMy", 0, -1); //Gen Projection

    // normalisation & condition
    TMatrixD m(Ny,Ny);
    for (int i = 1; i <= Ny; ++i) {
        double normalisation = RMy->GetBinContent(i);
               normalisation += miss->GetBinContent(i);
        if (normalisation > 0)
        for (int j = 1; j <= Nx; ++j) {
            double content = RM->GetBinContent(j,i);
            content /= normalisation;
            m((j-1)/2,i-1) += content;
        }
    }
    TDecompSVD svd(m);
    TVectorD v = svd.GetSig();
	cout << "Condition :"<<svd.Condition()<<endl;
    double Max = v[0];
    for(int k = 0; k < Ny; ++k) {
        if (abs(v[k]) < feps) break;
        cout << setw(5) << k << setw(15) << v[k] << setw(15) << Max/v[k] << '\n';
    }
}
*/
double Condition (TH2 * RM)
{
    cout << __LINE__ << "\tCalculating condition" << endl;
    const int Nx = RM->GetNbinsX(),
              Ny = RM->GetNbinsY();

    TH1 * RMy = RM->ProjectionY("RMy", 0, -1);

    // normalisation & condition
    TMatrixD m(Ny,Ny); // unfortunately, we have to swap the axes...
    for (int j = 1; j <= Ny; ++j) {
        double normalisation = RMy->GetBinContent(j);
        if (normalisation > 0)
        for (int i = 1; i <= Nx; ++i)
            m(j-1,i-1) = RM->GetBinContent(i,j) / normalisation;
    }
    TDecompSVD svd(m);
    TVectorD v = svd.GetSig();
    cout << "Condition :"<<svd.Condition()<<endl;
    double Max = v[0], Min = v[0];
    for(int k = 0; k < Ny; ++k) {
        if (abs(v[k]) < 1e-5) break;
        Min = v[k];
        cout << setw(5) << k << setw(15) << v[k] << setw(15) << Max/Min << '\n';
    }
    if (Min > 0) return Max/Min;
    else         return -numeric_limits<double>::infinity();
}
double ConditionV2(TH2* RM){
	const int nBinsX = RM->GetNbinsX(),
	          nBinsY = RM->GetNbinsY();	
    
	if(nBinsX!=nBinsY){
        	cout << endl;
        	cout << "x and y bins do not match" << endl;
        	cout << endl;
    	}

	TMatrixD matrix(nBinsY,nBinsX);
	for(int i=0;i<=nBinsX;i++){
		for(int j=0;j<=nBinsY;j++){
			matrix(j-1,i-1) = RM->GetBinContent(i,j);
		}
	}

	TDecompSVD decomp(matrix);
	cout << "The condition number: " << decomp.Condition() << endl;

	//double determinant;
	//TMatrixD mInverse = matrix.Invert(&determinant);
	//cout << "The determinant: " << determinant << endl;
    	//_determinant = determinant;
}
