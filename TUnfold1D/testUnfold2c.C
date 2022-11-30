//Version-2
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
#include <TProfile.h>
#include <THStack.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLine.h>
#include <TLorentzVector.h>
#include <iomanip>

#include "TUnfoldBinning.h"
#include "TUnfoldBinningXML.h"
#include "TUnfoldDensity.h"
#include "TUnfold.h"
#include "TUnfoldIterativeEM.h"
#include "TUnfoldSys.h"
#include <TVectorD.h>
#include <TDecompSVD.h>

#define CLOSURE
//#define BLTest

using namespace std;
static const auto feps = numeric_limits<float>::epsilon();

void testUnfold2c1D_v2()
{
  //switch on histogram errors
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  //MC -> 0,1,2,3 PY8, PY8_Flat, MG, HW7
  int const nmc=4;   //Number of MC used
  int const umc=3;   //MC will used for Unfold
  int const idata=3; //MC will used for Unfold
  int irbin = 1;     //Rebin

  //const TString Pyinput = "../Input/26Nov2022/PY8_bin.root"; // PY8 bin
  const TString Pyinput = "../Input/26Nov2022/PY8_flat.root"; // PY8 bin
  //const TString Pyinput = "../Input/26Nov2022/MG5_PY8_bin.root"; // PY8 bin
  //const TString Pyinput = "../Input/26Nov2022/HW7_flat.root"; // PY8 bin

  //const TString datainput = "../Input/26Nov2022/Data_UL2017.root"; // Data
  
  //const TString datainput = "../Input/26Nov2022/PY8_bin.root"; // PY8 bin for closure
  //const TString datainput = "../Input/26Nov2022/PY8_flat.root"; // PY8 bin for closure
  //const TString datainput = "../Input/26Nov2022/MG5_PY8_bin.root"; // PY8 bin for closure
  const TString datainput = "../Input/26Nov2022/HW7_flat.root"; // PY8 bin for closure

  //Input Data and MC histogram
  TFile *inputData=new TFile(datainput);
  TFile *RMinput=new TFile(Pyinput);

  TFile *inputMC[nmc];
  inputMC[0]=new TFile("../Input/26Nov2022/PY8_flat.root");     // PY8 bin
  inputMC[1]=new TFile("../Input/26Nov2022/PY8_flat.root");    // PY8 flat
  inputMC[2]=new TFile("../Input/26Nov2022/PY8_flat.root"); // MG5+PY8
  inputMC[3]=new TFile("../Input/26Nov2022/HW7_flat.root");    // HW7 flat

  TFile *outputFile=new TFile("../Unfolded/26Nov2022/1D/Unfolded_HW7_flat_PY8_flat.root","recreate");   //Unfolded Data and Covarince matrix, efficincy,fake rate, purity, stability
  
  string HistDir = "analyzeBasicPat";
  //const char* Dimtag = "2d";
  //const char* Histtag = "dd_";

  char histname[100], name[100];//, title[100];
  const int nHLTmx=10; //PT Range

  const int njetptmn = nHLTmx;
  double leadingPtThreshold[njetptmn+1] = {92, 119, 185, 251, 319, 388, 467, 518, 579, 669, 3000}; //Fit Value singlejet trigger

  static const int ndef = 3;
  static const int njet = 2;
  static const int nkappa = 10;

  const char* obs_def[ndef]={"Q","Q_{L}","Q_{T}"};
  const char* jet_num[njet]={"Leading-Jet","Sub-Leading-Jet"};
  const char* k_fact[nkappa]={"k=0.1","k=0.2","k=0.3","k=0.4","k=0.5","k=0.6","k=0.7","k=0.8","k=0.9","k=1.0"};

  TDirectoryFile *DirData = new TDirectoryFile("Data","Input Data");

  TDirectoryFile *outputDir[nmc];
  outputDir[0]=new TDirectoryFile("Pythia8"," Pythia8 , MC and Probability Matrix");
  outputDir[1]=new TDirectoryFile("Py8Flat"," Pythia8 Flat sample , MC and Probability Matrix");
  outputDir[2]=new TDirectoryFile("MG5","Madgraph, MC and Probability Matrix");
  outputDir[3]=new TDirectoryFile("HW7","Herwig7 MC and Probability Matrix");

  TDirectoryFile *DirRMinput = new TDirectoryFile("DirRMinput","MC used in Unfolding");
  TDirectoryFile *folddir = new TDirectoryFile("Folded","Gen fold with Probablility matrix");
  TDirectoryFile *Unfolddir = new TDirectoryFile("Unfold","Unfolded, Refold, correlation 1D");

//------------------------------------------------------------- 
  void setgstyle();
  TH1D* ReadHist1D(string name,TFile* root);
  TH2D* ReadHist2D(string name,TFile* root);
  void Integralhist(TH1 *hist);
  void Fold(TH2D* HistoMatrix, TH1D* HistReco, TH1D* HistoGen, TH1D* miss, TH1D* fake, TH1D* HistoCorrect);
  void BLT (TH1 * dataDist, TH2 * dataCov, TH1 * MC, int rebin = 1);
  void Chi2Root(TH1 * data, TH1 * MC, int rebin = 1);
  double Chi2(const TH1* hData, const TH2* covmat, const TH1* hGen, int skip/* = -1*/);
  //void Condition(TH2* RM, TH1* miss);
  double ConditionV2(TH2 * RM);
  double ConditionV3(TH2* RM);
  int subtract_background(TH2D* h2d_correl, TH1D* reco, TH1D* gen, TH1D* data, double* fakerate, double* effi, double* purity, double* stbl);
//-------------------------------------------------------------- 
  TH1D *Data_Reco[ndef][njet][nkappa][njetptmn];

#ifdef CLOSURE
  TH1D *Data_Gen[ndef][njet][nkappa][njetptmn];
  TH1D *Data_fake[ndef][njet][nkappa][njetptmn];
  TH1D *Data_miss[ndef][njet][nkappa][njetptmn];
  TH1D *Data_fakerate[ndef][njet][nkappa][njetptmn];
  TH1D *Data_fakerateInv[ndef][njet][nkappa][njetptmn];
  TH1D *Data_missrate[ndef][njet][nkappa][njetptmn];
  TH1D *Data_missrateInv[ndef][njet][nkappa][njetptmn];
  TH1D *Data_reco_fake[ndef][njet][nkappa][njetptmn];
  TH1D *Data_gen_miss[ndef][njet][nkappa][njetptmn];
#endif
  
  TH1D *RMinput_Reco[ndef][njet][nkappa][njetptmn];
  TH1D *RMinput_Gen[ndef][njet][nkappa][njetptmn];
  TH1D *RMinput_fakerate[ndef][njet][nkappa][njetptmn];
  TH1D *RMinput_fakerateInv[ndef][njet][nkappa][njetptmn];
  TH1D *RMinput_missrate[ndef][njet][nkappa][njetptmn];
  TH1D *RMinput_missrateInv[ndef][njet][nkappa][njetptmn];
  TH1D *RMinput_fake[ndef][njet][nkappa][njetptmn];
  TH1D *RMinput_reco_fake[ndef][njet][nkappa][njetptmn];
  TH1D *RMinput_miss[ndef][njet][nkappa][njetptmn];
  TH1D *RMinput_gen_miss[ndef][njet][nkappa][njetptmn];
  TH2D *RMinput_RM[ndef][njet][nkappa][njetptmn];

  TH1D *MC_Reco[nmc][ndef][njet][nkappa][njetptmn];
  TH1D *MC_fake[nmc][ndef][njet][nkappa][njetptmn]; 
  TH1D *MC_fakerate[nmc][ndef][njet][nkappa][njetptmn];
  TH1D *MC_reco_fake[nmc][ndef][njet][nkappa][njetptmn];

  TH1D *MC_Gen[nmc][ndef][njet][nkappa][njetptmn];  
  TH1D *MC_miss[nmc][ndef][njet][nkappa][njetptmn];      
  TH1D *MC_missrate[nmc][ndef][njet][nkappa][njetptmn];
  TH1D *MC_Gen_miss[nmc][ndef][njet][nkappa][njetptmn];

  TH2D *h2dGenDetMC[nmc][ndef][njet][nkappa][njetptmn]; 
  TH1D* hist_purity[nmc][ndef][njet][nkappa][njetptmn];
  TH1D* hist_stbl[nmc][ndef][njet][nkappa][njetptmn];

  //TH1D *fakerate[nmc][ndef][njet][nkappa][njetptmn];
  //TH1D *missrate[nmc][ndef][njet][nkappa][njetptmn];

//-------------------------Read Input Data/MC and Response matrix--------------------
for(int id=0; id<ndef; id++){
 	for (int ij=0; ij<njet; ij++){
		for (int ik=0; ik<nkappa; ik++){
			for(int ipt=0; ipt<njetptmn; ipt++){
       	 	DirData->cd();
		Data_Reco[id][ij][ik][ipt] =(TH1D*)ReadHist1D(HistDir+"/reco_jc_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_pt"+to_string(ipt)+"_eta0",inputData)->Clone();

       		for (int ibin =1 ; ibin <  Data_Reco[id][ij][ik][ipt]->GetNbinsX()+1; ibin++ ){
       			if(Data_Reco[id][ij][ik][ipt]->GetBinContent(ibin) == 0) {cout << "Data Reco Bin is Zero for bin number : ***** "<<  ibin  << endl; }
       			}
#ifdef CLOSURE
		Data_Gen[id][ij][ik][ipt] = (TH1D*)ReadHist1D(HistDir+"/gen_jc_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_pt"+to_string(ipt)+"_eta0",inputData)->Clone();
	
		TH1D *Datafake = (TH1D*)ReadHist1D(HistDir+"/recofake_jc_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_pt"+to_string(ipt)+"_eta0",inputData)->Clone();
		TH1D *Datamiss = (TH1D*)ReadHist1D(HistDir+"/genmiss_jc_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_pt"+to_string(ipt)+"_eta0",inputData)->Clone();

		Data_fake[id][ij][ik][ipt] = (TH1D*)Datafake->Clone();
		Data_miss[id][ij][ik][ipt] = (TH1D*)Datamiss->Clone();

		cout <<"Fake = "<<Data_fake[id][ij][ik][ipt]->GetEntries() <<" Reco-fake: "<<(Data_Reco[id][ij][ik][ipt]->GetEntries() - Data_fake[id][ij][ik][ipt]->GetEntries())<<endl;
                cout << "Miss = "<<Data_miss[id][ij][ik][ipt]->GetEntries() <<" Gen-Miss: " <<(Data_Gen[id][ij][ik][ipt]->GetEntries() - Data_miss[id][ij][ik][ipt]->GetEntries())<<endl;

		Data_fakerate[id][ij][ik][ipt] = (TH1D*)Datafake->Clone();
		Data_fakerate[id][ij][ik][ipt]->Divide(Data_fakerate[id][ij][ik][ipt], Data_Reco[id][ij][ik][ipt], 1, 1, "b");	
		sprintf(name,"Data_fakerate_d%i_j%i_k%i_pt%i_eta0", id, ij, ik, ipt);
		Data_fakerate[id][ij][ik][ipt]->SetNameTitle(name,name);

		Data_fakerateInv[id][ij][ik][ipt] = (TH1D*)Data_fakerate[id][ij][ik][ipt]->Clone(); Data_fakerateInv[id][ij][ik][ipt]->Reset();
		for(int i =1; i<=Data_fakerate[id][ij][ik][ipt]->GetNbinsX(); i++){
			double factor = Data_fakerate[id][ij][ik][ipt]->GetBinContent(i);
			double content = 1-factor; Data_fakerateInv[id][ij][ik][ipt]->SetBinContent(i, content);
		}

		Data_missrate[id][ij][ik][ipt] = (TH1D*)Datamiss->Clone();
		Data_missrate[id][ij][ik][ipt]->Divide(Data_missrate[id][ij][ik][ipt], Data_Gen[id][ij][ik][ipt], 1, 1, "b");
		sprintf(name,"Data_missrate_d%i_j%i_k%i_pt%i_eta0", id, ij, ik, ipt);
		Data_missrate[id][ij][ik][ipt]->SetNameTitle(name,name);

		Data_missrateInv[id][ij][ik][ipt] = (TH1D*)Data_missrate[id][ij][ik][ipt]->Clone(); Data_missrateInv[id][ij][ik][ipt]->Reset();
		for(int i =1; i<=Data_missrate[id][ij][ik][ipt]->GetNbinsX(); i++){
                        double factor = Data_missrate[id][ij][ik][ipt]->GetBinContent(i);
                        double content = 1-factor; Data_missrateInv[id][ij][ik][ipt]->SetBinContent(i, content);
                }

		sprintf(name,"DataRecominusfake_d%i_j%i_k%i_pt%i_eta0",id,ij,ik,ipt);
		Data_reco_fake[id][ij][ik][ipt] = (TH1D*)Data_fakerateInv[id][ij][ik][ipt]->Clone();
		Data_reco_fake[id][ij][ik][ipt]->SetNameTitle(name,name);
		Data_reco_fake[id][ij][ik][ipt]->Multiply(Data_Reco[id][ij][ik][ipt]);

		sprintf(name,"DataGenminusmiss_d%i_j%i_k%i_pt%i_eta0",id,ij,ik,ipt);
		Data_gen_miss[id][ij][ik][ipt] = (TH1D*)Data_missrateInv[id][ij][ik][ipt]->Clone();
		Data_gen_miss[id][ij][ik][ipt]->SetNameTitle(name,name);
		Data_gen_miss[id][ij][ik][ipt]->Multiply(Data_Gen[id][ij][ik][ipt]);

		//Data_fakerate[id][ij][ik][ipt]->Write();
		//Data_missrate[id][ij][ik][ipt]->Write();
		
		Data_reco_fake[id][ij][ik][ipt]->Write();
		Data_gen_miss[id][ij][ik][ipt]->Write();
#endif
//---------------------------------MC for Read RM------------------------------------
		DirRMinput->cd();
		RMinput_Reco[id][ij][ik][ipt] = (TH1D*)ReadHist1D(HistDir+"/reco_jc_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_pt"+to_string(ipt)+"_eta0",RMinput)->Clone();
		for(int i=1; i<RMinput_Reco[id][ij][ik][ipt]->GetNbinsX()+1; i++){if (RMinput_Reco[id][ij][ik][ipt]->GetBinContent(i) ==0) {cout <<"RMinput Bin entry Nil : bin no : "<< i <<endl;}}
		
		RMinput_Gen[id][ij][ik][ipt] = (TH1D*)ReadHist1D(HistDir+"/gen_jc_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_pt"+to_string(ipt)+"_eta0",RMinput)->Clone();
		
		RMinput_RM[id][ij][ik][ipt] = (TH2D*)ReadHist2D(HistDir+"/RM_jc_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_pt"+to_string(ipt)+"_eta0",RMinput)->Clone();
		
		RMinput_fake[id][ij][ik][ipt] = (TH1D*)ReadHist1D(HistDir+"/recofake_jc_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_pt"+to_string(ipt)+"_eta0",RMinput)->Clone();
		
		RMinput_miss[id][ij][ik][ipt] = (TH1D*)ReadHist1D(HistDir+"/genmiss_jc_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_pt"+to_string(ipt)+"_eta0",RMinput)->Clone();

        	RMinput_fakerate[id][ij][ik][ipt] = (TH1D*)RMinput_fake[id][ij][ik][ipt]->Clone(); 
		RMinput_fakerate[id][ij][ik][ipt]->Divide(RMinput_fakerate[id][ij][ik][ipt], RMinput_Reco[id][ij][ik][ipt], 1, 1, "b");
		sprintf(name, "RMinputfake_rate_d%i_j%i_k%i_pt%i_eta0",id,ij,ik,ipt);
		RMinput_fakerate[id][ij][ik][ipt]->SetNameTitle(name,name); 

        	RMinput_fakerateInv[id][ij][ik][ipt] =(TH1D*)RMinput_fakerate[id][ij][ik][ipt]->Clone(); 
		RMinput_fakerateInv[id][ij][ik][ipt]->Reset();
		for(int i =1; i<=RMinput_fakerate[id][ij][ik][ipt]->GetNbinsX(); i++){ // Check Bug ??
                        double factor = RMinput_fakerate[id][ij][ik][ipt]->GetBinContent(i);
                        double content = 1-factor; RMinput_fakerateInv[id][ij][ik][ipt]->SetBinContent(i, content);
                }

		RMinput_missrate[id][ij][ik][ipt] = (TH1D*)RMinput_miss[id][ij][ik][ipt]->Clone(); 
		RMinput_missrate[id][ij][ik][ipt]->Divide(RMinput_missrate[id][ij][ik][ipt], RMinput_Gen[id][ij][ik][ipt], 1, 1, "b");
        	sprintf(name, "RMinputmiss_rate_d%i_j%i_k%i_pt%i_eta0",id,ij,ik,ipt);
        	RMinput_missrate[id][ij][ik][ipt]->SetNameTitle(name,name);

        	RMinput_missrateInv[id][ij][ik][ipt] =(TH1D*)RMinput_missrate[id][ij][ik][ipt]->Clone(); RMinput_missrateInv[id][ij][ik][ipt]->Reset();
        	for(int i =1; i<=RMinput_missrate[id][ij][ik][ipt]->GetNbinsX(); i++){
                        double factor = RMinput_missrate[id][ij][ik][ipt]->GetBinContent(i);
                        double content = 1-factor; RMinput_missrateInv[id][ij][ik][ipt]->SetBinContent(i, content);
                }

		sprintf(name,"RMinputRecominusfake_d%i_j%i_k%i_pt%i_eta0",id,ij,ik,ipt);
        	RMinput_reco_fake[id][ij][ik][ipt] = (TH1D*)RMinput_fakerateInv[id][ij][ik][ipt]->Clone();
		RMinput_reco_fake[id][ij][ik][ipt]->SetNameTitle(name,name);
		RMinput_reco_fake[id][ij][ik][ipt]->Multiply(RMinput_Reco[id][ij][ik][ipt]);

		sprintf(name,"RMinputGenminusmiss_d%i_j%i_k%i_pt%i_eta0",id,ij,ik,ipt);
        	RMinput_gen_miss[id][ij][ik][ipt] = (TH1D*)RMinput_missrateInv[id][ij][ik][ipt]->Clone();
        	RMinput_gen_miss[id][ij][ik][ipt]->SetNameTitle(name,name);
        	RMinput_gen_miss[id][ij][ik][ipt]->Multiply(RMinput_Gen[id][ij][ik][ipt]);

		RMinput_fakerate[id][ij][ik][ipt]->Write();
		RMinput_missrate[id][ij][ik][ipt]->Write();
		RMinput_reco_fake[id][ij][ik][ipt]->Write();
		RMinput_gen_miss[id][ij][ik][ipt]->Write();
//----------------------------------Read MC------------------------------------------
	for(int imc=0; imc<nmc; imc++){
		cout << "MC number : "<<imc<<endl;
		outputDir[imc]->cd();
		MC_Reco[imc][id][ij][ik][ipt] = (TH1D*)ReadHist1D(HistDir+"/reco_jc_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_pt"+to_string(ipt)+"_eta0",inputMC[imc])->Clone();
		MC_Gen[imc][id][ij][ik][ipt] = (TH1D*)ReadHist1D(HistDir+"/gen_jc_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_pt"+to_string(ipt)+"_eta0",inputMC[imc])->Clone();
		MC_fake[imc][id][ij][ik][ipt] = (TH1D*)ReadHist1D(HistDir+"/recofake_jc_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_pt"+to_string(ipt)+"_eta0",inputMC[imc])->Clone();
		MC_miss[imc][id][ij][ik][ipt] = (TH1D*)ReadHist1D(HistDir+"/genmiss_jc_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_pt"+to_string(ipt)+"_eta0",inputMC[imc])->Clone();

		cout<<"Fake ="<< MC_fake[imc][id][ij][ik][ipt]->GetEntries() <<" Reco-fake: "<< (MC_Reco[imc][id][ij][ik][ipt]->GetEntries() - MC_fake[imc][id][ij][ik][ipt]->GetEntries()) <<endl;
		cout<<"Miss ="<< MC_miss[imc][id][ij][ik][ipt]->GetEntries() <<" Gen-miss: "<< (MC_Gen[imc][id][ij][ik][ipt]->GetEntries() - MC_miss[imc][id][ij][ik][ipt]->GetEntries()) <<endl;
		//Response Matrix
		h2dGenDetMC[imc][id][ij][ik][ipt] = (TH2D*)ReadHist2D(HistDir+"/RM_jc_d"+to_string(id)+"_j"+to_string(ij)+"_k"+ to_string(ik)+"_pt"+to_string(ipt)+"_eta0",inputMC[imc])->Clone();
		cout << "Corr = "<< h2dGenDetMC[imc][id][ij][ik][ipt]->GetEntries() << endl;

		//Calculate Fake rate and Miss rate
		TH1D* fakerate = (TH1D*)MC_fake[imc][id][ij][ik][ipt]->Clone();
		fakerate->Divide(fakerate, MC_Reco[imc][id][ij][ik][ipt], 1, 1, "b");
		sprintf(name,"fake_rate_d%i_j%i_k%i_pt%i_eta0", id, ij, ik, ipt);
		fakerate->SetNameTitle(name,name);

		TH1D* missrate = (TH1D*)MC_miss[imc][id][ij][ik][ipt]->Clone();
		missrate->Divide(missrate, MC_Gen[imc][id][ij][ik][ipt],1,1,"b");
		sprintf(name,"miss_rate_d%i_j%i_k%i_pt%i_eta0", id, ij, ik, ipt);
		missrate->SetNameTitle(name,name);

		fakerate->SetMinimum(-0.05); fakerate->SetMaximum(1.01);
		missrate->SetMinimum(-0.05); missrate->SetMaximum(1.01);

		MC_fakerate[imc][id][ij][ik][ipt] = (TH1D*)fakerate->Clone();
		MC_missrate[imc][id][ij][ik][ipt] = (TH1D*)missrate->Clone();

/*

		fakerate[imc][id][ij][ik][ipt] = (TH1D*)MC_fake[imc][id][ij][ik][ipt]->Clone();
                fakerate[imc][id][ij][ik][ipt]->Divide(fakerate[imc][id][ij][ik][ipt], MC_Reco[imc][id][ij][ik][ipt], 1, 1, "b");
                sprintf(name,"fakerate_d%i_j%i_k%i_pt%i_eta0", id, ij, ik, ipt);
                fakerate[imc][id][ij][ik][ipt]->SetNameTitle(name,name);

                missrate[imc][id][ij][ik][ipt] = (TH1D*)MC_miss[imc][id][ij][ik][ipt]->Clone();
                missrate[imc][id][ij][ik][ipt]->Divide(missrate[imc][id][ij][ik][ipt], MC_Gen[imc][id][ij][ik][ipt],1,1,"b");
                sprintf(name,"missrate_d%i_j%i_k%i_pt%i_eta0", id, ij, ik, ipt);
                missrate[imc][id][ij][ik][ipt]->SetNameTitle(name,name);

                fakerate[imc][id][ij][ik][ipt]->SetMinimum(-0.05); fakerate[imc][id][ij][ik][ipt]->SetMaximum(1.01);
                missrate[imc][id][ij][ik][ipt]->SetMinimum(-0.05); missrate[imc][id][ij][ik][ipt]->SetMaximum(1.01);

                MC_fakerate[imc][id][ij][ik][ipt] = (TH1D*)fakerate[imc][id][ij][ik][ipt]->Clone();
                MC_missrate[imc][id][ij][ik][ipt] = (TH1D*)missrate[imc][id][ij][ik][ipt]->Clone();

*/
		//Check RM Projection with Reco(gen) - Fake(miss) :
		TH1* RMx = h2dGenDetMC[imc][id][ij][ik][ipt]->ProjectionX(); 
		TH1* RMy = h2dGenDetMC[imc][id][ij][ik][ipt]->ProjectionY();
		
		sprintf(name,"ProjectX_d%i_j%i_k%i_pt%i_eta0", id, ij, ik, ipt); 
		RMx->SetNameTitle(name,name);
                
		sprintf(name,"ProjectY_d%i_j%i_k%i_pt%i_eta0", id, ij, ik, ipt); 
		RMy->SetNameTitle(name,name);

                sprintf(name,"Recominusfake_d%i_j%i_k%i_pt%i_eta0",id, ij, ik, ipt);
                MC_reco_fake[imc][id][ij][ik][ipt] = (TH1D*)MC_Reco[imc][id][ij][ik][ipt]->Clone(); MC_reco_fake[imc][id][ij][ik][ipt]->Reset();
                MC_reco_fake[imc][id][ij][ik][ipt]->SetNameTitle(name,name);

                sprintf(name,"Genminusmiss_d%i_j%i_k%i_pt%i_eta0",id, ij, ik, ipt);
                MC_Gen_miss[imc][id][ij][ik][ipt] = (TH1D*)MC_Gen[imc][id][ij][ik][ipt]->Clone();  MC_Gen_miss[imc][id][ij][ik][ipt]->Reset();
                MC_Gen_miss[imc][id][ij][ik][ipt]->SetNameTitle(name,name);

		for (int i = 1; i <= MC_reco_fake[imc][id][ij][ik][ipt]->GetNbinsX(); ++i) {
                double content =  MC_Reco[imc][id][ij][ik][ipt]->GetBinContent(i); 
		double factor =  MC_fake[imc][id][ij][ik][ipt]->GetBinContent(i);
                content -= factor;  MC_reco_fake[imc][id][ij][ik][ipt]->SetBinContent(i, content);
                MC_reco_fake[imc][id][ij][ik][ipt]->SetBinError(i, MC_fake[imc][id][ij][ik][ipt]->GetBinError(i)+ MC_Reco[imc][id][ij][ik][ipt]->GetBinError(i));
                }

                for (int i = 1; i <=  MC_Gen_miss[imc][id][ij][ik][ipt]->GetNbinsX(); ++i) {
                double content = MC_Gen[imc][id][ij][ik][ipt]->GetBinContent(i); 
		double factor = MC_miss[imc][id][ij][ik][ipt]->GetBinContent(i);
                content -= factor;   MC_Gen_miss[imc][id][ij][ik][ipt]->SetBinContent(i, content);
                MC_Gen_miss[imc][id][ij][ik][ipt]->SetBinError(i, MC_miss[imc][id][ij][ik][ipt]->GetBinError(i)+ MC_Gen[imc][id][ij][ik][ipt]->GetBinError(i));
                } // Bug in 2D ??
                
		MC_fakerate[imc][id][ij][ik][ipt]->Write();
		MC_missrate[imc][id][ij][ik][ipt]->Write();
		MC_reco_fake[imc][id][ij][ik][ipt]->Write();  
		MC_Gen_miss[imc][id][ij][ik][ipt]->Write();
                RMx->Write(); RMy->Write();

                //Stability and Purity
                TH1* hist_pu = (TH1D*)MC_miss[imc][id][ij][ik][ipt]->Clone(); hist_pu->Reset();
                TH1* hist_st = (TH1D*)MC_miss[imc][id][ij][ik][ipt]->Clone(); hist_st->Reset();
                TH2* RMcopy  = (TH2D*)h2dGenDetMC[imc][id][ij][ik][ipt]->Clone(); RMcopy->RebinX(2);
		for(int binRec=1; binRec<= hist_pu->GetNbinsX(); binRec++) {
                        double sum=0.;
                        for(int binGen=1; binGen<=hist_pu->GetNbinsX(); binGen++) {
                                //sum += RMcopy->GetBinContent(binGen,binRec);
                                sum += RMcopy->GetBinContent(binRec,binGen);
                                }
                        double p=0.;
                        if(sum>0.0) {
                                p = RMcopy->GetBinContent(binRec,binRec)/sum;
                                }
                                hist_pu->SetBinContent(binRec,p);
                        }
                        hist_pu->SetMinimum(-0.07); hist_pu->SetMaximum(1.01);
		//hist_pu->Write();
                //TH2* RMcopy = (TH2D*)RM_RecoGen->Clone(); RMcopy->RebinX(2);
                TH1* RMxcopy = (TH1D*)RMx->Clone(); RMxcopy->Rebin(2);
                //for(int ibin =1; ibin <= hist_pu->GetNbinsX(); ibin++){ hist_pu->SetBinContent(ibin, RMcopy->GetBinContent(ibin,ibin));};
                for(int ibin =1; ibin <= hist_pu->GetNbinsX(); ibin++){ hist_st->SetBinContent(ibin, RMcopy->GetBinContent(ibin,ibin));};
                //hist_pu->Divide(hist_pu,RMxcopy, 1, 1, "b");
                hist_st->Divide(hist_st,RMy, 1, 1, "b");

                //TH1* RMMx= (TH1D*)RMx->Clone(); RMMx->Rebin(2);
                //hist_pu->Divide(RMy,RMMx, 1, 1, "b");

		/*
                TH1D *h_pu = (TH1D*)fakerate->Clone(); h_pu->Reset();
                TH1D *h_st = (TH1D*)fakerate->Clone(); h_st->Reset();
                int ir = h_pu->GetNbinsX(); int ig = missrate->GetNbinsX();
                double fk[ir]; double ef[ig]; double pu[ir]; double st[ir];

                for(int i =1; i<h_pu->GetNbinsX()+1; i++){
                h_pu->SetBinContent(i,pu[i]);
                h_st->SetBinContent(i,st[i]);
                }
		*/

                hist_purity[imc][id][ij][ik][ipt]=(TH1D*)hist_pu->Clone();
                hist_stbl[imc][id][ij][ik][ipt]=(TH1D*)hist_st->Clone();
                
		sprintf(name,"Purity_d%i_j%i_k%i_pt%i_eta0",id, ij, ik, ipt);
		hist_purity[imc][id][ij][ik][ipt]->SetNameTitle(name,name);
		hist_purity[imc][id][ij][ik][ipt]->Write();
		
		sprintf(name,"stability_d%i_j%i_k%i_pt%i_eta0",id, ij, ik, ipt);
		hist_stbl[imc][id][ij][ik][ipt]->SetNameTitle(name,name); 
		hist_stbl[imc][id][ij][ik][ipt]->Write();

				}
			}
		}
	}
}
cout << "Histogram Read Data and MC Done "<<endl;
//--------------------------------------------------------------Read Madgraph
/*
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
////--------------------------------------------------------------HW7
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
*/
//------------------------Fold check : Patrick 1 Sep 20---------
//Get Probability Matrix (gen-miss)*probability = (reco-fake) Of course, don't forget to account for miss and fake entries (if applicable).
folddir->cd();
for(int id=0; id <ndef; id++){
     	for (int ij=0; ij<njet; ij++){
                for (int ik=0; ik<nkappa; ik++){
     			for(int ipt = 0; ipt < njetptmn ; ipt++){

      		TH2D* RM   = (TH2D*)h2dGenDetMC[umc][id][ij][ik][ipt]->Clone();
      		TH1D* Reco = (TH1D*)MC_Reco[umc][id][ij][ik][ipt]->Clone();
      		TH1D* Gen  = (TH1D*)MC_Gen[umc][id][ij][ik][ipt]->Clone();
      		TH1D* fake = (TH1D*)MC_fake[umc][id][ij][ik][ipt]->Clone();
      		TH1D* miss = (TH1D*)MC_miss[umc][id][ij][ik][ipt]->Clone();

		RM->RebinY(irbin);Gen->Rebin(irbin);miss->Rebin(irbin);

      		TH1D* Folded = (TH1D*)MC_Reco[umc][id][ij][ik][ipt]->Clone(); Folded->Reset();
		Fold(RM, Reco, Gen, miss, fake, Folded);

		sprintf(name,"Fold_d%i_j%i_k%i_pt%i_eta0", id, ij, ik, ipt);
      		Folded->SetNameTitle(name,name);
      		Folded->Write();
			}
     		}
   	}
}
//-------------Condition number of Probability Matrix-----------
for(int id=0; id <ndef; id++){
        for (int ij=0; ij<njet; ij++){
                for (int ik=0; ik<nkappa; ik++){
     			for(int ipt = 0 ; ipt < njetptmn ; ipt++){

		TH2D* RM  = (TH2D*)RMinput_RM[id][ij][ik][ipt]->Clone();
      		//TH1D* miss = (TH1D*)MC_miss2[id][ij][ik][ipt]->Clone();  // check with 2D ??
      		RM->RebinX(2); //miss->Rebin(2);
		//RM->RebinY(irbin); miss->Rebin(irbin);
		//cout << "Checking Condition Number: "<<endl;
      		cout <<setw(3) << id <<setw(3) << ij << setw(3) << ik << setw(3) << ipt <<'\n';
      		//Condition(RM, miss);
		ConditionV2(RM);
		//ConditionV3(RM);
			}
     		}
   	}
}
//--------------------------------------------------------------
/*
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
       
       		//cout << endl; 
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
//-----------------------cross check efficincy and purity calculation (Tunfold example code)
for(int id=0; id <ndef; id++){
	for (int ij=0; ij<njet; ij++){
        	for (int ik=0; ik<nkappa; ik++){
     			for(int ipt = 0 ; ipt < njetptmn ; ipt++){
     
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
//-----------------------------------------cross check efficincy and purity calculation
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
*/
//---------------------Start Unfolding---------------------------
Unfolddir->cd();
for(int id=0; id <ndef; id++){
	for (int ij=0; ij<njet; ij++){
        	for (int ik=0; ik<nkappa; ik++){
     			for(int ipt = 0 ; ipt < njetptmn ; ipt++){
		cout<<" TUnfolding : "<<" Jet Definition : "<<id<<" : Jet Number "<< ij <<" : Kappa "<<ik<<" "<<" : pt "<<ipt<<" ";
		//Rebin for match the condition of reco vs gen bin
		RMinput_RM[id][ij][ik][ipt]->RebinY(irbin);
		RMinput_Gen[id][ij][ik][ipt]->Rebin(irbin);

		// Define Inputs
		TH2* RMin = (TH2D*)RMinput_RM[id][ij][ik][ipt]->Clone();
		TH1* input = (TH1D*)Data_Reco[id][ij][ik][ipt]->Clone();
#ifdef CLOSURE
		TH1* mcgen = (TH1D*)Data_Gen[id][ij][ik][ipt]->Clone();
		TH1* mcbackground = (TH1D*)Data_fakerate[id][ij][ik][ipt]->Clone();
                TH1* mc_missrate = (TH1D*)Data_missrate[id][ij][ik][ipt]->Clone();
                TH1* mc_miss = (TH1D*)Data_miss[id][ij][ik][ipt]->Clone();
                cout << " Closure Test input " << endl;
#else
                TH1* mcgen = (TH1D*)RMinput_Gen[id][ij][ik][ipt]->Clone();
                TH1* mcbackground = (TH1D*)RMinput_fakerate[id][ij][ik][ipt]->Clone();
                TH1* mc_missrate = (TH1D*)RMinput_missrate[id][ij][ik][ipt]->Clone();
                TH1* mc_miss = (TH1D*)RMinput_miss[id][ij][ik][ipt]->Clone();
                cout << " Data Unfolding input " << endl;
#endif
                mcbackground->Multiply(input); //correction for Fake as background subtraction : Patrick  CHECK Multiplication factor ??

		//Get Reco bins
		double rxbins[RMinput_Reco[id][ij][ik][ipt]->GetNbinsX()+1]={0};
		for(int ix=0; ix<RMinput_Reco[id][ij][ik][ipt]->GetNbinsX()+1; ix++){
			rxbins[ix] = RMinput_Reco[id][ij][ik][ipt]->GetXaxis()->GetBinLowEdge(ix+1);
		}

		//Get Gen bins
		double gxbins[RMinput_Gen[id][ij][ik][ipt]->GetNbinsX()+1]={0};
		for(int ix=0; ix<RMinput_Gen[id][ij][ik][ipt]->GetNbinsX()+1; ix++){
			gxbins[ix] = RMinput_Gen[id][ij][ik][ipt]->GetXaxis()->GetBinLowEdge(ix+1);
		}

/*
		double rxbins[MC_Reco[id][ij][ik][ipt]->GetNbinsX()+1]={0};
                for (int ix=0; ix<MC_Reco[id][ij][ik][ipt]->GetNbinsX()+1; ix++) {
                        rxbins[ix] = MC_Reco[id][ij][ik][ipt]->GetXaxis()->GetBinLowEdge(ix+1);
                        }

                //Get Gen bins
                double gxbins[MC_Gen[id][ij][ik][ipt]->GetNbinsX()+1]={0};
                for (int ix=0; ix<MC_Gen[id][ij][ik][ipt]->GetNbinsX()+1; ix++) {
                        gxbins[ix] = MC_Gen[id][ij][ik][ipt]->GetXaxis()->GetBinLowEdge(ix+1);
                        }
*/
       		double biasScale = 1;
       		const char *REGULARISATION_DISTRIBUTION=0;
       		const char *REGULARISATION_AXISSTEERING="*[UOB]";

       		//https://root.cern.ch/doc/master/testUnfold5d_8C.html : this get input covariance matrix : Data covariance matrix
       		sprintf(name,"Data_covariance_d%i_j%i_k%i_pt%i_eta0", id, ij, ik, ipt);
       		TH2D* covM = new  TH2D(name,name, input->GetNbinsX(), rxbins, input->GetNbinsX(), rxbins); //?? check bins
       		covM->Sumw2();
       		for (int ix=1; ix<input->GetNbinsX()+1; ix++) {
          		double err = input->GetBinError(ix);
	  		covM->SetBinContent(ix,ix,err*err);
       			}
       			covM->Write();

       		//TUnfoldDensity class
       		//No regularisation --------------------------------
       		TUnfoldDensity tunfoldNoRegularisation(RMin, // Response matrix
					      TUnfold::kHistMapOutputVert,          // truth level on y-axis of response matrix
					      TUnfoldDensity::kRegModeNone,         // without regularisation
					      TUnfoldDensity::kEConstraintNone,     // no constrain on area
					      TUnfoldDensity::kDensityModeNone);//,0,0,REGULARISATION_DISTRIBUTION,REGULARISATION_AXISSTEERING);//,binningCoarseGen, binningFineReco);
       		//TUnfoldDensity::kDensityModeBinWidthAndUser);//,0,0,REGULARISATION_DISTRIBUTION,REGULARISATION_AXISSTEERING);//,binningCoarseGen, binningFineReco);

       		//tunfoldNoRegularisation.SubtractBackground(mcbackground, "Background", 1.0, 0.03); // hist,name,scale, scale error
		//tunfoldNoRegularisation.SubtractBackground(mcbackground, "Background", 1.0,0.04);  // hist,name,scale, scale error
		tunfoldNoRegularisation.SubtractBackground(mcbackground, "fake", 1.0, 0.00);

		//int status = tunfoldNoRegularisation.SetInput(input,biasScale,0,covarianceM);
		//int status = tunfoldNoRegularisation.SetInput(input,biasScale);  // coded in TUnfold2D code
		int status = tunfoldNoRegularisation.SetInput(input);

       		int nBadErrors = status%10000, nUnconstrOutBins = status/10000;
       		cout << nBadErrors << " bad errors and " << nUnconstrOutBins << " unconstrained output bins" << endl;

       		//Choose a value of tau to unfold with tau,0 means no regularization
       		tunfoldNoRegularisation.DoUnfold(0.0);//,input,biasScale);

       		char unfoldhist[100], title[100], NoReg_InEmatrix[100], NoReg_InEmatrix_tit[100], NoReg_RhoIJ[100], NoReg_RhoIJ_tit[100], NoReg_Ematrix[100], NoReg_Ematrix_tit[100], foldback[100], foldback_title[100], probMat[100], probmat_title[100];// NoReg_prob[100], NoReg_probtitle[100];

	       	//unfolding result, signal only
		sprintf(unfoldhist, "Tunfold_Noreg_d%i_j%i_k%i_pt%i_eta0", id, ij, ik, ipt);
		sprintf(title, "Tunfolded Noreg %i 2.5 %s %s %s", int(leadingPtThreshold[ipt]), obs_def[id], jet_num[ij], k_fact[ik] );
		TH1 *Unfolded_noRegularisation = tunfoldNoRegularisation.GetOutput(unfoldhist,title);//,0,"*[UO]" ,true);//,"","*[UO]");//,"signal");
       		//TH1 *hist_PTunfolded_noRegularisation = tunfoldNoRegularisation.GetOutput(unfoldhist,title);//,"","*[UO]");//,"signal");
       		//TH1 *hist_PTunfolded_noRegularisation = tunfoldNoRegularisation.GetOutput("hist_PTunfolded_noRegularisation", "P_{T,unfolded} [GeV]","signal");

		sprintf(NoReg_InEmatrix, "InEmat_Noreg_d%i_j%i_k%i_pt%i_eta0", id, ij, ik, ipt);
                sprintf(NoReg_InEmatrix_tit, "Input Ematrix No Regularisation %i 2.5 %s %s %s", int(leadingPtThreshold[ipt]), obs_def[id], jet_num[ij], k_fact[ik]);
		TH2 *hist_In_Emat_noRegularisation = tunfoldNoRegularisation.GetEmatrixInput(NoReg_InEmatrix, NoReg_InEmatrix_tit);

		sprintf(NoReg_RhoIJ, "Tunfold_Noreg_corr_d%i_j%i_k%i_pt%i_eta0", id, ij, ik, ipt);
		sprintf(NoReg_RhoIJ_tit, "RhoIJtotal No Regularisation %i 2.5 %s %s %s", int(leadingPtThreshold[ipt]), obs_def[id], jet_num[ij], k_fact[ik]);
       		TH2 *hist_RhoIJ_noRegularisation = tunfoldNoRegularisation.GetRhoIJtotal(NoReg_RhoIJ, NoReg_RhoIJ_tit);//,"signal");

		sprintf(NoReg_Ematrix, "Tunfold_Noreg_Emat_d%i_j%i_k%i_pt%i_eta0", id, ij, ik, ipt);
		sprintf(NoReg_Ematrix_tit, "EMatrix No Regularisation %i 2.5 %s %s %s", int(leadingPtThreshold[ipt]), obs_def[id], jet_num[ij], k_fact[ik]);
       		TH2 *hist_Emat_noRegularisation = tunfoldNoRegularisation.GetEmatrixTotal(NoReg_Ematrix, NoReg_Ematrix_tit);//,"signal");

		sprintf(foldback, "Tunfold_NoReg_Refold_d%i_j%i_k%i_pt%i_eta0", id, ij, ik, ipt);
		sprintf(foldback_title, "TunfoldFolded back  NoReg %i 2.5 %s %s %s", int(leadingPtThreshold[ipt]), obs_def[id], jet_num[ij], k_fact[ik]);
       		TH1 *hist_foldedback_NoReg = tunfoldNoRegularisation.GetFoldedOutput(foldback, foldback_title);//,"signal");

		sprintf(probMat, "Tunfold_Noreg_probM_d%i_j%i_k%i_pt%i_eta0", id, ij, ik, ipt);
		sprintf(probmat_title, "Probability matrix Noreg %i 2.5 %s %s %s", int(leadingPtThreshold[ipt]), obs_def[id], jet_num[ij], k_fact[ik]);
       		TH2 *hist_prob_Noreg = tunfoldNoRegularisation.GetProbabilityMatrix(probMat, probmat_title,TUnfold::kHistMapOutputVert);//,"signal");

		TH1 *BLTUnf = (TH1*)Unfolded_noRegularisation->Clone(); // For Bottom Line Test
       		//correction for Miss entries  : Partick 
       		for (int i = 1; i <= Unfolded_noRegularisation->GetNbinsX(); ++i) {
	 		double content = Unfolded_noRegularisation->GetBinContent(i);
	 		double factor = 1;
			//factor += (MC_miss[id][ij][ik][ipt]->GetBinContent(i)/(MC_Gen[id][ij][ik][ipt]->GetBinContent(i) - MC_miss[id][ij][ik][ipt]->GetBinContent(i)));
			factor += (mc_miss->GetBinContent(i)/(mcgen->GetBinContent(i) - mc_miss->GetBinContent(i)));
         		content *= factor;
	 		Unfolded_noRegularisation->SetBinContent(i, content);
       			}

		// Check ??
		//Unfolded_noRegularisation->Scale(1/(Unfolded_noRegularisation->Integral()));
		//mcgen->Scale(1/mcgen->Integral());

		for (int i =1; i <= Unfolded_noRegularisation->GetNbinsX(); i++){cout<<setprecision(4)<<"   "<<(Unfolded_noRegularisation->GetBinContent(i))/(MC_Gen[idata][id][ij][ik][ipt]->GetBinContent(i));}; 
		cout <<endl;

		for (int i =1; i<=Unfolded_noRegularisation->GetNbinsX(); i++){cout<<setprecision(4)<<"   "<<(Unfolded_noRegularisation->GetBinContent(i))/(mcgen->GetBinContent(i));}; 
		cout <<endl;
		
		TH1 *CTTest = (TH1*)Unfolded_noRegularisation->Clone(); CTTest->Reset();
		for (int i =1; i<=Unfolded_noRegularisation->GetNbinsX(); ++i) {CTTest->SetBinContent(i,(Unfolded_noRegularisation->GetBinContent(i)/mcgen->GetBinContent(i))); };
		cout<<endl;

		sprintf(name,"CTTest_d%i_j%i_k%i_pt%i_eta0",id,ij,ik,ipt);
		CTTest->SetNameTitle(name,name);
		CTTest->SetMinimum(0.8); CTTest->SetMaximum(1.2);
		CTTest->Write();

       		/*  //Quick check for closure Ratio Plots
       		hist_PTunfolded_noRegularisation->Scale(1/(hist_PTunfolded_noRegularisation->Integral()));
       		mcgen->Scale(1/mcgen->Integral());
       		hist_PTunfolded_noRegularisation->Divide(mcgen);
       		hist_PTunfolded_noRegularisation->SetMinimum(0.85); hist_PTunfolded_noRegularisation->SetMaximum(1.15);
       		*/

       		Unfolded_noRegularisation->Write();
		hist_In_Emat_noRegularisation->Write();
       		hist_Emat_noRegularisation->Write();
       		hist_RhoIJ_noRegularisation->Write();
       		hist_foldedback_NoReg->Write();
       		hist_prob_Noreg->Write();

#ifdef BLTest
		int ib =0;
		TH1 * BLTdata = (TH1*)Data_Reco[id][ij][ik][ipt]->Clone(); //Data Reco
                TH1 * BLTreco = (TH1*)MC_Reco[ib][id][ij][ik][ipt]->Clone(); //Data Reco
                TH1 * BLTdata1 = (TH1*)RMinput_fakerateInv[id][ij][ik][ipt]->Clone(); BLTdata1->Multiply(Data_Reco[id][ij][ik][ipt]);  //Data -background

                TH1 * BLTgen = (TH1*)MC_Gen[ib][id][ij][ik][ipt]->Clone(); //Data Reco
                TH1 * BLTunf = (TH1*)Unfolded_noRegularisation->Clone(); //Data Reco

                sprintf(name,"BLT_Data_covariance_d%i_j%i_k%i_pt%i_eta0",id, ij, ik, ipt);
                TH2D* covMBLT = RecoBin->CreateErrorMatrixHistogram(name,false); covMBLT->Sumw2();
                for (int ix=1; ix<BLTdata1->GetNbinsX()+1; ix++) {
                	double err = BLTdata1->GetBinError(ix);
                        covMBLT->SetBinContent(ix,ix,err*err);
                        }
                //https://root.cern.ch/doc/master/testUnfold5d_8C.html      : this get input covariance matrix : Data covariance matrix
                sprintf(name,"MC_covariance_d%i_j%i_k%i_pt%i_eta0",id, ij, ik,ipt);
                TH2D* covUf = GenBin->CreateErrorMatrixHistogram(name,false); covUf->Sumw2();
                for (int ix=1; ix<covUf->GetNbinsX()+1; ix++) {
                	double err = Unfolded_noRegularisation->GetBinError(ix);
                        covUf->SetBinContent(ix,ix,err*err);
                        }
        cout << " From Density Chi2A() : " << tunfoldNoRegularisation.GetChi2A() << " " << tunfoldNoRegularisation.GetChi2L()  << " " << tunfoldNoRegularisation.GetNdf() <<" Chi2/NDf : "  <<tunfoldNoRegularisation.GetChi2A()/(tunfoldNoRegularisation.GetNdf()+1) <<endl;
        
	cout << " chi2 det level 1: "; BLT(BLTdata, covM, MC_Reco[ib][id][ij][ik][ipt], 1);
        cout << " chi2 det level 2: "; BLT(BLTdata1, covM, MC_reco_fake[ib][id][ij][ik][ipt], 1);
        cout << " chi2 det level 3: "; BLT(hist_foldedback_NoReg, covMBLT, MC_reco_fake[ib][id][ij][ik][ipt], 1) ; cout << endl;

        cout << " Unfold before Miss input uncert: "; BLT(BLTUnf, hist_In_Emat_noRegularisation, MC_Gen_miss[ib][id][ij][ik][ipt]);
        cout << " Unfold before Miss Tolat: "; BLT(BLTUnf, hist_Emat_noRegularisation, MC_Gen_miss[ib][id][ij][ik][ipt]); cout << endl;

        cout << " Unfold after Miss 2: "; BLT(Unfolded_noRegularisation, hist_In_Emat_noRegularisation,MC_Gen[ib][id][ij][ik][ipt]);
        cout << " Unfold after Miss 1: "; BLT(Unfolded_noRegularisation, hist_Emat_noRegularisation, MC_Gen[ib][id][ij][ik][ipt]);

        BLTdata->Rebin(2);
        BLTreco->Rebin(2);
        Integralhist(BLTdata);
        Integralhist(BLTreco);
        
	TH1D* BLTDetRatio= (TH1D*)BLTreco->Clone(); BLTDetRatio->SetNameTitle(Form("BLTDet_d%i_j%i_k%i_pt%i_eta0",id, ij, ik, ipt),Form("BLTDet d%i j%i k%i pt%i eta0",id, ij, ik, ipt));
        BLTDetRatio->Divide(BLTdata);
        BLTDetRatio->Write();

        Integralhist(Unfolded_noRegularisation);
        Integralhist(MC_Gen[ib][id][ij][ik][ipt]);
        TH1D* BLTGenRatio= (TH1D*)MC_Gen[ib][id][ij][ik][ipt]->Clone(); BLTGenRatio->SetNameTitle(Form("BLTgen_d%i_j%i_k%i_pt%i_eta0",id,ij,ik),Form("BLTgen d%i j%i k%i pt%i eta0",id,ij,ik,ipt));
        BLTGenRatio->Divide(Unfolded_noRegularisation);
        BLTGenRatio->Write();
#endif
        //End of Variables loop
        cout << endl ;
       	//cout << "Unfolding without Regularization complete" << endl;
	     		}
   		}
	}
}
cout << "1D Unfolding without Regularization complete" << endl;
delete outputFile;
}
//-----------------------------------------------------------
//----------------------FUNCTIONS----------------------------
//-----------------------------------------------------------
void setgstyle(){
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
}
//-----------------------------------------------------------
TH1D* ReadHist1D(string name, TFile* root){
TString histname = name; cout << histname<< endl;
TH1D* hist=(TH1D*)root->Get(histname);
hist->Write();
return hist;
}
//-----------------------------------------------------------
TH2D* ReadHist2D(string name, TFile* root){
TString histname = name; cout << histname<< endl;
TH2D* hist=(TH2D*)root->Get(histname);
hist->Write();
return hist;
}
//-----------------------------------------------------------
void Integralhist(TH1 *hist){hist->Scale(1/(hist->Integral()));}
//-----------------------------------------------------------
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
                        sum = sum +(fake->GetBinContent(i)); //fake correction
                        HistoCorrect->SetBinContent(i,sum);
                        HistoCorrect->SetBinError(i,sqrt(Err));
                        }
}//end Fold
//---------------------BOTTOM LINE TEST----------------------
void BLT (TH1 * dataDistX, TH2 * dataCovX, TH1 * MCX, int rebin = 1)
{

    TH1 * dataDist = (TH1*)dataDistX->Clone();
    TH2 * dataCov = (TH2*)dataCovX->Clone();
    TH1 * MC = (TH1*)MCX->Clone();

    dataDist->Rebin(rebin);
    MC->Rebin(rebin);
    dataCov->Rebin2D(rebin, rebin);
    vector<int> indices;
    for (int i = 1; i <= dataCov->GetNbinsX(); ++i)
        if (dataCov->GetBinContent(i,i) > 0) indices.push_back(i);
    int ndf = indices.size();

    TVectorD v(ndf);
    TMatrixD M(ndf,ndf);
    for (int i = 0; i < ndf; ++i) {
        int index = indices.at(i);
        double iData = dataDist->GetBinContent(index),
               iMC   = MC      ->GetBinContent(index);
        v(i) = iData-iMC;
        for (int j = 0; j < ndf; ++j) {
            int jndex = indices.at(j);
            M(i,j) = dataCov->GetBinContent(index,jndex);
        }
    }
//    TMatrixD vT = v;
  //  vT.T();

    M.Invert();

    double chi2 = v*(M*v);

    cout << "chi2/ndf = " << chi2 << " / " << ndf << " = " << (chi2/ndf) << endl;
}
//-----------------------------------------------------------
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
//-----------------------------------------------------------
double Chi2(const TH1* hData, const TH2* covmat, const TH1* hGen, int skip/* = -1*/)
{
  //hGen->Print("all");
  //covmat->Print("all");
  //hData->Print("all");
  //throw;
  int n = hData->GetNbinsX();
  if(skip)
    n -= 1;
  TMatrixD res(1, n);
  for(int i = 0; i < n; i++)
    res(0, i) = hData->GetBinContent(i + 1) - ((hGen) ? hGen->GetBinContent(i + 1) : 0.0);
  //res.Print("all");
  //hData->Print("all");
  //hGen->Print("all");
  TMatrixD resT = res;
  resT.T();
  TMatrixD cov(n, n);
  for(int i = 0; i < n; i++)
    for(int j = 0; j < n; j++)
      cov(i, j) = covmat->GetBinContent(i + 1, j + 1);
  cov.Invert();
  double chi2 = (res * cov * resT)(0, 0);
  //hData->Print("all");
  //hGen->Print("all");
  //printf("CHI2 %.6f\n", chi2);
  return chi2;
}
//-----------------------------------------------------------
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
//-----------------Condition number calculation (Patrick code)----------
double ConditionV2(TH2 * RM)
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
//-----------------------------------------------------------Condition number calculation
double ConditionV3(TH2* RM){
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
//-----------------------------------------------------------
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
//-----------------------------------------------------------
