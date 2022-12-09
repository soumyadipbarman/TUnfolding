// Author: Stefan Schmitt
// DESY, 14.10.2008

//  Version 17.9, parallel to changes in TUnfold
//
//  History:
//    Version 17.8, parallel to changes in TUnfold
//    Version 17.7, in parallel to changes in TUnfold
//    Version 17.6, in parallel to changes in TUnfold
//    Version 17.5, in parallel to changes in TUnfold
//    Version 17.4, in parallel to changes in TUnfold
//    Version 17.3, in parallel to changes in TUnfold
//    Version 17.2, updated for writing out XML code
//    Version 17.1, in parallel to changes in TUnfold
//    Version 17.0 example for multi-dimensional unfolding
//

#include <iostream>
#include <fstream>
#include <TFile.h>
#include "TUnfoldBinningXML.h"
#include <TF2.h>

using namespace std;

/*
  This file is part of TUnfold.

  TUnfold is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  TUnfold is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with TUnfold.  If not, see <http://www.gnu.org/licenses/>.
*/

///////////////////////////////////////////////////////////////////////
// 
// Test program for the classes TUnfoldDensity and TUnfoldBinning
//
// A toy test of the TUnfold package
//
// This is an example of unfolding a two-dimensional distribution
// also using an auxillary measurement to constrain some background
//
// The example comprizes several macros
//   testUnfold5a.C   create root files with TTree objects for
//                      signal, background and data
//            -> write files  testUnfold5_signal.root
//                            testUnfold5_background.root
//                            testUnfold5_data.root
//
//   testUnfold5b.C   create a root file with the TUnfoldBinning objects
//            -> write file  testUnfold5_binning.root
//
//   testUnfold5c.C   loop over trees and fill histograms based on the
//                      TUnfoldBinning objects
//            -> read  testUnfold5_binning.root
//                     testUnfold5_signal.root
//                     testUnfold5_background.root
//                     testUnfold5_data.root
//
//            -> write testUnfold5_histograms.root
//
//   testUnfold5d.C   run the unfolding
//            -> read  testUnfold5_histograms.root
//            -> write testUnfold5_result.root
//                     testUnfold5_result.ps
// 
///////////////////////////////////////////////////////////////////////

void testUnfold5b()
{

  // write binning schemes to root file
  TFile *binningSchemes=new TFile("testUnfold5_binning_07Dec2022.root","recreate");

  // reconstructed pt, eta, discriminator
#define NBIN_PT_FINE 8
#define NBIN_ETA_FINE 10
#define NBIN_DISCR 4

  // generated pt, eta
#define NBIN_PT_COARSE 3
#define NBIN_ETA_COARSE 3

  // pt binning
  Double_t ptBinsFine[NBIN_PT_FINE+1]=
     {3.5,4.0,4.5,5.0,6.0,7.0,8.0,10.0,13.0};
  Double_t ptBinsCoarse[NBIN_PT_COARSE+1]=
     {    4.0,    5.0,    7.0,    10.0};
  // eta binning
  Double_t etaBinsFine[NBIN_ETA_FINE+1]=
     {-3.,-2.5,-2.0,-1.,-0.5,0.0,0.5,1.,2.,2.5,3.};
  Double_t etaBinsCoarse[NBIN_ETA_COARSE+1]=
     {         -2.0,    -0.5,    0.5,   2. };

  // discriminator bins
  Double_t discrBins[NBIN_DISCR+1]={0.,0.15,0.5,0.85,1.0};

  //=======================================================================
  // detector level binning scheme

  TUnfoldBinning *detectorBinning=new TUnfoldBinning("detector");
  // highest discriminator bin has fine binning
  TUnfoldBinning *detectorDistribution=
     detectorBinning->AddBinning("detectordistribution");
  detectorDistribution->AddAxis("pt",NBIN_PT_FINE,ptBinsFine,
                                false, // no underflow bin (not reconstructed)
                                true // overflow bin
                                );
  detectorDistribution->AddAxis("eta",NBIN_ETA_FINE,etaBinsFine,
                                false, // no underflow bin (not reconstructed)
                                false // no overflow bin (not reconstructed)
                                );
  detectorDistribution->AddAxis("discriminator",NBIN_DISCR,discrBins,
                                false, // no underflow bin (empty)
                                false // no overflow bin (empty)
                                );
  /* TUnfoldBinning *detectorExtra=
     detectorBinning->AddBinning("detectorextra",7,"one;zwei;three");
     detectorBinning->PrintStream(cout); */

  //=======================================================================
  // generator level binning
  TUnfoldBinning *generatorBinning=new TUnfoldBinning("generator");

  // signal distribution is measured with coarse binning
  // underflow and overflow bins are needed ot take care of
  // what happens outside the phase-space
  TUnfoldBinning *signalBinning = generatorBinning->AddBinning("signal");
  signalBinning->AddAxis("ptgen",NBIN_PT_COARSE,ptBinsCoarse,
                         true, // underflow bin
                         true // overflow bin
                         );
  signalBinning->AddAxis("etagen",NBIN_ETA_COARSE,etaBinsCoarse,
                         true, // underflow bin
                         true // overflow bin
                         );

  // this is just an example how to set bin-dependent factors
  // for the regularisation
  TF2 *userFunc=new TF2("userfunc","1./x+0.2*y^2",ptBinsCoarse[0],
                        ptBinsCoarse[NBIN_PT_COARSE],
                        etaBinsCoarse[0],etaBinsCoarse[NBIN_ETA_COARSE]);
  signalBinning->SetBinFactorFunction(1.0,userFunc);

  // background distribution is unfolded with fine binning
  // !!! in the reconstructed variable !!!
  //
  // This has the effect of "normalizing" the background in each
  // pt,eta bin to the low discriminator values
  // Only the shape of the discriminator in each (pt,eta) bin
  // is taken from Monte Carlo
  //
  // This method has been applied e.g. in
  //   H1 Collaboration, "Prompt photons in Photoproduction"
  //   Eur.Phys.J. C66 (2010) 17
  //
  TUnfoldBinning *bgrBinning = generatorBinning->AddBinning("background");
  bgrBinning->AddAxis("ptrec",NBIN_PT_FINE,ptBinsFine,
                      false, // no underflow bin (not reconstructed)
                      true // overflow bin
                      );
  bgrBinning->AddAxis("etarec",NBIN_ETA_FINE,etaBinsFine,
                      false, // no underflow bin (not reconstructed)
                      false // no overflow bin (not reconstructed)
                      );
  generatorBinning->PrintStream(cout);

  detectorBinning->Write();
  generatorBinning->Write();

  ofstream xmlOut("testUnfold5binning.xml");
  TUnfoldBinningXML::ExportXML(*detectorBinning,xmlOut,kTRUE,kFALSE);
  TUnfoldBinningXML::ExportXML(*generatorBinning,xmlOut,kFALSE,kTRUE);
  TUnfoldBinningXML::WriteDTD();
  xmlOut.close();

  const int nkappa=10;
  const int njetetamn=1;

/*
  int recojcd1bins[nkappa]={40,40,40,40,20,20,20,20,20,20};
  double recojcd1minran[nkappa]={-8.0,-5.0,-4.0,-3.0,-2.0,-2.0,-2.0,-2.0,-1.0,-1.0};
  double recojcd1maxran[nkappa]={8.0,5.0,4.0,3.0,2.0,2.0,2.0,2.0,1.0,1.0};

  int recojcd23bins[nkappa]={20,20,20,20,20,20,20,20,20,20};
  double recojcd23minran[nkappa]={-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0};
  double recojcd23maxran[nkappa]={1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};

  int genjcd1bins[nkappa]={20,20,20,20,10,10,10,10,10,10};
  double genjcd1minran[nkappa]={-8.0,-5.0,-4.0,-3.0,-2.0,-2.0,-2.0,-2.0,-1.0,-1.0};
  double genjcd1maxran[nkappa]={8.0,5.0,4.0,3.0,2.0,2.0,2.0,2.0,1.0,1.0};

  int genjcd23bins[nkappa]={10,10,10,10,10,10,10,10,10,10};
  double genjcd23minran[nkappa]={-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0};
  double genjcd23maxran[nkappa]={1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};

  //24Sept2022
  //Reco level
  int recojcd1bins[nkappa]={32,32,32,20,20,20,20,20,20,20};
  double recojcd1minran[nkappa]={-6.0,-4.0,-3.0,-2.0,-2.0,-1.0,-1.0,-1.0,-1.0,-1.0};
  double recojcd1maxran[nkappa]={6.0,4.0,3.0,2.0,2.0,1.0,1.0,1.0,1.0,1.0};

  int recojcd23bins[nkappa]={20,20,20,20,20,20,20,20,20,20};
  double recojcd23minran[nkappa]={-0.4,-0.6,-0.6,-0.8,-0.8,-1.0,-1.0,-1.0,-1.0,-1.0};
  double recojcd23maxran[nkappa]={0.4,0.6,0.6,0.8,0.8,1.0,1.0,1.0,1.0,1.0};

  //Gen level
  int genjcd1bins[nkappa]={16,16,16,10,10,10,10,10,10,10};
  double genjcd1minran[nkappa]={-6.0,-4.0,-3.0,-2.0,-2.0,-1.0,-1.0,-1.0,-1.0,-1.0};
  double genjcd1maxran[nkappa]={6.0,4.0,3.0,2.0,2.0,1.0,1.0,1.0,1.0,1.0};

  int genjcd23bins[nkappa]={10,10,10,10,10,10,10,10,10,10};
  double genjcd23minran[nkappa]={-0.4,-0.6,-0.6,-0.8,-0.8,-1.0,-1.0,-1.0,-1.0,-1.0};
  double genjcd23maxran[nkappa]={0.4,0.6,0.6,0.8,0.8,1.0,1.0,1.0,1.0,1.0};
*/

  //07Dec2022
  //Reco level
  int recojcd1bins[nkappa]={32,32,32,20,20,20,20,20,20,20};
  double recojcd1minran[nkappa]={-6.0,-4.0,-3.0,-2.0,-2.0,-1.0,-1.0,-1.0,-1.0,-1.0};
  double recojcd1maxran[nkappa]={6.0,4.0,3.0,2.0,2.0,1.0,1.0,1.0,1.0,1.0};

  int recojcd23bins[nkappa]={20,20,20,20,20,20,20,20,20,20};
  double recojcd23minran[nkappa]={-0.4,-0.6,-0.6,-0.8,-0.8,-1.0,-1.0,-1.0,-1.0,-1.0};
  double recojcd23maxran[nkappa]={0.4,0.6,0.6,0.8,0.8,1.0,1.0,1.0,1.0,1.0};

  //Gen level
  int genjcd1bins[nkappa]={16,16,16,10,10,10,10,10,10,10};
  double genjcd1minran[nkappa]={-6.0,-4.0,-3.0,-2.0,-2.0,-1.0,-1.0,-1.0,-1.0,-1.0};
  double genjcd1maxran[nkappa]={6.0,4.0,3.0,2.0,2.0,1.0,1.0,1.0,1.0,1.0};

  int genjcd23bins[nkappa]={10,10,10,10,10,10,10,10,10,10};
  double genjcd23minran[nkappa]={-0.4,-0.6,-0.6,-0.8,-0.8,-1.0,-1.0,-1.0,-1.0,-1.0};
  double genjcd23maxran[nkappa]={0.4,0.6,0.6,0.8,0.8,1.0,1.0,1.0,1.0,1.0};


  // Bin edges for Trigger Phase-space
  static const int nHLTmx=10;
  double recohtbins[nHLTmx+1] = {92, 119, 185, 251, 319, 388, 467, 518, 579, 669, 3000.0};

  // Equidistant binning for Trigger Phase-space
  int recohtnbins[nkappa]={10,10,10,10,10,10,10,10,10,10};
  double recohtbinsmin[nkappa] = {92, 119, 185, 251, 319, 388, 467, 518, 579, 669};
  double recohtbinsmax[nkappa] = {119, 185, 251, 319, 388, 467, 518, 579, 669, 3000};

  
  TUnfoldBinning *binsRec2D_D1J1[nkappa][njetetamn];
  TUnfoldBinning *binsRec2D_D1J2[nkappa][njetetamn];

  TUnfoldBinning *binsRec2D_D2J1[nkappa][njetetamn];
  TUnfoldBinning *binsRec2D_D2J2[nkappa][njetetamn];

  TUnfoldBinning *binsRec2D_D3J1[nkappa][njetetamn];
  TUnfoldBinning *binsRec2D_D3J2[nkappa][njetetamn];

  TUnfoldBinning *binsGen2D_D1J1[nkappa][njetetamn];
  TUnfoldBinning *binsGen2D_D1J2[nkappa][njetetamn];

  TUnfoldBinning *binsGen2D_D2J1[nkappa][njetetamn];
  TUnfoldBinning *binsGen2D_D2J2[nkappa][njetetamn];

  TUnfoldBinning *binsGen2D_D3J1[nkappa][njetetamn];
  TUnfoldBinning *binsGen2D_D3J2[nkappa][njetetamn];


  TUnfoldBinning *RecoBinning2D_D1J1[nkappa][njetetamn];
  TUnfoldBinning *RecoBinning2D_D1J2[nkappa][njetetamn];

  TUnfoldBinning *RecoBinning2D_D2J1[nkappa][njetetamn];
  TUnfoldBinning *RecoBinning2D_D2J2[nkappa][njetetamn];

  TUnfoldBinning *RecoBinning2D_D3J1[nkappa][njetetamn];
  TUnfoldBinning *RecoBinning2D_D3J2[nkappa][njetetamn];


  TUnfoldBinning *GenBinning2D_D1J1[nkappa][njetetamn];
  TUnfoldBinning *GenBinning2D_D1J2[nkappa][njetetamn];

  TUnfoldBinning *GenBinning2D_D2J1[nkappa][njetetamn];
  TUnfoldBinning *GenBinning2D_D2J2[nkappa][njetetamn];

  TUnfoldBinning *GenBinning2D_D3J1[nkappa][njetetamn];
  TUnfoldBinning *GenBinning2D_D3J2[nkappa][njetetamn];

  char RecoBinName[100], GenBinName[100], Axisname[100];

  for (int iet=0; iet<njetetamn; iet++) {                         //Default 1
        for (int ik=0; ik<10; ik++){
                sprintf(RecoBinName, "Detector2d_d0_j0_k%i_eta%i", ik, iet);
                binsRec2D_D1J1[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_d0_j0_k%i_eta%i",ik, iet);
                RecoBinning2D_D1J1[ik][iet]= binsRec2D_D1J1[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "d0_j0_k%i", ik);
                RecoBinning2D_D1J1[ik][iet]->AddAxis(Axisname, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik], false, false);
                sprintf(Axisname, "ht");
                //RecoBinning2D_D1J1[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
		RecoBinning2D_D1J1[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);


                sprintf(RecoBinName, "Detector2d_d0_j1_k%i_eta%i", ik, iet);
                binsRec2D_D1J2[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_d0_j1_k%i_eta%i",ik, iet);
                RecoBinning2D_D1J2[ik][iet]= binsRec2D_D1J2[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "d0_j1_k%i", ik);
                RecoBinning2D_D1J2[ik][iet]->AddAxis(Axisname, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik], false, false);
                sprintf(Axisname, "ht");
                //RecoBinning2D_D1J2[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
		RecoBinning2D_D1J2[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);


                sprintf(RecoBinName, "Detector2d_d1_j0_k%i_eta%i", ik, iet);
                binsRec2D_D2J1[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_d1_j0_k%i_eta%i",ik, iet);
                RecoBinning2D_D2J1[ik][iet]= binsRec2D_D2J1[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "d1_j0_k%i", ik);
                RecoBinning2D_D2J1[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "ht");
                //RecoBinning2D_D2J1[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
		RecoBinning2D_D2J1[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_d1_j1_k%i_eta%i", ik, iet);
                binsRec2D_D2J2[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_d1_j1_k%i_eta%i",ik, iet);
                RecoBinning2D_D2J2[ik][iet]= binsRec2D_D2J2[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "d1_j1_k%i", ik);
                RecoBinning2D_D2J2[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "ht");
                //RecoBinning2D_D2J2[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
		RecoBinning2D_D2J2[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_d2_j0_k%i_eta%i", ik, iet);
                binsRec2D_D3J1[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_d2_j0_k%i_eta%i",ik, iet);
                RecoBinning2D_D3J1[ik][iet]= binsRec2D_D3J1[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "d2_j0_k%i", ik);
                RecoBinning2D_D3J1[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "ht");
                //RecoBinning2D_D3J1[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
		RecoBinning2D_D3J1[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_d2_j1_k%i_eta%i", ik, iet);
                binsRec2D_D3J2[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_d2_j1_k%i_eta%i",ik, iet);
                RecoBinning2D_D3J2[ik][iet]= binsRec2D_D3J2[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "d2_j1_k%i", ik);
                RecoBinning2D_D3J2[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "ht");
                //RecoBinning2D_D3J2[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
		RecoBinning2D_D3J2[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

		sprintf(GenBinName, "Generator2d_d0_j0_k%i_eta%i", ik, iet);
                binsGen2D_D1J1[ik][iet] = new TUnfoldBinning(GenBinName);
                sprintf(GenBinName, "Genbin2d_d0_j0_k%i_eta%i",ik, iet);
                GenBinning2D_D1J1[ik][iet]= binsGen2D_D1J1[ik][iet]->AddBinning(GenBinName);
                sprintf(Axisname, "d0_j0_k%i", ik);
                GenBinning2D_D1J1[ik][iet]->AddAxis(Axisname, genjcd1bins[ik], genjcd1minran[ik], genjcd1maxran[ik], false, false);
                sprintf(Axisname, "ht");
                //GenBinning2D_D1J1[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
		GenBinning2D_D1J1[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(GenBinName, "Generator2d_d0_j1_k%i_eta%i", ik, iet);
                binsGen2D_D1J2[ik][iet] = new TUnfoldBinning(GenBinName);
                sprintf(GenBinName, "Genbin2d_d0_j1_k%i_eta%i",ik, iet);
                GenBinning2D_D1J2[ik][iet]= binsGen2D_D1J2[ik][iet]->AddBinning(GenBinName);
                sprintf(Axisname, "d0_j1_k%i", ik);
                GenBinning2D_D1J2[ik][iet]->AddAxis(Axisname, genjcd1bins[ik], genjcd1minran[ik], genjcd1maxran[ik], false, false);
                sprintf(Axisname, "ht");
                //GenBinning2D_D1J2[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
		GenBinning2D_D1J2[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(GenBinName, "Generator2d_d1_j0_k%i_eta%i", ik, iet);
                binsGen2D_D2J1[ik][iet] = new TUnfoldBinning(GenBinName);
                sprintf(GenBinName, "Genbin2d_d1_j0_k%i_eta%i",ik, iet);
                GenBinning2D_D2J1[ik][iet]= binsGen2D_D2J1[ik][iet]->AddBinning(GenBinName);
                sprintf(Axisname, "d1_j0_k%i", ik);
                GenBinning2D_D2J1[ik][iet]->AddAxis(Axisname, genjcd23bins[ik], genjcd23minran[ik], genjcd23maxran[ik], false, false);
                sprintf(Axisname, "ht");
                //GenBinning2D_D2J1[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
		GenBinning2D_D2J1[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(GenBinName, "Generator2d_d1_j1_k%i_eta%i", ik, iet);
                binsGen2D_D2J2[ik][iet] = new TUnfoldBinning(GenBinName);
                sprintf(GenBinName, "Genbin2d_d1_j1_k%i_eta%i",ik, iet);
                GenBinning2D_D2J2[ik][iet]= binsGen2D_D2J2[ik][iet]->AddBinning(GenBinName);
                sprintf(Axisname, "d1_j1_k%i", ik);
                GenBinning2D_D2J2[ik][iet]->AddAxis(Axisname, genjcd23bins[ik], genjcd23minran[ik], genjcd23maxran[ik], false, false);
                sprintf(Axisname, "ht");
                //GenBinning2D_D2J2[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
		GenBinning2D_D2J2[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(GenBinName, "Generator2d_d2_j0_k%i_eta%i", ik, iet);
                binsGen2D_D3J1[ik][iet] = new TUnfoldBinning(GenBinName);
                sprintf(GenBinName, "Genbin2d_d2_j0_k%i_eta%i",ik, iet);
                GenBinning2D_D3J1[ik][iet]= binsGen2D_D3J1[ik][iet]->AddBinning(GenBinName);
                sprintf(Axisname, "d2_j0_k%i", ik);
                GenBinning2D_D3J1[ik][iet]->AddAxis(Axisname, genjcd23bins[ik], genjcd23minran[ik], genjcd23maxran[ik], false, false);
                sprintf(Axisname, "ht");
                //GenBinning2D_D3J1[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
		GenBinning2D_D3J1[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(GenBinName, "Generator2d_d2_j1_k%i_eta%i", ik, iet);
                binsGen2D_D3J2[ik][iet] = new TUnfoldBinning(GenBinName);
                sprintf(GenBinName, "Genbin2d_d2_j1_k%i_eta%i",ik, iet);
                GenBinning2D_D3J2[ik][iet]= binsGen2D_D3J2[ik][iet]->AddBinning(GenBinName);
                sprintf(Axisname, "d2_j1_k%i", ik);
                GenBinning2D_D3J2[ik][iet]->AddAxis(Axisname, genjcd23bins[ik], genjcd23minran[ik], genjcd23maxran[ik], false, false);
                sprintf(Axisname, "ht");
                //GenBinning2D_D3J2[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
		GenBinning2D_D3J2[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);
                }
}

	for (int iet=0; iet<njetetamn; iet++) {
                for(int ik=0; ik<10; ik++){
                        binsRec2D_D1J1[ik][iet]->Write();
                        binsRec2D_D1J2[ik][iet]->Write();

                        binsRec2D_D2J1[ik][iet]->Write();
                        binsRec2D_D2J2[ik][iet]->Write();

                        binsRec2D_D3J1[ik][iet]->Write();
                        binsRec2D_D3J2[ik][iet]->Write();

                        binsGen2D_D1J1[ik][iet]->Write();
                        binsGen2D_D1J2[ik][iet]->Write();

                        binsGen2D_D2J1[ik][iet]->Write();
                        binsGen2D_D2J2[ik][iet]->Write();

                        binsGen2D_D3J1[ik][iet]->Write();
                        binsGen2D_D3J2[ik][iet]->Write();
                        }
                }

  delete binningSchemes;
}
