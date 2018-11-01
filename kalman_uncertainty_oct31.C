// Author: Molly Taylor
// Date: 10/02/18
// Goal: plot histograms of uncertainty in R, PHI, and Z for six INTT configurations.

// layout 0 = four layers, laddertypes 0-1-1-1
// layout 1 = four layers, laddertypes 1-1-0-1
// layout 2 = three outer layers, laddertypes 1-0-1
// layout 3 = three outer layers, laddertypes 1-1-1
// layout 4 = two outer layers, laddertypes 0-1
// layout 5 = two outer layers, laddertypes 1-1

#include <string>
#include "TFile.h"
#include "sPhenixStyle.C" //attempting to incorporate this successfully...3


// define constants common to all sets:
//const int n_layouts = 8;				// number of INTT layouts we are testing
const int n_layouts = 2;				// number of INTT layouts we are testing
//const int n_intt_layers[n_layouts] = {4,4,3,3,2,2,1,0};	// number of INTT layers in each layout
const string basepath = "~/sphenix/data/";	// base path to where data is located
//const string layout[n_layouts] = {"0111","1101","101","111","01","11","1","n"};	// convert layout# to layout configuration
const string layout[n_layouts] = {"11","01"};

	
const int n_datasets=6;
const string studypath="mom_scan";//doubles as a naming convention?
const string subpath[n_datasets]={"pi+_pt0.5GeV_phi-180-180d_z0cm_eta-1.2-1.2",
				  "pi+_pt0.6GeV_phi-180-180d_z0cm_eta-1.2-1.2",
				  "pi+_pt0.7GeV_phi-180-180d_z0cm_eta-1.2-1.2",
				  "pi+_pt0.8GeV_phi-180-180d_z0cm_eta-1.2-1.2",
				  "pi+_pt0.9GeV_phi-180-180d_z0cm_eta-1.2-1.2",
				  "pi+_pt1.0GeV_phi-180-180d_z0cm_eta-1.2-1.2"};
//				  "pi+_pt1.1GeV_phi-180-180d_z0cm_eta-1.2-1.2"};
const string setname[n_datasets]={"0.5GeV","0.6GeV","0.7GeV","0.8GeV","0.9GeV","1.0GeV"};
// */



//and define variables that it'll make our life easier to be able to access directly:
TFile *fin[n_datasets][n_layouts];
bool isZombie[n_datasets][n_layouts];
TTree *ntuple[n_datasets][n_layouts];	// pointers to trees in each of the respective data files

//generally, I'm organizing these as returned pointers, the canvas, then the stuff needed to generate the histogram, and the names/labels/limits associated with it.
void drawAndSaveSet2D(TH2F **histout, TCanvas *c, string drawcommand, string histname, string axislabels, int xbins,float xlow,float xhigh,int ybins,float ylow,float yhigh);
void drawAndFitAndSaveSet1D(TH1F **histout,TF1 **fitout, TCanvas *c, string drawcommand, string histname, string axislabels, int xbins,float xlow,float xhigh);
void drawFitParamAndSaveSet(TH1F **histout, TCanvas *c, TF1 **fit, int param, string histname, string axislabels, float minrange, float maxrange);
void drawFitErrorAndSaveSet(TH1F **histout, TCanvas *c, TF1 **fit, int param, string histname, string axislabels, float minrange, float maxrange);
void drawRMSAndSaveSet(TH1F **histout, TCanvas *c, TH1F **histin, string histname, string axislabels, float minrange, float maxrange);


void kalman_uncertainty_oct31() {
  //rcc says: don't need to do this just to read basic structures:  gSystem->Load("libg4hough.so");

  //TMacro sphstyle("sPhenixStyle.C");
  //gROOT->LoadMacro("sPhenixStyle.C+");

    //gROOT->ProcessLine(".L sPhenixStyle.C");
  
  SetsPhenixStyle();

  
	//  files, branches and their associated 'isZombie' (this could be called each time by fin[z][q]->isZombie() as well)
  	const string path=basepath+studypath+"/";



	//open all the files we need.  first index is subpath, second index is layour.
	for (int i = 0; i < n_datasets; i++){
	  for (int j = 0; j < n_layouts; j++){
	    isZombie[i][j]=false;
	    fin[i][j] = new TFile((path + subpath[i]+"/G4_sPHENIX_" + layout[j] + ".root_g4kalman_eval.root").c_str(),"READ");
	    isZombie[i][j]=fin[i][j]->IsZombie();
	    if (isZombie[i][j]) continue; //skip them if the file is broken.
	    fin[i][j]->GetObject("kalman_extrapolation_eval",ntuple[i][j]);
	  }
	}


	// leaf names to pass by reference
	float sigma_phi, sigma_z, sigma_r;
	



	//start drawing things:
	
	//create a canvas to draw all the 2D plots:
	TCanvas *c0=new TCanvas("c0","c0",800,600);
	c0->Divide(2,1);
	//	TCanvas *c0=new TCanvas("c0","c0",1000,800);
	//c0->Divide(4,2);
	gStyle->SetOptFit(111111);


	
	//looking at regions we're hitting:
	TH2F *hClustPhiZ[n_datasets][n_layouts];
	drawAndSaveSet2D((TH2F**)hClustPhiZ,c0, "z2t:phi2", "hClustPhiZ", ";phi clust (rad.);z clust(cm)", 50,-4,4,50,-100,100);

	//repeat for true clusters:
	TH2F *hTruePhiZ[n_datasets][n_layouts];
	drawAndSaveSet2D((TH2F**) hTruePhiZ,c0,"z2t:phi2t","hTruePhiZ",";phi g4 (rad.);z g4(cm)",50,-4,4,50,-100,100);

	//compare true and kalman extrapolated position
	TH2F *hTrueGuessPhiZ[n_datasets][n_layouts];
	drawAndSaveSet2D((TH2F**) hTrueGuessPhiZ,c0,"(z2te-z2t)*1e4:(phi2te-phi2t)*r2*1e4","hTrueGuessPhiZ",";phi guess-g4* (um);z guess-g4(um)",50,-4000,4000,50,-20000,20000);

	//compare true and nearest-cluster position
	TH2F *hTrueClustPhiZ[n_datasets][n_layouts];
	drawAndSaveSet2D((TH2F**) hTrueClustPhiZ,c0,"(z2c-z2t)*1e4:(phi2c-phi2t)*r2*1e4","hTrueClustPhiZ",";phi clust-g4* (um);z clust-g4(um)",50,-4000,4000,50,-20000,20000);
	
	//compare true and linear-fit position
	TH2F *hTrueLinPhiZ[n_datasets][n_layouts];
	drawAndSaveSet2D((TH2F**) hTrueLinPhiZ,c0,"(z2lin-z2t)*1e4:(phi2lin-phi2t)*r2*1e4","hTrueLinPhiZ",";phi lin-g4* (um);z lin-g4(um)",50,-2000,12000,50,-35000,35000);

	//compare phi of g4 and kalman vs number of hits:
	TH2F *hTrueNhitsPhi[n_datasets][n_layouts];
	drawAndSaveSet2D((TH2F**) hTrueNhitsPhi,c0,"nhits:(phi2te-phi2t)*r2*1e4","hTrueNhitsPhi",";phi guess-g4* (um);nhits",40,-6000,6000,10,-0.5,9.5);

	//compare phi of g4 and kalman vs pt:
	TH2F *hTruePtPhi[n_datasets][n_layouts];
	drawAndSaveSet2D((TH2F**) hTruePtPhi,c0,"pt:(phi2te-phi2t)*r2*1e4","hTruePtPhi",";phi guess-g4* (um);pt",80,-4000,4000,40,0,3);
	//compare nhits vs pt:
	TH2F *hNhitsPt[n_datasets][n_layouts];
	drawAndSaveSet2D((TH2F**) hNhitsPt,c0,"nhits:pt","hNhitsPt",";pt;nhits",50,0,2,12,-0.5,11.5);
	//compare nintt vs pt:
	TH2F *hNinttPt[n_datasets][n_layouts];
	drawAndSaveSet2D((TH2F**) hNinttPt,c0,"nintt:pt","hNinttPt",";pt;nintt",50,0,2,12,-0.5,11.5);
	//compare nmvtx vs pt:
	TH2F *hNmvtxPt[n_datasets][n_layouts];
	drawAndSaveSet2D((TH2F**) hNmvtxPt,c0,"nmvtx:pt","hNmvtxPt",";pt;nmvtx",50,0,2,7,-0.5,6.5);
	
	//compare phi of g4 and predicted position
	TH1F *hTrueGuessPhi[n_datasets][n_layouts];
	TF1 *fTrueGuessPhi[n_datasets][n_layouts];
	  drawAndFitAndSaveSet1D((TH1F**) hTrueGuessPhi,(TF1**) fTrueGuessPhi, c0,"(phi2te-phi2t)*r2t*1e4","hTrueGuessPhi",";phi guess-g4* (um)",50,-6000,6000);

	//compare z of g4 and predicted position
	TH1F *hTrueGuessZ[n_datasets][n_layouts];
	TF1 *fTrueGuessZ[n_datasets][n_layouts];
	  drawAndFitAndSaveSet1D((TH1F**) hTrueGuessZ,(TF1**) fTrueGuessZ, c0,"(z2te-z2t)*1e4","hTrueGuessZ",";z guess-g4 (um)",50,-3000,3000);
	  //was 10k for mom scan.

	//look at trends in true and predicted position
	//assumes z and phi histograms are packed the same way.

	//create a canvas to draw all the fit results plots:
	TCanvas *c1=new TCanvas("c1","c1",800,800);

	TH1F *hResolutionSigmaZ[n_layouts];
	drawFitParamAndSaveSet((TH1F**)hResolutionSigmaZ, c1, (TF1 **) fTrueGuessZ, 2, "hResolutionSigmaZ", ";Gaus. #sigma_z of g4-guess",0,5000);
	TH1F *hResolutionRMSZ[n_layouts];
	drawRMSAndSaveSet((TH1F**)hResolutionRMSZ, c1, (TH1F **) hTrueGuessZ, "hResolutionRMSZ", ";RMS #sigma_z of g4-guess",0,5000);
	TH1F *hResolutionMeanZ[n_layouts];
	drawFitParamAndSaveSet((TH1F**)hResolutionMeanZ, c1, (TF1 **) fTrueGuessZ, 1, "hResolutionMeanZ", ";Mean of #sigma_z of g4-guess",-500,500);
	TH1F *hResolutionMeanErrZ[n_layouts];
	drawFitErrorAndSaveSet((TH1F**)hResolutionMeanErrZ, c1, (TF1 **) fTrueGuessZ, 1, "hResolutionMeanErrZ", ";Err on mean #sigma_z of g4-guess",0,500);

	TH1F *hResolutionSigmaPhi[n_layouts];
	drawFitParamAndSaveSet((TH1F**)hResolutionSigmaPhi, c1, (TF1 **) fTrueGuessPhi, 2, "hResolutionSigmaPhi", ";Gaus. #sigma_ #phi of g4-guess",0,5000);
	TH1F *hResolutionRMSPhi[n_layouts];
	drawRMSAndSaveSet((TH1F**)hResolutionRMSPhi, c1, (TH1F **) hTrueGuessPhi, "hResolutionRMSPhi", ";RMS #sigma_ #phi of g4-guess",0,5000);
	TH1F *hResolutionMeanPhi[n_layouts];
	drawFitParamAndSaveSet((TH1F**)hResolutionMeanPhi, c1, (TF1 **) fTrueGuessPhi, 1, "hResolutionMeanPhi", ";Mean of #sigma_ #phi of g4-guess",-500,500);	TH1F *hResolutionMeanErrPhi[n_layouts];
	drawFitErrorAndSaveSet((TH1F**)hResolutionMeanErrPhi, c1, (TF1 **) fTrueGuessPhi, 1, "hResolutionMeanErrPhi", ";Err on mean #sigma_ #phi of g4-guess",0,500);

	return;
	
}

void drawAndSaveSet2D(TH2F **histout, TCanvas *c, string drawcommand, string histname, string axislabels, int xbins,float xlow,float xhigh,int ybins,float ylow,float yhigh)
{
  //assumes c has enough divisions to cover all of these.
  int nsets=n_datasets;
  int nsubs=n_layouts;
  for (int i = 0; i < nsets; i++){
    for (int j = 0; j < nsubs; j++){
      histout[i*nsubs+j]= new TH2F((histname + setname[i]+ layout[j]).c_str(),(layout[j]+axislabels).c_str(),xbins,xlow,xhigh,ybins,ylow,yhigh);

      c->cd(j+1);
      if (!isZombie[i][j])
	ntuple[i][j]->Draw((drawcommand+">>"+histname + setname[i]+ layout[j]).c_str());
      histout[i*nsubs+j]->Draw("colz");
    }
    c->SaveAs((histname+"_"+studypath+"_"+setname[i]+".pdf").c_str());

  }
  return;
}
void drawAndFitAndSaveSet1D(TH1F **histout,TF1 **fitout, TCanvas *c, string drawcommand, string histname, string axislabels, int xbins,float xlow,float xhigh)
{
  //assumes c has just enough divisions to cover all of these.
  int nsets=n_datasets;
  int nsubs=n_layouts;
  TH1F* ht;
  TF1* ft;
  for (int i = 0; i < nsets; i++){
    for (int j = 0; j < nsubs; j++){
      c->cd(j+1);
      ht=histout[i*nsubs+j]= new TH1F((histname + setname[i]+ layout[j]).c_str(),(layout[j]+axislabels).c_str(),xbins,xlow,xhigh);
      ft=fitout[i*nsubs+j]=new TF1(("fit"+histname + setname[i]+ layout[j]).c_str(),"gaus(0)",xlow,xhigh);

      ft->SetParNames("scale","x0","sigma");
      ft->SetParameters(ht->GetEntries(),ht->GetMean(),ht->GetRMS());
      if (!isZombie[i][j])
	ntuple[i][j]->Draw((drawcommand+">>"+histname + setname[i]+ layout[j]).c_str());

      histout[i*nsubs+j]->Draw();
      ht->Fit(ft);
    }
    c->SaveAs((histname+"_"+studypath+"_"+setname[i]+".pdf").c_str());

  }
  return;
}

void drawFitParamAndSaveSet(TH1F **histout, TCanvas *c, TF1 **fit, int param, string histname, string axislabels, float minrange, float maxrange)
{
  TLegend *leg;
  c->cd();
	
  leg = new TLegend(0.20,0.75,0.480,0.9);
  leg->SetNColumns(2);
    leg->AddEntry("","#it{#bf{sPHENIX}} Preliminary","");

    for (int i=0;i<n_layouts;i++){
    histout[i]=new TH1F((histname+ layout[i]).c_str(),axislabels.c_str(),n_datasets,-0.5,n_datasets-0.5);
    histout[i]->SetLineColor(i+1);
    histout[i]->SetLineWidth(2);
    for (int j=0;j<n_datasets;j++){
      int bin=histout[i]->Fill(setname[j].c_str(),fit[j*n_layouts+i]->GetParameter(param));
      histout[i]->SetBinError(bin,fit[j*n_layouts+i]->GetParError(param));
    }

    if (i==0){
      histout[i]->GetYaxis()->SetRangeUser(minrange,maxrange);
      histout[i]->Draw();
    }else {
      histout[i]->Draw("SAME");
    }
  }
  for (int i=0;i<n_layouts;i++){
    leg->AddEntry(histout[i],layout[i].c_str(),"l");
  }
  leg->Draw();
	  
  c->SaveAs((histname+studypath+".pdf").c_str());
  return;
}

void drawFitErrorAndSaveSet(TH1F **histout, TCanvas *c, TF1 **fit, int param, string histname, string axislabels, float minrange, float maxrange)
{
  TLegend *leg;
  c->cd();
	
  leg = new TLegend(0.20,0.75,0.480,0.9);
  leg->SetNColumns(2);
    leg->AddEntry("","#it{#bf{sPHENIX}} Preliminary","");

  for (int i=0;i<n_layouts;i++){
    histout[i]=new TH1F((histname+ layout[i]).c_str(),axislabels.c_str(),n_datasets,-0.5,n_datasets-0.5);
    histout[i]->SetLineColor(i+1);
    histout[i]->SetLineWidth(2);
    for (int j=0;j<n_datasets;j++){
      int bin=histout[i]->Fill(setname[j].c_str(),fit[j*n_layouts+i]->GetParError(param));
      //histout[i]->SetBinError(bin,fit[j*n_layouts+i]->GetParError(param));
    }

    if (i==0){
      histout[i]->GetYaxis()->SetRangeUser(minrange,maxrange);
      histout[i]->Draw();
    }else {
      histout[i]->Draw("SAME");
    }
  }
  for (int i=0;i<n_layouts;i++){
    leg->AddEntry(histout[i],layout[i].c_str(),"l");
  }
  leg->Draw();
	  
  c->SaveAs((histname+studypath+".pdf").c_str());
  return;
}

void drawRMSAndSaveSet(TH1F **histout, TCanvas *c, TH1F **histin, string histname, string axislabels, float minrange, float maxrange)
{
  TLegend *leg;
  c->cd();
	
  leg = new TLegend(0.20,0.15,0.480,0.3);
  leg->SetNColumns(2);
    leg->AddEntry("","#it{#bf{sPHENIX}} Preliminary","");

  for (int i=0;i<n_layouts;i++){
    histout[i]=new TH1F((histname+ layout[i]).c_str(),axislabels.c_str(),n_datasets,-0.5,n_datasets-0.5);
    histout[i]->SetLineColor(i+1);
    histout[i]->SetLineWidth(2);
    for (int j=0;j<n_datasets;j++){
      int bin=histout[i]->Fill(setname[j].c_str(),histin[j*n_layouts+i]->GetRMS());
      histout[i]->SetBinError(bin,histin[j*n_layouts+i]->GetRMSError());
    }

    if (i==0){
      histout[i]->GetYaxis()->SetRangeUser(minrange,maxrange);
      histout[i]->Draw();
    }else {
      histout[i]->Draw("SAME");
    }
  }
  for (int i=0;i<n_layouts;i++){
    leg->AddEntry(histout[i],layout[i].c_str(),"l");
  }
  leg->Draw();
	  
  c->SaveAs((histname+studypath+".pdf").c_str());
  return;
}

















