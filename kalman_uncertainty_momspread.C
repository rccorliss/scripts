// Author: Molly Taylor and Ross Corliss
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

//beautification details:
const string sPHENIXbanner="#it{#bf{sPHENIX}} Simulation";

// define constants common to all sets:
const string basepath = "~/sphenix/data/auto/";
//const string basepath = "/gpfs/mnt/gpfs04/sphenix/user/mitay/data/";

const string outpath="112818/";
const string particlename="pi+";

#define set2

#ifdef set1
const string studypath="mom_spread_112418";//doubles as a naming convention?
const string outname="spread";
const int n_layouts = 8;
const string particle="pi+";
const string layout[n_layouts]={"000p","00pp","00zp","0ppp","0pzp","ppzp","zppp","00pp_out"};
const string setname[n_layouts]={"000p","00pp","00zp","0ppp","0pzp","ppzp","zppp","00pp_out"};
const int color[n_layouts]={9,4,7,8,5,42,2,6};
const string fileprefix="";
const string filemiddle="_pt0.5-3.0GeV_";
const string filesuffix=".root_g4kalman_eval.root";
#endif

#ifdef set2
const string studypath="mom_scan_112418";//doubles as a naming convention?
const string outname="discrete";
const int n_layouts = 8;
const string particle="pi+";
const string layout[n_layouts]={"000p","00pp","00zp","0ppp","0pzp","ppzp","zppp","00pp_out"};
const string setname[n_layouts]={"000p","00pp","00zp","0ppp","0pzp","ppzp","zppp","00pp_out"};
const int color[n_layouts]={9,4,7,8,5,42,2,6};
const string fileprefix="";
const string filemiddle="_pt_discrete_";
const string filesuffix=".root_g4kalman_eval.root";
#endif
#ifdef set2_112118
const string studypath="mom_spread_112118";//doubles as a naming convention?
const string outname="all_7";
const int n_layouts = 5;
const string particle="pi+";
const string layout[n_layouts]={"000p","00pp","0ppp","zppp","00pp_out"};
const string setname[n_layouts]={"000p","00pp","0ppp","zppp","00pp_out"};
const int color[n_layouts]={9,4,7,8,5};
const string fileprefix="";
const string filemiddle="_pt0.5-3.0GeV_";
const string filesuffix=".root_g4kalman_eval.root";
#endif

#ifdef set3_112118
const string studypath="mom_scan_112118";//doubles as a naming convention?
const string outname="all_7";
const int n_layouts = 4;
const string particle="pi+";
const string layout[n_layouts]={"0ppp","0pzp","zppp","ppzp"};
const string setname[n_layouts]={"0ppp","0pzp","zppp","ppzp"};
const int color[n_layouts]={7,1,2,8};
const string fileprefix="";
const string filemiddle="_pt_discrete_";
const string filesuffix=".root_g4kalman_eval.root";
#endif


//and define variables that it'll make our life easier to be able to access directly:
int n_sample_bins=26;
const int n_datasets=26;
float sample_min;
float sample_max;
string sample_name;
TFile *fin[n_layouts];
bool isZombie[n_layouts];
TTree *ntuple[n_layouts];	// pointers to trees in each of the respective data files


//generally, I'm organizing these as returned pointers, the canvas, then the stuff needed to generate the histogram, and the names/labels/limits associated with it.
void drawAndSaveSet2D(TH2F **histout, TCanvas *c, string drawcommand, string histname, string axislabels, int xbins,float xlow,float xhigh,int ybins,float ylow,float yhigh);
void drawAndSaveSet2D(TH2F **histout, TCanvas *c, string drawcommand, string weightcommand, string histname, string axislabels, int xbins,float xlow,float xhigh,int ybins,float ylow,float yhigh);
void drawAndSaveSet1D(TH1F **histout, TCanvas *c, string drawcommand, string histname, string axislabels, int xbins,float xlow,float xhigh);
void drawAndFitAndSaveSet1D(TH1F **histout,TF1 **fitout, TCanvas *c, string drawcommand, string histname, string axislabels, int xbins,float xlow,float xhigh);
void drawAndFitAndSaveSet1D(TH1F **histout,TF1 **fitout, TCanvas *c, string drawcommand, string weightcommand, string histname, string axislabels, int xbins,float xlow,float xhigh); //with weighting so we can make a cut.
void drawFitParamAndSaveSet(TH1F **histout, TCanvas *c, TF1 **fit, int param, string histname, string axislabels, float minrange, float maxrange);
void drawFitErrorAndSaveSet(TH1F **histout, TCanvas *c, TF1 **fit, int param, string histname, string axislabels, float minrange, float maxrange);
void drawRMSAndSaveSet(TH1F **histout, TCanvas *c, TH1F **histin, string histname, string axislabels, float minrange, float maxrange);


void kalman_uncertainty_momspread() {
//rcc says: don't need to do this just to read basic structures:  gSystem->Load("libg4hough.so");

  SetsPhenixStyle();

  const string path=basepath+studypath+"/";

  //open all the files we need.  first index is subpath, second index is layout.
  for (int i = 0; i < n_layouts; i++){
    ntuple[i]=nullptr;
    isZombie[i]=false;
    fin[i] = new TFile((path +fileprefix+particle+filemiddle+layout[i]+filesuffix).c_str(),"READ");
    isZombie[i]=fin[i]->IsZombie();
    if (isZombie[i]) continue; //skip them if the file is broken.
    fin[i]->GetObject("kalman_eval",ntuple[i]);
  }


	//start drawing things:
	
	//create a canvas to draw all the 2D plots:
	//	TCanvas *c0=new TCanvas("c0","c0",800,600);
	//c0->Divide(2,1);
	TCanvas *c0;


	if (n_layouts==2 ){
	  c0=new TCanvas("c0","c0",1000,800);
	  c0->Divide(2,1);
	}
	else if (n_layouts==3 ){
	  c0=new TCanvas("c0","c0",1500,800);
	  c0->Divide(3,1);
	}
	else{
	  c0=new TCanvas("c0","c0",1000,800);
	  c0->Divide(4,2);
	}
//	gStyle->SetOptFit(111111);



n_sample_bins=26;
 sample_min=0.45;
 sample_max=3.05;
 sample_name="true_pti";	

	TCanvas *c2=new TCanvas("c2","c2",800,800);

	TH1F *hResolutionPhiCalc[n_layouts];
	drawDataAndSaveSet((TH1F**)hResolutionPhiCalc, c2, "hResPhiCalc",";pT [GeV];#phi res. at 30cm [mm]",0,2);

	return;
	
	//looking at regions we're hitting:
	TH2F *hTruePhiZ[n_datasets][n_layouts];
	drawAndSaveSet2D((TH2F**) hTruePhiZ,c0,"z30t:phi30t","ok30t","hTruePhiZ",";phi g4 [rad.];z g4[cm]",50,-4,4,50,-100,100);

	//compare true and kalman extrapolated position
	TH2F *hTrueGuessPhiZ[n_datasets][n_layouts];
	drawAndSaveSet2D((TH2F**) hTrueGuessPhiZ,c0,"(z30te-z30t)*10:(phi30te-phi30t)*r30t*10","ok30t&&ok30te","hDelta2D",";phi guess-g4* [mm];z guess-g4[mm]",50,-3,3,50,-15,15);

		//compare true and kalman extrapolated position

       	TH2F *hG4vsClusterTrack2D[n_datasets][n_layouts];
	drawAndSaveSet2D((TH2F**) hG4vsClusterTrack2D,c0,"(z30te-g4_z30te)*10:(phi30te-g4_phi30te)*r30te*10","ok30te&&g4_ok30te","hG4vsClusterTrack2D",";phi cltrack-g4track [mm];z cltrack-g4track[mm]",50,-3,3,50,-15,15);

	TH2F *hG4vsTruth2D[n_datasets][n_layouts];
	drawAndSaveSet2D((TH2F**) hG4vsClusterTrack2D,c0,"(g4_z30t-g4_z30te)*10:(g4_phi30t-g4_phi30te)*g4_r30t*10","g4_ok30t&&g4_ok30te","hG4vsTruth2D",";phi truth-g4track [mm];z truth-g4track[mm]",50,-3,3,50,-15,15);

	//check ng4hits
	TH1F *hG4Hits[n_datasets][n_layouts];
	drawAndSaveSet1D((TH1F**) hG4Hits, c0,"ng430t","hNG4Hits",";ng4hits",3,-0.5,2.5);

	//compare phi of g4 and predicted position
	TH1F *hTrueGuessPhi[n_datasets][n_layouts];
	TF1 *fTrueGuessPhi[n_datasets][n_layouts];
	  drawAndFitAndSaveSet1D((TH1F**) hTrueGuessPhi,(TF1**) fTrueGuessPhi, c0,"(phi30te-phi30t)*r30t*10","ok30t&&ok30te","hDeltaPhi",";phi kalman-g4* [mm]",80,-8,8);

	  //compare z of g4 and predicted position
	TH1F *hTrueGuessZ[n_datasets][n_layouts];
	TF1 *fTrueGuessZ[n_datasets][n_layouts];
	  drawAndFitAndSaveSet1D((TH1F**) hTrueGuessZ,(TF1**) fTrueGuessZ, c0,"(z30te-z30t)*10","ok30t&&ok30te","hDeltaZ",";z guess-g4 [mm]",50,-4,4);


	//the next two don't work, but the previous two do...?  g4_ok30t is not setting correctly.
	  /*
        //compare phi of g4 and predicted position
	TH1F *hTrueGuessPhig4[n_datasets][n_layouts];
	TF1 *fTrueGuessPhig4[n_datasets][n_layouts];
	  drawAndFitAndSaveSet1D((TH1F**) hTrueGuessPhig4,(TF1**) fTrueGuessPhig4, c0,"(g4_phi30te-g4_phi30t)*g4_r30t*10","g4_ok30t&&g4_ok30te","hDeltaPhiG4",";phi g4kalman-g4truth [mm]",80,-8,8);

	  //compare z of g4 and predicted position
	TH1F *hTrueGuessZg4[n_datasets][n_layouts];
	TF1 *fTrueGuessZg4[n_datasets][n_layouts];
	  drawAndFitAndSaveSet1D((TH1F**) hTrueGuessZg4,(TF1**) fTrueGuessZg4, c0,"(g4_z30te-g4_z30t)*10","g4_ok30t&&g4_ok30te","hDeltaZG4",";z g4kalman-g4truth [mm]",50,-3,3);
	  */
	  //was 10k for mom scan.


	TH1F *hTrueGuessZtemp[n_datasets][n_layouts];
	TF1 *fTrueGuessZtemp[n_datasets][n_layouts];

	//for (float ptlow=0;ptlow<ptmax;ptlow+=ptstep){
	//drawAndFitAndSaveSet1D((TH1F**) hTrueGuessZg4,(TF1**) fTrueGuessZg4, c0,"(g4_z30te-g4_z30t)*10",Form("true_pti>%f && true_pti<%f",ptlow,pthigh),"hDeltaZtemp",";z g4kalman-g4truth [mm]",50,-3,3);
	//}


	  

	//look at trends in true and predicted position
	//assumes z and phi histograms are packed the same way.

	//create a canvas to draw all the fit results plots:
	TCanvas *c1=new TCanvas("c1","c1",800,800);

	TH1F *hResolutionSigmaZ[n_layouts];
	drawFitParamAndSaveSet((TH1F**)hResolutionSigmaZ, c1, (TF1 **) fTrueGuessZ, 2, "hResZ", Form(";%s [GeV];z res. at 30cm [mm]",sample_name.c_str()),0,2);


	TH1F *hResolutionSigmaPhi[n_layouts];
	drawFitParamAndSaveSet((TH1F**)hResolutionSigmaPhi, c1, (TF1 **) fTrueGuessPhi, 2, "hResPhi",Form(";%s [GeV];#phi res. at 30cm [mm]",sample_name.c_str()),0,2);

	/*
	TH1F *hResolutionSigmaZg4[n_layouts];

	drawFitParamAndSaveSet((TH1F**)hResolutionSigmaZg4, c1, (TF1 **) fTrueGuessZg4, 2, "hResZg4",Form(";%s [GeV];z res. at 30cm w/ g4track[mm]",sample_name.c_str()),0,2);


	TH1F *hResolutionSigmaPhig4[n_layouts];
	drawFitParamAndSaveSet((TH1F**)hResolutionSigmaPhig4, c1, (TF1 **) fTrueGuessPhig4, 2, "hResPhig4",Form(";%s [GeV];#phi res. at 30cm w/ g4track[mm]",sample_name.c_str()),0,2);
	*/


	return; //don't bother doing the eta-dependent fits.  they just eat up time.


	//task list:
	// - compare gaussian widths at 30cm for 00pp and zppp as a function of pT true
	//do this by defining only those samples.

	// for later:  - compare gaussian widths at 30cm for 00pp and zppp as a function of pT reco, summed across slices

		//compare phi of g4 and predicted position
	TH1F *htempTrueGuessPhi[n_datasets][n_layouts];
	TF1 *ftempTrueGuessPhi[n_datasets][n_layouts];
	TH1F *htempResolution[n_layouts];
	float etabinmidpoint=0.0;
	float etabinwidth=0.8;
	drawAndFitAndSaveSet1D((TH1F**) htempTrueGuessPhi,(TF1**) ftempTrueGuessPhi, c0,"(phi30te-phi30t)*r30t*10",
			       Form("abs(abs(0.5*log( (sqrt(pt*pt+pz*pz)+pz)/(sqrt(pt*pt+pz*pz)-pz))) -%f)<=%f",etabinmidpoint,etabinwidth/2.0),
			       "hDeltaPhi_eta_0",";#Delta #phi [um]",50,-6,6);
	drawFitParamAndSaveSet((TH1F**)htempResolution, c1, (TF1 **) ftempTrueGuessPhi, 2, "hResolutionSigmaPhi_central", ";pT [GeV];#phi res. at 30cm [mm]",0,2);
	etabinwidth=0.4;//since we have two separate regions that fill this now.
	etabinmidpoint=0.6;
	drawAndFitAndSaveSet1D((TH1F**) htempTrueGuessPhi,(TF1**) ftempTrueGuessPhi, c0,"(phi30te-phi30t)*r30t*10",
			       Form("abs(abs(0.5*log( (sqrt(pt*pt+pz*pz)+pz)/(sqrt(pt*pt+pz*pz)-pz))) -%f)<=%f",etabinmidpoint,etabinwidth/2.0),
			       "hDeltaPhi_eta_0.6",";#Delta #phi [um]",50,-6,6);
	drawFitParamAndSaveSet((TH1F**)htempResolution, c1, (TF1 **) ftempTrueGuessPhi, 2, "hResolutionSigmaPhi_mid", ";pT [GeV];#phi res. at 30cm [mm]",0,2);
	etabinmidpoint=1.0;
	drawAndFitAndSaveSet1D((TH1F**) htempTrueGuessPhi,(TF1**) ftempTrueGuessPhi, c0,"(phi30te-phi30t)*r30t*10",
			       Form("abs(abs(0.5*log( (sqrt(pt*pt+pz*pz)+pz)/(sqrt(pt*pt+pz*pz)-pz))) -%f)<=%f",etabinmidpoint,etabinwidth/2.0),
			       "hDeltaPhi_eta_1.0",";#Delta #phi [um]",50,-6,6);
	drawFitParamAndSaveSet((TH1F**)htempResolution, c1, (TF1 **) ftempTrueGuessPhi, 2, "hResolutionSigmaPhi_high", ";pT [GeV];#phi res. at 30cm [mm]",0,2);
	//compare z of g4 and predicted position
	TH1F *htempTrueGuessZ[n_datasets][n_layouts];
	TF1 *ftempTrueGuessZ[n_datasets][n_layouts];
	//	TH1F *htempResolution[n_layouts];
	etabinwidth=0.8;
	etabinmidpoint=0.0;
	drawAndFitAndSaveSet1D((TH1F**) htempTrueGuessZ,(TF1**) ftempTrueGuessZ, c0,"(z30te-z30t)*10",
			       Form("abs(abs(0.5*log( (sqrt(pt*pt+pz*pz)+pz)/(sqrt(pt*pt+pz*pz)-pz))) -%f)<=%f",etabinmidpoint,etabinwidth/2.0),
			       "hDeltaZ_eta_0",";#Delta #phi [um]",50,-6,6);
	drawFitParamAndSaveSet((TH1F**)htempResolution, c1, (TF1 **) ftempTrueGuessZ, 2, "hResolutionSigmaZ_central", ";pT [GeV];z res. at 30cm [mm]",0,2);
	etabinwidth=0.4;//since we have two separate regions that fill this now.
	etabinmidpoint=0.4;
	drawAndFitAndSaveSet1D((TH1F**) htempTrueGuessZ,(TF1**) ftempTrueGuessZ, c0,"(z30te-z30t)*10",
			       Form("abs(abs(0.5*log( (sqrt(pt*pt+pz*pz)+pz)/(sqrt(pt*pt+pz*pz)-pz))) -%f)<=%f",etabinmidpoint,etabinwidth/2.0),
			       "hDeltaZ_eta_0.6",";#Delta #phi [um]",50,-6,6);
	drawFitParamAndSaveSet((TH1F**)htempResolution, c1, (TF1 **) ftempTrueGuessZ, 2, "hResolutionSigmaZ_mid", ";pT [GeV];z res. at 30cm [mm]",0,2);
	etabinmidpoint=0.8;
	drawAndFitAndSaveSet1D((TH1F**) htempTrueGuessZ,(TF1**) ftempTrueGuessZ, c0,"(z30te-z30t)*10",
			       Form("abs(abs(0.5*log( (sqrt(pt*pt+pz*pz)+pz)/(sqrt(pt*pt+pz*pz)-pz))) -%f)<=%f",etabinmidpoint,etabinwidth/2.0),
			       "hDeltaZ_eta_1.0",";#Delta #phi [um]",50,-6,6);
	drawFitParamAndSaveSet((TH1F**)htempResolution, c1, (TF1 **) ftempTrueGuessZ, 2, "hResolutionSigmaZ_high", ";pT [GeV];z res. at 30cm [mm]",0,2);
	return;
	
	// - do the above but cut into three eta bins: eta<0.4, eta<0.8, eta<1.2 for eta of hit (or polar angle?  or ...?)
	//something like:  drawAndSaveSet2D(hist
	// - repeat all that for the set of zppp and ppzp

	//presentation plots:
	
	return;
	
}

void drawAndSaveSet2D(TH2F **histout, TCanvas *c, string drawcommand, string histname, string axislabels, int xbins,float xlow,float xhigh,int ybins,float ylow,float yhigh)
{
  drawAndSaveSet2D(histout,c,drawcommand, string("1"),histname,axislabels,xbins,xlow,xhigh,ybins,ylow,yhigh);
    return;
}
void drawAndSaveSet2D(TH2F **histout, TCanvas *c, string drawcommand, string weightcommand, string histname, string axislabels, int xbins,float xlow,float xhigh,int ybins,float ylow,float yhigh)
{
  //assumes c has enough divisions to cover all of these.
  //now that we have ptspread, we have to generate one hist for each 'set':
  
  int nsets=n_sample_bins;
  int nsubs=n_layouts;
  float setwidth=(sample_max-sample_min)/n_sample_bins;
  for (int i = 0; i < nsets; i++){
    float setlow=setwidth*i+sample_min;
    float sethigh=setlow+setwidth;
    for (int j = 0; j < nsubs; j++){
      histout[i*nsubs+j]= new TH2F((string(histname) + string(to_string(i))+ string(layout[j])).c_str(),(string(layout[j])+string(axislabels)).c_str(),xbins,xlow,xhigh,ybins,ylow,yhigh);

      c->cd(j+1);
      if (!isZombie[j])
	ntuple[j]->Draw((drawcommand+">>"+histname + to_string(i)+ layout[j]).c_str(),Form("(%s)*(%f<%s && %s<%f)",weightcommand.c_str(),setlow,sample_name.c_str(),sample_name.c_str(),sethigh));
      histout[i*nsubs+j]->Draw("colz");
    }
    c->SaveAs((outpath+histname+"_"+outname+"_"+to_string(i)+".pdf").c_str());

  }
  return;
}

void drawAndSaveSet1D(TH1F **histout, TCanvas *c, string drawcommand, string histname, string axislabels, int xbins,float xlow,float xhigh)
{
  //assumes c has just enough divisions to cover all of these.
 int nsets=n_sample_bins;
  int nsubs=n_layouts;
  float setwidth=(sample_max-sample_min)/n_sample_bins;

  TH1F* ht;
  for (int i = 0; i < nsets; i++){
    float setlow=setwidth*i+sample_min;
    float sethigh=setlow+setwidth;
    for (int j = 0; j < nsubs; j++){
      c->cd(j+1);
      ht=histout[i*nsubs+j]= new TH1F((string(histname) + string(to_string(i))+ string(layout[j])).c_str(),(string(layout[j])+string(axislabels)).c_str(),xbins,xlow,xhigh);
      if (!isZombie[j])
	ntuple[j]->Draw((drawcommand+">>"+histname + to_string(i)+ layout[j]).c_str(),Form("%f<%s && %s<%f",setlow,sample_name.c_str(),sample_name.c_str(),sethigh));

      histout[i*nsubs+j]->Draw();
    }
    c->SaveAs((outpath+histname+"_"+outname+"_"+to_string(i)+".pdf").c_str());

  }
  return;
}


void drawAndFitAndSaveSet1D(TH1F **histout,TF1 **fitout, TCanvas *c, string drawcommand, string histname, string axislabels, int xbins,float xlow,float xhigh)
{
  drawAndFitAndSaveSet1D(histout,fitout,c,drawcommand, string("1"),histname,axislabels,xbins,xlow,xhigh);
  /* just route it to the weight version with weight of one:
  //assumes c has just enough divisions to cover all of these.
 int nsets=n_sample_bins;
  int nsubs=n_layouts;
  float setwidth=(sample_max-sample_min)/n_sample_bins;
  TH1F* ht;
  TF1* ft;
  for (int i = 0; i < nsets; i++){
    for (int j = 0; j < nsubs; j++){
      c->cd(j+1);
      ht=histout[i*nsubs+j]= new TH1F((string(histname) + to_string(i)+ string(layout[j])).c_str(),(string(layout[j])+string(axislabels)).c_str(),xbins,xlow,xhigh);
      ft=fitout[i*nsubs+j]=new TF1(("fit"+histname + to_string(i)+ layout[j]).c_str(),"gaus(0)",xlow,xhigh);


      if (!isZombie[j])
	ntuple[j]->Draw((drawcommand+">>"+histname + to_string(i)+ layout[j]).c_str(),Form("%f<%s && %s<%f",setwidth*i+sample_min,sample_name.c_str(),sample_name.c_str(),setwidth*(i+1)+sample_min));
      ft->SetParNames("scale","x0","sigma");
      ft->SetParameters(ht->GetEntries(),ht->GetMean(),ht->GetRMS());
      histout[i*nsubs+j]->Draw();
      ht->Fit(ft);
    }
    c->SaveAs((outpath+histname+"_"+outname+"_"+to_string(i)+".pdf").c_str());

  }
  */
  return;
}

void drawAndFitAndSaveSet1D(TH1F **histout,TF1 **fitout, TCanvas *c, string drawcommand, string weightcommand, string histname, string axislabels, int xbins,float xlow,float xhigh)
{
  //assumes c has just enough divisions to cover all of these.
  int nsets=n_sample_bins;
  int nsubs=n_layouts;
  float setwidth=(sample_max-sample_min)/n_sample_bins;
  TH1F* ht;
  TF1* ft;
  for (int i = 0; i < nsets; i++){
    for (int j = 0; j < nsubs; j++){
      c->cd(j+1);
      ht=histout[i*nsubs+j]= new TH1F((string(histname) + string(to_string(i))+ string(layout[j])).c_str(),(string(layout[j])+string(axislabels)).c_str(),xbins,xlow,xhigh);
      ft=fitout[i*nsubs+j]=new TF1(("fit"+histname + to_string(i)+ layout[j]).c_str(),"gaus(0)",xlow,xhigh);

      ft->SetParNames("scale","x0","sigma");
      if (!isZombie[j])
	ntuple[j]->Draw((drawcommand+">>"+histname + to_string(i)+ layout[j]).c_str(),Form("(%s)*(%f<%s && %s<%f)",weightcommand.c_str(),setwidth*i+sample_min,sample_name.c_str(),sample_name.c_str(),setwidth*(i+1)+sample_min));

      histout[i*nsubs+j]->Draw();
      ft->SetParameters(ht->GetEntries(),ht->GetMean(),ht->GetRMS());
      ht->Fit(ft);
    }
    c->SaveAs((outpath+histname+"_"+outname+"_"+to_string(i)+".pdf").c_str());

  }
  return;
}

void drawFitParamAndSaveSet(TH1F **histout, TCanvas *c, TF1 **fit, int param, string histname, string axislabels, float minrange, float maxrange)
{
  TLegend *leg;
  c->cd();
  int nsets=n_sample_bins;
  int nsubs=n_layouts;
  float setwidth=(sample_max-sample_min)/n_sample_bins;
  leg = new TLegend(0.46,0.75,0.72,0.9);
  leg->SetNColumns(2);
    leg->AddEntry("",sPHENIXbanner.c_str(),"");
    leg->AddEntry("","","");//for two columns, we need an additional spacer.

    for (int i=0;i<nsubs;i++){
      histout[i]=new TH1F((string(histname)+ string(layout[i])).c_str(),axislabels.c_str(),n_sample_bins,sample_min,sample_max);
      histout[i]->SetLineColor(color[i]);
      histout[i]->SetMarkerColor(color[i]);
      //    histout[i]->SetLineWidth(2);
      for (int j=0;j<nsets;j++){
	int bin=histout[i]->Fill(setwidth*(j+0.5)+sample_min,fit[j*n_layouts+i]->GetParameter(param));
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
	  
  c->SaveAs((outpath+histname+"_"+outname+".pdf").c_str());
  return;
}

void drawFitErrorAndSaveSet(TH1F **histout, TCanvas *c, TF1 **fit, int param, string histname, string axislabels, float minrange, float maxrange)
{
  TLegend *leg;
  c->cd();
  int nsets=n_sample_bins;
  int nsubs=n_layouts;
  float setwidth=(sample_max-sample_min)/n_sample_bins;
  leg = new TLegend(0.20,0.75,0.480,0.9);
  leg->SetNColumns(2);
  leg->AddEntry("",sPHENIXbanner.c_str(),"");
    leg->AddEntry("","","");//for two columns, we need an additional spacer.

  for (int i=0;i<n_layouts;i++){
    histout[i]=new TH1F((string(histname)+ string(layout[i])).c_str(),axislabels.c_str(),n_sample_bins,sample_min,sample_max);
    histout[i]->SetLineColor(color[i]);
        histout[i]->SetMarkerColor(color[i]);

    histout[i]->SetLineWidth(2);
    for (int j=0;j<n_datasets;j++){
      int bin=histout[i]->Fill(setwidth*(j+0.5)+sample_min,fit[j*n_layouts+i]->GetParError(param));
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
	  
  c->SaveAs((outpath+histname+"_"+outname+".pdf").c_str());
  return;
}

void drawRMSAndSaveSet(TH1F **histout, TCanvas *c, TH1F **histin, string histname, string axislabels, float minrange, float maxrange)
{
  TLegend *leg;
  c->cd();
	
  leg = new TLegend(0.20,0.15,0.480,0.3);
  leg->SetNColumns(2);
  leg->AddEntry("",sPHENIXbanner.c_str(),"");
    leg->AddEntry("","","");//for two columns, we need an additional spacer.

  for (int i=0;i<n_layouts;i++){
    histout[i]=new TH1F((string(histname)+ string(layout[i])).c_str(),axislabels.c_str(),n_sample_bins,sample_min,sample_max);
    histout[i]->SetLineColor(color[i]);
        histout[i]->SetMarkerColor(color[i]);

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
	  
  c->SaveAs((outpath+histname+"_"+outname+".pdf").c_str());
  return;
}











void drawStaticAndSave(TH1F **histout, TCanvas *c, string histname, string axislabels, float minrange, float maxrange)
{
  const int ndata=8;
   const float pt[8]={0.5,0.6,0.7,0.8,0.9,1.0,2.0,3.0};
 const float a[8][8] = {
	{0.16428, 0.13812, 0.12134, 0.10483, 0.09357, 0.08661, 0.05171, 0.03864},
	{0.16853, 0.13518, 0.12134, 0.10648, 0.09618, 0.08767, 0.05439, 0.03985},
	{0.12238, 0.10615, 0.08938, 0.08360, 0.07551, 0.06882, 0.04373, 0.03267},
	{0.13065, 0.10902, 0.09148, 0.08199, 0.07423, 0.06986, 0.04110, 0.03386},
	{0.11011, 0.09197, 0.07899, 0.06924, 0.06413, 0.05857, 0.03588, 0.02914},
	{0.10204, 0.08638, 0.07285, 0.06924, 0.06288, 0.05554, 0.03848, 0.03031},
	{0.10606, 0.09197, 0.07899, 0.06924, 0.06288, 0.05756, 0.03588, 0.02914},
	{0.08627, 0.07269, 0.06282, 0.05684, 0.05303, 0.04661, 0.03074, 0.02564} };


  
  TLegend *leg;
  c->cd();
  int nsets=n_sample_bins;
  int nsubs=n_layouts;
  float setwidth=(sample_max-sample_min)/n_sample_bins;
  leg = new TLegend(0.46,0.75,0.72,0.9);
  leg->SetNColumns(2);
    leg->AddEntry("",sPHENIXbanner.c_str(),"");
    leg->AddEntry("","","");//for two columns, we need an additional spacer.

    for (int i=0;i<nsubs;i++){
      histout[i]=new TH1F((string(histname)+ string(layout[i])).c_str(),axislabels.c_str(),n_sample_bins,sample_min,sample_max);
      histout[i]->SetLineColor(color[i]);
      histout[i]->SetMarkerColor(color[i]);
      //    histout[i]->SetLineWidth(2);
      for (int j=0;j<ndata;j++){
	int bin=histout[i]->Fill(pt[j],a[i][j]);
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
	  
  c->SaveAs((outpath+histname+"_"+outname+".pdf").c_str());
  return;
}




a[8][8] = {
	{0.16428, 0.13812, 0.12134, 0.10483, 0.09357, 0.08661, 0.05171, 0.03864},
	{0.16853, 0.13518, 0.12134, 0.10648, 0.09618, 0.08767, 0.05439, 0.03985},
	{0.12238, 0.10615, 0.08938, 0.08360, 0.07551, 0.06882, 0.04373, 0.03267},
	{0.13065, 0.10902, 0.09148, 0.08199, 0.07423, 0.06986, 0.04110, 0.03386},
	{0.11011, 0.09197, 0.07899, 0.06924, 0.06413, 0.05857, 0.03588, 0.02914},
	{0.10204, 0.08638, 0.07285, 0.06924, 0.06288, 0.05554, 0.03848, 0.03031},
	{0.10606, 0.09197, 0.07899, 0.06924, 0.06288, 0.05756, 0.03588, 0.02914},
	{0.08627, 0.07269, 0.06282, 0.05684, 0.05303, 0.04661, 0.03074, 0.02564} };
