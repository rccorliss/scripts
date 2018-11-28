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
const string studypath="mom_scan_112418_manual";//doubles as a naming convention?

const string outpath="112418/";
const string particlename="pi+";

#define set1

#ifdef set1
const string outname="classic";
const int n_layouts = 7;
const int n_datasets=8;
const string pathsuffix[n_datasets]={"_pt0.5GeV_phi-180-180d_z0cm_eta-1.2-1.2",
				  "_pt0.6GeV_phi-180-180d_z0cm_eta-1.2-1.2",
				  "_pt0.7GeV_phi-180-180d_z0cm_eta-1.2-1.2",
				  "_pt0.8GeV_phi-180-180d_z0cm_eta-1.2-1.2",
				  "_pt0.9GeV_phi-180-180d_z0cm_eta-1.2-1.2",
				  "_pt1.0GeV_phi-180-180d_z0cm_eta-1.2-1.2",
				  "_pt2.0GeV_phi-180-180d_z0cm_eta-1.2-1.2",
				  "_pt3.0GeV_phi-180-180d_z0cm_eta-1.2-1.2"};
const string setname[n_datasets]={"0.5","0.6","0.7","0.8","0.9","1.0","2.0","3.0"};
const string layout[n_layouts]={"000p","00pp","00zp","0ppp","0pzp","ppzp","zppp"};
#endif

#ifdef set2
const int n_layouts = 2;				// number of INTT layouts we are testing
const string outname="00pp_vs_00zp";
const string layout[n_layouts]={"00pp",  "00zp"};
const int n_datasets=8;
const string pathsuffix[n_datasets]={"_pt0.5GeV_phi-180-180d_z0cm_eta-1.2-1.2",
				  "_pt0.6GeV_phi-180-180d_z0cm_eta-1.2-1.2",
				  "_pt0.7GeV_phi-180-180d_z0cm_eta-1.2-1.2",
				  "_pt0.8GeV_phi-180-180d_z0cm_eta-1.2-1.2",
				  "_pt0.9GeV_phi-180-180d_z0cm_eta-1.2-1.2",
				  "_pt1.0GeV_phi-180-180d_z0cm_eta-1.2-1.2",
				  "_pt2.0GeV_phi-180-180d_z0cm_eta-1.2-1.2",
				  "_pt3.0GeV_phi-180-180d_z0cm_eta-1.2-1.2"};
const string setname[n_datasets]={"0.5","0.6","0.7","0.8","0.9","1.0","2.0","3.0"};
#endif
#ifdef set3
const int n_layouts = 3;				// number of INTT layouts we are testing
const string outname="pp_ppp_zppp";
const string layout[n_layouts]={"00pp", "0ppp","zppp"};
const int n_datasets=8;
const string pathsuffix[n_datasets]={"_pt0.5GeV_phi-180-180d_z0cm_eta-1.2-1.2",
				  "_pt0.6GeV_phi-180-180d_z0cm_eta-1.2-1.2",
				  "_pt0.7GeV_phi-180-180d_z0cm_eta-1.2-1.2",
				  "_pt0.8GeV_phi-180-180d_z0cm_eta-1.2-1.2",
				  "_pt0.9GeV_phi-180-180d_z0cm_eta-1.2-1.2",
				  "_pt1.0GeV_phi-180-180d_z0cm_eta-1.2-1.2",
				  "_pt2.0GeV_phi-180-180d_z0cm_eta-1.2-1.2",
				  "_pt3.0GeV_phi-180-180d_z0cm_eta-1.2-1.2"};
const string setname[n_datasets]={"0.5","0.6","0.7","0.8","0.9","1.0","2.0","3.0"};
#endif
#ifdef set4
const int n_layouts = 2;				// number of INTT layouts we are testing
const string outname="ppzp_vs_zppp";
const string layout[n_layouts]={"ppzp","zppp"};
const int n_datasets=8;
const string pathsuffix[n_datasets]={"_pt0.5GeV_phi-180-180d_z0cm_eta-1.2-1.2",
				  "_pt0.6GeV_phi-180-180d_z0cm_eta-1.2-1.2",
				  "_pt0.7GeV_phi-180-180d_z0cm_eta-1.2-1.2",
				  "_pt0.8GeV_phi-180-180d_z0cm_eta-1.2-1.2",
				  "_pt0.9GeV_phi-180-180d_z0cm_eta-1.2-1.2",
				  "_pt1.0GeV_phi-180-180d_z0cm_eta-1.2-1.2",
				  "_pt2.0GeV_phi-180-180d_z0cm_eta-1.2-1.2",
				  "_pt3.0GeV_phi-180-180d_z0cm_eta-1.2-1.2"};
const string setname[n_datasets]={"0.5","0.6","0.7","0.8","0.9","1.0","2.0","3.0"};
#endif
//const string layout[n_layouts]={"00pp",  "00zp","0ppp","0pzp","ppzp","zppp"};

#ifdef compare_particles
//this hijacks the existing structures, except that it changes how the initial files are loaded so they pull
// multiple particles rather than configurations.  Naming conventions should work fine.
const string single_layout="zppp";
const int n_layouts=4;
const string outname="particles";

const string layout[n_layouts]={"e-","mu-","pi-","pi+"}; 
const int n_datasets=8;
const string pathsuffix[n_datasets]={"_pt0.5GeV_phi-180-180d_z0cm_eta-1.2-1.2",
				  "_pt0.6GeV_phi-180-180d_z0cm_eta-1.2-1.2",
				  "_pt0.7GeV_phi-180-180d_z0cm_eta-1.2-1.2",
				  "_pt0.8GeV_phi-180-180d_z0cm_eta-1.2-1.2",
				  "_pt0.9GeV_phi-180-180d_z0cm_eta-1.2-1.2",
				  "_pt1.0GeV_phi-180-180d_z0cm_eta-1.2-1.2",
				  "_pt2.0GeV_phi-180-180d_z0cm_eta-1.2-1.2",
				  "_pt3.0GeV_phi-180-180d_z0cm_eta-1.2-1.2"};
const string setname[n_datasets]={"0.5","0.6","0.7","0.8","0.9","1.0","2.0","3.0"};
#endif


//and define variables that it'll make our life easier to be able to access directly:
TFile *fin[n_datasets][n_layouts];
bool isZombie[n_datasets][n_layouts];
TTree *ntuple[n_datasets][n_layouts];	// pointers to trees in each of the respective data files

//generally, I'm organizing these as returned pointers, the canvas, then the stuff needed to generate the histogram, and the names/labels/limits associated with it.
void drawAndSaveSet2D(TH2F **histout, TCanvas *c, string drawcommand, string histname, string axislabels, int xbins,float xlow,float xhigh,int ybins,float ylow,float yhigh);
void drawAndSaveSet1D(TH1F **histout, TCanvas *c, string drawcommand, string histname, string axislabels, int xbins,float xlow,float xhigh);
void drawAndFitAndSaveSet1D(TH1F **histout,TF1 **fitout, TCanvas *c, string drawcommand, string histname, string axislabels, int xbins,float xlow,float xhigh);
void drawAndFitAndSaveSet1D(TH1F **histout,TF1 **fitout, TCanvas *c, string drawcommand, string weightcommand, string histname, string axislabels, int xbins,float xlow,float xhigh); //with weighting so we can make a cut.
void drawFitParamAndSaveSet(TH1F **histout, TCanvas *c, TF1 **fit, int param, string histname, string axislabels, float minrange, float maxrange);
void drawFitErrorAndSaveSet(TH1F **histout, TCanvas *c, TF1 **fit, int param, string histname, string axislabels, float minrange, float maxrange);
void drawRMSAndSaveSet(TH1F **histout, TCanvas *c, TH1F **histin, string histname, string axislabels, float minrange, float maxrange);


void kalman_uncertainty_nov14() {
  //rcc says: don't need to do this just to read basic structures:  gSystem->Load("libg4hough.so");

  
  SetsPhenixStyle();

#ifndef compare_particles
  	const string path=basepath+studypath+"/";
	string subpath[n_datasets];
	for (int i = 0; i < n_datasets; i++){
	  subpath[i]=particlename+pathsuffix[i];
	}

	//open all the files we need.  first index is subpath, second index is layout.
	for (int i = 0; i < n_datasets; i++){
	  for (int j = 0; j < n_layouts; j++){
	    isZombie[i][j]=false;
	    fin[i][j] = new TFile((path + subpath[i]+"/G4_sPHENIX_" + layout[j] + ".root_g4kalman_eval.root").c_str(),"READ");
	    isZombie[i][j]=fin[i][j]->IsZombie();
	    if (isZombie[i][j]) continue; //skip them if the file is broken.
	    fin[i][j]->GetObject("kalman_eval",ntuple[i][j]);
	  }
	}
#endif
#ifdef compare_particles
	//open all the files we need.  first index is subpath, second index is layout -- in this case, particle species.
	const string path=basepath+studypath+"/";

	for (int i = 0; i < n_datasets; i++){
	  for (int j = 0; j < n_layouts; j++){
	    isZombie[i][j]=false;
	    fin[i][j] = new TFile((path + layout[j]+pathsuffix[i]+"/G4_sPHENIX_" + single_layout + ".root_g4kalman_eval.root").c_str(),"READ");
	    isZombie[i][j]=fin[i][j]->IsZombie();
	    if (isZombie[i][j]) continue; //skip them if the file is broken.
	    fin[i][j]->GetObject("kalman_eval",ntuple[i][j]);
	  }
	}
#endif

	// leaf names to pass by reference
	float sigma_phi, sigma_z, sigma_r;
	



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
	gStyle->SetOptFit(111111);



	

	
	//looking at regions we're hitting:
	TH2F *hTruePhiZ[n_datasets][n_layouts];
	drawAndSaveSet2D((TH2F**) hTruePhiZ,c0,"z30t:phi30t","hTruePhiZ",";phi g4 [rad.];z g4[cm]",50,-4,4,50,-100,100);

	//compare true and kalman extrapolated position
	TH2F *hTrueGuessPhiZ[n_datasets][n_layouts];
	drawAndSaveSet2D((TH2F**) hTrueGuessPhiZ,c0,"(z30te-z30t)*10:(phi30te-phi30t)*r30t*10","hDelta2D",";phi guess-g4* [mm];z guess-g4[mm]",50,-3,3,50,-15,15);

		//compare true and kalman extrapolated position

       	TH2F *hG4vsClusterTrack2D[n_datasets][n_layouts];
	drawAndSaveSet2D((TH2F**) hG4vsClusterTrack2D,c0,"(z30te-g4_z30te)*10:(phi30te-g4_phi30te)*r30te*10","hG4vsClusterTrack2D",";phi cltrack-g4track [mm];z cltrack-g4track[mm]",50,-3,3,50,-15,15);

	TH2F *hG4vsTruth2D[n_datasets][n_layouts];
	drawAndSaveSet2D((TH2F**) hG4vsClusterTrack2D,c0,"(g4_z30t-g4_z30te)*10:(g4_phi30t-g4_phi30te)*g4_r30t*10","hG4vsTruth2D",";phi truth-g4track [mm];z truth-g4track[mm]",50,-3,3,50,-15,15);

	//check ng4hits
	TH1F *hG4Hits[n_datasets][n_layouts];
	drawAndSaveSet1D((TH1F**) hG4Hits, c0,"ng430t","hNG4Hits",";ng4hits",3,-0.5,2.5);

	//compare phi of g4 and predicted position
	TH1F *hTrueGuessPhi[n_datasets][n_layouts];
	TF1 *fTrueGuessPhi[n_datasets][n_layouts];
	  drawAndFitAndSaveSet1D((TH1F**) hTrueGuessPhi,(TF1**) fTrueGuessPhi, c0,"(phi30te-phi30t)*r30t*10","hDeltaPhi",";phi kalman-g4* [mm]",80,-8,8);

	  //compare z of g4 and predicted position
	TH1F *hTrueGuessZ[n_datasets][n_layouts];
	TF1 *fTrueGuessZ[n_datasets][n_layouts];
	  drawAndFitAndSaveSet1D((TH1F**) hTrueGuessZ,(TF1**) fTrueGuessZ, c0,"(z30te-z30t)*10","hDeltaZ",";z guess-g4 [mm]",50,-3,3);

	//compare phi of g4 and predicted position
	TH1F *hTrueGuessPhig4[n_datasets][n_layouts];
	TF1 *fTrueGuessPhig4[n_datasets][n_layouts];
	  drawAndFitAndSaveSet1D((TH1F**) hTrueGuessPhig4,(TF1**) fTrueGuessPhig4, c0,"(g4_phi30te-g4_phi30t)*g4_r30t*10","hDeltaPhiG4",";phi g4kalman-g4truth [mm]",80,-8,8);

	  //compare z of g4 and predicted position
	TH1F *hTrueGuessZg4[n_datasets][n_layouts];
	TF1 *fTrueGuessZg4[n_datasets][n_layouts];
	  drawAndFitAndSaveSet1D((TH1F**) hTrueGuessZg4,(TF1**) fTrueGuessZg4, c0,"(g4_z30te-g4_z30t)*10","hDeltaZG4",";z g4kalman-g4truth [mm]",50,-3,3);

	  //was 10k for mom scan.

	  /*
	TH1F *hTrueGuessZtemp[n_datasets][n_layouts];
	TF1 *fTrueGuessZtemp[n_datasets][n_layouts];
	for (float ptlow=0;ptlow<ptmax;ptlow+=ptstep){
	drawAndFitAndSaveSet1D((TH1F**) hTrueGuessZg4,(TF1**) fTrueGuessZg4, c0,"(g4_z30te-g4_z30t)*10",Form("true_pti>%f && true_pti<%f",ptlow,pthigh),"hDeltaZtemp",";z g4kalman-g4truth [mm]",50,-3,3);
	}
	  */

	  

	//look at trends in true and predicted position
	//assumes z and phi histograms are packed the same way.

	//create a canvas to draw all the fit results plots:
	TCanvas *c1=new TCanvas("c1","c1",800,800);

	TH1F *hResolutionSigmaZ[n_layouts];
	drawFitParamAndSaveSet((TH1F**)hResolutionSigmaZ, c1, (TF1 **) fTrueGuessZ, 2, "hResZ", ";pT [GeV];z res. at 30cm [mm]",0,3);


	TH1F *hResolutionSigmaPhi[n_layouts];
	drawFitParamAndSaveSet((TH1F**)hResolutionSigmaPhi, c1, (TF1 **) fTrueGuessPhi, 2, "hResPhi",";pT [GeV];#phi res. at 30cm [mm]",0,3);

	TH1F *hResolutionSigmaZg4[n_layouts];
	drawFitParamAndSaveSet((TH1F**)hResolutionSigmaZg4, c1, (TF1 **) fTrueGuessZg4, 2, "hResZg4", ";pT [GeV];z res. at 30cm w/ g4track[mm]",0,3);


	TH1F *hResolutionSigmaPhig4[n_layouts];
	drawFitParamAndSaveSet((TH1F**)hResolutionSigmaPhig4, c1, (TF1 **) fTrueGuessPhig4, 2, "hResPhig4",";pT [GeV];#phi res. at 30cm w/ g4track[mm]",0,3);




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
  //printf(layout[1].c_str());
  //printf(axislabels.c_str());
  //printf("\n");
  //assumes c has enough divisions to cover all of these.
  int nsets=n_datasets;
  int nsubs=n_layouts;
  for (int i = 0; i < nsets; i++){
    for (int j = 0; j < nsubs; j++){
      histout[i*nsubs+j]= new TH2F((string(histname) + string(setname[i])+ string(layout[j])).c_str(),(string(layout[j])+string(axislabels)).c_str(),xbins,xlow,xhigh,ybins,ylow,yhigh);

      c->cd(j+1);
      if (!isZombie[i][j])
	ntuple[i][j]->Draw((drawcommand+">>"+histname + setname[i]+ layout[j]).c_str());
      histout[i*nsubs+j]->Draw("colz");
    }
    c->SaveAs((outpath+histname+"_"+outname+"_"+setname[i]+".pdf").c_str());

  }
  return;
}

void drawAndSaveSet1D(TH1F **histout, TCanvas *c, string drawcommand, string histname, string axislabels, int xbins,float xlow,float xhigh)
{
  //assumes c has just enough divisions to cover all of these.
  int nsets=n_datasets;
  int nsubs=n_layouts;
  TH1F* ht;
  for (int i = 0; i < nsets; i++){
    for (int j = 0; j < nsubs; j++){
      c->cd(j+1);
      ht=histout[i*nsubs+j]= new TH1F((string(histname) + string(setname[i])+ string(layout[j])).c_str(),(string(layout[j])+string(axislabels)).c_str(),xbins,xlow,xhigh);
      if (!isZombie[i][j])
	ntuple[i][j]->Draw((drawcommand+">>"+histname + setname[i]+ layout[j]).c_str());

      histout[i*nsubs+j]->Draw();
    }
    c->SaveAs((outpath+histname+"_"+outname+"_"+setname[i]+".pdf").c_str());

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
      ht=histout[i*nsubs+j]= new TH1F((string(histname) + string(setname[i])+ string(layout[j])).c_str(),(string(layout[j])+string(axislabels)).c_str(),xbins,xlow,xhigh);
      ft=fitout[i*nsubs+j]=new TF1(("fit"+histname + setname[i]+ layout[j]).c_str(),"gaus(0)",xlow,xhigh);

      ft->SetParNames("scale","x0","sigma");
      ft->SetParameters(ht->GetEntries(),ht->GetMean(),ht->GetRMS());
      if (!isZombie[i][j])
	ntuple[i][j]->Draw((drawcommand+">>"+histname + setname[i]+ layout[j]).c_str());

      histout[i*nsubs+j]->Draw();
      ht->Fit(ft);
    }
    c->SaveAs((outpath+histname+"_"+outname+"_"+setname[i]+".pdf").c_str());

  }
  return;
}

void drawAndFitAndSaveSet1D(TH1F **histout,TF1 **fitout, TCanvas *c, string drawcommand, string weightcommand, string histname, string axislabels, int xbins,float xlow,float xhigh)
{
  //assumes c has just enough divisions to cover all of these.
  int nsets=n_datasets;
  int nsubs=n_layouts;
  TH1F* ht;
  TF1* ft;
  for (int i = 0; i < nsets; i++){
    for (int j = 0; j < nsubs; j++){
      c->cd(j+1);
      ht=histout[i*nsubs+j]= new TH1F((string(histname) + string(setname[i])+ string(layout[j])).c_str(),(string(layout[j])+string(axislabels)).c_str(),xbins,xlow,xhigh);
      ft=fitout[i*nsubs+j]=new TF1(("fit"+histname + setname[i]+ layout[j]).c_str(),"gaus(0)",xlow,xhigh);

      ft->SetParNames("scale","x0","sigma");
      if (!isZombie[i][j])
	ntuple[i][j]->Draw((drawcommand+">>"+histname + setname[i]+ layout[j]).c_str(),weightcommand.c_str());

      histout[i*nsubs+j]->Draw();
      ft->SetParameters(ht->GetEntries(),ht->GetMean(),ht->GetRMS());
      ht->Fit(ft);
    }
    c->SaveAs((outpath+histname+"_"+outname+"_"+setname[i]+".pdf").c_str());

  }
  return;
}

void drawFitParamAndSaveSet(TH1F **histout, TCanvas *c, TF1 **fit, int param, string histname, string axislabels, float minrange, float maxrange)
{
  TLegend *leg;
  c->cd();
	
  leg = new TLegend(0.46,0.75,0.72,0.9);
  leg->SetNColumns(2);
    leg->AddEntry("",sPHENIXbanner.c_str(),"");
    leg->AddEntry("","","");//for two columns, we need an additional spacer.

    for (int i=0;i<n_layouts;i++){
      histout[i]=new TH1F((string(histname)+ string(layout[i])).c_str(),axislabels.c_str(),n_datasets,-0.5,n_datasets-0.5);
    histout[i]->SetLineColor(i+1);
    histout[i]->SetMarkerColor(i+1);
    //    histout[i]->SetLineWidth(2);
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
	  
  c->SaveAs((outpath+histname+"_"+outname+".pdf").c_str());
  return;
}

void drawFitErrorAndSaveSet(TH1F **histout, TCanvas *c, TF1 **fit, int param, string histname, string axislabels, float minrange, float maxrange)
{
  TLegend *leg;
  c->cd();
	
  leg = new TLegend(0.20,0.75,0.480,0.9);
  leg->SetNColumns(2);
  leg->AddEntry("",sPHENIXbanner.c_str(),"");
    leg->AddEntry("","","");//for two columns, we need an additional spacer.

  for (int i=0;i<n_layouts;i++){
    histout[i]=new TH1F((string(histname)+ string(layout[i])).c_str(),axislabels.c_str(),n_datasets,-0.5,n_datasets-0.5);
    histout[i]->SetLineColor(i+1);
        histout[i]->SetMarkerColor(i+1);

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
    histout[i]=new TH1F((string(histname)+ string(layout[i])).c_str(),axislabels.c_str(),n_datasets,-0.5,n_datasets-0.5);
    histout[i]->SetLineColor(i+1);
        histout[i]->SetMarkerColor(i+1);

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

















