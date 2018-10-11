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

void kalman_uncertainty_vs_config() {
  //rcc says: don't need to do this just to read basic structures:  gSystem->Load("libg4hough.so");
  //rcc says: coding conventions usually say 'h' stands for a histogram and 'c' for a canvas, 'f' for a function.  maybe 'g' for a graph?  I don't use them often enough.  To me 'g' sounds like a constant, like 'k'.

	// define constants
	const int n_layouts = 8;				// number of INTT layouts we are testing
	const int n_maps_layers = 3;				// number of MVTX layers is always the same
	const int n_intt_layers[n_layouts] = {4,4,3,3,2,2,1,0};	// number of INTT layers in each layout
	const string path = "~/sphenix/data/by_mom/";	// absolute path to where data is located
	const string layout[n_layouts] = {"0111","1101","101","111","01","11","1","n"};	// convert layout# to layout configuration
	bool isZombie[n_layouts];

	// open files
	TFile *fin[n_layouts];

	for (int i = 0; i < n_layouts; ++i){
	  isZombie[i]=false;
	  fin[i] = new TFile((path + "G4_sPHENIX_" + layout[i] + ".root_g4kalman_eval.root").c_str(),"READ");
	  isZombie[i]=fin[i]->IsZombie();
	}
	

	// declare branch pointers
	TTree *branch[n_layouts];	// pointers to branches in each of the respective data files

	// leaf names to pass by reference
	float sigma_phi, sigma_z, sigma_r;
	
	// create histograms
	TH1F *hSigmaZ[n_layouts];
	TH1F *hSigmaPHI[n_layouts];
	TH1F *hSigmaR[n_layouts];
	TH2F *hSigmaPhiSigmaZ[n_layouts];
	TH1F *hDeltaPhi[n_layouts];
	TH1F *hDeltaZ[n_layouts];
	TH2F *hDeltaPhiDeltaZ[n_layouts];
	
	
	for (int i = 0; i < n_layouts; ++i) {
	  hSigmaZ[i] = new TH1F(("sigmaZ" + layout[i]).c_str(),(layout[i]+";#sigma_Z (um)").c_str(),100,-10,800);
		//hSigmaZ[i]->SetMaximum(60);
		hSigmaZ[i]->SetLineColor(i+1);

		hSigmaPHI[i] = new TH1F(("sigmaPHI" + layout[i]).c_str(),(layout[i]+";#sigma_#phi (um)").c_str(),100,-10,600);
		//hSigmaPHI[i]->SetMaximum(300);
		hSigmaPHI[i]->SetLineColor(i+1);

		hSigmaR[i] = new TH1F(("sigmaR" + layout[i]).c_str(),(layout[i]+";#sigma_R (um)").c_str(),100,-10,100);
		//hSigmaR[i]->SetMaximum(200); setting the limits in the definition above fixes them anyway.
		hSigmaR[i]->SetLineColor(i+1);

		hSigmaPhiSigmaZ[i] = new TH2F(("sigmaPhiSigmaZ" + layout[i]).c_str(),(layout[i]+";#sigma_Z (um);#sigma_#phi (um)").c_str(),50,-10,800,50,-10,600);

		hDeltaPhi[i]= new TH1F(("hDeltaPhi" + layout[i]).c_str(),(layout[i]+";#phi_guess-#phi_true (um)").c_str(),100,-500,500);
		hDeltaZ[i]= new TH1F(("hDeltaZ" + layout[i]).c_str(),(layout[i]+";z_guess-z_true (um)").c_str(),200,-2500,2500);
		hDeltaPhiDeltaZ[i]= new TH2F(("hDeltaPhiDeltaZ" + layout[i]).c_str(),(layout[i]+";#phi_guess-#phi_true (um);z_guess-z_true (um)").c_str(),100,-500,500,200,-2500,2500);
		//hSigmPhiSigmaZa[i]->SetLineColor(i+1);
	}

	// fill histograms
	for (int j = 0; j < n_layouts; ++j) {
	  if (isZombie[j]) continue; //skip them if the file is broken.
		fin[j]->GetObject("kalman_extrapolation_eval",branch[j]);
		int n_entries = branch[j]->GetEntries();
		branch[j]->SetBranchAddress("sigma_z2",&sigma_z);
		branch[j]->SetBranchAddress("sigma_phi2",&sigma_phi);
		branch[j]->SetBranchAddress("sigma_r2",&sigma_r);
		branch[j]->Draw(("(z2-z2t)*1e4:(phi2-phi2t)*1e4>>hDeltaPhiDeltaZ" + layout[j]).c_str());
		branch[j]->Draw(("(phi2-phi2t)*1e4>>hDeltaPhi" + layout[j]).c_str());
		branch[j]->Draw(("(z2-z2t)*1e4>>hDeltaZ" + layout[j]).c_str());

		for (int i = 0; i < n_entries; ++i) {
			branch[j]->GetEvent(i);

			// get rid of strange outliers
			// since you've set the bounds of your plots, you don't need to do this check:  if ((sigma_z > -100) & (sigma_phi > -100) & (sigma_r > -100)) {
			//I added the conversion to the standard deviation, in microns:
			  hSigmaZ[j]->Fill(TMath::Sqrt(sigma_z)*1e4);
			  hSigmaPHI[j]->Fill(TMath::Sqrt(sigma_phi)*1e4);
			  hSigmaR[j]->Fill(TMath::Sqrt(sigma_r)*1e4);
			  hSigmaPhiSigmaZ[j]->Fill(TMath::Sqrt(sigma_z)*1e4,TMath::Sqrt(sigma_phi)*1e4);//2D histogram
			  //}
		}
				}

	//create a canvas to draw all the 2D plots:
	TCanvas *c0=new TCanvas("c0","c0",1000,800);
	c0->Divide(4,2);
	c0->SetTitle(" track fit minus true at ~30cm");
	for (int i=0;i<n_layouts;i++){
	  c0->cd(i+1);
	  hSigmaPhiSigmaZ[i]->Draw("colz");
	}
	c0->SaveAs("SigmaPhi_vs_SigmaZ_multi.pdf");

	c0->SetTitle(" track fit minus true at ~30cm");
	for (int i=0;i<n_layouts;i++){
	  c0->cd(i+1);
	  hDeltaPhiDeltaZ[i]->Draw("colz");
	}
	c0->SaveAs("DeltaPhi_vs_DeltaZ_multi.pdf");
	for (int i=0;i<n_layouts;i++){
	  c0->cd(i+1);
	  hDeltaPhi[i]->Draw();
	}
	c0->SaveAs("DeltaPhi_multi.pdf");
	for (int i=0;i<n_layouts;i++){
	  c0->cd(i+1);
	  hDeltaZ[i]->Draw();
	}
	c0->SaveAs("DeltaZ_multi.pdf");
	// create canvas and legend
	TCanvas *c = new TCanvas("c1","c1",1000,800);
	TLegend *legend;

	// draw and format Z histogram
	hSigmaZ[0]->Draw();
	for (int i = 1; i < n_layouts; ++i) {	
		hSigmaZ[i]->Draw("same");
	}

	legend = new TLegend(0.11,0.55,0.25,0.89);

	for (int i = 0; i < n_layouts; ++i) {
		legend->AddEntry(hSigmaZ[i],layout[i].c_str(),"l");
	}

	legend->Draw();
	c->SaveAs("SigmaZ_vs_config.png");
	delete legend;

	// draw and format PHI histogram
	hSigmaPHI[0]->Draw();

	for (int i = 1; i < n_layouts; ++i) {		
		hSigmaPHI[i]->Draw("same");
	}

	legend = new TLegend(0.11,0.55,0.25,0.89);

	for (int i = 0; i < n_layouts; ++i) {
		legend->AddEntry(hSigmaPHI[i],layout[i].c_str(),"l");
	}

	legend->Draw();
	c->SaveAs("SigmaPHI_vs_config.png");
	delete legend;

	// draw and format R histogram
	hSigmaR[0]->Draw();

	for (int i = 1; i < n_layouts; ++i) {		
		hSigmaR[i]->Draw("same");
	}

	legend = new TLegend(0.11,0.55,0.25,0.89);

	for (int i = 0; i < n_layouts; ++i) {
		legend->AddEntry(hSigmaR[i],layout[i].c_str(),"l");
	}

	legend->Draw();
	c->SaveAs("SigmaR_vs_config.png");
	delete legend;

	// plot RMS and error for Z
	double x[n_layouts], y[n_layouts], ex[n_layouts], ey[n_layouts];

	for (int i = 0; i < n_layouts; ++i) {
		x[i] = hSigmaZ[i]->GetMean();
		ex[i] = hSigmaZ[i]->GetMeanError();
		y[i] = i;
		ey[i] = 0;
	}
	c->Clear();
	//you can also do this with a histogram.  Make one bin per layout, and use a fill that also takes an error value (or setbincontent, setbinerror):


	// a short routine to compile a set of means from an array of histograms:
	TH1F *hRMSZ=new TH1F("hRMSDelta","RMS values of #Deltasigma_z;layout;RMS (um)",n_layouts,-0.5,n_layouts-0.5);//I do the -0.5 to make sure the integers are centered in the bins.  Mostly a style thing.
	TH1F **histarray=hDeltaZ;
	TH1F *histout=hRMSZ;
	histout->SetCanExtend(TH1::kAllAxes); //lets us keep adding bins even if we go over the number we defined.
	//histout->SetLine
	//	histout->SetOptStat(0);
	for (int i=0;i<n_layouts;i++){
	  int bin=histout->Fill(layout[i].c_str(),histarray[i]->GetRMS());
	  histout->SetBinError(bin,histarray[i]->GetRMSError());
	}
	histout->Draw();
	c->SaveAs("ZRMS_vs_config.png");


	TH1F *hRMSPhi=new TH1F("hRMSPhi","RMS values of #Deltasigma_#phi;layout;RMS (um)",n_layouts,-0.5,n_layouts-0.5);//I do the -0.5 to make sure the integers are centered in the bins.  Mostly a style thing.
	histarray=hDeltaPhi;
        histout=hRMSPhi;
	histout->SetCanExtend(TH1::kAllAxes); //lets us keep adding bins even if we go over the number we defined.
	//histout->SetLine
	//histout->SetOptStat(0);
	for (int i=0;i<n_layouts;i++){
	  int bin=histout->Fill(layout[i].c_str(),histarray[i]->GetRMS());
	  histout->SetBinError(bin,histarray[i]->GetRMSError());
	}
	histout->Draw();
	c->SaveAs("PhiRMS_vs_config.png");

	
	TGraphErrors *gRMSZ = new TGraphErrors(n_layouts,x,y,ex,ey);
	
	hRMSZ->Draw();
	c->SaveAs("ZRMS_vs_config.png");
	for (int i = 0; i < n_layouts; ++i) {
		cout << layout[i] << "\t" << hSigmaZ[i]->GetRMS() << "\t" << hSigmaZ[i]->GetRMSError() << endl;
	}

	// close files so we can open new ones
	for (int i = 0; i < n_layouts; ++i) {
		fin[i]->Close();
		delete fin[i];
	}

	// prove INTT layout by plotting error in Z versus R
	float ez, r;
	TH2F *hZResVsR[n_layouts];
	TCanvas *c2=new TCanvas("c2","c2",1000,800);
	c2->Divide(4,2);
	for (int i = 0; i < n_layouts; ++i) {
		c2->cd(i+1);		
		fin[i] = new TFile((path + "G4_sPHENIX_" + layout[i] + ".root_g4svtx_eval.root").c_str(),"READ");
		if (fin[i]->IsZombie()) continue;

		fin[i]->GetObject("ntp_cluster",branch[i]);
		hZResVsR[i] = new TH2F(("hZResVsR" + layout[i]).c_str(),(layout[i]+";R (cm);Z hitres (um?)").c_str(),60,0,16,10,0,0.6);
		branch[i]->Draw(("ez:r>>hZResVsR"+layout[i]).c_str());//,"1","colz");
		hZResVsR[i]->Draw("colz");
		//	c->SaveAs(("Zerror_vs_R_" + layout[i] + ".png").c_str());
	}
	c2->SaveAs("ZRes_vs_R_multi.pdf");

	return;
	// close files
	for (int i = 0; i < n_layouts; ++i) {
		fin[i]->Close();
		delete fin[i];
	}
	// delete variables
	// when histograms are deleted something crashes, so for now there is a memory leak
	delete c;
	delete legend;

	gSystem->Exit(0);
}
