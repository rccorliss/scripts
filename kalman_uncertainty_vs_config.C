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
  // don't need to do this just to read basic structures:  gSystem->Load("libg4hough.so");

	// define constants
	const int n_layouts = 6;				// number of INTT layouts we are testing
	const int n_maps_layers = 3;				// number of MVTX layers is always the same
	const int n_intt_layers[n_layouts] = {4,4,3,3,2,2};	// number of INTT layers in each layout
	const string path = "~/sphenix/data/";	// absolute path to where data is located
	const string layout[n_layouts] = {"0111","1101","101","111","01","11"};	// convert layout# to layout configuration

	// open files
	TFile *fin[n_layouts];

	for (int i = 0; i < n_layouts; ++i)
		fin[i] = new TFile((path + "G4_sPHENIX_" + layout[i] + ".root_g4kalman_eval.root").c_str(),"READ");

	// declare branch pointers
	TTree *branch[n_layouts];	// pointers to branches in each of the respective data files

	// leaf names to pass by reference
	float sigma_phi, sigma_z, sigma_r;
	
	// create histograms
	TH1F *hSigmaZ[n_layouts];
	TH1F *hSigmaPHI[n_layouts];
	TH1F *hSigmaR[n_layouts];
	
	for (int i = 0; i < n_layouts; ++i) {
		hSigmaZ[i] = new TH1F(("sigmaZ" + layout[i]).c_str(),layout[i].c_str(),100,-10,500);
		hSigmaZ[i]->SetMaximum(60);
		hSigmaZ[i]->SetLineColor(i+1);

		hSigmaPHI[i] = new TH1F(("sigmaPHI" + layout[i]).c_str(),layout[i].c_str(),100,-10,500);
		hSigmaPHI[i]->SetMaximum(300);
		hSigmaPHI[i]->SetLineColor(i+1);

		hSigmaR[i] = new TH1F(("sigmaR" + layout[i]).c_str(),layout[i].c_str(),100,-10,500);
		hSigmaR[i]->SetMaximum(200);
		hSigmaR[i]->SetLineColor(i+1);
	}

	// fill histograms
	for (int j = 0; j < n_layouts; ++j) {
		fin[j]->GetObject("kalman_extrapolation_eval",branch[j]);
		int n_entries = branch[j]->GetEntries();
		branch[j]->SetBranchAddress("sigma_z",&sigma_z);
		branch[j]->SetBranchAddress("sigma_phi",&sigma_phi);
		branch[j]->SetBranchAddress("sigma_r",&sigma_r);

		for (int i = 0; i < n_entries; ++i) {
			branch[j]->GetEvent(i);

			// get rid of strange outliers
			// since you've set the bounds of your plots, you don't need to do this check:  if ((sigma_z > -100) & (sigma_phi > -100) & (sigma_r > -100)) {
			//I added the conversion to the standard deviation, in microns:
			  hSigmaZ[j]->Fill(TMath::Sqrt(sigma_z)*1e4);
			  hSigmaPHI[j]->Fill(TMath::Sqrt(sigma_phi)*1e4);
			  hSigmaR[j]->Fill(TMath::Sqrt(sigma_r)*1e4);
			  //}
		}
	}
	// create canvas and legend
	TCanvas *h = new TCanvas("c1","c1",700,500);
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
	h->SaveAs("SigmaZ_vs_config.png");
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
	h->SaveAs("SigmaPHI_vs_config.png");
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
	h->SaveAs("SigmaR_vs_config.png");
	delete legend;

	// plot RMS and error for Z
	double x[n_layouts], y[n_layouts], ex[n_layouts], ey[n_layouts];

	for (int i = 0; i < n_layouts; ++i) {
		x[i] = hSigmaZ[i]->GetRMS();
		ex[i] = hSigmaZ[i]->GetRMSError();
		y[i] = i;
		ey[i] = 0;
	}
	h->Clear();
	TGraphErrors *hRMSZ = new TGraphErrors(n_layouts,x,y,ex,ey);
	hRMSZ->Draw();
	h->SaveAs("ZRMS_vs_config.png");

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
	for (int i = 0; i < n_layouts; ++i) {
		fin[i] = new TFile((path + "G4_sPHENIX_" + layout[i] + ".root_g4svtx_eval.root").c_str(),"READ");
		fin[i]->GetObject("ntp_cluster",branch[i]);

		branch[i]->Draw("ez:r","(r<20)","colz");
		h->SaveAs(("Zerror_vs_R_" + layout[i] + ".png").c_str());
	}

	// close files
	for (int i = 0; i < n_layouts; ++i) {
		fin[i]->Close();
		delete fin[i];
	}

	// delete variables
	// when histograms are deleted something crashes, so for now there is a memory leak
	delete h;
	delete legend;

	gSystem->Exit(0);
}
