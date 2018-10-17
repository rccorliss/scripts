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

void kalman_uncertainty_oct17() {
  //rcc says: don't need to do this just to read basic structures:  gSystem->Load("libg4hough.so");

	// define constants
	const int n_layouts = 8;				// number of INTT layouts we are testing
	const int n_maps_layers = 3;				// number of MVTX layers is always the same
	const int n_intt_layers[n_layouts] = {4,4,3,3,2,2,1,0};	// number of INTT layers in each layout
	const string basepath = "~/sphenix/data/";	// base path to where data is located


	/*
	// detector image sets -- large vertex range
	const int n_datasets=3;
	const string studypath="detector_image";//doubles as a naming convention?
	const string subpath[n_datasets]={"pi+_pt0.5GeV_phi-180-180d_z-10-10cm_eta-1.2-1.2",
					  "pi+_pt10GeV_phi-180-180d_z-10-10cm_eta-1.2-1.2",
					  "pi+_pt1GeV_phi-180-180d_z-10-10cm_eta-1.2-1.2"};
	const string setname[n_datasets]={"0.5GeV","1.0GeV","10GeV"};
	*/
	
	const int n_datasets=6;
	const string studypath="mom_scan";//doubles as a naming convention?
	const string subpath[n_datasets]={"pi+_pt0.5GeV_phi-180-180d_z0cm_eta-1.2-1.2",
					  "pi+_pt0.6GeV_phi-180-180d_z0cm_eta-1.2-1.2",
					  "pi+_pt0.7GeV_phi-180-180d_z0cm_eta-1.2-1.2",
					  "pi+_pt0.8GeV_phi-180-180d_z0cm_eta-1.2-1.2",
					  "pi+_pt0.9GeV_phi-180-180d_z0cm_eta-1.2-1.2",
					  "pi+_pt1GeV_phi-180-180d_z0cm_eta-1.2-1.2"};
	const string setname[n_datasets]={"0.5GeV","0.6GeV","0.7GeV","0.8GeV","0.9GeV","1.0GeV"};


	/*
	// high energy sets:
	// purpose:  prove that resolutions are sane at high energy
	// so:  plot z and phi deviations as a function of angle.  At fixed _pt_ they should be similar in phi and scaling in z.
	const int n_datasets=3;
	const string studypath="high_energy";//doubles as a naming convention?
	const string subpath[n_datasets]={"pi+_pt100GeV_phi0d_z0cm_eta0.2",
					  "pi+_pt100GeV_phi15d_z0cm_eta0.2",
					  "pi+_pt100GeV_phi30d_z0cm_eta0.2"};
	const string setname[n_datasets]={"0deg","15deg","30deg"};
	*/

	//common to all the sets:
	const string path=basepath+studypath+"/";//to limit how much we have to correct downstream.
	const string layout[n_layouts] = {"0111","1101","101","111","01","11","1","n"};	// convert layout# to layout configuration

	//  files, branches and their associated 'isZombie' (this could be called each time by fin[z][q]->isZombie() as well)
	TFile *fin[n_datasets][n_layouts];
	bool isZombie[n_datasets][n_layouts];
	TTree *tuple[n_datasets][n_layouts];	// pointers to trees in each of the respective data files

	//open all the files we need.  first index is subpath, second index is layour.
	for (int i = 0; i < n_datasets; i++){
	  for (int j = 0; j < n_layouts; j++){
	    isZombie[i][j]=false;
	    fin[i][j] = new TFile((path + subpath[i]+"/G4_sPHENIX_" + layout[j] + ".root_g4kalman_eval.root").c_str(),"READ");
	    isZombie[i][j]=fin[i][j]->IsZombie();
	    if (isZombie[i][j]) continue; //skip them if the file is broken.
	    fin[i][j]->GetObject("kalman_extrapolation_eval",tuple[i][j]);
	  }
	}


	// leaf names to pass by reference
	float sigma_phi, sigma_z, sigma_r;
	
	// create histograms
	TH1F *hSigmaZ[n_datasets][n_layouts];
	TH1F *hSigmaPhi[n_datasets][n_layouts];
	TH1F *hSigmaR[n_datasets][n_layouts];
	TH2F *hSigmaPhiSigmaZ[n_datasets][n_layouts];
	TH1F *hDeltaPhi[n_datasets][n_layouts];
	TH1F *hDeltaZ[n_datasets][n_layouts];
	TH2F *hDeltaPhiDeltaZ[n_datasets][n_layouts];
	
	
	for (int i = 0; i < n_datasets; i++){
	  for (int j = 0; j < n_layouts; j++){
	    hSigmaZ[i][j] = new TH1F(("sigmaZ" + subpath[i]+layout[j]).c_str(),(layout[j]+";#sigma_Z (um)").c_str(),100,-10,800);
	    hSigmaZ[i][j]->SetLineColor(i+1);

	    hSigmaPhi[i][j] = new TH1F(("sigmaPhi" + subpath[i]+ layout[j]).c_str(),(layout[j]+";#sigma_#phi (um)").c_str(),100,-10,600);
	    hSigmaPhi[i][j]->SetLineColor(i+1);
	    
	    hSigmaR[i][j] = new TH1F(("sigmaR" + subpath[i]+ layout[j]).c_str(),(layout[i]+";#sigma_R (um)").c_str(),100,-10,100);
	    //hSigmaR[i]->SetMaximum(200); setting the limits in the definition above fixes them anyway.
	    hSigmaR[i][j]->SetLineColor(i+1);
	    
	    hSigmaPhiSigmaZ[i][j] = new TH2F(("sigmaPhiSigmaZ" + subpath[i]+ layout[j]).c_str(),(layout[i]+";#sigma_Z (um);#sigma_#phi (um)").c_str(),50,-10,800,50,-10,600);
	    
	    hDeltaPhi[i][j]= new TH1F(("hDeltaPhi" + subpath[i]+ layout[j]).c_str(),(layout[i]+";#phi_guess-#phi_true (um)").c_str(),100,-500,500);
	    hDeltaZ[i][j]= new TH1F(("hDeltaZ" + subpath[i]+ layout[j]).c_str(),(layout[i]+";z_guess-z_true (um)").c_str(),200,-2500,2500);
	    hDeltaPhiDeltaZ[i][j]= new TH2F(("hDeltaPhiDeltaZ" + subpath[i]+ layout[j]).c_str(),(layout[i]+";#phi_guess-#phi_true (um);z_guess-z_true (um)").c_str(),100,-500,500,200,-2500,2500);
	    //hSigmPhiSigmaZa[i]->SetLineColor(i+1);
	  }
	}



	//start drawing things:
	
	//create a canvas to draw all the 2D plots:
	TCanvas *c0=new TCanvas("c0","c0",1000,800);

	
	//looking at regions we're hitting:
	TH2F *hClustPhiZ[n_datasets][n_layouts];
	c0->Divide(4,2);
	for (int i = 0; i < n_datasets; i++){
	  for (int j = 0; j < n_layouts; j++){
	    if (isZombie[i][j]) continue;
	    	    c0->cd(j+1);

	    hClustPhiZ[i][j]= new TH2F(("hClustPhiZ" + setname[i]+ layout[j]).c_str(),(layout[j]+";phi clust (rad.);z clust(cm)").c_str(),50,-4,4,50,-100,100);
	    tuple[i][j]->Draw(("z2t:phi2t>>hClustPhiZ" + setname[i]+ layout[j]).c_str());
	    hClustPhiZ[i][j]->Draw("colz");
	  }
	c0->SaveAs(("ClustPhiZ_multi_"+studypath+setname[i]+".pdf").c_str());
	}

	//repeat for true clusters:
	TH2F *hTruePhiZ[n_datasets][n_layouts];
	//c0->Divide(4,2);
	for (int i = 0; i < n_datasets; i++){
	  for (int j = 0; j < n_layouts; j++){
	    if (isZombie[i][j]) continue;
	    	    c0->cd(j+1);

	    hTruePhiZ[i][j]= new TH2F(("hTruePhiZ" + setname[i]+ layout[j]).c_str(),(layout[j]+";phi clust (rad.);z clust(cm)").c_str(),50,-4,4,50,-100,100);
	    tuple[i][j]->Draw(("z2t:phi2t>>hTruePhiZ" + setname[i]+ layout[j]).c_str());
	    hTruePhiZ[i][j]->Draw("colz");
	  }
	c0->SaveAs(("TruePhiZ_multi_"+studypath+setname[i]+".pdf").c_str());
	}

	//compare true and predicted position
	TH2F *hTrueGuessPhiZ[n_datasets][n_layouts];
	//c0->Divide(4,2);
	for (int i = 0; i < n_datasets; i++){
	  for (int j = 0; j < n_layouts; j++){
	    if (isZombie[i][j]) continue;
	    	    c0->cd(j+1);

	    hTrueGuessPhiZ[i][j]= new TH2F(("hTrueGuessPhiZ" + setname[i]+ layout[j]).c_str(),(layout[j]+";phi guess-clust* (um);z guess-clust(um)").c_str(),50,-4000,4000,50,-20000,20000);
	    tuple[i][j]->Draw(("(z2-z2t)*1e4:(phi2-phi2t)*r2*1e4>>hTrueGuessPhiZ" + setname[i]+ layout[j]).c_str());
	    hTrueGuessPhiZ[i][j]->Draw("colz");
	  }
	c0->SaveAs(("TrueGuessPhiZ_multi_"+studypath+setname[i]+".pdf").c_str());
	}

	//compare true and cluster position
	TH2F *hTrueClustPhiZ[n_datasets][n_layouts];
	//c0->Divide(4,2);
	for (int i = 0; i < n_datasets; i++){
	  for (int j = 0; j < n_layouts; j++){
	    if (isZombie[i][j]) continue;
	    c0->cd(j+1);

	    hTrueClustPhiZ[i][j]= new TH2F(("hTrueClustPhiZ" + setname[i]+ layout[j]).c_str(),(layout[j]+";phi clust-true* (um);z clust-true(um)").c_str(),50,-4000,4000,50,-20000,20000);
	    tuple[i][j]->Draw(("(z2c-z2t)*1e4:(phi2c-phi2t)*r2*1e4>>hTrueClustPhiZ" + setname[i]+ layout[j]).c_str());
	    hTrueClustPhiZ[i][j]->Draw("colz");
	  }
	c0->SaveAs(("TrueClustPhiZ_multi_"+studypath+setname[i]+".pdf").c_str());
	}

	
	//compare true and predicted position-- phi only
	TH1F *hTrueGuessPhi[n_datasets][n_layouts];
	//c0->Divide(4,2);
	for (int i = 0; i < n_datasets; i++){
	  for (int j = 0; j < n_layouts; j++){
	    if (isZombie[i][j]) continue;
	    c0->cd(j+1);
	    int bonus=1;
	    if (layout[j]=="n") bonus=10;
	    hTrueGuessPhi[i][j]= new TH1F(("hTrueGuessPhi" + setname[i]+ layout[j]).c_str(),(layout[j]+";phi guess-clust* (um)").c_str(),50,-4000*bonus,4000*bonus);
	    tuple[i][j]->Draw(("(phi2-phi2t)*r2*1e4>>hTrueGuessPhi" + setname[i]+ layout[j]).c_str());
	    hTrueGuessPhi[i][j]->Draw();
	  }
	c0->SaveAs(("TrueGuessPhi_multi_"+studypath+setname[i]+".pdf").c_str());
	}

		//compare true and predicted position-- z only
	TH1F *hTrueGuessZ[n_datasets][n_layouts];
	//c0->Divide(4,2);
	for (int i = 0; i < n_datasets; i++){
	  for (int j = 0; j < n_layouts; j++){
	    if (isZombie[i][j]) continue;
	    	    c0->cd(j+1);

	    hTrueGuessZ[i][j]= new TH1F(("hTrueGuessZ" + setname[i]+ layout[j]).c_str(),(layout[j]+";z guess-clust(um)").c_str(),50,-20000,20000);
	    tuple[i][j]->Draw(("(z2-z2t)*1e4>>hTrueGuessZ" + setname[i]+ layout[j]).c_str());
	    hTrueGuessZ[i][j]->Draw();
	  }
	c0->SaveAs(("TrueGuessZ_multi_"+studypath+setname[i]+".pdf").c_str());
	}


	//look at trends in true and predicted position
	//assumes z and phi histograms are packed the same way.

	//create a canvas to draw all the 2D plots:
	TCanvas *c1=new TCanvas("c1","c1",800,800);
	
	//plot one trend for each of the geometries, one histogram cell for each of the datasets
	TH1F *histout[n_layouts];
	TLegend *leg[2];
	c1->cd();
	  leg[0] = new TLegend(0.20,0.15,0.40,0.3);
	  leg[1] = new TLegend(0.60,0.15,0.80,0.3);

	for (int i=0;i<n_layouts;i++){
	  histout[i]=new TH1F(("hRMSDeltaZ"+ layout[i]).c_str(),"RMS values of #Delta_z;;RMS (um)",n_datasets,-0.5,n_datasets-0.5);
	  histout[i]->SetLineColor(i+1);
	  histout[i]->SetLineWidth(2);
	  for (int j=0;j<n_datasets;j++){
	    int bin=histout[i]->Fill(setname[j].c_str(),hTrueGuessZ[j][i]->GetRMS());
	    histout[i]->SetBinError(bin,hTrueGuessZ[j][i]->GetRMSError());
	  }

	  if (i==0){
	    histout[i]->GetYaxis()->SetRangeUser(0,5000);
	    histout[i]->Draw();
	  }else {
	    histout[i]->Draw("SAME");
	  }
	}
	for (int i=0;i<n_layouts;i++){
	  if (i>=n_layouts/2){
	    leg[1]->AddEntry(histout[i],layout[i].c_str(),"l");
	  }else{ leg[0]->AddEntry(histout[i],layout[i].c_str(),"l");}
	}
	leg[0]->Draw();
	leg[1]->Draw();
	  
	c1->SaveAs(("ZRMS_vs_config_"+studypath+".pdf").c_str());


	leg[0] = new TLegend(0.20,0.65,0.40,0.8);
	leg[1] = new TLegend(0.60,0.65,0.80,0.8);

	for (int i=0;i<n_layouts;i++){
	  histout[i]=new TH1F(("hRMSDeltaPhi"+ layout[i]).c_str(),"RMS values of #Delta_ #phi;;RMS (um)",n_datasets,-0.5,n_datasets-0.5);
	  histout[i]->SetLineColor(i+1);
	  histout[i]->SetLineWidth(2);
	  for (int j=0;j<n_datasets;j++){
	    int bin=histout[i]->Fill(setname[j].c_str(),hTrueGuessPhi[j][i]->GetRMS());
	    histout[i]->SetBinError(bin,hTrueGuessPhi[j][i]->GetRMSError());
	  }

	  if (i==0){
	    histout[i]->GetYaxis()->SetRangeUser(0,5000);
	    histout[i]->Draw();
	  }else {
	    histout[i]->Draw("SAME");
	  }
	}
	for (int i=0;i<n_layouts;i++){
	  if (i>=n_layouts/2){
	    leg[1]->AddEntry(histout[i],layout[i].c_str(),"l");
	  }else{ leg[0]->AddEntry(histout[i],layout[i].c_str(),"l");}
	}
	  leg[0]->Draw();
	  leg[1]->Draw();
	  	c1->SaveAs(("PhiRMS_vs_config_"+studypath+".pdf").c_str());


	return;


       
	return;
	
	
	// fill histograms
	for (int i = 0; i < n_datasets; i++){
	  
	  for (int j = 0; j < n_layouts; ++j) {

	    int n_entries = tuple[i][j]->GetEntries();
	    tuple[i][j]->SetBranchAddress("sigma_z2",&sigma_z);
	    tuple[i][j]->SetBranchAddress("sigma_phi2",&sigma_phi);
	    tuple[i][j]->SetBranchAddress("sigma_r2",&sigma_r);
	    tuple[i][j]->Draw(("(z2-z2t)*1e4:(phi2-phi2t)*1e4>>hDeltaPhiDeltaZ" + subpath[i]+layout[j]).c_str());
	    tuple[i][j]->Draw(("(phi2-phi2t)*1e4>>hDeltaPhi" + subpath[i]+layout[j]).c_str());
	    tuple[i][j]->Draw(("(z2-z2t)*1e4>>hDeltaZ" + subpath[i]+layout[j]).c_str());
	    
	    for (int k = 0; k < n_entries;k++) {
	      tuple[i][j]->GetEvent(k);
	      hSigmaZ[i][j]->Fill(TMath::Sqrt(sigma_z)*1e4);
	      hSigmaPhi[i][j]->Fill(TMath::Sqrt(sigma_phi)*1e4);
	      hSigmaR[i][j]->Fill(TMath::Sqrt(sigma_r)*1e4);
	      hSigmaPhiSigmaZ[i][j]->Fill(TMath::Sqrt(sigma_z)*1e4,TMath::Sqrt(sigma_phi)*1e4);//2D histogram
		}
	  }
	
	  //create a canvas to draw all the 2D plots:
	  c0->Divide(4,2);
	  c0->SetTitle(" track fit minus true at ~30cm");
	  for (int j=0;j<n_layouts;j++){
	    c0->cd(j+1);
	    hSigmaPhiSigmaZ[i][j]->Draw("colz");
	  }
	  c0->SaveAs(("SigmaPhi_vs_SigmaZ_multi_"+subpath[i]+".pdf").c_str());

	  c0->SetTitle(" track fit minus true at ~30cm");
	  for (int j=0;j<n_layouts;j++){
	    c0->cd(j+1);
	    hDeltaPhiDeltaZ[i][j]->Draw("colz");
	  }
	  c0->SaveAs(("DeltaPhi_vs_DeltaZ_multi_"+subpath[i]+".pdf").c_str());
	  for (int j=0;j<n_layouts;j++){
	  c0->cd(j+1);
	  hDeltaPhi[i][j]->Draw();
	}
	  c0->SaveAs(("DeltaPhi_multi_"+subpath[i]+".pdf").c_str());
	
	for (int j=0;j<n_layouts;j++){
	  c0->cd(j+1);
	  hDeltaZ[i][j]->Draw();
	}
	c0->SaveAs(("DeltaZ_multi_"+subpath[i]+".pdf").c_str());
	// create canvas and legend
	TCanvas *c = new TCanvas("c1","c1",1000,800);
	TLegend *legend;

	// draw and format Z histogram
	hSigmaZ[i][0]->Draw();
	for (int j = 1; j < n_layouts; j++) {	
		hSigmaZ[i][j]->Draw("same");
	}

	legend = new TLegend(0.11,0.55,0.25,0.89);

	for (int j = 0; j < n_layouts; j++) {
		legend->AddEntry(hSigmaZ[i][j],layout[j].c_str(),"l");
	}

	legend->Draw();
	c->SaveAs(("SigmaZ_vs_config_"+subpath[i]+".pdf").c_str());
	delete legend;

	// draw and format PHI histogram
	hSigmaPhi[i][0]->Draw();

	for (int j = 1; j < n_layouts;j++) {		
		hSigmaPhi[i][j]->Draw("same");
	}

	legend = new TLegend(0.11,0.55,0.25,0.89);

	for (int j = 0; j < n_layouts; j++) {
		legend->AddEntry(hSigmaPhi[i][j],layout[i].c_str(),"l");
	}

	legend->Draw();
	c->SaveAs("SigmaPHI_vs_config.png");
	delete legend;

	// plot RMS and error for Z
	double x[n_layouts], y[n_layouts], ex[n_layouts], ey[n_layouts];

	//you can also do this with a histogram.  Make one bin per layout, and use a fill that also takes an error value (or setbincontent, setbinerror):


	// a short routine to compile a set of means from an array of histograms:
	TH1F *hRMSZ=new TH1F("hRMSDelta","RMS values of #Deltasigma_z;layout;RMS (um)",n_layouts,-0.5,n_layouts-0.5);//I do the -0.5 to make sure the integers are centered in the bins.  Mostly a style thing.
	TH1F **histarray=hDeltaZ[i];
	TH1F *histout=hRMSZ;
	histout->SetCanExtend(TH1::kAllAxes); //lets us keep adding bins even if we go over the number we defined.
	//histout->SetLine
	//	histout->SetOptStat(0);
	for (int j=0;j<n_layouts;j++){
	  int bin=histout->Fill(layout[j].c_str(),histarray[j]->GetRMS());
	  histout->SetBinError(bin,histarray[j]->GetRMSError());
	}
	histout->Draw();
	c->SaveAs(("ZRMS_vs_config_"+subpath[i]+".pdf").c_str());


	TH1F *hRMSPhi=new TH1F("hRMSPhi","RMS values of #Deltasigma_#phi;layout;RMS (um)",n_layouts,-0.5,n_layouts-0.5);//I do the -0.5 to make sure the integers are centered in the bins.  Mostly a style thing.
	histarray=hDeltaPhi[i];
        histout=hRMSPhi;
	histout->SetCanExtend(TH1::kAllAxes); //lets us keep adding bins even if we go over the number we defined.
	//histout->SetLine
	//histout->SetOptStat(0);
	for (int j=0;j<n_layouts;j++){
	  int bin=histout->Fill(layout[j].c_str(),histarray[j]->GetRMS());
	  histout->SetBinError(bin,histarray[j]->GetRMSError());
	}
	histout->Draw();
	c->SaveAs(("PhiRMS_vs_config_"+subpath[i]+".pdf").c_str());

	
	TGraphErrors *gRMSZ = new TGraphErrors(n_layouts,x,y,ex,ey);
	
	hRMSZ->Draw();
	c->SaveAs(("ZRMS_vs_config_"+subpath[i]+".pdf").c_str());
	for (int j = 0; j < n_layouts; ++j) {
		cout << layout[j] << "\t" << hSigmaZ[i][j]->GetRMS() << "\t" << hSigmaZ[i][j]->GetRMSError() << endl;
	}
	}
	return;

	/*
	// close files so we can open new ones
	for (int i = 0; i < n_layouts; ++i) {
		fin[i][j]->Close();
		delete fin[i];
	}
	*/

	/*
	// prove INTT layout by plotting error in Z versus R
	float ez, r;
	TH2F *hZResVsR[n_layouts];
	TCanvas *c2=new TCanvas("c2","c2",1000,800);
	c2->Divide(4,2);
	for (int i = 0; i < n_layouts; ++i) {
		c2->cd(i+1);		
		fin[i] = new TFile((path + "G4_sPHENIX_" + subpath[i]+ layout[j] + ".root_g4svtx_eval.root").c_str(),"READ");
		if (fin[i]->IsZombie()) continue;

		fin[i]->GetObject("ntp_cluster",tuple[i]);
		hZResVsR[i] = new TH2F(("hZResVsR" + subpath[i]+ layout[j]).c_str(),(layout[i]+";R (cm);Z hitres (um?)").c_str(),60,0,16,10,0,0.6);
		tuple[i]->Draw(("ez:r>>hZResVsR"+layout[i]).c_str());//,"1","colz");
		hZResVsR[i]->Draw("colz");
		//	c->SaveAs(("Zerror_vs_R_" + layout[i] + ".png").c_str());
	}
	c2->SaveAs("ZRes_vs_R_multi.pdf");

	return;
	// close files
	for (int i = 0; i < n_layouts; ++i) {
		fin[i][j]->Close();
		delete fin[i];
	}
	// delete variables
	// when histograms are deleted something crashes, so for now there is a memory leak
	delete c;
	delete legend;

	gSystem->Exit(0);
	*/
}
