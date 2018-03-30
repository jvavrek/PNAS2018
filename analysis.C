//---------------------------------------------------------------------------------------
// Analysis code for HVRL 2017 campaign
// Jayson Vavrek, MIT, 2017

#include "TH1.h"
#include "TF1.h"
#include "TKDE.h"
#include "TString.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TLine.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TSpectrum.h"
#include "TLatex.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
using std::vector;
using std::cout;
using std::endl;

//---------------------------------------------------------------------------------------
// Data(, data, data! I can't make bricks without clay!)
//---------------------------------------------------------------------------------------

TTree *dataTree = new TTree("dataTree","tree of concise HVRL data");

vector<double> dead_time_fraction0;
vector<double> dead_time_fraction2;
vector<double> dead_time_fraction3;
vector<double> real_time;
vector<double> current;
vector<TString> lyso_filenames;
int nRuns;

// hacky way to actually create the data vectors at initialization: hide the function call
// inside a dummy function call, since ROOT allows variable defs but not function calls in
// this scope; anyway, it's often easier to work with std::vector's than TTree's, so convert
// the data format here
//
// Eventually might be nice to move to TTree structure, though...
bool createDataVectors(TString filename = "conciseData.txt")
{
	dataTree->ReadFile(filename,"run/I:time/D:DTfrac0/D:DTfrac2/D:DTfrac3/D:current/D:lysoFile/C",',');
	nRuns = dataTree->GetEntries();
	if (nRuns == 0) {
		cout << "Error reading " << filename << endl;
		cout << "Exiting..." << endl;
		return false;
	}
	else {
		TLeaf *run_leaf  = dataTree->GetBranch("run")->GetLeaf("run");
		TLeaf *time_leaf = dataTree->GetBranch("time")->GetLeaf("time");
		TLeaf *curr_leaf = dataTree->GetBranch("current")->GetLeaf("current");
		TLeaf *dtf0_leaf = dataTree->GetBranch("DTfrac0")->GetLeaf("DTfrac0");
		TLeaf *dtf2_leaf = dataTree->GetBranch("DTfrac2")->GetLeaf("DTfrac2");
		TLeaf *dtf3_leaf = dataTree->GetBranch("DTfrac3")->GetLeaf("DTfrac3");
		TLeaf *lyso_leaf = dataTree->GetBranch("lysoFile")->GetLeaf("lysoFile");

		for (size_t i = 0; i < nRuns; ++i) {
			dataTree->GetEntry(i);
			real_time.push_back( 60.0*time_leaf->GetValue() );
			dead_time_fraction0.push_back( dtf0_leaf->GetValue() );
			dead_time_fraction2.push_back( dtf2_leaf->GetValue() );
			dead_time_fraction3.push_back( dtf3_leaf->GetValue() );
			current.push_back( curr_leaf->GetValue() );
			lyso_filenames.push_back( lyso_leaf->GetValuePointer() );
		}
	}
	cout << "Attempting to read data for " << nRuns << " runs..." << endl;
	return true;
}
bool done = createDataVectors();

vector<double> compute_live_time(int det)
{
	vector<double> tl;
	vector<double> dtf;
	if (det == 0) dtf = dead_time_fraction0;
	else if (det == 2) dtf = dead_time_fraction2;
	else if (det == 3) dtf = dead_time_fraction3; 
	for (size_t i = 0; i < nRuns; ++i)
		tl.push_back(real_time[i]*(1.0-dtf[i]));

	return tl;
}
const vector<double> live_time0 = compute_live_time(0);
const vector<double> live_time2 = compute_live_time(2);
const vector<double> live_time3 = compute_live_time(3);


//---------------------------------------------------------------------------------------
// Functions to convert raw data into useful histograms
//---------------------------------------------------------------------------------------

TH1D *plot_tka(TString filename, int nbins = 32768, bool verbose = true)
{
	TString runname = filename;
	runname.ReplaceAll(".TKA","");
	runname.ReplaceAll("hpge/","");
	runname.ReplaceAll("lyso/","");

	ifstream file(filename);
	if (file) {
		// initialize the histogram, correcting for the first two entries being
		// overwritten with time information; keep the number of bins a power of
		// two so we can rebin by factors of two later to polish the data
		TH1D *h = new TH1D(runname, runname, nbins, 0, nbins);

		string line;
		int counter = 0;
		while (getline(file,line))
		{
			++counter;
			int val = atoi(line.c_str());
			if (counter <= 2) val = 0;
			h->SetBinContent(counter,val);
		}
		h->Sumw2();
		return h;
	}
	else {
		if (verbose) cout << "Warning! File " << filename << " not found." << endl;
		return 0;
	}
}


size_t nDetectors(size_t i) {
	if (i >= 71) return 3;
	else return 1;
}


// assumes Cs-137 and Co-60
// takes adc channels in ascending order of energy
// hardcoded adc's unfortunately
TF1 *calibrateHPGe(int detNumber, TString opt = "goff", int adc0 = 0, int adc1 = 0, int adc2 = 0){
	TString filename;
	if (detNumber==0){
		filename = "hpge/cal_sept12_morning_det0+3_Co60+Cs137_0.TKA";
		adc0=6263; adc1=11090; adc2=12602;
	}
	else if (detNumber==2){
		filename = "hpge/cal_sept12_morning_det2_Co60+Cs137_2.TKA";
		adc0=6366; adc1=11308; adc2=12865;
	}
	else if (detNumber==3){
		filename = "hpge/cal_sept12_morning_det0+3_Co60+Cs137_3.TKA";
		adc0=6917; adc1=12263; adc2=13920;
	}

	TH1D *h = plot_tka(filename);
	h->Rebin(8);
	if (opt != "goff") h->Draw();

	TF1 *fit;
	if (adc0>0 && adc1>0 && adc2>0){
		if (opt != "goff") {
			TCanvas *c2 = new TCanvas("c2","c2");
			c2->cd();
		}
		TGraph *g = new TGraph(3);
		g->SetPoint(0,adc0,0.662);
		g->SetPoint(1,adc1,1.173);
		g->SetPoint(2,adc2,1.333);
		g->SetMarkerStyle(2);
		g->SetMarkerColor(kBlack);
		g->SetMarkerSize(1.2);
		if (opt != "goff") g->Draw("EPA");
	
		fit = new TF1("fit","x++1",h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
		g->Fit("fit");

		if (opt != "goff") fit->Draw("same");
	}

	cout << "slope = " << fit->GetParameter(0) << endl;
	cout << "incpt = " << fit->GetParameter(1) << endl;

	return fit;
}


// the "standard" linear calibrations
TF1 *standard_cal(TString det = "HPGe") {
	TF1 *f = new TF1("f", "x++1", 0, 30000);

	double slope, incpt;
	// the "standard" calibration for HPGe
	if (det == "HPGe"){
		slope = 0.0959e-3;
		incpt = -1.834e-3;
	}
	// the "standard" calibration for LYSO, using 0=0 and the 662 keV peak at ch 2887 in 08/02 5pm run
	else if (det == "LYSO"){
		slope = 0.2293e-3;
		incpt = 0.0;
	}
	else if (det == "LaBr"){
		TF1 *fitLaBr = LaBr_Th232_cal();
		slope = fitLaBr->GetParameter(1);
		incpt = fitLaBr->GetParameter(0);
	}
	else if (det == "HPGe0"){
		TF1 *fit = calibrateHPGe(0);
		slope = fit->GetParameter(0);
		incpt = fit->GetParameter(1);
	}
	else if (det == "HPGe2"){
		TF1 *fit = calibrateHPGe(2);
		slope = fit->GetParameter(0);
		incpt = fit->GetParameter(1);
	}
	else if (det == "HPGe3"){
		TF1 *fit = calibrateHPGe(3);
		slope = fit->GetParameter(0);
		incpt = fit->GetParameter(1);
	}
	f->SetParameters(slope, incpt);
	return f;
}

TGraph *vector_plot(vector<double> x, vector<double> y) {
	int n = x.size();
	TGraph *g = new TGraphErrors(n, &(x[0]), &(y[0]));
	return g;
}

// Takes an ostensibly-calibrated HPGe spectrum h and recalibrates it assuming there's
// a 511 keV, 1001 keV, and 2212 keV line.
TF1 *autoCalibration(TH1 *h){
	vector<double> evec, fvec;
	evec.push_back(0.511);
	evec.push_back(1.001);
	evec.push_back(2.212);
	double etol = 0.100;

	for (size_t i = 0; i < evec.size(); ++i){
		h->GetXaxis()->SetRangeUser(evec[i]-etol, evec[i]+etol);
		TSpectrum *s = new TSpectrum(1); // find the strongest peak
		s->Search(h,2,"goff");
		fvec.push_back((s->GetPositionX())[0]);
	}

	TGraph *newCal = vector_plot(fvec, evec);
	TF1 *newFit = new TF1("newFit","x++1",0.0,3.0);
	newCal->Fit("newFit");
	return newFit;
}

// function to apply the calibrations to histograms
// works for calibrating both TH1D spectra and TH2F PSD plots through templating
template <typename T>
T *calibrate(T *h, TF1 *calibration_fit) {
	int nbinsX = h->GetNbinsX();
	int nbinsY = h->GetNbinsY();

	double xmin = h->GetXaxis()->GetXmin();
	double xmax = h->GetXaxis()->GetXmax();
	double ymin = h->GetYaxis()->GetXmin();
	double ymax = h->GetYaxis()->GetXmax(); // bizarre ROOT convention...

	// calibrated energy limits
	double Emin = calibration_fit->Eval(xmin);
	double Emax = calibration_fit->Eval(xmax);
	// calibrated PSD limits; note: not actually calibrating here
	double pmin = ymin;
	double pmax = ymax;

	TString runname = h->GetName();
	runname += "_cal";

	// create the new histogram
	T *hc = h->Clone();
	hc->Reset();
	hc->SetName(runname);
	hc->GetXaxis()->SetTitle("energy #it{E} [MeV]");
	hc->GetYaxis()->SetTitle(h->GetYaxis()->GetTitle());
	hc->SetTitle(h->GetTitle());

	hc->GetXaxis()->SetLabelSize(0.04);
	hc->GetYaxis()->SetLabelSize(0.04);
	hc->GetXaxis()->SetTitleSize(0.04);
	hc->GetYaxis()->SetTitleSize(0.04);
	hc->GetXaxis()->SetTitleOffset(1.0);
	hc->GetYaxis()->SetTitleOffset(1.0);


	// set new histogram limits
	if (nbinsY == 1) {
		hc->SetBins(h->GetNbinsX(), Emin, Emax);
	}
	else if (nbinsY > 1) {
		hc->SetBins(h->GetNbinsX(), Emin, Emax, h->GetNbinsY(), pmin, pmax);
	}

	for (size_t i = 1; i <= nbinsX*nbinsY; ++i) {
		hc->SetBinContent(i, h->GetBinContent(i));
		hc->SetBinError(  i, h->GetBinError(i));
	}

	return hc;
}


template <typename T>
T *applyAutoCalibration(T *h) {
	return calibrate(h, autoCalibration(h));
}


// Make a copy of a vector of pointers with new pointers to the same objects.
// Specifically used for TH1D* here but could be generalized in the future.
// Useful for copying spectra to new histograms safely (without aliasing).
vector< TH1D* > deepVectorCopy(vector< TH1D* > v) {
	vector< TH1D* > vnew = v;
	for (size_t i = 0; i < v.size(); ++i){
		if (v[i] == NULL) vnew[i] = NULL;
		else vnew[i] = v[i]->Clone();
	}
	return vnew;
} 


// this is the vector of the data in its rawest form (at least after conversion
// to a TH1D from raw text files)
vector< TH1D* > build_hv_uncal(int detNumber) {
	if (detNumber == 1) {cout << "Error! No Det1. Aborting..." << endl; exit(1);}

	vector< TH1D* > v;
	bool vOpt = true;
	for (size_t i = 0; i < nRuns; ++i)
	{
		TString fname = Form("hpge/run%i",i);
		size_t nDets = nDetectors(i);
		if ((nDets > 1 && detNumber == 3) || detNumber != 3) fname.Append(Form("_%i",detNumber));
		fname.Append(".TKA");
		if (detNumber == 2 || (detNumber != 3 && i <= 70)) vOpt = false; // manual override
		TH1D *h = plot_tka(fname, 32768, vOpt);
		if (h != NULL) v.push_back(h);
		else v.push_back(0);
	}

	int nonZeroCounter = 0;
	for (size_t i = 0; i < v.size(); ++i)
		if (v[i] != 0) ++nonZeroCounter;

	cout << "Found " << nonZeroCounter << " / " << v.size() << " spectra for Detector " << detNumber << endl;
	return v;
}
const vector< TH1D* > hv_uncal3 = build_hv_uncal(3);
const vector< TH1D* > hv_uncal2 = build_hv_uncal(2);
const vector< TH1D* > hv_uncal0 = build_hv_uncal(0);


vector<double> getLiveTimeVector(int detNumber){
	if (detNumber == 0) return live_time0;
	else if (detNumber == 2) return live_time2;
	else if (detNumber == 3) return live_time3;
	else {cout << "getLiveTimeVector detNumber " << detNumber << " error!" << endl; exit(1);}
}


vector< TH1D* > getUncalHistVector(int detNumber){
	if (detNumber == 0) return deepVectorCopy(hv_uncal0);
	else if (detNumber == 2) return deepVectorCopy(hv_uncal2);
	else if (detNumber == 3) return deepVectorCopy(hv_uncal3);
	else {cout << "getUncalHistVector detNumber " << detNumber << " error!" << endl; exit(1);}
}


// and this is a vector using the "standard" calibration taken on day 1 of experiments
// this could be generalized...
vector< TH1D* > build_hv_cal(int detNumber) {
	cout << endl << "Calibrating detector " << detNumber << "..." << endl;
	vector< TH1D* > hv_uncal = getUncalHistVector(detNumber);

	TF1 *fit  = NULL; TF1 *fit0 = NULL; TF1 *fit2 = NULL; TF1 *fit3 = NULL;
	if (detNumber == 3){
		fit  = standard_cal("HPGe");
		fit3 = standard_cal("HPGe3");
	}
	else if (detNumber == 2) fit2 = standard_cal("HPGe2");
	else if (detNumber == 0) fit0 = standard_cal("HPGe0");
	else {cout << "Error! no calibration for Det " << detNumber << endl; exit(1);}

	vector< TH1D* > v;
	for (size_t i = 0; i < hv_uncal.size(); ++i)
	{
		TH1D *h = hv_uncal[i];
		TF1 *theFit;
		if (h != NULL){
			// adaptive calibration
			if (detNumber == 3){
				if (i > 70) theFit = fit;
				else theFit = fit3;
			} 
			else if (detNumber == 2) theFit = fit2;
			else if (detNumber == 0) theFit = fit0;

			TH1D *hnew = calibrate(h, theFit);
	
			// set up at least some initial axes labels
			hnew->GetXaxis()->SetTitle("energy #it{E} [MeV]");
			hnew->GetYaxis()->SetTitle("counts per bin");
			hnew->SetTitle(hnew->GetName());

			v.push_back(hnew);
		}
		else v.push_back(0);
	}
	return v;
}
const vector< TH1D* > hv_cal3 = build_hv_cal(3);
const vector< TH1D* > hv_cal2 = build_hv_cal(2);
const vector< TH1D* > hv_cal0 = build_hv_cal(0);


vector< TH1D* > getCalHistVector(int detNumber){
	if (detNumber == 0) return deepVectorCopy(hv_cal0);
	else if (detNumber == 2) return deepVectorCopy(hv_cal2);
	else if (detNumber == 3) return deepVectorCopy(hv_cal3);
	else {cout << "getCalHistVector detNumber error!" << endl; exit(1);}
}


// It's also useful to rebin all the histograms in a vector as needed.
// This should only be called on the local copies of data so we can keep the const.
void VRebin(vector< TH1D* > v, int rebinFactor = 8) {
	for (size_t i = 0; i < v.size(); ++i)
		if (v[i] != NULL)
			v[i]->Rebin(rebinFactor);
}


// A general-purpose histogram sum routine with different normalization options.
// Takes a vector v of TH1D*, and appropriately sums between iMin and iMax under normalization option.
// Currently only supports one detector at at time.
TH1D *norm_sum(TString name, vector< TH1D* > v, int detNumber, size_t iMin, size_t iMax, TString option = "lc") {
	// currently accepts only live charge, live time, or none
	// for LaBr it might be useful to expand to real time and charge
	if (option != "lc" && option != "lt" && option != "n") {
		cout << "Error! Invalid option " << option << ". Aborting..." << endl; exit(1);
	}

	TH1D *hnorm = v[iMin]->Clone(name);
	hnorm->SetTitle(name);
	hnorm->SetName(name);
	hnorm->Reset();
	hnorm->Sumw2();

	for (size_t i = iMin; i <= iMax; ++i){
		if (is_problematic_NRF(i)) continue;
		hnorm->Add(v[i]);
	}

	double norm_factor = 0.0;
	vector<double> lt = getLiveTimeVector(detNumber);
	for (size_t i = iMin; i <= iMax; ++i){
		if (is_problematic_NRF(i)) continue;

		if (option == "lt") norm_factor += lt[i];
		else if (option == "lc") norm_factor += (lt[i] * current[i]);
		else continue;
	}

	if (option != "n") hnorm->Scale(1.0/norm_factor);
	return hnorm;
}


//---------------------------------------------------------------------------------------
// Fitting functions
//---------------------------------------------------------------------------------------

// a more intelligent guess at the peak height (constant) parameter
// find the maximum bin content difference within +/- nbinsLR of the peak
double constant_search(TH1 *h, double peak_energy, int nbinsLR = 3)
{
	int peak_bin = h->FindBin(peak_energy);
	double constant = 0.0;
	for (int i = -nbinsLR; i <= nbinsLR; ++i)
	{
		for (int j = i+1; j <= nbinsLR; ++j)
		{
			double dy = fabs(h->GetBinContent(peak_bin + i) - h->GetBinContent(peak_bin + j));
			//cout << peak_bin + i << " " << peak_bin + j << " " << dy << endl;
			if (dy > constant) constant = dy;
		}
	}
	if (constant == 0.0) {cout << "Error: constant still 0.0" << endl; exit(1);}
	//cout << "constant = " << constant << endl;
	return constant;
}


// function to generate a fit string that uses the non-default form of the Gaussian
TString generate_function_string(int nPeaks = 1, TString bkgd_string = "expo(0)"){
	TString result = bkgd_string;
	for (size_t i = 0; i < nPeaks; ++i){
		size_t parNo = 2 + 3*i;
		result += Form("+[%i]/[%i]/sqrt(2.0*TMath::Pi())*exp(-(x-[%i])*(x-[%i])/(2.0*[%i]*[%i]))",
			parNo, parNo+2, parNo+1, parNo+1, parNo+2, parNo+2);
	}
	return result;
}


// Single Gaussian plus linear background
TF1 *fit_gaus_pol1(TString name, TH1 *h, double peak_energy, double region_width = 0.020, double sigma = -1)
{
	TF1 *fit = new TF1(name, "pol1(0)+gaus(2)", peak_energy - region_width/2.0, peak_energy + region_width/2.0);

	// starting estimates of peak parameters
	if (sigma < 0) sigma = 0.25*region_width;
	double e_lo = peak_energy - 2.0*sigma;
	double e_hi = peak_energy + 2.0*sigma;
	double c_lo = h->GetBinContent(h->FindBin(e_lo));
	double c_hi = h->GetBinContent(h->FindBin(e_hi));
	double slope = (c_hi-c_lo)/(e_hi-e_lo);
	double incpt = c_hi - slope * e_hi;
	//double constant = h->GetBinContent(h->FindBin(peak_energy)) - h->GetBinContent(h->FindBin(e_hi)); // OLD
	double constant = constant_search(h, peak_energy);
	

	fit->SetParameters(incpt, slope, constant, peak_energy, sigma);
	// plus some extra bounds to stabilize the more complicated fit
	fit->SetParLimits(2,1e-2*constant,1e2*constant);
	fit->SetParLimits(3,peak_energy-2*sigma,peak_energy+2*sigma);
	fit->SetParLimits(4,sigma/10.0,10.0*sigma);

	fit->SetLineColor(h->GetLineColor());
	fit->SetNpx(1000);

	h->Fit(fit, "BRM+");

  return fit;
}


// Single Gaussian plus constant background
TF1 *fit_gaus_pol0(TString name, TH1 *h, double peak_energy, double region_width = 0.020, double sigma = -1)
{
	TF1 *fit = new TF1(name, "pol0(0)+gaus(1)", peak_energy - region_width/2.0, peak_energy + region_width/2.0);

	// starting estimates of peak parameters
	if (sigma < 0) sigma = 0.25*region_width;
	double e_lo = peak_energy - 2.0*sigma;
	double e_hi = peak_energy + 2.0*sigma;
	double c_lo = h->GetBinContent(h->FindBin(e_lo));
	double c_hi = h->GetBinContent(h->FindBin(e_hi));
	double c_avg = 0.5*(c_lo+c_hi);
	//double constant = h->GetBinContent(h->FindBin(peak_energy)) - h->GetBinContent(h->FindBin(e_hi)); // OLD
	double constant = constant_search(h, peak_energy);

	fit->SetParameters(c_avg, constant, peak_energy, sigma);
	// plus some extra bounds to stabilize the more complicated fit
	fit->SetParLimits(1,1e-3*constant,1e2*constant);
	fit->SetParLimits(2,peak_energy-2*sigma,peak_energy+2*sigma);
	fit->SetParLimits(3,sigma/2.0,2.0*sigma);

	fit->SetLineColor(h->GetLineColor());
	fit->SetNpx(1000);

	h->Fit(fit, "BRM+");

  return fit;
}


// Single Gaussian plus exponential background
TF1 *fit_gaus_expo(TString name, TH1D *h, double peak_energy, double region_width = 0.020, double sigma = -1)
{
	TF1 *fit = new TF1(name, "expo(0)+gaus(2)", peak_energy - region_width/2.0, peak_energy + region_width/2.0);

	// starting estimates of peak parameters
	if (sigma < 0) sigma = 0.25*region_width;
	double e_lo = peak_energy - 2.0*sigma;
	double e_hi = peak_energy + 2.0*sigma;
	double c_lo = h->GetBinContent(h->FindBin(e_lo));
	double c_hi = h->GetBinContent(h->FindBin(e_hi));
	double ex_decay = TMath::Log(c_hi/c_lo)/(e_hi-e_lo);
	double ex_constant = TMath::Log(c_hi)/(ex_decay * e_hi);
	//double constant = h->GetBinContent(h->FindBin(peak_energy)) - h->GetBinContent(h->FindBin(e_hi)); // OLD
	double constant = constant_search(h, peak_energy);

	fit->SetParameters(ex_constant, ex_decay, constant, peak_energy, sigma);
	// plus some extra bounds to stabilize the more complicated fit
	fit->SetParLimits(2,1e-3*constant,1e2*constant);
	fit->SetParLimits(3,peak_energy-2*sigma,peak_energy+2*sigma);
	fit->SetParLimits(4,sigma/2.0,2.0*sigma);

	fit->SetLineColor(h->GetLineColor());
	fit->SetNpx(1000);

	h->Fit(fit, "BRM+");

	return fit;
}


// Two Gaussians plus linear background
// See also fit_uranium_multiplet
TF1 *fit_doublet_pol1(TH1 *h, double e1, double e2, double region_width = 0.010)
{
	TF1 *doublet_fit = new TF1("doublet_fit", "pol1(0)+gaus(2)+gaus(5)", e1 - region_width/2.0, e2 + region_width/2.0);

	// starting estimates of peak parameters
	double sigma = 0.25*region_width;
	double e_lo = e1 - 2.0*sigma;
	double e_hi = e2 + 2.0*sigma;
	double c_lo = h->GetBinContent(h->FindBin(e_lo));
	double c_hi = h->GetBinContent(h->FindBin(e_hi));
	double slope = (c_hi-c_lo)/(e_hi-e_lo);
	double incpt = c_hi - slope * e_hi;
	double constant1 = h->GetBinContent(h->FindBin(e1)) - h->GetBinContent(h->FindBin(e_lo)); // could add constant_search here
	double constant2 = h->GetBinContent(h->FindBin(e2)) - h->GetBinContent(h->FindBin(e_hi));
	if (constant1 < 0 && constant2 > 0) constant1 = constant2;
	if (constant2 < 0 && constant1 > 0) constant2 = constant1;

	doublet_fit->SetParameters(incpt, slope, constant1, e1, sigma, constant2, e2, sigma);
	// plus some extra bounds to stabilize the more complicated fit
	doublet_fit->SetParLimits(2,1e-3*constant1,1e2*constant1);
	doublet_fit->SetParLimits(3,e1-2*sigma,e1+2*sigma);
	doublet_fit->SetParLimits(4,sigma/2.0,2.0*sigma);
	doublet_fit->SetParLimits(5,1e-3*constant2,1e2*constant2);
	doublet_fit->SetParLimits(6,e2-2*sigma,e2+2*sigma);
	doublet_fit->SetParLimits(7,sigma/2.0,2.0*sigma);

	doublet_fit->SetLineColor(h->GetLineColor());

	h->Fit(doublet_fit, "BRM");
	doublet_fit->SetNpx(1000);

	return doublet_fit;
}


// Three Gaussians plus linear background
// See also fit_uranium_multiplet
TF1 *fit_triplet_pol1(TH1 *h, double e1, double e2, double e3, double region_width = 0.020)
{
	TF1 *triplet_fit = new TF1("triplet_fit", "pol1(0)+gaus(2)+gaus(5)+gaus(8)", e1 - region_width/2.0, e3 + region_width/2.0);

	// starting estimates of peak parameters
	double sigma = 0.25*region_width;
	double e_lo = e1 - 2.0*sigma;
	double e_hi = e3 + 2.0*sigma;
	double c_lo = h->GetBinContent(h->FindBin(e_lo));
	double c_hi = h->GetBinContent(h->FindBin(e_hi));
	double slope = (c_hi-c_lo)/(e_hi-e_lo);
	double incpt = c_hi - slope * e_hi;
	double constant1 = h->GetBinContent(h->FindBin(e1)) - h->GetBinContent(h->FindBin(e_lo)); // could add constant_search here
	double constant2 = h->GetBinContent(h->FindBin(e3)) - h->GetBinContent(h->FindBin(e_hi));
	if (constant1 < 0 && constant2 > 0) constant1 = constant2;
	if (constant2 < 0 && constant1 > 0) constant2 = constant1;
	if (constant3 < 0 && constant1 > 0) constant3 = constant1;

	triplet_fit->SetParameters(incpt, slope, constant1, e1, sigma, constant2, e2, sigma, constant3, e3, sigma);
	// plus some extra bounds to stabilize the more complicated fit
	triplet_fit->SetParLimits(2,1e-3*constant1,1e2*constant1);
	triplet_fit->SetParLimits(3,e1-2*sigma,e1+2*sigma);
	triplet_fit->SetParLimits(4,sigma/2.0,2.0*sigma);
	triplet_fit->SetParLimits(5,1e-3*constant2,1e2*constant2);
	triplet_fit->SetParLimits(6,e2-2*sigma,e2+2*sigma);
	triplet_fit->SetParLimits(7,sigma/2.0,2.0*sigma);
	triplet_fit->SetParLimits(8,1e-3*constant3,1e2*constant3);
	triplet_fit->SetParLimits(9,e3-2*sigma,e3+2*sigma);
	triplet_fit->SetParLimits(10,sigma/2.0,2.0*sigma);

	triplet_fit->SetLineColor(h->GetLineColor());

	h->Fit(triplet_fit, "BRM");
	triplet_fit->SetNpx(1000);

	return triplet_fit;
}


// Two Gaussians plus exponential background
TF1 *fit_doublet_expo(TH1D *h, double e1, double e2, double region_width = 0.010)
{
	TF1 *doublet_fit = new TF1("doublet_fit", "expo(0)+gaus(2)+gaus(5)", e1 - region_width/2.0, e2 + region_width/2.0);

	// starting estimates of peak parameters
	double sigma = 0.25*region_width;
	double e_lo = e1 - 2.0*sigma;
	double e_hi = e2 + 2.0*sigma;
	double c_lo = h->GetBinContent(h->FindBin(e_lo));
	double c_hi = h->GetBinContent(h->FindBin(e_hi));
	double constant1 = h->GetBinContent(h->FindBin(e1)) - h->GetBinContent(h->FindBin(e_lo)); // could add constant_search here
	double constant2 = h->GetBinContent(h->FindBin(e2)) - h->GetBinContent(h->FindBin(e_hi));
	double ex_decay = TMath::Log(c_hi/c_lo)/(e_hi-e_lo);
	double ex_constant = TMath::Log(c_hi)/(ex_decay * e_hi);

	doublet_fit->SetParameters(ex_constant, ex_decay, constant1, e1, sigma, constant2, e2, sigma);
	// plus some extra bounds to stabilize the more complicated fit
	doublet_fit->SetParLimits(2,1e-3*constant1,1e2*constant1);
	doublet_fit->SetParLimits(3,e1-2*sigma,e1+2*sigma);
	doublet_fit->SetParLimits(4,sigma/2.0,2.0*sigma);
	doublet_fit->SetParLimits(5,1e-3*constant2,1e2*constant2);
	doublet_fit->SetParLimits(6,e2-2*sigma,e2+2*sigma);
	doublet_fit->SetParLimits(7,sigma/2.0,2.0*sigma);

	doublet_fit->SetLineColor(h->GetLineColor());
	doublet_fit->SetNpx(1000);

	h->Fit(doublet_fit, "BRM");

	return doublet_fit;
}


// Main function to fit the U-238 NRF peaks (2.176, 2.209, 2.245 MeV multiplet)
// peakCode values:
//  - -2: (experimental) -1 with doublet fit near the Al-27 peak at 2.212 MeV
//  - -1: 2.176 MeV and its branch plus 2.245 MeV and the branch from 2.209 MeV
//  -  3: 2.176, 2.209, 2.245 MeV (default)
//  -  5: 2.176, 2.209, 2.245 MeV + branches of former two
//  -  6: 2.176, 2.209, 2.245 MeV + branches of all three
TF1 *fit_uranium_multiplet(TH1 *h, int peakCode = 3, double emin = NULL, double emax = NULL, double eShift = 0.0)
{
	if (peakCode != -1 && peakCode != -2 && peakCode != 3 && peakCode != 5 && peakCode != 6)
	{
		cout << "Error! peakCode " << peakCode << " not supported. Exiting..." << endl;
		exit(1);
	}
	size_t nPeaksTmp;
	if (peakCode > 0) {nPeaksTmp = peakCode;}
	else {
		if (peakCode == -1) nPeaksTmp = 4;
		else if (peakCode == -2) nPeaksTmp = 8;
	}
	const size_t nPeaks = nPeaksTmp;
	const size_t nParameters = 2 + 3*nPeaks;

	TString fitname = "multiplet_fit_";
	fitname += h->GetName();

	// start by building vectors of initial peak energies
	vector<double> ePeak_vec;
	double branch_diff = 0.045;
	if (peakCode == -1){
		ePeak_vec.push_back(2.176);
		ePeak_vec.push_back(ePeak_vec[0]-branch_diff);
		ePeak_vec.push_back(2.245);
		ePeak_vec.push_back(2.164);
	}
	else if (peakCode == -2){
		ePeak_vec.push_back(2.176);
		ePeak_vec.push_back(ePeak_vec[0]-branch_diff);
		ePeak_vec.push_back(2.245);
		ePeak_vec.push_back(2.164);
		ePeak_vec.push_back(2.145);
		ePeak_vec.push_back(2.200);
		ePeak_vec.push_back(2.209);
		ePeak_vec.push_back(2.212);
	}
	else {
		ePeak_vec.push_back(2.176 + eShift);
		ePeak_vec.push_back(2.209 + eShift);
		ePeak_vec.push_back(2.2465 + eShift);
		if (nPeaks >= 5){
			ePeak_vec.push_back(ePeak_vec[0]-branch_diff);
			ePeak_vec.push_back(ePeak_vec[1]-branch_diff);
		}
		if (nPeaks == 6) ePeak_vec.push_back(ePeak_vec[2]-branch_diff+0.001);
	}
	double eRangeMin = (emin == NULL ? vector_min(ePeak_vec) - 0.010 : emin);
	double eRangeMax = (emax == NULL ? vector_max(ePeak_vec) + 0.020 : emax);

	// initialize the TF1 for the final fit
	TString function_string = generate_function_string(nPeaks);
	TF1 *multiplet_fit = new TF1(fitname,                // name
														   function_string.Data(), // fit function
														   eRangeMin, eRangeMax);  // range
	multiplet_fit->SetLineColor(h->GetLineColor());
	multiplet_fit->SetLineWidth(1);
	multiplet_fit->SetNpx(2000);

	// some HAX to make the fit parameters human-readable;
	for (size_t i = 0; i < nParameters; ++i){
		TString name;
		if (i == 0) name = "ExpConst";
		else if (i == 1) name = "ExpSlope";
		else
		{
			size_t j = (i-2)%3;
			if (j == 0) name = "Area";
			else if (j == 1) name = "Mean";
			else if (j == 2) name = "Sigma";
		}
		size_t k = ((i-2)/3)%nPeaks + 1;
		if (i > 1) name += Form("%i",k);
		multiplet_fit->SetParName(i,name);
	}

	// start with preliminary independent fits to each of the peaks
	// to get parameter estimates
	vector<TF1*> fit_vec;
	for (size_t i = 0; i < ePeak_vec.size(); ++i)
	{
		TF1 *fit_tmp = fit_gaus_pol1("fit_tmp", h, ePeak_vec[i], 0.015);
		fit_vec.push_back(fit_tmp);
	}
	cout << "created preliminary peak fits" << endl;

	// create a preliminary estimate of the background
	TF1* bkd_prelim = new TF1("bkd_prelim","expo",eRangeMin,eRangeMax);
	h->Fit("bkd_prelim","BRM");
	double expConst = bkd_prelim->GetParameter(0);
	double expSlope = bkd_prelim->GetParameter(1);
	multiplet_fit->SetParameter(0, expConst);
	multiplet_fit->SetParameter(1, expSlope);
	multiplet_fit->SetParLimits(0, 0.5*expConst, 2.0*expConst);
	multiplet_fit->SetParLimits(1, 2.0*expSlope, 0.5*expSlope);
	cout << "created preliminary background fit" << endl;

	for (size_t i = 0; i < ePeak_vec.size(); ++i){
		size_t j = 3*i + 2;
		TF1 *theFit;

		// manual override for closely-spaced peaks near 2.210 MeV
		// hopefully your appetite is good because this is a big plate of spaghetti
		if (peakCode > 0 && (i == 1 || i == 5)) theFit = fit_doublet_pol1(h, 2.2025 + eShift, 2.210 + eShift);
		else if (peakCode == -2 && i >= ePeak_vec.size()-2) theFit = fit_doublet_pol1(h, 2.2095, 2.2125, 0.008);
		else theFit = fit_vec[i];

		int parNoCorr = 0;
		if (peakCode > 0 && i == 1) parNoCorr = 3;
		else if (peakCode == -2 && i == 7) parNoCorr = 3;

		double constant = theFit->GetParameter(2 + parNoCorr);
		double mean     = theFit->GetParameter(3 + parNoCorr);
		double sigma    = theFit->GetParameter(4 + parNoCorr);

		// Here we need to convert from the {constant, mean, sigma} parameter set that the "gaus"
		// functions use to the {area, mean, sigma} parameter set in the custom gaussian we're
		// using to get around the correlated uncertainties issue.
		double area = TMath::Sqrt(2.0*TMath::Pi()) * constant * sigma;
		multiplet_fit->SetParameter(j, area);
		multiplet_fit->SetParameter(j+1, mean);
		multiplet_fit->SetParameter(j+2, sigma);

		multiplet_fit->SetParLimits(j,   0.10*area, 6.0*area);
		multiplet_fit->SetParLimits(j+1, mean - sigma/2.0, mean + sigma/2.0);
		multiplet_fit->SetParLimits(j+2, 0.30*sigma, 3.0*sigma);
	}
	
	h->Fit(fitname, "BRM");
	cout << "created final multiplet fit: " << multiplet_fit->GetName() << endl;
	cout << "chi2/ndf = " << multiplet_fit->GetChisquare()/multiplet_fit->GetNDF() << endl;
	cout << "parameter relative errors:" << endl;
	for (size_t i = 0; i < nParameters; ++i){
		cout << "  " << multiplet_fit->GetParError(i)/multiplet_fit->GetParameter(i) << endl;
	}

	cout << "peak resolutions (sigma):" << endl;
	for (size_t i = 0; i < ePeak_vec.size(); ++i){
		size_t j = 3*i + 2;
		double en = multiplet_fit->GetParameter(j+1);
		double resSigma = multiplet_fit->GetParameter(j+2);
		cout << "  " << en << " MeV: " << resSigma*100.0/en << " perCent" << endl;
	}
	cout << endl;

	return multiplet_fit;
}


// Extract the area and, optionally, the uncertainty of a gaussian from
// a single TF1 *myFit of a TH1 *h that is of the form "<2-par background> + gaus"
// after already having called the relevant fitting function.
double gaussian_area(TF1 *myFit, TH1 *hist, double &error, bool altForm = false){
	double h = myFit->GetParameter(2);
	double dh = myFit->GetParError(2);
	double s = myFit->GetParameter(4);
	double ds = myFit->GetParError(4);

	double counts = h * s * sqrt(2.0*TMath::Pi())/hist->GetBinWidth(1);
	double rel_err = sqrt(dh*dh/h/h + ds*ds/s/s);
	error = rel_err * counts;

	return counts;
}

// Similar function designed to store the gaus areas + errors from fit_uranium_multiplet.
// Note here the convention is "<2-par background> + gaus + ... + gaus".
// Optionally, set altForm = true if using the custom Gaussian function
vector<double> gaussian_areas(TF1 *myFit, TH1 *hist, vector<double> &errors, vector<double> &energies, bool altForm = false){
	vector<double> counts;
	const size_t nParameters = myFit->GetNpar();
	const size_t nPeaks = (nParameters - 2)/3;

	for (size_t i = 0; i < nPeaks; ++i)
	{
		size_t j = 2+3*i;
		double h = myFit->GetParameter(j)/hist->GetBinWidth(1);
		double e = myFit->GetParameter(j+1);
		double s = myFit->GetParameter(j+2);
		double dh = myFit->GetParError(j)/hist->GetBinWidth(1);
		double ds = myFit->GetParError(j+2);

		double c = h * s * sqrt(2.0*TMath::Pi())/hist->GetBinWidth(1);
		double dc = c * frac_ratio_error(h,dh,s,ds);

		// if altForm = true, h represents the area already; else it is the height
		if (altForm){
			counts.push_back(h);
			errors.push_back(dh);
		}
		else {
			counts.push_back(c);
			errors.push_back(dc);
		}
		energies.push_back(e);
	}

	return counts;
}


//---------------------------------------------------------------------------------------
// Misc utility functions
//---------------------------------------------------------------------------------------

// quick function to compute error in a ratio a/b
double frac_ratio_error(double a, double da, double b, double db)
{ return sqrt( da*da/a/a + db*db/b/b ); }


// throw error function for invalid option
void validate_option(TString option) {
	if ( option != "08/01"  && option != "08/02"  && option != "07/26a" && option != "07/26b"
		&& option != "04/27"  && option != "04/28"  && option != "08sum"  && option != "04sum"
		&& option != "09/12"  && option != "09/13"  && option != "09/13a" && option != "09/13b"
		&& option != "09/14"  && option != "09/14a" && option != "09/14b" && option != "09/15a"
		&& option != "09/15b" && option != "09/15c" && option != "09/15d") {
		cout << ">>> Invalid option " << option << " ! Aborting..." << endl;
		exit(1);
	}
}

// function to get iMin, iMax loop limits (INclusive) given an option string
void get_option_limits(TString option, size_t &iMin, size_t &iMax) {
	validate_option(option);
	if      (option == "07/26a") {iMin = 13;  iMax = 17;}
	else if (option == "07/26b") {iMin = 18;  iMax = 23;}
	else if (option == "08/01" ) {iMin = 28;  iMax = 41;}
	else if (option == "08/02" ) {iMin = 42;  iMax = 50;}
	else if (option == "04/27" ) {iMin = 52;  iMax = 55;}
	else if (option == "04/28" ) {iMin = 57;  iMax = 69;}
	else if (option == "08sum" ) {iMin = 28;  iMax = 50;}
	else if (option == "04sum" ) {iMin = 52;  iMax = 69;}
	else if (option == "09/12" ) {iMin = 76;  iMax = 80;}
	else if (option == "09/13" ) {iMin = 85;  iMax = 107;}
	else if (option == "09/13a") {iMin = 85;  iMax = 93;}
	else if (option == "09/13b") {iMin = 94;  iMax = 107;}
	else if (option == "09/14" ) {iMin = 115; iMax = 139;}
	else if (option == "09/14a") {iMin = 115; iMax = 125;}
	else if (option == "09/14b") {iMin = 127; iMax = 138;}
	else if (option == "09/15a") {iMin = 144; iMax = 153;} // 154 is significantly high in Det0
	else if (option == "09/15b") {iMin = 156; iMax = 160;} // this is the set of runs that had a bumped ionization chamber
	else if (option == "09/15c") {iMin = 161; iMax = 168;}
	else if (option == "09/15d") {iMin = 170; iMax = 183;}
}


bool is_problematic_NRF(size_t i) {
	// 23: lots of pileup
	// 42: gain issues
	// 43: gain issues
	// 44: gain issues
	// 54: unknown problem
	// 56: beam off
	// 58: beam blip
	if (i == 23 || i == 42 || i == 43 || i == 44 || i == 54 || i == 56 || i == 58) return true;
	else return false;
}

bool is_problematic_LYSO(size_t i) {
	// 23: current spikes
	// 31: no LYSO data
	// 43: no LYSO data
	// 53: no LYSO data
	// 56: beam off
	// 58: beam blip
	if (i == 23 || i == 31 || i == 43 || i == 53 || i == 56 || i == 58) return true;
	else return false;
}

bool is_problematic(size_t i) {
	if (is_problematic_NRF(i) || is_problematic_LYSO(i)) return true;
	else return false;
}

template <typename T>
void vector_print(vector<T> v)
{ for (size_t i = 0; i < v.size(); ++i) cout << v[i] << endl; cout << endl; }

void table_print(vector< vector<double> > v) {
	for (size_t j = 0; j < v[0].size(); ++j)
	{
		for (size_t i = 0; i < v.size(); ++i)
			cout << v[i][j] << "\t";
		cout << endl;
	}
}

template <typename T>
T vector_sum(vector<T> v) {
	T sum = 0.0;
	for (size_t i = 0; i < v.size(); ++i)
		sum += v[i];
	return sum;
}

template <typename T>
T vector_quadsum(vector<T> v) {
	T sum = 0.0;
	for (size_t i = 0; i < v.size(); ++i)
		sum += v[i]*v[i];
	return sqrt(sum);
}

// vector min/max, since iterators aren't working (namespace issue?)
// templates also not working
double vector_min(vector<double> v) {
	double m = v[0];
	for (size_t i = 0; i < v.size(); ++i)
		if (v[i] < m) m = v[i];
	return m;
}

double vector_max(vector<double> v) {
	double m = v[0];
	for (size_t i = 0; i < v.size(); ++i)
		if (v[i] > m) m = v[i];
	return m;
}

double vector_avg(vector<double> v) {
	double m = 0.0;
	for (size_t i = 0; i < v.size(); ++i)
		m += v[i];
	m /= (1.0*v.size());
	return m;
}

template <typename T>
vector<T> vector_scale(vector<T> v, double x) {
	vector<T> vnew;
	for (size_t i = 0; i < v.size(); ++i) vnew.push_back(v[i]*x);
	return vnew;
}

TGraphErrors *vector_err_plot(vector<double> x, vector<double> dx, vector<double> y, vector<double> dy) {
	int n = x.size();
	TGraphErrors *g = new TGraphErrors(n, &(x[0]), &(y[0]), &(dx[0]), &(dy[0]));
	return g;
}


TLatex *draw_label(TH1D *h, TString text = "label", double height = 0.8, double width = 0.6, Color_t color = NULL, bool ndc = true)
{
	TLatex *t = new TLatex(width, height, "#bf{"+text+"}");
	t->SetNDC(ndc);
	if (color == NULL) t->SetTextColor(h->GetLineColor());
	else t->SetTextColor(color);
	t->SetTextSize(h->GetXaxis()->GetLabelSize());
	t->Draw();
	return t;
}

// reimplementation for TGraph's rather than TH1D's; template function would be nicer but
// it gets persnickety about string <--> TString conversions
TLatex *draw_label(TGraph *g, TString text = "label", double height = 0.8, double width = 0.6, Color_t color = NULL, bool ndc = true)
{
	TLatex *t = new TLatex(width, height, "#bf{"+text+"}");
	t->SetNDC(ndc);
	if (color == NULL) t->SetTextColor(g->GetMarkerColor());
	else t->SetTextColor(color);
	t->SetTextSize(g->GetXaxis()->GetLabelSize());
	t->Draw();
	return t;
}

// draw arrows to show peaks and their branches
void overlay_arrows(TCanvas *c, TH1 *h1, TH1 *h2 = NULL)
{
	c->cd();
	vector<double> Evec;
	Evec.push_back(2.176); Evec.push_back(2.176-0.045);
	Evec.push_back(2.209); Evec.push_back(2.209-0.045);
	Evec.push_back(2.245); Evec.push_back(2.245-0.045);
	Evec.push_back(2.146); Evec.push_back(2.212);

	double x1_last = 0.0;
	double y1_last = 0.0;
	for (int i = 0; i < Evec.size(); i++)
	{
		// compute the arrow coordinates (x1,y1) -> (x2,y2)
		double x1 = h1->GetBinCenter(h1->FindBin(Evec[i]));
		double x2 = x1;

		double range = c->GetFrame()->GetY2() - c->GetFrame()->GetY1();

		double offset = 0.05*range;

		double y2_h1 = 0;
		double y2_h2 = 0;
		int x1_bin = h1->FindBin(x1);
		y2_h1 = h1->GetBinContent(x1_bin) + h1->GetBinError(x1_bin) + offset;
		if (h2 != NULL) y2_h2 = h2->GetBinContent(x1_bin) + h2->GetBinError(x1_bin) + offset;
		double y2 = max(y2_h1, y2_h2);
		double y1 = y2+1.5*offset;

		// draw the arrow
		TArrow *arrow = new TArrow(x1,y1,x2,y2,0.015,"|>");
		arrow->SetLineWidth(1);
		arrow->SetAngle(30);
		arrow->SetLineColor(kBlack);
		arrow->Draw();

		// draw lines between the branched decay arrows
		if (i%2 == 1 && i < 6){
			TLine *line = new TLine(x1_last, y1_last, x1, y1);
			line->SetLineWidth(1);
			line->SetLineColor(kBlack);
			line->Draw();
			double xlab = x1_last - 0.013;
			double ylab = y1_last + 0.5*offset;
			draw_label(hfinal0, "U-238", ylab, xlab, kBlack, false);
		}

		if (i == 6) draw_label(hfinal0, "U-238", y1 + offset/2.0, x1-0.0055, kBlack, false);
		if (i == 7) draw_label(hfinal0, "Al-27", y1 + offset/2.0, x1-0.0055, kBlack, false);

		x1_last = x1;
		y1_last = y1;
	}

}


// quick calibration from Th-232 gamma spectrum
TF1 *LaBr_Th232_cal()
{
	double energies[] = {0.238, 0.583, 0.911, 0.969, 2.614}; // MeV
	double adcs[]     = {690,   1660,  2600,  2750,  7200};
	TGraph *g = new TGraph(5, adcs, energies);
	g->SetMarkerSize(1);
	g->SetMarkerStyle(8);
	g->Draw("AP");
	TF1 *fit = new TF1("fit", "pol1");
	g->Fit("fit","M");
	return fit;
}


void draw_arrow_and_label(TCanvas *c, TH1 *h, double x, TString text)
{
	c->cd();
	c->Update();

	int xbin = h->FindBin(x);
	double range = c->GetFrame()->GetY2() - c->GetFrame()->GetY1();
	double offset = 0.05*range;
	double y2 = h->GetBinContent(xbin) + h->GetBinError(xbin) + offset;
	double y1 = y2 + 2*offset;
	TArrow *t = new TArrow(x,y1,x,y2,0.01,">");
	t->SetLineWidth(2);
	t->SetAngle(30);
	t->Draw();

	TLatex *label = new TLatex(x, y1 + offset/2.0, "#bf{"+text+"}");
	label->SetTextSize(0.75*h->GetXaxis()->GetLabelSize());
	label->SetTextAngle(90);
	label->SetTextAlign(12);
	label->Draw();
}


//---------------------------------------------------------------------------------------
// The main analysis code
//---------------------------------------------------------------------------------------

// Convention: 0 template, 1 candidate
// Useful pairs:
//  - 09/13a, 09/13b: genuine DU vs hoax Pb, Sept 13
//  - 09/14b, 09/14a: genuine DU vs hoax Pb, Sept 14 (reversed order)
//  - 09/15b, 09/15c: half-hoax with bumped ion chamber vs half-hoax
//  - 09/15a, 09/15c: double-thick DU target vs half-hoax, Sept 15
//  - 09/15a, 09/15d: double-thick DU target vs double-thick Pb target, Sept 15
//  - 09/15c, 09/15d: half-hoax vs double-thick Pb target, Sept 15
// Style codes:
//  - 0: raw histogram with arrow overlay, no fits (Areg's specs)
//  - 1: fits overlay, no arrows (Areg's specs)
//  - 2: data-rich style with text and fit overlay
TCanvas *comparisonPlot(TString cname, int detNumber, TString option0 = "09/13a", TString option1 = "09/13b", int styleCode = 2, bool rebin = true) {
	cout << "\n#### Beginning new analysis!\n####\n####" << endl;
	vector<int> detVec;
	if (detNumber >= 0) detVec.push_back(detNumber);
	else {
		cout << "Using combined datasets!" << endl;
		detVec.push_back(0); detVec.push_back(2); detVec.push_back(3);
	}

	if (detNumber < 0 && !rebin){
		cout << "Error! rebin = false only valid for one detector at a time! Exiting..." << endl;
		exit(1);
	}

	validate_option(option0);
	validate_option(option1);
	size_t iMin0, iMax0, iMin1, iMax1;
	get_option_limits(option0, iMin0, iMax0);
	get_option_limits(option1, iMin1, iMax1);

	// Output file to save rates, since this function needs to get called multiple times
	ofstream rfile;

	// First setup the histograms with arbitrary first Det, then overwrite them in the for loop.
	// This is just to establish the bins we will fill in the for loop + axis titles, etc. and
	// realistically could be done with just an empty histo created from a vector of bin low edges.
	vector< TH1D* > v = getCalHistVector(detVec[0]);
	VRebin(v,8);
	TString name0 = "hfinal0";
	TString name1 = "hfinal1";
	name0 += cname;
	name1 += cname;
	TH1D *hfinal0 = applyAutoCalibration(norm_sum(name0, v, detVec[0], iMin0, iMax0));
	TH1D *hfinal1 = applyAutoCalibration(norm_sum(name1, v, detVec[0], iMin1, iMax1));
	if (rebin){
		hfinal0 = rebin_interp(hfinal0);
		hfinal1 = rebin_interp(hfinal1);
	}

	hfinal0->Reset();
	hfinal0->Sumw2();
	hfinal1->Reset();
	hfinal1->Sumw2();

	for (size_t j = 0; j < detVec.size(); ++j) {
		vector< TH1D* > vp = getCalHistVector(detVec[j]);
		vector< double> t  = getLiveTimeVector(detVec[j]);
		vector< double> a  = current;
		VRebin(vp,8);

		// generate test rebinned (via VRebin above) histograms in order to get a good estimate of the autoCalibration's
		TH1D *htest0   = norm_sum("htest0", vp, detVec[j], iMin0, iMax0);
		TH1D *htest1   = norm_sum("htest1", vp, detVec[j], iMin1, iMax1);
		TF1  *autoCal0 = autoCalibration(htest0);
		TF1  *autoCal1 = autoCalibration(htest1);

		// Apply the autoCalibration to the roughly-calibrated histogram and rebin by 8x then by interpolations
		// Bin errors are recomputed for each detector, then propagated automatically in TH1::Add().
		vector< TH1D* > vpp = getCalHistVector(detVec[j]);
		VRebin(vpp,8);

		TH1D *hsum0 = calibrate(norm_sum(Form("hsum0_det%i",detVec[j]), vpp, detVec[j], iMin0, iMax0), autoCal0);
		if (rebin) hsum0 = rebin_interp(hsum0);
		double liveCharge0 = 0.0;
		for (size_t i = iMin0; i <= iMax0; ++i) liveCharge0 += t[i]*a[i];
		recompute_bin_errors(hsum0, liveCharge0);
		hfinal0->Add(hsum0);
	
		TH1D *hsum1 = calibrate(norm_sum(Form("hsum1_det%i",detVec[j]), vpp, detVec[j], iMin1, iMax1), autoCal1);
		if (rebin) hsum1 = rebin_interp(hsum1);
		double liveCharge1 = 0.0;
		for (size_t i = iMin1; i <= iMax1; ++i) liveCharge1 += t[i]*a[i];
		recompute_bin_errors(hsum1, liveCharge1);
		hfinal1->Add(hsum1);
	}

	hfinal0->SetTitle("hfinal0");
	hfinal0->SetName("hfinal0");
	hfinal1->SetTitle("hfinal1");
	hfinal1->SetName("hfinal1");
	hfinal0->SetStats(0);

	hfinal0->SetLineColor(kBlack);
	hfinal1->SetLineColor(kRed);

	double yrange = (detNumber == -1 ? 0.015 : 0.005);
	hfinal0->GetYaxis()->SetRangeUser(0,yrange);
	hfinal0->GetXaxis()->SetRangeUser(2.1,2.28);
	hfinal0->GetYaxis()->SetTitleOffset(1.3);
	hfinal0->GetYaxis()->SetTitle(Form("interp. counts per %2.2f keV per #muA#upoints (live)",hfinal0->GetBinWidth(1)*1.0e3));
	hfinal0->GetXaxis()->SetTitle("energy #it{E} [MeV]");

	TString detString = (detNumber >= 0 ? Form("%i", detNumber) : "ALL");

	TCanvas *c = new TCanvas(cname, cname, 900, 600);
	c->cd();

	hfinal0->Draw();
	hfinal1->Draw("same");

	gPad->SetLogy(0);
	gPad->SetTicks(1,1);

	if (styleCode > 0){
		int fitCode = -2;
	
		TF1 *fit0 = fit_uranium_multiplet(hfinal0, fitCode);
		vector<double> areas0, dareas0, energies0;
		areas0 = gaussian_areas(fit0, hfinal0, dareas0, energies0, true);
	
		TF1 *fit1 = fit_uranium_multiplet(hfinal1, fitCode);
		vector<double> areas1, dareas1, energies1;
		areas1 = gaussian_areas(fit1, hfinal1, dareas1, energies1, true);
	
		// get the Al-27 areas
		double  area0_2212 =  areas0[areas0.size()-1];
		double darea0_2212 = dareas0[areas0.size()-1];
		double  area1_2212 =  areas1[areas1.size()-1];
		double darea1_2212 = dareas1[areas1.size()-1];
	
		if (fitCode == -2){
			// remove the 2.212 MeV Al-27 contribution and 2.209 MeV U-238 contribution, possibly others
			// sigmas should be ~8,8,5,10 with peaksToCut=2, and similar otherwise
			int peaksToCut = 2;
			for (size_t i = 0; i < peaksToCut; ++i)
			{
				areas0.pop_back();
				areas1.pop_back();
				dareas0.pop_back();
				dareas1.pop_back();
			}
		}
		for (size_t i = 0; i < areas0.size(); ++i)
			cout << "  using E = " << energies0[i] << " MeV for fit statistic" << endl;
		cout << endl;

		cout << Form("2.176 MeV (genu) area: (%1.2f ± %1.2f) / mC", 1e3*areas0[0], 1e3*dareas0[0]) << endl;
		cout << Form("2.176 MeV (hoax) area: (%1.2f ± %1.2f) / mC", 1e3*areas1[0], 1e3*dareas1[0]) << endl;
		cout << Form("2.245 MeV (genu) area: (%1.2f ± %1.2f) / mC", 1e3*areas0[2], 1e3*dareas0[2]) << endl;
		cout << Form("2.245 MeV (hoax) area: (%1.2f ± %1.2f) / mC", 1e3*areas1[2], 1e3*dareas1[2]) << endl;
		cout << Form("2.212 MeV (genu) area: (%1.2f ± %1.2f) / mC", 1e3*area0_2212, 1e3*darea0_2212) << endl;
		cout << Form("2.212 MeV (hoax) area: (%1.2f ± %1.2f) / mC", 1e3*area1_2212, 1e3*darea1_2212) << endl;
		cout << endl;

		TString rebin_str = (rebin ? "interp" : "raw");
		rfile.open(Form("rates_file_%s_%s_det%s.txt", cname.Data(), rebin_str.Data(), detString.Data()));
		rfile << Form("%.6f, %.6f, 2.176,%s", areas0[0], dareas0[0], option0.Data()) << endl;
		rfile << Form("%.6f, %.6f, 2.176,%s", areas1[0], dareas1[0], option1.Data()) << endl;
		rfile << Form("%.6f, %.6f, 2.245,%s", areas0[2], dareas0[2], option0.Data()) << endl;
		rfile << Form("%.6f, %.6f, 2.245,%s", areas1[2], dareas1[2], option1.Data()) << endl;
		rfile << Form("%.6f, %.6f, 2.212,%s", area0_2212, darea0_2212, option0.Data()) << endl;
		rfile << Form("%.6f, %.6f, 2.212,%s", area1_2212, darea1_2212, option1.Data()) << endl;
		rfile.close();
	
		// for later...
		double rc0 = fit0->GetChisquare()/(1.0*fit0->GetNDF());
		double rc1 = fit1->GetChisquare()/(1.0*fit1->GetNDF());
	
		// print the individual peak discrepancies and their quadrature sum
		vector<double> indivDiscs;
		for (size_t i = 0; i < areas0.size(); ++i){
			double 	a0 = areas0[i];
			double 	a1 = areas1[i];
			double da0 = dareas0[i];
			double da1 = dareas1[i];
			double disc = (a1-a0)/sqrt(da1*da1 + da0*da0);
			indivDiscs.push_back(disc);
			cout << energies0[i] << " MeV individual discrepancy: " << disc << " sigma." << endl;
		}
		cout << "quadrature sum of sigmas: " << vector_quadsum(indivDiscs) << endl;
	
		double  totalArea0 = 3600.0 * vector_sum(areas0);
		double  totalArea1 = 3600.0 * vector_sum(areas1);
		double dtotalArea0 = 3600.0 * vector_quadsum(dareas0);
		double dtotalArea1 = 3600.0 * vector_quadsum(dareas1);
	
		cout << option0 << " vs " << option1 << ", Det" << detString << endl;
		cout << "total area (genu) = (" << totalArea0 << " ± " << dtotalArea0 << ") / uA.h (live)" << endl;
		cout << "total area (hoax) = (" << totalArea1 << " ± " << dtotalArea1 << ") / uA.h (live)" << endl;
	
		double ratio = totalArea1/totalArea0;
		double dratio = ratio * frac_ratio_error(totalArea0, dtotalArea0, totalArea1, dtotalArea1);
		double discrepancy = (totalArea1 - totalArea0)/sqrt(dtotalArea1*dtotalArea1 + dtotalArea0*dtotalArea0);
		cout << "ratio = (" << ratio << " ± " << dratio << ")" << endl;
		cout << "discrepancy = " << discrepancy << " sigma" << endl;
	
		if (styleCode == 1){
			TLegend *legend = new TLegend(0.1,0.7,0.45,0.9); // x1, y1, x2, y2; for top-left legend, only change y1, x2
   		legend->AddEntry(hfinal0,"template (DU) spectrum","lpe");
   		TString cand_spec_text = (cname == "c1" ? "candidate (DU) spectrum" : "hoax (Pb) spectrum");
   		legend->AddEntry(hfinal1,cand_spec_text,"lpe");
   		legend->AddEntry(hfinal0->GetFunction("multiplet_fit_hfinal0"),"template (DU) fit","l");
   		TString cand_fit_text = (cname == "c1" ? "candidate (DU) fit" : "hoax (Pb) fit");
   		legend->AddEntry(hfinal1->GetFunction("multiplet_fit_hfinal1"),cand_fit_text,"l");
   		legend->AddEntry((TObject*)0, Form("discrepancy = %2.2f #sigma", discrepancy),"");
   		legend->Draw();
			hfinal0->SetTitle("");
		}
		else if (styleCode == 2){
			draw_label(hfinal0, Form("net fit area = (%2.2f #pm %2.2f)/#muA#upointh (#chi^{2}/#nu = %2.2f)", totalArea0, dtotalArea0, rc0), 0.85, 0.135);
			draw_label(hfinal1, Form("net fit area = (%2.2f #pm %2.2f)/#muA#upointh (#chi^{2}/#nu = %2.2f)", totalArea1, dtotalArea1, rc1), 0.81, 0.135);
			draw_label(hfinal0, Form("#rightarrow ratio = %2.2f #pm %2.2f", ratio, dratio), 0.77, 0.135);
			draw_label(hfinal0, Form("#rightarrow discrepancy = %2.2f #sigma", discrepancy), 0.73, 0.135);
			hfinal0->SetTitle(Form("HVRL Det%s genuine (black, %s) vs hoax (red, %s)",detString.Data(),option0.Data(),option1.Data()));
		}

		cout << "Al-27 areas:" << endl;
		cout << "  genu: " << 3600.0*area0_2212 << " ± " << 3600.0*darea0_2212 << " per uA.h" << endl;
		cout << "  hoax: " << 3600.0*area1_2212 << " ± " << 3600.0*darea1_2212 << " per uA.h" << endl;

		cout << "1-2 MeV integral difference = "
		     << hfinal1->Integral(hfinal1->FindBin(1.0),hfinal1->FindBin(2.0))/hfinal0->Integral(hfinal0->FindBin(1.0),hfinal0->FindBin(2.0)) << endl;
	}
	else if (styleCode == 0){
		c->Update();
		hfinal0->SetTitle("");
		overlay_arrows(c, hfinal0, hfinal1);
	}

	c->Update();
	return c;
}


void plotAll(int styleCode = 1) {
	int detCode = -1;
	bool rebin = true;
	TCanvas *c0 = comparisonPlot("c0", detCode, "09/13a", "09/13b", styleCode, rebin); c0->SaveAs("c0.eps"); // template I vs hoax Ia
	TCanvas *c1 = comparisonPlot("c1", detCode, "09/13a", "09/14b", styleCode, rebin); c1->SaveAs("c1.eps"); // template I vs genuine candidate Ig
	TCanvas *c2 = comparisonPlot("c2", detCode, "09/13a", "09/14a", styleCode, rebin); c2->SaveAs("c2.eps"); // template I vs hoax Ib
	TCanvas *c3 = comparisonPlot("c3", detCode, "09/15a", "09/15d", styleCode, rebin); c3->SaveAs("c3.eps"); // template II vs hoax IIc
	TCanvas *c4 = comparisonPlot("c4", detCode, "09/15a", "09/15c", styleCode, rebin); c4->SaveAs("c4.eps"); // template II vs hoax IId
}


void recompute_bin_errors(TH1 *h, double scaleFactor){
	for (size_t i = 1; i <= h->GetNbinsX(); ++i){
		h->SetBinError(i, TMath::Sqrt(h->GetBinContent(i)*scaleFactor)/scaleFactor);
	}
}


// Generic rebinning of h0 to a new histogram with bin_centers via linear interpolation.
// Sets new bin error to sqrt of bin content, so make sure to call before scaling by live charge.
// Defaults to [0.0, 3.0) MeV in 1 keV bins.
TH1D *rebin_interp(TH1D *h0, vector<double> bin_low_edges = NULL){
	if (bin_low_edges == NULL){
		vector<double> be;
		double emin = 0.0;
		double emax = 3.0;
		double ebin = 0.001;
		for (double e = emin; e < emax; e += ebin){
			be.push_back(e);
		}
		bin_low_edges = be;
	}

	TString name = h0->GetName();
	name += "_rb";
	TH1D *hr = new TH1D(name, name, bin_low_edges.size()-1, &(bin_low_edges[0]));
	hr->Sumw2();

	for (size_t i = 1; i <= hr->GetNbinsX(); ++i) {
		double x = hr->GetBinCenter(i);
		double y = h0->Interpolate(x);
		hr->SetBinContent(i, y);
		hr->SetBinError(i, 0); // set to 0 here, then call recompute_bin_errors
	}

	// Subtle point: make sure to account for the change in bin width.
	//
	// E.g. initial two bins at 500 and 1000 counts with bin widths of 1 keV. Integral is simply 1500 counts.
	// Now interpolate a new bin between them: 750 counts in the bin. Unscaled bins are now 500,750,1000, and
	// it looks like we've created counts out of nowhere because the integral is 2250. But our bin width is now
	// only two thirds what it was before (same histogram range is covered by three bins instead of two), so
	// multiply through: 2250*2/3 = 1500, and our new bin contents are 333,500,667. While the "counts per bin"
	// values have changed, if we're careful about how we tally counts and instead compute "counts per x keV",
	// where x is the binwidth, we recover the 500,750,1000 as estimates of the 'differential counts'.
	//
	// The remaining question is to what do we set the uncertainty of the new bins? Probably just sqrt(new counts), right?
	double scale = hr->GetBinWidth(1)/h0->GetBinWidth(1);
	//cout << "SCALE = " << scale << endl;
	hr->Scale(scale);

	return hr;
}


