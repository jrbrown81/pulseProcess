/*
 pulseProcess

Version 5.0
 Looks for header data in file to define record length.
 If not found, issues warning and reads from settings.txt
 
 Version 4.0
 Now reads in pulse shaping settings at run time from a file ('settings.txt').

 Version 3.0
 Suspect there are some other 2.x versions around so jump to version 3.0
 Corrected bug which wasn't using baseline subtracted data for gaussian shaping.
 Now stores pulses by default. -p options doesn't store pulses.
 baseline subtracted data now written to file as well.

 Version 2.1
 Messed about with version, e.g. defaults to saving pulses again, can be used with -n option...

 Version 2.0
 Made further additions to gaussianfilter. Can now define number of poles,
 specify the pole zero correction and gain factor scales shaped pulse to comparible
 amplitude to un shaped pulse.

 Version 1.1
 Doesn't save pulses by default, can save them using -p option.
 Can also  only process some hits with -n option. Note, can't use multiple options.

 e.g. ./pulseProcess -n 10 wave0.txt wave0.root

 to only process first 10 pulses.

 Version 1.0
 Some old code removed, and gaussianfilter parameters are now defined near the top
 for ease of adjustment.

 Based on txt2tree.C (see below).

-------------------
  Developed from:

  txt2treeJB.C

  Originally created by Paul Davies on 14/09/2015.
Modified by Stefanos Paschalis 02/2017 for MPhys projects (Harry PSA with V1730Band Ewan - Muons with V1740)
 Nuclear Physics Group, Department of Physics, University of York

 Further modified by J. Brown for postion senistivity analysis, 23/05/17
-------------------
*/
/*

To compile this code needs to be linked to ROOT libraries, use:

g++ `root-config --cflags --libs` -Wl,--no-as-needed -lTree -lRIO -lCore -lHist pulseProcess.C -o pulseProcess5

--no-as-needed, -lTree -lRIO -lCore -lHist is used to work around some linking issues in Ubuntu. Probably not needed on less pedantic systems

data is saved into a tree, 1 histogram per pulse. To look at them in root use:

tree->Draw("hist.Draw()","","goff",1,12) which will plot one event starting at event 1

*/


#include "TMath.h"	
#include "TF1.h"
#include "TH1.h"
#include "TFitResult.h"
#include "TMatrixDSym.h"
#include "TRandom.h"
#include "TStyle.h"
#include "TTree.h"
#include "TFile.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TGraph.h"

#include <fstream>

//switch on debugging print statements
// 0 nothing
// 1 header data being read in
// 2 header data being written to file
// 4 both of the above
// 8 all data - pulse and header
#define debug 0

//using namespace std;

std::string version="5.0";
// New options to store pulses (or not) and process a limited number of pulses

struct headerdata {
//for info Wavedump txt files have the following header info

/*
Record Length: 8200
BoardID:  0
Channel: 0
Event Number: 4122
Pattern: 0x0000
Trigger Time Stamp: 2109637025
DC offset (DAC): 0xF5C1
*/

  UInt_t recordlengthI;
  double recordlength;
  double boardid;
  double channel;
  double eventnumber;
  double pattern;
  double triggertimestamp;
  double dcoffset;

} ;


// Some user parameters.
const int array_length=70000;
// These are overwritten by values read in from 'settings.txt'
// range to be used for baseline integration
double bslMin=0;
double bslMax=1000;
int polarity= 1;  // 1 for positive pulses, -1 for negative
int timeConst=1000;  // samples
double poleZ=0.0015;
int n_poles=4;
int headers=0;
int traceLength=8200;

//	Gain = (n_poles!)/(n_poles^n_poles * exp(-n_poles));
	double gain=5.11858; // for n_poles=4
//	double gain=3.69453; // for n_poles=2
//	double gain=6.22575; // for n_poles=6

struct pulsedata {
  double data[array_length];  // pulse
  double data_bsl[array_length];  // baseline subtracted pulse
  double data_cfd[array_length];
  double data_t0_aligned[array_length];
  double data_gaussianfilter[array_length];
  double bsl;
//   double pulseMax;
  double Q1;
  double Q2;
  double Q1pp;
  double Q2pp;
  double t0;
  double gaussmax;

} ;

// int timeStep[array_length];

//declaring functions which are defined below
//int txt2tree(char * infile, char * outfile);
int txt2tree(char * infile, char * outfile, bool storePulses, int nPulses=-100);
void Usage();
void readSettings(TString settingsFile);
void printheaderdata(headerdata hdrdata);
void clearheaderdata(headerdata hdrdata);
void clearpulsedata(pulsedata &plsdata);
void postprocess(pulsedata &plsdata);
void baseline(pulsedata &plsdata);
void CFD(pulsedata &plsdata, Int_t samples, Int_t delay, Float_t fraction);
void GetCFD(pulsedata &plsdata, Int_t samples);
void t0_align(pulsedata &plsdata, Int_t samples);
void Q1Q2pp(pulsedata &plsdata, Int_t samples);
void gaussianfilter(pulsedata &plsdata, Int_t samples, Double_t tau, Double_t pz, int n_poles);
void GetGaussPeak(pulsedata &plsdata, Int_t samples);


int main(int argc, char **argv)
{

   std::cout << std::endl << "------------------------" << std::endl
   << "pulseProcess Version " << version << std::endl
   << "------------------------" << std::endl << std::endl;

if(argc == 3){

//  std::cout << argv[0] << std::endl;
  std::cout << "input filename:  " << argv[1] << std::endl;
  std::cout << "output filename: " << argv[2] << std::endl;
	
	txt2tree(argv[argc-2], argv[argc-1],1);
}
else if(argc > 3){
	if(argv[1][0]=='-') {
		if(argv[1][1]=='p') {
  		   std::cout << "input filename:  " << argv[argc-2] << std::endl;
  		   std::cout << "output filename: " << argv[argc-1] << std::endl;
			std::cout << std::endl << "******************************" << std::endl
         << "Not storing pulses to rootfile" << std::endl
         << "******************************" << std::endl << std::endl;
			txt2tree(argv[argc-2], argv[argc-1],0);
		}
		else if(argv[1][1]=='n' && argc==5 && isdigit(argv[2][0])) {
         int nPulses=atoi(argv[2]);
  			std::cout << "input filename:  " << argv[argc-2] << std::endl;
  			std::cout << "output filename: " << argv[argc-1] << std::endl;
			std::cout << "Processing " << nPulses << " pulses" 	<< std::endl;
			txt2tree(argv[argc-2], argv[argc-1], 1, nPulses);
		}
		else {
			Usage();
			return 0;
		}
	}
	else {
		Usage();
		return 0;
	}
}
else{
  Usage();
  return 0;
}

return 1;

}

int txt2tree(char * infile, char * outfile, bool storePulses, int nPulses){

//   settings
   std::ifstream settings("settings.txt");
   if (!settings){
      std::cout << "'settings.txt' file not found! Exiting." << std::endl;
      return 0;
   } else readSettings("settings.txt");


std::ifstream in(infile);

std::string data;

size_t first;
size_t last;

headerdata hdrdata;
pulsedata plsdata;

//some counters to help with sanity checks
//number of points read in from a pulse
int datapointcounter=0;
//number of waveforms found in file
int waveformcounter = 0;
int processedCounter = 0;
// for(int i=0; i<array_length; i++) timeStep[i]=i;

//the output file and tree
TFile * out = new TFile(outfile,"RECREATE");
TTree * tree = new TTree("tree","tree");

//create a branch to store the histogram pulses
// tree->Branch("hist","TH1I",&hist,32000,0);

//for info this is where we will store the header info, this is declared at the start of the file
/*
struct headerdata {
  double recordlength;
  double boardid;
  double channel;
  double eventnumber;
  double pattern;
  double triggertimestamp;
  double dcoffset;

} ;
*/
//create a branch for the header information
tree->Branch("headerinfo",&hdrdata,"recordlength/I:recordlength/d:boardid/d:channel/d:eventnumber/d:pattern/d:triggertimestamp/d:dcoffset/d");

tree->Branch("recordlength",&hdrdata.recordlengthI,"recordlengthI/I");

if(storePulses) {
   std::cout << "******* Storing pulses to root tree *******" << std::endl;
   
   if(!headers) hdrdata.recordlengthI=traceLength;
   
	tree->Branch("pulsedata",plsdata.data,"pulsedata[recordlengthI]/d");
// tree->Branch("pulsedata",plsdata.data,"pulsedata[10000]/d");
 	tree->Branch("pulsedata_bsl",plsdata.data_bsl,"pulsedata_bsl[recordlengthI]/d");
//  tree->Branch("timeStep",&timeStep,"timeStep[10000]/I");
//  tree->Branch("pulsedata_cfd",plsdata.data_cfd,"pulsedata_cfd[10000]/d");
//  tree->Branch("pulsedata_t0_aligned",plsdata.data_t0_aligned,"pulsedata_t0_aligned[10000]/d");
	tree->Branch("pulsedata_gauss",plsdata.data_gaussianfilter,"pulsedata_gauss[recordlengthI]/d");
}

//tree->Branch("bsl",&plsdata.bsl,"bsl/d");
// tree->Branch("pulseMax",&plsdata.pulseMax,"pulseMax/d");
tree->Branch("Q1",&plsdata.Q1,"Q1/d");
tree->Branch("Q2",&plsdata.Q2,"Q2/d");
tree->Branch("Q1pp",&plsdata.Q1pp,"Q1/d");
tree->Branch("Q2pp",&plsdata.Q2pp,"Q2/d");
tree->Branch("t0",&plsdata.t0,"t0/d");
tree->Branch("gaussmax",&plsdata.gaussmax,"gaussmax/d");

//lets start with a clean header data structure
if(headers) clearheaderdata(hdrdata);
clearpulsedata(plsdata);

while( !in.eof() && waveformcounter!=nPulses ){

	//get the data and save as a string
	getline(in,data);
	//reading in as a string so lets use the first item in the header to denote a new event
	//then check for each of the header titles, if it doesn't contain a header title it must be data
	if( data.find("Record Length: ") != std::string::npos ) {
      if(datapointcounter==0 && headers==0) {
         std::cout << "******* Headers detected *******" << std::endl;
         headers=1;
      }
		if(debug) std::cout << "*******Beginning of new event header********" << std::endl;
		if(debug == 1 || debug == 4 || debug == 8) std::cout << data << "  " << std::endl;
		//some string maths which will delete the title and leave just the numerical value
		last = data.find(":") + 2; // don't want the : or the whitespace 
		data.erase(0,last);
	   if(headers) {
         hdrdata.recordlength = atof(data.data());
	      hdrdata.recordlengthI = atoi(data.data());
         traceLength=hdrdata.recordlengthI;
         if(datapointcounter==0) std::cout << "Using Trace Length: " << traceLength << std::endl;
      }
	}
	else if( data.find("BoardID:") != std::string::npos ) {
		if(debug == 1 || debug == 4 || debug == 8) std::cout << data << std::endl;
		last = data.find(":") + 2; // don't want the : or the whitespace 
		data.erase(0,last);
	   if(headers) hdrdata.boardid = atof(data.data());
	}
	else if( data.find("Channel:") != std::string::npos ) {
		if(debug == 1 || debug == 4 || debug == 8) std::cout << data << std::endl;
		last = data.find(":") + 2; // don't want the : or the whitespace 
		data.erase(0,last);
	   if(headers) hdrdata.channel = atof(data.data());
	}
	else if( data.find("Event Number:") != std::string::npos ) {
		if(debug == 1 || debug == 4 || debug == 8) std::cout << data << std::endl;
		last = data.find(":") + 2; // don't want the : or the whitespace 
		data.erase(0,last);
	   if(headers) hdrdata.eventnumber = atof(data.data());
	}
	else if( data.find("Pattern:") != std::string::npos ) {
		if(debug == 1 || debug == 4 || debug == 8) std::cout << data << std::endl;
		last = data.find(":") + 2; // don't want the : or the whitespace 
		data.erase(0,last);
	   if(headers) hdrdata.pattern = atof(data.data());
	}
	else if( data.find("Trigger Time Stamp:") != std::string::npos ) {
		if(debug == 1 || debug == 4 || debug == 8) std::cout << data << std::endl;
		last = data.find(":") + 2; // don't want the : or the whitespace 
		data.erase(0,last);
	   if(headers) hdrdata.triggertimestamp = atof(data.data());
	}
	else if( data.find("DC offset (DAC):") != std::string::npos ) {
		if(debug == 1 || debug == 4 || debug == 8) std::cout << data << std::endl;
		last = data.find(":") + 2; // don't want the : or the whitespace 
		data.erase(0,last);
	   if(headers) {
         hdrdata.dcoffset = atof(data.data());
	      if(debug == 2 || debug == 4 || debug == 8) printheaderdata(hdrdata);
      }
	}
	else { // This must be data as not header info.
      if(headers==0 && waveformcounter==0 && datapointcounter==0) {
         std::cout << std::endl
                     << "|*******************************************|" << std::endl
                     << "|  Warning: No headers found in this file   |" << std::endl
                     << "|  Is traceLength in settings file correct? |" << std::endl
                     << "|*******************************************|" << std::endl << std::endl;
      }
 
      // Check if complete pulse has been read, and process
      if(datapointcounter==traceLength) {
         waveformcounter++;
         if(debug) std::cout << "*******Processing event********" << std::endl;
         postprocess(plsdata);
         processedCounter++;
         if(debug) std::cout << "*******Filling tree********" << std::endl;
         tree->Fill();
         if (processedCounter%10000 == 1) {
            if(debug) std::cout << "*******Writing tree********" << std::endl;
            tree->Write("",TObject::kOverwrite);
         }
      // now we collect data from the next event so we should reset the datapointcounter
         if(headers) clearheaderdata(hdrdata);
         clearpulsedata(plsdata);
         datapointcounter = 0;
      //print something so we know the whole thing is working
         if(processedCounter%1000 == 0) std::cout << "Pulses processed: " << processedCounter << std::endl;
      }
		datapointcounter++;
		double tmpdata = atof(data.data());
		plsdata.data[datapointcounter-1] = atof(data.data());
      if(debug == 8) std::cout << waveformcounter << " " << datapointcounter << " " << plsdata.data[datapointcounter-1] << std::endl;

// note, need the -1 here as datapointcounter starts from 1
	   if(datapointcounter-1>=bslMin && datapointcounter-1<bslMax){
		   plsdata.bsl+=atof(data.data());}
	}
   
} // End of "while( !in.eof() && waveformcounter!=nPulses )"

std::cout << std::endl << "Total pulses processed: " << processedCounter << std::endl;
//////////

out->cd();
tree->Write("",TObject::kOverwrite);
out->Close();

return 1;

}

//some helper functions
void Usage(){

std::cout << "Usage: \n" << "'./pulseProcess <inputfile> <outputfile>'" << std::endl;
std::cout << "'./pulseProcess [-p] <inputfile> <outputfile>'   to NOT store pulses" << std::endl;
std::cout << "'./pulseProcess [-n <N_pulses>] <inputfile> <outputfile>'   to process a limited number of pulses" << std::endl;

}

void printheaderdata(headerdata hdrdata){

std::cout << "Record Length: " << hdrdata.recordlength <<std::endl;
std::cout << "Board ID: " << hdrdata.boardid << std::endl;
std::cout << "Channel: " << hdrdata.channel << std::endl;
std::cout << "Event Number: " << hdrdata.eventnumber << std::endl;
std::cout << "Pattern: " << hdrdata.pattern << std::endl;
std::cout << "Trigger Time Stamp: " << hdrdata.triggertimestamp << std::endl;
std::cout << "DC Offset (DAQ): " << hdrdata.dcoffset << std::endl;

}

void clearheaderdata(headerdata hdrdata){

hdrdata.recordlength = 0.0;
hdrdata.boardid = 0.0;
hdrdata.channel = 0.0;
hdrdata.eventnumber = 0.0;
hdrdata.pattern = 0.0;
hdrdata.triggertimestamp = 0.0;
hdrdata.dcoffset = 0.0;

}


void clearpulsedata(pulsedata &plsdata){

	memset(plsdata.data,0,sizeof(plsdata.data));
	memset(plsdata.data_bsl,0,sizeof(plsdata.data_bsl));
	memset(plsdata.data_cfd,0,sizeof(plsdata.data_cfd));
	memset(plsdata.data_t0_aligned,0,sizeof(plsdata.data_t0_aligned));
	memset(plsdata.data_gaussianfilter,0,sizeof(plsdata.data_gaussianfilter));
	plsdata.Q1 = -1;
	plsdata.Q2 = -1;
	plsdata.Q1pp = -1;
	plsdata.Q2pp = -1;
	plsdata.bsl = 0;
// 	plsdata.pulseMax = -1;
	plsdata.gaussmax = -1;
}

void postprocess(pulsedata &plsdata){

  baseline(plsdata);
//  CFD(plsdata,array_length,80,0.30);
//  GetCFD(plsdata,array_length);
// t0_align(plsdata,array_length);
//  Q1Q2pp(plsdata,array_length);
//gaussianfilter(plsdata,array_length,10,10);
//	if(poleZ!=0) poleZ=1/poleZ;
//	if(poleZ==0) poleZ=1/poleZ;
  gaussianfilter(plsdata,array_length,timeConst,poleZ,n_poles);
  GetGaussPeak(plsdata, array_length);
}

void baseline(pulsedata &plsdata){

  for(Int_t i=0 ; i<array_length; i++){
//    plsdata.data_bsl[i] = plsdata.data[i]-plsdata.bsl/50.;
      if(plsdata.data[i]!=0) plsdata.data_bsl[i] = polarity*(plsdata.data[i]-plsdata.bsl/(bslMax-bslMin));
  }

}


void CFD(pulsedata &plsdata, Int_t samples, Int_t delay, Float_t fraction){

  for(Int_t i=0 ; i<samples; i++){
    if(i+delay<samples)
      plsdata.data_cfd[i] = plsdata.data_bsl[i]-fraction*plsdata.data_bsl[i+delay];
    else
      plsdata.data_cfd[i] = plsdata.data_bsl[i]-fraction*plsdata.data_bsl[i];
  }
}


void GetCFD(pulsedata &plsdata, Int_t samples){

  TGraph gr;
  TH1F hi("t0","",array_length,10,200);
  gr.SetName("gr");
  //  TF1 *fexp = new TF1("fexp","expo");
  TF1 *f1 = new TF1("f1","pol1");

  Int_t nmax, nmin, n;
  Double_t max, min, zero_crossing, slope, t0;
  nmax = 0; max = 0.;
  nmin = 0; min = 0.;
  n=0;
  zero_crossing = 0.; slope = 0.; t0 = 0.;

  //  for(Int_t i=0; i<samples; i++){
  for(Int_t i=400; i<500; i++){
    //    if(i<600) zero_crossing+=plsdata.data_cfd[i];
    if(plsdata.data_cfd[i] > max) {nmax = i; max = plsdata.data_cfd[i];}
    if(plsdata.data_cfd[i] < min) {nmin = i; min = plsdata.data_cfd[i];}
  }
  zero_crossing = zero_crossing/600.;
  //  std::cout << nmin << "\t" << min << "\t" << nmax << "\t" << max << std::endl;
  for(Int_t i=nmin; i<nmax; i++){
    gr.SetPoint(n,i,plsdata.data_cfd[i]);
    n++;
    //    if(plsdata.data_cfd[i]<=zero_crossing && plsdata.data_cfd[i+1]>=zero_crossing){
    if(plsdata.data_cfd[i]<=0 && plsdata.data_cfd[i+1]>=0){
      t0 = (2*i+1)/2.;
      hi.Fill(i);
    }
  }
  t0 = hi.GetMean();
  //  std::cout << "t0 from hi\t" << t0 << std::endl;
  plsdata.t0 = t0;
  if(n>1){
    //   gr.Fit("f1","Q");
    //    std::cout << "t0\t" << gr.GetFunction("f1")->GetX(zero_crossing,nmin,nmax) << std::endl;
  }
  gr.Clear();
  delete f1;
  hi.Clear();
}


void t0_align(pulsedata &plsdata, Int_t samples){

  for(Int_t i=0; i<samples; i++){
    //    if(i+(1500-plsdata.t0))
    if(i-(80-plsdata.t0)<samples)
      plsdata.data_t0_aligned[i] = plsdata.data_bsl[i-(80-Int_t(plsdata.t0))];
    else
      plsdata.data_t0_aligned[i] = -100;

      }
}

void Q1Q2pp(pulsedata &plsdata, Int_t samples){
  for(Int_t i=0; i<samples; i++){
    if(i>79 && i<=125){
      plsdata.Q1pp+=plsdata.data_t0_aligned[i];}
    else if(i>500 && i<=1700){
      plsdata.Q2pp+=plsdata.data_t0_aligned[i];}
  }
}

void gaussianfilter(pulsedata &plsdata, Int_t samples, Double_t tau, Double_t pz, int n_poles) {

	int i, j;
	double a0, a1, b1;
	double *y, *z;
	double x1, x2;

	y = new double[samples];
	z = new double[samples];

	for(i = 0; i < samples; i++){
		y[i] = 0.0;
		z[i] = 0.0;
	}

//single pole high pass filter
//with pole-zero correction
	b1 = TMath::Exp(-1.0 / tau);
	a0 = (1.0 + b1) / 2.0;
	a1 = -1*(1.0+b1) / 2.0;

	for(i = 1; i < samples; i++){
		x2 = plsdata.data_bsl[i];
		x1 = plsdata.data_bsl[i - 1];
//		y[i] = b1*y[i-1] + a0 * x2 + a1*x1+ x1 / pz;
		y[i] = b1*y[i-1] + a0 * x2 + a1*x1 + x1*pz;
	}
     
// n-pole low pass filter
	b1 = TMath::Exp(-1.0 / tau);
	a0 = 1.0 - b1;
	
	for(j = 0; j < n_poles; j++){
		for(i = 1; i < samples; i++)
			z[i] = b1 * z[i-1] + a0*y[i];
		for(i = 1; i < samples; i++)
			y[i] = z[i];
	}
	for(i = 1; i < samples; i++) plsdata.data_gaussianfilter[i]=gain*z[i];

//// average of every 10 samples, to shorten shaped pulse length
//	double sum=0;
//	int counter=1;
//	int index=1;
//	for(i = 1; i < samples; i++) {
//		sum+=z[i];
//		counter++;
//		if(counter==10) {
//			plsdata.data_gaussianfilter[index] = gain*sum/counter;
//			index++;
//			sum=0;
//			counter=1;
//		}
//
//// 		plsdata.data_gaussianfilter[i] = /*gain**/ z[i];
//	}

	delete[] y;
	delete[] z;

}


void GetGaussPeak(pulsedata &plsdata, Int_t samples){

 Double_t max=0;

 for(Int_t i = 1; i < samples; i++)
   if(max < plsdata.data_gaussianfilter[i])
     max = plsdata.data_gaussianfilter[i];

 plsdata.gaussmax = max;
}

void readSettings(TString settingsFile){

   std::ifstream settings(settingsFile);
   std::cout << std::endl << "Reading pulse shaping settings from file: '" << settingsFile << "'" << std::endl;
   std::string line;
   while(!settings.eof() ){
      getline(settings,line);
      TString sline = TString(line);
      if (sline.BeginsWith("/")) continue;
      if (sline.Contains("bslMin: ")) {
         sline.ReplaceAll("bslMin: ","");
         bslMin = sline.Atoi();
         std::cout << "bslMin = " << bslMin << std::endl;
      }
      if (sline.Contains("bslMax: ")) {
         sline.ReplaceAll("bslMax: ","");
         bslMax = sline.Atoi();
         std::cout << "bslMax = " << bslMax << std::endl;
      }
      if (sline.Contains("polarity: ")) {
         sline.ReplaceAll("polarity: ","");
         polarity = sline.Atoi();
         std::cout << "polarity = " << polarity << std::endl;
      }
      if (sline.Contains("timeConst: ")) {
         sline.ReplaceAll("timeConst: ","");
         timeConst = sline.Atoi();
         std::cout << "timeConst = " << timeConst << std::endl;
      }
      if (sline.Contains("poleZ: ")) {
         sline.ReplaceAll("poleZ: ","");
         poleZ = sline.Atof();
         std::cout << "poleZ = " << poleZ << std::endl;
      }
      if (sline.Contains("n_poles: ")) {
         sline.ReplaceAll("n_poles: ","");
         n_poles = sline.Atoi();
         std::cout << "n_poles = " << n_poles << std::endl;
      }
//      if (sline.Contains("headers: ")) {
//         sline.ReplaceAll("headers: ","");
//         headers = sline.Atoi();
//         std::cout << "headers = " << headers << std::endl;
//      }
      if (sline.Contains("traceLength: ")) {
         sline.ReplaceAll("traceLength: ","");
         traceLength = sline.Atoi();
         std::cout << "traceLength = " << traceLength << std::endl;
      }
   }
   std::cout << std::endl;

}


