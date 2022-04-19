/*======================================================================
 LO, NLO, NNLO TOTAL CROSS SECTION distributions for 
      p p(bar) -> q qbar (q) -> (Z,gamma^*) + X,
                  q g -> (Z,gamma^*) + q + X
                  g g -> (Z,gamma^*) + X
      p p(bar) -> q qbar -> DY + X,
 NNLO formulae from Hamberg, van Neerven, & Matsuura
 ====================================================================*/

#include <iomanip>
#include <fstream>
#include "QCDbasics.h"
#include "NAClasses.h"
#include "DClasses.h"
#include "dilog.h"
#include "pdf.h"  // pdf distributions
#include "Vlumifns.h"  // luminosity functions required for Drell-Yan
#include "Vlumifns_LHApdf.h"  // luminosity functions required for Drell-Yan
#include "totfns.h"
#include "mode.h"

using namespace std;

double E_CM;         // hadron-hadron center-of-mass energy
double Q, muR, muF;  // DY-mass; renormalization & factorization scales
double muRrel, muFrel;  // global muR/Q, muF/Q (needed to do M (Q) integral)
double alphat;  // local value of alpha_QED(Q)
exchange exchM;  // value of "exchange" for doing M (Q) integral
double Nf;  // Nf = number of light quarks flavors; should ~ depend on Q.
collider coll;
int ranseed;
bool NNLO_only;
int f_NNLO_only;  // fudge parameter to isolate NNLO terms
                  // normally should be 1
                  // Set to 0 for only the NNLO hard cross section terms
int f_quiet = 0; // Set to 1 to suppress intermediate printing by setV().


std::string pdfMode;
std::string pdfFile;
bool useMyAlphaRunning=false;
bool useOtherPDF;
int pdfSet;
int parton_flag; // 1: all, 2:qqbar 3: qg 4: gg 5:qq 6: qqbar_plus_qq 

int direction; int nbrYPnts;




/*================================================================
 The full (NLO, NNLO) DY total calculation divides up into 2 pieces:
 1) A Born kinematics piece (Born), with only DY in final state, and
  proportional to delta(1-z).
 2) A boosted Born kinematics pieces (boost).
=================================================================*/

/*=================================================================
  The overall q qbar -> DY normalization:
  The "Born-level" partonic cross section in pb is sigma_DY_0.
  It omits the quark charges, which will be supplied elsewhere.
  (Now sigma_DY_0 is factored out front of everything.)
  The "2/M" normalization gives d sigma/dM.                     */


#include "integration_tot.h"

#include "options.h"


class VtotOptionsHandler : public OptionsHandler {
public:
VtotOptionsHandler();
	virtual ~VtotOptionsHandler(){};
};

VtotOptionsHandler::VtotOptionsHandler(){
	add(new ValueSettingOption<double>("E_CM",E_CM,"Sets the hadron-hadron center of mass energy in GeV units."));
	add(new ValueSettingOption<double>("Q",Q,"Sets the vector boson mass in GeV units."));
	add(new ValueSettingOption<double>("muFoverQ",muFrel,"Sets the ratio muF/Q."));
	add(new ValueSettingOption<double>("muRoverQ",muRrel,"Sets the the ratio muR/Q."));
	add(new ValueSettingOption<double>("Nf",Nf,"Sets the number of quark flavors."));
	add(new ValueSettingOption<double>("Alphat",alphat,"Sets the value of alpha."));
	add(new ValueSettingOption<int>("RandomSeed",ranseed,"Sets the random seed, must be a positive integer."));
	add(new multipleValueOption<collider>("Collider",coll,"pp",pp,"ppbar",ppbar,"Sets the type of collider. ") );
	add(new yesOrNoOption("NNLO_only",NNLO_only,"Sets whether to compute NNLO terms only."));
	add(new yesOrNoOption("UseMyAlphaQEDRunning",useMyAlphaRunning,"Sets whether to use my running QED coupling, better for low-mass"));
	add(new ValueSettingOption<std::string>("PDF_mode",pdfMode,"Sets the PDF mode."));
	add(new yesOrNoOption("UseOtherPDF",useOtherPDF,"Sets whether to use a different PDF setting"));
	add(new ValueSettingOption<double>("AlphasZ",alpha_s_Z,"Sets the value of alphas at the Z pole."));
	add(new multipleValueOption<int>("Order",order_flag,"LO",0,"NLO",1,"NNLO",2,"Sets the order in the hard cross section."));
	add(new ValueSettingOption<std::string>("PDFfile",pdfFile,"Sets the PDF file."));
	add(new ValueSettingOption<int>("PDFset",pdfSet,"Sets the PDF set."));
	add(new multipleValueOption<int>("InitialPartonChannels",parton_flag,"All",1,"qqbar",2,"qg",3,"gg",4,"qq",5,"qqbar_plus_qq",6,"Sets which parton channel(s) to take into account."));
	add(new multipleValueOption<exchange>("VectorBoson",exchM,"gamma_only",gamma_only,"Zgamma_interf",Zgamma_interf,"Z_only",Z_only,"Zgamma",Zgamma,"Wplus",Wplus,"Wminus",Wminus,"Sets which vector boson to use."));
	add(new ValueSettingOption<int>("NumberOfYPoints",nbrYPnts,"Sets the number of rapidity values at which to compute."));
	add(new multipleValueOption<int>("PrintDirection",direction,"Forward",+1,"Reverse",-1,"Sets the order in which output is printed (for td). ") );
	//enableDebug();
}


void defaultSettings(){
		
Nf = 5.;  // Default number of light quark flavors in hard cross section.
NNLO_only=false;


// Default value of alpha QED.
 alphat = 1./128.;   // USE FOR ALL Z/W COMPARISONS WITH FRANK!!!
// alphat = 1./132.1;  // Used for fxd. tgt. M=8 GeV in  hep-ph/0306192

// Uncomment for fixed-target Drell-Yan, M vs. y. table of K factors
// (overwrites previous!):
// coll = pp;   E_CM = 38.757;  Q = 8.;  Nf = 5.;   muF = Q;  muR = muF; 
// f_NNLO_only = 0;  alphat = 1./132.1;   setV(gamma_only,Q,alphat,Nf,0);   


	useOtherPDF=false;
	pdfSet=-1;
	parton_flag=1;
	nbrYPnts=19;
	direction=1;
	}


void printBanner(){
 cout << "======================================================== " << endl;
 cout << "=                       Vrap-0.8                       = " << endl;
 cout << "=                                                      = " << endl;
 cout << "=       Computes rapidity distributions for            = " << endl;
 cout << "=   production of electroweak vector bosons in         = " << endl;
 cout << "=     hadronic collisions through NNLO in QCD          = " << endl;
 cout << "=                                                      = " << endl;
 cout << "=                 Based on the article                 = " << endl;
 cout << "=    `High precision QCD at hadron colliders:          = " << endl;
 cout << "=  Electroweak gauge boson rapidity distributions      = " << endl;
 cout << "=  at NNLO', by Charalampos Anastasiou, Lance Dixon,   = " << endl;
 cout << "=  Kirill Melnikov and Frank Petriello,                = " << endl;
 cout << "=  Phys. Rev. D69, 094008 (2004) [hep-ph/0312266]      = " << endl;
 cout << "=                                                      = " << endl;
 cout << "= This version includes an interface to parton         = " << endl; 
 cout << "= distribution functions via the Les Houches Accord    = " << endl;
 cout << "= (LHAPDF). This interface, and several other program  = " << endl;
 cout << "= improvements, provided by Daniel Maitre              = " << endl;
 cout << "=                                                      = " << endl;
 cout << "=  Questions?  lance at slac dot stanford dot edu      = " << endl;
 cout << "======================================================== " << endl; 
 cout << endl;

	}
	
void printParamInfo(){	
 cout << "-------------------------------------------------------- " << endl; 
 cout << " ranseed =  " << ranseed << endl; 
 cout << " order_flag =  " << order_flag << endl; 
 cout << " E_CM =  " << E_CM 
      << " GeV;   collider = " << coll << "  (1 = pp;  2 = ppbar)"  << endl; 
 cout << " alpha_s(M_Z) =  " << alpha_s_Z << endl; 
 cout << " Q =  " << Q << "   muR/Q = " << muRrel 
                       << "   muF/Q = " << muFrel << endl; 
 cout << " alpha_QED used =  " << alphat << endl; 
 cout << " Nf =  " << Nf << "     NFf = " << NFf(int(Nf)) << endl;
 cout << " alpha_s(muR)/Pi =  " << alpha_s(muR)/PI << endl;
 cout << " DY prefactor = " << DY_prefactor(Q,alphat) << endl;
 if (useOtherPDF){
	cout << " PDF file = " << pdfFile << endl;
	cout << " PDF set = " << pdfSet << endl;
} else {
	cout << " PDF mode = " << pdfMode << endl;
}

cout << "-------------------------------------------------------- " << endl; 
};

//===========  MAIN PROGRAM  ==========================================


int main(int argc,char* argv[]){

	printBanner();

	defaultSettings();

		VtotOptionsHandler VOH;
	string filename;
	if (argc ==2 ){
		filename = string(argv[1]);
		if (filename == "help" || filename[0]=='-'){
			cout << "Here is the list of options that can be set in the input file\n " <<endl;
			VOH.printHelp(cout);
			return 0;
		}
		cout << "Reading input from: " << filename << endl;
	} else {
		cerr << "Usage: Vrap InputFile\n" << endl;
		return 1;
	}
	
	ifstream inputfile(filename.c_str());

	if (!inputfile){
		cerr << "Could not open " << filename <<"." << endl;
		return 2;
	}

	string message;
	bool success = VOH.process_file(inputfile,message);
	if (!success){
		cerr << "Parsing of input file " << filename << " failed: " << message << endl;
		return 3;
	}


switch(parton_flag){
	case 1: compute_all(); break;	
	case 2: compute_qqbar(); break;	
	case 3: compute_qg(); break;	
	case 4: compute_gg(); break;	
	case 5: compute_qq(); break;	
	case 6: compute_qqbar_plus_qq(); break;	
}


// alphat = alpha_QED(Q);   // my running QED coupling, better for low-mass
if ( useMyAlphaRunning){
	alphat = alpha_QED(Q);
}

 // Set f_NNLO_only = 0 (1)  for NNLO terms only (LO+NLO+NNLO)
if (NNLO_only){
 f_NNLO_only = 0;
} else {
 f_NNLO_only = 1;	
}

setV(exchM,Q,alphat,Nf,0);

if (useOtherPDF){
	set_mode(pdfFile,pdfSet);
} else {
	set_mode(pdfMode);
}

   muR = muRrel*Q;    muF = muFrel*Q;

printParamInfo();

	
	
// Also for Joey --- NNLO cross section with MRST2004 NLO partons 
 // but NNLO alpha_s:
 // order_flag = 2;  alpha_s_Z = 0.1167;  pdf_init(mrst,15);  
 // Alekhin distributions and alpha_s:  
 // Alekhin02 [LO,NLO,NNLO](kord=1,2,3), variable-flavor (kschem=1), 
 // nominal (kset=0), are obtained by setting order_flag = 0,1,2.
 //  order_flag = 2; pdf_init(alekhin,order_flag); alpha_s_Z = Alekhin_alpha_s(m_Z); 




 // cross-section at one value of muF, muR:
 DVector total_cross_section = int_y();

// scan_mu(mu_r_lower,mu_r_upper,n_points)
// scans through Q*mu_r_lower < (muR=muF) < Q*mu_r_upper  with n_points steps:
 // DMatrix muresultMatrix = scan_mu(0.5,2.0,2);

 return 0;
}
