/*======================================================================
 LO, NLO, NNLO rapidity distributions for 
      p p(bar) -> q qbar (q) -> (Z,gamma^*) + X,
                  q g -> (Z,gamma^*) + q + X
                  g g -> (Z,gamma^*) + X
 BOOST & REAL INTEGRANDS ARE COMBINED INTO ONE.
 ====================================================================*/

#include <iomanip>
#include "QCDbasics.h"
#include "NAClasses.h"
#include "DClasses.h"
#include "dilog.h"
#include "BornNLOfns.h"  // Born kinematics & NLO cross section fns.
#include "qbarqfns.h"  // NNLO hard cross section fns. for q-\bar{q} channel
#include "qgfns.h"  // NNLO hard cross section fns. for q-g channel
#include "ggfns.h"  // NNLO hard cross section fns. for g-g channel
#include "qiqjfns.h"  // NNLO hard cross section fns. for q_i \neq q_j channel
#include <fstream>
#include "settings.h"
#include "Vlumifns_LHApdf.h" // luminosity functions required for Drell-Yan
#include "LHApdf.h"
#include "mode.h"
#include "integration.h"

#include "options.h"


class VrapOptionsHandler : public OptionsHandler {
    public:
        VrapOptionsHandler();
            virtual ~VrapOptionsHandler(){};
    };

VrapOptionsHandler::VrapOptionsHandler(){
	add(new ValueSettingOption<double>("E_CM",E_CM,"Sets the hadron-hadron center of mass energy in GeV units."));
	add(new ValueSettingOption<double>("Q",Q,"Sets the vector boson mass in GeV units."));
	add(new ValueSettingOption<double>("muFoverQ",muFrel,"Sets the ratio muF/Q."));
	add(new ValueSettingOption<double>("muRoverQ",muRrel,"Sets the the ratio muR/Q."));
	add(new ValueSettingOption<double>("Nf",Nf,"Sets the number of quark flavors."));
	add(new ValueSettingOption<double>("Alphat",alphat,"Sets the value of alpha."));
	add(new ValueSettingOption<int>("RandomSeed",ranseed,"Sets the random seed, must be a positive integer."));
	add(new multipleValueOption<collider>("Collider",coll,"pp",pp,"ppbar",ppbar,"piso",piso,"Sets the type of collider. ") );
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
        add(new multipleValueOption<int>("OutputFormat",o_f,"TopDrawStyle",0,"ListValues",1,"Sets the style output is printed (for td or just list dsig/dy). ") );
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
    std::cout << "======================================================== " << std::endl;
    std::cout << "=                       Vrap-0.9                       = " << std::endl;
    std::cout << "=                                                      = " << std::endl;
    std::cout << "=       Computes rapidity distributions for            = " << std::endl;
    std::cout << "=   production of electroweak vector bosons in         = " << std::endl;
    std::cout << "=     hadronic collisions through NNLO in QCD          = " << std::endl;
    std::cout << "=                                                      = " << std::endl;
    std::cout << "=                 Based on the article                 = " << std::endl;
    std::cout << "=    `High precision QCD at hadron colliders:          = " << std::endl;
    std::cout << "=  Electroweak gauge boson rapidity distributions      = " << std::endl;
    std::cout << "=  at NNLO', by Charalampos Anastasiou, Lance Dixon,   = " << std::endl;
    std::cout << "=  Kirill Melnikov and Frank Petriello,                = " << std::endl;
    std::cout << "=  Phys. Rev. D69, 094008 (2004) [hep-ph/0312266]      = " << std::endl;
    std::cout << "=                                                      = " << std::endl;
    std::cout << "= This version includes an interface to parton         = " << std::endl; 
    std::cout << "= distribution functions via the Les Houches Accord    = " << std::endl;
    std::cout << "= (LHAPDF). This interface, and several other program  = " << std::endl;
    std::cout << "= improvements, provided by Daniel Maitre              = " << std::endl;
    std::cout << "=                                                      = " << std::endl;
    std::cout << "=  Questions?  lance at slac dot stanford dot edu      = " << std::endl;
    std::cout << "======================================================== " << std::endl; 
    std::cout << std::endl;
}
	
void printParamInfo(){	
    std::cout << "-------------------------------------------------------- " << std::endl; 
    std::cout << " ranseed =  " << ranseed << std::endl; 
    std::cout << " order_flag =  " << order_flag << std::endl; 
    std::cout << " E_CM =  " << E_CM 
        << " GeV;   collider = " << coll << "  (1 = pp;  2 = ppbar)"  << std::endl; 
    std::cout << " alpha_s(M_Z) =  " << alpha_s_Z << std::endl; 
    std::cout << " Q =  " << Q << "   muR/Q = " << muRrel 
                        << "   muF/Q = " << muFrel << std::endl; 
    std::cout << " alpha_QED used =  " << alphat << std::endl; 
    std::cout << " Nf =  " << Nf << "     NFf = " << NFf(int(Nf)) << std::endl;
    std::cout << " alpha_s(muR)/Pi =  " << alpha_s(muR)/PI << std::endl;
    std::cout << " DY prefactor = " << DY_prefactor(Q,alphat) << std::endl;
    if (useOtherPDF){
        std::cout << " PDF file = " << pdfFile << std::endl;
        std::cout << " PDF set = " << pdfSet << std::endl;
    } else {
        std::cout << " PDF mode = " << pdfMode << std::endl;
    }

    std::cout << "-------------------------------------------------------- " << std::endl; 
};

//===========  MAIN PROGRAM  ==========================================

int main(int argc,char* argv[]){

	printBanner();

	defaultSettings();
	
	VrapOptionsHandler VOH;
    std::string filename;
	if (argc >=2 ){
		filename = std::string(argv[1]);
		if (filename == "help" || filename[0]=='-'){
			std::cout << "Here is the list of options that can be set in the input file\n " <<std::endl;
			VOH.printHelp(std::cout);
			return 0;
		}
		std::cout << "Reading input from: " << filename << std::endl;
	} else {
        std::cerr << "Usage: Vrap InputFile\n" << std::endl;
		return 1;
	}
	
    std::ifstream inputfile(filename.c_str());

	if (!inputfile){
        std::cerr << "Could not open " << filename <<"." << std::endl;
		return 2;
	}

    std::string message;
	bool success = VOH.process_file(inputfile,message);
	if (!success){
        std::cerr << "Parsing of input file " << filename << " failed: " << message << std::endl;
		return 3;
	}

	if(argc>2) {
	  Q = strtod(argv[2], NULL);
	}
	double rapY = 0;
	if(argc>3) {
	  rapY = strtod(argv[3], NULL);
	}

	if(argc>4) {
		pdfSet = strtod(argv[4], NULL);
	}


    switch(parton_flag){
        // From `Vlumifns.h`, sets up which parton channel are active
        // flags: f_qqbar, f_qg, f_gg, f_qq;
        case 1: compute_all(); break;	
        case 2: compute_qqbar(); break;	
        case 3: compute_qg(); break;	
        case 4: compute_gg(); break;	
        case 5: compute_qq(); break;	
        case 6: compute_qqbar_plus_qq(); break;	
    }

    // alphat = alpha_QED(Q);   // my running QED coupling, better for low-mass
    if (useMyAlphaRunning) {
        alphat = alpha_QED(Q);
    }

    // Set f_NNLO_only = 0 (1)  for NNLO terms only (LO+NLO+NNLO)
    if (NNLO_only){
        f_NNLO_only = 0;
    } else {
        f_NNLO_only = 1;	
    }

    // Compute all couplings for all luminosity channels
    // these are then used by the different functions inside `Vlumifns.C`
    setV(exchM,Q,alphat,Nf,0);


// Uncomment for fixed-target Drell-Yan, M vs. y. table of K factors
// (overwrites previous!):
// coll = pp;   E_CM = 38.757;  Q = 8.;  Nf = 5.;   muF = Q;  muR = muF; 
// f_NNLO_only = 0;  alphat = 1./132.1;   setV(gamma_only,Q,alphat,Nf,0);   


    if (useOtherPDF){
        // initializes LHAPDF
        set_mode(pdfFile, pdfSet);
    } else {
        set_mode(pdfMode);
    }

   muR = muRrel*Q;
   muF = muFrel*Q;

   printParamInfo();

//double alpha_s_Z_PDF;
// LHAPDF::initPDFSetByName("abkm09_5_nlo.LHgrid");
// alpha_s_Z_PDF = alphasPDF(m_Z); 
// std::cout << " alpha_s(M_Z) =  " << alpha_s_Z_PDF << std::endl; 

// To do the total cross section integral, for normalization/checks:
// DVector norm_factor = int_y(0.,1.);

// To scan the total cross section in (muR = muF), through
//    Q*mu_r_lower < (muR=muF) < Q*mu_r_upper     with n_points steps:
// DMatrix mutotresultMatrix = scan_mu_tot(0.2,5.0,10);

// Compute and print unnormalized rapidity distributions, 
// with arguments (ymax,n_points)

// "sym_scan_rap_y(n_points,forwrev,o_f)"
// scans through 0 < y < ymax with n_points steps,
// for the case of a symmetric distribution.
// forwrev = 1  -> print output from -ymax to +ymax
// forwrev = -1  -> print output from +ymax to -ymax
// o_f = output_format = 0  -> print for topdraw: y, d sigma/dy  ( error
// o_f = output_format = 1  -> just list d sigma/dy values.
//
//DMatrix resultMatrix = sym_scan_rap_y(nbrYPnts,direction,o_f); // "normally" 19   /// <- this was uncommented MB
// Case of an asymmetric distribution, scans from -ymax to +ymax:
// DMatrix resultMatrix = asym_scan_rap_y(nbrYPnts,direction,o_f); // "normally" 38

// Scan all curves at once (symmetric case):
// sym_scan_all(alekhin,nbrYPnts,0.5,2.0);
// Scan all curves at once (asymmetric case):
// asym_scan_all(alekhin,nbrYPnts,0.5,2.0);
// Scan all W+/W- ratio curves at once (symmetric case):
// Wpm_ratio_scan_all(alekhin,nbrYPnts,0.5,2.0);

// Just one (two?) rapidity point, y:
// order_flag = 2;

    y = rapY;
    std::cout << "\n > Starting calculation:\n\n";
    DVector temp_ans = rap_y();
    std::cout << "\nFinal result: " << temp_ans[0] << " +/- " << temp_ans[1] << std::endl;
    std::cout << temp_ans[0] << std::endl;

// At fixed rapidity y, scan d^2sigma/dM/dY through
//     Q*mu_r_lower < (muR=muF) < Q*mu_r_upper     with n_points steps:
// arguments are: (y, mu_r_lower,mu_r_upper,n_points).
// scan_mu(0.0,0.2,5.0,10);

// Same as "scan_mu" except for W+/W- ratio:
// Wpm_ratio_scan_mu(2.0,0.2,5.0,10);

// Make table of K factors for J. Stirling:
// makeKtable(1.0, 5, 4.0, 16.0, 6);
 
// Make table of "C" factors for J. Stirling:
// makeCtable(1.0, 5, 4.0, 16.0, 6);
// makeCtable(1.0, 20, 14.0, 20.0, 3);

// Make table of "D" and "E" factors for J. Stirling:
//  makeDEtable(1.0, 5, 4.0, 16.0, 6);

// Integral over M (Simpson's rule; scan in y):
// int_M_sym_scan_y(Zgamma,66.,116.,200,3,0);

// Monte Carlo integral over M from Ml to Mu, at fixed y:
 // y = 0.;   Ml = 66.;  Mu = 116.;  exchM = Zgamma;  
 // DVector temp_y_M = rap_y_M();

// Monte Carlo integral over M from Ml to Mu, symmetric y scan (n_points):
// sym_scan_rap_y_M(Zgamma,66.,116.,nbrYPnts,direction);

// Monte Carlo integral over M from Ml to Mu, asymmetric y scan (n_points):
// asym_scan_rap_y_M(Wplus,60.,100.,nbrYPnts,direction);

// hard_fn_test();
// lumi_test();

/*
// Use int_y(xi_l,xi_u) to integrate over a specified rapidity interval
 y_lower = 0.0; 
 y_upper = 1.0; 
 xi_l = 0.5 * (1. + tanh(y_lower));  // range in terms of xi variable.
 xi_u = 0.5 * (1. + tanh(y_upper));
 DVector temp_int_y = int_y(xi_l, xi_u);
*/

 return 0;
}
