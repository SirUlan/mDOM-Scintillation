#include "mdomDetectorConstruction.hh"
#include "mdomPhysicsList.hh"
#include "mdomPrimaryGeneratorAction.hh"
#include "mdomRunAction.hh"
#include "mdomEventAction.hh"
#include "mdomTrackingAction.hh"
#include "mdomSteppingAction.hh"
#include "mdomSteppingVerbose.hh"
#include "mdomAnalysisManager.hh"
#include "mdomStackingAction.hh"
#include "mdomScintillation.hh"
#include "mdomPMT.hh"

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4ThreeVector.hh"

#include "G4UItcsh.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"	//xxx

#include "argtable2.h"
#include <ctime>
#include <sys/time.h>
#include <time.h>

#include "Randomize.hh"
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>

#include <cmath>	// for abs() of doubles
// #include "G4SystemOfUnits.hh"	// if compiling with Geant4.10

unsigned int	stats_buffer_max_size = 10;	// how many hits to keep in memory before purging to file in EndOfEventAction
unsigned int	stats_buffer_last_purge_at = 0;	// at what hits count was the hits file last written to
std::vector<G4int>	stats_PMT_hit;
std::vector<G4int>	stats_OM_hit;

G4double	gworldsize;
G4int		gsimevents;
G4String	ggunfilename;
G4String	ghitsfilename;
G4int		gPMT;
G4int		gGlass;
G4int		gGel;
G4double	gRefCone_angle;
G4int		gConeMat;
G4int		gHolderColor;
G4int		gDOM;
G4int		gEnvironment;
G4bool		gVisual;
G4bool		gInteractive;
G4String	gHittype;

G4double	gscintYield;
G4double	gscintTimeConst;
G4double	gscintSpectrum;
G4double	gTemperature;
G4double	gdistanceToSource;
G4double	gSourceRadius;
G4double	gElectronFactor;
G4int		gDecayConditional;
G4int		gQE;
G4int		gNucleus;
// G4String	greffilename;

G4bool		gKillAll;
G4long		current_event_id;
G4double gSimulatedTime;

struct timeval	gTime_Run_Start;
struct timeval	gTime_Run_End;
long randseed;

MdomAnalysisManager gAnalysisManager;
mdomPMT gPMTAnalysis;


G4UImanager* UI;

void clearstats() {
	stats_PMT_hit.clear();
	stats_OM_hit.clear();
}

std::vector<G4String> explode(G4String s, char d) {
	std::vector<G4String> o;
	int i,j;
	i = s.find_first_of("#");
	if (i == 0) return o;
 	while (s.size() > 0) {
		i = s.find_first_of(d);
		j = s.find_last_of(d);
		o.push_back(s.substr(0, i));
		if (i == j) {
			o.push_back(s.substr(j+1));
			break;
		}
		s.erase(0,i+1);
 	}
	return o;// o beinhaltet s ohne d
}
std::vector<G4String> explode(char* cs, char d) {
	std::vector<G4String> o;
	G4String s = cs;
	return explode(s,d);
}
std::vector<double> readColumnDouble (G4String fn, int col) {
	std::vector<double>	values;
	unsigned int c;
	double	a;
	c = col;
	std::ifstream	infile;
	std::vector<G4String> n;
	char l[256];
	G4String l2;
	infile.open(fn);
	while (infile.good() && !infile.eof()) {
		infile.getline(l,255);
		l2 = l;
		n = explode(l2,'\t');
		if (n.size()>=c) {
			a = atof(n.at(c-1));
			values.push_back(a);
		}
	}
	infile.close();

	return values;//values enthält den c. Wert aus fn (aus jeder Spalte,welche  nach 255 zeichen oder durch \n beendet wird?)
}


void Isotope_GPS(G4String Isotope) {
  UI = G4UImanager::GetUIpointer();
//       if ( Isotope == "U238" ){
//       UI->ApplyCommand("/control/execute ../DecayFiles/IonGui/NoDecayTime/U238FC.gui");
//     }
//     if ( Isotope == "U235" ){
//       UI->ApplyCommand("/control/execute ../DecayFiles/IonGui/NoDecayTime/UFC.gui");
//     }
//     if ( Isotope == "Th232" ){
//       UI->ApplyCommand("/control/execute ../DecayFiles/IonGui/NoDecayTime/Th232FC.gui");
//     }
//     if ( Isotope == "K40" ){
//       UI->ApplyCommand("/control/execute ../DecayFiles/IonGui/NoDecayTime/K.gui");
//     }

    UI->ApplyCommand("/control/verbose  0");
    UI->ApplyCommand("/run/verbose 0");
    UI->ApplyCommand("/event/verbose    0");
    UI->ApplyCommand("/tracking/verbose 0");
    UI->ApplyCommand("/run/initialize");
    UI->ApplyCommand("/process/em/fluo true");
    UI->ApplyCommand("/process/em/auger true");
    UI->ApplyCommand("/process/em/pixe true");
    UI->ApplyCommand("/gps/particle GenericIon");

    if ( Isotope == "U238" ){
      UI->ApplyCommand("/gps/ion 92 238 0");
      
    }
    if ( Isotope == "U235" ){
      UI->ApplyCommand("/gps/ion 92 235 0");
      
    }
    if ( Isotope == "Th232" ){
      UI->ApplyCommand("/gps/ion 90 232 0");
      
    }
    if ( Isotope == "K40" ){
      UI->ApplyCommand("/gps/ion 19 40 0");
      
    }
    UI->ApplyCommand("/mdomDECAYcmd/fullChain 1");
    
     //std::ostringstream osstring;
    //osstring << -0.5*gdistanceToSource;
    
    //G4String halfz = osstring.str();
    //G4cout << "Origin point for emission: 0 0 " << halfz << "\n\n\n" << G4endl;
    UI->ApplyCommand("/gps/pos/centre 0 0 0 m");
    //UI->ApplyCommand("/gps/particle e-");
    UI->ApplyCommand("/gps/ene/mono 0 eV");
    UI->ApplyCommand("/gps/pos/type Volume");
    UI->ApplyCommand("/gps/pos/shape Sphere");
    UI->ApplyCommand("/gps/pos/radius 0.5 m");// 240 mm");
    // UI->ApplyCommand("/process/inactivate Scintillation");
// 
//     UI->ApplyCommand("/process/inactivate Cerenkov");
    UI->ApplyCommand("/gps/ang/type iso");
    UI->ApplyCommand("/gps/pos/confine World_phys"); //Glass_phys");
}

int mdom_K40() {
        
	
	struct timeval time_for_randy;
	gettimeofday(&time_for_randy, NULL);

	randseed = time_for_randy.tv_sec+4294*time_for_randy.tv_usec;

	//CLHEP::HepRandom::setTheEngine(new CLHEP::RanluxEngine(randseed,4));
	//CLHEP::HepRandom::setTheEngine(new CLHEP::MixMaxRng);
	CLHEP::HepRandom::setTheEngine(new CLHEP::MTwistEngine(randseed));
	
	std::stringstream command;

	G4RunManager* runManager = new G4RunManager;

	mdomDetectorConstruction* detector;
	detector = new mdomDetectorConstruction();
	runManager->SetUserInitialization(detector);
	
	G4VUserPhysicsList* physics = new mdomPhysicsList;
	runManager->SetUserInitialization(physics);

	  
	#ifdef G4VIS_USE
 		G4VisManager* visManager = new G4VisExecutive;
 		visManager->Initialize();
 		visManager->SetVerboseLevel(0);
	#endif

	G4VUserPrimaryGeneratorAction* gen_action = new mdomPrimaryGeneratorAction();
	runManager->SetUserAction(gen_action);
	
	G4UserStackingAction* stacking_action = new mdomStackingAction();
	runManager->SetUserAction(stacking_action);

	G4UserRunAction* run_action = new mdomRunAction();
	runManager->SetUserAction(run_action);

	G4UserEventAction* event_action = new mdomEventAction();
	runManager->SetUserAction(event_action);

 	G4UserTrackingAction* tracking_action = new mdomTrackingAction();
 	runManager->SetUserAction(tracking_action);

	G4UserSteppingAction* stepping_action = new mdomSteppingAction();
	runManager->SetUserAction(stepping_action);

	runManager->Initialize();
	

	UI = G4UImanager::GetUIpointer();

	command.str("");
	command << "/control/execute " << ggunfilename;
	UI->ApplyCommand(command.str());
	
	G4double Activities[] = {0,61*13.,4.61*13,0.59*13,1.28*13};
	G4int NrOfDecays[] = {0,0,0,0,0};
	G4double TimeInterval = 60;
	gSimulatedTime = TimeInterval*s;

	G4String Decay_Isotopes[5] = {"none","K40","U238","U235","Th232"};
	G4String Isotope;
        
        double startingtime= clock() / CLOCKS_PER_SEC;
	for ( int i = 0; i < (int) gsimevents; i++) {
	  
// 	for ( int k = 1; k < (int) 5; k++) {
// 	 NrOfDecays[k] = G4int(G4Poisson(TimeInterval*Activities[k]));
// 	 G4cout << NrOfDecays[k] << G4endl;
// 	}
// 	
// 	
// 	for ( int k = 1; k < (int) 5; k++) {
// 	  Isotope = Decay_Isotopes[k];
// 	  //G4cout << "Simulating " <<  Isotope << " " << NrOfDecays[k] << G4endl;
// 	  Isotope_GPS(Isotope);
//           command.str("");
// 	  command << "/run/beamOn " << NrOfDecays[k];
// 	  UI->ApplyCommand(command.str());
//           
// //           for (int l = 0; l <  NrOfDecays[k]; l++){
// //               G4double tiii = (G4UniformRand()*gSimulatedTime)/ns;
// // 
// //           Isotope_GPS(Isotope);
// //           //G4cout << "/gps/time "<< std::setprecision(20) << tiii/s << " s";//std::setprecision(20) << tiii/s  << "\t";
// //           //command.str("");
// // 	  //command << "/gps/time "<< std::setprecision(20) << tiii/ns << " ns" ;
// //           //G4cout << command.str() << G4endl;
// //           //UI->ApplyCommand(command.str());
// //           command.str("");
// // 	  command << "/run/beamOn 1";
// // 	  UI->ApplyCommand(command.str());
// //           }
//           
// 	  
// 	}
            
        Isotope_GPS("K40");
        command.str("");
        command << "/run/beamOn 1";
        UI->ApplyCommand(command.str());
        G4String name1= ghitsfilename;
	gAnalysisManager.datafileTest.open(name1.c_str(), std::ios::out|std::ios::app);
	gAnalysisManager.WriteAccept();
	gAnalysisManager.hits_all_events.clear();
        gAnalysisManager.Reset();
	gAnalysisManager.datafileTest.close();
	}
	double finishtime=clock() / CLOCKS_PER_SEC;
	G4cout << "Computation time: " << finishtime-startingtime << " seconds." << G4endl;
	// Opens new user interface prompt and visualization after simulation was run
	if (gInteractive){
		int argumc = 1;
		char* argumv[] = {"dummy", NULL};
		G4UIExecutive* UIEx = new G4UIExecutive(argumc, argumv);
		if (gVisual){
			UI->ApplyCommand("/control/execute ../aux/init_vis.mac");
		}
		UIEx->SessionStart();
		delete UIEx;
	}

#ifdef G4VIS_USE
	delete visManager;
#endif

	delete runManager;
	return 0;
}

int main(int argc,char *argv[])
{
	struct arg_dbl  *worldsize	= arg_dbl0("wW", "world","<n>","\t\tradius of world sphere in m");
	struct arg_int  *events		= arg_int0("nN", "numevents,nevents","<n>","\tnumber of decays per run");
	struct arg_file *gunfile	= arg_file0("gG","gun","<file.txt>","\t\tfile containing GPS parameters");
	struct arg_int  *pmt		= arg_int0("pP", "pmt,PMT","<n>","\t\tPMT type [12199S, etel, 12199e]");
	
	struct arg_int  *glass		= arg_int0("uU", "glass","<n>","\t\t\tglass type [VITROVEX, Chiba, Kopp, myVitroVex, myChiba, WOMQuartz, fusedSilica]");
	struct arg_int	*gel 		= arg_int0("jJ", "gel", "<n>", "\t\t\tgel type [WACKER, Chiba, IceCube, Wacker_company]");
	struct arg_dbl  *cone_ang   = arg_dbl0("aA", "cone_ang","<n>","\t\t\topening semi-angle of cone; (45 deg)");	
	struct arg_int	*conemat 	= arg_int0("kK", "conemat", "<n>", "\t\t\tcone material [V95, v98, aluminium, total98]");
	struct arg_int	*holdercol 	= arg_int0("cC", "holdercol", "<n>", "\t\t\tcone color [BLACK, white (Lambertian R = 98%)]");
	struct arg_int	*dom 		= arg_int0("mM", "om, dom", "<n>", "\t\t\tmodule type [MDOM, PDOM,HalfMDOM]");
	
	struct arg_int  *environment= arg_int0("eE", "environment","<n>","\t\tmedium in which the setup is emmersed [AIR=0, ice=1, spice=2]");
	struct arg_file *outputfile	= arg_file0("oO","output","<file.txt>","\t\tfilename for hits data");
	struct arg_int  *hittype	= arg_int0("hH", "hits","<n>","\t\thit collection [onlyHits, detailed]");
	struct arg_lit	*interactive = arg_lit0("iI","interact","\t\topens user interface after run");
	struct arg_lit	*visual		= arg_lit0("vV","visual","\t\tshows visualization of module after run (also calls interactive)");
	struct arg_dbl	*scintYield	= arg_dbl0("yY", "scintYield", "<n>", "\t\tScintillation Yield of the glass (only Vitrovex). Default 20/MeV");
	struct arg_dbl	*scintTimeConst	= arg_dbl0("tT", "scintTimeConst", "<n>", "\t\tScintillation's Time constant of the glass (only Vitrovex) in ns. Default 300000.");
	struct arg_dbl	*scintSpectrum	= arg_dbl0("sS", "scintSpectrum", "<n>", "\t\tMove the scintillation's spectrum by # nm. Default 0 nm.");
	struct arg_dbl	*Temperature 	= arg_dbl0("bB", "Temperature", "<n>", "\t\t Temperature for material property selection");
	struct arg_dbl	*distanceToProbe	= arg_dbl0("dD", "distanceToProbe", "<n>", "\t\t Distance from PMT to source");
	struct arg_dbl	*SourceRadius	= arg_dbl0("rR", "SourceRadius", "<n>", "\t\t SourceRadius");
	struct arg_dbl	*ElectronFactor	= arg_dbl0("fF", "ElectronFactor", "<n>", "\t\t ElectronFactor");
	struct arg_lit	*help		= arg_lit0(NULL,"help","\t\tprint this help and exit");
	struct arg_int	*nucleus	= arg_int0(NULL,"nucleus","<n>","\t\t decay chain or isotope [1=K40, 2=U238, 3=U235, 4=Th232]");
	struct arg_int *QE = arg_int0("qQ","QE","<n>","\t\tQuantum efficiency ON = 1 or off =0. Killing events. Default 0.");
	struct arg_end  *end		= arg_end(22);
	
	void* argtable[] = {worldsize, 
						events, 
						gunfile, 
						pmt,
						
						glass,
						gel,
						cone_ang,
						conemat,
						holdercol,
						dom,
						
						environment, 
						outputfile, 
						hittype, 
						interactive, 
						visual, 
						scintYield,
						scintTimeConst,
						scintSpectrum,
						Temperature,
						distanceToProbe,
						SourceRadius,
						ElectronFactor,
						help,
						nucleus,
						QE,						
						end};
						
	const char* progname = "mdom_K40";
	int nerrors;
	int exitcode=0;

	// verify the argtable[] entries were allocated sucessfully
	if (arg_nullcheck(argtable) != 0) {
		/* NULL entries were detected, some allocations must have failed */
		printf("%s: insufficient memory\n",progname);
		exitcode=1;
		arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));
		return exitcode;
	}

	// set any command line default values prior to parsing
	worldsize->dval[0] = 10.0;
	events->ival[0] = 0;
	gunfile->filename[0] = "mdom_K40.gps";
	pmt->ival[0] = 0;	// use new R12199 as default
	
	glass->ival[0] = 0;	// use VITROVEX as default
	gel->ival[0] = 0;	// use Wacker SilGel 612 A/B as default
	cone_ang->dval[0] = 45.0; // [degrees]	
	conemat->ival[0] = 0;	// use Alemco V95 as default
	holdercol->ival[0] = 0;	// use classic black holder as default
	dom->ival[0] = 0;	// use mDOM as default
	
	environment->ival[0] = 0;	// use air as default
	outputfile->filename[0] = "K40";
	hittype->ival[0] = 0;	// store only hits as default
	scintYield->dval[0] = 20.0;
	scintTimeConst->dval[0] = 100.0;
	scintSpectrum->dval[0] = 0.0;
	Temperature->dval[0] = -35;
	distanceToProbe->dval[0] = 7.3;
	SourceRadius->dval[0] = 0.0;
	ElectronFactor->dval[0] = 9.5;
	QE->ival[0]=1;
	nucleus->ival[0] = 0;
	/* Parse the command line as defined by argtable[] */
    nerrors = arg_parse(argc,argv,argtable);

    /* special case: '--help' takes precedence over error reporting */
    if (help->count > 0)
	{
        printf("\nGEANT4 simulation of the mDOM: K40 decays\n");
        printf("\nUsage: %s", progname);
        arg_print_syntax(stdout,argtable,"\n");
        arg_print_glossary(stdout,argtable,"  %-25s %s\n");
        printf("\n");
        exitcode=0;
        goto hell;
	}

    /* If the parser returned any errors then display them and exit */
    if (nerrors > 0)
	{
        /* Display the error details contained in the arg_end struct.*/
        arg_print_errors(stdout,end,progname);
        printf("Try '%s --help' for more information.\n",progname);
        exitcode=1;
        goto hell;
	}

    /* special case: uname with no command line options induces brief help */
    if (argc==1)
	{
        printf("Try '%s --help' for more information.\n",progname);
        exitcode=0;
        goto hell;
	}

//	assign command-line arguments to variables:
	gworldsize = worldsize->dval[0];
	gsimevents = events->ival[0];
	ggunfilename = gunfile->filename[0];
	gPMT = pmt->ival[0];
	
	gGlass = glass->ival[0];
	gGel = gel->ival[0];
	gRefCone_angle = cone_ang->dval[0];		
	gConeMat = conemat->ival[0];
	gHolderColor = holdercol->ival[0];
	gDOM = dom->ival[0];
	
	gEnvironment = environment->ival[0];
	ghitsfilename = outputfile->filename[0];
	gscintYield = scintYield->dval[0];
	gscintTimeConst = scintTimeConst->dval[0];
	gscintSpectrum = scintSpectrum->dval[0];
	gTemperature = Temperature->dval[0];
	gdistanceToSource = distanceToProbe->dval[0];
	gSourceRadius = SourceRadius->dval[0];
	gElectronFactor = ElectronFactor->dval[0];
	gQE = QE->ival[0];
	
	gNucleus = nucleus->ival[0];
	
	if (hittype->ival[0]==0){
		gHittype = "onlyHits";
	}
	if (hittype->ival[0]==1){
		gHittype = "detailed";
	}
	if (interactive->count > 0) gInteractive = true; else gInteractive = false;
	if (visual->count > 0) {
		gVisual = true;
		gInteractive = true;
	}
	else {
		gVisual = false;
	}

	//	check params for sanity
	gPMTAnalysis.loadThePMTInfo();
	mdom_K40();

hell:
    /* deallocate each non-null entry in argtable[] */
	arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));
	//	return exitcode;

}

