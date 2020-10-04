#include "ColAHelL.h"
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <iomanip>

using namespace std;

#include "../RapidXML/rapidxml.hpp"
using namespace rapidxml;

static int SUCCESSX = 58;

ColAHelL::ColAHelL()
{
	this->Spect = 0;
	this->Proc = 0;
	this->CTof = 0;
	this->HGrams = 0;
	this->evt = 0;
//	this->Cfg = 0;
	this->PrsHist = 0;
	for(int j=0;j<64;j++) {
		this->presorters[j]=0; 
		this->e[j] = 0;
		this->i[j] = 0;
	}
	this->mol = 0;
	this->big_mol = 0;
}
ColAHelL::~ColAHelL()
{
	// delete all CH objects...
	if(this->Spect)
	{ delete this->Spect; this->Spect=0; } 
	if(this->Proc)
	{ delete this->Proc; this->Proc=0; }
	if(this->CTof)
	{ delete this->CTof; this->CTof=0; }
	if(this->HGrams)
	{ delete this->HGrams; this->HGrams=0; }
	if(this->evt)
	{ delete this->evt; this->evt=0; }
	if(this->PrsHist)
	{ delete this->PrsHist; this->PrsHist=0; }
	for(int i=0;i<num_presorters;i++) {
		if(this->presorters[i])
			{ delete this->presorters[i]; this->presorters[i]=0; }
	}
	for(int i=0;i<64;i++) {
		if(this->e[i])
			{ delete this->e[i]; this->e[i]=0; }
		if(this->i[i])
			{ delete this->i[i]; this->i[i]=0; }
	}
	if(this->mol)
	{ delete this->mol; this->mol=0; }
	if(this->big_mol)
	{ delete this->big_mol; this->big_mol=0; }
}

int ColAHelL::InitAll()
{
	return InitAll("ColAHelL.cfg", false);
}

int ColAHelL::InitAll(bool CrLogFile)
{
	return InitAll("ColAHelL.cfg", CrLogFile);
}

int ColAHelL::InitAll(const char* FName, bool CrLogFile=false)
{
	//switch cout to file if file log requested
	if(CrLogFile) {
		filestr.open ("CH-log.txt");
		backup = cout.rdbuf();     // backup cout's streambuf
		psbuf = filestr.rdbuf();   // get file's streambuf
		cout.rdbuf(psbuf);         // assign streambuf to cout
	}

	this->Status = 1;
	
	cout << setprecision(3);
	White(true);
	cout << "\nSetting up...\n";
	cout << "---------------------------------------------------\n";
	Red(true);
	cout << "                             )          (     \n";
	cout << "   (         (     (      ( /(       (  )\\ )  \n";
	cout << "   )\\        )\\    )\\     )\\())   (  )\\(()/(  \n";
	cout << " (((_)   (  ((_)((((_)(  ((_)\\   ))\\((_)/(_)) \n";
	cout << " )\\___   )\\  _   )\\ _ )\\  _((_) /((_)_ (_))   \n";
	cout << "((/ __| ((_)| |  (_)_\\(_)| || |(_)) | || |    \n";
	White(true);
	cout << " | (__ / _ \\| |   / _ \\  | __ |/ -_)| || |__  \n";
	Yellow(true);
	cout << "  \\___|\\___/|_|  /_/ \\_\\ |_||_|\\___||_||____| \n";
	cout << "                                              \n";
	White(false);
	printf("                                      Version: %.2f\n",CH_VERSION);
	White(true);
	cout << "---------------------------------------------------\n";
	White(false);

	SetupEventStruct();
	SetupParticleClasses();
	if(!SetupHistogramming()) {
		this->Status = -30;
	}

	if(!OpenXML((char*)FName)) {
		
		
		// restore cout
		if(CrLogFile) {
			cout.rdbuf(backup);
			filestr.close();
		}		

		this->Status = -99;
		return this->Status;
	}
	
	cout << "\n\n\n--====================================================================";
	if(this->Status>0) {
		cout << "\nSetting up ColAHelL: ";
		Green(true);
		gotoX(SUCCESSX); 
		cout << " ...success!\n";
		White(false);	
	}	
	if(this->Status<0 && this->Status>-90) {
		Yellow(true);
		cout << "\nTHERE WERE WARNINGS while setting up ColAHelL!\n Please check the log above. Some features might now work!\n";
		White(false);	
	}
	if(this->Status<-89) {
		Red(true);
		cout << "\nTHERE WERE ERRORS while setting up ColAHelL!\n Please check the log above. ColAHelL will not work!\n";
		White(false);	
	}
	cout << "====================================================================--\n\n\n";

	srand(unsigned int(time(NULL)));
	int r = rand() % 20;
	if(r==7)
		DisplayConfucius();
	else
		DisplayFortune();

	// restore cout
	if(CrLogFile) {
		cout.rdbuf(backup);
		filestr.close();
	}

	return this->Status;
}

int ColAHelL::ReplaceCHCfg(const char* FName)
{

	this->Status = 1;

	White(true);
	cout << "\nReplacing ColAHelL configuration...\n";
	cout << "---------------------------------------------------\n";
	
	White(false);
	cout << "Deleting ColAHelL objects...\n";
	// delete all CH objects...
	if(this->Spect)
	{ delete this->Spect; this->Spect=0; } 
	if(this->Proc)
	{ delete this->Proc; this->Proc=0; }
	if(this->CTof)
	{ delete this->CTof; this->CTof=0; }

	for(int i=0;i<num_presorters;i++) {
		if(this->presorters[i])
			{ delete this->presorters[i]; this->presorters[i]=0; }
	}

	cout << "Opening new configuration file...\n";
	if(!OpenXML((char*)FName)) {
		
		this->Status = -99;
		return this->Status;
	}
	
	cout << "\n\n\n--====================================================================";
	if(this->Status>0) {
		cout << "\nReplacing ColAHelL configuration: ";
		Green(true);
		gotoX(SUCCESSX); 
		cout << " ...success!\n";
		White(false);	
	}	
	if(this->Status<0 && this->Status>-90) {
		Yellow(true);
		cout << "\nTHERE WERE WARNINGS while setting up ColAHelL!\n Please check the log above. Some features might now work!\n";
		White(false);	
	}
	if(this->Status<-89) {
		Red(true);
		cout << "\nTHERE WERE ERRORS while setting up ColAHelL!\n Please check the log above. ColAHelL will not work!\n";
		White(false);	
	}
	cout << "====================================================================--\n\n\n";

	return this->Status;
}

void ColAHelL::DisplayFortune() 
{
    ifstream file("ColAHelL/Src/fcd.dat");
    string Fortune; 
	srand(unsigned int(time(NULL)));
	int r = rand() % 899 + 1;
	for(int i=0;i<r;i++)
		getline(file, Fortune);

	cout << "ColAHelL has some advice for you, my friend:\n";
	White(true);
	cout << "'";
	// Take care of word wrapping...!
	int cwidth = getwidth()-3;
	if(Fortune.length() > cwidth) {
		int ln = 1;
		while (Fortune.length()>0) {
			if(Fortune.length() > cwidth) {
				int cut = Fortune.find_last_of(" ",cwidth*ln);
				cout << Fortune.substr(0,cut) << endl;
				Fortune = Fortune.substr(cut+1);
			} else {
				cout << Fortune << "'" << endl << endl;
				break;
			}
		}
	} else {
		White(true);
		cout << Fortune << "'" << endl << endl;
	}
	White(false);
}

void ColAHelL::DisplayConfucius() 
{
    ifstream file("ColAHelL/Src/cfw.dat");
    string Fortune; 
	srand(unsigned int(time(NULL)));
	int r = rand() % 296 + 1;
	for(int i=0;i<r;i++)
		getline(file, Fortune);

	cout << "It is your lucky day!\nColAHelL's master, the great ";
	Blue(true);
	cout << "Confucius";
	White(false);
	cout << " is willing to teach you, my friend:\n";
	White(true);
	cout << "'";
	// Take care of word wrapping...!
	int cwidth = getwidth()-3;
	if(Fortune.length() > cwidth) {
		int ln = 1;
		while (Fortune.length()>0) {
			if(Fortune.length() > cwidth) {
				int cut = Fortune.find_last_of(" ",cwidth*ln);
				cout << Fortune.substr(0,cut) << endl;
				Fortune = Fortune.substr(cut+1);
			} else {
				cout << Fortune << "'" << endl << endl;
				break;
			}
		}
	} else {
		White(true);
		cout << Fortune << "'" << endl << endl;
	}
	White(false);
}

void ColAHelL::Error(string s) {
    ifstream file("ColAHelL/Src/demons.dat");
    string Er; 
	srand(unsigned int(time(NULL)));
	int r = rand() % 222 + 1;
	for(int i=0;i<r;i++)
		getline(file, Er);

	White(true);

	cout << endl << "     _.--\"\"--._	"<< endl;
	cout << "    /  _    _  \\	"<< endl;
	cout << " _  ( (_\\  /_) )  _"<< endl;
	cout << "{ \\._\\   /\\   /_./ }"<< endl;
	cout << "/_\"=-.}______{.-=\"_\\"<< endl;
	cout << " _  _.=(\"\"\"\")=._  _ "<< endl;
	cout << "(_'\"_.-\"`~~`\"-._\"'_)"<< endl;
	cout << " {_\"            \"_} " << endl;
	White(true);
	cout << "\nYour Demons SCREAM:\n";
	Red(true);
	cout << "'" << Er << "'" << endl << endl;
	cout << s;
	White(false);
	cout << "Press any key to continue!";
	_getch(); // wait for keypress
}

bool ColAHelL::OpenXML(char* FName="ColAHelL.cfg")
{
	ifstream theFile(FName);
	if(!theFile.is_open()) {
		Error("Error: Cannot open ColAHelL config file!\n");
		return false;
	}

	cout << "Opening XML file.\n";
	
	cout << "Accessing XML file...\n";

	// Read the xml file into a vector	
	vector<char> buffer((istreambuf_iterator<char>(theFile)), istreambuf_iterator<char>());
	buffer.push_back('\0');
	// Parse the buffer using the xml file parsing library into doc 
	try {
		doc.parse<0>(&buffer[0]);
	}
    catch (const std::runtime_error& e)
    {
        Red(true);
		cout << "\nA runtime error occurred: '" << e.what() << "'" << std::endl << std::endl;
		White(false);
		return false;
    }
    catch (const rapidxml::parse_error& e)
    {
		Red(true);
		cout << "\nA parse error occurred, massage is: ";
		White(true);
		cout << e.what() << endl;
		Red(true);
		cout << "Error is located roughly here:\n"; 
		White(false);		
		cout << e.where<char>() << std::endl << std::endl;
		White(false);
		return false;
    }
    catch (const std::exception& e)
    {
        Red(true);
		cout << "\nAn error occurred: '" << e.what() << "'" << std::endl << std::endl;
		White(false);
		return false;
    }
    catch (...)
    {
        Red(true);
		cout << "\nAn unknown error occurred..." << std::endl << std::endl;
		White(false);
		return false;
    }

	xml_node<> * root_node = doc.first_node("ColAHelL_cfg");
	if(root_node) {
		if(root_node->first_attribute("Version")) {
			double Ver = atof(root_node->first_attribute("Version")->value());
			cout << "Config file version: ";
			Green(true);
			cout << Ver;
			White(false);
		}
	} else {
		Error("'ColAHelL_cfg' node is missing in config-file!\n");
		return false;	
	}


	Green(true);
	gotoX(SUCCESSX);
	cout << " ...success!\n";
	White(false);
	cout << "====================================================================--\n\n";	
	
	if(!SetupTOFCalc()) {
		this->Status = -98;
	}

	if(!SetupSpectrometer()) {
		this->Status = -80;
	}

	if(!SetupDetectors()) {
		this->Status = -20;
	}

	if(!SetupMomentumParameters()) {
		this->Status = -15;
	}
	
	if(!SetupPresorters()) {
		this->Status = -50;
	}	
		
	if(!SetupReactions()) {
		this->Status = -40;
	}
	return true;
}

bool ColAHelL::SetupTOFCalc() 
{
	// Find our root node
	xml_node<> * root_node = doc.first_node("ColAHelL_cfg");
	
	cout << "\nSetting up time-of-flight calculator:\n";
	cout << "=====================================\n";

	// Look for TOF-Calc-node
	xml_node<> *node = root_node->first_node("TOF_calc");
	if(node) {
		if(node->first_attribute("Type")) {

			// calc TOF from BM
			if(ToUp(node->first_attribute("Type")->value()) == "BM") 
			{
				cout << "Times-of-flight will be calculated using a bunchmarker.\n";

				if(node->first_attribute("TDC_channel")) {
					int Ch = atoi(node->first_attribute("TDC_channel")->value());
					if(Ch<0 || Ch>64) {
						Error("Error: No valid TDC-channel for TOF calculations!\n");
						return false;
					}
					if(node->first_attribute("BM_period")) {
						// create a tof_calc with BM period
						CTof = new tof_calc_class(0,Ch,atof(node->first_attribute("BM_period")->value()));
						cout << "BM-channel is channel " << this->CTof->GetBMchannel() <<", BM-period is " << this->CTof->GetBMspacing() << "ns.\n";
					}  else {
						// create a tof_calc without BM period
						CTof = new tof_calc_class(0,Ch);
						cout << "BM-channel is channel " << this->CTof->GetBMchannel() <<".\n";	
					}
				
					if(node->first_attribute("Shift")) {
						double Sh = atof(node->first_attribute("Shift")->value());
						cout << "All TOFs are shifted by " << Sh <<" ns prior to modulo operation.\n";			
						CTof->set_mod_shift_e(Sh);
					}

					cout << "\nSetup of TOF calculator:"; 
					Green(true);
					gotoX(SUCCESSX);
					cout << " ...success!\n";
					White(false);
					cout << "====================================================================--\n\n";
					return true;
				}
			}

			// calc TOF from projectile
			if(ToUp(node->first_attribute("Type")->value()) == "PROJ") 
			{
				cout << "Times-of-flight will be calculated using the projectile detector.\n";

				CTof = new tof_calc_class(1);

				if(node->first_attribute("Shift")) {
					double Sh = atof(node->first_attribute("Shift")->value());
					cout << "All TOFs are shifted by " << Sh <<" ns.\n";			
					CTof->set_mod_shift_p(Sh);
				}

				if(node->first_attribute("TDC_Channel")) {
					Yellow(true);
					int Ch = atoi(node->first_attribute("TDC_Channel")->value());
					cout << "Warning: A bunchmarker channel (ch "<< Ch <<") was specified. Ignoring this setting...\n";
					White(false);
					this->Status = 0;
				}

				if(node->first_attribute("BM_period")) {
					Yellow(true);
					cout << "Warning: A bunchmarker period ("<< atof(node->first_attribute("BM_period")->value()) <<" ns) was specified. Ignoring this setting...\n";
					White(false);
					this->Status = 0;
				}

				cout << "\nSetup of TOF calculator:";
				Green(true);
				gotoX(SUCCESSX);
				cout << " ...success!\n";
				White(false);
				cout << "====================================================================--\n\n";
				return true;
			}

			// calc TOF from electron
			if(ToUp(node->first_attribute("Type")->value()) == "ELEC") 
			{
				cout << "Times-of-flight will be calculated using the electron detector.\n";

				CTof = new tof_calc_class(2);

				if(node->first_attribute("Shift")) {
					double Sh = atof(node->first_attribute("Shift")->value());
					cout << "All TOFs are shifted by "<< Sh <<" ns.\n";			
					CTof->set_mod_shift_e(Sh);
				}

				if(node->first_attribute("TDC_Channel")) {
					Yellow(true);
					int Ch = atoi(node->first_attribute("TDC_Channel")->value());
					cout << "Warning: A bunchmarker channel (ch "<< Ch <<") was specified. Ignoring this setting...\n";
					White(false);
					this->Status = 0;
				}

				if(node->first_attribute("BM_period")) {
					Yellow(true);
					cout << "Warning: A bunchmarker period ("<< atof(node->first_attribute("BM_period")->value()) <<" ns) was specified. Ignoring this setting...\n";
					White(false);
					this->Status = 0;
				}

				cout << "\nSetup of TOF calculator:";
				Green(true);
				gotoX(SUCCESSX);
				cout << " ...success!\n";
				White(false);
				cout << "====================================================================--\n\n";
				return true;
			}

			// did not find any of the three supported types...
			Error("Error: Invalid TYPE for time-of-flight calculation specified.\n");
			return false;
		}
	}
	Error("Error: No valid definition for TOF calculations were found!\n");
	return false;
}

bool ColAHelL::SetupSpectrometer() 
{
	// Find our root node
	xml_node<> * root_node = doc.first_node("ColAHelL_cfg");

	cout << "\nSetting up spectrometer + experiment properties:\n";
	cout << "===================================================\n";

	this->Spect = new spectrometer_class();
	this->Proc = new processor_class();

	// Look for COLTRIMS-node
	xml_node<> *node = root_node->first_node("COLTRIMS");

	if(!node) {
		Yellow(true);
		cout << "Warning: Did not find any COLTRIMS definitions. Is this on purpose?\n";
		White(false);
		this->Status = 0;
		return false;
	}

	xml_node<> *enode = node->first_node("Electron_arm");
	xml_node<> *inode = node->first_node("Ion_arm");
	xml_node<> *Bnode = node->first_node("Bfield");
	xml_node<> *Jnode = node->first_node("Jet");


	if((enode==0) || (inode==0)) {
		Yellow(true);
		cout << "Warning: No valid spectrometer definition found. Is this on purpose?\n";
		White(false);
		this->Status = 0;
		return true;
	}

	bool elec_arm = false;
	bool ion_arm = false;

	double Len_mm[3];
	double EVpcm[3];

	// now set up electron side
	if(enode) {
		if(enode->first_attribute("Type")) {
			// linear approx.? Maybe in the future.. Future is now!
			if(ToUp(enode->first_attribute("Type")->value()) == "LINEAR") {
				// linear approx.
				cout << "Electron arm defined, using the linear approximation.\n";
				Yellow(true);
				cout << "(This is rather unusual.. Are you sure...? ...Sure, sure..?).\n";
				White(false);
				if(enode->first_attribute("E_Vpcm")) {
					double fieldE = atof(enode->first_attribute("E_Vpcm")->value());
					if(fieldE > 10000.0 || fieldE<0.0) {
						Yellow(true);
						cout << "Warning: Reading a huge or negative E-field... ";
						Red(true);
						cout << "That doesn't seem right.\n";
						White(false);
						this->Status = 0;
					}
					cout << "E-field at target region is " << fieldE <<" V/cm.\n";
					this->Spect->set_electron_arm_lin(fieldE);
					elec_arm = true;
				} else {
					Error("No E-field value specified. This won't work...\n");
					return false;
				}
				if(enode->first_attribute("Mean_e_TOF")) {
					double meanTOF = atof(enode->first_attribute("Mean_e_TOF")->value());
					cout << "Mean electron time of flight is " << meanTOF <<" ns.\n";
					this->Spect->set_mean_tof_e(meanTOF);
					elec_arm = true;
				} else {
					Error("No mean electron flight time specified. This won't work...\n");
					return false;
				}
			}
			// Regions definitions
			if(ToUp(enode->first_attribute("Type")->value()) == "REGIONS") {
				int i=0;
				for(xml_node<> * reg_node = enode->first_node("Region"); reg_node; reg_node = reg_node->next_sibling("Region")) {
					if(i>2) {
						Yellow(true);
						cout << "Warning: ColAHelL supports up to 3 E-regions only! Further regions are ignored.\n";		
						White(false);
						break;
					}
					if(reg_node->first_attribute("L_mm"))
						Len_mm[i] = atof(reg_node->first_attribute("L_mm")->value());
					else {	
						Error("Found a region definition without length value. This won't work...\n");
						return false;
					}
					if(reg_node->first_attribute("E_Vpcm"))
						EVpcm[i] = atof(reg_node->first_attribute("E_Vpcm")->value());
					else {	
						Error("Found a region definition without E-field value. This won't work...\n");
						return false;
					}
					i++;
				}
				if(i==0) {
					Error("No a regions specified. This won't work...\n");
					return false;				
				}

				int num_reg =  i;
			
				cout << "Electron arm consisting of " << num_reg << " regions:\n";
				for(int i=0;i<num_reg;i++) {
					cout << "	R"<<i+1<<": "<< Len_mm[i]<<" mm, "<< EVpcm[i] <<" V/cm\n";	
				}
			
				this->Spect->set_electron_arm(num_reg, Len_mm, EVpcm);
				elec_arm = true;			
			}	
		}
	}

	// now set up ion side
	if(inode) {
		if(inode->first_attribute("Type")) {
			if(ToUp(inode->first_attribute("Type")->value())=="LINEAR") {
				// linear approx.
				cout << "Ion arm defined, using the linear approximation.\n";
				if(inode->first_attribute("E_Vpcm")) {
					double fieldE = atof(inode->first_attribute("E_Vpcm")->value());
					if(fieldE > 10000.0 || fieldE<0.0) {
						Yellow(true);
						cout << "Warning: Reading a huge or negative E-field... ";
						Red(true);
						cout << "That doesn't seem right.\n";
						White(false);
						this->Status = 0;
					}
					cout << "E-field at target region is " << fieldE <<" V/cm.\n";
					this->Spect->set_ion_arm_lin(fieldE);
					ion_arm = true;
				} else {
					Error("No E-field value specified. This won't work...\n");
					return false;
				}
			}

			if(ToUp(inode->first_attribute("Type")->value())=="REGIONS") {
			// Regions definitions


				int i=0;
				for(xml_node<> * reg_node = inode->first_node("REGION"); reg_node; reg_node = reg_node->next_sibling("REGION")) {
					if(i>2) {
						Yellow(true);
						cout << "Warning: ColAHelL supports up to 3 E-regions only! Further regions are ignored.\n";		
						White(false);
						break;
					}
					if(reg_node->first_attribute("L_mm"))
						Len_mm[i] = atof(reg_node->first_attribute("L_mm")->value());
					else {	
						Error("Found a region definition without length value. This won't work...\n");
						return false;
					}
					if(reg_node->first_attribute("E_Vpcm"))
						EVpcm[i] = atof(reg_node->first_attribute("E_Vpcm")->value());
					else {	
						Error("Found a region definition without E-field value. This won't work...\n");
						return false;
					}
					i++;
				}
				if(i==0) {
					Error("No a regions specified. This won't work...\n");
					return false;				
				}

				int num_reg =  i;

				cout << "Ion arm consisting of "<< num_reg <<" regions:\n";
				for(int i=0;i<num_reg;i++) {
					cout << "	R" << i+1 << ": " << Len_mm[i] << " mm, " << EVpcm[i] << " V/cm\n";	
				}
			
				this->Spect->set_ion_arm(num_reg, Len_mm, EVpcm);
				ion_arm = true;
			}
		}
	}

	// B-field
	if(elec_arm) {
		if(Bnode) {
			bool B_clockwise = false;
			if(Bnode->first_attribute("DIRECTION")) {
				if(ToUp(Bnode->first_attribute("DIRECTION")->value()) == "CLOCKWISE") {
					B_clockwise = true;
					cout << "(gyration is clockwise) ";
				} else {
					cout << "(gyration is counter clockwise) ";		
				}
			}

			bool B_ns = false;
			if(Bnode->first_attribute("UNIT")) {
				if(ToUp(Bnode->first_attribute("UNIT")->value()) == "NS")
					B_ns = true;
			}
			double St;
			if(Bnode->first_attribute("STRENGTH")) {
				St = atof(Bnode->first_attribute("STRENGTH")->value());
				cout << "Setting B-field to ";
				Green(true);
				cout << St << " " << Bnode->first_attribute("UNIT")->value() <<".\n";		
				White(false);
				this->Spect->set_Bfield(St,B_clockwise,B_ns);	
			} else {
				Yellow(true);
				cout << "Warning, no B-field strength defined. Setting B field to 0 G.\n";
				White(false);
				this->Spect->set_Bfield(0.0,true,false);
			}
		} else {
			White(true);
			cout << "No B-Field defined, setting B to 0 G.\n";
			White(false);
			this->Spect->set_Bfield(0.0,true,false);
		} 
	} else {
		if(Bnode) {
			Error("B-field is specified, but no electron arm definiton found!\nThat does not seem right...\n");
			return false;
		}	
	}

	if(Jnode) {
		bool mmpns=false;
		double ang = 0.0;

		if(Jnode->first_attribute("Unit")) {
			if(ToUp(Jnode->first_attribute("Unit")->value())=="MMPNS")
				mmpns = true;
		}

		if(Jnode->first_attribute("Angle")) 
			ang = atof(Jnode->first_attribute("Angle")->value());

		if(Jnode->first_attribute("Velocity")) {
			double vel = atof(Jnode->first_attribute("Velocity")->value()); 
			// convert to m/s
			if(mmpns)
				vel = vel*1e+9/1000.0;
			if(fabs(vel) > 1e+6) {
				Error("Jet velocity is extremely high. That doesn't seem right...\n");
				return false;
			}
			this->Spect->set_VJet(vel,ang);	
			cout << "Jet velocity set to " << vel << " m/s, direction is " << ang << " deg.\n";
		} else {
			Error("No jet velocity defined. That does not seem right...\n");
			return false;
		}
	}

	if(elec_arm || ion_arm) {
		this->Proc->Spect = this->Spect;

		cout << "\nSetup of spectrometer:";
		Green(true);
		gotoX(SUCCESSX);
		cout << " ...success!\n";
		White(false);
		cout << "====================================================================--\n\n";
		return true;
	} else {
		Error("Error: Did not find valid ion or electron arm definition!\n");
		return false;
	}
}

bool ColAHelL::SetupDetectors() {
	// Find our root node
	xml_node<> * root_node = doc.first_node("ColAHelL_cfg");

	// electron det
	xml_node<> *node = root_node->first_node("Electron_det");
	if(node) {
		double size = 80.0;
		if(node->first_attribute("Size"))
			size = atof(node->first_attribute("Size")->value());

		cout << "\nSetting parameters for electron detector:\n";
		cout << "=========================================\n";
		cout << "Detector size is: " << size <<" mm.\n";
		this->HGrams->edet_size = size;
		CopyDetPara(node,this->Proc->e_fac);
	}

	// ion det
	node = root_node->first_node("Ion_det");
	if(node) {
		double size = 80.0;
		if(node->first_attribute("Size"))
			size = atof(node->first_attribute("Size")->value());

		cout << "\nSetting parameters for ion detector:\n";
		cout << "=========================================\n";
		cout << "Detector size is: " << size <<" mm.\n";
		this->HGrams->rdet_size = size;
		CopyDetPara(node,this->Proc->r_fac);
	}

	// projectile det
	node = root_node->first_node("Projectile_det");
	if(node) {
		double size = 80.0;
		if(node->first_attribute("Size"))
			size = atof(node->first_attribute("Size")->value());

		cout << "\nSetting parameters for projectile detector:\n";
		cout << "===========================================\n";
		cout << "Detector size is: " << size <<" mm.\n";
		this->HGrams->pdet_size = size;
		CopyDetPara(node,this->Proc->p_fac);
	}
	return true;
}

void ColAHelL::CopyDetPara(xml_node<> *node, cor_param *fac) {
	xml_node<> *Mnode = node->first_node("Mirror");
	xml_node<> *dTnode = node->first_node("dTOffset_mmpns");
	xml_node<> *Shnode = node->first_node("Shift");
	xml_node<> *Stnode = node->first_node("Stretch");

	if(Mnode) {
		if(Mnode->first_attribute("Axis")) {
			if(ToUp(Mnode->first_attribute("Axis")->value())=="X") {
				fac->mir_x = true;
				cout << "Mirroring detector by inverting x-axis.\n";
			}
			if(ToUp(Mnode->first_attribute("Axis")->value())=="Y") {
				fac->mir_y = true;
				cout << "Mirroring detector by inverting y-axis.\n";
			}
			if((ToUp(Mnode->first_attribute("Axis")->value())=="XY") || (ToUp(Mnode->first_attribute("Axis")->value())=="YX")) { 
				fac->mir_x = true;
				fac->mir_y = true;
				cout << "Mirroring detector by inverting x-axis and y-axis.\n";
			}
		}
	}
		
	if(dTnode) {
		if(Mnode->first_attribute("X"))
			fac->EBx = atof(Mnode->first_attribute("X")->value());
		if(Mnode->first_attribute("Y"))
			fac->EBy = atof(Mnode->first_attribute("Y")->value());

		cout << "Setting up E/B drift of x="<< fac->EBx <<" mm/ns and y="<< fac->EBy <<" mm/ns.\n";
	}

	if(node->first_attribute("ProcessingOrder")) {
		if(ToUp(node->first_attribute("ProcessingOrder")->value())=="RS") {
			fac->raw_order = false;
			cout << "Raw detector data will be rotated first, then shifted + stretched...\n";
		} else {
			fac->raw_order = true;
			cout << "Raw detector data will be shifted + stretched first, then rotated...\n";
		}
	} 

	if(node->first_attribute("Rotate")) {
		fac->rot_ang = atof(node->first_attribute("Rotate")->value());
	}
	cout << "Rotating detector by "<< fac->rot_ang << " deg.\n";

	if(Shnode) {
		if(Shnode->first_attribute("X"))
			fac->dx = atof(Shnode->first_attribute("X")->value());
		if(Shnode->first_attribute("Y"))
			fac->dy = atof(Shnode->first_attribute("Y")->value());
		if(Shnode->first_attribute("T"))
			fac->dt = atof(Shnode->first_attribute("T")->value());
		cout << "Shifting detector:\n x: "<< fac->dx <<" mm\n y: "<< fac->dy <<" mm\n t: "<< fac->dt <<" ns\n";
	}

	if(Stnode) {
		if(Stnode->first_attribute("X"))
			fac->x_stretch = atof(Stnode->first_attribute("X")->value());
		if(Stnode->first_attribute("Y"))
			fac->y_stretch = atof(Stnode->first_attribute("Y")->value());
		if(Stnode->first_attribute("Total"))
			fac->overall_stretch = atof(Stnode->first_attribute("Total")->value());
		cout << "Stretching detector by a factor of:\n x: "<< fac->x_stretch <<"\n y: "<< fac->y_stretch <<"\n overall: "<< fac->overall_stretch <<"\n";
	}


}

bool ColAHelL::SetupMomentumParameters() {
	// Find our root node
	xml_node<> * root_node = doc.first_node("ColAHelL_cfg");

	// electron momenta
	xml_node<> *node = root_node->first_node("Electron_momenta");
	if(node) {
		cout << "\nSetting additional parameters for electron momenta:\n";
		cout << "=====================================================\n";
		CopyMomPara(node,this->Proc->e_mom_fac);
	}

	// ion momenta
	node = root_node->first_node("Ion_momenta");
	if(node) {
		cout << "\nSetting additional parameters for ion momenta:\n";
		cout << "=====================================================\n";
		CopyMomPara(node,this->Proc->r_mom_fac);
	}

	// projectile momenta
	node = root_node->first_node("Projectile_momenta");
	if(node) {
		cout << "\nSetting additional parameters for projectile momenta:\n";
		cout << "=======================================================\n";
		CopyMomPara(node,this->Proc->p_mom_fac);
	}
	return true;
}

void ColAHelL::CopyMomPara(xml_node<> *node, mom_cor_param *fac) {

	xml_node<> *Mnode = node->first_node("Mirror");
	xml_node<> *Shnode = node->first_node("Shift");
	xml_node<> *Stnode = node->first_node("Stretch");

	if(Shnode) {
		if(Shnode->first_attribute("X"))
			fac->dx = atof(Shnode->first_attribute("X")->value());
		if(Shnode->first_attribute("Y"))
			fac->dy = atof(Shnode->first_attribute("Y")->value());
		if(Shnode->first_attribute("Z"))
			fac->dz = atof(Shnode->first_attribute("Z")->value());

		cout << "Shifting momenta:\n x: "<< fac->dx <<" a.u.\n y: "<< fac->dy <<" a.u.\n z: "<< fac->dz <<" a.u.\n";
	}

	if(node->first_attribute("Rotate")) {
		fac->rot_ang = atof(node->first_attribute("Rotate")->value());
		cout << "Rotating momenta (around spectrometer axis) by "<< fac->rot_ang <<" deg.\n";
	}

	if(Stnode) {
		if(Stnode->first_attribute("X"))
			fac->x_stretch = atof(Stnode->first_attribute("X")->value());
		if(Stnode->first_attribute("Y"))
			fac->y_stretch = atof(Stnode->first_attribute("Y")->value());
		if(Stnode->first_attribute("Z"))
			fac->z_stretch = atof(Stnode->first_attribute("Z")->value());
		if(Stnode->first_attribute("Total"))
			fac->overall_stretch = atof(Stnode->first_attribute("Total")->value());

		cout << "Stretching momenta by a factor of:\n x: "<<fac->x_stretch <<"\n y: "<< fac->y_stretch <<"\n z: "<< fac->z_stretch <<"\n overall: "<< fac->overall_stretch <<"\n";
	}

	if(Mnode) {
		if(Mnode->first_attribute("Axis")) {
			string s=ToUp(Mnode->first_attribute("Axis")->value());
			std::size_t found = s.find("Z");
			if(found!=std::string::npos) {
				fac->mir_z = true;
				cout << "Mirroring momenta at xy-plane (inverting z).\n";
			}
			found = s.find("Y");
			if(found!=std::string::npos) {
				fac->mir_y = true;
				cout << "Mirroring momenta at xz-plane (inverting y).\n";
			}
			found = s.find("X");
			if(found!=std::string::npos) {
				fac->mir_x = true;
				cout << "Mirroring momenta at yz-plane (inverting x).\n";
			}			
		}
	}
}


bool ColAHelL::SetupPresorters() 
{
	this->maxNumPrsHists = 0;
	CH_event_struct e;
	this->prs_evt_out.push_back(e);
	this->prs_evt_out.clear();

	// Find our root node
	xml_node<> * root_node = doc.first_node("ColAHelL_cfg");

	// Presorters
	xml_node<> *node = root_node->first_node("Presorter");
	if(node) {

		cout << "\nSetting up presorters:\n";
		cout << "======================\n";

		
		int ListOfChannels[256];
		int ElementsInListOfChannels = 0;
		//ListOfChannels.clear();
		this->num_presorters = 0;

		int NumDifferentPRS = 0;

		//loop over all presorters and create them..
		for(node = root_node->first_node("Presorter"); node; node = node->next_sibling("Presorter")) {
			
			// this is not necessary, but makes sorting of the sequence
			// the presorters are applied easier in the end:
			bool prs_e_tof = false;
			bool prs_e_pos = false;
			bool prs_e_hit = false;
			bool prs_r_tof = false;
			bool prs_r_pos = false;
			bool prs_r_hit = false;
			bool prs_p_tof = false;
			bool prs_p_pos = false;
			bool prs_p_hit = false;
			bool prs_pipico = false;
			bool prs_sumdiff = false;
			bool prs_polyatomic = false;

			bool prs_usr[8];
			for(int j=0;j<8;j++)
				prs_usr[j] = false;

			int ch = 0;
			if(node->first_attribute("Channel"))
				ch = atoi(node->first_attribute("Channel")->value());
			else {
				Yellow(true);
				cout << "Warning: you did not provide a channel number!" << endl;
				White(false);
			}

			White(true);
			cout << "\nCreating new presorter,";
			Green(true);
			cout << " channel "<< ch <<".\n";
			White(false);
			cout << "--==============================================\n";
			presorters[num_presorters] = new presorter_class(ch,this->PrsHist,this->Spect,this->CTof);
			num_presorters++;
			// check if we have a presorter with same channel number already
			bool foundChannel = false;
			for(int j=0;j<ElementsInListOfChannels;j++)
			{
				if(ch == ListOfChannels[j]) {
					foundChannel = true;
					break;
				}
			}
			
			// if it is a new channel assign new histo block to presorter, if 
			// not use old one (=do nothing as it has already been assigned...)

			if(!foundChannel) {
				this->presorters[num_presorters - 1]->SetHistoBlock(NumDifferentPRS++);
				ListOfChannels[ElementsInListOfChannels++] = ch;
			}

			xml_node<> *Enode = node->first_node("Elec");
			xml_node<> *Inode = node->first_node("Ion");
			xml_node<> *Pnode = node->first_node("Proj");
			xml_node<> *Pinode = node->first_node("PIPICO");
			xml_node<> *SDnode = node->first_node("SumDiff");
			xml_node<> *Polnode = node->first_node("POLYATOMIC");

			// electron presorters
			if(Enode) {
				xml_node<> *TRnode = Enode->first_node("TOF_Range");
				xml_node<> *TPnode = Enode->first_node("TOF_Peak");
				xml_node<> *POSnode = Enode->first_node("Pos");
				xml_node<> *Hnode = Enode->first_node("Hits");

				cout << "Electrons: \n";
				if(TRnode) {
					double tmin = -100.0;
					double tmax = 100.0;
					if(TRnode->first_attribute("Min"))
						tmin = atof(TRnode->first_attribute("Min")->value());
					if(TRnode->first_attribute("Max"))
						tmax = atof(TRnode->first_attribute("Max")->value());
					prs_e_tof = true;
					this->presorters[num_presorters - 1]->set_tof_e(tmin, tmax);
					cout << " Times-of-flight: " << tmin << " ns < TOF < "<< tmax <<" ns.\n";
				}

				if(TPnode) {
					double tmin = -100.0;
					double tmax = 100.0;
					if(TPnode->first_attribute("Width") && TPnode->first_attribute("Center")) {
						tmin = atof(TPnode->first_attribute("Center")->value()) - atof(TPnode->first_attribute("Width")->value());
						tmax = atof(TPnode->first_attribute("Center")->value()) + atof(TPnode->first_attribute("Width")->value());
						prs_e_tof = true;
						this->presorters[num_presorters - 1]->set_tof_e(tmin, tmax);
						cout << " Times-of-flight: " << tmin << " ns < TOF < "<< tmax <<" ns.\n";
					}
				}

				if(POSnode) {
					if(POSnode->first_attribute("Center_X") && POSnode->first_attribute("Center_Y") && POSnode->first_attribute("Radius")) {
						double posx = atof(POSnode->first_attribute("Center_X")->value());
						double posy = atof(POSnode->first_attribute("Center_Y")->value());
						double posrad = atof(POSnode->first_attribute("Radius")->value());	
						prs_e_pos = true;
						this->presorters[num_presorters - 1]->set_pos_e(posx, posy, posrad);
						cout << " Hits are valid around ("<< posx <<"/"<< posy <<") within "<< posrad <<" mm.\n";
					}
				}

				if(Hnode) {
					int hmin = 0;
					int hmax = 64;
					if(Hnode->first_attribute("Min"))
						hmin = atoi(Hnode->first_attribute("Min")->value());
					if(Hnode->first_attribute("Max"))
						hmax = atoi(Hnode->first_attribute("Max")->value());
					prs_e_hit = true;
					this->presorters[num_presorters - 1]->set_hits_e(hmin, hmax);
					cout << " Number of hits on electron detector is restricted to "<< hmin <<" <= hits <= "<< hmax <<".\n";		
				}	

				cout << "------------------------------\n";
			}					

			// ion presorters
			if(Inode) {
				xml_node<> *TRnode = Inode->first_node("TOF_Range");
				xml_node<> *TPnode = Inode->first_node("TOF_Peak");
				xml_node<> *POSnode = Inode->first_node("Pos");
				xml_node<> *Hnode = Inode->first_node("Hits");

				cout << "Ions: \n";
				if(TRnode) {
					double tmin = -100.0;
					double tmax = 100.0;
					if(TRnode->first_attribute("Min"))
						tmin = atof(TRnode->first_attribute("Min")->value());
					if(TRnode->first_attribute("Max"))
						tmax = atof(TRnode->first_attribute("Max")->value());
					prs_r_tof = true;
					this->presorters[num_presorters - 1]->set_tof_r(tmin, tmax);
					cout << " Times-of-flight: " << tmin << " ns < TOF < "<< tmax <<" ns.\n";
				}

				if(TPnode) {
					double tmin = -100.0;
					double tmax = 100.0;
					if(TPnode->first_attribute("Width") && TPnode->first_attribute("Center")) {
						tmin = atof(TPnode->first_attribute("Center")->value()) - atof(TPnode->first_attribute("Width")->value());
						tmax = atof(TPnode->first_attribute("Center")->value()) + atof(TPnode->first_attribute("Width")->value());
						prs_r_tof = true;
						this->presorters[num_presorters - 1]->set_tof_r(tmin, tmax);
						cout << " Times-of-flight: " << tmin << " ns < TOF < "<< tmax <<" ns.\n";
					}
				}

				if(POSnode) {
					if(POSnode->first_attribute("Center_X") && POSnode->first_attribute("Center_Y") && POSnode->first_attribute("Radius")) {
						double posx = atof(POSnode->first_attribute("Center_X")->value());
						double posy = atof(POSnode->first_attribute("Center_Y")->value());
						double posrad = atof(POSnode->first_attribute("Radius")->value());	
						prs_r_pos = true;
						this->presorters[num_presorters - 1]->set_pos_r(posx, posy, posrad);
						cout << " Hits are valid around ("<< posx <<"/"<< posy <<") within "<< posrad <<" mm.\n";
					}
				}

				if(Hnode) {
					int hmin = 0;
					int hmax = 64;
					if(Hnode->first_attribute("Min"))
						hmin = atoi(Hnode->first_attribute("Min")->value());
					if(Hnode->first_attribute("Max"))
						hmax = atoi(Hnode->first_attribute("Max")->value());
					prs_r_hit = true;
					this->presorters[num_presorters - 1]->set_hits_r(hmin, hmax);
					cout << " Number of hits on ion detector is restricted to "<< hmin <<" <= hits <= "<< hmax <<".\n";		
				}	
				cout << "------------------------------\n";
			}					

			// projectile presorters
			if(Pnode) {
				xml_node<> *TRnode = Pnode->first_node("TOF_Range");
				xml_node<> *TPnode = Pnode->first_node("TOF_Peak");
				xml_node<> *POSnode = Pnode->first_node("Pos");
				xml_node<> *Hnode = Pnode->first_node("Hits");

				cout << "Projectiles: \n";
				if(TRnode) {
					double tmin = -100.0;
					double tmax = 100.0;
					if(TRnode->first_attribute("Min"))
						tmin = atof(TRnode->first_attribute("Min")->value());
					if(TRnode->first_attribute("Max"))
						tmax = atof(TRnode->first_attribute("Max")->value());
					prs_p_tof = true;
					this->presorters[num_presorters - 1]->set_tof_p(tmin, tmax);
					cout << " Times-of-flight: " << tmin << " ns < TOF < "<< tmax <<" ns.\n";
				}

				if(TPnode) {
					double tmin = -100.0;
					double tmax = 100.0;
					if(TPnode->first_attribute("Width") && TPnode->first_attribute("Center")) {
						tmin = atof(TPnode->first_attribute("Center")->value()) - atof(TPnode->first_attribute("Width")->value());
						tmax = atof(TPnode->first_attribute("Center")->value()) + atof(TPnode->first_attribute("Width")->value());
						prs_p_tof = true;
						this->presorters[num_presorters - 1]->set_tof_p(tmin, tmax);
						cout << " Times-of-flight: " << tmin << " ns < TOF < "<< tmax <<" ns.\n";
					}
				}

				if(POSnode) {
					if(POSnode->first_attribute("Center_X") && POSnode->first_attribute("Center_Y") && POSnode->first_attribute("Radius")) {
						double posx = atof(POSnode->first_attribute("Center_X")->value());
						double posy = atof(POSnode->first_attribute("Center_Y")->value());
						double posrad = atof(POSnode->first_attribute("Radius")->value());	
						prs_p_pos = true;
						this->presorters[num_presorters - 1]->set_pos_p(posx, posy, posrad);
						cout << " Hits are valid around ("<< posx <<"/"<< posy <<") within "<< posrad <<" mm.\n";
					}
				}

				if(Hnode) {
					int hmin = 0;
					int hmax = 64;
					if(Hnode->first_attribute("Min"))
						hmin = atoi(Hnode->first_attribute("Min")->value());
					if(Hnode->first_attribute("Max"))
						hmax = atoi(Hnode->first_attribute("Max")->value());
					prs_p_hit = true;
					this->presorters[num_presorters - 1]->set_hits_p(hmin, hmax);
					cout << " Number of hits on projectile detector is restricted to "<< hmin <<" <= hits <= "<< hmax <<".\n";		
				}	

				cout << "------------------------------\n";
			}

			// sumdiff presorters
			if(SDnode) {
				int num_ions = 0;
				double CSum = 0.0;
				double WSum = 0.0;
				double CDiff = 0.0;
				double WDiff = 0.0;

				cout << "Summ/Diff presorter enabled for ";

				if(SDnode->first_attribute("Ions"))
					num_ions = atoi(SDnode->first_attribute("Ions")->value());
				else {

					Error("\nPlease specify how many ions to include!\n");
					return false;
				}

				if(SDnode->first_attribute("Center_Sum"))
					CSum = atof(SDnode->first_attribute("Center_Sum")->value());
				else {
					Error("\nPlease specify the center of the sum-TOF!\n");
					return false;
				}
				if(SDnode->first_attribute("Width_Sum"))
					WSum = atof(SDnode->first_attribute("Width_Sum")->value());
				else {
					Error("\nPlease specify the width of the sum-TOF!\n");
					return false;
				}

				if(SDnode->first_attribute("Center_Diff"))
					CDiff = atof(SDnode->first_attribute("Center_Diff")->value());
				else {
					Error("\nPlease specify the center of the difference-TOF!\n");
					return false;
				}
				if(SDnode->first_attribute("Width_Diff"))
					WDiff = atof(SDnode->first_attribute("Width_Diff")->value());
				else {
					Error("\nPlease specify the width of the difference-TOF!\n");
					return false;
				}

				prs_sumdiff = true;
				this->presorters[num_presorters - 1]->set_sumdiff(num_ions,CSum,WSum,CDiff,WDiff);
				
				Green(true);
				cout << num_ions << endl;
				White(false);
				cout << "Sum ion TOF is centered at ";
				Green(true);
				cout << CSum;
				White(false);
				cout << " ns, with a width of ";
				Green(true);
				cout << WSum;
				White(false);
				cout << " ns.\n";
				cout << "Ion TOF difference is centered at ";
				Green(true);
				cout << CDiff;
				White(false);
				cout << " ns, width a width of ";
				Green(true);
				cout << WDiff;
				White(false);
				cout << " ns.\n";

				cout << "------------------------------\n";
			}

			// PIPICO presorter
			if(Pinode) {
				xml_node<> *IDnode = Pinode->first_node("Ion_def");
				xml_node<> *IDnode2 = IDnode->next_sibling("Ion_def");
				if(IDnode && IDnode2) {
					double width = 200.0;
					double m1,m2 = 14.0;
					double q1,q2 = 1.0;
					if(IDnode->first_attribute("Mass")) 
						m1 = atof(IDnode->first_attribute("Mass")->value());
					if(IDnode2->first_attribute("Mass")) 
						m2 = atof(IDnode2->first_attribute("Mass")->value());
					if(IDnode->first_attribute("Charge")) 
						q1 = atof(IDnode->first_attribute("Charge")->value());
					if(IDnode->first_attribute("Charge")) 
						q2 = atof(IDnode2->first_attribute("Charge")->value());
					if(Pinode->first_attribute("Width"))
						width = atof(Pinode->first_attribute("Width")->value());

					cout << "PIPICO presorter enabled, width is: ";
					White(true);
					cout << width <<" ns.\n";
					Green(true);
					cout << "m1/q1 = " << m1 <<"/"<< q1 <<", m2/q2 = " << m2 <<"/"<< q2; 
					White(false);
					this->presorters[num_presorters - 1]->set_pipico(m1,q1,m2,q2,width);
					prs_pipico = true;
					cout << "\n------------------------------\n";
				} else {
					Error("Need at least two ion definitions for a PIPICO presorter!\n");
					return false;
				}
			}

			// POLYATOMIC presorter
			if(Polnode) {
				int num_ions = -1;
				if(Polnode->first_attribute("Ions"))
					num_ions = atoi(Polnode->first_attribute("Ions")->value());
				else {
					Error("Please define how many ions are expected ('Ions = ..') ...\n");					
					return false;
				}
				if(num_ions>16) {
					Error("Sorry, only up to 16 ions are supported...\n");				
					return false;
				}

				cout << "POLYATOMIC presorter enabled, looking for ";
				Green(true);
				cout << num_ions; 
				White(false);
				cout << " ions.\n";
				this->presorters[num_presorters - 1]->set_polyatomic(num_ions);
				prs_polyatomic = true;

				int ions_defined = 0;
				for(xml_node<> * ID_node = Polnode->first_node("Ion_def"); ID_node; ID_node = ID_node->next_sibling("Ion_def")) {
					cout << "Ion definitions:\n";

					cout << "	Ion ";
					White(true);
					cout << ions_defined + 1;
					White(false);
					cout << ": ";

					double tmin = -1.0;
					double tmax = -1.0;
					if(ID_node->first_attribute("TOF_min")) 
						tmin = atof(ID_node->first_attribute("TOF_min")->value());
					if(ID_node->first_attribute("TOF_max")) 
						tmax = atof(ID_node->first_attribute("TOF_max")->value());

					if(tmin < 0.0 || tmax < 0.0) {
						Error("You need to set pairs of TOF-min and TOF-max, not just one of the two!\n");				
						return false;
					}
					cout << "TOF from "; 
					Green(true);
					cout << tmin; 
					White(false);
					cout << " ns to ";
					Green(true);
					cout << tmax;
					White(false);
					cout << " ns.\n";
					this->presorters[num_presorters - 1]->polyatomic_add_ion(tmin,tmax);
					ions_defined++;
					if(ions_defined>num_ions) 
						break;
				}
				
				if(ions_defined<num_ions) {
					Error("You defined a few ions, but not enough!\n");
					return false;
				}
				cout << "------------------------------\n";
			}

			// user defined presorters
			int un = 0;
			for(xml_node<> * Unode = node->first_node("User_defined"); Unode; Unode = Unode->next_sibling("User_defined")) {

				cout << "User defined Presorter: \n";
				if(un>8) {
					Error("Too many user presorters defined! (more than 8)..");
					return false;
				}
				int uprs_num = 0;
				if(Unode->first_attribute("Num"))
					uprs_num = atoi(Unode->first_attribute("Num")->value());
					cout << " Presorter function: ";
					Yellow(true);
					cout << uprs_num << endl;
					White(false);

					prs_usr[un] = true;

					// get all user parameters

					int pn = 0;
					for(xml_node<> * Pparnode = Unode->first_node("Prs_par"); Pparnode; Pparnode = Pparnode->next_sibling("Prs_par")) {
						if(pn>32) {
							Error("Too many user parameters (more than 32)..");
							return false;
						}
						double user_par = 0.0;
						if(Pparnode->first_attribute("Val"))
							user_par = atof(Pparnode->first_attribute("Val")->value());
						cout << " - User parameter " << pn << " set to: ";
						Yellow(true);
						cout << user_par << endl;
						White(false);
						this->presorters[num_presorters - 1]->add_user_parameter(uprs_num, user_par);
						pn++;
					}
				
		
				cout << "------------------------------\n";
				un++;
			}

			// set up presorter function pointers
			// The sequence listed here is the same sequence the presorters are
			// run (if more than one presorter is defined for a given channel)!

			// 1st: all presorters, that require some kind of cioncidence condition 
			if(prs_pipico) 
				this->presorters[num_presorters - 1]->append_presorter_to_list(&CH::presorter_class::pipico);
			if(prs_sumdiff) 
				this->presorters[num_presorters - 1]->append_presorter_to_list(&CH::presorter_class::sumdiff);
			if(prs_polyatomic)
				this->presorters[num_presorters - 1]->append_presorter_to_list(&CH::presorter_class::polyatomic);

			if(prs_usr[0]) 
				this->presorters[num_presorters - 1]->append_presorter_to_list(&CH::presorter_class::user_prs_0);
			if(prs_usr[1])
				this->presorters[num_presorters - 1]->append_presorter_to_list(&CH::presorter_class::user_prs_1);
			if(prs_usr[2]) 
				this->presorters[num_presorters - 1]->append_presorter_to_list(&CH::presorter_class::user_prs_2);
			if(prs_usr[3]) 
				this->presorters[num_presorters - 1]->append_presorter_to_list(&CH::presorter_class::user_prs_3);
			if(prs_usr[4]) 
				this->presorters[num_presorters - 1]->append_presorter_to_list(&CH::presorter_class::user_prs_4);
			if(prs_usr[5]) 
				this->presorters[num_presorters - 1]->append_presorter_to_list(&CH::presorter_class::user_prs_5);
			if(prs_usr[6]) 
				this->presorters[num_presorters - 1]->append_presorter_to_list(&CH::presorter_class::user_prs_6);
			if(prs_usr[7]) 
				this->presorters[num_presorters - 1]->append_presorter_to_list(&CH::presorter_class::user_prs_7);

			// now: disable single particles that are still left...
			if(prs_e_tof) 
				this->presorters[num_presorters - 1]->append_presorter_to_list(&CH::presorter_class::tof_e);
			if(prs_r_tof) 
				this->presorters[num_presorters - 1]->append_presorter_to_list(&CH::presorter_class::tof_r);
			if(prs_p_tof) 
				this->presorters[num_presorters - 1]->append_presorter_to_list(&CH::presorter_class::tof_p);
			if(prs_e_pos) 
				this->presorters[num_presorters - 1]->append_presorter_to_list(&CH::presorter_class::pos_e);
			if(prs_r_pos) 
				this->presorters[num_presorters - 1]->append_presorter_to_list(&CH::presorter_class::pos_r);
			if(prs_p_pos) 
				this->presorters[num_presorters - 1]->append_presorter_to_list(&CH::presorter_class::pos_p);
			if(prs_e_hit) 
				this->presorters[num_presorters - 1]->append_presorter_to_list(&CH::presorter_class::hits_e);
			if(prs_r_hit) 
				this->presorters[num_presorters - 1]->append_presorter_to_list(&CH::presorter_class::hits_r);
			if(prs_p_hit) 
				this->presorters[num_presorters - 1]->append_presorter_to_list(&CH::presorter_class::hits_p);	
		}

		this->maxNumPrsHists = NumDifferentPRS * NHISTS;

		cout << "\nSetup of presorters:";
		Green(true);
		gotoX(SUCCESSX);
		cout << " ...success!\n";
		White(false);
		cout << "====================================================================--\n\n";
	}

	return true;
}

void ColAHelL::prepare_tag(int type, rtag_struct *tag, xml_node<> *node, bool negate) {
	std::string tags_ascii[27] = {"X","Y","TOF","PX","PY","PZ","PXY","PXZ","PYZ","P","PHI_DEG","PHIPOS","THETA_DEG","ENERGY",
								  "PSUMX","PSUMY","PSUMZ","PSUM","PSUM_PHI_DEG","PSUM_THETA_DEG",
								  "PRELX","PRELY","PRELZ","PREL","PREL_PHI_DEG","PREL_THETA_DEG","KER"};
	tag->type = type;
	if(node->first_attribute("From"))
		tag->min = atof(node->first_attribute("From")->value());
	else {
		Error("Error: Value 'From' is missing! Sorry, this won't work.\n");
		this->Status = -1;
		return;
	}
	if(node->first_attribute("To"))
		tag->max = atof(node->first_attribute("To")->value());
	else {
		Error("Error: Value 'To' is missing! Sorry, this won't work.\n");
		this->Status = -1;
		return;
	}

	tag->negate = negate;

	tag->particle_num = -1;
	if(node->first_attribute("Particle")) {
		tag->particle_num = (short)atof(node->first_attribute("Particle")->value());
		cout << "   Applying the following condition to particle number ";
		Green(true);
		cout << tag->particle_num;
		White(false);
		cout <<":\n";		
	}

	if(node->first_attribute("Var")) {
		std::string var = node->first_attribute("Var")->value();
		for(int i=0;i<NUMBER_OF_TAGS;i++) {
			if(var == tags_ascii[i]) {
				// single particle tag in diatomic or polyatomic (i.e. not electron)..
				if((type >2) && (i<FIRST_DIATOMIC_TAG)) {
					tag->type = RTYPE_ION;
				}
				// diatomic tag for single particle.. Won't work..
				if(type == RTYPE_ELEC || type == RTYPE_ION) {
					if(i >= FIRST_DIATOMIC_TAG) { 
						Error("Error: A molecule tag was used to invalidate a single particle. Sorry, this won't work.\n");
						this->Status = -1;
						return;
					}
				}
				tag->tag_type = i;

				if(negate) {
					cout << "   Valid if '";
					White(true);	
					cout << var.c_str();
					White(false);
					cout <<"' is within ";
					Green(true);
					cout << tag->min <<" to "<< tag->max;
					White(false);
					cout <<".\n";
				} else {
					cout << "   Not valid if '";
					White(true);	
					cout << var.c_str();
					White(false);
					cout <<"' is from ";
					Red(true);
					cout << tag->min <<" to "<< tag->max;
					White(false);
					cout <<".\n";
				}
				return;
			}
		}
		string e = "Error: Invalid tag definition '" + var + "'...! This won't work.\n";
		Error(e);
		this->Status = -1;
	} else {
		Error("Error: No tag provided...! This won't work.\n");
		this->Status = -1;
	}
}

bool ColAHelL::SetupReactions() 
{
	this->maxNumCHHists = 0;
	// Find our root node
	xml_node<> * root_node = doc.first_node("ColAHelL_cfg");

	//loop over all reaction definitions
	int i=0;
	for(xml_node<> * node = root_node->first_node("Reaction"); node; node = node->next_sibling("Reaction")) {
	
		if(i==0) {
			cout << "\nSetting up reaction definitions:\n";
			cout << "===================================\n";
		}

		if(this->Proc->number_of_reactions_defined>MAX_NUM_REACTIONS-1) {
			Error("Error: You defined too many reactions.\n");
			this->Status = -1;
			return false;									
		}

		// create new reaction
		this->Proc->new_reaction();
		int rnum = this->Proc->number_of_reactions_defined - 1;
		this->Proc->reaction_list[rnum]->type = -1;

		// apply reaction definition
		int ch = 0;
		if(node->first_attribute("Channel"))
			ch = atoi(node->first_attribute("Channel")->value());
		else {
			Yellow(true);
			cout << "Warning you did not specify a presorter channel to connect to!" << endl;
			White(false);
		}

		this->Proc->reaction_list[rnum]->channel = ch;
		
		if(node->first_attribute("Name"))
			sprintf(this->Proc->reaction_list[rnum]->name, node->first_attribute("Name")->value());
		else
			sprintf(this->Proc->reaction_list[rnum]->name, "Dummy_" + ch);

		if(node->first_attribute("ID"))
			this->Proc->reaction_list[rnum]->ID = atoi(node->first_attribute("ID")->value());
		else
			this->Proc->reaction_list[rnum]->ID = ch;

		cout << "\n\nSetting up reaction '"<< this->Proc->reaction_list[rnum]->name <<"' connected to presorter channel "<< this->Proc->reaction_list[rnum]->channel;
		cout << "\nReaction ID is " << this->Proc->reaction_list[rnum]->ID << ":\n";
		cout << "===============================================================---\n|\n";

		// electron det
		xml_node<> * Tnode = node->first_node("Electron_det");
		if(Tnode) {
			cout << "\nSetting additional parameters for electron detector:\n";
			cout << "=====================================================\n";
			CopyDetPara(Tnode,this->Proc->reaction_list[rnum]->e_fac);
		}

		// ion det
		Tnode = node->first_node("Ion_det");
		if(Tnode) {
			cout << "\nSetting additional parameters for ion detector:\n";
			cout << "=====================================================\n";
			CopyDetPara(Tnode,this->Proc->reaction_list[rnum]->r_fac);
		}

		// projectile det
		Tnode = node->first_node("Projectile_det");
		if(Tnode) {
			cout << "\nSetting additional parameters for projectile detector:\n";
			cout << "=======================================================\n";
			CopyDetPara(Tnode,this->Proc->reaction_list[rnum]->p_fac);
		}

		// electron momenta
		Tnode = node->first_node("Electron_momenta");
		if(Tnode) {
			cout << "\nSetting additional parameters for electron momenta:\n";
			cout << "=====================================================\n";
			CopyMomPara(Tnode,this->Proc->reaction_list[rnum]->e_mom_fac);
		}

		// ion momenta
		Tnode = node->first_node("Ion_momenta");
		if(Tnode) {
			cout << "\nSetting additional parameters for ion momenta:\n";
			cout << "=====================================================\n";
			CopyMomPara(Tnode,this->Proc->reaction_list[rnum]->r_mom_fac);
		}

		// projectile momenta
		Tnode = node->first_node("Projectile_momenta");
		if(Tnode) {
			cout << "\nSetting additional parameters for projectile momenta:\n";
			cout << "=======================================================\n";
			CopyMomPara(Tnode,this->Proc->reaction_list[rnum]->p_mom_fac);
		}

		Tnode = node->first_node("SINGLEION");
		if(Tnode) {
			this->Proc->reaction_list[rnum]->type = RTYPE_ION;
			cout << "\nType is: ";
			Green(true);
			cout << "single ion\n";
			White(false);

			xml_node<> *IDnode = Tnode->first_node("Ion_def");
			if(IDnode) {
				double m,q,mean_tof = 0.0;
				if(IDnode->first_attribute("Mass"))
					m = atof(IDnode->first_attribute("Mass")->value());
				if(IDnode->first_attribute("Charge"))
					q = atof(IDnode->first_attribute("Charge")->value());

				this->Proc->reaction_list[rnum]->mass.push_back(m);
				this->Proc->reaction_list[rnum]->charge.push_back(q);
				cout << "Ion defined as m = "<< (int)this->Proc->reaction_list[rnum]->mass[0] <<" amu, q = "<< (int)this->Proc->reaction_list[rnum]->charge[0] <<" au.\n";

				if(IDnode->first_attribute("Mean_tof")) {
					mean_tof = atof(IDnode->first_attribute("Mean_tof")->value());
					this->Proc->reaction_list[rnum]->t_mean.push_back(mean_tof);
					cout << "Mean TOF is "<< this->Proc->reaction_list[rnum]->t_mean[0] <<" ns.\n";
				} else
					this->Proc->reaction_list[rnum]->t_mean.push_back(mean_tof);

			} else {
				Error("Error: no ion definitions found...! This won't work.\n");
				this->Status = -1;
				return false;
			}

			cout << "\nApplying conditions: \n";
			int added_tags=0;
			//iterate over all "if"
			for(xml_node<> *Invnode = Tnode->first_node("Invalidate_if"); Invnode; Invnode = Invnode->next_sibling("Invalidate_if")) {
				this->Proc->new_rtag(rnum);
				prepare_tag(RTYPE_ION, this->Proc->reaction_list[rnum]->tag[this->Proc->reaction_list[rnum]->num_tags-1], Invnode, false);	
				added_tags++;
			}
			//iterate over all "if_not"
			for(xml_node<> *Invnode = Tnode->first_node("Invalidate_ifnot"); Invnode; Invnode = Invnode->next_sibling("Invalidate_ifnot")) {
				this->Proc->new_rtag(rnum);
				prepare_tag(RTYPE_ION, this->Proc->reaction_list[rnum]->tag[this->Proc->reaction_list[rnum]->num_tags-1], Invnode, true);
				added_tags++;
			}
			if(added_tags == 0) {
				White(true);
				cout << "   none.\n";
				White(false);
			}
		}
		
		Tnode = node->first_node("DIATOMIC");
		if(Tnode) {
			this->Proc->reaction_list[rnum]->type = RTYPE_DIATOMIC;
			cout << "\nType is: ";
			Green(true);
			cout << "diatomic molecule\n";
			White(false);
			
			this->Proc->reaction_list[rnum]->diatomic_prel = RELMOMSTD;
			if(Tnode->first_attribute("RELMOM")) {
				if(ToUp(Tnode->first_attribute("RELMOM")->value())=="NO_CM_P") {
					this->Proc->reaction_list[rnum]->diatomic_prel = RELMOMNOP;					
					cout << "Will calculate relative momenta assuming";
					Green(true);
					cout <<	" sum momentum";
					White(false);
					cout <<  " is zero.\n";
				}
				if(ToUp(Tnode->first_attribute("RELMOM")->value())=="NO_R") {
					this->Proc->reaction_list[rnum]->diatomic_prel = RELMOMNOR;
					Green(true);
					cout << "Will calculate relative momenta assuming";
					Green(true);
					cout <<	" target size";
					White(false);
					cout <<  " is zero.\n";
					White(false);
				}
			}

			int ions_defined = 0;
			//loop over all ion definitions..
			for(xml_node<> *IDnode = Tnode->first_node("Ion_def"); IDnode; IDnode = IDnode->next_sibling("Ion_def")) {

				double m = 0.0;
				double q = -1.0;
				double mean_tof = 0.0;
				
				// retrieve all parameters from cfg
				if(IDnode->first_attribute("Charge"))
					q = atof(IDnode->first_attribute("Charge")->value());
				if(IDnode->first_attribute("Mass"))
					m = atof(IDnode->first_attribute("Mass")->value());
				if(IDnode->first_attribute("Mean_tof")) 
					mean_tof = atof(IDnode->first_attribute("Mean_tof")->value());

				this->Proc->reaction_list[rnum]->mass.push_back(m);
				this->Proc->reaction_list[rnum]->charge.push_back(q);
				this->Proc->reaction_list[rnum]->t_mean.push_back(mean_tof);

				if(q < 0.5) {
					White(true);
					cout << "Found neutral particle in reaction defintion, assuming molecular dissociation:\n";
					White(false);
					this->Proc->reaction_list[rnum]->type = RTYPE_DIATOMIC_DISS;
				}

				if(q < -0.5) {
					Error("Error: Ion definition without specifying a charge. This won't work.\n");
					this->Status = -1;
					return false;					
				}
				if(m < 0.5) {
					Error("Error: Ion definition without specifying a mass. This won't work.\n");
					this->Status = -1;
					return false;					
				}
				if(mean_tof < 0.5 && this->Spect->ion_side->linear_approximation) {
					Error("Error: Ion definition without specifying mean time-of-flight. This won't work.\n");
					this->Status = -1;
					return false;					
				}
				if(m > 0.5 && !this->Spect->ion_side->linear_approximation) {
					cout << "Ion ";
					Green(true);
					cout << ions_defined + 1;
					White(false);
					cout << " defined as ";
					White(true);
					cout << "m = "<< (int)this->Proc->reaction_list[rnum]->mass[ions_defined] <<" amu, q = "<< (int)this->Proc->reaction_list[rnum]->charge[ions_defined] << " au";
					White(false);
					cout << ".\n";
				}
				if(mean_tof > 0.5 && this->Spect->ion_side->linear_approximation) {
					cout << "Ion ";
					Green(true);
					cout << ions_defined + 1;
					White(false);
					cout << " defined as ";
					White(true);
					cout << "mean tof = "<< (int)this->Proc->reaction_list[rnum]->t_mean[ions_defined] <<" ns, m = "<< (int)this->Proc->reaction_list[rnum]->mass[ions_defined] <<" amu, q = "<< (int)this->Proc->reaction_list[rnum]->charge[ions_defined] << " au";
					White(false);
					cout << ".\n";
				}
				ions_defined++;
			}
			if(ions_defined != 2)  {
				Error("Error: Wrong number of ion definitions found...! (You need two, eh?) This won't work.\n");
				this->Status = -1;
				return false;
			}

			if(Tnode->first_node("Randomize")) {
				this->Proc->reaction_list[rnum]->randomize_ions = true;
				cout << "Ions will be ";
				Green(true);
				cout <<"randomly swapped";
				White(false);
				cout << ".. (Is this molecule homonuclear?)\n";
			}

			cout << "\nApplying conditions: \n";
			int added_tags = 0;
			//iterate over all "if"
			for(xml_node<> *Invnode = Tnode->first_node("Invalidate_if"); Invnode; Invnode = Invnode->next_sibling("Invalidate_if")) {
				this->Proc->new_rtag(rnum);
				prepare_tag(RTYPE_DIATOMIC, this->Proc->reaction_list[rnum]->tag[this->Proc->reaction_list[rnum]->num_tags-1], Invnode, false);	
				added_tags++;
			}
			//iterate over all "if_not"
			for(xml_node<> *Invnode = Tnode->first_node("Invalidate_ifnot"); Invnode; Invnode = Invnode->next_sibling("Invalidate_ifnot")) {
				this->Proc->new_rtag(rnum);
				prepare_tag(RTYPE_DIATOMIC, this->Proc->reaction_list[rnum]->tag[this->Proc->reaction_list[rnum]->num_tags-1], Invnode, true);	
				added_tags++;
			}	
			if(added_tags == 0) {
				White(true);
				cout << "   none.\n";
				White(false);
			}
		}

		Tnode = node->first_node("POLYATOMIC");
		if(Tnode) {
			this->Proc->reaction_list[rnum]->type = RTYPE_POLYATOMIC;
			this->Proc->reaction_list[rnum]->incomplete = false;

			cout << "\nType is: ";
			Green(true);
			cout << "polyatomic molecule";
			White(false);			
			cout << ", consisting of ";

			int num_ions = 0;
			if(Tnode->first_attribute("Ions")) 
				num_ions = atoi(Tnode->first_attribute("Ions")->value());
			if(num_ions > 0) {
				Green(true);
				cout << num_ions;
				White(false);	
				cout << " ions.\n";
				this->Proc->reaction_list[rnum]->number_of_ions = num_ions;
			} else {
				Error("..wait! Number of ions has not been specified.. This won't work.\n");
				this->Status = -1;
				return false;				
			}

			bool use_ion_matrix = false;
			double max_target_val = 0.0;
			double ambig_par = 1000.0;
			if(Tnode->first_attribute("UseIonMatrix")) {
				if(ToUp(Tnode->first_attribute("UseIonMatrix")->value()) == "TRUE")
					use_ion_matrix = true;			
			}

			this->Proc->reaction_list[rnum]->use_ion_matrix = use_ion_matrix;

			//if we want to use the ion matrix, we need the other parameters, as well...
			if(use_ion_matrix) {
				Green(true);
				cout << "Using the ion matrix to identify the ions..!\n";
				White(false);
				if(Tnode->first_attribute("MaxTargetValue")) 
					max_target_val = atof(Tnode->first_attribute("MaxTargetValue")->value());
				else {
					Error("Error: no max. target value provided. (You need that one when using the ion matrix..) This won't work.\n");
					this->Status = -1;
					return false;					
				}
				if(Tnode->first_attribute("Ambiguity")) 
					ambig_par = atof(Tnode->first_attribute("Ambiguity")->value());
				else {
					Error("Error: no ambiguity value provided. (You need that one when using the ion matrix..) This won't work.\n");
					this->Status = -1;
					return false;					
				}
				this->Proc->reaction_list[rnum]->max_target_value = max_target_val;
				this->Proc->reaction_list[rnum]->ambiguity_parameter = ambig_par;

				cout << "Max. target value is set to: ";
				Green(true);
				cout << this->Proc->reaction_list[rnum]->max_target_value << endl;
				White(false);
				cout << "Ambiguity parameter is set to: ";
				Green(true);
				cout << this->Proc->reaction_list[rnum]->ambiguity_parameter << endl;
				White(false);
			} else {
				White(true);
				cout << "Ion m/q assignment will be done by time-of-flight.\n";
				White(false);
			}

			int ions_defined = 0;
			//loop over all ion definitions..
			for(xml_node<> *IDnode = Tnode->first_node("Ion_def"); IDnode; IDnode = IDnode->next_sibling("Ion_def")) {

				double m = 0.0;
				double q = -1.0;
				double mean_tof = 0.0;
				double tmin = 0.0;
				double tmax = 0.0;
				double twidth = 0.0;

				// retrieve all parameters from cfg
				if(IDnode->first_attribute("Charge"))
					q = atof(IDnode->first_attribute("Charge")->value());
				if(IDnode->first_attribute("Mass"))
					m = atof(IDnode->first_attribute("Mass")->value());
				if(IDnode->first_attribute("Mean_tof")) 
					mean_tof = atof(IDnode->first_attribute("Mean_tof")->value());
				if(IDnode->first_attribute("TOF_width")) 
					twidth = atof(IDnode->first_attribute("TOF_width")->value());
				if(IDnode->first_attribute("TOF_min")) 
					tmin = atof(IDnode->first_attribute("TOF_min")->value());
				if(IDnode->first_attribute("TOF_max")) 
					tmax = atof(IDnode->first_attribute("TOF_max")->value());

				if(twidth < 1.0) {
					twidth = (tmax - tmin);
				}
				if(mean_tof < 1.0) {
					mean_tof = tmin + twidth/2.0;
				}
									
				this->Proc->reaction_list[rnum]->mass.push_back(m);
				this->Proc->reaction_list[rnum]->charge.push_back(q);
				this->Proc->reaction_list[rnum]->t_mean.push_back(mean_tof);
				this->Proc->reaction_list[rnum]->t_width[ions_defined] = twidth;

				if(q < -0.5) {
					Error("Error: Ion definition without specifying a charge. This won't work.\n");
					this->Status = -1;
					return false;					
				}
				if(q < 0.5) {
					if(this->Proc->reaction_list[rnum]->incomplete) {
						Error("Error: Found more than one neutral particle! This won't work.\n");
						this->Status = -1;
						return false;					
					} else {
						this->Proc->reaction_list[rnum]->incomplete = true;
						this->Proc->reaction_list[rnum]->number_of_ions--;
						White(true);
						cout << "Found ion with charge 0, assuming incomplete breakup.\n";
						White(false);
					}
				}
				if(m < 0.5 && !this->Spect->ion_side->linear_approximation) {
					Error("Error: Ion definition without specifying a mass. This won't work.\n");
					this->Status = -1;
					return false;					
				}
				if(mean_tof < 0.5 && this->Spect->ion_side->linear_approximation && q > 0.5) {
					Error("Error: Ion definition without specifying mean time-of-flight. This won't work.\n");
					this->Status = -1;
					return false;					
				}
				if(m > 0.5 && !this->Spect->ion_side->linear_approximation) {
					cout << "Ion ";
					Green(true);
					cout << ions_defined + 1;
					White(false);
					cout << " defined as ";
					White(true);
					cout << "m = "<< (int)this->Proc->reaction_list[rnum]->mass[ions_defined] <<" amu, q = "<< (int)this->Proc->reaction_list[rnum]->charge[ions_defined] << " au";
					White(false);
					cout << ".\n";
				}
				if(mean_tof > 0.5 && this->Spect->ion_side->linear_approximation) {
					cout << "Ion ";
					Green(true);
					cout << ions_defined + 1;
					White(false);
					cout << " defined as ";
					White(true);
					cout << "mean tof = "<< (int)this->Proc->reaction_list[rnum]->t_mean[ions_defined] <<" ns, q = "<< (int)this->Proc->reaction_list[rnum]->charge[ions_defined] << " au";
					White(false);
					cout << ".\n";
				}
				ions_defined++;
			}
			if(ions_defined != num_ions)  {
				Error("Error: Wrong number of ion definitions found...! This won't work.\n");
				this->Status = -1;
				return false;
			}

			if(Tnode->first_node("Randomize")) {
				this->Proc->reaction_list[rnum]->randomize_ions = true;
				xml_node<> *IDnode = Tnode->first_node("Randomize");
				if(IDnode->first_attribute("Ions")) {
					// Extract ion numbers from list..			
					std::string str = IDnode->first_attribute("Ions")->value(); 
					std::stringstream ss(str);
					for(int i = 0; ss >> i; ) {
				        this->Proc->reaction_list[rnum]->rnd_array[this->Proc->reaction_list[rnum]->rnd_count]=i;
						this->Proc->reaction_list[rnum]->rnd_count++;
						if(this->Proc->reaction_list[rnum]->rnd_count > this->Proc->reaction_list[rnum]->number_of_ions) {
							Error("Too many ions provided for swapping...! This won't work.\n");
							this->Status = -1;
							return false;							
						}
					}
				} else {
					this->Proc->reaction_list[rnum]->rnd_count = this->Proc->reaction_list[rnum]->number_of_ions;
					for(int i=0;i>this->Proc->reaction_list[rnum]->number_of_ions;i++) 		
						this->Proc->reaction_list[rnum]->rnd_array[i]=i;
				}
				cout << "Ions ";
				White(true);
				for(int i=0;i < this->Proc->reaction_list[rnum]->rnd_count-1;i++) 
					cout << this->Proc->reaction_list[rnum]->rnd_array[i] << ",";
				cout << this->Proc->reaction_list[rnum]->rnd_array[this->Proc->reaction_list[rnum]->rnd_count-1];
				White(false);
				cout <<" will be ";
				Green(true);
				cout <<"randomly swapped";
				White(false);
				cout << "..\n";
			}

			cout << "\nApplying conditions: \n";
			int added_tags = 0;
			//iterate over all "if"
			for(xml_node<> *Invnode = Tnode->first_node("Invalidate_if"); Invnode; Invnode = Invnode->next_sibling("Invalidate_if")) {
				this->Proc->new_rtag(rnum);
				prepare_tag(RTYPE_POLYATOMIC, this->Proc->reaction_list[rnum]->tag[this->Proc->reaction_list[rnum]->num_tags-1], Invnode, false);
				added_tags++;
			}
			//iterate over all "if_not"
			for(xml_node<> *Invnode = Tnode->first_node("Invalidate_ifnot"); Invnode; Invnode = Invnode->next_sibling("Invalidate_ifnot")) {
				this->Proc->new_rtag(rnum);
				prepare_tag(RTYPE_POLYATOMIC, this->Proc->reaction_list[rnum]->tag[this->Proc->reaction_list[rnum]->num_tags-1], Invnode, true);
				added_tags++;
			}	
			if(added_tags == 0) {
				White(true);
				cout << "   none.\n";
				White(false);
			}
		}

		Tnode = node->first_node("ELECTRON");
		if(Tnode) {
			this->Proc->reaction_list[rnum]->type++;
			cout << "\nReaction includes electrons...\n";

			cout << "\nApplying conditions to electrons: \n";
			int added_tags = 0;
			//iterate over all "if"
			for(xml_node<> *Invnode = Tnode->first_node("Invalidate_if"); Invnode; Invnode = Invnode->next_sibling("Invalidate_if")) {
				this->Proc->new_rtag(rnum);
				prepare_tag(RTYPE_ELEC, this->Proc->reaction_list[rnum]->tag[this->Proc->reaction_list[rnum]->num_tags-1], Invnode, false);	
				added_tags++;
			}
			//iterate over all "if_not"
			for(xml_node<> *Invnode = Tnode->first_node("Invalidate_ifnot"); Invnode; Invnode = Invnode->next_sibling("Invalidate_ifnot")) {
				this->Proc->new_rtag(rnum);
				prepare_tag(RTYPE_ELEC, this->Proc->reaction_list[rnum]->tag[this->Proc->reaction_list[rnum]->num_tags-1], Invnode, true);	
				added_tags++;
			}
			if(added_tags == 0) {
				White(true);
				cout << "   none.\n";
				White(false);
			}
		}

		// No reaction type definition was given, so abort!
		if(this->Proc->reaction_list[rnum]->type<0) {
			Error("You did not define a reaction type! This won't work.\n");
			this->Status = -1;
			return false;	
		}

		Tnode = node->first_node("Plot");
		if(Tnode) {
			cout << "\nSpecial histograms will be created:\n";
			cout << "===================================\n";

			xml_node<> *IDnode = Tnode->first_node("MFPAD");
			if(IDnode) {
				cout << "\nWill plot MFPADs:\n";
				if(IDnode->first_attribute("Type")) {
					if(ToUp(IDnode->first_attribute("Type")->value())=="LINEAR") {
						this->Proc->reaction_list[rnum]->MF_cond->type = 0;
						White(true);
						cout << "linearly polarized light\n";
						White(false);
					}
					if(ToUp(IDnode->first_attribute("Type")->value())=="CIRCULAR") {
						this->Proc->reaction_list[rnum]->MF_cond->type = 1;
						Green(true);
						cout << "circularly polarized light\n";
						White(false);
					}
					if(ToUp(IDnode->first_attribute("Type")->value())=="AUGER") {
						this->Proc->reaction_list[rnum]->MF_cond->type = 2;
						Green(true);
						cout << "Auger electron (no polarization dependence)\n";
						White(false);
					}
					if(this->Proc->reaction_list[rnum]->MF_cond->type == -1) {
						Error("...invalid MFPAD type! This won't work.\n");
						this->Status = -1;
						return false;					
					}
				} else {
					Error("...no MFPAD type provided! This won't work.\n");
					this->Status = -1;
					return false;	
				}
				
				if(IDnode->first_attribute("EE_min")) {
					this->Proc->reaction_list[rnum]->MF_cond->emin = atof(IDnode->first_attribute("EE_min")->value());
				}
				if(IDnode->first_attribute("EE_max")) {
					this->Proc->reaction_list[rnum]->MF_cond->emax = atof(IDnode->first_attribute("EE_max")->value());
				}
				if(IDnode->first_attribute("KER_min")) {
					this->Proc->reaction_list[rnum]->MF_cond->rmin = atof(IDnode->first_attribute("KER_min")->value());
				}
				if(IDnode->first_attribute("KER_max")) {
					this->Proc->reaction_list[rnum]->MF_cond->rmax = atof(IDnode->first_attribute("KER_max")->value());
				}
				cout << "Restricting plot to:\n";
				Green(true);
				cout << this->Proc->reaction_list[rnum]->MF_cond->emin;
				cout << " eV";
				White(false);
				cout << " << electron energy << ";
				Green(true);
				cout << this->Proc->reaction_list[rnum]->MF_cond->emax;
				cout << " eV\n";
				cout << this->Proc->reaction_list[rnum]->MF_cond->rmin;	
				cout << " eV";
				White(false);
				cout << " << KER << ";
				Green(true);
				cout << this->Proc->reaction_list[rnum]->MF_cond->rmax;
				cout << " eV\n";
				White(false);
			}

			IDnode = Tnode->first_node("LFPAD");
			if(IDnode) {
				cout << "\nWill plot labframe angular distributions for ";
				if(IDnode->first_attribute("Type")) {
					if(ToUp(IDnode->first_attribute("Type")->value())=="LINEAR") {
						this->Proc->reaction_list[rnum]->LF_cond->type = 0;
						White(true);
						cout << "linearly polarized light\n";
						White(false);
					}
					if(ToUp(IDnode->first_attribute("Type")->value())=="CIRCULAR") {
						this->Proc->reaction_list[rnum]->LF_cond->type = 1;
						Green(true);
						cout << "circularly polarized light\n";
						White(false);
					}
					if(this->Proc->reaction_list[rnum]->LF_cond->type == -1) {
						Error("...invalid LFPAD type! This won't work.\n");
						this->Status = -1;
						return false;					
					}
				} else {
					Error("...no LFPAD type provided! This won't work.\n");
					this->Status = -1;
					return false;	
				}
				
				if(IDnode->first_attribute("EE_min")) {
					this->Proc->reaction_list[rnum]->LF_cond->emin = atof(IDnode->first_attribute("EE_min")->value());
				}
				if(IDnode->first_attribute("EE_max")) {
					this->Proc->reaction_list[rnum]->LF_cond->emax = atof(IDnode->first_attribute("EE_max")->value());
				}

				cout << "Restricting plot to:\n";
				Green(true);
				cout << this->Proc->reaction_list[rnum]->LF_cond->emin;
				cout << " eV";
				White(false);
				cout << " << electron energy << ";
				Green(true);
				cout << this->Proc->reaction_list[rnum]->LF_cond->emax;
				cout << " eV\n";
				White(false);
			}

			IDnode = Tnode->first_node("Dalitz");
			if(IDnode) {
				if(IDnode->first_attribute("Ions")) {
					// Extract ion numbers from list..			
					std::string str = IDnode->first_attribute("Ions")->value(); 
					std::stringstream ss(str);
					int n=0;
					for(int i = 0; ss >> i; ) {
						if(n>2) {
							Error("Too many ions provided for Dalitz plot...! This won't work.\n");
							this->Status = -1;
							return false;							
						}
						this->Proc->reaction_list[rnum]->Dalitz_array[n++]=i;						
					}
					if(n<3) {
						Error("Warning: Less than 3 ions provided for Dalitz plot...! This doesn't seem right.\n");
						this->Status = -1;
						return false;	
					}
				}
				
				cout << "\nWill create a Dalitz plot of ions: ";
				White(true);
				for(int i=0;i < 2;i++) 
					cout << this->Proc->reaction_list[rnum]->Dalitz_array[i] << ",";
				cout << this->Proc->reaction_list[rnum]->Dalitz_array[2] << endl ;
				White(false);
			}

			IDnode = Tnode->first_node("Newton");
			if(IDnode) {
				if(IDnode->first_attribute("Ions")) {
					// Extract ion numbers from list..			
					std::string str = IDnode->first_attribute("Ions")->value(); 
					std::stringstream ss(str);
					int n=0;
					for(int i = 0; ss >> i; ) {
						if(n>2) {
							Error("Too many ions provided for Newton plot...! This won't work.\n");
							this->Status = -1;
							return false;							
						}
						this->Proc->reaction_list[rnum]->Newton_array[n++]=i;						
					}
					if(n<3) {
						Error("Warning: Less than 3 ions provided for Newton plot...! This doesn't seem right.\n");
						this->Status = -1;
						return false;
					}
				}
				
				cout << "\nWill create a Newton plot of ions: ";
				White(true);
				for(int i=0;i < 2;i++) 
					cout << this->Proc->reaction_list[rnum]->Newton_array[i] << ",";
				cout << this->Proc->reaction_list[rnum]->Newton_array[2] << endl ;
				White(false);			
			}

			IDnode = Tnode->first_node("Two_elec");
			if(IDnode) {		
				cout << "\nWill create two-electron plots.";
				this->Proc->reaction_list[rnum]->two_elec = true;
				if(IDnode->first_attribute("Master_EE_min")) {
					this->Proc->reaction_list[rnum]->e_master_min = atof(IDnode->first_attribute("Master_EE_min")->value());
				}
				if(IDnode->first_attribute("Master_EE_max")) {
					this->Proc->reaction_list[rnum]->e_master_max = atof(IDnode->first_attribute("Master_EE_max")->value());
				}
				if(this->Proc->reaction_list[rnum]->e_master_min<0.0 && this->Proc->reaction_list[rnum]->e_master_max>1000000.0) {
					cout << "'First' electron is master...\n";
				} else {
					cout << "Master electron properties:\n";
					Green(true);
					cout << this->Proc->reaction_list[rnum]->e_master_min;
					cout << " eV";
					White(false);
					cout << " << electron energy << ";
					Green(true);
					cout << this->Proc->reaction_list[rnum]->e_master_max;
					cout << " eV\n";
					White(false);			
				}
			}

			IDnode = Tnode->first_node("Photon_scan");
			if(IDnode) {		
				cout << "\nWill create plots for photon energy scans.\n";
				this->Proc->reaction_list[rnum]->photon_scan = true;

				if(IDnode->first_attribute("Scan_channel")) {
					this->Proc->reaction_list[rnum]->ph_scan_channel = atoi(IDnode->first_attribute("Scan_channel")->value());
				}
				if(IDnode->first_attribute("hv_min")) {
					this->Proc->reaction_list[rnum]->ph_min = atof(IDnode->first_attribute("hv_min")->value());
				}
				if(IDnode->first_attribute("hv_max")) {
					this->Proc->reaction_list[rnum]->ph_max = atof(IDnode->first_attribute("hv_max")->value());
				}
				if(IDnode->first_attribute("hv_step")) {
					this->Proc->reaction_list[rnum]->ph_step = atof(IDnode->first_attribute("hv_step")->value());
				}

				if(this->Proc->reaction_list[rnum]->ph_scan_channel <0) {
					Error("Error: Please specify channel containing scan data.\n");
					this->Status = -1;
					return false;
				}
				if(this->Proc->reaction_list[rnum]->ph_min <0) {
					Error("Error: Please specify lower histogram boundary.\n");
					this->Status = -1;
					return false;
				}
				if(this->Proc->reaction_list[rnum]->ph_max <0) {
					Error("Error: Please specify upper histogram boundary.\n");
					this->Status = -1;
					return false;
				}
				if(this->Proc->reaction_list[rnum]->ph_step <0) {
					Error("Error: Please specify histogram stepsize.\n");
					this->Status = -1;
					return false;
				}

				cout << "Histograms will be defined from ";
				Green(true);
				cout << this->Proc->reaction_list[rnum]->ph_min;
				White(false);
				cout << " eV to ";
				Green(true);
				cout << this->Proc->reaction_list[rnum]->ph_max;
				White(false);		
				cout << " eV, stepsize is ";
				Green(true);
				cout << this->Proc->reaction_list[rnum]->ph_step;
				White(false);
				cout << " eV.\n\n";
			}
		}

		cout << "\n| End of definition of reaction '"<< this->Proc->reaction_list[rnum]->name <<"'.\n";
		cout << "===============================================================---";
		i++;
	}

	this->maxNumCHHists = this->Proc->number_of_reactions_defined * HISTS_PER_CHANNEL * 16;

	if(this->Proc->number_of_reactions_defined>0) {
		cout << "\nSetup of reactions:";
		Green(true);
		gotoX(SUCCESSX);
		cout << " ...success!\n";
		White(false);
		cout << "====================================================================--\n\n";
	}
	return true;
}

bool ColAHelL::SetupParticleClasses() {
	for(int j=0;j<64;j++) {
		this->e[j] = new electron_class();
		this->i[j]  = new ion_class();
	}
	this->mol = new diatomic_class();
	this->big_mol = new polyatomic_class();
	return true;
}

bool ColAHelL::SetupHistogramming() {
	this->HGrams = new histograms_class();
	this->PrsHist = new histo_handler();
	return true;
}

void ColAHelL::CalculateTOFs() {
	this->CTof->calc(evt);
}

std::vector<CH_event_struct> ColAHelL::PresortEvent() {
	
	prs_evt_out.clear();
	for(int i=0;i<num_presorters;i++) {
		// reset output event
		ResetEventStruct(&this->presorters[i]->evto, this->evt->event_number);
		// give current event to presorter
		this->presorters[i]->evti = *this->evt;
		// run presorter + append out event if successful
		if(this->presorters[i]->Run()) {
			prs_evt_out.push_back(this->presorters[i]->evto);
		}
	}
	return prs_evt_out;
}

bool ColAHelL::CheckEvent(CH_event_struct * evt, int channel) {
	bool result = true;
	for(int i=0;i<num_presorters;i++) {

		ResetEventStruct(&this->presorters[i]->evto, this->evt->event_number);
		this->presorters[i]->evti = *this->evt;
		
		if(this->presorters[i]->Get_flag() == channel)
		 result &= this->presorters[i]->Run();
	}
	return result;
}

void ColAHelL::ProcessEvent(CH_event_struct * evt, int ReactionDefNum) {
	this->Proc->process_all(evt, this->e, this->i, this->mol, this->big_mol, ReactionDefNum);
	this->current_reaction = this->Proc->get_current_reaction();
}
void ColAHelL::ProcessEvent(int ReactionDefNum) {
	this->Proc->process_all(this->evt, this->e, this->i, this->mol, this->big_mol, ReactionDefNum);
	this->current_reaction = this->Proc->get_current_reaction();
}

void ColAHelL::ResetEventStruct(CH_event_struct * evt, __int64 event_number) {
	evt->event_number = event_number;
	evt->bunchmarker = 0.0;
	evt->reaction = 0;
	evt->scanval.clear();
	
	evt->e.x.clear();
	evt->e.y.clear();
	evt->e.time.clear();
	evt->e.tof.clear();
	evt->e.method.clear();
	evt->e.num_hits = 0;

	evt->r.x.clear();
	evt->r.y.clear();
	evt->r.time.clear();
	evt->r.tof.clear();
	evt->r.method.clear();
	evt->r.num_hits = 0;

	evt->p.x.clear();
	evt->p.y.clear();
	evt->p.time.clear();
	evt->p.tof.clear();
	evt->p.method.clear();
	evt->p.num_hits = 0;
}

void ColAHelL::SetupEventStruct() {
	this->evt = new CH_event_struct();
}

double ColAHelL::GetEDetSize() {
	if(this->HGrams->edet_size>1.0)
		return this->HGrams->edet_size;
	return -1.0;
}
double ColAHelL::GetIDetSize() {
	if(this->HGrams->rdet_size>1.0)
		return this->HGrams->edet_size;
	return -1.0;
}
double ColAHelL::GetPDetSize() {
	if(this->HGrams->pdet_size>1.0)
		return this->HGrams->edet_size;
	return -1.0;
}
int ColAHelL::GetMaxPrsHists() {
	return this->maxNumPrsHists;
}
histo_handler* ColAHelL::GetPrsHists() {
	return this->PrsHist;
}

int ColAHelL::GetMaxUserHists() {
	return this->HGrams->MaxUserHists;
}
int ColAHelL::GetMaxHists() {
	return this->maxNumCHHists;
}
histo_handler* ColAHelL::GetHists(int Type) {
	switch(Type) {
		case 0:
			return this->HGrams->Hist_e;
		case 1:
			return this->HGrams->Hist_ions;
		case 3:
			return this->HGrams->Hist_User;
		default:
			cout << endl << "Critical error detected! Unknown histohandler requested!";
			return 0;
	}
}

bool ColAHelL::IsReactionDefined(int reaction) {
	int reac = this->Proc->find_reaction(reaction);
	if(reac<0)
		return false;
	else 
		return true;
}

int ColAHelL::EvtBelongsToReactions(int channel) {
	return this->Proc->EvtBelongsToReactions(channel);
}

void ColAHelL::SetMasterFolder(const char* folder) {
	this->HGrams->SetMasterFolder(folder);
}

void ColAHelL::FillCHHists() {
	//
	this->HGrams->electron_valid = this->Proc->at_least_one_e_valid;
	this->HGrams->ion_valid = this->Proc->at_least_one_r_valid;
	// plot everything...
	// The processor updated global things like "current_reaction"!
	this->HGrams->UserAnalysis(this->current_reaction, this->e, this->i, this->mol, this->big_mol, this->evt);
	this->HGrams->plot(this->current_reaction, this->e, this->i, this->mol, this->big_mol, this->evt);
}

mom_cor_param* ColAHelL::GetMomCorParam(int Det) {
	switch(Det) {
		case 0:
			return this->Proc->e_mom_fac;
			break;
		case 1:
			return this->Proc->r_mom_fac;
			break;
		case 2:
			return this->Proc->p_mom_fac;
			break;
		default:
			return this->Proc->e_mom_fac;
	}
}

cor_param* ColAHelL::GetCorParam(int Det) {
	switch(Det) {
		case 0:
			return this->Proc->e_fac;
			break;
		case 1:
			return this->Proc->r_fac;
			break;
		case 2:
			return this->Proc->p_fac;
			break;
		default:
			return this->Proc->e_fac;
	}
}


