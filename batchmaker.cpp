#include <iostream>
#include <cstring>
#include <fstream>
#include <math.h>
#define PI 3.14159265
using namespace std;

//Makes a file "18rings.txt" that represents an initial condition of 18 concentric circles

int main(int argc, char* argv[]) {

	ofstream ofs;
	char test[100] = "";
	char cores[100] = "";
	char p[100] = ""; 
	char hours[100] = "";
	char mins[100] = "";
	char nodes[100] = "";
	char mem[100] = "";
	char DPM[100] = "";
	char WUS_CF[100]= "";
	char CK_CF[100] = "";
	int divDataCutoff;

	for (int i = 1; i < argc; i++) { 
		if (!strcmp(argv[i], "-p")) { 
			strcpy(p,argv[i+1]);
		} else if (!strcmp(argv[i], "-test")) { 
			strcpy(test,argv[i+1]);
		} else if (!strcmp(argv[i], "-hours")) { 
			strcpy(hours,argv[i+1]);
		} else if (!strcmp(argv[i], "-cores")) { 
			strcpy(cores,argv[i+1]);
		} else if (!strcmp(argv[i], "-mins")) { 
			strcpy(mins,argv[i+1]);
		} else if (!strcmp(argv[i], "-nodes")) { 
			strcpy(nodes, argv[i+1]);
		} else if (!strcmp(argv[i], "-mem")) { 
			strcpy(mem, argv[i+1]);
		} else if (!strcmp(argv[i], "-help")) { 
			goto helplabel;
		} else if (!strcmp(argv[i], "-WR")) { 
			strcpy(WUS_CF, argv[i+1]);
		} else if (!strcmp(argv[i], "-CKR")) { 
			strcpy(CK_CF, argv[i+1]);
		} else if (!strcmp(argv[i], "-div")) { 
			strcpy(DPM, argv[i+1]);
		}
	}

	ofs.open("AUTO_BATCH.sh");

	if (!strlen(test)) {
		cout << "Test name: ";
		cin >> test; 
		if (cin.fail()) { 
			cout << "Integer only." << endl;
			ofs.close();
			return 1;
		}
		cin.clear();
	}
	if (!strlen(cores)) { 
		cout << "Number of cores: ";
		cin >> cores; 
		if (cin.fail()) { 
			cout << "Integer only." << endl;
			ofs.close();
			return 1;
		}
	}
	if (!strlen(hours) && !strlen(mins)) { 
		cout << "Hours: ##";
		cin >> hours; 
		if (cin.fail()) { 
			ofs.close();
			return 1;
		}
		cout << "Minutes: ##";
		cin >> mins; 
		if (cin.fail()) { 
			ofs.close();
			return 1;
		}
	} else {
		if (!strlen(hours)) strcpy(hours, "00");
		if (!strlen(mins)) strcpy(mins, "00");
	}
	if (!strlen(p)) {
		cout << "Partition: ";
		cin >> p; 
		if (cin.fail()) { 
			ofs.close();
			return 1;
		}
	}
	if (!strlen(mem)) strcpy(mem, "2");
	if (!strlen(nodes)) strcpy(nodes,"1");
	ofs << "#!/bin/bash -l\n";
	ofs << "#SBATCH --nodes=" << nodes << "\n";
	ofs << "#SBATCH --ntasks=1\n";
	ofs << "#SBATCH --cpus-per-task=" << cores << "\n";
	ofs << "#SBATCH --mem-per-cpu=" << mem << "G\n";
	ofs << "#SBATCH --time=0-" << hours << ":" << mins << ":00\n";
	ofs << "#SBATCH --output=my"<< test << ".stdout\n";
	ofs << "#SBATCH --job-name=\"" << test << "\"\n";
	ofs << "#SBATCH -p " << p << " \n";

	ofs << "export OMP_NUM_THREADS " << cores << "\n";
	ofs << "mkdir Animate_Cyt_" << test << "\n";
	ofs << "mkdir Nematic_test_1 \n";
	ofs << "mkdir Locations_test_1\n";
	ofs << "mkdir Animate_No_Cyt_" << test << "\n";
	ofs << "./program Animate_Cyt_" << test << " Locations_test_" << test << " Nematic_test_" 
		<< test << " Animate_No_Cyt_" << test << 
		 " -WR " << WUS_CF << " -CKR " << CK_CF << " -div "
		<< DPM << "\n";
	ofs.close();


	return 0;


	helplabel:
	cout << "This file generates a SBATCH script called AUTO_BATCH.sh.\n";
	cout << "Options and their defaults are:" << endl;
	cout << "-p <partition>.  Omitting this will give a console prompt." << endl;
	cout << "-test <name>.  Ommitting this will give a console prompt." << endl;
	cout << "-hours <##> and -mins <##>.  Two digits only.  Filling only of these" << 
		" will default the other to 00." << endl;
	cout << "-cores <#>. Ommitting this will give a console prompt to define it." <<endl;
	cout << "-nodes <#>. Defaults to 1." << endl;
	cout << "-mem <#>.  Specifies number of gigabytes of memory per CPU. Defaults to 2." << endl;
	cout << "Example:  ./a.out -p short -hours 02 -cores 12 -test mech_div_5" << endl;

	return 0;
}
