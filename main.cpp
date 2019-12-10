#ifndef _MAIN_INCLUDE_
#define _MAIN_INCLUDE_
// Main.cpp
//============================

//===========================
// Include Dependencies
#include <iostream>
#include <math.h>
#include <string>
#include <vector>
#include <fstream>
#include <ctime>
#include <stdio.h>
#include <memory>
//#include <random>
#include <functional>
#include <chrono>
#include "phys.h"
#include "coord.h"
#include "node.h"
#include "cell.h"
#include "tissue.h"
#include "rand.h"
//==============================

using namespace std;

//FREQUENTLY CHANGED PARAMETERS
bool OUT_OF_PLANE_GROWTH = true;
//EXPERIMENTAL PARAMTERS
int DIV_MECHANISM = 0;
double WUS_RAD_CONTRACTION_FACTOR = 1;
double CK_RAD_CONTRACTION_FACTOR = 1;
int TENSILE_CALC = 1;
//Must be declared in externs.h
//For clarity, listed as comments in phys.h

int main(int argc, char* argv[]) {

	// reads in name of folder that 
	// stores vtk files, given in run.sh
	//usually called Animate1
	string anim_folder = argv[1];
	//Reads in the name of folder that stores VTK files for
	//visualization without cytoplasm nodes.
	string no_cyt_folder = argv[4];
	//reads in name of folder that 
	//stores data output, given in run.sh
	string locations_cyt_folder = argv[3];
	//this is a folder that holds data output
	//without cytoplasm node info
	string locations_no_cyt_folder = argv[2];
	//keep track of time
	for (int i = 1; i < argc; i++) { 
		if (!strcmp(argv[i], "-WR")) { 
			WUS_RAD_CONTRACTION_FACTOR = stod(argv[i+1]);
		} else if (!strcmp(argv[i], "-CKR")) { 
			CK_RAD_CONTRACTION_FACTOR = stod(argv[i+1]);
		} else if (!strcmp(argv[i], "-div")) { 
			DIV_MECHANISM = stoi(argv[i+1]);
		} else if (!strcmp(argv[i], "-TC")) { 
			TENSILE_CALC = stoi(argv[i+1]);
		} else if (!strcmp(argv[i], "-OOP_off")) { 
			OUT_OF_PLANE_GROWTH = false;
		}
	}
	if (DIV_MECHANISM == 0) { 
		cout << "DIV_MECHANISM not set.  Exiting..." << endl;
		exit(1);
	} else if (DIV_MECHANISM > 3 || DIV_MECHANISM < 1) { 
		cout << "DIV_MECHANISM Input failed. Exiting..." << endl;
		exit(1);
	}
	int start = clock();	
	mt19937::result_type seed = time(0);
	
	//distrubution for different cell cycle lengths
	auto normal_rand1 = bind(normal_distribution<double> (15800,2300),mt19937(seed));
	auto normal_rand2 = bind(normal_distribution<double> (19400,1800),mt19937(seed));
	auto normal_rand3 = bind(normal_distribution<double> (24800,3600), mt19937(seed));
	auto normal_rand4 = bind(normal_distribution<double> (44600,18200), mt19937(seed));
	vector<int> dist1;
	vector<int> dist2;
	vector<int> dist3;
	vector<int> dist4;
	for (unsigned int i=0;i<600;i++){
		dist1.push_back((int) normal_rand1());
		dist2.push_back((int) normal_rand2());
		dist3.push_back((int) normal_rand3());
		dist4.push_back((int) normal_rand4());
	}
	//.txt file that tells initial
	//cell configuration 
	//cout << "before cell file is read in" << endl;
	string init_tissue = "staggered_generated.txt";
	//cout << "Read in cell starter" << endl;	


	//instantiate tissue
	//new cell and node objects
	//are made in this call
	Tissue growing_Tissue(init_tissue);
	//TIssue growing_Tissue_experiment(init_tissue);
	growing_Tissue.assign_dist_vecs(dist1, dist2, dist3, dist4);
	//cout << "Finished creating Cells" << endl;
	growing_Tissue.update_Signal();
	//cout << "Signal" << endl;
	//growing_Tissue.update_growth_direction();
	//cout << "growth direction" << endl;
	//parameters for time step
	double numSteps = 500;

	// Variable for dataoutput
	int digits;
	string format = ".vtk";
	string Number;
	string initial = "/Plant_Cell_";
	string Filename;
	ofstream ofs_anim;
	int out = 0;

	string noCyt_Filename;
	ofstream ofs_noCyt;
	string no_cyt_initial = "/Plant_Cell_NC_";

	//int digits2;
	string Number2;
	string Locations_no_cyt;
	ofstream ofs_loc_no_cyt;
	int out2 = 0;

	//int digits3;
	string Number3;
	string Locations_cyt;
	ofstream ofs_loc_cyt;
	string locations_initial = "/Locations_";
	int out3 = 0;

	//(No longer use this loop, instead flag for terminal.)
	//loop for time steps
	//which matlab file tells you how many
	//seconds each time step reprsents (2.5?)
	//for(int Ti = 0; Ti*dt < numSteps; Ti++) {

	int terminal_timer = 0;
	bool is_terminal = false;
	while(terminal_timer < 3000) {
		//keep track of simulation runs
		if (Ti %1000 == 0 && !is_terminal) {
			is_terminal = growing_Tissue.terminal_Tissue();
			if (is_terminal) {

			}
		} else { 
			terminal_timer++;
		}

		// Tissue Growth

		//fills vector of neighbor cells for each cell
		//in tissue class this goes through each cell and calls
		//updated neighbor on each cell
		//growing_Tissue.BAD_CATCH(1,Ti);
		if(Ti%5000==0) {
			//cout << "Find Neighbors" << endl;
			growing_Tissue.update_Neighbor_Cells();
		}	

		if(Ti == 10000){
			growing_Tissue.update_Signal();
		}
		if(Ti % 30000 == 0){
			//cout << "update signal" << endl;
			growing_Tissue.update_Signal();
			//growing_Tissue.update_growth_direction();
		}
		//adds one new cell wall node per cell everytime it is called
		//dont call it right away to give cell time to find initial configuration
		if(Ti >= 10000){
			if(Ti%1000 == 0){	
				//cout << "add new cell wall nodes if needed" << endl;
				growing_Tissue.add_Wall(Ti);
			}
		}
		//was used previously to help stability of cells
		//deletes a cell wall node if too close together
		if(Ti > 10000){
			//if(Ti%1000 == 500){}
			if (Ti % 1000 == 500) { 
				//cout << "delete wall" << endl;
				growing_Tissue.delete_Wall(Ti);
				//growing_Tissue.delete_Wall_Node_Check(Ti);
			}
		}

		//make adhesion pairs for each cell
		if(Ti < 10000){
			if(Ti%1000 == 0) {
				//cout << "adhesion early" << endl;
				growing_Tissue.update_Adhesion();
			}
		} else {	
			if(Ti%50000 == 0) {
				//cout << "adhesion"<< endl;
				growing_Tissue.update_Adhesion();
			}
		}
		//adds internal node according to 
		//individual cell growth rate
		if(Ti >= 10000){
			//cout << "cell cycle" << endl;
			growing_Tissue.update_Cell_Cycle(Ti);
		}
		//will divide cell if time
		cout << "divide necessary cells" << endl;
		if(Ti >= 10000){
			growing_Tissue.division_check();
		}

		//Calculate new forces on cells and nodes
		cout << "forces" << endl;
		growing_Tissue.calc_New_Forces(Ti);

		//Update node positions
		//cout << "locations" << endl;
		growing_Tissue.update_Cell_Locations(Ti);	

		//cout << "Finished" << endl;

		// print to dataOutput and vtk files
		if(Ti%2500==0) {
			digits = ceil(log10(out + 1));
			if (digits == 1 || digits == 0) {
				Number = "0000" + to_string(out);
			}	
			else if (digits == 2) {
				Number = "000" + to_string(out);
			}	
			else if (digits == 3) {
				Number = "00" + to_string(out);
			}
			else if (digits == 4) {
				Number = "0" + to_string(out);
			}

			Filename = anim_folder + initial + Number + format;

			ofs_anim.open(Filename.c_str());
			//true is in reference to printing cytoplasm
			growing_Tissue.print_VTK_File(ofs_anim,true);
			ofs_anim.close();	

			noCyt_Filename = no_cyt_folder + no_cyt_initial + Number + format;

			ofs_noCyt.open(noCyt_Filename.c_str());
			//false is in reference to printing cytoplasm
			growing_Tissue.print_VTK_File(ofs_noCyt,false);
			ofs_noCyt.close();	

			out++;
		}
		//growing_Tissue.BAD_CATCH(11,Ti);
		/*if(Ti%1000==0) {
			digits2 = ceil(log10(out2 + 1));
			if (digits2 == 1 || digits2 == 0) {
				Number2 = "0000" + to_string(out2);
			}	
			else if (digits2 == 2) {
				Number2 = "000" + to_string(out2);
			}	
			else if (digits2 == 3) {
				Number2 = "00" + to_string(out2);
			}
			else if (digits2 == 4) {
				Number2 = "0" + to_string(out2);
			}

			nem_Filename = nem_folder+ initial2 + Number2 + format;

			ofs_nem.open(nem_Filename.c_str());
			growing_Tissue.print_VTK_Direction_File(ofs_nem);
			ofs_nem.close();	
			out2++;
		}*/
		//data output from simulations
		//for cell center etc
		/*if(Ti%10000 == 1){
		  nem_Filename = nematic_folder + nem_initial + to_string(Number2) + ".txt";
		  ofs_nem.open(nem_Filename.c_str());
		  growing_Tissue.nematic_output(ofs_nem);
		  ofs_nem.close();
		  Number2++;
		  }
		  if(Ti%1000 == 1){
		  locations_Filename = locations_folder + locations_initial + to_string(Number3) + ".txt";
		  ofs_loc.open(locations_Filename.c_str());
		  growing_Tissue.locations_output(ofs_loc);
		  ofs_loc.close();
		  Number3++;
		  }

			nem_Filename = nematic_folder + nem_initial + to_string(Number2) + ".txt";
			ofs_nem.open(nem_Filename.c_str());
			growing_Tissue.nematic_output(ofs_nem);
			ofs_nem.close();
			Number2++;
		}*/
		//locations with cyt nodes
		if(Ti%10000 == 0){
			Locations_no_cyt = locations_no_cyt_folder + locations_initial + to_string(out2) + ".txt";
			ofs_loc_no_cyt.open(Locations_no_cyt.c_str());
			growing_Tissue.locations_output(ofs_loc_no_cyt,false);
			ofs_loc_no_cyt.close();
			out2++;
			Locations_cyt = locations_cyt_folder + locations_initial + to_string(out3) + ".txt";
			ofs_loc_cyt.open(Locations_cyt.c_str());
			growing_Tissue.locations_output(ofs_loc_cyt,true);
			ofs_loc_cyt.close();
			out3++;
		}
		//growing_Tissue.BAD_CATCH(12,Ti);

	}
	/*nem_Filename = nematic_folder + nem_initial + to_string(Number2) + ".txt";
	  ofs_nem.open(nem_Filename.c_str());
	  growing_Tissue.nematic_output(ofs_nem);
	  ofs_nem.close();*/

	int stop = clock();

	cout << "Time: " << (stop - start) / double(CLOCKS_PER_SEC) * 1000 << endl;

	return 0;

}










#endif
