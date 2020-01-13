//tissue.h
//=================
//Include Guards
#ifndef _TISSUE_H_INCLUDED_
#define _TISSUE_H_INCLUDED_
//=========================
//forward declarations

//=======================
//Include dependencies
#include <string>
#include <vector>
#include <fstream>
#include <memory>
#include <random>
#include "phys.h"
#include "coord.h"
#include "node.h"
#include "cell.h"
#include "externs.h"
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
//=========================
// Tissue Class Declaration

//class Tissue: public enable_shared_from_this<Tissue> {
class Tissue {
	private:
		// We'll need a better data structure for later
		vector<shared_ptr<Cell>> cells;
		int num_cells;
		int num_deleted;
		mt19937 gen;
		vector<int> dist1;
		vector<int> dist2;
		vector<int> dist3;
		vector<int> dist4;
		vector<int> counts;	
	public:
		Tissue(string filename, mt19937 gen);
		void get_Cells(vector<shared_ptr<Cell>>& cells);
		//set/get the number of cells in the tissue
		void update_Num_Cells(shared_ptr<Cell>& new_Cell);
		int  get_num_cells() {return num_cells;}
		Coord Compute_L1_AVG();
		int return_counts(int index);
		void set_up_counts();
		void set_counts(int index);
		void assign_dist_vecs(vector<int> dist1, vector<int>dist2, vector<int> dist3, vector<int> dist4);
		int get_next_random(int dist, int count);
		void update_Signal();
		void update_growth_direction();
		void update_Neighbor_Cells();
		void add_Wall(int Ti);
		void delete_Wall(int Ti);
		void update_Adhesion();
		double get_normal_number(double mean, double sigma); //mt19937 get_Random_Generator(){return gen;}
		//not in use
		void update_Linear_Bending_Springs();
		
		void update_Cell_Cycle(int Ti);
		void division_check();
		void calc_New_Forces(int Ti);
		void update_Cell_Locations(int Ti);
		
		//For tissue halting purposes

		bool terminal_Tissue();
		//stuff for data output
		void plot_direction_vec(ofstream& ofs);
		void print_Data_Output(ofstream& ofs);

		void locations_output(ofstream& ofs,bool cytoplasm);
		//int update_VTK_Indices();
		//void print_VTK_File(ofstream& ofs);

		//void locations_output(ofstream& ofs);
		int update_VTK_Indices(bool cytoplasm);
		void print_VTK_File(ofstream& ofs, bool cytoplasm);
		void print_VTK_Direction_File(ofstream& ofs);	
		void inc_Num_Deleted();
		//Destructor
		~Tissue();
		//Debugging
		void NAN_CATCH(int Ti);
		void BAD_CATCH(int call, int Ti);
};


//===========================
//End of file

#endif
