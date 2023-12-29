#ifndef DG_UNIT_H
#define DG_UNIT_H

#include <vector>
#include <unordered_map>
//#include <array>
#include "dg_mapped_geometry_class.h"


/// @brief
/// Structure to deal with the mortar.
struct Mortar{

	int n_max;	// maximum poly order between left and right element
	
	int l_max;	// maximum h-refinement level between left and rigth elemnt. 

	double a_l, a_r;	// offset of the mortar with respect to the larger element.
	double b_l, b_r;	// scale of the mortar with respect to the larger element. 

	std::vector<double> psi_l;	// L2 projection on the left interface of the mortar
	std::vector<double> psi_r; 	// L2 projection on the right interface of the mortar

	std::vector<double> nflux;	// numerical flux on the right element interface
};

/// @brief
/// Element unit. All the needed information is stored as a unit.
/// @param child_position ith-child
class Unit{

public:

	int n, m; 	// polynomial orders (x, y)	

	int index[3]{0, 0, 0}; 	// element index[i, j, k]
 
	char status; 	// Hilbert status

	int child_position{};	// ith-child
	
	struct Face;	
	std::vector<std::vector<Face>> facen;	// a structure for the neighbour on each face, initialized in dg_unit.cpp

	double xcoords[2]{0.0, 0.0};	// x coordinates
	double ycoords[2]{0.0, 0.0};	// y coordinates

	std::unordered_map<int, std::vector<double>> solution;	// solutions <equ, (n + 1) * (m + 1)>

	std::vector<double> solution_int_l;	// solution on the element left interface (porder + 1)
	std::vector<double> solution_int_r;	// solution on the element right interface
	
	// ghost space to store the neighbours' solutions on the interface
	// hash <neighbour's key, solution_int>
	std::unordered_map<long long int, std::vector<double>> ghost; 

	std::vector<double> nflux_l;	// numerical flux on the left interface
	std::vector<double> nflux_r;	// numerical flux on the right interface

	std::unordered_map<int, std::vector<double>> solution_time_der;	// solution time derivative <equ, (n + 1) * (m + 1)>

	std::unordered_map<int, std::vector<double>> G; 	// Runge-Kutta 3rd <equ, (n + 1) * (m + 1)>
	
	std::vector<double> ref_x;		// element boundary in x direction
	std::vector<double> ref_y;		// element boundary in y direction

	Unit* next = nullptr;	// pointer to the next Unit

	bool hrefine = false;	// h-refinemnt 
	bool prefine = false;	// p-refinemnt 
	bool coarsen = false; 

	// mortar ---------------------------------------------
	Mortar mortar;	// mortar struct instance	
	//-----------------------------------------------------

	// constructor (default)
	Unit();

	MappedGeometry holdmetrics;
		
};


/// @param face_type interface type of current direction.
/// @param hlevel hlevel of neighbour element.
/// @param porderx polynomial order of neighbour element in x direction.
/// @param pordery polynomial order of neighbour element in y direction.
/// @param key neighbour's key.
/// @param rank if face_type == 'M' rank == facing process's MPI.
struct Unit::Face{
	
	char face_type;
	
	int hlevel;
	
	int porderx;

	int pordery;
	
	long long int key;
	
	int rank;

	std::vector<double> ref_x; 
	std::vector<double> ref_y;


	// default construtor1
	Face(); 

	// constructor2
	Face(char c, int h, int nx, int ny, long long int k, int r, double* ref1, double* ref2);

	// constructor3
	Face(char c, int h, int nx, int ny, long long int k, int r, std::vector<double>& ref1, std::vector<double>& ref2);
	
	// copy constructor ------------------------------------------------
	Face(const Face& face);	// copy another instance

	Face(const std::vector<Face>::iterator p);	// copy by pointer
	//-------------------------------------------------------------------
	
};
	


#endif
