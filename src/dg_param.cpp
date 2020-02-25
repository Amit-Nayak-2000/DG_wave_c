#include <iostream>
#include <mpi.h>
#include <string>
#include "dg_param.h"

//variable you must change--------------------------------------------

/// @brief mesh file information
/// @param fileplace The path of mesh file and mesh file name
namespace fileinfo{
	const std::string fileplace = "../gmsh_files/4_elements.msh";
}

/// @brief Domain size

/// @param exp_x element number in x direction(exponential order), i.e. 2^(exp_x)

/// @param exp_y element number in y direction(exponential order), i.e. 2^(exp_y)

/// @param gx_l domain left boundary (x direction, must >= 0) 
/// @param gx_r domain right boundary (x direction, must >= 0) 
/// @param gy_l domain left boundary (y direction, must >= 0) 
/// @param gy_r domain right boundary (y direction, must >= 0) 
/// @param nmin minimum polynomial degree in x and y direction
/// @param nmax maximum polynomial degree in x and y direction
/// @param hlevel_max maximum h-refinement level. 
namespace grid{
	const int exp_x = 1; 
	const int exp_y = 1; 
	
	const double gx_l = 0.0;
	const double gx_r = 1.0; 
	const double gy_l = 0.0;
	const double gy_r = 1.0; 

	const int nmin = 2;	// x direction
	const int nmax = 4;

	const int hlevel_max = 2;	
};
//---------------------------------------------------------------------


// variables you could change------------------------------------------

/// @brief
/// Time related variables
/// @param t_total total time integal
/// @param nt time step number
namespace dg_time{
	const double t_total = 1.0;

	const int nt = 2;

};

/// @brief
/// Function related parameters
/// @param num_of_equation  number of equations
/// @param C wave speed
namespace dg_fun{

	const int num_of_equation = 3;

	const double C = 1.0;
};

//----------------------------------------------------------------------

// variables you do not need to change----------------------------------

/// @brief mpi variables
/// @param rank process rank
/// @param num_proc total number of processor
namespace mpi{
	int rank;    
	int num_proc;    
};

//-----------------------------------------------------------------------
