#include <string>
#include "dg_param.h"
#include <cmath>

//variable you must change--------------------------------------------

/// @brief mesh file information
/// @param fileplace The path of mesh file and mesh file name
/// @param output_place output directory. 
namespace fileinfo{
	 std::string fileplace;
//	const std::string fileplace = "../gmsh_files/strong/1024.msh";

	const std::string output_place = "../outputs/";

	const std::string eff_filename = "../efficiency/eff_no_LB";	// efficiency output

	const std::string crosssection_filename = "../cross_section/cs_4elem_hp.csv";	// write the result on the cross-setion

	const std::string exact_error_filename = "../exact_error/error_refine.dat";	// exact error 
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
	int exp_x; 
	int exp_y; 
	
	double gx_l;
	double gx_r; 
	double gy_l;
	double gy_r; 

	int nmin;	
	int nmax;

	int hlevel_max;	
};
//---------------------------------------------------------------------


// variables you could change------------------------------------------

/// @brief
/// Time related variables
/// @param t_total total time integal
/// @param nt time step number
namespace dg_time{

	double t_total;
//	const double t_total = 0.5;

	int nt;

};

/// @brief
/// Function related parameters
/// @param num_of_equation  number of equations
/// @param C wave speed
namespace dg_fun{

	const int num_of_equation = 3;

	double C;
//	const double C =  1.0 / (4.0 * std::atan(1.0));
};

/// @brief
/// refinement (hp) refinement swtich.
/// @param adapt refinement switch. 
/// @param refine_frequency time step interval of refinement. 
/// @param fit_point_num The number of points that are used to compute least square fit
/// @param tolerance_min the minimum discretization tolerance . If the estimated error exceeds this, refine. 
/// @param tolerance_max the maximum discretization tolerance. If the estimated error smaller than this, coarsen. 
/// @param load_balaning Repartitioning switch. 
namespace dg_refine{

	bool adapt;

	int refine_frequency;	// every time step refine once

	int fit_point_num;

	double tolerance_min;

	double tolerance_max;

	bool load_balancing;

};

//----------------------------------------------------------------------

/// @brief
/// variable to control the outputs
/// @param
namespace dg_io{

	const bool io = true;

	const int output_frequency = 1;
};


// variables you do not need to change----------------------------------

/// @brief mpi variables
/// @param rank process rank
/// @param num_proc total number of processor
namespace mpi{
	int rank;    
	int num_proc;    
};

//-----------------------------------------------------------------------
