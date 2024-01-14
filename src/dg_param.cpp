#include <string>
#include "dg_param.h"
#include <cmath>

//variable you must change--------------------------------------------

/// @brief mesh file information
/// @param fileplace The path of mesh file and mesh file name
/// @param output_place output directory. 
namespace fileinfo{
	const std::string fileplace = "../gmsh_files/256_elements.msh";
	// const std::string fileplace = "../gmsh_files/1element.msh";
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
	const int exp_x = 4; 
	const int exp_y = 4; 
	
	const double gx_l = 0.0;
	const double gx_r = 1.0; 
	const double gy_l = 0.0;
	const double gy_r = 1.0; 

	const int nmin = 8;	
	const int nmax = 12;

	const int hlevel_max = 2;	
};
//---------------------------------------------------------------------


// variables you could change------------------------------------------

/// @brief
/// Time related variables
/// @param t_total total time integal
/// @param nt time step number
namespace dg_time{

	//wave
	// const double t_total = 1.0e-2* 5;
	// const int nt = 5000;
	// const double t_total = 0.5;
	// const int nt = 50000;


	//curved channel
	// const double t_total = 200*1.0e-2;
	// const int nt = 200000;

	// //acoustic scatter time
	const double t_total = 3;
	const int nt = 60000;
	// const double t_total = 1.0e-3;
	// const int nt = 20;

};

/// @brief
/// Function related parameters
/// @param num_of_equation  number of equations
/// @param C wave speed
namespace dg_fun{

	const int num_of_equation = 3;

	const double C = 1.0;
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

	const bool adapt = false;

	const int refine_frequency = 1;	// every time step refine once

	const int fit_point_num = 4;

	const double tolerance_min = 1.0e-4; // originally e-6

	const double tolerance_max = 1.0e-9;  // originally e-10

	const bool load_balancing = false;

};

//e-4, e-8 for curved annulus? mabye e-6?
//----------------------------------------------------------------------

/// @brief
/// variable to control the outputs
/// @param
namespace dg_io{

	const bool io = true;

	const int output_frequency = 500;
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
