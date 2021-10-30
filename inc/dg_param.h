#ifndef DG_PARAM_H
#define DG_PARAM_H

#include <string>

// mesh file----------------------------------------
namespace fileinfo{
	extern  std::string fileplace;
	extern const std::string output_place;
	extern const std::string eff_filename;
	extern const std::string crosssection_filename;
	extern const std::string exact_error_filename;
}
// -------------------------------------------------

// Grid size ---------------------------------------
namespace grid{
	extern  int exp_x;
	extern  int exp_y;

	extern  double gx_l;
	extern  double gx_r;
	extern  double gy_l;
	extern  double gy_r;
	
	extern  int nmin;
	extern  int nmax;

	extern  int hlevel_max;
};
//--------------------------------------------------

// time---------------------------------------------
namespace dg_time{

	extern  double t_total;
	extern  int nt;

};
//--------------------------------------------------


// dg function related parameters--------------------
namespace dg_fun{

	extern const int num_of_equation;
	
	extern  double C;
};

//---------------------------------------------------

// refinemnt ----------------------------------------
namespace dg_refine{

	extern  bool adapt;

	extern  int refine_frequency;

	extern  int fit_point_num;
	
	extern  double tolerance_min;

	extern  double tolerance_max;

	extern  bool load_balancing;
};
//---------------------------------------------------

// output variables --------------------------------
namespace dg_io{

	extern const bool io;

	extern const int output_frequency;
};
//--------------------------------------------------

// mpi variables-------------------------------------
namespace mpi{
	extern int rank;
	extern int num_proc;
};
//---------------------------------------------------


#endif
