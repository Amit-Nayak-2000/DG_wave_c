#include <mpi.h>
#include "dg_param.h"
#include "dg_nodal_2d_storage.h"
#include "dg_constructor.h"
/// @brief
/// Driver for DG approxiation. Algorithm 51. 
/// First, get DG basis parameters, such as collocation points and weights.
/// Then marches by each time step. Using explicit 3rd order Runge-Kutta methods.
void Driver_for_DG_approximation(){

	// construct basis
	Construct_basis();

	// time step
	double delta_t = time::t_total / time::nt;

	// current time
	double tn{};

	



}
