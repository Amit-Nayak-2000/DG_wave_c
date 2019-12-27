#include <mpi.h>
#include <cmath>

/// @brief Generate the dual coordinate of the current element.
/// @param x element k's node 1 x coordinate
/// @param y element k's node 1 y coordinate
/// @param delta_x current element size
/// @param delta_y current element size
/// @param coord return current element dual (integer) coordinate (i, j)
void Get_dual_coord_2d(double& x, double& y, double& delta_x, double& delta_y, int* coord){

	*coord = std::round(x / delta_x); 
	*(coord+1) = std::round(y / delta_y); 


}
