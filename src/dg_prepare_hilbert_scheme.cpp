#include <mpi.h>
#include "dg_param.h"
#include "dg_read_mesh_2d.h"
#include "dg_gen_dual_graph.h"
#include "dg_hilbert_sort.h"
#include "dg_prepare_hilbert_scheme.h"
#include "dg_gen_status.h"

/// @brief
/// 1. Read the .msh file generated by GMSH. 
/// 2. Sort the local element node numbering sequence. 
/// 3. Generate dual graph of the mesh. 
/// 4. Re-numbering all the elements by using Hilbert curve.
void Hilbert_numbering(){
	
	if(mpi::rank == 0){
		
		// read .msh file
		Read_mesh_2d();		
		
		// generate dual graph
		Gen_dual_graph_2d();
		
		// renumbering elements by hilbert curve
		Hilbert_sort_2d();
		
		// generate element status
		Gen_status_all();
	}


}
