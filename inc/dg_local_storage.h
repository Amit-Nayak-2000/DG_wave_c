#ifndef DG_LOCAL_STORAGE_H
#define DG_LOCAL_STORAGE_H

#include "dg_unit.h"
#include <unordered_map>
#include <vector>
#include "dg_boundary_table.h" 
#include "dg_load_struct.h"

namespace local{
 
	extern int local_elem_num;
 
        extern int original_elem_num;

        extern int* elem_range;
	extern int* rank_indicator;
	
	extern double* x_local; 
	extern double* y_local; 
	
	extern char* status;

	extern Unit* head;
	
	extern std::unordered_map<int, Unit*> Hash_elem;

	extern std::vector<std::vector<double>> solution_int_l;
	extern std::vector<std::vector<double>> solution_int_r;
	
};


namespace hrefinement{

	
	extern std::vector<table_elem> north;
	extern std::vector<table_elem> south;
	extern std::vector<table_elem> west;
	extern std::vector<table_elem> east;

	extern std::vector<accum_elem> south_accum;
	extern std::vector<accum_elem> north_accum;
	extern std::vector<accum_elem> west_accum;
	extern std::vector<accum_elem> east_accum;
};

namespace LB{

	extern std::vector<pmap> proc_mapping_table;

	extern int elem_accum;

	extern struct sending_envelope Send;
	
	extern Unit* end;
	extern Unit* my_rank_last;
};

namespace result{

	extern double* exact;
	
	extern double* error;

	extern double* L2_norm;

};

#endif
