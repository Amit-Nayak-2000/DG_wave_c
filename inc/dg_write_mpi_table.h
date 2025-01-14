#ifndef DG_WRITE_MPI_TABLE_H
#define DG_WRITEMPI_TABLE_H

#include <vector>
#include <unordered_map>
#include "dg_boundary_table.h"

//void Write_mpi_table(std::unordered_map<int, std::vector<mpi_table>>& table);

void Write_table_all(std::unordered_map<int, std::vector<mpi_table>>& table1, 
			std::unordered_map<int, std::vector<mpi_table>>& table2);

#endif
