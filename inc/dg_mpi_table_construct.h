#ifndef DG_MPI_TABLE_CONSTRUCT_H
#define DG_MPI_TABLE_CONSTRUCT_H

#include "dg_boundary_table.h"
#include <vector>

void Construct_mpi_table(std::vector<table_elem>& north, std::vector<table_elem>& south);

#endif