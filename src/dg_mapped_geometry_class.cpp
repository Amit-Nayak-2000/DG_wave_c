#include "dg_mapped_geometry_class.h"
#include "dg_single_index.h"
#include "dg_transfinite_quad_map.h"
#include "dg_transfinite_quad_metrics.h"
#include <cstdlib>
#include <cmath>
#include <iostream>
#include "dg_param.h"

MappedGeometry::MappedGeometry(){

    this->N = grid::nmin;
    
    //Initialize the Mapped Geometry Arrays
    x_node = std::vector<double>((N + 1) * (N + 1));
    y_node = std::vector<double>((N + 1) * (N + 1));

    x_boundary = std::vector<double>((N + 1)*4);
    y_boundary = std::vector<double>((N + 1)*4);

    delx_delxi = std::vector<double>((N + 1) * (N + 1));
    delx_deleta = std::vector<double>((N + 1) * (N + 1));
    dely_delxi = std::vector<double>((N + 1) * (N + 1));
    dely_deleta = std::vector<double>((N + 1) * (N + 1));

    jacobian = std::vector<double>((N + 1) * (N + 1));

    boundary_normal = std::vector<double>((N + 1)*4);
    boundary_normal_y = std::vector<double>((N + 1)*4);
    scaling_factor = std::vector<double>((N + 1)*4);

}

MappedGeometry::MappedGeometry(CurveInterpolant boundaries[4], std::vector<double> gl, int porder){

    this->N = porder;
    std::cout << "New constructor being called" << std::endl;
    
    //Initialize the Mapped Geometry Arrays
    x_node = std::vector<double>((N + 1) * (N + 1));
    y_node = std::vector<double>((N + 1) * (N + 1));

    x_boundary = std::vector<double>((N + 1)*4);
    y_boundary = std::vector<double>((N + 1)*4);

    delx_delxi = std::vector<double>((N + 1) * (N + 1));
    delx_deleta = std::vector<double>((N + 1) * (N + 1));
    dely_delxi = std::vector<double>((N + 1) * (N + 1));
    dely_deleta = std::vector<double>((N + 1) * (N + 1));

    jacobian = std::vector<double>((N + 1) * (N + 1));

    boundary_normal = std::vector<double>((N + 1)*4);
    boundary_normal_y = std::vector<double>((N + 1)*4);
    scaling_factor = std::vector<double>((N + 1)*4);


    for(int j = 0; j <= N; j++){
        for(int i = 0; i <= N; i++){

            int index = Get_single_index(i, j, N+1);

            TransfiniteQuadMap(gl[i], gl[j], boundaries, x_node[index], y_node[index]);
            
            TransfiniteQuadMetrics(gl[i], gl[j], boundaries, delx_delxi[index], delx_deleta[index], dely_delxi[index], dely_deleta[index]);

            jacobian[index] = (delx_delxi[index] * dely_deleta[index]) - (delx_deleta[index] * dely_delxi[index]);
        }
    }


    for(int j = 0; j <= N; j++){
        //right boundary reference square
        int boundaryindex1 = Get_single_index(j, 1, 4);
        int metricindex1 = Get_single_index(N, j, N+1);

        TransfiniteQuadMap(1, gl[j], boundaries, x_boundary[boundaryindex1], y_boundary[boundaryindex1]);

        TransfiniteQuadMetrics(1, gl[j], boundaries, delx_delxi[metricindex1], delx_deleta[metricindex1], dely_delxi[metricindex1], dely_deleta[metricindex1]);

        jacobian[metricindex1] =  (delx_delxi[metricindex1] * dely_deleta[metricindex1]) - (delx_deleta[metricindex1] * dely_delxi[metricindex1]);

        scaling_factor[boundaryindex1] = std::sqrt((dely_deleta[metricindex1] * dely_deleta[metricindex1]) + (delx_deleta[metricindex1] * delx_deleta[metricindex1]));

        boundary_normal[boundaryindex1] = (jacobian[metricindex1] / std::abs(jacobian[metricindex1])) * ((dely_deleta[metricindex1]) / (scaling_factor[boundaryindex1]));
        boundary_normal_y[boundaryindex1] = (jacobian[metricindex1] / std::abs(jacobian[metricindex1])) * (( - delx_deleta[metricindex1]) / (scaling_factor[boundaryindex1]));
        

        //left boundary reference square
        int boundaryindex2 = Get_single_index(j, 3, 4);
        int metricindex2 = Get_single_index(0, j, N+1);

        TransfiniteQuadMap(-1, gl[j], boundaries, x_boundary[boundaryindex2], y_boundary[boundaryindex2]);

        // std::cout << x_boundary[boundaryindex2] << ", " << y_boundary[boundaryindex2] << std::endl;
        
        TransfiniteQuadMetrics(-1, gl[j], boundaries, delx_delxi[metricindex2], delx_deleta[metricindex2], dely_delxi[metricindex2], dely_deleta[metricindex2]);

        jacobian[metricindex2] =  (delx_delxi[metricindex2] * dely_deleta[metricindex2]) - (delx_deleta[metricindex2] * dely_delxi[metricindex2]);

        scaling_factor[boundaryindex2] = std::sqrt((dely_deleta[metricindex2] * dely_deleta[metricindex2]) + (delx_deleta[metricindex2] * delx_deleta[metricindex2]));

        boundary_normal[boundaryindex2] = -(jacobian[metricindex2] / std::abs(jacobian[metricindex2])) * ((dely_deleta[metricindex2]) / (scaling_factor[boundaryindex2]));
        boundary_normal_y[boundaryindex2] = -(jacobian[metricindex2] / std::abs(jacobian[metricindex2])) * ((- delx_deleta[metricindex2]) / (scaling_factor[boundaryindex2]));
        
    }

    //double check this
    for(int i = 0; i <= N; i++){
        //bottom boundary reference square
        int boundaryindex3 = Get_single_index(i, 0, 4);
        int metricindex3 = Get_single_index(i, 0, N+1);

        TransfiniteQuadMap(gl[i], -1, boundaries, x_boundary[boundaryindex3], y_boundary[boundaryindex3]);

        TransfiniteQuadMetrics(gl[i], -1, boundaries, delx_delxi[metricindex3], delx_deleta[metricindex3], dely_delxi[metricindex3], dely_deleta[metricindex3]);

        jacobian[metricindex3] =  (delx_delxi[metricindex3] * dely_deleta[metricindex3]) - (delx_deleta[metricindex3] * dely_delxi[metricindex3]);

        scaling_factor[boundaryindex3] = std::sqrt((dely_delxi[metricindex3] * dely_delxi[metricindex3]) + (delx_delxi[metricindex3] * delx_delxi[metricindex3]));

        boundary_normal[boundaryindex3] = -(jacobian[metricindex3] / std::abs(jacobian[metricindex3])) * ((-dely_delxi[metricindex3]) / (scaling_factor[boundaryindex3]));
        boundary_normal_y[boundaryindex3] = -(jacobian[metricindex3] / std::abs(jacobian[metricindex3])) * (( delx_delxi[metricindex3]) / (scaling_factor[boundaryindex3]));

        //top boundary reference square
        int boundaryindex4 = Get_single_index(i, 2, 4);
        int metricindex4 = Get_single_index(i, N, N+1);

        TransfiniteQuadMap(gl[i], 1, boundaries, x_boundary[boundaryindex4], y_boundary[boundaryindex4]);

        TransfiniteQuadMetrics(gl[i], 1, boundaries, delx_delxi[metricindex4], delx_deleta[metricindex4], dely_delxi[metricindex4], dely_deleta[metricindex4]);

        jacobian[metricindex4] =  (delx_delxi[metricindex4] * dely_deleta[metricindex4]) - (delx_deleta[metricindex4] * dely_delxi[metricindex4]);

        scaling_factor[boundaryindex4] = std::sqrt((dely_delxi[metricindex4] * dely_delxi[metricindex4]) + (delx_delxi[metricindex4] * delx_delxi[metricindex4]));

        boundary_normal[boundaryindex4] = (jacobian[metricindex4] / std::abs(jacobian[metricindex4])) * ((-dely_delxi[metricindex4]) / (scaling_factor[boundaryindex4]));
        boundary_normal_y[boundaryindex4] = (jacobian[metricindex4] / std::abs(jacobian[metricindex4])) * ((delx_delxi[metricindex4]) / (scaling_factor[boundaryindex4]));

        // std::cout << x_boundary[boundaryindex4] << ", " << y_boundary[boundaryindex4] << std::endl;
        
    }

}

void MappedGeometry::updateOrder(int N){
    this->N = N;
    
    //Initialize the Mapped Geometry Arrays
    x_node = std::vector<double>((N + 1) * (N + 1));
    y_node = std::vector<double>((N + 1) * (N + 1));

    x_boundary = std::vector<double>((N + 1)*4);
    y_boundary = std::vector<double>((N + 1)*4);

    delx_delxi = std::vector<double>((N + 1) * (N + 1));
    delx_deleta = std::vector<double>((N + 1) * (N + 1));
    dely_delxi = std::vector<double>((N + 1) * (N + 1));
    dely_deleta = std::vector<double>((N + 1) * (N + 1));

    jacobian = std::vector<double>((N + 1) * (N + 1));

    boundary_normal = std::vector<double>((N + 1)*4);
    boundary_normal_y = std::vector<double>((N + 1)*4);
    scaling_factor = std::vector<double>((N + 1)*4);

}

void MappedGeometry::initMetrics(CurveInterpolant boundaries[4], std::vector<double> gl){

    // std::cout << "being initalized" << std::endl;

    for(int j = 0; j <= N; j++){
        for(int i = 0; i <= N; i++){

            int index = Get_single_index(i, j, N+1);

            TransfiniteQuadMap(gl[i], gl[j], boundaries, x_node[index], y_node[index]);
            
            TransfiniteQuadMetrics(gl[i], gl[j], boundaries, delx_delxi[index], delx_deleta[index], dely_delxi[index], dely_deleta[index]);

            jacobian[index] = (delx_delxi[index] * dely_deleta[index]) - (delx_deleta[index] * dely_delxi[index]);
        }
    }


    for(int j = 0; j <= N; j++){
        //right boundary reference square
        int boundaryindex1 = Get_single_index(j, 1, 4);
        int metricindex1 = Get_single_index(N, j, N+1);

        TransfiniteQuadMap(1, gl[j], boundaries, x_boundary[boundaryindex1], y_boundary[boundaryindex1]);

        TransfiniteQuadMetrics(1, gl[j], boundaries, delx_delxi[metricindex1], delx_deleta[metricindex1], dely_delxi[metricindex1], dely_deleta[metricindex1]);

        jacobian[metricindex1] =  (delx_delxi[metricindex1] * dely_deleta[metricindex1]) - (delx_deleta[metricindex1] * dely_delxi[metricindex1]);

        scaling_factor[boundaryindex1] = std::sqrt((dely_deleta[metricindex1] * dely_deleta[metricindex1]) + (delx_deleta[metricindex1] * delx_deleta[metricindex1]));

        boundary_normal[boundaryindex1] = (jacobian[metricindex1] / std::abs(jacobian[metricindex1])) * ((dely_deleta[metricindex1]) / (scaling_factor[boundaryindex1]));
        boundary_normal_y[boundaryindex1] = (jacobian[metricindex1] / std::abs(jacobian[metricindex1])) * (( - delx_deleta[metricindex1]) / (scaling_factor[boundaryindex1]));

        //left boundary reference square
        int boundaryindex2 = Get_single_index(j, 3, 4);
        int metricindex2 = Get_single_index(0, j, N+1);

        TransfiniteQuadMap(-1, gl[j], boundaries, x_boundary[boundaryindex2], y_boundary[boundaryindex2]);

        // std::cout << x_boundary[boundaryindex2] << ", " << y_boundary[boundaryindex2] << std::endl;
        
        TransfiniteQuadMetrics(-1, gl[j], boundaries, delx_delxi[metricindex2], delx_deleta[metricindex2], dely_delxi[metricindex2], dely_deleta[metricindex2]);

        jacobian[metricindex2] =  (delx_delxi[metricindex2] * dely_deleta[metricindex2]) - (delx_deleta[metricindex2] * dely_delxi[metricindex2]);

        scaling_factor[boundaryindex2] = std::sqrt((dely_deleta[metricindex2] * dely_deleta[metricindex2]) + (delx_deleta[metricindex2] * delx_deleta[metricindex2]));

        boundary_normal[boundaryindex2] = -(jacobian[metricindex2] / std::abs(jacobian[metricindex2])) * ((dely_deleta[metricindex2]) / (scaling_factor[boundaryindex2]));
        boundary_normal_y[boundaryindex2] = -(jacobian[metricindex2] / std::abs(jacobian[metricindex2])) * ((- delx_deleta[metricindex2]) / (scaling_factor[boundaryindex2]));
        
    }

    //double check this
    for(int i = 0; i <= N; i++){
        //bottom boundary reference square
        int boundaryindex3 = Get_single_index(i, 0, 4);
        int metricindex3 = Get_single_index(i, 0, N+1);

        TransfiniteQuadMap(gl[i], -1, boundaries, x_boundary[boundaryindex3], y_boundary[boundaryindex3]);

        TransfiniteQuadMetrics(gl[i], -1, boundaries, delx_delxi[metricindex3], delx_deleta[metricindex3], dely_delxi[metricindex3], dely_deleta[metricindex3]);

        jacobian[metricindex3] =  (delx_delxi[metricindex3] * dely_deleta[metricindex3]) - (delx_deleta[metricindex3] * dely_delxi[metricindex3]);

        scaling_factor[boundaryindex3] = std::sqrt((dely_delxi[metricindex3] * dely_delxi[metricindex3]) + (delx_delxi[metricindex3] * delx_delxi[metricindex3]));

        boundary_normal[boundaryindex3] = -(jacobian[metricindex3] / std::abs(jacobian[metricindex3])) * ((-dely_delxi[metricindex3]) / (scaling_factor[boundaryindex3]));
        boundary_normal_y[boundaryindex3] = -(jacobian[metricindex3] / std::abs(jacobian[metricindex3])) * (( delx_delxi[metricindex3]) / (scaling_factor[boundaryindex3]));

        //top boundary reference square
        int boundaryindex4 = Get_single_index(i, 2, 4);
        int metricindex4 = Get_single_index(i, N, N+1);

        TransfiniteQuadMap(gl[i], 1, boundaries, x_boundary[boundaryindex4], y_boundary[boundaryindex4]);

        TransfiniteQuadMetrics(gl[i], 1, boundaries, delx_delxi[metricindex4], delx_deleta[metricindex4], dely_delxi[metricindex4], dely_deleta[metricindex4]);

        jacobian[metricindex4] =  (delx_delxi[metricindex4] * dely_deleta[metricindex4]) - (delx_deleta[metricindex4] * dely_delxi[metricindex4]);

        scaling_factor[boundaryindex4] = std::sqrt((dely_delxi[metricindex4] * dely_delxi[metricindex4]) + (delx_delxi[metricindex4] * delx_delxi[metricindex4]));

        boundary_normal[boundaryindex4] = (jacobian[metricindex4] / std::abs(jacobian[metricindex4])) * ((-dely_delxi[metricindex4]) / (scaling_factor[boundaryindex4]));
        boundary_normal_y[boundaryindex4] = (jacobian[metricindex4] / std::abs(jacobian[metricindex4])) * ((delx_delxi[metricindex4]) / (scaling_factor[boundaryindex4]));

        // std::cout << x_boundary[boundaryindex4] << ", " << y_boundary[boundaryindex4] << std::endl;
        
    }

}

void MappedGeometry::recalcscalenorms(CurveInterpolant boundaries[4], std::vector<double> gl){
    // std::cout << "being initalized" << std::endl;

    for(int j = 0; j <= N; j++){
        for(int i = 0; i <= N; i++){

            int index = Get_single_index(i, j, N+1);
            
            TransfiniteQuadMetrics(gl[i], gl[j], boundaries, delx_delxi[index], delx_deleta[index], dely_delxi[index], dely_deleta[index]);

            jacobian[index] = (delx_delxi[index] * dely_deleta[index]) - (delx_deleta[index] * dely_delxi[index]);
        }
    }


    for(int j = 0; j <= N; j++){
        //right boundary reference square
        int boundaryindex1 = Get_single_index(j, 1, 4);
        int metricindex1 = Get_single_index(N, j, N+1);

        TransfiniteQuadMetrics(1, gl[j], boundaries, delx_delxi[metricindex1], delx_deleta[metricindex1], dely_delxi[metricindex1], dely_deleta[metricindex1]);

        jacobian[metricindex1] =  (delx_delxi[metricindex1] * dely_deleta[metricindex1]) - (delx_deleta[metricindex1] * dely_delxi[metricindex1]);

        scaling_factor[boundaryindex1] = std::sqrt((dely_deleta[metricindex1] * dely_deleta[metricindex1]) + (delx_deleta[metricindex1] * delx_deleta[metricindex1]));

        boundary_normal[boundaryindex1] = (jacobian[metricindex1] / std::abs(jacobian[metricindex1])) * ((dely_deleta[metricindex1]) / (scaling_factor[boundaryindex1]));
        boundary_normal_y[boundaryindex1] = (jacobian[metricindex1] / std::abs(jacobian[metricindex1])) * (( - delx_deleta[metricindex1]) / (scaling_factor[boundaryindex1]));

        //left boundary reference square
        int boundaryindex2 = Get_single_index(j, 3, 4);
        int metricindex2 = Get_single_index(0, j, N+1);
        
        TransfiniteQuadMetrics(-1, gl[j], boundaries, delx_delxi[metricindex2], delx_deleta[metricindex2], dely_delxi[metricindex2], dely_deleta[metricindex2]);

        jacobian[metricindex2] =  (delx_delxi[metricindex2] * dely_deleta[metricindex2]) - (delx_deleta[metricindex2] * dely_delxi[metricindex2]);

        scaling_factor[boundaryindex2] = std::sqrt((dely_deleta[metricindex2] * dely_deleta[metricindex2]) + (delx_deleta[metricindex2] * delx_deleta[metricindex2]));

        boundary_normal[boundaryindex2] = -(jacobian[metricindex2] / std::abs(jacobian[metricindex2])) * ((dely_deleta[metricindex2]) / (scaling_factor[boundaryindex2]));
        boundary_normal_y[boundaryindex2] = -(jacobian[metricindex2] / std::abs(jacobian[metricindex2])) * ((- delx_deleta[metricindex2]) / (scaling_factor[boundaryindex2]));
        
    }

    //double check this
    for(int i = 0; i <= N; i++){
        //bottom boundary reference square
        int boundaryindex3 = Get_single_index(i, 0, 4);
        int metricindex3 = Get_single_index(i, 0, N+1);

        TransfiniteQuadMetrics(gl[i], -1, boundaries, delx_delxi[metricindex3], delx_deleta[metricindex3], dely_delxi[metricindex3], dely_deleta[metricindex3]);

        jacobian[metricindex3] =  (delx_delxi[metricindex3] * dely_deleta[metricindex3]) - (delx_deleta[metricindex3] * dely_delxi[metricindex3]);

        scaling_factor[boundaryindex3] = std::sqrt((dely_delxi[metricindex3] * dely_delxi[metricindex3]) + (delx_delxi[metricindex3] * delx_delxi[metricindex3]));

        boundary_normal[boundaryindex3] = -(jacobian[metricindex3] / std::abs(jacobian[metricindex3])) * ((-dely_delxi[metricindex3]) / (scaling_factor[boundaryindex3]));
        boundary_normal_y[boundaryindex3] = -(jacobian[metricindex3] / std::abs(jacobian[metricindex3])) * (( delx_delxi[metricindex3]) / (scaling_factor[boundaryindex3]));

        //top boundary reference square
        int boundaryindex4 = Get_single_index(i, 2, 4);
        int metricindex4 = Get_single_index(i, N, N+1);

        TransfiniteQuadMetrics(gl[i], 1, boundaries, delx_delxi[metricindex4], delx_deleta[metricindex4], dely_delxi[metricindex4], dely_deleta[metricindex4]);

        jacobian[metricindex4] =  (delx_delxi[metricindex4] * dely_deleta[metricindex4]) - (delx_deleta[metricindex4] * dely_delxi[metricindex4]);

        scaling_factor[boundaryindex4] = std::sqrt((dely_delxi[metricindex4] * dely_delxi[metricindex4]) + (delx_delxi[metricindex4] * delx_delxi[metricindex4]));

        boundary_normal[boundaryindex4] = (jacobian[metricindex4] / std::abs(jacobian[metricindex4])) * ((-dely_delxi[metricindex4]) / (scaling_factor[boundaryindex4]));
        boundary_normal_y[boundaryindex4] = (jacobian[metricindex4] / std::abs(jacobian[metricindex4])) * ((delx_delxi[metricindex4]) / (scaling_factor[boundaryindex4]));

        
    }
}

void MappedGeometry::outputMetrics(){
    std::cout << "Nodes" << std::endl;
		for(int j = 0; j <= this->N; ++j){
			for(int i = 0; i <= this->N; ++i){
				int num_p = Get_single_index(i, j, this->N + 1);
				std::cout << "{" << x_node[num_p] << ", " << y_node[num_p] << "} ";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;

		std::cout << "del(x,y)xi" << std::endl;
		for(int j = 0; j <= this->N; ++j){
			for(int i = 0; i <= this->N; ++i){
				int num_p = Get_single_index(i, j, this->N + 1);
				std::cout << "{" << delx_delxi[num_p] << ", " << dely_delxi[num_p] << "} ";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;

		std::cout << "del(x,y)eta" << std::endl;
		for(int j = 0; j <= this->N; ++j){
			for(int i = 0; i <= this->N; ++i){
				int num_p = Get_single_index(i, j, this->N + 1);
				std::cout << "{" << delx_deleta[num_p] << ", " << dely_deleta[num_p] << "} ";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;

		std::cout << "Jacobian" << std::endl;
		for(int j = 0; j <= this->N; ++j){
			for(int i = 0; i <= this->N; ++i){
				int num_p = Get_single_index(i, j, this->N + 1);
				std::cout << "{" << jacobian[num_p] << "} ";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;

		std::cout << "Boundaries" << std::endl;
		for(int j = 0; j <= 3; j++){
			for(int i = 0; i <= this->N; i++){
				int num_p = Get_single_index(i, j, 4);
				std::cout << "{" << x_boundary[num_p] << ", " << y_boundary[num_p] << "}" << std::endl;
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;

		std::cout << "Normals" << std::endl;
		for(int j = 0; j <= 3; j++){
			for(int i = 0; i <= this->N; i++){
				int num_p = Get_single_index(i, j, 4);
				std::cout << "{" << boundary_normal[num_p] << ", " << boundary_normal_y[num_p] << "}" << std::endl;
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;

		std::cout << "Scaling Factor" << std::endl;
		for(int j = 0; j <= 3; j++){
			for(int i = 0; i <= this->N; i++){
				int num_p = Get_single_index(i, j, 4);
				std::cout << "{" << scaling_factor[num_p] << "}" << std::endl;
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
}
