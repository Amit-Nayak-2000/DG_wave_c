#include "dg_user_defined.h"
#include "dg_param.h"
#include "dg_read_params.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "exprtk.hpp"

void readParameters(double values[24], std::string &filelocation);
void initializeParameters();

void initializeParameters(){
    double values[24];
    std::string filelocation;
    // int BCval;
    readParameters(values, filelocation);

    //dg_user_defined.cpp
    user::kx = values[0];
    user::ky = values[1];
    user::D = values[2];
    user::xx0 = values[3];
    user::yy0 = values[4];
    user::BC = values[23];
    // user::BC = BCval;

    //dg_param.cpp
    //namespace grid:
    grid::exp_x = values[5];
    grid::exp_y = values[6];
    grid::gx_l = values[7];
    grid::gx_r = values[8];
    grid::gy_l = values[9];
    grid::gy_r = values[10];
    grid::nmin = values[11];
    grid::nmax = values[12];
    grid::hlevel_max = values[13];
    //namespace dg_time:
    dg_time::t_total = values[14];
    dg_time::nt = values[15];
    //namespace dg_fun:
    dg_fun::C = values[16];
    //namespace dg_refine
    if(values[17] == 1){
        dg_refine::adapt = true;
    }
    else{
        dg_refine::adapt = false;
    }
    dg_refine::refine_frequency = values[18];
    dg_refine::fit_point_num = values[19];
    dg_refine::tolerance_min = values[20];
    dg_refine::tolerance_max = values[21];

    if(values[22] == 1){
        dg_refine::load_balancing = true;
    }
    else{
        dg_refine::load_balancing = false;
    }

    //namespace fileinfo
    fileinfo::fileplace = filelocation;



    
}

void readParameters(double values[24], std::string &filelocation){
    //variables used for parsing literals from string
    std::string parse;
    std::string delimiter = ": ";
    size_t pos = 0;

    //variables used for parsing expressions
    typedef exprtk::expression<double>   expression_t;
    typedef exprtk::parser<double>       parser_t;

    expression_t expression;
    parser_t parser;

    std::ifstream MyFile("userparameters.txt");

    for(int i = 0; i < 24; i++){
        getline(MyFile, parse);
        pos = parse.find(delimiter);
        parse = parse.substr(pos+2);

        if (!parser.compile(parse,expression)){
            printf("Compilation error...\n");
            return;
        }

        values[i] = expression.value();

    }

    getline(MyFile, parse);
    pos = parse.find(delimiter);
    parse = parse.substr(pos+2);
    filelocation = parse;

    // getline(MyFile, parse);
    // pos = parse.find(delimiter);
    // parse = parse.substr(pos+2);
    // if (!parser.compile(parse,expression)){
    //         printf("Compilation error...\n");
    //         return;
    // }

    // BCval = expression.value();

    // std::cout << "BCval is: " << BCval << std::endl;

    MyFile.close();

}

