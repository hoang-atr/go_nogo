//
//  utils.hpp
//  cerebellum
//
//  Created by hoang on 2024/03/13.
//

#ifndef utils_hpp
#define utils_hpp

#include <stdio.h>

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>

#define MAX_LICK_COUNT 6.0
#define NEGATIVE_ERROR_START_RATE 7.0
#define NEGATIVE_ERROR_END_RATE 2.0
#define POSITIVE_ERROR_START_RATE 3.0
#define POSITIVE_ERROR_END_RATE 6.0
#define BASE_ERROR_RATE 2.0

using namespace std;

int read_integers_from_file(std::string filename, vector<int> &numbers);
int write_floats_to_file(std::string filename, vector<float> lick_rate);

int set_output_file_name(int run, string &fnIO1, string &fnIO2, string &fnPC1, string &fnPC2, string &fnPfPC1, string &fnPfPC2, string &fnLick, string &fnCN1, string &fnCN2, string &fnIOcoupling1, string &fnIOcoupling2);

float spike_rate_to_lick_rate(float rate);

float error_lick_to_poisson_rate(float error_lick);
float mean_omitnan(std::vector< std::vector< float > > vector);

#endif /* utils_hpp */
