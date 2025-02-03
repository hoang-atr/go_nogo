//
//  utils.cpp
//  cerebellum
//
//  Created by hoang on 2024/03/13.
//
#include <vector>
#include "utils.hpp"

float mean_omitnan(std::vector< std::vector< float > > vector) {
    int count = 0;
    float sum = 0;
    for (int i=0; i<vector.size(); i++) {
        for (int j=0; j<vector[i].size(); j++) {
            if (!std::isnan(vector[i][j])) {
                count = count + 1;
                sum = sum + vector[i][j];
            }
        }
    }
    if (count>0) {
        return (float)sum/count;
    } else {
        return -1;
    }
}

int read_integers_from_file(std::string filename, vector<int> &numbers) {  
  int number;

  ifstream input_file(filename);

  if (!input_file.is_open()) {
    cerr << "Could not open the file - '" << filename << "'" << endl;
    return 0;
  }

  while (input_file >> number) {
    numbers.push_back(number);
  }

  input_file.close();

  return (int)numbers.size();
}

int write_floats_to_file(std::string filename, vector<float> lick_rate) {
    
    ofstream output_file(filename, std::ios::binary);
    
    if (!lick_rate.empty())
        output_file.write(reinterpret_cast<char*>(&lick_rate[0]), lick_rate.size() * sizeof(lick_rate[0]));
    
    output_file.close();
    
    return 1;
}

int set_output_file_name(int run, string &fnIO1, string &fnIO2, string &fnPC1, string &fnPC2, string &fnPfPC1, string &fnPfPC2, string &fnLick, string &fnCN1, string &fnCN2, string &fnIOcoupling1, string &fnIOcoupling2) {
    std::ostringstream oss;
    
    oss << "results/spk_IO1_r" << run << ".dat";
    fnIO1 = oss.str();
    
    oss.str("");
    oss << "results/spk_PC1_r" << run << ".dat";
    fnPC1 = oss.str();
    
    oss.str("");
    oss << "results/wt_PfPC1_r" << run << ".dat";
    fnPfPC1 = oss.str();
    
    oss.str("");
    oss << "results/spk_IO2_r" << run << ".dat";
    fnIO2 = oss.str();
    
    oss.str("");
    oss << "results/spk_PC2_r" << run << ".dat";
    fnPC2 = oss.str();
    
    oss.str("");
    oss << "results/wt_PfPC2_r" << run << ".dat";
    fnPfPC2 = oss.str();
    
    oss.str("");
    oss << "results/spk_lick_r" << run << ".dat";
    fnLick = oss.str();
    
    oss.str("");
    oss << "results/spk_CN1_r" << run << ".dat";
    fnCN1 = oss.str();
    
    oss.str("");
    oss << "results/spk_CN2_r" << run << ".dat";
    fnCN2 = oss.str();
    
    oss.str("");
    oss << "results/IOcoupling1_r" << run << ".dat";
    fnIOcoupling1 = oss.str();
    
    oss.str("");
    oss << "results/IOcoupling2_r" << run << ".dat";
    fnIOcoupling2 = oss.str();
    
    return run;
}

float spike_rate_to_lick_rate(float rate) {
    float const a = 0.3f;
    float const b = 16.0f;
    
    return  (MAX_LICK_COUNT - MAX_LICK_COUNT / (1+std::exp(-a*rate+b))); // sigmoid function
}

float error_lick_to_poisson_rate(float error_lick) {
    
    float rate = BASE_ERROR_RATE;
    float slope = 0;
    if (error_lick < 0) {
        slope = -(NEGATIVE_ERROR_START_RATE-NEGATIVE_ERROR_END_RATE)/MAX_LICK_COUNT;
        rate = (float) slope*error_lick + NEGATIVE_ERROR_END_RATE;
    } else if (error_lick > 0) {
        slope = (POSITIVE_ERROR_START_RATE-POSITIVE_ERROR_END_RATE)/(MAX_LICK_COUNT/2);
        rate = (float) slope*error_lick + POSITIVE_ERROR_END_RATE;
        
        if (rate<POSITIVE_ERROR_START_RATE) {
            rate = POSITIVE_ERROR_START_RATE;
        }
    }
    
    return rate;
}


