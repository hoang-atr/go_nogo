// include CARLsim user interface
#include <carlsim.h>
#include <periodic_spikegen.h>
#include <string>
#include <sstream>
#include <iostream>
#include <time.h>

#include "utils.hpp"

#define WT_ERR_IO 70.0
#define WT_IO_IO 10.0
#define WT_IO_PC 50.0

#define WT_PF_PC 10.0
#define WT_PF_PC_MIN 0.0
#define WT_PF_PC_MAX 25.0

#define WT_PC_CN 50.0

#define RATE_MF 4.0
#define WT_MF_GC 10.0
#define PERC_MF_GC 0.05

#define WT_MF_CN 10.0
#define WT_CN_LICK 10.0

#define BASELINE_RATE_IO_LTP 8
#define BASELINE_RATE_IO_LTD 1

#define INITIAL_ERROR_RATE 6

#define DELAY_MF_GC 1.0
#define DELAY_MF_CN 1.0
#define DELAY_ERR_IO 1.0
#define DELAY_IO_IO 1.0
#define DELAY_PF_PC 1.0
#define DELAY_PC_CN 1.0
#define DELAY_CN_LICK 1.0

#define CHANGE_RATE_IO1 0.8
#define CHANGE_RATE_IO2 0.5
#define INIT_RATE_CN 26
#define INIT_COUPLING 4.7

#define GO 1
#define NO_GO 2

#define ALPHA 0.5

#define LEARNING_RATE_LTP 0.01
#define LEARNING_RATE_LTD 0.05

#define NUM_GC 2000
#define NUM_MF 100

int main(int argc, char *argv[]) {
    std::vector <int> trials;
    int nTrials = read_integers_from_file(argv[1], trials);
    int nRun = atoi(argv[2]);
//    int nTrials = read_integers_from_file("trial.txt", trials);
//    int nRun = 1;
    
    std::cout << "Number of runs: " << nRun << " - Number of trials: " << nTrials << std::endl;
        
    int const trial_duration = 500;     // duration of a trial in ms
    
    for (int run = 1; run<=nRun; run++) {
        long t1 = clock();
        std::cout << "=== Run " << run << ": ";
        
        std::string fnIO1, fnIO2;
        std::string fnPC1, fnPC2;
        std::string fnPfPC1, fnPfPC2;
        std::string fnLick;
        std::string fnCN1, fnCN2;
        std::string fnIOcoupling1, fnIOcoupling2;
        
        set_output_file_name(run, fnIO1, fnIO2, fnPC1, fnPC2, fnPfPC1, fnPfPC2, fnLick, fnCN1, fnCN2, fnIOcoupling1, fnIOcoupling2);
        
        // ---------------- CONFIG STATE -------------------
        // create a network on CPU
        CARLsim sim("cerebellum", CPU_MODE, SILENT);
        
        // configure the network
        Grid3D grid_structure(10,10,1); // post is on a 3x3 grid
        
        // Mossy fibers
        int gMF1 = sim.createSpikeGeneratorGroup("MF1", NUM_MF, EXCITATORY_NEURON);
        int gMF2 = sim.createSpikeGeneratorGroup("MF2", NUM_MF, EXCITATORY_NEURON);
        PoissonRate stiRate1(NUM_MF);
        PoissonRate stiRate2(NUM_MF);
        
        // Granule cells
        int gGC1 = sim.createGroup("GC1", NUM_GC, EXCITATORY_NEURON);
        sim.setNeuronParameters(gGC1, 0.02f, 0.2f, -65.0f, 8.0f);
        
        int gGC2 = sim.createGroup("GC2", NUM_GC, EXCITATORY_NEURON);
        sim.setNeuronParameters(gGC2, 0.02f, 0.2f, -65.0f, 8.0f);
                    
        // Error neurons
        int gErr1 = sim.createSpikeGeneratorGroup("error1", grid_structure, EXCITATORY_NEURON);
        int gErr2 = sim.createSpikeGeneratorGroup("error2", grid_structure, EXCITATORY_NEURON);
        PoissonRate errRate1(grid_structure.N);
        PoissonRate errRate2(grid_structure.N);
        
        // Purkinje cells
        Grid3D gridPC(10,10,1);
        int gPC1=sim.createGroupLIF("PC1", grid_structure, INHIBITORY_NEURON);
        sim.setNeuronParametersLIF(gPC1, 10, 5, -67.0f, -70.0f, RangeRmem(0.5f, 1.0f));
        
        int gPC2=sim.createGroupLIF("PC2", grid_structure, INHIBITORY_NEURON);
        sim.setNeuronParametersLIF(gPC2, 10, 5, -67.0f, -70.0f, RangeRmem(0.5f, 1.0f));
                
        // Inferior olive
        int gIO1=sim.createGroupLIF("IO1", grid_structure, EXCITATORY_NEURON);
        sim.setNeuronParametersLIF(gIO1, 10, 50, -65.0f, -70.0f, RangeRmem(0.5f, 1.0f));
        
        int gIO2=sim.createGroupLIF("IO2", grid_structure, EXCITATORY_NEURON);
        sim.setNeuronParametersLIF(gIO2, 10, 50, -65.0f, -70.0f, RangeRmem(0.5f, 1.0f));
        
        // DCN
        int gCN1=sim.createGroupLIF("CN1", grid_structure, INHIBITORY_NEURON);
        sim.setNeuronParametersLIF(gCN1, 10, 5, -67.0f, -70.0f, RangeRmem(0.5f, 1.0f));
        
        int gCN2=sim.createGroupLIF("CN2", grid_structure, INHIBITORY_NEURON);
        sim.setNeuronParametersLIF(gCN2, 10, 5, -67.0f, -70.0f, RangeRmem(0.5f, 1.0f));
        
        
        // Connections
        // group #1
        int cnnMF_GC1 = sim.connect(gMF1, gGC1, "random", RangeWeight(WT_MF_GC), PERC_MF_GC,  RangeDelay(DELAY_MF_GC), RadiusRF(-1.0), SYN_FIXED);
        int cnnErr_IO1 = sim.connect(gErr1, gIO1, "one-to-one", RangeWeight(WT_ERR_IO), 1.0f, RangeDelay(DELAY_ERR_IO));
        int cnnIO1 = sim.connect(gIO1, gIO1, "gaussian", RangeWeight(WT_IO_IO), 1.0f, RangeDelay(DELAY_IO_IO), RadiusRF(5,5,0));
        int cnnGC1_PC1 = sim.connect(gGC1, gPC1, "full", RangeWeight(WT_PF_PC_MIN,WT_PF_PC,WT_PF_PC_MAX), 1.0,  RangeDelay(DELAY_PF_PC), RadiusRF(-1.0), SYN_PLASTIC);
        int cnnGC2_PC1 = sim.connect(gGC2, gPC1, "full", RangeWeight(WT_PF_PC_MIN,WT_PF_PC,WT_PF_PC_MAX), 1.0,  RangeDelay(DELAY_PF_PC), RadiusRF(-1.0), SYN_PLASTIC);
        int cnnIO_PC1 = sim.connect(gIO1, gPC1, "one-to-one", RangeWeight(WT_IO_PC), 1.0f);
        
        int cnnPC_CN1 = sim.connect(gPC1, gCN1, "one-to-one", RangeWeight(WT_PC_CN), 1.0f, RangeDelay(DELAY_PC_CN));
        
        int cnnMF1_CN1 = sim.connect(gMF1, gCN1, "full", RangeWeight(WT_MF_CN), 1.0f,  RangeDelay(DELAY_MF_CN));
        int cnnMF2_CN1 = sim.connect(gMF2, gCN1, "full", RangeWeight(WT_MF_CN), 1.0f,  RangeDelay(DELAY_MF_CN));
                
        // group #2
        int cnnMF_GC2 = sim.connect(gMF2, gGC2, "random", RangeWeight(WT_MF_GC), PERC_MF_GC,  RangeDelay(DELAY_MF_GC), RadiusRF(-1.0), SYN_FIXED);
        int cnnErr_IO2 = sim.connect(gErr2, gIO2, "one-to-one", RangeWeight(WT_ERR_IO), 1.0f, RangeDelay(DELAY_ERR_IO));
        int cnnIO2 = sim.connect(gIO2, gIO2, "gaussian", RangeWeight(WT_IO_IO), 1.0f, RangeDelay(DELAY_IO_IO), RadiusRF(5,5,0));
        int cnnGC1_PC2 = sim.connect(gGC1, gPC2, "full", RangeWeight(WT_PF_PC_MIN,WT_PF_PC,WT_PF_PC_MAX), 1.0,  RangeDelay(DELAY_PF_PC), RadiusRF(-1.0), SYN_PLASTIC);
        int cnnGC2_PC2 = sim.connect(gGC2, gPC2, "full", RangeWeight(WT_PF_PC_MIN,WT_PF_PC,WT_PF_PC_MAX), 1.0,  RangeDelay(DELAY_PF_PC), RadiusRF(-1.0), SYN_PLASTIC);
        int cnnIO_PC2 = sim.connect(gIO2, gPC2, "one-to-one", RangeWeight(WT_IO_PC), 1.0f);
        
        int cnnPC_CN2 = sim.connect(gPC2, gCN2, "one-to-one", RangeWeight(WT_PC_CN), 1.0f, RangeDelay(DELAY_PC_CN));
        
        int cnnMF1_CN2 = sim.connect(gMF1, gCN2, "full", RangeWeight(WT_MF_CN), 1.0f,  RangeDelay(DELAY_MF_CN));
        int cnnMF2_CN2 = sim.connect(gMF2, gCN2, "full", RangeWeight(WT_MF_CN), 1.0f,  RangeDelay(DELAY_MF_CN));
        
        sim.setConductances(false);
        
        // ---------------- SETUP STATE -------------------
        // build the network
        sim.setupNetwork();
        
        // set some monitors
        SpikeMonitor* smIO1 = sim.setSpikeMonitor(gIO1,fnIO1);
        SpikeMonitor* smPC1 = sim.setSpikeMonitor(gPC1,fnPC1);
        SpikeMonitor* smCN1 = sim.setSpikeMonitor(gCN1,fnCN1);
        
        SpikeMonitor* smIO2 = sim.setSpikeMonitor(gIO2,fnIO2);
        SpikeMonitor* smPC2 = sim.setSpikeMonitor(gPC2,fnPC2);
        SpikeMonitor* smCN2 = sim.setSpikeMonitor(gCN2,fnCN2);
                
        ConnectionMonitor *cmGC1_PC1 = sim.setConnectionMonitor(gGC1, gPC1, fnPfPC1);
        ConnectionMonitor *cmGC2_PC2 = sim.setConnectionMonitor(gGC2, gPC2, fnPfPC2);
        
        ConnectionMonitor *cmCoupling1 = sim.setConnectionMonitor(gIO1, gIO1, fnIOcoupling1);
        ConnectionMonitor *cmCoupling2 = sim.setConnectionMonitor(gIO2, gIO2, fnIOcoupling2);
        
        float curr_err_rate_1 = INITIAL_ERROR_RATE;
        float curr_err_rate_2 = INITIAL_ERROR_RATE;
            
        std::vector<float> Lick;
                
        float prev_spkrateCN1_Go = INIT_RATE_CN;
        float prev_spkrateCN2_Go = INIT_RATE_CN;
        float prev_spkrateCN1_NoGo = INIT_RATE_CN;
        float prev_spkrateCN2_NoGo = INIT_RATE_CN;
                        
        // ---------------- RUN STATE -------------------
        for (int trial=0; trial<nTrials; trial++) {
            int progress = (int)trial*100/nTrials;
            if ( progress % 10 == 0 ) {
                std::cout << progress << "% ...";
                std::cout.flush();
            }
                        
            int cue_type = trials[trial];
            if (cue_type==GO) {
                stiRate1.setRates(RATE_MF);
                sim.setSpikeRate(gMF1, &stiRate1);
                                
                stiRate2.setRates(0.0f);
                sim.setSpikeRate(gMF2, &stiRate2);
                
                errRate1.setRates(curr_err_rate_1);
                sim.setSpikeRate(gErr1, &errRate1);
                                
                errRate2.setRates(BASE_ERROR_RATE);
                sim.setSpikeRate(gErr2, &errRate2);
                
                float bias_coupling = -CHANGE_RATE_IO1*(prev_spkrateCN1_Go-INIT_RATE_CN);
                sim.biasWeights(cnnIO1, bias_coupling);
                
                bias_coupling = -CHANGE_RATE_IO2*(prev_spkrateCN2_Go-INIT_RATE_CN);
                sim.biasWeights(cnnIO2, bias_coupling);
                                    
            }else {
                stiRate1.setRates(0.0f);
                sim.setSpikeRate(gMF1, &stiRate1);
                
                stiRate2.setRates(RATE_MF);
                sim.setSpikeRate(gMF2, &stiRate2);
                                
                errRate1.setRates(BASE_ERROR_RATE);
                sim.setSpikeRate(gErr1, &errRate1);
                
                errRate2.setRates(curr_err_rate_2);
                sim.setSpikeRate(gErr2, &errRate2);
                                                
                float bias_coupling = -CHANGE_RATE_IO2*(prev_spkrateCN2_NoGo-INIT_RATE_CN);
                sim.biasWeights(cnnIO2, bias_coupling);
                
                bias_coupling = -CHANGE_RATE_IO1*(prev_spkrateCN1_NoGo-INIT_RATE_CN);
                sim.biasWeights(cnnIO1, bias_coupling);
            }            
                                    
            smIO1->startRecording();
            smPC1->startRecording();
            smCN1->startRecording();
            
            smIO2->startRecording();
            smPC2->startRecording();
            smCN2->startRecording();
                        
            sim.runNetwork(0,trial_duration);
            
            smPC1->stopRecording();
            smIO1->stopRecording();
            smCN1->stopRecording();
            cmGC1_PC1->takeSnapshot();
            
            smPC2->stopRecording();
            smIO2->stopRecording();
            smCN2->stopRecording();
            cmGC2_PC2->takeSnapshot();
            
            std::vector<std::vector<float>> coupling1 = cmCoupling1->takeSnapshot();
            float curr_coupling1 = mean_omitnan(coupling1);
            std::vector<std::vector<float>> coupling2 = cmCoupling2->takeSnapshot();
            float curr_coupling2 = mean_omitnan(coupling2);
                        
            //reset coupling strength
            sim.biasWeights(cnnIO1, INIT_COUPLING-curr_coupling1);
            sim.biasWeights(cnnIO2, INIT_COUPLING-curr_coupling2);
                                          
            float spkrateCN1 = smCN1->getPopMeanFiringRate();
            float spkrateCN2 = smCN2->getPopMeanFiringRate();
                        
            float lick_rate = spike_rate_to_lick_rate(spkrateCN1+spkrateCN2);
            if (lick_rate>MAX_LICK_COUNT) {lick_rate=MAX_LICK_COUNT;}
            Lick.push_back(lick_rate);
                                                
            if (cue_type==GO) {
                prev_spkrateCN1_Go = spkrateCN1;
                prev_spkrateCN2_Go = spkrateCN2;
            } else {
                prev_spkrateCN1_NoGo = spkrateCN1;
                prev_spkrateCN2_NoGo = spkrateCN2;
            }
                                    
            // compute the error rate
            float err_lick = 0;
            if (cue_type==GO) {
                // error rate
                err_lick = MAX_LICK_COUNT-lick_rate;
                curr_err_rate_1 = error_lick_to_poisson_rate(err_lick);
            } else {
                // error rate
                err_lick = -lick_rate;
                curr_err_rate_2 = error_lick_to_poisson_rate(err_lick);
            }
            
            float bias = 0;
            // updating pf-PC weight
            float pop_spkrate_IO1 = smIO1->getPopMeanFiringRate();
            bias = -LEARNING_RATE_LTP*(pop_spkrate_IO1 - BASELINE_RATE_IO_LTP);
            sim.biasWeights(cnnGC1_PC1, bias);
            
            float pop_spkrate_IO2 = smIO2->getPopMeanFiringRate();
            bias = -LEARNING_RATE_LTD*(pop_spkrate_IO2 - BASELINE_RATE_IO_LTD);
            sim.biasWeights(cnnGC2_PC2, bias);
        }
        
        write_floats_to_file(fnLick, Lick);
        
        long t2 = clock();
        double time = (double)(t2 - t1) / CLOCKS_PER_SEC;
        
        std::cout << "done! Running time: " << time << " sec ===" << std::endl;
    }
    
    return 0;
}

