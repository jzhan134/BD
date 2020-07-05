/*
 * starting configuration:
 *      file_name - use a specific file
 *      library - use one of the configurations from library folder
 *      random - randomize particles into dilute phase
 */
#include "BasicDataStructures.h"
#include "BrownianDynamicEngine.hpp"
#include "BrownianDynamics3D.hpp"
#include "BrownianDynamicsWithEP.hpp"

using namespace std;

void controlScheme(BrownianDynamicEngine*, int);
void controlOnCurvedSurface(BrownianDynamicEngine*);
void isotropicControlScheme(BrownianDynamicEngine*);
void anisotropicControlScheme(BrownianDynamicEngine*);
void openLoopAnisotropicControlScheme(BrownianDynamicEngine*);
void controlWithElectrodeArray(BrownianDynamicEngine* model);
void equalizeClusterSize(BrownianDynamicEngine* model);
void getLibraryConfig(BrownianDynamicEngine* model);
void getRawTrajectory(BrownianDynamicEngine* model);

const int totalSimTime(500);
const int epochTime(20);
const int num_thread(1);
int totExecCycle(1), execCycle(0);
std::default_random_engine generator;
std::uniform_real_distribution<double> distribution(0.0,1.0);

int main() {
    thread *threads = new thread[num_thread];
    for (int idx = 0; idx < min(num_thread, totExecCycle); idx++){
//        BrownianDynamicEngine* model = new BrownianDynamicEngine();
//        threads[idx] = thread(anisotropicControlScheme, model);
//        
//        BrownianDynamicEngine* model = new BrownianDynamics3D();
//        threads[idx] = thread(controlOnCurvedSurface, model);
        
        BrownianDynamicEngine* model = new BrownianDynamicsWithEP();
        // threads[idx] = thread(controlWithElectrodeArray, model);
       threads[idx] = thread(equalizeClusterSize, model);
    }
    
    for (int idx = 0; idx < num_thread && threads[idx].joinable(); idx++){
        threads[idx].join();
    }
    delete[] threads;
    return 0;
}

void getLibraryConfig (BrownianDynamicEngine* model){
    while (execCycle++ < totExecCycle){
        model->getDefectedState("./library/fields/Assembly/600/10v10.txt", 2.1);
        model->initialization("library");
    }
}

void controlScheme(BrownianDynamicEngine* model, int idx){
    while (execCycle++ < totExecCycle){
        stringstream num;
        num<< execCycle%4*2+4;
        model->initialization("random");
        model->coreBD("./library/fields/Assembly/300/10v"+num.str()+".txt", 2, 20);
    }
//    while (execCycle++ < totExecCycle){
//        model->initialization("300_Crystal.txt");
//    model->getDefectedState("./library/fields/Assembly/10v10.txt", 1);
//        model->coreBD("./library/fields/Assembly/10v3.txt", 5, 50);
//        model->coreBD("./library/fields/Assembly/300/10v10.txt", 1.4, 50);
//    }

//    model->coreBD("./library/fields/Assembly/10v10.txt", 1, 20);
//    model->coreBD("./library/fields/Assembly/10v4.txt", 1, 100);
    
//    model->initialization("900_loose.txt");
//    model->coreBD("./library/fields/Assembly/900/10v10.txt", 2.4, 20);
//    model->coreBD("./library/fields/Assembly/900/10v4.txt", 2.4, 100);
}

void isotropicControlScheme(BrownianDynamicEngine* model){
    double fieldStrength(1);
    while (execCycle++ < totExecCycle){
        model->initialization("library");
        for (int t = 0; t < totalSimTime && model->getPsi6() < 0.97; t++){
            if (t%epochTime == 0) {
                fieldStrength = 1;
            }
            if (t%(2*epochTime) == 0){ // relax mode
//                 closed loop
//                if (model->getPsi6() < 0.6){
//                    fieldStrength = 0.1;
//                } else if (model->getPsi6() < 0.8) {
//                    fieldStrength = 0.2;
//                } else {
//                    fieldStrength = 0.3;
//                }
//                
                // open loop
                fieldStrength = 0.1;
            }
            model->coreBD("./library/fields/Assembly/10v10.txt", fieldStrength);
        }
    }
}

//2.7: 900 1.4: 300, 2.3: 600
void anisotropicControlScheme(BrownianDynamicEngine* model){
    FIELDROTATION currRot;
    double V = 2.7;
    while (execCycle++ < totExecCycle){
        int num = execCycle-1;
        
        // initial configuration
//        model->initialization("library300");
        model->getDefectedState("./library/fields/Assembly/900/10v10.txt", 1.2*2.7);
        
        // step 1: form defect free structure
        for (int period = 0; period < totalSimTime/epochTime; period++){
//            currRot = DEFAULT;
//            currRot = (distribution(generator) < 0.5) ? DEFAULT : ROT90;
            currRot = model->getPrefFieldByGB();
            currRot = (currRot == DEFAULT) ? ROT90 : DEFAULT;
            for (int t = 0; t < epochTime && model->getPsi6() < 0.97; t++){
                model->coreBD("./library/fields/Assembly/900/10v4.txt", V, 1, currRot);
            }
            for (int t = 0; t < epochTime && model->getPsi6() < 0.97; t++){
                model->coreBD("./library/fields/Assembly/900/10v10.txt", V, 1);
            }
        }
        
        if (model->getPsi6() < 0.97) { continue; } // exit from an failed cycle
        
        // step 2: restore circular morphology
        double currMorph = model->getMorphology();
        if (currMorph > 0){
            while (model->getMorphology() > 0.05){
                model->coreBD("./library/fields/Assembly/900/10v4.txt", V, 0.1, ROT90);
            }
        } else {
            while (model->getMorphology() < 0.05){
                model->coreBD("./library/fields/Assembly/900/10v4.txt", V, 0.1);
            }
        }
        
//        // stage 3: hold with isotropic field
        for(int t = 0; t < 50 && model->getPsi6() < 0.98; t++) {
            model->coreBD("./library/fields/Assembly/900/10v10.txt", V);
        }
    }
}

void getRawTrajectory(BrownianDynamicEngine* model){
    FIELDROTATION currRot;
    while (execCycle++ < totExecCycle){
        
        // create an random starting configuration
//        model->initialization("library300");
        model->getDefectedState("./library/fields/Assembly/300/10v10.txt", 1.2);
        
        // control period
        for (int period = 0; period < 20; period++){
            
            // state
            double state1 = model->getGB();
            double state2 = model->getPsi6();
            
            
            // action
            double randomNum = distribution(generator);
            string randomCircularity;
            if (randomNum<0.25){
                randomCircularity = "4";
            } else if (randomNum<0.5){
                randomCircularity = "5";
            } else if (randomNum<0.75){
                randomCircularity = "6";
            } else {
                randomCircularity = "7";
            }
            currRot = (distribution(generator) < 0.5) ? DEFAULT : ROT90;
            
            // simulate
            for (int t = 0; t < epochTime && model->getPsi6() < 0.97; t++){
                model->coreBD("./library/fields/Assembly/300/10v" + randomCircularity + ".txt", 1.2, 1, currRot);
            } 
            for (int t = 0; t < epochTime && model->getPsi6() < 0.97; t++){
                model->coreBD("./library/fields/Assembly/300/10v10.txt", 1.2, 1);
            }
            
            // reward
            double reward1 = model-> getPsi6();
            double reward2 = model->getMorphology();
            
            // output the experience
            cout << state1 << "\t" << state2 << "\t" << randomCircularity << "\t" << currRot << "\t" << reward1 << "\t" << reward2 << endl;
            
            // end if the assembly is done
            if (model->getPsi6() >= 0.97) break;
            
            // restore to  isotropic morphology
            double currMorph = model->getMorphology();
            if (currMorph > 0){
                while (model->getMorphology() > 0.05){
                    model->coreBD("./library/fields/Assembly/300/10v4.txt", 1.2, 0.1, ROT90);
                }
            } else {
                while (model->getMorphology() < 0.05){
                    model->coreBD("./library/fields/Assembly/300/10v4.txt", 1.2, 0.1);
                }
            }
        }
    }
}

void openLoopAnisotropicControlScheme(BrownianDynamicEngine* model){
    FIELDROTATION currRot;
    string currFile; 
    int t1;
    while (execCycle++ < totExecCycle){
        int num = execCycle-1;
        currRot = (distribution(generator)>0.5) ? DEFAULT : ROT90;
        model->initialization("library");
        // step 1: form defect free structure
        for (int period = 0; period < totalSimTime/epochTime; period++){
            currRot = (currRot == ROT90) ? DEFAULT : ROT90;
            for (int t = 0; t < epochTime && (model->getPsi6() < 0.95 || fabs(model->getMorphology()) > 0.07); t++){
                model->coreBD("./library/fields/Assembly/10v4.txt", 1.2, 0.1, currRot);
            }
            for (int t = 0; t < epochTime && (model->getPsi6() < 0.95 || fabs(model->getMorphology()) > 0.07); t++){
                model->coreBD("./library/fields/Assembly/10v10.txt", 1.2, 0.1);
            }
        }
        model->coreBD("./library/fields/Assembly/10v10.txt", 1.5, 10);
    }
}

void controlOnCurvedSurface(BrownianDynamicEngine* model){
//    model->initialization("900_R40.txt");
    model->initialization("random");
    model->coreBD("./library/fields/3D/900_10v10.txt", 10, 10);
    
    // 300 particles conditions
   model->coreBD("./library/fields/3D/10v10.txt", 3, 10);
   model->coreBD("./library/fields/3D/10v10.txt", 3, 15);
   model->coreBD("./library/fields/3D/10v6.txt", 3, 15, FIELDROTATION::ROT90);
   model->coreBD("./library/fields/3D/10v6.txt", 3, 15);
}

void controlWithElectrodeArray(BrownianDynamicEngine* model){
    model->initialization("random");
    model->coreBD("./library/fields/array/3x3_frame.txt", 0, 50);
    model->coreBD("./library/fields/array/3x3_frame.txt", 0.5, 2);
    model->coreBD("./library/fields/array/3x3_frame.txt", 1, 2);
    model->coreBD("./library/fields/array/3x3_frame.txt", 1.5, 2);
    model->coreBD("./library/fields/array/3x3_frame.txt", 2, 2);
    model->coreBD("./library/fields/array/3x3_frame.txt", 5, 12);
}


void equalizeClusterSize(BrownianDynamicEngine* model){
    model->initialization("3600_quench.txt");
    // model->coreBD("./library/fields/array/3x3_quad.txt", 0.5, 10);
    model->migrateParticlesUnderEP(400);
    model->coreBD("./library/fields/array/3x3_quad.txt", 5, 10);
}