/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   BrownianDynamicsWithEP.hpp
 * Author: jzhan134
 *
 * Created on September 26, 2019, 10:17 AM
 */

#ifndef BROWNIANDYNAMICSWITHEP_HPP
#define BROWNIANDYNAMICSWITHEP_HPP
#include "BrownianDynamicEngine.hpp"
using namespace std;
class BrownianDynamicsWithEP : public BrownianDynamicEngine {
public:

    void initialization(std::string startingConfig);

    void coreBD(string fieldFileName, double fieldStrength_,
        double execTime = 1, FIELDROTATION fieldDirc_ = FIELDROTATION::DEFAULT);
    
    void migrateParticlesUnderEP(int);
    
    void ReadE(const std::string filename, int EPFlag);
    
    void barycentricInterpolation();

    int countNum(int gridX, int gridY){
        int num = 0;
        for (int i = 0; i < np; i++){
            if (fabs(p[i].x - (gridX-1)*60000) <= 30000 && fabs(p[i].y- (gridY-1)*60000) <= 30000){
                num++;
            }
        }
        return num;
    }

};
#endif /* BROWNIANDYNAMICSWITHEP_HPP */

