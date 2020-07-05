/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   BrownianDynamics3D.hpp
 * Author: jzhan134
 *
 * Created on September 25, 2019, 4:53 PM
 */

#ifndef BROWNIANDYNAMICS3D_HPP
#define BROWNIANDYNAMICS3D_HPP
#pragma once
#include "BrownianDynamicEngine.hpp"

using namespace std;
extern const int total_time;
extern ofstream summary;
class BrownianDynamics3D : public BrownianDynamicEngine {
public:    
    void initialization(std::string startingConfig);
        
    void coreBD(
        string fieldFileName, 
        double fieldStrength_,
        double execTime = 1,
        FIELDROTATION fieldDirc_ = FIELDROTATION::DEFAULT
    );
    
    void ReadE(std::string filename, int);
    
    void CalOp();
    
    void RandomizeConfig();

    inline void ProjectOntoSurf(BD_Particle& p){
        if (sqrt(p.x * p.x + p.y * p.y) <= rim) {
            double r = sqrt(p.x * p.x + p.y * p.y + (p.z + zCenter) * (p.z + zCenter));
            p.x -= p.x * (r - R) / r;
            p.y -= p.y * (r - R) / r;
            p.z -= (p.z + zCenter) * (r - R) / r ;
        } else { 
            p.z = 0;
        }
    }
};

#endif /* BROWNIANDYNAMICS3D_HPP */

