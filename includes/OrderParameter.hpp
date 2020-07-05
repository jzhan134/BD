/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   OrderParameter.hpp
 * Author: jzhan134
 *
 * Created on September 26, 2019, 10:03 AM
 */

#ifndef ORDERPARAMETER_HPP
#define ORDERPARAMETER_HPP
#include "BasicDataStructures.h"

class OrderParameter {
public:
    OrderParameter(BD_Particle* p_){
        for (int i = 0; i < np; i++){
            p[i] = p_++;
        }
    };
protected:
    BD_Particle* p[np];
};


#endif /* ORDERPARAMETER_HPP */

