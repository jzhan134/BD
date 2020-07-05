#pragma once
#include <cstdlib>
#include <chrono>
#include <iostream>
#include <thread>
#include <cstdlib>
#include <ctime>
#include <string>
#include <cmath>
#include <memory>
#include <random>
#include <vector>
#include <random>
#include <chrono>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <math.h>
#include <stdlib.h>
#include <queue>
#include "omp.h"
#include <cassert>
#include <map>

using namespace std;

enum CONTROLCATEGORY {
    CONST,
    ISO,
    ANISO,
    MANUAL,
    CURVEDSURF,
    ARRAY
};

enum PARTICLETYPE {
    INNER_RIM,
    PERIPHERAL,
    GRAIN,
    CRYSTAL
};

enum FIELDROTATION {
    DEFAULT,
    ROT90
};


enum STARTINGCONFIG {
    SPECIFIC,
    LIBRARY,
    RANDOM
};

struct BD_Particle{
    double x;
    double y;
    double z = 0;
    double D;
    double Fx;
    double Fy;
    double Fz;
    double Psi;
    double c6;
    double Theta;
    vector<int> op_neighbor;
    PARTICLETYPE type;
    int dirc;
    int coordinate_number;
};

struct OP_para{
    double Psi;
    double Theta;
    vector<int> op_neighbor;
    int type;
};

template<class T = double>
struct State{
    T psi6;
    T rg;
    T c6=0;
    double GB_fit; // [0 90] w.r.t. WE anisotropic field
    T Ix;
    T Iy;
};

const int np = 100;      // particle number
const int E_tot = 801; // size of field look-up table
const double periodicWindow = 90000;