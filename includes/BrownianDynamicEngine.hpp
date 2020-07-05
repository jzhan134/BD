#ifndef BROWNIANDYNAMICENGINE_HPP
#define BROWNIANDYNAMICENGINE_HPP
#pragma once

#include "BasicDataStructures.h"

using namespace std;
extern const int total_time;
extern ofstream summary;
class BrownianDynamicEngine{
public:
    BrownianDynamicEngine(){
        IMAGEANALYSISNEEDED = true;
        prevCheckSum = 0;
        ReadDiffusivity();
    }
    ~BrownianDynamicEngine(){
        if (trajOs.is_open()){ trajOs.close();}
        if (opOs.is_open()){ opOs.close();}
    }
    
    virtual void initialization(std::string startingConfig);
    
    virtual void coreBD( string fieldFileName, double fieldStrength_,
        double execTime = 1,FIELDROTATION fieldDirc_ = FIELDROTATION::DEFAULT);
    
    virtual void migrateParticlesUnderEP(int){};
    
    virtual void RandomizeConfig();
    
    void getDefectedState(string fieldFileName, double fieldStrength_);
    
    void CalDss(); 

    virtual void barycentricInterpolation();    

    // Image Analysis core
    virtual void CalOp();
    
    void ImageAna();
    
    vector<vector<int>> Connectivity(vector<int>&);
    
    virtual int countNum(int gridX, int gridY){return 0;}
    
    void ReadDiffusivity();
    
    // overwrite for 3D field
    virtual void ReadE(std::string filename, int);
    
    
    virtual void ProjectOntoSurf(BD_Particle& p){};
    
    void Readxyz(const string);
    
    void OutputTrajectory(ostream&);
    
    void OutputOrderParameter(ostream&);
    
    inline double getPsi6(){ return ensemble.psi6;}
    
    inline FIELDROTATION getPrefFieldByGB(){
        return (fabs(ensemble.GB_fit-120) > 45)? 
            FIELDROTATION::ROT90 : FIELDROTATION::DEFAULT;
    } 
    
    inline double getGB(){ return ensemble.GB_fit;} 
    
    inline FIELDROTATION getPrefFieldByMorph(){
        return (ensemble.Ix - ensemble.Iy < 0) ? 
            FIELDROTATION::DEFAULT : FIELDROTATION::ROT90;
    }
    
    inline double getMorphology(){
        return (ensemble.Ix  > ensemble.Iy) ? 
            sqrt(fabs(ensemble.Ix - ensemble.Iy))/(ensemble.rg/a) :
            -sqrt(fabs(ensemble.Ix - ensemble.Iy))/(ensemble.rg/a);
    }
    
    inline FIELDROTATION getFieldDirection() {return fieldDirection;}
    
    inline State<> GetCurrState() const {return ensemble;} 
    
    inline double Distance2D(BD_Particle& p1, BD_Particle& p2) const {
        return sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2));
    }
    
    inline double Distance3D(BD_Particle& p1, BD_Particle& p2) const {
        return sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2) + pow(p1.z - p2.z, 2));
    }
    
    inline bool isOutsideOfTable(int idx){
        return (idx < 0 || idx > E_tot);
    }

protected:
    // constants
    const double pi = 3.14159265;
    const int nstep = 10000;        // steps to simulate 1s
    const double dt = 0.1;          // (ms) simulation time of each step
    const double a = 800;        // particle radius (nm)
    const double rim = 20000;       // intersection of hemisphere on xy plain
    const double R = 20000;        // radius of hemisphere surface
    const double zCenter = sqrt(R * R - rim * rim);
    const double tempr = 20.0;      // temperature (C)
    const double fcm = -0.4667;     //Clausius–Mossotti Factor
    const double kb = 1.380658e-23; // Boltzmann constant
    const double kappa = a/10;      // Debye length
    const double pfpp = 2.2975*a;   // electrostatic repulsion
    const double fac1 = 5.9582e7/a;
    const double fac2 = 40.5622*sqrt((273 + tempr) / a) / sqrt(dt);
    const double gravFactor = 4/3*pi*pow(a*1e-9,3) * 1.65 * 9.8 * 1e-9 * 1e21;
    const double FppFactor = 1e18*kb*(tempr+273)*kappa*pfpp/a;
    const double Fhw = 0.417; // hard-sphere repulsion 
    const double rcut = 5.0*a; // range of neighbor for force calculation
    const double rmin = 2.64*a; // range of neighbor for order parameter calculation
    const double ppPreFactor = 1e9 * (273 + tempr) * kb * 21.86 * 6.25 * 8 * pow(a,3);
    const double pfPreFactor = 1e18 * 2 * pi * 80 * 8.85e-12 * pow(a*1e-9, 3) * fcm;
    const double epPreFactor = 1.5 * 80 * 8.85e-12 * 1e-4 / (a*1e-9);
    
    const double d_idx = 250; // (nm) resolution of field loop-up table
    
    // diffusivity maps
    static const int rgdssbin = 25;
    static const int distdssbin = 50;
    int dsscount[rgdssbin][distdssbin];
    double dssarray[rgdssbin][distdssbin]; 
    
    // path of particle starting configurations
    std::string configPath = "./library/configurations/";

    bool IMAGEANALYSISNEEDED;
    BD_Particle p[np]; // coordinate and force of each particle
    OP_para op[np]; // local order parameters
    vector<double> GB_queue; // time series of grain boundary in last 5s
    
    // data structure for lookup tables
    double ETable[E_tot][E_tot][2];
    double dETable[E_tot][E_tot][4];
    double dE2Table[E_tot][E_tot][3];
    double dE2TableII[E_tot][E_tot][2];
    
    double E[np][2]; 
    double dE[np][4];
    double dE2[np][3];
    double EP_E[np][2];
    
    State<> ensemble; // current state of ensemble
    double fieldStrength; // applied field strength with respect to the table
    FIELDROTATION fieldDirection; 
    int elapsedTime; // elapsed time in a complete simulation cycle
    int currCycle;
    static int repCycle;
    int prevCheckSum;
    int prevEPCheckSum;
    ofstream trajOs, opOs; // output file names

    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution<> rand_normal{0.0,1.0};
//    std::normal_distribution<> rand_normal2{-1.0,1.0};
    std::uniform_real_distribution<> rand_uniform{0.0,1.0};
};
#endif /* BROWNIANDYNAMICS_HPP */