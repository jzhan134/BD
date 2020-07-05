#include "BrownianDynamicsWithEP.hpp"

void BrownianDynamicsWithEP::initialization(std::string startingConfig)
{
    // reset timer for the current cycle
    elapsedTime = 0;
    
    // update finished number of cycles
    stringstream repCycleOS;
    repCycleOS << repCycle++;
    
    // read initial configuration
    Readxyz(startingConfig);
//    ReadE("./library/fields/array/3x3_frame.txt", 0);
//    ReadE("./library/fields/array/3x3_quad.txt", 0);
    // create new trajectory output file
    if (trajOs.is_open()) trajOs.close();
    trajOs.open("./traj/xyz" + repCycleOS.str() + ".dat");  
    OutputTrajectory(trajOs);
}

// time step: 1e-4s
// total time: 0.1s
// output freq: 0.1s
void BrownianDynamicsWithEP::coreBD(
        string fieldFileName,
        double fieldStrength_, 
        double execTime,
        FIELDROTATION fieldDirc_)  
{
    double rijsep,xij, yij, distIJ, distOffSet;
    double dx, dy;
    execTime *= 10;
    double pfFactor = pfPreFactor * pow(fieldStrength_,2);
    int cluster[np] = {0};
    ReadE(fieldFileName, 0);
    
    for (int t = 0; t < execTime; t++){
        CalDss();
        barycentricInterpolation();
        
        // combine particles into clusters
        for (int i = 0; i < np; i++){
            cluster[i] = round(p[i].x / 60000)*3 + round(p[i].y / 60000);
        }
        
        // 0.1s of simulation is divided in 2000 steps of 0.05ms each
        // 0.1s of simulation is divided in 1000 steps of 0.01ms each
        for (int step = 0; step < nstep/10; step++){ 
            for (int i = 0; i < np; i++){
                p[i].Fx = 0.0;
                p[i].Fy = 0.0;
            }
            #pragma omp parallel for
            for (int i = 0; i < np; i++) {
                for (int j = i+1; j < np; j++){
                    // if (cluster[i] == cluster[j]){
                        rijsep = Distance2D(p[i],p[j]);
                        if (rijsep < rcut){
                            xij = p[j].x - p[i].x;
                            yij = p[j].y - p[i].y;
                            if (rijsep < rcut){
                                xij = p[j].x - p[i].x;
                                yij = p[j].y - p[i].y;
                                if (rijsep < 2*a + 40){
                                    dx = xij * (2*a + 40 - rijsep)  /sqrt(xij * xij + yij * yij);
                                    dy = yij * (2*a + 40 - rijsep)  /sqrt(xij * xij + yij * yij);
                                    p[i].x -= dx;
                                    p[i].y -= dy;
                                    p[j].x += dx;
                                    p[j].y += dy;
                                } else {
                                    double Fpp = (rijsep < 2*a +20) ? 
                                        Fhw * (1 - (rijsep - 2*a) / 20) :
                                        1e18*kb*(tempr+273)*kappa*pfpp*
                                            exp(-kappa*(rijsep-2.0*a)/a)/a;
                                    p[i].Fx -= Fpp*xij/rijsep;
                                    p[i].Fy -= Fpp*yij/rijsep;
                                    p[j].Fx += Fpp*xij/rijsep;
                                    p[j].Fy += Fpp*yij/rijsep;
                                }
                            }
                        }
                    // }
                }
                // DEP particle-field interaction
                p[i].Fx += pfFactor * dE2[i][0];
                p[i].Fy += pfFactor * dE2[i][1];

                double randx(0), randy(0);
                randx =  rand_normal(gen) * sqrt(1.0 / p[i].D);
                randy =  rand_normal(gen) * sqrt(1.0 / p[i].D);

                p[i].x += p[i].D * (p[i].Fx*fac1 + randx*fac2) * dt;
                p[i].y += p[i].D * (p[i].Fy*fac1 + randy*fac2) * dt;

                // periodic boundary
                if (p[i].x > periodicWindow) p[i].x = (p[i].x - 2*periodicWindow);
                if (p[i].x < -periodicWindow) p[i].x = (p[i].x + 2*periodicWindow);
                if (p[i].y > periodicWindow) p[i].y = (p[i].y - 2*periodicWindow);
                if (p[i].y < -periodicWindow) p[i].y = (p[i].y + 2*periodicWindow);
            }
        }
        
        // track overall time stamp
        cout << (double) (elapsedTime++)/10 << endl; // increment every 0.1s
        OutputTrajectory(trajOs);
    }
}

void BrownianDynamicsWithEP::migrateParticlesUnderEP(int targetNum) 
{
    double rijsep,xij, yij, distIJ, distOffSet;
    double EPStrength;
    double pfFactor = pfPreFactor * 4;
    int iniNum = countNum(1, 1);
    int timeCount = 0;
    
    ReadE("./library/fields/array/3x3_quad.txt", 0); // read DEP field
    ReadE("./library/fields/array/3x3_EP.txt", 1); // read EP field
    
    while (countNum(1, 1) > targetNum){
        // check endpoint every 0.01s
        CalDss();
        barycentricInterpolation();
        EPStrength = (countNum(1, 1) == iniNum) ? 0.4 : 0.05;
        
        // 0.01s of simulation is divided in 500 steps of 0.02ms each,
        for (int step = 0; step < 500; step++){ 
            // check for exit condition every time step
            for (int i = 0; i < np; i++){
                p[i].Fx = 0.0;
                p[i].Fy = 0.0;
            }
            #pragma omp parallel for
            for (int i = 0; i < np; i++) {
                for (int j = i+1; j < np; j++){
                    rijsep = Distance2D(p[i],p[j]);
                    if (rijsep < rcut){
                        xij = p[j].x - p[i].x;
                        yij = p[j].y - p[i].y;
                        if (rijsep < rcut){
                            xij = p[j].x - p[i].x;
                            yij = p[j].y - p[i].y;
                            if (rijsep < 2*a + 40){
                                p[i].x -= xij * (2*a + 40 - rijsep)  /sqrt(xij * xij + yij * yij);
                                p[i].y -= yij * (2*a + 40 - rijsep)  /sqrt(xij * xij + yij * yij);
                                p[j].x += xij * (2*a + 40 - rijsep)  /sqrt(xij * xij + yij * yij);
                                p[j].y += yij * (2*a + 40 - rijsep)  /sqrt(xij * xij + yij * yij);
                            } else {
                                double Fpp = (rijsep < 2*a +20) ? 
                                    Fhw * (1 - (rijsep - 2*a) / 20) :
                                    1e18*kb*(tempr+273)*kappa*pfpp*
                                        exp(-kappa*(rijsep-2.0*a)/a)/a;
                                p[i].Fx -= Fpp*xij/rijsep;
                                p[i].Fy -= Fpp*yij/rijsep;
                                p[j].Fx += Fpp*xij/rijsep;
                                p[j].Fy += Fpp*yij/rijsep;
                            }
                        }
                    }
                }
                
                // DEP particle-field interaction
                p[i].Fx += pfFactor * dE2[i][0];
                p[i].Fy += pfFactor * dE2[i][1];
                
                // EP particle-field interaction
                p[i].Fx += epPreFactor * EPStrength * EP_E[i][0];
                p[i].Fy += epPreFactor * EPStrength * EP_E[i][1];

                double randx(0), randy(0);
                randx =  rand_normal(gen) * sqrt(1.0 / p[i].D);
                randy =  rand_normal(gen) * sqrt(1.0 / p[i].D);

                p[i].x += p[i].D * (p[i].Fx*fac1 + randx*fac2) * dt/5;
                p[i].y += p[i].D * (p[i].Fy*fac1 + randy*fac2) * dt/5;
                
                // periodic boundary
                if (p[i].x > periodicWindow) p[i].x = (p[i].x - 2*periodicWindow);
                if (p[i].x < -periodicWindow) p[i].x = (p[i].x + 2*periodicWindow);
                if (p[i].y > periodicWindow) p[i].y = (p[i].y - 2*periodicWindow);
                if (p[i].y < -periodicWindow) p[i].y = (p[i].y + 2*periodicWindow);
            }
        }
        timeCount++;
        
        // track overall time stamp
        if (timeCount%10 == 0){
            cout << (double) (elapsedTime++)/10 << endl; // increment every 0.1s
            OutputTrajectory(trajOs);
        }
    }
}

void BrownianDynamicsWithEP::barycentricInterpolation()
// interpolate electric field for each particles based on lookup tables.
// DEP force (dE2) interpolated from dE2Table
// EP force (E) interpolated from ETable
// both information are calculated even if only only would be used in the
// dynamic engine.

{
    int x[3],y[3];
    double dx,dy;
    int idx_x, idx_y;
    double det, lambda1, lambda2, lambda3;

    for (int i = 0; i < np; i++){  
        dx = p[i].x/d_idx + (E_tot-1)/2;
        dy = p[i].y/d_idx + (E_tot-1)/2;
        idx_x = floor (dx);
        idx_y = floor (dy);
        
        if (dx-idx_x >= 0.5) {
            x[0] = (double)idx_x + 1;
            x[1] = (double)idx_x;
            x[2] = (double)idx_x + 1;
        } else {
            x[0] = (double)idx_x;
            x[1] = (double)idx_x + 1;
            x[2] = (double)idx_x;
        }
        
        if (dy-idx_y >= (0.5)) {
            y[0] = (double)idx_y + 1;
            y[1] = (double)idx_y + 1;
            y[2] = (double)idx_y;
        } else{
            y[0] = (double)idx_y;
            y[1] = (double)idx_y;
            y[2] = (double)idx_y + 1;
        }
        
        // grid point cannot exceed table limit
        
        for (int i = 0; i < 3; i++){
            if (x[i] > E_tot) x[i] = E_tot;
            if (x[i] < 0) x[i] = 0;
            if (y[i] > E_tot) y[i] = E_tot;
            if (y[i] < 0) y[i] = 0;
        }
        
        det = (y[1] - y[2])*(x[0] - x[2]) + (x[2] - x[1])*(y[0] - y[2]);
        lambda1 = ((y[1] - y[2]) * (dx - x[2]) + (x[2] - x[1]) * (dy - y[2]))/det;
        lambda2 = ((y[2] - y[0]) * (dx - x[2]) + (x[0] - x[2]) * (dy - y[2]))/det;
        lambda3 = 1.0 - lambda1 - lambda2;
        
        // find weighted local field magnitudes
        E[i][0] =   lambda1*ETable[x[0]][y[0]][0]   + lambda2*ETable[x[1]][y[1]][0]   + lambda3*ETable[x[2]][y[2]][0];
        E[i][1] =   lambda1*ETable[x[0]][y[0]][1]   + lambda2*ETable[x[1]][y[1]][1]   + lambda3*ETable[x[2]][y[2]][1];
        dE2[i][0] = lambda1*dE2Table[x[0]][y[0]][0] + lambda2*dE2Table[x[1]][y[1]][0] + lambda3*dE2Table[x[2]][y[2]][0];
        dE2[i][1] = lambda1*dE2Table[x[0]][y[0]][1] + lambda2*dE2Table[x[1]][y[1]][1] + lambda3*dE2Table[x[2]][y[2]][1];
        dE2[i][0] /= 1e9;
        dE2[i][1] /= 1e9;
        
    }
}

void BrownianDynamicsWithEP::ReadE(const std::string filename, int EPFlag) {
    // if the same lookup table is used in the previous run, skip the reading process
    int chksum = 0;
    for (int i = 0; i < filename.length(); i++) chksum += (int) filename[i];
    if (prevCheckSum == chksum || prevEPCheckSum == chksum) {
        return;
    } else {
        if (EPFlag == 0) prevCheckSum = chksum;
        else if (EPFlag == 1) prevEPCheckSum = chksum;
    }
    ifstream is;
    is.open(filename.c_str());
    assert(is.is_open());
    string line;
    double dum;
    for (int i = 0; i < E_tot; i++) {
        for (int j = 0; j < E_tot; j++) {
            getline(is, line);
            stringstream linestream(line);
            linestream >> dum;
            linestream >> dum;
            if (EPFlag == 0){
                linestream >> dum;
                linestream >> dum;
                linestream >> dE2Table[j][i][0];
                linestream >> dE2Table[j][i][1];
            } else if (EPFlag == 1){
                linestream >> ETable[j][i][0];
                linestream >> ETable[j][i][1];
            }
        }
    }
    is.close();
}