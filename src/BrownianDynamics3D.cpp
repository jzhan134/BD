#include "BrownianDynamics3D.hpp"

/***************************** DYNAMICS SIMULATION ****************************/


void BrownianDynamics3D::initialization(std::string startingConfig)
{
    // reset timer for the current cycle
    elapsedTime = 0;
    
    // update finished number of cycles
    stringstream repCycleOS;
    repCycleOS << repCycle++;
    
    // read initial configuration
    Readxyz(startingConfig);
    
    // calculate the correct height such that the particles are on the surface
    for (int i = 0; i < np; i++){
        if (sqrt(p[i].x*p[i].x + p[i].y*p[i].y) < rim){
            p[i].z = sqrt(R*R - (p[i].x*p[i].x + p[i].y*p[i].y)) - zCenter;
        } else {
            p[i].z = 0;
        }
    }
    
    // create new trajectory output file
    if (trajOs.is_open()) trajOs.close();
    trajOs.open("./traj/xyz" + repCycleOS.str() + ".dat");  
    CalOp();
    OutputTrajectory(trajOs);
}

void BrownianDynamics3D::coreBD(
        string fieldFileName, 
        double fieldStrength_,
        double execTime,
        FIELDROTATION fieldDirc_) 
{
    // the dynamics of particle on a curved surface with perfect hemisphere
    // shape. The center is at (0,0,0), with radius R.
    // The algorithm is to find the normal displacement of particles in 3
    // dimensional space as normal Brownian Dynamics, and then project the 
    // displacement onto the hemisphere.
    // The force in x, y, and z axes are considered. Gravity is considered
    // assume the energy landscape is independent of z
    double rijsep,xij, yij, zij;
    double Fdepx, Fdepy;
    
    execTime *= 10;
    fieldStrength = fieldStrength_;
    fieldDirection = fieldDirc_;
    double ppFactor = ppPreFactor * pow(fieldStrength,2);
    double pfFactor = pfPreFactor * pow(fieldStrength,2);
    
    ReadE(fieldFileName, 0);
    
    for (int t = 0; t < execTime; t++){
        
        // update local particle parameters every 0.1s
        CalDss();
        barycentricInterpolation();

        // calculate the dynamics in 0.1s with step size of 0.1ms
        for (int step = 0; step < nstep; step++){        
            for (int i = 0; i < np; i++){
                p[i].Fx = 0.0;
                p[i].Fy = 0.0;
                p[i].Fz = 0.0;
            }
            #pragma omp parallel for
            for (int i = 0; i < np; i++) {

                // particle-particle pair-wise forces
                for (int j = i+1; j < np; j++){
                    rijsep = Distance3D(p[i],p[j]);
                    if (rijsep < rcut){
                        xij = p[j].x - p[i].x;
                        yij = p[j].y - p[i].y;
                        zij = p[j].z - p[i].z;
                        if (rijsep < 2 * a + 40){
                            p[i].x -= xij * (2*a + 40 - rijsep) /sqrt(xij * xij + yij * yij + zij * zij);
                            p[i].y -= yij * (2*a + 40 - rijsep) /sqrt(xij * xij + yij * yij + zij * zij);
                            p[i].z -= zij * (2*a + 40 - rijsep) /sqrt(xij * xij + yij * yij + zij * zij);
                            p[j].x += xij * (2*a + 40 - rijsep) /sqrt(xij * xij + yij * yij + zij * zij);
                            p[j].y += yij * (2*a + 40 - rijsep) /sqrt(xij * xij + yij * yij + zij * zij);
                            p[j].z += zij * (2*a + 40 - rijsep) /sqrt(xij * xij + yij * yij + zij * zij);
                            ProjectOntoSurf(p[i]);
                            ProjectOntoSurf(p[j]);
                        } else {
                            double Fpp = (rijsep < 2*a +20) ? Fhw * (1 - (rijsep - 2*a) / 20) :
                                1e18*kb*(tempr+273)*kappa*pfpp*
                                    exp(-kappa*(rijsep-2.0*a)/a)/a;
                            p[i].Fx -= Fpp*xij/rijsep;
                            p[i].Fy -= Fpp*yij/rijsep;
                            p[i].Fz -= Fpp*zij/rijsep;
                            p[j].Fx += Fpp*xij/rijsep;
                            p[j].Fy += Fpp*yij/rijsep;
                            p[j].Fz += Fpp*zij/rijsep;
                        }
                        

                        // field-induced dipole-dipole interaction
                        if (rijsep > 2*a){

                            double Er = xij*E[i][0] + yij*E[i][1];
                            double Emod2 = (E[i][0]*E[i][0] + E[i][1]*E[i][1]);
                            double cos2 = (Er/(rijsep*sqrt(Emod2)))*(Er/(rijsep*sqrt(Emod2)));
                            double rij5 = rijsep*rijsep*rijsep*rijsep*rijsep;

                            double F1 = -(E[i][0]*dE[i][0] + E[i][1]*dE[i][2])*Er*Er/Emod2
                                    + (xij*dE[i][0] + yij*dE[i][2] - E[i][0])*Er
                                    + xij*Er*Er/rijsep/rijsep;
                            F1 *= (3.0/rij5);
                            double F2 = 0.5*(3*cos2-1)*(3*xij/rij5)*Emod2;
                            double F3 = 0.5*(3*cos2-1)/(rijsep*rijsep*rijsep)*dE2[i][0];
                            Fdepx = ppFactor*(F1 + F2 + F3);

                            F1 = -(E[i][0]*dE[i][1] + E[i][1]*dE[i][3])*Er*Er/Emod2
                                    + (xij*dE[i][1] + yij*dE[i][3] - E[i][1])*Er
                                    + yij*Er*Er/rijsep/rijsep;
                            F1 *= (3.0/rij5);
                            F2 = 0.5*(3*cos2-1)*(3*yij/rij5)*Emod2;
                            F3 = 0.5*(3*cos2-1)/(rijsep*rijsep*rijsep)*dE2[i][1];
                            Fdepy = ppFactor*(F1 + F2 + F3);
                        } else {
                            Fdepx = 0.0;
                            Fdepy = 0.0;
                        }
//                        p[i].Fx += Fdepx;
//                        p[i].Fy += Fdepy;
//                        p[j].Fx -= Fdepx;
//                        p[j].Fy -= Fdepy;
                    }
                }

                // particle-field interaction
                p[i].Fx += pfFactor * dE2[i][0];
                p[i].Fy += pfFactor * dE2[i][1];
                p[i].Fz += pfFactor * dE2[i][2];

                // electrostatic repulsion between particle and substrate
                if (p[i].z < a) {
                    p[i].Fz += Fhw;
                } else {
                    p[i].Fz += 2*1e18*kb*(tempr+273)*kappa*pfpp*exp(-kappa*(p[i].z-a)/a)/a;
                }

                // gravity
//                p[i].Fz -= gravFactor;

                // random force
                double randx(0), randy(0), randz(0);
                randx =  rand_normal(gen) * sqrt(1.0 / p[i].D);
                randy =  rand_normal(gen) * sqrt(1.0 / p[i].D);
                randz =  rand_normal(gen) * sqrt(1.0 / p[i].D);

                // displacement update
                p[i].x += p[i].D * (p[i].Fx*fac1 + randx*fac2) * dt;
                p[i].y += p[i].D * (p[i].Fy*fac1 + randy*fac2) * dt;
                p[i].z += p[i].D * (p[i].Fz*fac1 + randz*fac2) * dt;
                ProjectOntoSurf(p[i]);
            }
        }
        elapsedTime++;
        if (elapsedTime%10 == 0){ 
            cout << (double) elapsedTime/10 << endl;
            CalOp();
            OutputTrajectory(trajOs);
        }
    }
}

void BrownianDynamics3D::ReadE(const std::string filename, int EPFlag) {
    int chksum = 0;
    for (int i = 0; i < filename.length(); i++){
        chksum += (int) filename[i];
    }
    if (prevCheckSum == chksum){
        return;
    } else {
        prevCheckSum = chksum;
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
            linestream >> ETable[j][i][0];
            linestream >> ETable[j][i][1];
            linestream >> dE2Table[j][i][0];
            linestream >> dE2Table[j][i][1];
            linestream >> dE2Table[j][i][2];
            linestream >> dETable[j][i][0];
            linestream >> dETable[j][i][1];
            linestream >> dETable[j][i][2];
            linestream >> dETable[j][i][3];
            }
        }
        is.close();
    }
}

void BrownianDynamics3D::CalOp() {
    for (int i = 0; i < np; i++){
        p[i].coordinate_number = 0;
    }
    for (int i = 0; i < np; i++) {
        for (int j = i + 1; j < np; j++){
            if (Distance3D(p[i],p[j]) < 2*a*1.366){
                (p[i].coordinate_number)++;
                (p[j].coordinate_number)++;
            }
        }
    }
}

void BrownianDynamics3D::RandomizeConfig(){
    double theta, psi;
    for (int i = 0; i < np; i++){
        theta = (double) (rand()%3600)/10/180*pi;
        psi = (double) (rand()%600)/10/180*pi;
//        cout << theta << " " << psi << endl;
        p[i].x = R*sin(psi)*cos(theta);
        p[i].y = R*sin(psi)*sin(theta);
        p[i].z = R*cos(psi) - zCenter;
    }
}