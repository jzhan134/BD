#include "BrownianDynamicEngine.hpp"
int BrownianDynamicEngine::repCycle(0);
/***************************** DYNAMICS SIMULATION ****************************/
void BrownianDynamicEngine::initialization(
        std::string startingConfig)
{
    // reset timer for the current cycle
    elapsedTime = 0;
    
    // update finished number of cycles
    stringstream repCycleOS;
    repCycleOS<< repCycle++;
    
    // read initial configuration
    Readxyz(startingConfig);
    
    // create new trajectory output file
    if (trajOs.is_open()) trajOs.close();
    trajOs.open("./traj/xyz" + repCycleOS.str() + ".dat");  
    
    
    // create new order parameter output file
    GB_queue.clear();
    if (opOs.is_open()) opOs.close();
    opOs.open("./traj/op" + repCycleOS.str() + ".dat");
    opOs << fixed;
    opOs << setprecision(2);
    opOs << setw(4);
    CalOp();
    ImageAna();
    OutputOrderParameter(opOs);
    OutputTrajectory(trajOs);
}

// runs BD core engine for a period of time during which the same type of field
// is used. The simulation outcome is saved every 0.1s.
void BrownianDynamicEngine::coreBD(
        string fieldFileName, 
        double fieldStrength_, 
        double execTime,
        FIELDROTATION fieldDirc_) 
{
    double rijsep,xij, yij;
    double Fdepx, Fdepy;
    
    execTime *= 10;
    fieldStrength = fieldStrength_;
    fieldDirection = fieldDirc_;
    double ppFactor = ppPreFactor * pow(fieldStrength,2);
    double pfFactor = pfPreFactor * pow(fieldStrength,2);
    
    ReadE(fieldFileName,0);
    for (int t = 0; t < execTime; t++){
        
        // update local particle parameters every 0.1s
        CalDss();
        barycentricInterpolation();
        
        // calculate the dynamics in 0.1s with step size of 0.1ms
        for (int step = 0; step < nstep; step++){ 
            
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
                        if (rijsep < 2 * a + 20 ){
                            p[i].x -= xij * (2*a + 20 - rijsep) /sqrt(xij * xij + yij * yij);
                            p[i].y -= yij * (2*a + 20 - rijsep) /sqrt(xij * xij + yij * yij);
                            p[j].x += xij * (2*a + 20 - rijsep) /sqrt(xij * xij + yij * yij);
                            p[j].y += yij * (2*a + 20 - rijsep) /sqrt(xij * xij + yij * yij);
                        } else {
                            double Fpp = (rijsep < 2*a +20) ? Fhw * (1 - (rijsep - 2*a) / 20) :
                                1e18*kb*(tempr+273)*kappa*pfpp*
                                    exp(-kappa*(rijsep-2.0*a)/a)/a;
                            p[i].Fx -= Fpp*xij/rijsep;
                            p[i].Fy -= Fpp*yij/rijsep;
                            p[j].Fx += Fpp*xij/rijsep;
                            p[j].Fy += Fpp*yij/rijsep;
                        }

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
                        p[i].Fx += Fdepx;
                        p[i].Fy += Fdepy;
                        p[j].Fx -= Fdepx;
                        p[j].Fy -= Fdepy;
                    }
                }

                p[i].Fx += pfFactor * dE2[i][0];
                p[i].Fy += pfFactor * dE2[i][1];

                double randx(0), randy(0);
                randx =  rand_normal(gen) * sqrt(1.0 / p[i].D);
                randy =  rand_normal(gen) * sqrt(1.0 / p[i].D);

                p[i].x += p[i].D * (p[i].Fx*fac1 + randx*fac2) * dt;
                p[i].y += p[i].D * (p[i].Fy*fac1 + randy*fac2) * dt;
            }
        }
        
        // store trajectory and order parameter every 1s
        elapsedTime++;
        if (elapsedTime%10 == 0){
//            cout << (double) elapsedTime/10 << endl;
            CalOp();
            ImageAna();
            OutputOrderParameter(opOs);
            OutputTrajectory(trajOs);
        }
    }
}

// use Barycentric coordinate system to approximate the local magnitude of field
// at specific position. 
// The first derivative of E is listed as: dExdx[], dExdy[], dEydx[], dEydy[]
// The field can be used after rotating by 90 degrees, which requires a rotation
// of particle coordinate
void BrownianDynamicEngine::barycentricInterpolation()
{
    int x[3],y[3];
    double dx,dy;
    int idx_x, idx_y;
    double det, lambda1, lambda2, lambda3;
    // lookup table can be scaled with respect to particle number with sqrt(N)

    for (int i = 0; i < np; i++){  
        dx = (fieldDirection == FIELDROTATION::ROT90) ? -p[i].y/d_idx : p[i].x/d_idx;
        dy = (fieldDirection == FIELDROTATION::ROT90) ? p[i].x/d_idx : p[i].y/d_idx;
        dx += (E_tot-1)/2;
        dy += (E_tot-1)/2;
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
        
        // weighted-average fields based on Barycentric coordinates
        if (isOutsideOfTable(x[0]) || isOutsideOfTable(x[1]) || isOutsideOfTable(x[2]) ||
                isOutsideOfTable(y[0]) || isOutsideOfTable(y[1]) || isOutsideOfTable(y[2])){
            lambda1 = 1/3;
            lambda2 = 1/3;
            lambda3 = 1/3;
        } else {
            det = (y[1] - y[2])*(x[0] - x[2]) + (x[2] - x[1])*(y[0] - y[2]);
            lambda1 = ((y[1] - y[2]) * (dx - x[2]) + (x[2] - x[1]) * (dy - y[2]))/det;
            lambda2 = ((y[2] - y[0]) * (dx - x[2]) + (x[0] - x[2]) * (dy - y[2]))/det;
            lambda3 = 1.0 - lambda1 - lambda2;
        }
        
        // find weighted local field magnitudes
        E[i][0] =   lambda1*ETable[x[0]][y[0]][0]   + lambda2*ETable[x[1]][y[1]][0]   + lambda3*ETable[x[2]][y[2]][0];
        E[i][1] =   lambda1*ETable[x[0]][y[0]][1]   + lambda2*ETable[x[1]][y[1]][1]   + lambda3*ETable[x[2]][y[2]][1];
        dE2[i][0] = lambda1*dE2Table[x[0]][y[0]][0] + lambda2*dE2Table[x[1]][y[1]][0] + lambda3*dE2Table[x[2]][y[2]][0];
        dE2[i][1] = lambda1*dE2Table[x[0]][y[0]][1] + lambda2*dE2Table[x[1]][y[1]][1] + lambda3*dE2Table[x[2]][y[2]][1];
        dE2[i][2] = lambda1*dE2Table[x[0]][y[0]][2] + lambda2*dE2Table[x[1]][y[1]][2] + lambda3*dE2Table[x[2]][y[2]][2];
        dE[i][0] =  lambda1*dETable[x[0]][y[0]][0]  + lambda2*dETable[x[1]][y[1]][0]  + lambda3*dETable[x[2]][y[2]][0];
        dE[i][1] =  lambda1*dETable[x[0]][y[0]][1]  + lambda2*dETable[x[1]][y[1]][1]  + lambda3*dETable[x[2]][y[2]][1];
        dE[i][2] =  lambda1*dETable[x[0]][y[0]][2]  + lambda2*dETable[x[1]][y[1]][2]  + lambda3*dETable[x[2]][y[2]][2];
        dE[i][3] =  lambda1*dETable[x[0]][y[0]][3]  + lambda2*dETable[x[1]][y[1]][3]  + lambda3*dETable[x[2]][y[2]][3];
        dE2[i][0] /= 1e9;
        dE2[i][1] /= 1e9;
        dE2[i][2] /= 1e9;
        dE[i][0] /= 1e9;
        dE[i][1] /= 1e9;
        dE[i][2] /= 1e9;
        dE[i][3] /= 1e9;
        
        // additional operation is needed if a rotated field is used
        if (fieldDirection == FIELDROTATION::ROT90){
            double temp;
            temp = E[i][0];
            E[i][0] = E[i][1];
            E[i][1] = -temp;
            
            temp = dE2[i][0];
            dE2[i][0] = dE2[i][1];
            dE2[i][1] = -temp;
            
            temp = dE[i][0];
            dE[i][0] = dE[i][1];
            dE[i][1] = -temp;
            
            temp = dE[i][2];
            dE[i][2] = dE[i][3];
            dE[i][3] = -temp;
        }
    }
    
}

/********************* ORDER PARAMETER AND IMAGE ANALYSIS *********************/
void BrownianDynamicEngine::CalOp() {
    
    double psir[np]{0.0}, psii[np]{0.0}, con[np]{0.0}; 
    double xmean(0.0), ymean(0.0), rgmean(0.0);
    double accumpsi6r(0.0), accumpsi6i(0.0);
    double degree;
    
    for (int i = 0; i < np; i++) {
        
        p[i].Psi = 0.0;
        p[i].Theta = 0.0;
        xmean += p[i].x;
        ymean += p[i].y;
    }
    xmean /= np;
    ymean /= np;
    ensemble.c6 = 0.0;

    for (int i = 0; i < np; i++) {
        if (!p[i].op_neighbor.empty()){
            p[i].op_neighbor.clear();
        }
        for (int j = 0; j < np; j++) {
            if (i != j && Distance2D(p[i],p[j]) < rmin) {
                p[i].op_neighbor.push_back(j);
                degree = atan((p[j].y-p[i].y)/(p[j].x-p[i].x));
                p[i].Theta += degree;
                psir[i] += cos(6 * degree);
                psii[i] += sin(6 * degree);
            }  
        } 
    }
    
    // normalize local psi6 and theta
    for (int i = 0; i < np; i++){
        if (!p[i].op_neighbor.empty()) {
            psir[i] /=  p[i].op_neighbor.size();
            psii[i] /=  p[i].op_neighbor.size();
            p[i].Psi = sqrt(psir[i]*psir[i] + psii[i]*psii[i]);
            p[i].Theta *= (180.0/pi/p[i].op_neighbor.size());
            p[i].Theta += 30;
        }
        p[i].op_neighbor.push_back(i);
        sort(p[i].op_neighbor.begin(),p[i].op_neighbor.end());
    }
    
    // calculate local c6
    for (int i = 0; i < np; i++) {
        for (int j = 0; j < np; j++) {
            if (i != j && Distance2D(p[i],p[j]) < rmin) {
                double numer = psir[i] * psir[j] + psii[i] * psii[j];
                double temp = psii[i] * psir[j] - psii[j] * psir[i];
                double denom = sqrt(numer * numer + temp*temp);
                double testv = numer / denom;
                if (testv >= 0.32) {
                    con[i] += 1;
                }
            }
        }
    }
    
    // ensemble psi6, c6, and Rg
    for (int i = 0; i < np; i++) {
        accumpsi6r += psir[i];
        accumpsi6i += psii[i];
        p[i].c6 = con[i]/6;
        ensemble.c6 += con[i];
        rgmean += (p[i].x - xmean)*(p[i].x - xmean);
        rgmean += (p[i].y - ymean)*(p[i].y - ymean);
    }
    ensemble.psi6 = sqrt(accumpsi6r * accumpsi6r + accumpsi6i * accumpsi6i)/np;
    ensemble.c6 /= (np*5.6);
    ensemble.rg = sqrt(rgmean/np);
    
    // the position of particle is projected on the two axes directons of field
    double ang1[2] = {cos(22.5*pi/180),sin(22.5*pi/180)};
    double ang2[2] = {cos(112.5*pi/180),sin(112.5*pi/180)};
    
    // principle moments of particles along two axes
    ensemble.Ix = 0;
    ensemble.Iy = 0;
    for (int i = 0; i < np; i++){
        double mo = (p[i].x - xmean)*ang1[0] + (p[i].y - ymean)*ang1[1];
        ensemble.Iy += mo*mo/np/a/a;
        mo = (p[i].x - xmean)*ang2[0] + (p[i].y - ymean)*ang2[1];
        ensemble.Ix += mo*mo/np/a/a;
    }
}

void BrownianDynamicEngine::ImageAna(){
    vector<vector<int>> GB;
    vector<int> Check_list;
    
    // categorize particles
    for (int i = 0; i < np; i++){
        if (p[i].op_neighbor.size() < 5){
            p[i].type = PERIPHERAL;// peripheral
        } else if(p[i].Psi < 0.95 || p[i].c6 < 0.9){
//        } else if(p[i].Psi < 0.97 || p[i].op_neighbor.size() < 7){
            p[i].type = GRAIN;// peripheral
        } else {
            p[i].type = CRYSTAL;// peripheral
        }
    }
    
    // identify particles inside of peripheral particles
    for (int i = 0; i < np; i++){
        if(p[i].type != PERIPHERAL){
            for (int pt:p[i].op_neighbor){
                if (p[pt].type == PERIPHERAL){
                    p[i].type = INNER_RIM;
                }
            }
        }
    }
    
    // if a peripheral particle is not connected to other peripheral particles,
    // redefine it as grain particle
//    Check_list.clear();
//    for (int i = 0; i < np; i++){
//        if (p[i].type == PERIPHERAL || p[i].type == INNER_RIM){
//            Check_list.push_back(i);
//        }
//    }
//    if (!Check_list.empty()){
//        vector<vector<int>> rim = Connectivity(Check_list);
//        for (int i = 1; i < rim.size(); i++){
//            for (int pt: rim[i]){
//                p[pt].type = GRAIN;
//            }
//        }
//    }
    
    // group grain particles by their connectivity, define the largest cluster 
    // as the major (only) grain boundary, and redefine others into crystal
    Check_list.clear();
    for (int i = 0; i < np; i++){
        if (p[i].type == GRAIN){Check_list.push_back(i);}
    }
    if (!Check_list.empty()){
        GB = Connectivity(Check_list);
        for (int i = 1; i < GB.size(); i++){
//            if (GB[i].size() <= 5) {
                for (int pt : GB[i]){
                    p[pt].type = CRYSTAL;
                }
//            }
        }
    }
    
    // fit the grain boundary by linear regression of major grain particles.
    // the range of fit is [0,179], and 180 if no grain exists
    double fac1, fac2;
    double x_GB(0.0), y_GB(0.0);
    double currError(0);
    vector <double> fullError;
    if (GB.empty()) {
        ensemble.GB_fit = pi;
    } else {
        for (int pt: GB[0]){
            p[pt].type == GRAIN;
            x_GB += p[pt].x/GB[0].size();
            y_GB += p[pt].y/GB[0].size();
        }
        
        // fit the orientation by least square method
        fullError.clear();
        for (int i = 0; i < 180; i++){
            currError = 0.0;
            for (int pt : GB[0]) {
                fac1 = fabs(p[pt].x*tan(i*pi/180)-p[pt].y - tan(i*pi/180)*x_GB + y_GB);
                fac2 = sqrt(tan(i*pi/180)*tan(i*pi/180) + 1);
                currError += pow((fac1/fac2),2.0);
            }
            fullError.push_back(sqrt(currError));
        }
        vector<double>::iterator itr = min_element(fullError.begin(), fullError.end());
        
        // the grain boundary is averaged over past 5 values
        if (GB_queue.size()>=5) GB_queue.erase(GB_queue.begin());       
        GB_queue.push_back(distance(fullError.begin(), itr));
        vector<double> GB_sort(GB_queue.begin(), GB_queue.end());
        sort(GB_sort.begin(), GB_sort.end());
        ensemble.GB_fit = GB_sort[floor(GB_sort.size()/2)];
    }
}

vector<vector<int>> BrownianDynamicEngine::Connectivity(vector<int>& wait_list){
    vector<int> temp_list;
    int size_new, size_old;
    vector<vector<int>> group, Ordered_group;

    while (!wait_list.empty()){
        // if any particle remained in the wait list, create a new group for 
        // the first particle
        if (!temp_list.empty()) temp_list.clear();
        temp_list.push_back(wait_list[0]);
        do {
            size_old = temp_list.size();
            for (int pt:temp_list){
                for (int j:p[pt].op_neighbor){
                    if (find(wait_list.begin(), wait_list.end(), j) != wait_list.end()
                            && find(temp_list.begin(), temp_list.end(), j) == 
                            temp_list.end()){
                        temp_list.push_back(j);
                    }
                }
            }
            size_new = temp_list.size();
        } while(size_new != size_old);
        sort(temp_list.begin(),temp_list.end());
        group.push_back (temp_list);
        for (int pt : temp_list){
            wait_list.erase(remove(wait_list.begin(),wait_list.end(),pt),wait_list.end());
        }
    }
    
    
    // sort final group list
    int largest_group;
    while (!group.empty()){
        int group_size = 0;
        for (int i = 0; i < group.size(); i++){
            if (group[i].size() >= group_size){
                group_size = group[i].size();
                largest_group = i;
            }
        }
        Ordered_group.push_back(group[largest_group]);
        group.erase (group.begin()+largest_group);
    }
    return Ordered_group;
}



/******************************* OTHER OPERATIONS *******************************/

void BrownianDynamicEngine::OutputOrderParameter(ostream& os) {
    os << fixed;
    os << setprecision(2);
    os << setw(4);
    os << (double) elapsedTime/10 << "\t";
    os << fieldStrength << "\t";
    os << fieldDirection << "\t";
    os << ensemble.psi6 << "\t";
    os << ensemble.c6 << "\t";
    os << ensemble.rg/a << "\t";
    os << ensemble.GB_fit << "\t";
    os << ensemble.Ix << "\t";
    os << ensemble.Iy << endl;
}

void BrownianDynamicEngine::OutputTrajectory(ostream& os) {
    os << fixed;
    os << setprecision(2);
    os << setw(7);
    for (int i = 0; i < np; i++) {
        os << i << "\t";
        os << p[i].x/a << "\t";
        os << p[i].y/a << "\t";
        os << p[i].z/a << "\t";
        os << p[i].type << "\t";
        os << p[i].coordinate_number << endl;
    }
}

void BrownianDynamicEngine::ReadE(const std::string filename, int EPFlag) {
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
            linestream >> dETable[j][i][0];
            linestream >> dETable[j][i][1];
            linestream >> dETable[j][i][2];
            linestream >> dETable[j][i][3];
            }
        }
        is.close();
    }
}

void BrownianDynamicEngine::Readxyz(string filename) {
    if (filename.compare("random") == 0){ // randomize particle positions
        RandomizeConfig();
        return;
    } else if (filename.compare("library600") == 0 || 
            filename.compare("library300") == 0 || 
            filename.compare("library900") == 0){ // select from a set of configs
        stringstream configID;
        int config_file = repCycle%150;
        configID << config_file;
        filename = configPath + filename + "/sphere" + configID.str() + ".txt";
    } else { // read a specific configuration given by filename
        filename = configPath + filename;
    }
    
    ifstream is;
    is.open(filename.c_str());
    assert(is.is_open());
    string line;
    double dum;
    for (int i = 0; i < np; i++) {
        getline(is, line);
        std::stringstream linestream(line);
        linestream >> dum;
        linestream >> p[i].x;
        linestream >> p[i].y;
        linestream >> p[i].z;
        p[i].x *= a;
        p[i].y *= a;
        p[i].z *= a;
    }
    is.close();
}

void BrownianDynamicEngine::ReadDiffusivity() {
    const int rgbin(25), rbin(50); // rg bin step size
    double dum;
    string line, fileName("./library/Diffusivity.txt");
    ifstream is;
    
    is.open(fileName.c_str());
    assert(is.is_open());
    for (int i = 0; i < rbin * rgbin; i++) {
        getline(is, line);
        std::stringstream linestream(line);
        linestream >> dum;
        linestream >> dum;
        linestream >> dum;
        linestream >> dum;
        linestream >> *(&(dssarray[0][0])+i);
        linestream >> *(&(dsscount[0][0])+i);
    }
    is.close();
}

void BrownianDynamicEngine::getDefectedState(string fieldFileName, double fieldStrength_){
    stringstream repCycleOS;
    currCycle = repCycle++;
    repCycleOS<< currCycle;
    do {
        elapsedTime = 0;
        if (trajOs.is_open()) trajOs.close();
        trajOs.open("./traj/xyz" + repCycleOS.str() + ".dat"); 
        if (opOs.is_open()) opOs.close();
        opOs.open("./traj/op" + repCycleOS.str() + ".dat");
        RandomizeConfig();
        for (int t = 0; t < 50; t++){
            coreBD(fieldFileName, fieldStrength_, 1);
            if (getPsi6() > 0.6) {
                break; 
            }
        }
    } while (getPsi6() > 0.6);
    trajOs.close();
    trajOs.open("./traj/xyz" + repCycleOS.str() + ".dat"); 
    opOs.close();
    opOs.open("./traj/op" + repCycleOS.str() + ".dat");
//    OutputTrajectory(trajOs);
}

void BrownianDynamicEngine::RandomizeConfig(){ 
    vector<int> idxList;
    int randMeshGrid = 180;
    int idx;
    for (int i = 0; i < np; i++){
        do {
            idx = (int) floor(rand_uniform(gen) * randMeshGrid * randMeshGrid);
        } while (find(idxList.begin(), idxList.end(), idx) != idxList.end());
        idxList.push_back(idx);
        p[i].x = (double) ((int) idx/randMeshGrid) - randMeshGrid / 2;
        p[i].y = (double) ((int)idx % randMeshGrid - randMeshGrid / 2);
        p[i].x += (double) rand_uniform(gen) * 2 - 1;
        p[i].y += (double) rand_uniform(gen) * 2 - 1;
        p[i].x *= 1000;
        p[i].y *= 1000;
    }
    
    for (int i = 0; i < np; i++) {
        for (int j = i+1; j < np; j++){
            double rijsep = Distance2D(p[i],p[j]);
            if (rijsep < rcut){
                double xij = p[j].x - p[i].x;
                double yij = p[j].y - p[i].y;
                if (rijsep < 2 * a + 40 ){
                    p[i].x -= xij * (2*a + 40 - rijsep) /sqrt(xij * xij + yij * yij);
                    p[i].y -= yij * (2*a + 40 - rijsep) /sqrt(xij * xij + yij * yij);
                    p[j].x += xij * (2*a + 40 - rijsep) /sqrt(xij * xij + yij * yij);
                    p[j].y += yij * (2*a + 40 - rijsep) /sqrt(xij * xij + yij * yij);
                }
            }
        }
    }
}

void BrownianDynamicEngine::CalDss() {
    const double rgdsmin = 22250;
    const double delrgdsmin = -250;
    const double dssmin = 0.15;
    const double dssmax = 0.50;
    const double distmin = 0.0;
    const double deldist = 1400;
    // calculate rg
    double xmean(0), ymean(0), rgmean(0), rg;
    for (int i = 0; i < np; i++) {
        xmean += p[i].x/np;
        ymean += p[i].y/np;
    }

    for (int i = 0; i < np; i++) {
        rgmean += (p[i].x - xmean) * (p[i].x - xmean);
        rgmean += (p[i].y - ymean) * (p[i].y - ymean);
    }
    rg = sqrt(rgmean/np);
    int rgbinindex = (int) ((rg - rgdsmin) / delrgdsmin) + 1;
    if (rgbinindex <= 0) {rgbinindex = 1;}
    
    // approximate diffusivity by rg and partilce position
    for (int i = 0; i < np; i++) {
        double disttemp = sqrt(pow(p[i].x - xmean, 2) + pow(p[i].y - ymean, 2));
        int distbinindex = (int) ((disttemp - distmin) / deldist) + 1;
        if (rgbinindex >= 1 && rgbinindex <= rgdssbin) {
            if (distbinindex >= 1 && distbinindex <= distdssbin) {
                if (dsscount[rgbinindex-1][distbinindex-1] >= 1) {
                    p[i].D  = dssarray[rgbinindex - 1][distbinindex - 1];
                } else {
                    p[i].D  = dssmax;
                }
            } else {
                p[i].D = dssmax;
            } 
        } else if (rgbinindex >= rgdssbin) {
            p[i].D = dssmin;
        }
    }
}