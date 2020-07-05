/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
//
//const int L = 287.0; // radius of depletants
//
//switch (OS_opt){
//    case 0:
//        Os_pressure = 0.00000;
//        break;
//    case 1:
//        Os_pressure = 0.00767;
//        break;
//    case 2:
//        Os_pressure = 0.01490;
//        break;
//    case 3:
//        Os_pressure = 0.02198;
//        break;
//    case 4:
//        Os_pressure = 0.02896;
//        break;
//}

//const double DG = 100000;
//void BrownianDynamicEngine::analyticalField(double Ex[], double Ey[], double dExdx[], 
//        double dExdy[], double dEydx[], double dEydy[], double dE2x[], double dE2y[]){
//    double STEP = 1e-3;
//    double EMAG0, EMAGx, EMAGy;
//    for (int i = 0 ; i < np; i++) {
//        Ex[i] = -4*p[i].x/DG;
//        Ey[i] = 4*p[i].y/DG;
//        dExdx[i] = 4/DG;
//        dExdy[i] = 0;
//        dEydx[i] = 0;
//        dEydy[i] = -4/DG;
//        
//        EMAG0 = EMAG(p[i].x, p[i].y);
//        EMAGx = EMAG(p[i].x + STEP, p[i].y);
//        EMAGy = EMAG(p[i].x, p[i].y + STEP);
//        dE2x[i] = (EMAGx * EMAGx - EMAG0 * EMAG0) / STEP;
//        dE2y[i] = (EMAGy  * EMAGy  - EMAG0  * EMAG0 ) / STEP;
//
//    }
//}

//        ensemble.GB_fit = fabs(90 - GB_sort[floor(GB_sort.size()/2)]);


    // calculate domain angle with respect to grain boundary orientation
//    double DM_temp[] = {ensemble.Domain_ori[0]-60, ensemble.Domain_ori[0],
//        ensemble.Domain_ori[0]+60, ensemble.Domain_ori[0] + 120,
//        ensemble.Domain_ori[0]+180};
//    double min_Theta = 31;
//    for (double DM : DM_temp){
//        if (abs(DM - ensemble.GB_fit) <= abs(min_Theta)){
//            min_Theta = DM - ensemble.GB_fit;
//        }
//    }
//    ensemble.Domain_ori[0] = min_Theta;
//    double DM_temp2[] = {ensemble.Domain_ori[1]-60, ensemble.Domain_ori[1],
//        ensemble.Domain_ori[1]+60, ensemble.Domain_ori[1] + 120,
//        ensemble.Domain_ori[1]+180};
//    min_Theta = 31;
//    for (double DM : DM_temp2){
//        if (abs(DM - ensemble.GB_fit) <= abs(min_Theta)){
//            min_Theta = DM - ensemble.GB_fit;
//        }
//    }
//    ensemble.Domain_ori[1] = min_Theta;

    /*
     * calculate the relative domain size
     * relative domain size is defined as follow: first, a vector is defined as
     * from the grain boundary center to one direction of the grain boundary;
     * next another vector is defined as from grain boundary center to each 
     * particle. Last, the cross product between the two vectors is calculated
     * the particles with opposite signs are grouped as different groups.
     */
     
//    int sz = 0;
//    for (int pt = 0; pt < np; pt++){
//        if  (cos(ensemble.GB_fit*pi/180)*(p[pt].y - y_GB) 
//                    - sin(ensemble.GB_fit*pi/180)*(p[pt].x - x_GB) > 0) sz++;
//    }
//    sz = (sz > (300-sz))? sz : 300 - sz;
//    if (SZ_queue.size()>=5) SZ_queue.erase(SZ_queue.begin());
//    SZ_queue.push_back(sz);
//    vector<double> SZ_sort(SZ_queue.begin(), SZ_queue.end());
//    sort(SZ_sort.begin(), SZ_sort.end());
//    ensemble.sz = SZ_sort[floor(SZ_sort.size()/2)];

//vector<int> SZ_queue;
//double Os_pressure;

//if (rijsep > 2*a && (rijsep < (2*a + 2*L))){
//                        FOS = (4.0/3.0)* Os_pressure*M_PI*
//                                (-0.75*(a+L)*(a+L)*1e-18 + 
//                                0.1875*rijsep*rijsep*1e-18)*1e9;
//                    }

//double BrownianDynamicEngine::EMAG(double RX, double RY) const {        
//// 1. Old expression of isotropic field  
////    double RT = sqrt(RX * RX + RY * RY);
////    double result = 4 * RT / DG;
////    double correctfactor = 2.081e-7*pow(RT/1000.0,4) - 1.539e-9*pow(RT/1000.0,3) + 
////            8.341e-5*pow(RT/1000, 2) + 1.961e-5*(RT/1000) + 1.028;
////    result *= correctfactor;
////    return result;
//    
//// 2. New version of anisotropic field
//    double result = 4.0*sqrt((RX/DG)*(RX/DG) + (RY/DG)*(RY/DG));
//    return result;
//}

//                    if (rijsep < 2*(a+10)){
//                        Fpp = Fhw*(2*(a+10)-rijsep)/20;
//                    } else {
//                        Fpp = 1e18*kb*(tempr+273)*kappa*pfpp*exp(-kappa*(rijsep-2.0*a)/a)/a;
//                    }
//                    Fpp = (rijsep < 2*a) ? Fhw : 1e18*kb*(tempr+273)*kappa*
//                            pfpp*exp(-kappa*(rijsep-2.0*a)/a)/a;