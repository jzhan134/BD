/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include "OrderParameter.hpp"
void OrderParameter::CalOp() {
    const double rmin = 2.64*a; // threshold for OP calculation
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
            if (i != j && Distance(p[i],p[j]) < rmin) {
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
//            cout << p[i].Theta << endl;
        }
        p[i].op_neighbor.push_back(i);
        sort(p[i].op_neighbor.begin(),p[i].op_neighbor.end());
    }
    
    // calculate local c6
    for (int i = 0; i < np; i++) {
        for (int j = 0; j < np; j++) {
            if (i != j && Distance(p[i],p[j]) < rmin) {
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
    double ang1[2] = {cos(22.5*pi/180),sin(22.5*pi/180)};
    double ang2[2] = {cos(112.5*pi/180),sin(112.5*pi/180)};
    double Iy, Ix;
    
    // ensemble long/short axis ratio
    Iy = 0;
    Ix = 0;
    for (int i = 0; i < np; i++){
        double mo = (p[i].x - xmean)*ang1[0] + (p[i].y - ymean)*ang1[1];
        Iy += mo*mo/np/a/a;
        mo = (p[i].x - xmean)*ang2[0] + (p[i].y - ymean)*ang2[1];
        Ix += mo*mo/np/a/a;
    }
    ensemble.Ix = Ix;
    ensemble.Iy = Iy;
}

void OrderParameter::ImageAna(){
    vector<vector<int>> GB;
    vector<int> Check_list;
    
    
    // categorize particles
    for (int i = 0; i < np; i++){
        if (p[i].op_neighbor.size() < 6){
            p[i].type = -1;// rim
        } else if(p[i].Psi < 0.95 || p[i].op_neighbor.size() < 7){
            p[i].type = 0;// grain boundary
        } else {
            p[i].type = 1;// crystalline
        }
    }
    
    // inner rim particles
    for (int i = 0; i < np; i++){
        if(p[i].type != -1){
            for (int pt:p[i].op_neighbor){
                if (p[pt].type == -1){
                    p[i].type = -2;
                }
            }
        }
    }
    
    // group unconnected rim particles into grain boundary
    Check_list.clear();
    for (int i = 0; i < np; i++){
        if (p[i].type == -1 || p[i].type == -2){
            Check_list.push_back(i);
        }
    }
    if (!Check_list.empty()){
        vector<vector<int>> rim = Connectivity(Check_list);
        for (int i = 1; i < rim.size(); i++){
            for (int pt: rim[i]){
                p[pt].type = 0;
            }
        }
    }
    
    
    Check_list.clear();
    for (int i = 0; i < np; i++){
        if (p[i].type == 0){Check_list.push_back(i);}
    }
    if (!Check_list.empty()){
        GB = Connectivity(Check_list);
    }

    // Group crystalline particles into domains
    Check_list.clear();
    for (int i = 0; i < np; i++){
        if (p[i].type == 1){
            Check_list.push_back(i);
        }
    }
    vector<double> domain_temp;
    if (!Check_list.empty()){
        vector<vector<int>> Domain = Connectivity(Check_list);
        for (int i = 0; i < Domain.size() && i < 2; i++){
            domain_temp.clear();
            for (int pt : Domain[i]){
                p[pt].type = i+1;
                domain_temp.push_back(p[pt].Theta);
            }
            std::sort(domain_temp.begin(),domain_temp.end());
            ensemble.Domain_ori[i] = domain_temp[round(domain_temp.size()/2)];
        }
    }
    
    // fit the grain boundary, the default is 180, otherwise it's within [0,179]
    double fac1, fac2;
    double x_GB(0.0), y_GB(0.0);
    vector <double> dist_sum;
    double local_dist_sum;
    if (!GB.empty()){
        // find center of grain particles
        for (int pt: GB[0]){
            p[pt].type = -3;
            x_GB += p[pt].x/GB[0].size();
            y_GB += p[pt].y/GB[0].size();
        }
        // fit the orientation by least square method
        dist_sum.clear();
        for (int i = -60; i < 120; i++){
            local_dist_sum = 0.0;
            for (int pt : GB[0]) {
                fac1 = fabs(p[pt].x*tan(i*pi/180)-p[pt].y - tan(i*pi/180)*x_GB + y_GB);
                fac2 = sqrt(tan(i*pi/180)*tan(i*pi/180) + 1);
                local_dist_sum += pow((fac1/fac2),2.0);
            }
            dist_sum.push_back(sqrt(local_dist_sum));
        }
        vector<double>::iterator itr = min_element(dist_sum.begin(), dist_sum.end());
        if (GB_queue.size()>=5) GB_queue.erase(GB_queue.begin());
        GB_queue.push_back(distance(dist_sum.begin(), itr));
        vector<double> GB_sort(GB_queue.begin(), GB_queue.end());
        sort(GB_sort.begin(), GB_sort.end());
        ensemble.GB_fit = GB_sort[floor(GB_sort.size()/2)]-60;

    }
}

vector<vector<int>> OrderParameter::Connectivity(vector<int>& wait_list){
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