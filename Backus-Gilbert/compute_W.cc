#include <iostream> 
#include <math.h> 
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <H5Cpp.h>
//#include <HDF5pp.h>
using namespace H5;
using namespace std; 
#ifndef M_PI
   #define M_PI 3.14159265358979323846
#endif

const int N=32;
const int M=32;
const double T=2.0;
const double min_t=-0.25;
const double max_t=0.25;
const double dt=(max_t-min_t)/N;
const double Gamma=100.0;
const double Beta=max_t-min_t;
const double mass=300.0;
const double dy=1.0/(M-1);


//definition of the kernel (equation 1)
double Kernel(const double& omega_i, const double& t_i){
    return (exp(omega_i*(abs(t_i)-Beta))+exp(-omega_i*abs(t_i)))/(1.0-exp(-omega_i*Beta));
}



//try to compute W with rescaled kernel but it doesn't change the fact that W_reg is singular
double Rescaled_Kernel(const double& omega_i, const double& t_i){
    return omega_i*Kernel(omega_i,t_i);
}

//Simpson method to approximate the integral
double simpsons_(const vector<double>& fx, const double& dx){ 
    double res = 0; 
    for (int i = 0; i<M; i++) { 
        if (i == 0 || i == N) 
            res += fx[i]; 
        else if (i % 2 != 0) 
            res += 4 * fx[i]; 
        else
            res += 2 * fx[i]; 
    } 
    res = res * (dx/ 3); 
    return res; 
}


//Computation of W (equation 16) using y as integration variable (integration from 0 to 1)
void Compute_W(vector<vector<double>>& W, const double& omega_i, const vector<double>& t, const vector<double>& y){
    //analogy with (equation 16): omega_i->w_0, i->tau_i, j->tau_j, k->omega(integration variable)
    for(unsigned int i=0;i<t.size();i++){
        for(unsigned int j=0;j<t.size();j++){
            //the matrix is symetric so we only compute the upper trianglar value, the others are set to zero
        	if(j<i+1){
            	vector<double> f_W(N);
            	for(int k=0;k<M;k++){
            		if(y[k]==1.0){
               			f_W[k]=pow(omega_i,2)/pow(Beta/2,2);
            		}else if (y[k]==0){
               			f_W[k]=0;
            		}else{
               			//integrant (equation 16)
           				f_W[k]=1.0/pow(y[k],2)*Rescaled_Kernel((1.0-y[k])/y[k],t[i])*Rescaled_Kernel((1.0-y[k])/y[k],t[j])*pow(omega_i-(1.0-y[k])/y[k],2);
            		}
            	}
            	//integration over omega (index k)
            	W[i][j]=simpsons_(f_W,dy);
            }else{
                W[i][j]=0.0;
            }  
        }
    }
}


void WriteW(const H5std_string FILE_NAME, const vector<vector<double>>& W, const int & taille){
    // Create HDF5 file and dataset
    H5File file(FILE_NAME, H5F_ACC_TRUNC);
    // dataset dimensions
    hsize_t dimsf[2];
    dimsf[0] = taille;
    dimsf[1] = taille;
    DataSpace dataspace(2, dimsf);


  	DataType datatype(H5::PredType::NATIVE_DOUBLE);
  	DataSet dataset = file.createDataSet("data", datatype, dataspace);

 
  	dataset.write(&W[0][0], H5::PredType::NATIVE_DOUBLE);


    dataset.close();
    dataspace.close();
    file.close();

}


int main(int argc, char *argv[]){
	istringstream ss(argv[1]);
	double omega_i;
	ss>>omega_i;
	vector <double> t(N);
    for(int i=0; i<N;i++){
        t[i]=min_t+i*dt;
    }
    /*
    for(int i=0; i<N;i++){
        if (t[i]==0){
            t.erase(t.begin()+i);
        }
    }
    */
    vector<double> y_integrate(M);
    for (int i=0;i<M;i++){
        y_integrate[i]=i*dy;
    }

    vector<vector<double>> W (t.size(),vector<double>(t.size()));
    Compute_W(W,omega_i,t,y_integrate);
    H5std_string filename="data_W_"+to_string(N)+"_"+to_string(omega_i)+".h5";
    WriteW(filename,W,t.size());
}