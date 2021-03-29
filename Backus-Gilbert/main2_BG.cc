#include <iostream> 
#include <math.h> 
#include <fstream>
#include <string>
#include <vector>
#include "H5Cpp.h"
using namespace H5;
using namespace std; 
#ifndef M_PI
   #define M_PI 3.14159265358979323846
#endif

//definition of the kernel (equation 1)
double Kernel(const double& omega_i, const double& t_i){
    return cosh(omega_i*(abs(t_i)-beta/2))/sinh(omega_i*beta/2);
}



//try to compute W with rescaled kernel but it doesn't change the fact that W_reg is singular
double Rescaled_Kernel(const double& omega_i, const double& t_i){
    return omega_i*cosh(omega_i*(abs(t_i)-beta/2))/sinh(omega_i*beta/2);
    //return pow((omega_i/T),2)/tanh(omega_i/(2*T))*cosh(omega_i*(abs(t_i)-beta/2))/sinh(omega_i*beta/2);
    //tanh(omega_i/(2*T))*cosh(omega_i*(abs(t_i)-beta/2))/sinh(omega_i*beta/2);
}



//Analytic defintion of the spectral function (equation 28)
double Analytic_rho(const double& omega_i){
    return 2*omega_i*Gamma/(M_PI*pow((omega_i*omega_i-Gamma*Gamma-mass*mass),2)+4*omega_i*omega_i*Gamma*Gamma);
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

//computation of the analytic correlator using the definition of rho (equation 28)
void analytic_G_E(vector<double>& G_E, const vector<double>& y, const vector<double>& t){
    for(int i=0;i<t.size();i++){
        vector<double>f_G_E(M);
        for (int j=0;j<M;j++){
        	if(y[j]==1){
            	f_G_E[j]=2.0*Gamma/(beta/2*M_PI*pow((pow(mass,2)+pow(Gamma,2)),2));
          }else if(y[j]==0){
          	f_G_E[j]=0;
          }else{
          	//integrant (equation 1)
          	f_G_E[j]=1.0/pow(y[j],2)*Kernel((1-y[j])/y[j],t[i])*Analytic_rho((1-y[j])/y[j]);
          }
        }
        G_E[i]=simpsons_(f_G_E,dy);
    }
}




//Computation of W (equation 16) using y as integration variable (integration from 0 to 1)
void Compute_W(vector<vector<double>>& W, const double& omega_i, const vector<double>& t, const vector<double>& y){
    //analogy with (equation 16): omega_i->w_0, i->tau_i, j->tau_j, k->omega(integration variable)
    for(int i=0;i<t.size();i++){
        for(int j=0;j<t.size();j++){
            //the matrix is symetric so we only compute the upper trianglar value, the others are set to zero
        	if(j<i+1){
            	vector<double> f_W(N);
            	for(int k=0;k<M;k++){
            		if(y[k]==1.0){
               			f_W[k]=pow(omega_i,2)/pow(beta/2,2);
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


//compuation of R (equation 17)
void Compute_R(vector<double>& R, const vector<double>& t, const vector<double>& y){
    //analogy with (equation 7):i->tau_i, j->omega (integration variable)
    for(int i=0;i<t.size();i++){
        vector<double> f_R(N);
        for(int j=0;j<M;j++){
        	if(y[j]==1){
            	f_R[j]=2.0/beta;
        	}else if (y[j]==0){
        		f_R[j]=0;        		
        	}else{
        		f_R[j]=1.0/pow(y[j],2)*Rescaled_Kernel((1-y[j])/y[j],t[i]);
        	}
        }
        R[i]=simpsons_(f_R,dy);
    }
}

void WriteInt(const int& N, const string& filename){
    string const N_data(filename);
    ofstream N_Flux(N_data.c_str());
    N_Flux<<N;
}

void WriteVec(const vector<double>& vec, const string& filename){
    string const data_vec(filename);
    ofstream vec_Flux(data_vec.c_str());
    for(int i=0; i<vec.size();i++){
        if(i==vec.size()-1){
            vec_Flux<<vec[i];
        }else{
            vec_Flux<<vec[i]<<endl;
        }
    }
}

int WriteW(const H5std_string FILE_NAME, const vector<vector<double>>& W, const int & size){
    const int   RANK = 2;
    const H5std_string  DATASET_NAME( "DoubleArray" );
    try{
      /*
       * Turn off the auto-printing when failure occurs so that we can
       * handle the errors appropriately
       */
      Exception::dontPrint();

      /*
       * Create a new file using H5F_ACC_TRUNC access,
       * default file creation properties, and default file
       * access properties.
       */
      H5File file( FILE_NAME, H5F_ACC_TRUNC );

      /*
       * Define the size of the array and create the data space for fixed
       * size dataset.
       */
      hsize_t     dimsf[2];              // dataset dimensions
      dimsf[0] = size;
      dimsf[1] = size;
      DataSpace dataspace( RANK, dimsf );

      /*
       * Define datatype for the data in the file.
       * We will store little endian INT numbers.
       */
      DoubleType datatype( PredType::NATIVE_DOUBLE );
      datatype.setOrder( H5T_ORDER_LE );

      /*
       * Create a new dataset within the file using defined dataspace and
       * datatype and default dataset creation properties.
       */
      DataSet dataset = file.createDataSet( DATASET_NAME, datatype, dataspace );

      /*
       * Write the data to the dataset using default memory space, file
       * space, and transfer properties.
       */
      dataset.write( W, PredType::NATIVE_DOUBLE );
   }  // end of try block

   // catch failure caused by the H5File operations
   catch( FileIException error ){
      error.printErrorStack();
      return -1;
   }

   // catch failure caused by the DataSet operations
   catch( DataSetIException error ){
      error.printErrorStack();
      return -1;
   }

   // catch failure caused by the DataSpace operations
   catch( DataSpaceIException error ){
      error.printErrorStack();
      return -1;
   }

   // catch failure caused by the DataSpace operations
   catch( DataTypeIException error ){
      error.printErrorStack();
      return -1;
   }

   return 0;  // successfully terminated
}



int main(int argc, char *argv[]){
	double omega_i=atoi(argv[1]);
	vector <double> t(N);
    for(int i=0; i<N;i++){
        t[i]=min_t+i*dt;
    }
    for(int i=0; i<N;i++){
        if (t[i]==0){
            t.erase(t.begin()+i);
        }
    }
    vector<double> y_integrate(M);
    vector<double> omega(M);
    for (int i=0;i<M;i++){
    	omega[i]=i*domega;
        y_integrate[i]=i*dy;
    }


    //intilisation of all physical quantities we want
    vector<double>G_E(t.size());
    //double* G_E=new double[N];
    analytic_G_E(G_E,y_integrate,t);
    vector<double>rho_analytic(M);
    for(int alpha=0;alpha<M;alpha++){
        rho_analytic[alpha]=Analytic_rho(omega[alpha]);
        vector<vector<double>> W (t.size(),vector<double>(t.size()));
        //double** W=create_W();
        Compute_W(W,omega[alpha],t,y_integrate);
        H5std_string filename="data_W_"+to_string(N)+"_"+to_string(alpha)+".txt";
        error=WriteW(filename,W,Nt)
        /*
        string filename="data_W_"+to_string(N)+"_"+to_string(alpha)+".txt";
        string const data_W(filename);
    	ofstream W_Flux(data_W.c_str());
    	for(int i=0;i<t.size();i++){
            for(int j=0;j<t.size();j++){
                if (i==t.size()-1 and j==t.size()-1){
                    W_Flux<<W[i][j];
                }else{
                    W_Flux<<W[i][j]<<endl;
                }               
            }
        }
        */
    }
    vector<double> R(t.size());
    Compute_R(R,t,y_integrate);
    
    //Write all data in a file
    string N_string=to_string(N);
    string M_string=to_string(M);
    WriteInt(N,"N.txt");
    WriteInt(M,"M.txt");
    WriteVec(rho_analytic,"data_rho_analytic_"+M_string+".txt");
    WriteVec(G_E,"data_G_E_"+N_string+".txt");
    WriteVec(t,"data_t_"+N_string+".txt");
    WriteVec(omega,"data_omega_"+M_string+".txt");
    WriteVec(R,"data_R_"+N_string+".txt");

}