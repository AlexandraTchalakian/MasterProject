#coding:utf-8
from math import *
import random
import numpy as np
from scipy import integrate as integrate
from scipy.integrate import simps
from numpy import trapz
from decimal import localcontext
import matplotlib.pyplot as plt
import os.path
from matplotlib import rc


def Read_N():
	with open("N.txt",'r') as data_N:
		N=int(data_N.readline())
	return N

N=Read_N()
Nt=N
#N=650
M=300
gamma=100
T=2
beta=1/T
Lambda=1e-9



def Read_omega():
	omega=np.zeros(N)
	file_name="data_omega_"+str(N)+".txt"
	with open(file_name,'r') as data_omega:
		for i in range(0,N):
			omega[i]= float(data_omega.readline())

	return omega;

def Read_omega_integrate():
	omega_integrate=np.zeros(N)
	file_name="data_omega_integrate_"+str(N)+".txt"
	with open(file_name,'r') as data_omega_integrate:
		for i in range(0,N):
			omega_integrate[i]= float(data_omega_integrate.readline())

	return omega_integrate;

def Read_t():
	t=np.zeros(Nt)
	file_name="data_t_"+str(N)+".txt"
	with open(file_name,'r') as data_t:
		for i in range(0,Nt):
			t[i]= float(data_t.readline())

	return t;

def Read_G_E():
	G_E_analytic=np.zeros(Nt)
	file_name="data_G_E_"+str(N)+".txt"
	with open(file_name,'r') as data_G_E:
		for i in range(0,Nt):
			G_E_analytic[i]= float(data_G_E.readline())

	return G_E_analytic;

def Read_analytic_rho():
	rho_analytic=np.zeros(N)
	file_name="data_rho_analytic_"+str(N)+".txt"
	with open(file_name,'r') as data_rho_analytic:
		for i in range(0,N):
			rho_analytic[i]= float(data_rho_analytic.readline())

	return rho_analytic;



def Read_W():
	v=np.zeros(N*Nt**2)
	W=np.zeros((N,Nt,Nt))
	file_name="data_W_"+str(N)+".txt"
	with open(file_name,'r') as data_W:
		for  line in range (0,N*Nt**2):
			#read all columuns of the file
			v[line]=float(data_W.readline())
			#print(v[line])

	#Here we have to match the right matrix element with the right columun number, 
	#that why we add N at the end of each i loop. 
	n=0
	#Read the upper triangular part
	for alpha in range (0,N):#omega
		for i in range (0,Nt):#t_i
			for j in range (0,i+1):#t_j
				W[alpha][i][j]=v[j+n]	
			n+=Nt	

	#Computation of the other elements by symmetry		
	for alpha in range (0,N):#omega
		for i in range (0,Nt):#t_i
			for j in range (i+1,Nt):#t_j
				W[alpha][i][j]=W[alpha][j][i]


	"""
	#To verify the symmetry of W
	for i in range (0,N):
		for j in range (0,N):
			print(W[0][i][j],W[0][j][i])
	"""

	return W

def Read_R():
	R=np.zeros(Nt)
	file_name="data_R_"+str(N)+".txt"
	with open(file_name,'r') as data_R:
		for i in range(0,Nt):
			R[i]= float(data_R.readline())

	return R;


def Inverse_W(W):
	(U, S, VT) = np.linalg.svd(W)
	Sigma=np.diag(S)
	D=np.zeros((Nt,Nt))
	for i in range (0,Nt):
		D[i][i]=1.0/Sigma[i][i]
	D_tild=np.zeros((Nt,Nt))
	for i in range (0,Nt):
		D_tild[i][i]=Sigma[i][i]/(Sigma[i][i]**2+Lambda**2)
	UT=U.transpose()
	V=VT.transpose()
	return V.dot(D_tild.dot(UT));





#computation of the numerical spectral function
def Spectral_function_construction(W,R,G_E_analytic):

	spectral_function=np.zeros(N)
	#The spectral function is a function of omega, 
	#in the discret version this is traduced by a vector evaluted in each omega value (index alpha).
	for alpha in range (0,N):
		W_omega=np.zeros((Nt,Nt))
		for i in range (0,Nt):
			for j in range (0,Nt):
				W_omega[i][j]=W[alpha][i][j]
		
		(eigenvalues_W,eigenvectors_W)=np.linalg.eig(W_omega)
		W_Inverse=Inverse_W(W_omega)
		#Print of the determinant
		#print(np.linalg.det(W_Inverse))
		#Computation of q (equation 15)
		q=W_Inverse.dot(R)/(R.dot(W_Inverse.dot(R)))
		#computation of the spectral function (equation 13)
		spectral_function[alpha]=q.dot(G_E_analytic)

	return spectral_function



def main():

	#Read all data 
	G_E_analytic=Read_G_E()
	rho_analytic=Read_analytic_rho()
	omega=Read_omega()
	omega_integrate=Read_omega_integrate()
	t=Read_t()
	W=Read_W()
	R=Read_R()

	#computation of the spectral function
	spectral_function=Spectral_function_construction(W,R,G_E_analytic)


	rescaled_spectral_function=np.zeros(N)
	for i in range(0,N):
		#rescaled_spectral_function[i]=spectral_function[i]*omega_integrate[i]
		rescaled_spectral_function[i]=spectral_function[i]*omega[i]
		#rescaled_spectral_function[i]=spectral_function[i]*(omega[i]/T)**2/tanh(omega[i]/(2*T))


	plt.figure(1)
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	plt.plot(omega,rho_analytic)
	#plt.plot(omega,spectral_function)
	plt.plot(omega,rescaled_spectral_function)
	plt.ylabel(r'$\rho(\omega)$')
	plt.xlabel(r'$\omega$')
	plt.legend(('analytic','numeric'), loc='upper right')
	plt.show()


	plt.figure(2)
	plt.yscale("log")
	plt.plot(t,G_E_analytic)
	plt.show()




if __name__ == '__main__':
	main()	