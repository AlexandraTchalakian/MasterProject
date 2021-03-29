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

N=512
M=300
gamma=100
T=2
beta=1/T
Lambda=0.5


def coviante_matrix(x):
	cov_x=np.zeros((N,N))
	mean=0
	for i in range (0,N):
		mean+=x[i]
	mean/=N
	for i in range(0,N):
		for j in range (0,N):
			cov_x[i][j]=(x[i]-mean)*(x[j]-mean)

	return cov_x

def Integrate_G_E(omega_i,t_i):
	return Analytic_rho(omega_i)*Kernel(omega_i,t_i)



def Kernel(omega_i,t_i):
	return cosh(omega_i*(abs(t_i)-beta/2))/sinh(omega_i*beta/2)
	#return (omega_i/T)**2/tanh(omega_i/(2*T))*cosh(omega_i*(abs(t_i)-beta/2))/sinh(omega_i*beta/2)
	#return tanh(omega_i/(2*T))*cosh(omega_i*(abs(t_i)-beta/2))/sinh(omega_i*beta/2)
	#return omega_i*cosh(omega_i*(abs(t_i)-beta/2))/sinh(omega_i*beta/2)
	#return (omega_i/T)**4/(tanh(omega_i/(4*T))*tanh(omega_i*beta/T))*cosh(omega_i*(t_i-beta/2))/sinh(omega_i*beta/2)

def Analytic_rho(omega_i):
	return 2*omega_i*gamma/(pi*(omega_i**2-gamma**2-M**2)**2+4*omega_i**2*gamma**2)


def Analytic_G_E(omega,t):
	G_E_analytic=np.zeros(N)
	err_G_analytic=np.zeros(N)
	for i in range (0,N):
		f=np.zeros(N)
		for j in range (0,N):
			f[j]=Integrate_G_E(omega[j],t[i])
		G_E_analytic[i]= simps(f, omega)
		"""

		try:			
			(G_analytic[i],err_G_analytic[i])=integrate.quad(lambda x:Integrate_G_E(x,t[i]),0,np.inf)
			
			
		except OverflowError:
			print(omega[i],'----',t[i],'----',i)
		"""
	return G_E_analytic


def Backus_Gilbert(omega,t,G_E_analytic):

	
	nomfichier_W="data_W_"+str(N)+".txt"
	exist_W=True
	if not os.path.isfile(nomfichier_W):
		exist_W=False

	v=np.zeros(N**3)
	if exist_W:
		with open(nomfichier_W,'r') as fichier:
			for  line in range (0,N**3):
				v[line]=float(fichier.readline())
				#print(v[line])
				



	W=np.zeros((N,N,N))
	if not exist_W:
		with open(nomfichier_W,'w') as fichier:
			for alpha in range (0,N):#omega
				for i in range (0,N):#t_i
					for j in range (0,N):#t_j	
						if j<i+1:
							f_W=np.zeros(N)
							for k in range (0,N):#integration
								f_W[k]=Kernel(omega[k],t[i])*Kernel(omega[k],t[j])*(omega[alpha]-omega[k])**2
							W[alpha][i][j]=simps(f_W,omega)
						else:
							W[alpha][i][j]=0	
						if i==N-1 and alpha==N-1 and j==N-1:
							fichier.write(str(W[alpha][i][j]))
						else:
							fichier.write(str(W[alpha][i][j])+"\n")

		with open(nomfichier_W,'r') as fichier:
			for  line in range (0,N**3):
				v[line]=float(fichier.readline())
	
	
	n=0
	for alpha in range (0,N):#omega
		for i in range (0,N):#t_i
			for j in range (0,i+1):#t_j
				W[alpha][i][j]=v[j+n]	
			n+=N	

		
	for alpha in range (0,N):#omega
		for i in range (0,N):#t_i
			for j in range (i+1,N):#t_j
				W[alpha][i][j]=W[alpha][j][i]


	nomfichier_R="data_R_"+str(N)+".txt"
	exist_R=True
	if not os.path.isfile(nomfichier_R):
		exist_R=False

	R=np.zeros(N)
	if exist_R:
		with open(nomfichier_R,'r') as fichier:
			for i in range (0,N):
				R[i]=float(fichier.readline())

	if not exist_R:
		with open(nomfichier_R,'w') as fichier:
			for i in range (0,N):#t_i
				f_R=np.zeros(N)
				for j in range (0,N):#t_j	
					f_R[j]=Kernel(omega[j],t[i])
				R[i]=simps(f_R,omega)
				fichier.write(str(R[i])+"\n")			
		with open(nomfichier_R,'r') as fichier:
			for i in range (0,N):
				R[i]=float(fichier.readline())
			

	spectral_function=np.zeros(N)
	for alpha in range (0,N):
		W_omega=np.zeros((N,N))
		for i in range (0,N):
			for j in range (0,N):
				W_omega[i][j]=W[alpha][i][j]
		

		(eigenvalues,eigenvectors)=np.linalg.eig(W_omega)
		for i in range (0,N):
			print(eigenvalues[i])
		print(np.linalg.det((1-Lambda)*W_omega+Lambda*coviante_matrix(G_E_analytic)))
		W_reg_inverse=np.linalg.inv((1-Lambda)*W_omega+Lambda*coviante_matrix(G_E_analytic))
		q=W_reg_inverse.dot(R)/(R.dot(W_reg_inverse.dot(R)))
		spectral_function[alpha]=q.dot(G_E_analytic)

	return spectral_function






"""
def Backus_Gilbert_2(omega,t,G_E_analytic):
#Does not converge
	W=[[[0]*N]*N]*N
	R=np.zeros(N)
	for alpha in range (0,N):
		for i in range (0,N):
			for j in range (0,N):
				W[alpha][i][j]=integrate.quad(lambda x: Kernel(x,t[i])*Kernel(x,t[j])*(omega[alpha]-x)**2,0.1,np.inf)
			R[i]=integrate.quad(lambda x:Kernel(x,t[i]),0,np.inf)


	spectral_function=np.zeros(N)

	for alpha in range (0,N):
		W_omega=np.zeros((N,N))
		for i in range (0,N):
			for j in range (0,N):
				W_omega[i][j]=W[alpha][i][j]
		print(np.linalg.det((1-Lambda)*W_omega+Lambda*coviante_matrix(G_E_analytic)))
		W_reg_inverse=np.linalg.inv((1-Lambda)*W_omega+Lambda*coviante_matrix(G_E_analytic))
		q=W_reg_inverse.dot(R)/(R.dot(W_reg_inverse.dot(R)))
		spectral_function[alpha]=q.dot(G_E_analytic)
"""


def main():
	omega=np.linspace(0.001,1000,N)
	t=np.linspace(-0.25,0.25,N)
	rho_analytic=np.zeros(N)
	ker=np.zeros(N)
	inte=np.zeros(N)

	
	for i in range (0,N):	
		rho_analytic[i]+=Analytic_rho(omega[i])
		ker[i]+=Kernel(omega[i],t[i])
		inte[i]+=Integrate_G_E(omega[i],t[i])

	G_E_analytic=Analytic_G_E(omega,t)
	spectral_function=Backus_Gilbert(omega,t,G_E_analytic)

	plt.figure(1)
	plt.plot(omega,rho_analytic)
	plt.plot(omega,spectral_function)
	plt.show()

	plt.figure(2)
	plt.plot(omega,ker)
	plt.show()

	plt.figure(3)
	plt.plot(omega,inte)
	plt.show()

	plt.figure(100)
	plt.yscale("log")
	plt.plot(t,G_E_analytic)
	plt.show()

if __name__ == '__main__':
	main()		