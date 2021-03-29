#coding:utf-8
from math import *
import random
import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib import cm
import MarkovSineGordonFinal as MC


N_str=sys.argv[1]
N=int(N_str)
a_str=sys.argv[2]
a=float(a_str)
Ncor=7000
Ncf=2000
A=0.25


def ComputeG(x,n):
	g1=0
	g2=0
	g3=0
	for j in range(0,N):
		g1+=x[(j+n)%N]*x[j]
		g2+=x[j]
		g3+=x[(j+n)%N]

	return (g1/N,g2/N,g3/N)


def MonteCarloAverage(x):
	G1=np.zeros((Ncf,N))
	G2=np.zeros((Ncf,N))
	G3=np.zeros((Ncf,N))
	for alpha in range(0,Ncf):
		for n in range(0,N):
			(G1[alpha][n],G2[alpha][n],G3[alpha][n])=ComputeG(x[alpha],n)
	return (G1,G2,G3)


def AverageG(G1,G2,G3,type):
	averageG1=np.zeros(N)
	averageG2=np.zeros(N)
	averageG3=np.zeros(N)
	averageG=np.zeros(N)
	for n in range (0,N):
		for alpha in range (0,Ncf):
			averageG1[n]+=G1[alpha][n]/Ncf
			averageG2[n]+=G2[alpha][n]/Ncf
			averageG3[n]+=G3[alpha][n]/Ncf
	for n in range(0,N):
		if type!="zero":
			averageG[n]=averageG1[n]-averageG2[n]*averageG3[n]
		else:
			averageG[n]=averageG1[n]
	return averageG


def Bootstrap(G):
	G_bootstrap=np.zeros((Ncf,N))
	for i in range(0,Ncf):
		for j in range(0,N):
			alpha = int(random.uniform(0,Ncf))
			G_bootstrap[i][j]=G[alpha][j]
	return G_bootstrap



def BootstrapAverage(G,N_bootstrap):
	avg=np.zeros((N_bootstrap,N))
	G_bootstrap=np.zeros((N_bootstrap,Ncf,N))
	for nb in range (0,N_bootstrap):
		G_bootstrap[nb]=Bootstrap(G)

	for nb in range (0,N_bootstrap):
		for j in range (0,N):
			for alpha in range (0,Ncf):
				avg[nb][j]+=G_bootstrap[nb][alpha][j]/Ncf
	return avg


def AverageErrorBootstrapG(G,N_bootstrap):
	avg=BootstrapAverage(G,N_bootstrap)
	average_G=np.zeros(N)
	var_G=np.zeros(N)
	for j in range (0,N):
		for nb in range (0,N_bootstrap):
			average_G[j]+=avg[nb][j]/N_bootstrap
			var_G[j]+=avg[nb][j]**2/(N_bootstrap-1)
		var_G[j]-=average_G[j]**2
		var_G[j]=sqrt(var_G[j])
	return(average_G,var_G)




def AverageErrorTwoPtFunction(G1,G2,G3,N_bootstrap):
	(average_G1,var_G1)=AverageErrorBootstrapG(G1,N_bootstrap)
	(average_G2,var_G2)=AverageErrorBootstrapG(G2,N_bootstrap)
	(average_G3,var_G3)=AverageErrorBootstrapG(G3,N_bootstrap)
	var_G23=np.zeros(N)
	cov_G=np.zeros(N)
	var_G=np.zeros(N)
	average_G=np.zeros(N)
	average_G23=np.zeros(N)
	error_G=np.zeros(N)
	for j in range (0,N):
		average_G[j]=average_G1[j]-average_G2[j]*average_G3[j]
		error_G[j]=sqrt(var_G1[j])-sqrt(var_G2[j]*var_G3[j])

	return (average_G,error_G)
	"""
	avg1=BootstrapAverage(G1,N_bootstrap)
	avg2=BootstrapAverage(G2,N_bootstrap)
	avg3=BootstrapAverage(G3,N_bootstrap)
	for j in range(0,N):
		var_G23[j]=var_G2[j]*var_G3[j]+var_G2[j]*average_G3[j]**2+var_G3[j]*average_G2[j]**2
		for nb in range(0,N_bootstrap):
			average_G23[j]+=avg2[nb][j]*avg3[nb][j]/(N_bootstrap)
			cov_G[j]+=(avg1[nb][j]-average_G1[j])*(avg2[nb][j]*avg3[nb][j]-average_G23[j])/(N_bootstrap-1)
		var_G[j]=sqrt(abs(var_G1[j]+var_G23[j]-2*cov_G[j]))
	return(average_G1-average_G23,var_G)
	"""

def time():
	t0=0
	t=np.zeros(N)
	for i in range(0,N):
		t[i]=t0+i*a
	return t



def main():
	t=time()
	x=MC.markov_steps(N,Ncf,Ncor,a)
	(G1,G2,G3)=MonteCarloAverage(x)
	two_pt_function=AverageG(G1,G2,G3,"zero")
	(two_pt_function_bootstrap,two_pt_function_error)=AverageErrorBootstrapG(G1,100)
	two_pt_function1=AverageG(G1,G2,G3,"no_zero")
	(two_pt_function1_bootstrap,two_pt_function1_error)=AverageErrorTwoPtFunction(G1,G2,G3,100)


	plt.figure(1)
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	#plt.errorbar(t,two_pt_function,two_pt_function_error,marker='.',linestyle='none')
	#plt.errorbar(t,two_pt_function_bootstrap,two_pt_function_error,marker='.',linestyle='none')
	plt.plot(t,two_pt_function,marker='.',linestyle='none')
	plt.ylabel(r'$\langle x(t+\Delta t)x(t)\rangle$')
	plt.xlabel(r'$t$')

	plt.figure(2)
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	#plt.errorbar(t,two_pt_function1,two_pt_function1_error,marker='.',linestyle='none')
	#plt.errorbar(t,two_pt_function1_bootstrap,two_pt_function1_error,marker='.',linestyle='none')
	plt.plot(t,two_pt_function1,marker='.',linestyle='none')
	plt.ylabel(r'$\langle x(t+\Delta t)x(t)\rangle$')
	plt.xlabel(r'$t$')
	plt.show()


if __name__ == '__main__':
	main()	