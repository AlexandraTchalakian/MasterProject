#coding:utf-8
from math import *
import random
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter


N=20
a=0.05
m=1.
epsilon=0.2
Ncor=5000
Ncf=1000
omega=1.


def Potential(y):
	return 0.5*m*omega**2*y**2

def Local_S(j,x):
	return m*x[j]*(x[j]-x[(j+1)%N]-x[(j-1)%N])/a+a*Potential(x[j])

def Total_S(x):
	S_tot=0
	x_dot=np.gradient(x)
	for j in range (0,N):
		S_tot+=a*(0.5*m*x_dot[j]**2+Potential(x[j]))

	return S_tot


def Update(x,N_hits):
	count_rejected=0
	total_count=0
	for i in range(0,N):
		for n_hits in range (0,N_hits):
			x_before=x[i]
			S_before=Local_S(i,x)
			x[i]=x[i]+random.uniform(-epsilon,epsilon)
			delta_S=Local_S(i,x)-S_before
			#print(exp(-delta_S),delta_S)
			total_count+=1
			if delta_S>0 and exp(-delta_S)<random.uniform(0,1):
				x[i]=x_before
				count_rejected+=1
	#print(count_rejected,total_count)



def ComputeG(x,n):
	Gn1=0
	Gn2=0
	Gn3=0
	for j in range (0,N):
		Gn1+=x[(j+n)%N]*x[j]/N
		Gn2+=x[(j+n)%N]/N
		Gn2+=x[j]/N
	return (Gn1,Gn2,Gn3)



def MonteCarloAverage():
	x=np.zeros(N)
	G1=np.zeros((Ncf,N))
	G2=np.zeros((Ncf,N))
	G3=np.zeros((Ncf,N))
	for k in range (0,5*Ncor):
		Update(x,15)
	for alpha in range(0,Ncf):
		for k in range(0,Ncor):
			Update(x,15)
		for n in range(0,N):
			(G1[alpha][n],G2[alpha][n],G3[alpha][n])=ComputeG(x,n)
	return (G1,G2,G3)



def AverageG(G1,G2,G3):
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
		averageG[n]=averageG1[n]-averageG2[n]*averageG3[n]
	return averageG



def Bootstrap(G):
	G_bootstrap=np.zeros((Ncf,N))
	for i in range(0,Ncf):
		for j in range(0,N):
			alpha = int(random.uniform(0,Ncf))
			G_bootstrap[i][j]=G[alpha][j]
	return G_bootstrap


def AverageErrorBootstrapG1(G,N_bootstrap):
	avg=np.zeros((N_bootstrap,N))
	G_bootstrap=np.zeros((N_bootstrap,Ncf,N))
	for nb in range (0,N_bootstrap):
		G_bootstrap[nb]=Bootstrap(G)

	for nb in range (0,N_bootstrap):
		for j in range (0,N):
			for alpha in range (0,Ncf):
				avg[nb][j]+=G_bootstrap[nb][alpha][j]/Ncf

	average_G=np.zeros(N)
	var_G=np.zeros(N)
	for j in range (0,N):
		for nb in range (0,N_bootstrap):
			average_G[j]+=avg[nb][j]/N_bootstrap
			var_G[j]+=avg[nb][j]**2/N_bootstrap
		var_G[j]-=average_G[j]**2
		var_G[j]=sqrt(var_G[j])
	return(average_G,var_G)


def AverageErrorBootstrapG123(G1,G2,G3,N_bootstrap):
	avg1=np.zeros((N_bootstrap,N))
	avg2=np.zeros((N_bootstrap,N))
	avg3=np.zeros((N_bootstrap,N))
	G_bootstrap1=np.zeros((N_bootstrap,Ncf,N))
	G_bootstrap2=np.zeros((N_bootstrap,Ncf,N))
	G_bootstrap3=np.zeros((N_bootstrap,Ncf,N))
	for nb in range (0,N_bootstrap):
		G_bootstrap1[nb]=Bootstrap(G1)
		G_bootstrap2[nb]=Bootstrap(G2)
		G_bootstrap3[nb]=Bootstrap(G3)

	for nb in range (0,N_bootstrap):
		for j in range (0,N):
			for alpha in range (0,Ncf):
				avg1[nb][j]+=G_bootstrap1[nb][alpha][j]/Ncf
				avg2[nb][j]+=G_bootstrap2[nb][alpha][j]/Ncf
				avg3[nb][j]+=G_bootstrap3[nb][alpha][j]/Ncf

	average_G=np.zeros(N)
	average_G23=np.zeros(N)
	average_G1=np.zeros(N)
	average_G2=np.zeros(N)
	average_G3=np.zeros(N)
	var_G23=np.zeros(N)
	var_G1=np.zeros(N)
	var_G2=np.zeros(N)
	var_G3=np.zeros(N)
	for j in range (0,N):
		for nb in range (0,N_bootstrap):
			average_G[j]+=(avg1[nb][j]+avg2[nb][j]*avg3[nb][j])/N_bootstrap
			average_G23[j]+=(avg2[nb][j]*avg3[nb][j])/N_bootstrap
			average_G1[j]+=avg2[nb][j]/N_bootstrap
			average_G2[j]+=avg2[nb][j]/N_bootstrap
			average_G3[j]+=avg3[nb][j]/N_bootstrap
			var_G1[j]+=avg1[nb][j]**2/N_bootstrap
			var_G2[j]+=avg2[nb][j]**2/N_bootstrap
			var_G3[j]+=avg3[nb][j]**2/N_bootstrap
		var_G1[j]-=average_G1[j]**2
		var_G2[j]-=average_G2[j]**2
		var_G3[j]-=average_G3[j]**2
		var_G23[j]=var_G2[j]*var_G3[j]+var_G2[j]*average_G3[j]**2+var_G3[j]*average_G2[j]**2
	cov_G=np.zeros(N)
	error_G=np.zeros(N)
	for j in range(0,N):
		for nb in range(0,N_bootstrap):
			cov_G[j]+=(avg1[nb][j]-average_G1[j])*(avg2[nb][j]*avg3[nb][j]-average_G23[j])/N_bootstrap
		error_G[j]=sqrt(abs(var_G1[j]+var_G23[j]-2*cov_G[j]))

	return(average_G,error_G)




def time():
	t0=0
	t=np.zeros(N)
	for i in range(0,N):
		t[i]=t0+i*a
	return t

def analyticG(t):
	analytic_G=np.zeros(len(t))
	analytic_G_dot=np.zeros(len(t))
	for i in range (0, len(t)):
		analytic_G[i]=1/(2*omega*m*sinh(omega*N*a/2.))*cosh(omega*(N*a*0.5-t[i]))
		analytic_G_dot[i]=omega/(2*m*sinh(omega*N*a/2.))*cosh(omega*(N*a*0.5-t[i]))
	return (analytic_G, analytic_G_dot)	


def main():
	t=time()
	
	two_pt_function=AverageG(G1,G2,G3)
	(two_pt_function_bootstrap1,two_pt_function_error1)=AverageErrorBootstrapG1(G1,100)
	#(two_pt_function_bootstrap2,two_pt_function_error2)=AverageErrorBootstrapG123(G1,G2,G3,100)

	
	(one_pt_function_bootstrap2,one_pt_function_error2)=AverageErrorBootstrapG1(G2,100)
	(one_pt_function_bootstrap3,one_pt_function_error3)=AverageErrorBootstrapG1(G3,100)
	two_pt_function_error=np.zeros(N)
	two_pt_function_bootstrap=np.zeros(N)
	for n in range (0,N):
		two_pt_function_error[n]=two_pt_function_error1[n]-one_pt_function_error2[n]*one_pt_function_error3[n]
		two_pt_function_bootstrap[n]=two_pt_function_bootstrap1[n]-one_pt_function_bootstrap2[n]*one_pt_function_bootstrap3[n]

	(two_pt_function_analytic,two_pt_function_dot_analytic)=analyticG(t)

	plt.figure()
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	plt.errorbar(t,two_pt_function,two_pt_function_error,marker='.',linestyle='none')
	plt.errorbar(t,two_pt_function_bootstrap1,two_pt_function_error1,marker='.',linestyle='none')
	plt.plot(t,two_pt_function_analytic)
	plt.ylabel(r'$\langle x(t+\Delta t)x(t)\rangle$')
	plt.xlabel(r'$t$')
	plt.show()


if __name__ == '__main__':
	main()	