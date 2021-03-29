#coding:utf-8
from math import *
import random
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib import cm

N=20
a=0.01
m=1.
epsilon=1.4
Ncor=500
Ncf=250
omega=1.




def update(x):
	for i in range(0,N):
		x_before=x[i]
		S_before=S(i,x)
		x[i]=x[i]+random.uniform(-epsilon,epsilon)
		delta_S=S(i,x)-S_before
		if delta_S>0 and exp(-delta_S)<random.uniform(0,1):
			x[i]=x_before

def update_2(x):
	for i in range(0,N):
		x_before=x[i]
		S_before=S_2(i,x)
		x[i]=x[i]+random.uniform(-epsilon,epsilon)
		delta_S=S_2(i,x)-S_before
		if delta_S>0 and exp(-delta_S)<random.uniform(0,1):
			x[i]=x_before



def potential(y,type):
	if type=="HO":
		return m*omega**2*y**2*0.5
	elif type=="AHO":
		return y**4*0.5
	else:
		return 0

def S_total(x):
	x_dot=np.gradient(x,a)
	S=0
	for i in range(0,N):
		S+=m*0.5*x_dot[i]**2*a+m*0.5*omega**2*x[i]**2*a
	return S

	

def S(j,x):
	jp=(j+1)%N
	jm=(j-1)%N
	return m*x[j]*(x[j]-x[jm]-x[jp])/a+a*potential(x[j],"HO")

def S_2(j,x):
	jp=(j+1)%N
	jm=(j-1)%N
	return (x[j]*(x[j]-x[jm]-x[jp]))**2/a**3+a*potential(x[j],"zeros")


def computeG(x,n,ndot):
	g1=0
	g2=0
	g3=0
	g4=0
	if ndot=="dot":
		for j in range(0,N):
			g1+=x[(j+n+1)%N]*x[(j+1)%N]/(4*a**2)
			g2+=x[(j+n+1)%N]*x[(j-1)%N]/(4*a**2)
			g3+=x[(j+n-1)%N]*x[(j+1)%N]/(4*a**2)
			g4+=x[(j+n-1)%N]*x[(j-1)%N]/(4*a**2)
		return (g1/N,g2/N,g3/N,g4/N)
		"""
			g1+=x_dot[(j+n)%N]*x_dot[j]
			g2+=x_dot[j]
			g3+=x_dot[(j+n)%N]

			g1+=(x[(n+j+1)%N]-x[((n+j))%N])*(x[(j+1)%N]-x[(j)%N])/(a**2)
			g2+=(x[(j+1)%N]-x[j])/(a**2)
			g3+=(x[(j+n+1)%N]-x[(j+n)%N])/(a**2)
		"""
	else:
		for j in range(0,N):
			g1+=x[(j+n)%N]*x[j]
			g2+=x[j]
			g3+=x[(j+n)%N]
		return (g1/N,g2/N,g3/N)



def computeG3(x,n):
	g=0
	for j in range(0,N):
		g+=x[(j+n)%N]**3*x[j]**3

	return g/N


def bootstrapG(G):
	G_bootstrap=np.zeros((len(G),N))
	for i in range(0,len(G)):
		for j in range(0,N):
			alpha = int(random.uniform(0,len(G)))
			G_bootstrap[i][j]=G[alpha][j]

	return G_bootstrap

""""
def averageG(G):
	average_G=np.zeros(N)
	
	for n in range (0,N):
		for alpha in range(0,Ncf):
			average_G[n]+=G[alpha][n]	

		average_G[n]/=Ncf

	return average_G
"""

def averageGdot(G1,G2,G3,G4):
	average_G1=np.zeros(N)
	average_G2=np.zeros(N)
	average_G3=np.zeros(N)
	average_G4=np.zeros(N)
	average_G=np.zeros(N)
	
	for n in range (0,N):
		for alpha in range(0,Ncf):
			average_G1[n]+=G1[alpha][n]	
			average_G2[n]+=G2[alpha][n]
			average_G3[n]+=G3[alpha][n]
			average_G4[n]+=G4[alpha][n]
		average_G[n]=(average_G1[n]-average_G2[n]-average_G3[n]+average_G4[n])/Ncf
	

	return average_G




def averageG(G1,G2,G3):
	average_G1=np.zeros(N)
	average_G2=np.zeros(N)
	average_G3=np.zeros(N)
	average_G=np.zeros(N)
	for n in range (0,N):
		for alpha in range(0,Ncf):
			average_G1[n]+=G1[alpha][n]	
			average_G2[n]+=G2[alpha][n]
			average_G3[n]+=G3[alpha][n]
		average_G[n]=average_G1[n]/Ncf-average_G2[n]*average_G3[n]/Ncf**2
			
	return average_G






def AverageErrorG(G,N_bootstrap):
	print(len(G))
	avg=np.zeros((len(G),N))
	for i in range (0,N_bootstrap):
		bootstrap=bootstrapG(G)
		for j in range (0,len(G)):
			for k in range (0,N):
				avg[j][k]+=bootstrap[j][k]/N_bootstrap

	average_G=np.zeros(N)

	var_G=np.zeros(N)
	for n in range(0,N):
		for alpha in range (0,len(G)):
			average_G[n]+=avg[alpha][n]
			var_G[n]+=avg[alpha][n]**2	

		average_G[n]/=len(G)
		var_G[n]/=len(G)
		var_G[n]-=average_G[n]**2
		var_G[n]/=len(G)

	error_G=np.zeros(N)
	for n in range(0,N):
		error_G[n]=sqrt(var_G[n])

	return (average_G,error_G)	

def computeGamma(x):
	G1 = np.zeros((Ncf,N))
	G1_dot = np.zeros((Ncf,N))
	G2 = np.zeros((Ncf,N))
	G2_dot = np.zeros((Ncf,N))
	G3 = np.zeros((Ncf,N))
	G3_dot = np.zeros((Ncf,N))
	G4_dot = np.zeros((Ncf,N))
	Splot=np.zeros(5*Ncor)
	x_axis=np.zeros(5*Ncor)		
	for i in range(0,5*Ncor):
		update(x)
		x_axis[i]=i
		Splot[i]=S_total(x)
	plt.figure()
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	plt.plot(x_axis,Splot)
	plt.ylabel(r'$S_{cl}$')
	plt.xlabel(r'$N_{cor}$')
	plt.show()
	for i in range(0,Ncf):
		for j in range (0,Ncor):
			update(x)
			#update_2(x)
		for n in range(0,N):
			(G1[i][n],G2[i][n],G3[i][n])=computeG(x,n,"ndot")
			(G1_dot[i][n],G2_dot[i][n],G3_dot[i][n],G4_dot[i][n])=computeG(x,n,"dot")
		


	return (G1, G1_dot,G2, G2_dot,G3, G3_dot,G4_dot)


def time():
	t=[]
	t_0=0
	for i in range(0,N):
		t.append(t_0+a*i)
	return t



def analyticG(t):
	analytic_G=np.zeros(len(t))
	analytic_G_dot=np.zeros(len(t))
	for i in range (0, len(t)):
		analytic_G[i]=1/(2*omega*m*sinh(omega*N*a/2.))*cosh(omega*(N*a*0.5-t[i]))
		analytic_G_dot[i]=omega/(2*m*sinh(omega*N*a/2.))*cosh(omega*(N*a*0.5-t[i]))

	return (analytic_G, analytic_G_dot)	

def main():
	x=np.zeros(N)
	(G1, G1_dot,G2, G2_dot,G3, G3_dot,G4_dot)=computeGamma(x)
	#average_G=averageG(G1)
	average_G=averageG(G1,G2,G3)
	#average_G_dot=averageG(G_dot)
	average_G_dot=averageGdot(G1_dot,G2_dot,G3_dot,G4_dot)
	abs_average_G_dot=np.zeros(N)
	for i in range (0,N):
		abs_average_G_dot[i]=abs(average_G_dot[i])
	(average_G1,error_G)=AverageErrorG(G1,100)
	#error_G_dot=errorG(G1_dot,10)
	t=time()
	(analytic_G, analytic_G_dot)=analyticG(t)



	plt.figure()
	plt.errorbar(t, average_G, error_G, marker='.',linestyle='none')
	plt.errorbar(t, average_G1, error_G, marker='.',linestyle='none')
	#plt.plot(t, average_G, marker='.',linestyle='none')
	plt.plot(t,analytic_G)
	plt.ylabel('G(t)')
	plt.xlabel('t')
	plt.legend(('analytic','bootstrap','numeric'), loc='upper center')
	plt.show()

	"""
	plt.figure(2)
	#plt.errorbar(t,abs_average_G_dot, error_G_dot, marker='.',linestyle='none')
	plt.plot(t, average_G_dot, marker='.',linestyle='none')
	plt.plot(t,analytic_G_dot)
	plt.ylabel('G(t)')
	plt.xlabel('t'),
	plt.legend(('analytic','numeric'), loc='upper center')
	plt.show()
	"""

if __name__ == '__main__':
	main()	