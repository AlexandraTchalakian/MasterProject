#coding:utf-8
import math
import random
import numpy as np
import matplotlib.pyplot as plt

N=20
a=0.5
m=1
epsilon=1.4
Ncor=20
Ncf=100

def initX():
	x=[]
	for i in range(0,N):
		x.append(0)

	return x



def update(x):
	for i in range(0,N):
		x_before=x[i]
		S_before=S(i,x)
		x[i]=x[i]+random.uniform(-epsilon,epsilon)
		delta_S=S(i,x)-S_before
		if delta_S>0 and math.exp(-delta_S)<random.uniform(0,1):
			x[i]=x_before


def potential(y,type):
	if type=="HO":
		return y**2*0.5
	elif type=="AHO":
		return y**4*0.5
	else:
		return 0

def S(j,x):
	jp=(j+1)%N
	jm=(j-1)%N
	return m*x[j]*(x[j]-x[jm]-x[jp])/a+a*potential(x[j],"HO")




def computeG(x,n):
	g=0
	for j in range(0,N):
		g+=x[(j+n)%N]*x[j]

	return g/N

def computeG3(x,n):
	g=0
	for j in range(0,N):
		g+=x[(j+n)%N]**3*x[j]**3

	return g/N


def bootstrap(G):
	G_bootstrap=np.zeros((len(G),N))
	for i in range(0,len(G)):
		for j in range(0,N):
			alpha = int(random.uniform(0,len(G)))
			G_bootstrap[i][j]=G[alpha][j]

	return G_bootstrap


def average_var_G(G):
	average_G=np.zeros(N)
	var_G=np.zeros(N)
	for n in range (0,N):
		for alpha in range(0,Ncf):
			average_G[n]+=G[alpha][n]
			var_G[n]+=G[alpha][n]**2	

		average_G[n]/=Ncf
		var_G[n]/=Ncf
		var_G[n]-=average_G[n]**2
		var_G[n]/=Ncf

	return (average_G, var_G)


def computeGamma(x):
	G = np.zeros((Ncf,N))

	for i in range(0,5*Ncor):
		update(x)
	for i in range(0,Ncf):
		for j in range (0,Ncor):
			update(x)
		for n in range(0,N):
			G[i][n]=computeG(x,n)

	G_bootstrap1=bootstrap(G)
	(average_G_bootstrap1, var_G_bootstrap1)=average_var_G(G_bootstrap1)
	G_bootstrap2=bootstrap(G)
	(average_G_bootstrap2, var_G_bootstrap2)=average_var_G(G_bootstrap2)
	G_bootstrap3=bootstrap(G)
	(average_G_bootstrap3, var_G_bootstrap3)=average_var_G(G_bootstrap3)
	G_bootstrap4=bootstrap(G)
	(average_G_bootstrap4, var_G_bootstrap4)=average_var_G(G_bootstrap4)

	(average_G, var_G) = average_var_G(G)
		

	return (average_G, var_G, average_G_bootstrap1, var_G_bootstrap1, average_G_bootstrap2, var_G_bootstrap2, average_G_bootstrap3, var_G_bootstrap3, average_G_bootstrap4, var_G_bootstrap4)	


def Delta_E(average_G, var_G):

	error=[]
	for i in range(0,N):
		error.append(math.sqrt(var_G[i]))
	delta_E=[]
	error_E=[]
	for i in range(0,N):
		error.append(math.sqrt(var_G[i]))
		if i==N-1:
			delta_E.append(math.log(math.fabs(average_G[N-1]/average_G[0]))/a)
			error_E.append(math.fabs(error[N-1]/average_G[N-1])+math.fabs(-error[0]/average_G[0]))
		else:
			delta_E.append(math.log(math.fabs(average_G[i]/average_G[i+1]))/a)
			error_E.append(math.fabs(error[i]/average_G[i])+math.fabs(-error[i+1]/average_G[i+1]))

	return (delta_E, error_E)


def time(N):
	t=[]
	t_0=0
	for i in range(0,N):
		t.append(t_0+a*i)
	return t


x=initX()
(average_G, var_G, average_G_bootstrap1, var_G_bootstrap1,average_G_bootstrap2, var_G_bootstrap2,average_G_bootstrap3, var_G_bootstrap3,average_G_bootstrap4, var_G_bootstrap4)=computeGamma(x)
(delta_E, error_E)=Delta_E(average_G,var_G)
(delta_E_B1, error_E_B1)=Delta_E(average_G_bootstrap1,var_G_bootstrap1)
(delta_E_B2, error_E_B2)=Delta_E(average_G_bootstrap2,var_G_bootstrap2)
(delta_E_B3, error_E_B3)=Delta_E(average_G_bootstrap3,var_G_bootstrap3)
(delta_E_B4, error_E_B4)=Delta_E(average_G_bootstrap4,var_G_bootstrap4)


t_E=time(N)



plt.figure(1)
plt.errorbar(t_E, delta_E, error_E, marker='o',linestyle='none',color='black',label="original")
plt.errorbar(t_E, delta_E_B1, error_E_B1, marker='.',linestyle='none',label="bootstrap1")
plt.errorbar(t_E, delta_E_B2, error_E_B2, marker='.',linestyle='none',label="bootstrap2")
plt.errorbar(t_E, delta_E_B3, error_E_B3, marker='.',linestyle='none',label="bootstrap3")
plt.errorbar(t_E, delta_E_B4, error_E_B4, marker='.',linestyle='none',label="bootstrap4")
#plt.axhline(y=1.5,color='black',linestyle="--",linewidth=0.3)
plt.axhline(y=1,color='black',linestyle="--",linewidth=0.3)
#plt.axhline(y=-1,color='black',linestyle="--",linewidth=0.3)
plt.xlim(-0.5,3.4)
plt.ylim(-3,3)
plt.xlabel("t")
plt.ylabel('$\Delta$E(t)')
plt.legend()
plt.title("Bootstrap test for HO")
plt.show()

