#coding:utf-8
from math import *
import random
import numpy as np
import os.path
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter


m=1
epsilon=0.5
A=0.25


def Potential(y):
	return m*A**2*(1-cos(y))


def Local_S(j,x,N,a):
	return m*x[j]*(x[j]-x[(j+1)%N]-x[(j-1)%N])/a+a*Potential(x[j])



def Total_S(x,N,a):
	S_tot=0
	x_dot=np.gradient(x)
	for j in range (0,N):
		S_tot+=a*(0.5*m*x_dot[j]**2+Potential(x[j]))

	return S_tot



def Update(x,N_hits,N,a):
	count_rejected=0
	total_count=0
	for i in range(0,N):
		for n_hits in range (0,N_hits):
			x_before=x[i]
			S_before=Local_S(i,x,N,a)
			x[i]=x[i]+random.uniform(-epsilon,epsilon)
			delta_S=Local_S(i,x,N,a)-S_before
			total_count+=1
			if delta_S>0 and exp(-delta_S)<random.uniform(0,1):
				x[i]=x_before
				count_rejected+=1
	#print(count_rejected,total_count)


def markov_steps(N,Ncf,Ncor,a):

	x=[0.0]*N
	exist=True
	nomfichier=[]
	N_str=str(N)
	Ncor_str=str(Ncor)
	a_str=str(a)
	y=np.zeros((Ncf,N))

	for alpha in range(0,Ncf):
		Ncf_string=str(alpha)
		nomfichier.append('SineGordonMarkov_'+Ncor_str+'_'+Ncf_string+'_'+N_str+'_'+a_str+'.txt')
		if not os.path.isfile(nomfichier[alpha]):
			exist=False


	if exist:
		for alpha in range(0,Ncf):
			with open(nomfichier[alpha],'r') as fichier:
				for k in range(0,len(x)):
					y[alpha][k]=float(fichier.readline())

					
	if not exist:
		Splot=np.zeros(5*Ncor)
		x_axis=np.zeros(5*Ncor)		
		for i in range(0,5*Ncor):
			Update(x,20,N,a)
			x_axis[i]=i
			Splot[i]=Total_S(x,N,a)
		plt.figure()
		plt.rc('text', usetex=True)
		plt.rc('font', family='serif')
		plt.plot(x_axis,Splot)
		plt.ylabel(r'$S_{cl}$')
		plt.xlabel(r'$N_{cor}$')
		plt.show()
		for alpha in range(0,Ncf):
			for j in range (0,Ncor):
				Update(x,20,N,a)
			with open(nomfichier[alpha],'w') as fichier:				
				for k in range(0,len(x)):
					fichier.write(str(x[k])+"\n")
			with open(nomfichier[alpha],'r') as fichier:
				for k in range(0,len(x)):
					y[alpha][k]=float(fichier.readline())
		


	return y
