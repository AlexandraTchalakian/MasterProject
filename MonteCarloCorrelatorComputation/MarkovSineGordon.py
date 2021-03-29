#coding:utf-8
from math import *
import random
import numpy as np
import sys
import os.path
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

m=1
epsilon=0.1
omega=1
A=0.5







def update(x,N,a):
	for i in range(0,N):
		x_before=x[i]
		S_before=S(i,x,N,a)
		x[i]=x[i]+random.uniform(-epsilon,epsilon)
		delta_S=S(i,x,N,a)-S_before
		#print(exp(-delta_S))
		if delta_S>0 and exp(-delta_S)<random.uniform(0,1):
			x[i]=x_before


def potential(y):
	return m*A**2*(1-cos(y))
	#return m*omega**2*y**2/2.0 

def S(j,x,N,a):
	jp=(j+1)%N
	jm=(j-1)%N
	return m*x[j]*(x[j]-x[jm]-x[jp])/a+a*potential(x[j])

def S_total(x,N,a):
	x_dot=np.gradient(x,a)
	S=0
	for i in range(0,N):
		S+=m*0.5*x_dot[i]**2*a+a*potential(x[i])
	return S



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
			update(x,N,a)
			x_axis[i]=i
			Splot[i]=S_total(x,N,a)
		plt.figure()
		plt.rc('text', usetex=True)
		plt.rc('font', family='serif')
		plt.plot(x_axis,Splot)
		plt.ylabel(r'$S_{cl}$')
		plt.xlabel(r'$N_{cor}$')
		plt.show()
		for alpha in range(0,Ncf):
			for j in range (0,Ncor):
				update(x,N,a)
			with open(nomfichier[alpha],'w') as fichier:				
				for k in range(0,len(x)):
					fichier.write(str(x[k])+"\n")
			with open(nomfichier[alpha],'r') as fichier:
				for k in range(0,len(x)):
					y[alpha][k]=float(fichier.readline())
		


	return y

"""
for j in range (0,N):
	print(S(j,x,N,a))
print("-------------------------------")
"""