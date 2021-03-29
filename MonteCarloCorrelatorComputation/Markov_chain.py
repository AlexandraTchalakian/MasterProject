from math import *
import random 
import numpy as np
import os.path
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter



def modulo(y,n):
	if abs(y)<=n:
		return y
	if y > n:
		return y-1
	if y <-n:
		return y+1


def potential(y,Type,Omega):
	if(Type=="paper"):
		return Omega**2*cos(2*pi*y)
	if(Type=="HO"):
		return 0.5*y**2
	if (Type=="SineGordon"):
		return Omega**2*(1-cos(2*pi*y))
	else:
		return 0


def Total_S(x,N,a,m,Type,Omega):
	S_tot=0
	for j in range (0,N):
		S_tot+=0.5*m*(modulo(x[(j+1)%N]-x[j],0.5))**2/a+a*potential(x[j],Type,Omega)

	return S_tot


def S(j,x,a,N,m,Omega,Type):
	jp=(j+1)%N
	jm=(j-1)%N
	return 0.5*m/a*(modulo((x[jp]-x[j]),0.5)**2+modulo((x[j]-x[jm]),0.5)**2)+a*potential(x[j],Type,Omega)





def update(x,a,N,Nhits, Delta,m,Omega,Type):
	count_rejected=0
	total_count=0
	for i in range(0,N):
		for nhits in range(0,Nhits):
			x_before=x[i]
			S_before=S(i,x,a,N,m,Omega,Type)
			x[i]=(x[i]+(1-2*random.uniform(0,1))*Delta)%1
			delta_S=S(i,x,a,N,m,Omega,Type)-S_before
			total_count+=1
			if delta_S>0 and exp(-delta_S)<random.uniform(0,1):
				x[i]=x_before
				count_rejected+=1
	#print(count_rejected,total_count)






def markov_steps(a,N,Ncf,Ncor,Nhits,Delta,m,Omega,Type):

	N_string=str(N)
	a_string=str(a)
	Ncor_string=str(Ncor)
	Omega_string=str(Omega)
	x=[0.0]*N
	exist=True
	nomfichier=[]
	y=np.zeros((Ncf,N))

	for alpha in range(0,Ncf):
		Ncf_string=str(alpha)
		nomfichier.append('Markov_'+Ncf_string+"_"+Ncor_string+"_"+N_string+"_"+Omega_string+"_"+a_string+"_"+Type+"_"+'a0.txt')
		if not os.path.isfile(nomfichier[alpha]):
			exist=False


	if exist:
		for alpha in range(0,Ncf):
			with open(nomfichier[alpha],'r') as fichier:
				for k in range(0,len(x)):
					y[alpha][k]=float(fichier.readline())

						
	if not exist:	
		for i in range(0,5*Ncor):
			update(x,a,N,Nhits,Delta,m,Omega,Type)
		for alpha in range(0,Ncf):
			for j in range (0,Ncor):
				update(x,a,N,Nhits,Delta,m,Omega,Type)
			with open(nomfichier[alpha],'w') as fichier:				
				for k in range(0,len(x)):
					fichier.write(str(x[k])+"\n")
			with open(nomfichier[alpha],'r') as fichier:
				for k in range(0,len(x)):
					y[alpha][k]=float(fichier.readline())	



	return y


	""""

		Splot=np.zeros(5*Ncor)
		x_axis=np.zeros(5*Ncor)	

			x_axis[i]=i
			Splot[i]=Total_S(x,N,a,m,Type,Omega)
		plt.figure()
		plt.rc('text', usetex=True)
		plt.rc('font', family='serif')Ncf_string
		plt.plot(x_axis,Splot)
		plt.ylabel(r'$S_{cl}$')
		plt.xlabel(r'$N_{cor}$')
		plt.show()

	"""