from math import *
import numpy as np
import Markov_chain as MC
import random 
import CompQ




def AverageG(G,N,Ncf):
	average_G=np.zeros(N)
	for n in range(0,N):
		for alpha in range (0,Ncf):
			average_G[n]+=G[alpha][n]/Ncf

	return average_G




def bootstrapG(G,N):
	G_bootstrap=np.zeros((len(G),N))
	for i in range(0,len(G)):
		for j in range(0,N):
			alpha = int(random.uniform(0,len(G)))
			G_bootstrap[i][j]=G[alpha][j]

	return G_bootstrap

def BootstrapAverage(G,N,Ncf,N_bootstrap):
	avg=np.zeros((N_bootstrap,N))
	G_bootstrap=np.zeros((N_bootstrap,Ncf,N))
	for nb in range (0,N_bootstrap):
		G_bootstrap[nb]=bootstrapG(G,N)

	for nb in range (0,N_bootstrap):
		for j in range (0,N):
			for alpha in range (0,Ncf):
				avg[nb][j]+=G_bootstrap[nb][alpha][j]/Ncf
	return avg
	


def AverageErrorG(G,N,Ncf,N_bootstrap):
	avg=BootstrapAverage(G,N,Ncf,N_bootstrap)
	average_G=np.zeros(N)
	var_G=np.zeros(N)
	for j in range (0,N):
		for nb in range (0,N_bootstrap):
			average_G[j]+=avg[nb][j]/N_bootstrap
			var_G[j]+=avg[nb][j]**2/N_bootstrap
		var_G[j]-=average_G[j]**2
		var_G[j]=sqrt(var_G[j])
	return(average_G,var_G)




def TwoPtFunction(x,a,N,Ncf,N_bootstrap):
	Gfrozen=[]
	G=np.zeros((Ncf,N))
	for alpha in range(0,Ncf):
		y=x[alpha]
		Qst=np.zeros(N)
		distances=np.zeros(len(y))
		for j in range(0,len(y)):
			distances[j]=MC.modulo(y[(j+1)%N]-y[j],0.5)
		for t in range(0,N):
			Qst[t]=np.sum(distances[0:t+1]) # Qst est le x pas modulo 2pi
		for t in range(0,len(y)):
			temp=0
			#for j in range(0,len(y)):
				#temp+=(Qst[(j+t)%N])*(Qst[j])
			G[alpha][t]=Qst[t]**2
		#if(np.abs(Qst[-1]) < 1e-5):
		if True:
			Gfrozen.append(G[alpha])
	(G_mean,G_err)=AverageErrorG(Gfrozen,N,len(Gfrozen),N_bootstrap)
	return (AverageG(Gfrozen,N,len(Gfrozen)),G_mean,G_err)
	#(G_mean,G_err)=AverageErrorG(G,N,Ncf,N_bootstrap)
	#return (AverageG(G,N,Ncf),G_mean,G_err)





"""
def TwoPtFunction(x,a,N,Ncf,N_bootstrap):
	G=np.zeros((Ncf,N))
	for alpha in range(0,Ncf):
		y=x[alpha]
		for t in range(0,len(y)):
			distances1=0
			distances2=0
			for j in range(0,len(y)):
				distances1+=MC.modulo(y[(j+t+1)%N]-y[(j+t)%N],0.5)
				distances2+=MC.modulo(y[(j+1)%N]-y[(j)%N],0.5)
			G[alpha][t]=distances1*distances2/N
	(G_mean,G_err)=AverageErrorG(G,N,Ncf,N_bootstrap)
	return (G_mean,G_err)
"""