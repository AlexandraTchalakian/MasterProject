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
	G1=np.zeros((Ncf,N))
	G2=np.zeros((Ncf,N))
	G3=np.zeros((Ncf,N))
	for alpha in range(0,Ncf):
		y=x[alpha]
		distances=np.zeros(len(y))
		for j in range(0,len(y)):
			distances[j]=MC.modulo(y[(j+1)%N]-y[j],0.5)
		Qst=np.cumsum(distances) # Qst est le x pas modulo 2pi
		for t in range(0,len(y)):
			temp1=0
			temp2=0
			temp3=0
			for j in range(0,len(y)):
				temp1+=Qst[(j+t)%N]*Qst[j]
				temp2+=Qst[(j+t)%N]
				temp3+=Qst[j]
			G1[alpha][t]=temp1/N
			G2[alpha][t]=temp2/N
			G3[alpha][t]=temp3/N
	(G1_mean,G1_err)=AverageErrorG(G1,N,Ncf,N_bootstrap)
	(G2_mean,G2_err)=AverageErrorG(G2,N,Ncf,N_bootstrap)
	(G3_mean,G3_err)=AverageErrorG(G3,N,Ncf,N_bootstrap)
	G_mean=np.zeros(N)
	G_err=np.zeros(N)
	for j in range(0,N):
		G_mean[j]=G1_mean[j]-G2_mean[j]*G3_mean[j]
		G_err[j]=G1_err[j]+(G2_err[j]/G2_mean[j]+G3_err[j]/G3_err[j])*G2_mean[j]*G3_mean[j]
	return (G_mean,G_err)



