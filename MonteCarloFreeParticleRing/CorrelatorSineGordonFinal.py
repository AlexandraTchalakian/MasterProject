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
	Gfrozen1=[]
	Gfrozen2=[]
	G1=np.zeros((Ncf,N))
	G2=np.zeros((Ncf,N))
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
			G1[alpha][t]=Qst[t]**2
			G2[alpha][t]=Qst[t]
		if(np.abs(Qst[-1]) < 1e-5):
		#if(np.abs(Qst[-1]) < 2.1 and np.abs(Qst[-1])>1.9):
		#if True:
			Gfrozen1.append(G1[alpha])
			Gfrozen2.append(G2[alpha])
	(G1_mean,G1_err)=AverageErrorG(Gfrozen1,N,len(Gfrozen1),N_bootstrap)
	(G2_mean,G2_err)=AverageErrorG(Gfrozen2,N,len(Gfrozen2),N_bootstrap)
	G_mean=np.zeros(N)
	G_err=np.zeros(N)

	for i in range (0,N):
		G_mean[i]=G1_mean[i]-G2_mean[i]**2
		G_err[i]=G1_err[i]+2*G2_err[i]*G2_mean[i]
	return(G_mean,G_err)

	#return (AverageG(Gfrozen,N,len(Gfrozen)),G_mean,G_err)
	#(G_mean,G_err)=AverageErrorG(G,N,Ncf,N_bootstrap)
	#return (AverageG(G,N,Ncf),G_mean,G_err)



def computeG2(x,n,N):
	g1=0
	g2=0
	g3=0
	for j in range(0,N):
		g1+=x[(j+n)%N]*x[j]
		g2+=x[(j+n)%N]
		g3+=x[j]
	return (g1/N,g2/N,g3/N)



def TwoPtFunctionVelocities(x,a,N,Ncf,N_bootstrap):
	x_dot=np.zeros((Ncf,N))
	G_dot=np.zeros((Ncf,N))
	G1_dot=np.zeros((Ncf,N))
	G2_dot=np.zeros((Ncf,N))
	G3_dot=np.zeros((Ncf,N))
	for alpha in range (0,Ncf):
		for i in range(0,N):
			if i==N-1:
				x_dot[alpha][N-1]=MC.modulo((x[alpha][0]-x[alpha][N-1]),0.5)/a
			else:
				x_dot[alpha][i]=MC.modulo((x[alpha][i+1]-x[alpha][i]),0.5)/a
		for n in range(0,N):
			(G1_dot[alpha][n],G2_dot[alpha][n],G3_dot[alpha][n])=computeG2(x_dot[alpha],n,N)

	(bootstrap_average1_dot,error_G1_dot)=AverageErrorG(G1_dot,N,Ncf,N_bootstrap)
	(bootstrap_average2_dot,error_G2_dot)=AverageErrorG(G2_dot,N,Ncf,N_bootstrap)
	(bootstrap_average3_dot,error_G3_dot)=AverageErrorG(G3_dot,N,Ncf,N_bootstrap)
	average=np.zeros(N)
	error=np.zeros(N)
	for j in range(0,N):
		average[j]=bootstrap_average1_dot[j]-bootstrap_average3_dot[j]*bootstrap_average2_dot[j]
		error[j]=error_G1_dot[j]+(error_G2_dot[j]/bootstrap_average2_dot[j]+error_G3_dot[j]/bootstrap_average3_dot[j])*bootstrap_average2_dot[j]*bootstrap_average3_dot[j]
	return (average,error)	
