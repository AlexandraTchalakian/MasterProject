from math import *
import numpy as np
import Markov_chain as MC
import random 
import CompQ

def computeG3(x,n,N):
	g1=0
	g2=0
	g3=0
	g4=0
	g5=0
	for j in range(0,N):
		g1+=x[(j+n)%N]*x[j]
		g2+=x[(j+n)%N]
		g3+=x[j]
		g4+=x[j]*x[j]
		g5+=x[(j+n)%N]*x[(j+n)%N]
		#g5+=(x[(j+n)%N]-x[j])**2
	return (g1/N,g2/N,g3/N,g4/N,g5/N)

def computeG2(x,n,N):
	g1=0
	g2=0
	g3=0
	for j in range(0,N):
		g1+=x[(j+n)%N]*x[j]
		g2+=x[(j+n)%N]
		g3+=x[j]

	return (g1/N,g2/N,g3/N)

def computeG(x,n,N):
	g=0
	for j in range(0,N):
		g+=x[(j+n)%N]*x[j]

	return g/N



def averageG2(G1,G2,G3,N,Ncf):
	average_G=np.zeros(N)
	average_G1=np.zeros(N)
	average_G2=np.zeros(N)
	average_G3=np.zeros(N)
	for n in range(0,N):
		for alpha in range (0,Ncf):
			average_G1[n]+=G1[alpha][n]/Ncf
			average_G2[n]+=G2[alpha][n]/Ncf
			average_G3[n]+=G3[alpha][n]/Ncf
		average_G[n]=average_G1[n]-average_G2[n]*average_G3[n]

	return average_G

def averageG3(G1,G2,G3,G4,G5,N,Ncf):
	average_G=np.zeros(N)
	average_G1=np.zeros(N)
	average_G2=np.zeros(N)
	average_G3=np.zeros(N)
	average_G4=np.zeros(N)
	average_G5=np.zeros(N)
	for n in range(0,N):
		for alpha in range (0,Ncf):
			average_G1[n]+=G1[alpha][n]/Ncf
			average_G2[n]+=G2[alpha][n]/Ncf
			average_G3[n]+=G3[alpha][n]/Ncf
			average_G4[n]+=G4[alpha][n]/Ncf
			average_G5[n]+=G5[alpha][n]/Ncf
		average_G[n]=(average_G5[n]-average_G2[n]**2)-2*(average_G1[n]-average_G2[n]*average_G3[n])+(average_G4[n]-average_G3[n]**2)
		#average_G[n]=average_G5[n]
	return average_G



def averageG(G,N,Ncf):
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



def TwoPtFunction(x,a,N,Ncf,N_bootstrap,ndot,terms):
	x_dot=np.zeros((Ncf,N))
	G=np.zeros((Ncf,N))
	G1=np.zeros((Ncf,N))
	G2=np.zeros((Ncf,N))
	G3=np.zeros((Ncf,N))
	G4=np.zeros((Ncf,N))
	G5=np.zeros((Ncf,N))
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
			if terms==1:
				G[alpha][n]=computeG(x[alpha],n,N)
				G_dot[alpha][n]=computeG(x_dot[alpha],n,N)
			if terms==2:
				(G1[alpha][n],G2[alpha][n],G3[alpha][n])=computeG2(x[alpha],n,N)
				(G1_dot[alpha][n],G2_dot[alpha][n],G3_dot[alpha][n])=computeG2(x_dot[alpha],n,N)
			if terms==3:
				(G1[alpha][n],G2[alpha][n],G3[alpha][n],G4[alpha][n],G5[alpha][n])=computeG3(x[alpha],n,N)


	if (ndot=="dot" and terms==1):
		average_G_dot = averageG(G_dot,N,Ncf)
		(bootstrap_average_dot,error_G_dot)=AverageErrorG(G_dot,N,Ncf,N_bootstrap)
		return (average_G_dot,bootstrap_average_dot,error_G_dot,G_dot)
	elif(ndot=="ndot" and terms==1):
		average_G = averageG(G,N,Ncf)
		(bootstrap_average,error_G)=AverageErrorG(G,N,Ncf,N_bootstrap)
		return (average_G,bootstrap_average,error_G,G)
	elif(ndot=="ndot" and terms==2):
		average_G=averageG2(G1,G2,G3,N,Ncf)
		(bootstrap_average1,error_G1)=AverageErrorG(G1,N,Ncf,N_bootstrap)
		(bootstrap_average2,error_G2)=AverageErrorG(G2,N,Ncf,N_bootstrap)
		(bootstrap_average3,error_G3)=AverageErrorG(G3,N,Ncf,N_bootstrap)
		error_G=np.zeros(N)
		bootstrap_average=np.zeros(N)
		for j in range(0,N):
			bootstrap_average[j]=bootstrap_average1[j]-bootstrap_average2[j]*bootstrap_average3[j]
			error_G[j]=error_G1[j]+(error_G2[j]/bootstrap_average2[j]+error_G3[j]/bootstrap_average3[j])*bootstrap_average2[j]*bootstrap_average3[j]
		return (average_G,bootstrap_average,error_G)
	elif(ndot=="dot" and terms==2):
		average_G_dot=averageG2(G1_dot,G2_dot,G3_dot,N,Ncf)
		return average_G_dot
	elif(ndot=="ndot" and terms==3):
		average_G=averageG3(G1,G2,G3,G4,G5,N,Ncf)
		(bootstrap_average1,error_G1)=AverageErrorG(G1,N,Ncf,N_bootstrap)
		(bootstrap_average2,error_G2)=AverageErrorG(G2,N,Ncf,N_bootstrap)
		(bootstrap_average3,error_G3)=AverageErrorG(G3,N,Ncf,N_bootstrap)
		(bootstrap_average4,error_G4)=AverageErrorG(G4,N,Ncf,N_bootstrap)
		(bootstrap_average5,error_G5)=AverageErrorG(G5,N,Ncf,N_bootstrap)
		error_G=np.zeros(N)
		bootstrap_average=np.zeros(N)
		for j in range(0,N):
			bootstrap_average[j]=(bootstrap_average5[j]-bootstrap_average2[j]**2)-2*(bootstrap_average1[j]-bootstrap_average2[j]*bootstrap_average3[j])+(bootstrap_average4[j]-bootstrap_average3[j]**2)
			error_G[j]=error_G5[j]+2*error_G2[j]*bootstrap_average2[j]+error_G1[j]+(error_G2[j]/bootstrap_average2[j]+error_G3[j]/bootstrap_average3[j])*bootstrap_average2[j]*bootstrap_average3[j]+error_G4[j]+2*error_G3[j]*bootstrap_average3[j]
		return (average_G,bootstrap_average,error_G)





