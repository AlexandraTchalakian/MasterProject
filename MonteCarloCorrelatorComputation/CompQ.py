from math import *
import numpy as np
import Markov_chain as MC

def computeQ(x,a,N):

	Q=0
	for i in range(0,N):
		if i==N-1:
			Q+=MC.modulo(x[0]-x[N-1],0.5)
		else:
			Q+=MC.modulo(x[i+1]-x[i],0.5)
	return Q



def averageM(A,M,Ncf):
	average=0
	for alpha in range (0,Ncf):
		average+=A[alpha]**M

	average/=float(Ncf)
	return average



def computeAverageQ(x,a,N,Ncf):
	Q=np.zeros(Ncf)	
	for alpha in range(0,Ncf):			
		Q[alpha]=computeQ(x[alpha],a,N)
		#intergral of dx/dt*dt
	average2_Q=averageM(Q,2,Ncf)
	print(averageM(Q,2,Ncf))
	average4_Q=averageM(Q,4,Ncf)

	return (average2_Q, average4_Q)