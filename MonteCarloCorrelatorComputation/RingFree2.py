from math import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib import cm
import Markov_chain as MC
import CorrelatorFree as CF
import CompQ
import sys

m=1
Ncor=500
Ncf=500
Delta=0.5
Omega=0
N=int(sys.argv[2])
a=float(sys.argv[1])
N_bootstrap=100
Beta=N*a
R=1/sqrt(4*pi**2*m)


def time(N,a):
	t=[]
	t_0=0
	for i in range(0,N):
		t.append(t_0+a*i)
	return t




def PlotTwoPtFunctionSineGordon(x,t,N):
	analyticTPF=AnaliticTwoPtFunction(t,N)
	(two_pt_function,two_pt_function_b,two_pt_function_error)=CF.TwoPtFunction(x,a,N,Ncf,N_bootstrap)

	plt.figure()
	plt.errorbar(np.asarray(t)+a,two_pt_function,two_pt_function_error,marker='.',linestyle='none')
	plt.plot(t,analyticTPF)
	plt.show()




def AnaliticTwoPtFunction(t,N):
	Q=np.linspace(-10,10,num=21)
	analyticTPF=np.zeros(N)
	for j in range (0,N):
		analyticTPF[j]=t[j]/m-1/(m*Beta)*t[j]**2
		for alpha in range (0,21):
			analyticTPF[j]+=0
			#analyticTPF[j]+=(t[j]/m-(1/(m*Beta)-4*pi**2*Q[alpha]**2/(Beta)**2)*t[j]**2)*exp(-2*pi**2*m*Q[alpha]**2/(Beta))
	return (analyticTPF)



def main():
	t=time(N,a)
	x=MC.markov_steps(a,N,Ncf,Ncor,10,Delta,m, Omega,"free")
	PlotTwoPtFunctionSineGordon(x,t,N)
	
	

if __name__ == '__main__':
	main()