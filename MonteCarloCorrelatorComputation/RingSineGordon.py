from math import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib import cm
import Markov_chain as MC
import CorrelatorSineGordonFinal as CS
import CompQ
import sys

m=1
Ncor=500
Ncf=500
Delta=0.5
Omega=0.5
N=int(sys.argv[2])
a=float(sys.argv[1])
N_bootstrap=100




def time(N,a):
	t=[]
	t_0=0
	for i in range(0,N):
		t.append(t_0+a*i)
	return t




def PlotTwoPtFunctionSineGordon(x,t,N):
	(two_pt_function,two_pt_function_error)=CS.TwoPtFunction(x,a,N,Ncf,N_bootstrap)

	coeffs1=np.polyfit(t,two_pt_function,deg=2,w=1/two_pt_function_error)
	poly1=np.poly1d(coeffs1)
	yfit1= lambda t: poly1(t)
	print(coeffs1)

	plt.figure()
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	plt.rc('font',**font) 
	plt.errorbar(t,two_pt_function,two_pt_function_error,marker='.',linestyle='none')
	plt.plot(t,yfit1(t))
	plt.ylabel(r'$\langle Q^2(t)\rangle$')
	plt.xlabel(r'$t$[s]')
	plt.show()




def main():
	t=time(N,a)
	x=MC.markov_steps(a,N,Ncf,Ncor,10,Delta,m, Omega,"SineGordon")
	PlotTwoPtFunctionSineGordon(x,t,N)
	
	

if __name__ == '__main__':
	main()
