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
Ncor=1000
Ncf=1000
Delta=0.5
Omega=0.5
N=int(sys.argv[2])
a=float(sys.argv[1])
N_bootstrap=100
Beta=N*a
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 16}




def time(N,a):
	t=[]
	t_0=0
	for i in range(0,N):
		t.append(t_0+a*i)
	return t

def AnaliticTwoPtFunction(t,N):
	Q=np.linspace(-10,10,num=21)
	analyticTPF=np.zeros(N)
	for j in range (0,N):
		analyticTPF[j]=t[j]/m-1/(m*Beta)*t[j]**2
		for alpha in range (0,21):
			analyticTPF[j]+=0
			#analyticTPF[j]+=(t[j]/m-(1/(m*Beta)-Q[alpha]**2/(Beta)**2)*t[j]**2)*exp(-m*Q[alpha]**2/(2*Beta))
	return (analyticTPF)



def PlotTwoPtFunctionSineGordon(x,t,N):
	(two_pt_function,two_pt_function_error)=CS.TwoPtFunction(x,a,N,Ncf,N_bootstrap)
	analyticTPF=AnaliticTwoPtFunction(t,N)
	coeffs1=np.polyfit(t,two_pt_function,deg=2)
	poly1=np.poly1d(coeffs1)
	yfit1= lambda t: poly1(t)
	print(coeffs1)
	plt.figure()
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	plt.rc('font',**font) 
	line1=plt.errorbar(t,two_pt_function,two_pt_function_error,marker='.',linestyle='none',label='Numerical solution')
	line2=plt.plot(t,yfit1(t),label=str(round(coeffs1[0],2))+"t²+"+str(round(coeffs1[1],2))+"t")
	line3=plt.plot(t,analyticTPF,label='Analytic solution for V=0')
	plt.ylabel(r'$\langle Q^2(t)\rangle$')
	plt.xlabel(r'$t$[s]')
	plt.legend(handles=[line1,line2[0],line3[0]])
	plt.title('N='+str(N))
	plt.show()



def PlotTwoPtFunctionSineGordonDot(x,t,N):
	(two_pt_function,two_pt_function_error)=CS.TwoPtFunctionVelocities(x,a,N,Ncf,N_bootstrap)
	analyticTPF=AnaliticTwoPtFunction(t,N)

	plt.figure()
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	plt.rc('font',**font) 
	line1=plt.errorbar(t,two_pt_function,two_pt_function_error,marker='.',linestyle='none',label='Numerical solution')

	line3=plt.plot(t,analyticTPF,label='Analytic solution for V=0')
	plt.ylabel(r'$\langle Q^2(t)\rangle$')
	plt.xlabel(r'$t$[s]')
	plt.title('N='+str(N))
	plt.show()


def main():
	t=time(N,a)
	x=MC.markov_steps(a,N,Ncf,Ncor,10,Delta,m, Omega,"SineGordon")
	#x=MC.markov_steps(a,N,Ncf,Ncor,10,Delta,m, Omega,"paper")
	PlotTwoPtFunctionSineGordon(x,t,N)
	PlotTwoPtFunctionSineGordonDot(x,t,N)
	
	

if __name__ == '__main__':
	main()
