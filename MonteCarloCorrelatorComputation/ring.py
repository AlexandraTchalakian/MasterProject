from math import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib import cm
import Markov_chain as MC
import Correlator
import CompQ
import sys

m=1
Ncor=500
Ncf=1000000
Delta=0.5
Omega=0
#Omega=0.5
N=int(sys.argv[2])
a=float(sys.argv[1])





def time(N,a):
	t=[]
	t_0=0
	for i in range(0,N):
		t.append(t_0+a*i)
	return t



def beta(N,a):
	return N*a




def susceptibility(average2_Q,beta):
	return average2_Q/float(beta)



def b2(average2_Q,average4_Q):
	return -(average4_Q - 3*average2_Q**2)/(12*average2_Q)



def PlotTwoPtFunctionSineGordon(x,t,N):
	(two_pt_function1,two_pt_function1_bootstrap,two_pt_function1_error,G1)=Correlator.TwoPtFunction(x,a,N,Ncf,N_bootstrap=100,ndot="ndot",terms=1)
	(two_pt_function1_dot,two_pt_function1_bootstrap_dot,two_pt_function1_error_dot,G1_dot)=Correlator.TwoPtFunction(x,a,N,Ncf,N_bootstrap=100,ndot="dot",terms=1)
	(two_pt_function2,two_pt_function2_bootstrap,two_pt_function2_error)=Correlator.TwoPtFunction(x,a,N,Ncf,N_bootstrap=100,ndot="ndot",terms=2)
	(two_pt_function3,two_pt_function3_bootstrap,two_pt_function3_error)=Correlator.TwoPtFunction(x,a,N,Ncf,N_bootstrap=100,ndot="ndot",terms=3)
	two_pt_function2_dot=Correlator.TwoPtFunction(x,a,N,Ncf,N_bootstrap=100,ndot="dot",terms=2)
	"""
	Q=np.zeros(Ncf)	
	for alpha in range(0,Ncf):			
		Q[alpha]=CompQ.computeQ(x[alpha],a,N)
		print(Q[alpha])

	#logx = np.log(t[1:int(N/2)])
	#logy1 = np.log(two_pt_function3[1:int(N/2)])
	#log_error=np.log(two_pt_function3_error[1:int(N/2)])
	logx = np.log(t[1:10])
	logy1 = np.log(two_pt_function3[1:10])
	log_error=np.log(two_pt_function3_error[1:10])
	coeffs1_log= np.polyfit(logx,logy1,deg=1,w=1/log_error)
	poly1_log= np.poly1d(coeffs1_log)
	yfit1_log= lambda t: np.exp(poly1_log(np.log(t)))
	coeffs1=np.polyfit(t[1:int(N/2)],two_pt_function3[1:int(N/2)],deg=2,w=1/two_pt_function3_error[1:int(N/2)])
	poly1=np.poly1d(coeffs1)
	yfit1= lambda t: poly1(t)
	#print(coeffs1_log)
	logx = np.log(t[N-21:N-1])
	logy1 = np.log(two_pt_function2[N-21:N-1])
	logy2 = np.log(two_pt_function2[N-21:N-1]+two_pt_function2_error[N-21:N-1])
	logy3 = np.log(two_pt_function2[N-21:N-1]-two_pt_function2_error[N-21:N-1])
	coeffs1= np.polyfit(logx,logy1,deg=1)
	coeffs2= np.polyfit(logx,logy2,deg=1)
	coeffs3= np.polyfit(logx,logy3,deg=1)
	poly1 = np.poly1d(coeffs1)
	yfit1 = lambda t: np.exp(poly1(np.log(t)))
	poly2 = np.poly1d(coeffs2)
	yfit2 = lambda t: np.exp(poly2(np.log(t)))
	poly3 = np.poly1d(coeffs3)
	yfit3 = lambda t: np.exp(poly3(np.log(t)))
	print(coeffs1,coeffs2,coeffs3)

	plt.figure(1)
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	#plt.errorbar(t,two_pt_function1,two_pt_function1_error,marker='.',linestyle='none')
	#plt.errorbar(t,two_pt_function1_bootstrap,two_pt_function1_error,marker='.',linestyle='none')
	#line1=plt.errorbar(t[N-21:N-1],two_pt_function2[N-21:N-1],two_pt_function2_error[N-21:N-1],marker='.',linestyle='none')
	#plt.errorbar(t,two_pt_function2,two_pt_function2_error,marker='.',linestyle='none')
	#plt.errorbar(t,two_pt_function3,two_pt_function3_error,marker='.',linestyle='none')
	plt.errorbar(t,two_pt_function3_bootstrap,two_pt_function3_error,marker='.',linestyle='none')
	#plt.errorbar(t[1:20],two_pt_function2[1:20],two_pt_function2_error[1:20],marker='.',linestyle='none')
	#plt.errorbar(t[1:int(N/2)],two_pt_function3[1:int(N/2)],two_pt_function2_error[1:int(N/2)],marker='.',linestyle='none')
	#line2=plt.plot(t[N-21:N-1],yfit1(t[N-21:N-1]),label="g(t)="+str(coeffs1[0])+"t"+str(coeffs1[1]))
	#line2=plt.plot(t[1:int(N/2)],yfit1(t[1:int(N/2)]),label="g(t)="+str(coeffs1[0])+"t²"+str(coeffs1[1])+"t"+str(coeffs1[2]))
	#line2=plt.plot(t[1:int(N/2)],yfit1_log(t[1:int(N/2)]))
	#line2=plt.plot(t[1:10],yfit1_log(t[1:10]))
	#line2=plt.plot(t[1:20],yfit1(t[1:20]),label="g(t)="+str(coeffs1[0])+"t²"+str(coeffs1[1])+"t"+str(coeffs1[2]))
	#line3=plt.plot(t[N-21:N-1],yfit2(t[N-21:N-1]),label="g(t)="+str(coeffs2[0])+"t"+str(coeffs2[1]))
	#line4=plt.plot(t[N-21:N-1],yfit3(t[N-21:N-1]),label="g-(t)="+str(int(coeffs3[0]))+"t"+str(int(coeffs3[1])))
	#plt.errorbar(t,two_pt_functtion2_bootstrap,two_pt_function2_error,marker='.',linestyle='none')
	#plt.loglog(t[N-1-30:N-1],two_pt_function2[N-1-30:N-1],marker='x',linestyle='none')
	#ax=plt.axes()
	#ax.set_yscale('log')
	#ax.set_xscale('log')
	#first_legend=plt.legend(handles=line2,loc='lower right')
	#plt.gca().add_artist(first_legend)
	#plt.legend(handles=line3)
	plt.ylabel(r'$\langle x(t+\Delta t)x(t)\rangle$')
	plt.xlabel(r'$t$')
	plt.title("N="+str(N)+" "+"a="+str(a))
	plt.show()
	
	plt.figure(2)
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	plt.errorbar(t,two_pt_function1_dot,two_pt_function1_error_dot,marker='.',linestyle='none')
	plt.errorbar(t,two_pt_function1_bootstrap_dot,two_pt_function1_error_dot,marker='.',linestyle='none')
	plt.plot(t,two_pt_function2_dot,marker='x',linestyle='none')
	plt.ylabel(r'$\langle x(t+\Delta t)x(t)\rangle$')
	plt.xlabel(r'$t$')
	"""
	return (two_pt_function3,two_pt_function3_bootstrap,two_pt_function3_error)



def AnaliticTwoPtFunction(t,N):
	Q=np.linspace(-10,10,num=21)
	analyticTPF=np.zeros(N)
	for j in range (0,N):
		for alpha in range (0,21):
			analyticTPF[j]+=(t[j]/m-(1/(m*N*a)-4*pi**2*Q[alpha]**2/(N*a)**2)*t[j]**2)*exp(-2*pi**2*m*Q[alpha]**2/(N*a))
	return (analyticTPF)



def main():
	t=time(N,a)
	#x=MC.markov_steps(a,N,Ncf,Ncor,10,Delta,m, Omega,"SineGordon")
	x=MC.markov_steps(a,N,Ncf,Ncor,10,Delta,m, Omega,"free")
	(two_pt_function,two_pt_function_bootstrap,two_pt_function_error)=PlotTwoPtFunctionSineGordon(x,t,N)
	print("ok")
	"""
	N=[100,200]
	Na=np.size(N)

	for i in range(0,Na):
		plt.figure(i)
		t=np.zeros(N[i])
		#x=MC.markov_steps(a[i],N[i],Ncf,Ncor,Delta,m, Omega,"HO");
		x=MC.markov_steps(a[i],N[i],Ncf,Ncor,10,Delta,m, Omega,"free")
		#x=MC.markov_steps(a,N[i],Ncf,Ncor,10,Delta,m, Omega,"SineGordon")
		t=time(N[i],a)
		analyticTPF=AnaliticTwoPtFunction(t,N[i])
		#(two_pt_function3,two_pt_function3_bootstrap,two_pt_function3_error)=PlotTwoPtFunctionSineGordon(x,t,N[i])
		(two_pt_function,two_pt_function_bootstrap,two_pt_function_error)=PlotTwoPtFunctionSineGordon(x,t,N[i])
		plt.rc('text', usetex=True)
		plt.rc('font', family='serif')
		#plt.errorbar(t,two_pt_function3,two_pt_function3_error,marker='.',linestyle='none')
		plt.errorbar(t,two_pt_function,two_pt_function_error,marker='.',linestyle='none')
		plt.plot(t,analyticTPF)
		#plt.plot(t,two_pt_function,marker='.',linestyle='none')
		plt.ylabel(r'$\langle x(t+\Delta t)x(t)\rangle$')
		plt.xlabel(r'$t$')
		plt.show()
	"""


if __name__ == '__main__':
	main()