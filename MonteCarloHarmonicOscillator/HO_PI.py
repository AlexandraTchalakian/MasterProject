#coding:utf-8
import math
import matplotlib.pyplot as plt

m=1
N=10
a=0.5
x_0=0
x_N=2
x=[]

n=0
while n<=N:
	j=n/N
	x.append(x_0+(x_N-x_0)*j)
	n+=1



def potential(x_temp,potential_type):
	if potential_type=="HO":
		return x_temp*x_temp*0.5
	elif potential_type=="AHO":
		return x_temp*x_temp*x_temp*x_temp*0.25
	else:
		return 0



def K_j(x0,x1,x2):
	return math.sqrt(m*0.5/(math.pi*a))*math.exp(-(m*a*(x1*x1-x1*(x2+x0))+a*potential(x1,"HO")))


def integal_j(lim_inf,lim_sup, j):
	return (lim_sup-lim_inf)/6*(K_j(x[j-1],lim_inf,x[j+1])+4*K_j(x[j-1],(lim_inf+lim_sup)*0.5,x[j+1])+K_j(x[j-1],lim_sup,x[j+1]))
		

def Graph():
	K_vec=[]
	x_plot=[]
	j=0
	
	while j<N:
		K_vec.append(integal_j(-5,5,j))
		x_plot.append(x[j])
		j+=1

	plt.plot(x_plot,K_vec,"o:")
	plt.show()	


Graph()
