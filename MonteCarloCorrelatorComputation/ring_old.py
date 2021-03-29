from math import *
import random 
import numpy as np
import matplotlib.pyplot as plt



m=1
Ncor=700
Ncf=3000
Delta=0.5
R=sqrt(1/(4*pi**2*m))
Omega=0


def modulo(x,n):
	if abs(x)<=n:
		return x
	if x > n:
		return x-1
	if x <-n:
		return x+1

def initX(N):
	x=[]
	for i in range(0,N):
		x.append(0)

	return x



def potential(y,type):
	if(type==0):
		return Omega**2*cos(2*pi*y)
	if(type==1):
		return 0.5*y**2
	else:
		return 0


# def S(j,x,a,N):
# 	jp=(j+1)%N
# 	jm=(j-1)%N
# 	x_HO=x[j]*(x[j]-x[jp]-x[jm])

# 	if abs(x[jp]-x[j])<=0.5 and abs(x[j]-x[jm])<=0.5:
# 		x_kin=2*x_HO

# 	if ((x[jp]-x[j]>0.5 and x[j]-x[jm]>0.5) or (x[jp]-x[j]<-0.5 and x[j]-x[jm]<-0.5)):
# 		x_kin=2*(x_HO+1)

# 	if((abs(x[jp]-x[j])<=0.5 and x[j]-x[jm]>0.5) or (x[jp]-x[j]<-0.5 and abs(x[j]-x[jm])<=0.5)):
# 		x_kin=2*(x_HO-x[j])+1	

# 	if ((x[jp]-x[j]>0.5 and abs(x[j]-x[jm])<=0.5) or (abs(x[jp]-x[j])<=0.5 and x[j]-x[jm]<-0.5)):
# 		x_kin=2*(x_HO+x[j])+1

# 	if x[jp]-x[j]>0.5 and x[j]-x[jm]<-0.5:
# 		x_kin=2*(x_HO+2*x[j]+1)

# 	if x[jp]-x[j]<-0.5 and x[j]-x[jm]>0.5:
# 		x_kin=2*(x_HO-2*x[j]+1)
		
# 	return 0.5*m*x_kin/a + a*potential(x[j],0)



def S(j,x,a,N):
	jp=(j+1)%N
	jm=(j-1)%N
	return 0.5*m/a*(modulo((x[jp]-x[j]),0.5)**2+modulo((x[j]-x[jm]),0.5)**2)+a*potential(x[j],0)


def update(x,a,N):
	for i in range(0,N):
		x_before=x[i]
		S_before=S(i,x,a,N)
		x[i]=(x[i]+(1-2*random.uniform(0,1))*Delta)%1
		delta_S=S(i,x,a,N)-S_before
		if delta_S>0 and exp(-delta_S)<random.uniform(0,1):
			x[i]=x_before




def computeQ(x,a,N):

	Q=0
	for i in range(0,N):
		if i==N-1:
			Q+=modulo(x[0]-x[N-1],0.5)
		else:
			Q+=modulo(x[i+1]-x[i],0.5)	
	return Q


def averageM(A,M):
	average=0
	for alpha in range (0,Ncf):
		average+=A[alpha]**M

	average/=float(Ncf)
	return average


def computeAverageQ(x,a,N):
	Q=np.zeros(Ncf)
	t=np.zeros(N)

	for i in range(0,N):
		t[i]=t[i]+a*i
	for i in range(0,5*Ncor):
		update(x,a,N)
	for alpha in range(0,Ncf):
		for j in range (0,Ncor):
			update(x,a,N)	
		#sauver les configuration			
		Q[alpha]=computeQ(x,a,N)
		#intergral of dx/dt*dt
	average2_Q=averageM(Q,2)
	print(averageM(Q,2))
	average4_Q=averageM(Q,4)

	for alpha in range(0,Ncf):
		print(Q[alpha])

	x_0=[]
	x_1=[]
	x_1_minus=[]
	x_2=[]
	x_2_minus=[]
	t_0=[]
	t_1=[]
	t_1_minus=[]
	t_2=[]
	t_2_minus=[]

	for i in range(0,N):
		for alpha in range(0,Ncf):
			if round(Q[alpha])==0:
				print("0")
				x_0.append(x[i])
				t_0.append(t[i])
			if round(Q[alpha])==1:
				print("1")
				x_1.append(Q[alpha])
				t_1.append(t[i])
			if round(Q[alpha])==-1:
				print("-1")
				x_1_minus.append(Q[alpha])
				t_1_minus.append(t[i])
			if round(Q[alpha])==2:
				print("2")
				x_2.append(Q[alpha])
				t_2.append(t[i])
			if round(Q[alpha])==-2:
				print("-2")
				x_2_minus.append(Q[alpha])
				t_2_minus.append(t[i])

	plt.figure(1)
	plt.plot(t_0,x_0, marker='o',linestyle='none')
	plt.plot(t_1,x_1, marker='o',linestyle='none')
	plt.plot(t_1_minus,x_1_minus, marker='o',linestyle='none')
	plt.xlabel("t")
	plt.ylabel('x(t)')
	plt.show()

	plt.figure(2)
	plt.hist(Q, bins=50)
	plt.xlabel("Q")
	plt.show()

	return (average2_Q, average4_Q)

			




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



def main():
	N=[10,20]
	Na=np.size(N)
	a=np.zeros(Na)
	average2_Q=np.zeros(Na)
	average4_Q=np.zeros(Na)
	Q=np.zeros(Na)
	Xi=np.zeros(Na)
	b_2=np.zeros(Na)
	for i in range(0,Na):
		x=initX(N[i]);
		a[i]=1/(5*float(N[i]))
		(average2_Q[i], average4_Q[i])=computeAverageQ(x,a[i],N[i])
		Xi[i]=susceptibility(average2_Q[i],beta(N[i],a[i]))
		b_2[i]=b2(average2_Q[i],average4_Q[i])

	# plt.figure(1)
	# plt.plot(a,Xi,marker='.',linestyle='none')
	# plt.show()

	

if __name__ == '__main__':
	main()