#coding:utf-8
from math import *
import random
import numpy as np
from scipy import integrate as integrate
from scipy.integrate import simps
from numpy import trapz
from decimal import localcontext
import matplotlib.pyplot as plt
import os.path
from matplotlib import rc
import h5py

N=32
Nt=N
def read_W(filename):
	f = h5py.File(filename, 'r')
	dset=f['data']
	for i in range (0,Nt):
		for j in range (0,Nt):
			print(dset[i][j])


def main():
	read_W("data_W_32_100.000000.h5")


if __name__ == '__main__':
	main()	