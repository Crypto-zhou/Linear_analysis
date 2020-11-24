# -*- coding: utf-8 -*-
# @Time    : 2020/10/24 14:09
# @Author  : 周瑞生
# @Email   : 201800122031@mail.sdu.edu.cn
# @File    : LAT_table.py
# @Software: PyCharm
import numpy as np
import pandas as pd


def LinearTestS(alpha,beta,x):
	Alpha='{:04b}'.format(alpha)
	Beta='{:04b}'.format(beta)
	X='{:04b}'.format(x)
	Sx='{:04b}'.format(int(S[x],16))
	Res=0
	for i in range(4):
		if Alpha[i]=='1':
			Res=Res^int(X[i],2)
		if Beta[i]=='1':
			Res=Res^int(Sx[i],2)
	return Res^1


def Create_Lat_Table(m,n):
	Lat=np.zeros((2**m,2**n),dtype=np.int)
	for alpha in range(2**m):
		for beta in  range(2**n):
			for x in range(2**m):
				if LinearTestS(alpha,beta,x):
					Lat[alpha][beta]+=1
	return Lat-(2**(m-1))

S=['6','4','c','5','0','7','2','e','1','f','3','d','8','a','9','b']
m=4
n=4
if __name__=='__main__':
	Lat=Create_Lat_Table(m,n)
	filename = 'LatTable'
	with open(filename+'.txt','w',encoding='utf-8',newline='') as f:
		f.write(str(pd.DataFrame(Lat)))
	np.save(filename+'.npy',Lat)
	L=np.load(filename+'.npy')
	print(L)
	pass