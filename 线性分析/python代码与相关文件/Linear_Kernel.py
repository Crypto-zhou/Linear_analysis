# -*- coding: utf-8 -*-
# @Time    : 2020/10/25 11:10
# @Author  : 周瑞生
# @Email   : 201800122031@mail.sdu.edu.cn
# @File    : Linear_Kernel.py
# @Software: PyCharm

import numpy as np
import pandas as pd
pd.set_option('display.max_columns', None)#显示所有行
pd.set_option('display.max_rows', None)#设置value的显示长度为100，默认为50
pd.set_option('max_colwidth',100)
import time
from functools import wraps
import random as rd
from treelib import Node, Tree
import copy
import math
import scipy.stats as st





def Hex2Bin(str_h):
	return '{:016b}'.format(int(str_h,16))

def Bin2Hex(str_b):
	return '{:04x}'.format(int(str_b,2))
# 计时
def runtime(func):
	@wraps(func) #修饰一下，保证函数名不改变
	def with_time(*args, **kwargs):
		t = time.time()
		res = func(*args, **kwargs) #返回原函数结果
		print('Run-Time: %f s\n' % (time.time()-t))
		return res
	return with_time



#线性分析代码部分,此主函数仅用到此函数中堆积引理计算方法.
class Linear_Analyse(object):
	def __init__(self):
		self.S_box_iv = ['4','8','6','a','1','3','0','5','c','e','d','f','2','b','7','9']
		self.CryM=self.ImportCryText()
		self.L = np.load('LatTable.npy')
		self.l = 4  # 要恢复的密钥比特数
		self.r = 2  # 选择前2**r个可能的密钥
		self.Plain=[] #进行测试需要的明文编号
		self.CheckKey=[2,5,13,12]
		self.hull=[]#[alpha,beta,bias]
		self.K5num = None
		self.LastCount=np.zeros((4,16))
		self.K5TestHullNum = np.zeros(4)
	def ImportCryText(self):#导入明密文对数组
		with open('CryText.txt', 'r') as f:  # 处理输入文件为一个np数组
			reader = f.read()
			reader = (reader.replace('\"', '')).split(',')
			CryM = np.array(reader)
		return CryM

	def Find_live_S(self,str_t):
		S = list(str_t)
		for i, v in enumerate(S): S[i] = int(v, 16)
		Sn = np.array(S)
		return np.nonzero(Sn)[0][0]

	@runtime
	def Linear_Main(self,LinearHull,Count):
		'''
		线性分析的主要思路部分
		:param LinearHull: 所有线性壳
		:param Count: 选择明文数
		:return: 无
		'''
		self.Plain= rd.sample(range(0, 2 ** 16), Count)
		for hull in LinearHull:
			self.K5num = self.Find_live_S(hull[1])
			if math.fabs(float(hull[2]))>=0.02 or (self.K5num==2 and math.fabs(float(hull[2]))>=0.015):
				self.Counter_Construct(hull,Count)
		print('K5的4个密钥的线性壳使用次数:',self.K5TestHullNum)
		with open('LinearRes_test.txt', 'a+', encoding='utf-8', newline='') as f:
			f.write("\n使用明文数为"+str(Count)+'时的计数结果:\n'+str(pd.DataFrame(self.LastCount)))
			K5=np.argmax(self.LastCount, axis=1)
			f.write("\n恢复出来的密钥为:"+str(K5))
			print("\n恢复出来的密钥为:"+str(K5))

	def StackingLemma(self, Eps):
		'堆积引理；eps:偏差列表；return 偏差'
		if Eps == []:
			return -1
		n = len(Eps)
		P = np.prod(Eps)
		# print(n,P)# 快速连乘积
		return (2 ** (n - 1))* P

	def int2arr(self,Num,Len=16):
		'''
		:param Num: 需要进行转换的十进制数字
		:param Len: 需要转换成数组的长度
		:return: 转换后的数组
		'''
		Arr=np.zeros(Len,dtype=int)
		Str=bin(Num)[2:]
		L=len(Str)
		for i in range(L):
			Arr[Len-1-i]=Str[L-1-i]
		return Arr
	def strH2arr(self,Str,Len=16):
		'''
			:param Str: 需要进行转换的十六进制字符串,长度必须为4
			:param Len: 需要转换成数组的长度
			:return: 转换后的数组
		'''
		return self.int2arr(int(Str,16),16)
	def Counter_Construct(self,hull,DtSize):
		'''
		:param hull: 线性壳['输入掩码','输出掩码','偏差']
		:param DtSize: 选择的明文个数
		:return: 偏差的计数器,取绝对值(此方法进行对对应线性壳进行计数)
		'''
		self.hull=hull
		if self.K5num==None:
			self.K5num = self.Find_live_S(hull[1])
		Counter0 = np.zeros(2 ** self.l)
		Alpha=self.strH2arr(self.hull[0],16)
		Beta = self.strH2arr(self.hull[1], 16)
		for plainText in self.Plain:
			plainText_arr=self.int2arr(plainText,16)
			CryText = self.CryM[plainText]
			for K5_ in range(16):
				Temp = int(CryText[self.K5num], 16) ^ K5_  # 由密文异或K50
				OutArr=self.int2arr(int('0'*self.K5num+self.S_box_iv[Temp]+'0'*(3-self.K5num), 16),16)
				if (Alpha@plainText_arr)%2==(Beta@OutArr)%2:
					Counter0[K5_] += 1
		Counter0=Counter0-DtSize//2
		Counter0=abs(Counter0)
		if self.IsHullValid(Counter0,DtSize):
			self.K5TestHullNum[self.K5num] += 1
			Sorted = np.argsort(-Counter0)[:4]
			for index in Sorted:
				self.LastCount[self.K5num][index]+=1
			# print(hull,"\n此时恢复K5中第",self.K5num,'密钥,排名前4位的是:',Sorted)
		return Counter0


	def Computer_CDF(self,hull,Count):
		'''
		:param hull: 线性壳['输入掩码','输出掩码','偏差']
		:param Count: 选择的明文个数
		:return: 无(此方法用于利用大数定理(伯努利分布)计算CDF函数值)
		'''
		self.Plain= rd.sample(range(0, 2 ** 16), Count)
		Counter=self.Counter_Construct(hull,Count)
		Counter_Gauss=np.zeros(16)
		for i in range(16):
			Counter_Gauss[i]=st.norm.cdf(Count*0.5-Counter[i],loc=Count*0.5,scale=1/2*math.sqrt(Count))
		print('线性壳为:',hull,'时,利用大数定律计算对应偏差出现的概率:\n',Counter_Gauss)


	def IsHullValid(self,Counter,Count):
		'''
		此方法用于判断线性壳是否可以使用,对应于第二个筛选条件
		:param Counter: 偏差数的计数器
		:param Count: 选择明密文的数目
		:return: 是(T)否(F)可以
		'''
		SelectedIndex = np.argsort(-Counter)[0]
		Counter2=copy.deepcopy(Counter)
		Counter22=np.delete(Counter2,SelectedIndex)
		for i in range(len(Counter22)):#如果有两个或以上最大值重复,直接弃用
			if Counter22[i]>=(Counter[SelectedIndex]-0.0001*Count):
				return False
		return True

Res_hull=[]
#寻找线性壳与概率函数
class FindLinearHull(Linear_Analyse):
	def __init__(self,alpha,beta):
		self.alpha=alpha
		self.beta=beta
		super().__init__() #父类初始化函数继承
		self.P = ['0', '4', '8', 'c', '1', '5', '9', 'd', '2', '6', 'a', 'e', '3', '7', 'b', 'f']
		self.FMaxRound = 2 #可以通过选择调整这两个(+self.BMaxRound)参数进行调整寻找线性壳的轮数
		self.BMaxRound = 2
		self.ID_ALL=[]
		self.LAT=np.load('LatTable.npy')
		self.threshold=1/128#选择的阈值,为经过测试后权衡准确率与速度的结果
		self.LeavesID=[]
		self.FLeavesID = []
	def Main_Function(self):
		#函数功能:
		#1.根据传入参数中输入掩码,输出掩码进行创建多叉树
		#2.对两个树进行中间相遇处理,找到四轮线性传递路线
		#3.将路线与树拓扑图保存到两个输出文件中,并计算线性壳偏差
		treeF = Tree()
		ID0 = '{:06x}'.format(rd.randint(0, 1000))  # 创建随机数ID
		self.ID_ALL.append(ID0)  # 添加到ID_ALL中,为了与之后创建的ID作区分
		treeF.create_node(self.alpha, ID0, data=[])  # 创建前向树的父节点
		self.FLeavesID.append(ID0)
		treeF = self.ForwardTree(treeF)
		self.ID_ALL.clear()
		treeB = Tree()  # 实例化后向树
		ID0 = '{:06x}'.format(rd.randint(0, 1000))  # 创建随机数ID
		self.ID_ALL.append(ID0)  # 添加到ID_ALL中,为了与之后创建的ID作区分
		treeB.create_node(self.beta, ID0, data=[])  # 创建后向树的父节点
		self.LeavesID.append(ID0)
		treeB=self.BackwardTree(treeB)  # 创建后向树函数
		Res=np.array(self.PathInformation(treeB, treeF))
		Sum=sum(Res)
		print("总共有",len(Res),'条线路')
		print(self.alpha, '-->', self.beta, "的线性壳对应的概率为:", Sum)
		Res_hull.append([self.alpha,self.beta,Sum])
		# treeF.save2file('treeForward.txt')
		# treeB.save2file('treeBackward.txt')
	def Perm(self,str_h):
		output = ['0'] * 16
		tmpS2 = Hex2Bin(str_h)  # 转换成二进制数
		for i in range(16):
			output[int(self.P[i], 16)] = tmpS2[i]  # 二进制数的P置换
		return Bin2Hex(''.join(output))  # 16长的二进制数的列表转换成4长16进制字符串


	#根据根节点创建前向线性多叉树函数
	def ForwardTree(self,tree):
		for round in range(self.FMaxRound):
			TempLeaves=copy.deepcopy(self.FLeavesID)
			self.FLeavesID.clear()
			for LeavesID_tmp in TempLeaves:
				NodeName = tree[LeavesID_tmp].tag
				Index = []
				OutMask=[]
				for i in NodeName:
					Index.append(int(i, 16))
					OutMask.append(np.where(self.LAT[int(i, 16),:] != 0))
				for i1 in OutMask[0][0]:
					for i2 in OutMask[1][0]:
						for i3 in OutMask[2][0]:
							for i4 in OutMask[3][0]:
								NewData = copy.deepcopy(tree[LeavesID_tmp].data)
								if self.LAT[Index[0]][i1] != 8:
									NewData.append(self.LAT[Index[0]][i1] / 16)
								if self.LAT[Index[1]][i2] != 8:
									# print(self.LAT[Index[1]][i2] / 16)
									NewData.append(self.LAT[Index[1]][i2] / 16)
								if self.LAT[Index[2]][i3] != 8:
									# print(self.LAT[Index[2]][i3] / 16)
									NewData.append(self.LAT[Index[2]][i3] / 16)
								if self.LAT[Index[3]][i4] != 8:
									# print(self.LAT[Index[3]][i4] / 16)
									NewData.append(self.LAT[Index[3]][i4]/16)
								Lemma = self.StackingLemma(NewData)
								if max(Lemma, -Lemma) < self.threshold:
									# print('终止',NewData,Lemma)
									continue
								while True:
									IDt = '{:06x}'.format(rd.randint(0, 100000))
									if IDt not in self.ID_ALL:
										self.ID_ALL.append(IDt)
										break
								NodeName_t = hex(i1)[2:] + hex(i2)[2:] + hex(i3)[2:] + hex(i4)[2:]
								NodeName_t = self.Perm(NodeName_t)  # 进入P盒

								tree.create_node(NodeName_t, IDt, parent=LeavesID_tmp, data=NewData)

								self.FLeavesID.append(IDt)
		return tree

	# 根据根节点创建后向线性多叉树函数
	def BackwardTree(self,tree):
		for round in range(self.BMaxRound):
			TempLeaves=copy.deepcopy(self.LeavesID)
			self.LeavesID.clear()
			for LeavesID_tmp in TempLeaves:
				NodeName = tree[LeavesID_tmp].tag
				Index = []
				OutMask=[]
				P_NodeName = self.Perm(NodeName)  # 进入P盒
				for i in P_NodeName:
					Index.append(int(i, 16))
					OutMask.append(np.where(self.LAT[:, int(i, 16)] != 0))
				for i1 in OutMask[0][0]:
					for i2 in OutMask[1][0]:
						for i3 in OutMask[2][0]:
							for i4 in OutMask[3][0]:
								NewData = copy.deepcopy(tree[LeavesID_tmp].data)
								if self.LAT[i1][Index[0]] != 8:
									NewData.append(self.LAT[i1][Index[0]] / 16)
								if self.LAT[i2][Index[1]] != 8:
									NewData.append(self.LAT[i2][Index[1]] / 16)
								if self.LAT[i3][Index[2]] != 8:
									NewData.append(self.LAT[i3][Index[2]] / 16)
								if self.LAT[i4][Index[3]] != 8:
									NewData.append(self.LAT[i4][Index[3]] / 16)
								Lemma=self.StackingLemma(NewData)
								if max(Lemma,-Lemma) < self.threshold:
									continue
								while True:
									IDt = '{:06x}'.format(rd.randint(0, 100000))
									if IDt not in self.ID_ALL:
										self.ID_ALL.append(IDt)
										break
								NodeName_t = hex(i1)[2:] + hex(i2)[2:] + hex(i3)[2:] + hex(i4)[2:]
								tree.create_node(NodeName_t, IDt, parent=LeavesID_tmp, data=NewData)
								self.LeavesID.append(IDt)
		return tree

	#对两个树进行中间相遇处理,得到相遇的路径,以及相关信息的处理
	def PathInformation(self,treeB, treeF):
		Res=[]
		PathsF = np.array(treeF.paths_to_leaves(),dtype=object)  # 将到叶节点的ID路径转换成np.array存储起来
		EndFID=[]
		EndIndexF=[]
		EndFTag=[]
		for i in range(len(PathsF)):
			if len(PathsF[i])!=self.FMaxRound+1:
				continue
			else:
				EndFID.append(PathsF[i][self.FMaxRound])
				EndIndexF.append(i)
				EndFTag.append(treeF[PathsF[i][self.FMaxRound]].tag)
		EndBID = []
		EndIndexB=[]
		EndBTag=[]
		PathsB = np.array(treeB.paths_to_leaves(),dtype=object)
		for i in range(len(PathsB)):
			if len(PathsB[i]) != self.BMaxRound + 1:
				continue
			else:
				EndBID.append(PathsB[i][self.BMaxRound])
				EndIndexB.append(i)
				EndBTag.append(treeB[PathsB[i][self.BMaxRound]].tag)
		for i in EndFTag:
			if i in EndBTag:
				# print("Find it!,相遇结点的tag为:", i, "此相遇的结点对应路径数:",
				#       EndBTag.count(i) * EndFTag.count(i))
				IndexB = []
				IndexF = []
				for iNum in range(len(EndBTag)):
					if i == EndBTag[iNum]:
						IndexB.append(iNum)
				for iNum in range(len(EndFTag)):
					if i == EndFTag[iNum]:
						IndexF.append(iNum)
				for i1 in IndexF:
					for j1 in IndexB:
						L1=treeF[EndFID[i1]].data
						L2=treeB[EndBID[j1]].data
						L2.extend(L1)
						Pr=self.StackingLemma(L2)
						L1.clear()
						L2.clear()
						Res.append(Pr)
		return Res


VarDiff_list=[]
if __name__=='__main__':
	#=========线性壳生成部分BEGIN
	T1=time.time()#设计失误,导致计时功能不能放到装饰器函数中.....写代码前,没考虑好
	for i in range(16):
		alpha='{:04x}'.format(1<<i)
		for j in range(16):
			beta = '{:04x}'.format(1 << j)
			L = FindLinearHull(alpha, beta)
			L.Main_Function()
	T2=time.time()
	with open('LinearHull.txt','w',encoding='utf-8',newline='') as f:
		f.write(str(np.array(Res_hull)))#方便观察,后续线性分析代码中不再使用
	print("测试256对线性壳寻找算法,总共用时:", T2-T1 ,"s")
	np.save("LinearHull.npy",np.array(Res_hull)) #方便线性分析继续使用
	# =========线性壳生成部分代码END


	#============线性分析
	 # p=0.03+0.52  2|p-1/2|-2

	# Count_Sample=[2222,8888]
	# LinearHull=np.load("LinearHull.npy")
	# for i in range(2):
	# 	VarDiff_list = []
	# 	l=Linear_Analyse()#实例化类

	# 	l.Linear_Main(LinearHull,Count_Sample[i])

	#============线性分析

