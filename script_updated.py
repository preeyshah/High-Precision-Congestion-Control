import os
import sys
import math
from numpy import sort
from numpy import array
filename = sys.argv[1]

# os.system('cat '+filename+' | grep DATA_RATE > temp1')
# tfile = open('temp1', 'r')
# filedata = ''
# for input in tfile:
# 	filedata += input

# filedata = filedata.split('\n')
# dr = filedata[0].split('\t')[-1]
# datarate = int(dr[0:-4])
# os.system('rm temp1')
# print(datarate)

datarate = 100.0

Incast_size = -1
os.system('cat '+filename+' | grep \'Incast info\' > temp1')
tfile = open('temp1', 'r')
filedata = ''
for input in tfile:
	filedata += input

if filedata=='':
	Incast_size = -1
else:
	Incast_size = int((filedata.split('Num')[0]).split(' ')[-1])

print(Incast_size)

is_w3 = False
os.system('cat '+filename+' | grep WORKLOAD > temp1')
tfile = open('temp1', 'r')
filedata = ''
for input in tfile:
	filedata += input

if filedata=='':
	is_w3= False
else:
	#print(filedata)
	workload = filedata.split('\t')[-1]
	if workload[0:2]=='W3':
		is_w3 = True
is_w3 = True
print(is_w3)

# FLOW_SIZE_IN_INCAST
# filedata
os.system('rm temp1')
os.system('cat '+filename+' | grep port > temp1')
tfile = open('temp1', 'r')
filedata = ''
for input in tfile:
	filedata += input

num_switches = 209
filedata = filedata.split('\n')
arr = []
arr2 = []
port_dict = [{} for _ in range(num_switches+1)]
dc_dict = [{} for _ in range(num_switches+1)]
size_dict = [{} for _ in range(num_switches+1)]
started = 0
finished = 0
failed = 0
completed =0
for i in range(len(filedata)):
	if filedata[i][0:8]=='Finished':
		try:
			finished += 1
			sep = filedata[i].split(' ')
			src =int((sep[3].split('.'))[1])
			#print('src')
			dst = int((sep[5].split('.'))[1])
			#print('dst')
			time = int(((sep[8].split('.'))[0]).split('+')[1])
			#print('time')
			size = int(sep[12])
			#print('size')
			ld = 8000
			ex = 3
			np = int(size/1000)
			if size < 1000:
				np = 1
			if size > 1000:
				size = 1000
			if int((src)/16)==int((dst)/16):
				ld = 4000
				ex =1

			base_time = ld+(size*8.0/datarate)*(np+ex)
			ratio = time*1.0/base_time
			
			arr.append([np,time, ratio])
		except:
			continue
	elif filedata[i][0:7]=='Started':
		started += 1

n = 10
nums = [0.0 for _ in range(n)]
sums = [0.0 for _ in range(n)]
maxims = [0.0 for _ in range(n)]
medians = [0.0 for _ in range(n)]
avgs = [0.0 for _ in range(n)]
percentiles = [0.0 for _ in range(n)]
percentiles_95 = [0.0 for _ in range(n)]
vals = [[] for _ in range(n)]

for i in range(len(arr)):
	j = -1
	if not is_w3:
		if arr[i][0]==Incast_size:
			j = 9
		elif arr[i][0] <= 3:
			j = 0
		elif arr[i][0] <= 12:
			j = 1
		elif arr[i][0] <= 48:
			j = 2
		elif arr[i][0] <= 192:
			j = 3
		elif arr[i][0] <= 768:
			j = 4
		elif arr[i][0] <= 3072:
			j = 5
		elif arr[i][0] <= 12288:
			j = 6
		elif arr[i][0] <= 49152:
			j = 7
		else:
			j = 8

	else:
		if arr[i][0]==Incast_size:
			j = 9
		elif arr[i][0] <= 1:
			j = 0
		elif arr[i][0] <= 2:
			j = 1
		elif arr[i][0] <= 4:
			j = 2
		elif arr[i][0] <= 8:
			j = 3
		elif arr[i][0] <= 16:
			j = 4
		elif arr[i][0] <= 64:
			j = 5
		elif arr[i][0] <= 256:
			j = 6
		elif arr[i][0] <= 1024:
			j = 7
		else:
			j = 8		
	nums[j] +=1
	sums[j] +=arr[i][2]
	vals[j].append(arr[i][2])
	# if arr[i][2] >= 10000:
	# 	print(j)
	# 	print(arr[i][1])
	# 	print("-----")

for i in range(len(medians)):
	if nums[i] != 0:
		s = sort(array(vals[i]))
		percentiles[i] = s[int(99.0/100.0*nums[i])]
		percentiles_95[i] = s[int(95.0/100.0*nums[i])]
		medians[i] = s[int(nums[i]/2)]
		avgs[i] = sums[i]/nums[i]
		maxims[i] = s[int(nums[i]-1)]
		#print(s)
		# print(len(s))


print("Flows started", end=" ")
print(started)
print("Flows completed", end=" ")
print(finished)
print("nums =", end=" ")
print(nums)
print("averages =", end=" ")
print(avgs)
print("medians =", end = " ")
print(medians)
print("95 percentile =", end=" ")
print(percentiles_95)
print("99 percentile =", end=" ")
print(percentiles)
print("Maximum =", end=" ")
print(maxims)
os.system('rm temp1')
	



n = 10
nums = [0.0 for _ in range(n)]
sums = [0.0 for _ in range(n)]
maxims = [0.0 for _ in range(n)]
medians = [0.0 for _ in range(n)]
avgs = [0.0 for _ in range(n)]
percentiles = [0.0 for _ in range(n)]
percentiles_95 = [0.0 for _ in range(n)]
vals = [[] for _ in range(n)]

