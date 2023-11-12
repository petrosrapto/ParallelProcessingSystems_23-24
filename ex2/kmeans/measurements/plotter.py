import matplotlib.pyplot as plt
import sys
import numpy as np
import matplotlib.pyplot as plt

fp = open(sys.argv[1])
line = fp.readline()
count = 0

clusters_ser = []
clusters_1_thread = []
clusters_2_thread = []
clusters_4_thread = []
clusters_8_thread = []
clusters_16_thread = []
clusters_32_thread = []
clusters_64_thread = []

time = []
x_axis = ["Sequential" , "1 thr" , "2 thr" , "4 thr" , "8 thr" , "16 thr" , "32 thr" , "64 thr"]

while line:
    tokens = line.split()
    if tokens == []:
        line = fp.readline()
        continue
    if tokens[0] == 'Final':
        for i in range(16):
            l = fp.readline()
            t = l.split()
            if count == 0:
                clusters_ser.append(list(map(float, t[2:18])))
            elif count == 1:
                clusters_1_thread.append(list(map(float, t[2:18])))
            elif count == 2:
                clusters_2_thread.append(list(map(float, t[2:18])))
            elif count == 3:
                clusters_4_thread.append(list(map(float, t[2:18])))
            elif count == 4:
                clusters_8_thread.append(list(map(float, t[2:18])))
            elif count == 5:
                clusters_16_thread.append(list(map(float, t[2:18])))
            elif count == 6:
                clusters_32_thread.append(list(map(float, t[2:18])))
            elif count == 7:
                clusters_64_thread.append(list(map(float, t[2:18])))
        count += 1
    if tokens[0] == 'nloops':
        s = tokens[5] ### it is of the form '7.7225s)' we need to get rid of the s and the )
        time.append(float(s[:-2]))


    line = fp.readline()

fp.close()

checks_passed = True
for i in range(16):
    for j in range(16):
        if(clusters_ser[i][j] != clusters_1_thread[i][j]):
            print("FOUND A MISTAKE IN 1 THREAD in position" , i , j , "supposed to be :" , clusters_ser[i][j])
            checks_passed = False
        if(clusters_ser[i][j] != clusters_2_thread[i][j]):
            print("FOUND A MISTAKE IN 2 THREAD in position" , i , j , "supposed to be :" , clusters_ser[i][j])
            checks_passed = False
        if(clusters_ser[i][j] != clusters_4_thread[i][j]):
            print("FOUND A MISTAKE IN 4 THREAD in position" , i , j , "supposed to be :" , clusters_ser[i][j])
            checks_passed = False
        if(clusters_ser[i][j] != clusters_8_thread[i][j]):
            print("FOUND A MISTAKE IN 8 THREAD in position" , i , j , "supposed to be :" , clusters_ser[i][j])
            checks_passed = False
        if(clusters_ser[i][j] != clusters_16_thread[i][j]):
            print("FOUND A MISTAKE IN 16 THREAD in position" , i , j , "supposed to be :" , clusters_ser[i][j])
            checks_passed = False
        if(clusters_ser[i][j] != clusters_32_thread[i][j]):
            print("FOUND A MISTAKE IN 32 THREAD in position" , i , j , "supposed to be :" , clusters_ser[i][j])
            checks_passed = False
        if(clusters_ser[i][j] != clusters_64_thread[i][j]):
            print("FOUND A MISTAKE IN 64 THREAD in position" , i , j , "supposed to be :" , clusters_ser[i][j])
            checks_passed = False

if(checks_passed):
    print("Checks passed")
else:
    exit()

print(time)

plt.figure(1)
plt.bar(x_axis, time)
plt.title('Total time')
plt.xlabel('')
plt.ylabel('Time in seconds', fontsize=14)
plt.savefig('total_time.png')

seq_time = time[0]
speedup = []
for x in time:
    speedup.append(seq_time/x)

print(speedup)

plt.figure(2)
plt.bar(x_axis, speedup)
plt.title('Speedup')
plt.xlabel('')
plt.ylabel('Speedup compared to sequential algorithm' , fontsize=14)
plt.savefig('speedup.png')
