import matplotlib.pyplot as plt
import numpy as np
import math
import argparse
from pdb import set_trace
import re

#define parser for arguments
parser=argparse.ArgumentParser(description='Plot Teensy Data')
parser.add_argument('logDir', type=str, help='Path to log file for plotting')
args=parser.parse_args()

#display input path
print(f'Your log path: {args.logDir}')

#load in data
with open(args.logDir, 'r') as file:
    teensyData=file.read()
    
#print
printData=False
if printData:
    print(teensyData)

splitData=teensyData.splitlines()
omegaBias=[float(i) for i in re.findall(r'[+-]?\d+\.\d+', splitData[0])]
accelBias=[float(i) for i in re.findall(r'[+-]?\d+\.\d+', splitData[1])]
splitData=splitData[2:]

#2 spaces every 5 break up time points
timeStride=6
nTm=math.floor(len(splitData)/timeStride)
times=np.zeros((nTm,1))
omegas=np.zeros((nTm,3))
accels=np.zeros((nTm,3))
thetas=np.zeros((nTm,2))

#print('Len Data='+str(nTm))

for iData in range(nTm):
    tmIx=iData*timeStride
    #print('TmIx='+str(tmIx))
    times[iData]=[float(i) for i in re.findall(r'[+-]?\d+\.\d+', splitData[tmIx])]
    omegas[iData,:]=[float(i) for i in re.findall(r'[+-]?\d+\.\d+', splitData[tmIx+1])]
    accels[iData,:]=[float(i) for i in re.findall(r'[+-]?\d+\.\d+', splitData[tmIx+2])]
    thetas[iData,:]=[float(i) for i in re.findall(r'[+-]?\d+\.\d+', splitData[tmIx+3])]

#set_trace()
print('')
print('Gyro Bias=' + str(omegaBias))
print('Accel Bias=' + str(accelBias))

#make plots
plt.figure(1)
plt.plot(times, omegas[:,0], color='blue', linewidth=2, label='omega_x')
plt.plot(times, omegas[:,1], color='black', linewidth=2, label='omega_y')
plt.plot(times, omegas[:,2], color='magenta', linewidth=2, label='omega_z')
plt.xlabel('Time [sec]')
plt.ylabel(r'$\omega [deg/sec]$')
plt.legend(loc='best')
#plt.show()

plt.figure(2)
plt.plot(times, accels[:,0], color='blue', linewidth=2, label='accel_x')
plt.plot(times, accels[:,1], color='black', linewidth=2, label='accel_y')
plt.plot(times, accels[:,2], color='magenta', linewidth=2, label='accel_z')
plt.xlabel('Time [sec]')
plt.ylabel(r'$accel [g]$')
plt.legend(loc='best')
#plt.show()

plt.figure(3)
plt.plot(times, thetas[:,0], color='blue', linewidth=2, label='theta_x')
plt.plot(times, thetas[:,1], color='black', linewidth=2, label='theta_y')
plt.xlabel('Time [sec]')
plt.ylabel(r'$\theta [deg]$')
plt.legend(loc='best')
plt.show()