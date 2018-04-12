# -*- coding: utf-8 -*-
"""
Created on Wed Apr 11 16:14:27 2018

@author: gghysels
"""


import os
import sys
import matplotlib.pyplot as plt

        
def readData():
    temp = []
    head = []
    seconds = []

    filename = raw_input('Enter file name of data file: ')
    try:
        with open(filename,'r') as f:
            for line in f:
                
                if line != 'END OF DATA FILE OF DATALOGGER FOR WINDOWS\t\t\n':
                    temp = line.split('\t')
                    temp[2] = temp[2].rstrip('\n')
                    seconds.append(float(temp[0]))
                    head.append(float(temp[3]))
    
                else:
                    print 'Data imported succesfully!'

    except Exception, e:
        print "Exception: %s" % str(e)
        sys.exit(1)
        
    return seconds, head, filename


#######################################################################
################################ MAIN #################################
######################################################################
 
#Change working directory to location datafile
os.chdir('C:\Users\gertg\OneDrive - Vrije Universiteit Brussel\Python\SlugTestAnalysis\Voorbeeld Cas')

#Read data, convert pressure to waterpressure, convert date to seconds, reads filename
data = readData()
seconds = data[0]
head = data[1]
filename = data[2]

#convert mm to cm
head[:] = [x/10 for x in head]

minHead = min(head)
t_0 = seconds[head.index(minHead)]

print 'Minimum head is', "%.2f" % minHead, 'cm at time t =', t_0, 's'


#Plot H vs T of entire dataset
plt.close("all")
plt.clf()
fig0 = plt.figure(0)
plt.plot(seconds,head, 'b-')
plt.xlabel('Time [s]')
plt.ylabel('H [cm $H_2O$]')
plt.title('H vs t')
plt.grid()
plt.show()
fig0.savefig(filename[:-4] + '_fig0.pdf', bbox_inches='tight')