# -*- coding: utf-8 -*-
"""
Created on Wed Apr 11 16:34:26 2018

@author: gghysels
"""


#importing packages
import os
import sys
import time
import numpy as np
import matplotlib.pyplot as plt
import math
from matplotlib.ticker import FormatStrFormatter



    

#asks for input file, opens inputfile, removes header and stores date, pressure and temperature       
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
    
# Polynomial Regression
def polyFit(x, y, degree):
    results = {}

    coeffs = np.polyfit(x, y, degree)
     # Polynomial Coefficients
    results['polynomial'] = coeffs.tolist()

    correlation = np.corrcoef(x, y)[0,1]

     # r
    results['correlation'] = correlation
     # r-squared
    results['determination'] = correlation**2

    return results


#function that does a linear regression of the data    
def linearRegression(head, seconds):
    
    results = {}

    poly = polyFit(seconds, head,1)

    a = poly['polynomial'][0]
    a_str = "%.9f" % a

    b = poly['polynomial'][1]
    b_str = "%.8f" % b

    x_fit = np.array(range(0,len(slugTestSeconds)+1))
    y_fit = eval('a*x_fit+b')
    y_exp = []
    y_exp[:] = [math.exp(y) for y in y_fit]
    #plt.plot(x_fit,y_exp, '--k', label='Regression (All Data)')

    r2 = poly['determination']
    r2_str = "%.5f" % r2

    regForm = 'y = ', a_str, ' x + ', b_str
    regForm_str = ''.join(regForm)
    

    results['a'] = a
    results['b'] = b
    results['r2'] = r2_str
    results['regForm'] = regForm_str
    results['x'] = x_fit
    results['y'] = y_exp

    
    return results


#Function that returns all values in a list between a max and min value
def sliceList(head, seconds, begin, end):
    n = 0
    r = 0
    t = 0
    partHead = []
    partLnHead = []
    results = {}

    while n < lenList:
        if head[n] < begin:
            if head[n] >end:
                partHead.append(head[n])        
        n = n+1

    while r < len(partHead):
        partLnHead.append(math.log(partHead[r]))
        r = r+1

    while t < len(head):

        if head[t] > begin:
            startTest = 0
            t = t+1
        else:
            startTest = t
            t = len(head)
    
    partSeconds = seconds[0:len(partHead)]
    partSeconds2 = [sec+startTest/2 for sec in partSeconds]
                
    results['partHead'] = partHead
    results['partLnHead'] = partLnHead    
    results['partSeconds'] = partSeconds2

    return results

    

def calculateKhBouwerRice(rc, rw, Le, T0, T0plus, T0first, d, Bt, aniso):
    
    results = {}
    
    b = Le
    r_asterisk = rw*math.sqrt(aniso)
    #print 'rw*: ', r_asterisk
    b_r = b/r_asterisk
    #print 'b/rw*: ', b_r
    
    A = 1.4720 + 0.03537*(b_r)-0.00008148*(b_r)**2+0.0000001028*(b_r)**3-0.00000000006484*(b_r)**4+0.00000000000001573*(b_r)**5
    B = 0.2372 + 0.005151*(b_r) - 0.000002682*(b_r)**2 - 0.0000000003491*(b_r)**3 + 0.0000000000004738*(b_r)**4

    #print 'A: ', A, 'B: ',B
    
    #test = (Bthick-(d+b))/r_asterisk
    
    coeff = math.log((Bt-(d+b))/r_asterisk)
    
    #print 'coeff: ', coeff
    
    if coeff > 6.0:
        coeff = 6.0
    
    coeff2 = ((1.1)/(math.log((d+b)/r_asterisk))+((A+(B*coeff))/b_r))**-1
    #print 'coeff2: ', coeff2    
    
    Kh_cms = (rc*rc*coeff2)/(2*b*T0)
    Kh_md = Kh_cms*864

    Kh_cms_plus = (rc*rc*coeff2)/(2*b*T0plus)
    Kh_md_plus = Kh_cms_plus*864

    Kh_cms_first = (rc*rc*coeff2)/(2*b*T0first)
    Kh_md_first = Kh_cms_first*864
    
    #print 'Kh_cms_plus:', Kh_cms_plus
    #print 'Kh_md_plus:', Kh_md_plus
    
    results['Kh_cms'] = Kh_cms
    results['Kh_md'] = Kh_md

    results['Kh_cms_plus'] = Kh_cms_plus
    results['Kh_md_plus'] = Kh_md_plus

    results['Kh_cms_first'] = Kh_cms_first
    results['Kh_md_first'] = Kh_md_first
    
    results['r_asterisk'] = r_asterisk
    results['A'] = A
    results['B'] = B
    results['coeff'] = coeff
    results['coeff2'] = coeff2
    
    return results

        
#calculates Kh for a confined aquifer, Hvorslev method, overdamped response (log), partially penetrating   
def calculateKhHvorslev(rc, rw, Le, Re, T0, T0plus, T0first, aniso):
    results = {}
    
    psi = math.sqrt(aniso)/((Le/rw))
    
    #print 'PSI ', psi

    F = ((1/(2*psi)) + math.sqrt(1 + (1/(2*psi)**2)))
    
    #print 'F ', F
    
    Kh_cms_full = ((rc*rc)*math.log(F)) / (2*Le*T0)
    Kh_md_full = Kh_cms_full * 864

    Kh_cms_part = ((rc*rc)*math.log(F)) / (2*Le*T0plus)
    Kh_md_part = Kh_cms_part * 864

    Kh_cms_first = ((rc*rc)*math.log(F)) / (2*Le*T0first)
    Kh_md_first = Kh_cms_first * 864    
    
    results['Kh_cms_full'] = Kh_cms_full
    results['Kh_md_full'] = Kh_md_full

    results['Kh_cms_part'] = Kh_cms_part
    results['Kh_md_part'] = Kh_md_part

    results['Kh_cms_first'] = Kh_cms_first
    results['Kh_md_first'] = Kh_md_first
    
    return results



#Variables
startTime = time.time()
i = 0 #counters
j = 0
k = 0
l = 0
m = 0
n = 0
o = 0
r = 0
t = 0
u = 0
head = [] #empty list

#######################################################################
################################ MAIN #################################
######################################################################


### PARAMETERS ###

rc = 1.658 #effective radius well casing (1.685), corrected for radius transponder (0.3)
Le = 21 #effective screen length = b
Re = 1.535 #effective radius slug test = effective radius well screen (see fetter p. 194)
#Re = 200*rc (this is only for fully penetrating??)
rw = 1.535 #effective radius well screen
d = 48 #z-position of top of screen (positive direction downward)
Bt = 1000 #aquifer thickness (taken from model Aa)
aniso = 1 #anisotropy ratio Kv/Kh
aniso2 = 0.1

os.chdir('C:\Users\gghysels\OneDrive - Vrije Universiteit Brussel\Python\SlugTestAnalysis\Voorbeeld Cas')

#Read data, convert pressure to waterpressure, convert date to seconds, reads filename
data = readData()
seconds = data[0]
head = data[1]
filename = data[2]

#convert mm to cm
head[:] = [x/10 for x in head]

minHead = min(head)
t_0 = seconds[head.index(minHead)]

print 'Minimum head is', "%.2f" % minHead, 'at time t=', t_0, 'seconds'


#Calculates base level (before initiation of slug)
print 'Identify moment of slug test initiation (last head value at base lvl):'

while j < 25:
    print '#',j, ":", head[head.index(minHead)-j]
    j = j + 1
    
stableLvl = input('Enter index of last head value before initiation slug test: ')
baseLvl = np.mean(head[head.index(minHead)-stableLvl-10:head.index(minHead)-stableLvl])
print 'Average base level is:',  "%.2f" % baseLvl, 'cm'

H_0 = baseLvl - minHead
print 'H_0 =', "%.2f" %H_0, 'cm'


#Plot H vs T of entire dataset
plt.close("all")
plt.clf()
fig0 = plt.figure(0)
plt.plot(seconds,head)
plt.xlabel('Time [s]')
plt.ylabel('H [cm $H_2O$]')
plt.title('H vs t')
plt.grid()
plt.show()
fig0.savefig(filename[:-4] + '_fig0.pdf', bbox_inches='tight')


#Localise part of curve that represents slug test
tempHead = head[head.index(minHead):]


while k < len(tempHead):
    if tempHead[k] < baseLvl:
        endTest = 0
        k = k+1

        if k == len(tempHead):
            endTest = k
    else:
        endTest = k
        k = len(tempHead)

print 'End of slugtest is located at time t =', seconds[endTest], 'seconds'



#Draw slug test curve 
slugTestHeadNorm = []  
slugTestHead = head[head.index(minHead):head.index(minHead)+endTest]
slugTestHead[:] = [baseLvl - x for x in slugTestHead]
slugTestHeadNorm[:] = [(y/H_0) for y in slugTestHead]
slugTestSeconds = seconds[0:endTest]
slugTestSeconds[:] = [x-slugTestSeconds[0] for x in slugTestSeconds]



fig1 = plt.figure(1)
plt.plot(slugTestSeconds,slugTestHeadNorm)
plt.xlabel('Time [s]')
plt.ylabel('Normalized Head')
plt.title('Normalized Head vs t')
plt.grid(True)
plt.show() 
fig1.savefig(filename[:-4] + '_fig1.pdf', bbox_inches='tight')


slugTestSecondsShift = [] 
slugTestSecondsShift[:] = [x+1 for x in slugTestSeconds]
    

fig2=plt.figure(2)
plt.plot(slugTestSecondsShift, slugTestHeadNorm)
plt.xscale('log')
plt.xlabel('Time [s]')
plt.ylabel('Normalized Head')
plt.title('Normalized Head vs t')
plt.grid(True)
fig2.savefig(filename[:-4] + '_fig2.pdf', bbox_inches='tight')


lnSlugTestHead=[]

#taking ln of NORMALIZED head ln(H/H0)
while l < len(slugTestHead):
    lnSlugTestHead.append(math.log(slugTestHeadNorm[l]))
    l = l+1

#Gives time at which normalized head is equal to 0.035
#Is used to get only the nice looking first part of the curve in fig3
while o < len(slugTestHeadNorm):
    if slugTestHeadNorm[o] > 0.035: 
        x_cutoff = 0
        o = o+1
    else:
        x_cutoff = o-1
        o = len(slugTestHeadNorm)
    
#plots Hvorslev semilog plot with normalized head vs t for all data    
fig3 = plt.figure(3)
plt.plot(slugTestSeconds,slugTestHeadNorm, label='Normalized Head')
plt.ylim([0.05,1])
plt.xlim([0,slugTestSeconds[x_cutoff]]) 
plt.yscale('log')
fig3.axes[0].get_yaxis().set_major_formatter(FormatStrFormatter('%.1f'))
plt.xlabel('Time [s]')
plt.ylabel('Normalized Head')
plt.grid(True) 


#linear regression of all data, plotted on fig2
fit_all = linearRegression(lnSlugTestHead, slugTestSeconds)

x_all = fit_all['x']
y_all = fit_all['y']

r2_all = fit_all['r2']
regForm_all = fit_all['regForm']
a_all = fit_all['a']
b_all = fit_all['b']

plt.plot(x_all,y_all, '--k', label='Regression (All Data)')


tempHeadNorm = list(slugTestHeadNorm) #copies list
lenList = len(tempHeadNorm)

#Linear regression of data between 0.15 and 0.25, plotted on fig2
slice_Hvorslev_part = sliceList(tempHeadNorm, slugTestSeconds, 0.25, 0.15)
partHead_Hvorslev = slice_Hvorslev_part['partHead']
partLnHead_Hvorslev = slice_Hvorslev_part['partLnHead']
partSeconds_Hvorslev = slice_Hvorslev_part['partSeconds']

fit_Hvorslev = linearRegression(partLnHead_Hvorslev, partSeconds_Hvorslev)
x_Hvorslev = fit_Hvorslev['x']
y_Hvorslev = fit_Hvorslev['y']
r2_Hvorslev = fit_Hvorslev['r2']
regForm_Hvorslev = fit_Hvorslev['regForm']
a_Hvorslev = fit_Hvorslev['a']
b_Hvorslev = fit_Hvorslev['b']

plt.plot(x_all,y_Hvorslev, '--g', label='Regression (0.15-0.25)')


#Linear regression of data between 0.20 and 0.30, plotted on fig2
slice_BR_part = sliceList(tempHeadNorm, slugTestSeconds, 0.30, 0.20)
partHead_BR = slice_BR_part['partHead']
partLnHead_BR = slice_BR_part['partLnHead']
partSeconds_BR = slice_BR_part['partSeconds']

fit_BR = linearRegression(partLnHead_BR, partSeconds_BR)
x_BR = fit_BR['x']
y_BR = fit_BR['y']
r2_BR = fit_BR['r2']
regForm_BR = fit_BR['regForm']
a_BR = fit_BR['a']
b_BR = fit_BR['b']

plt.plot(x_all,y_BR, '--r', label='Regression (0.20-0.30)')




T_0 = math.exp(a_all*(0.368*H_0)+b_all)
#Linear regression of data between 1 and head correspond with T_0*1.5, plotted on fig2
#first = tempHeadNorm[int(round(T_0*2*1.5))]
slice_first = sliceList(tempHeadNorm, slugTestSeconds, 1.01, 0.3)
headFirst = slice_first['partHead']
LnHeadFirst = slice_first['partLnHead']
secondsFirst = slice_first['partSeconds']

fit_first = linearRegression(LnHeadFirst, secondsFirst)
x_first = fit_first['x']
y_first = fit_first['y']
r2_first = fit_first['r2']
regForm_first = fit_first['regForm']
a_first = fit_first['a']
b_first = fit_first['b']

plt.plot(x_all,y_first, '--m', label='Regression (0.3-1)')
legend = plt.legend(loc='upper right')
plt.show()


fig3.savefig(filename[:-4] + '_fig3.pdf', bbox_inches='tight')


#Calculate H_0+ interval 0.25-0.15
H_0plusNorm = math.exp(b_Hvorslev)*H_0
H_0plus = H_0plusNorm

#Calculate H_0+ interval 0.30-0.20
H_0plusNorm2 = math.exp(b_BR)*H_0
H_0plus2 = H_0plusNorm2

#Calculate H_0+ interval 1-H(1.5*T_0)
H_0plusNorm3 = math.exp(b_first)*H_0
H_0plus3 = H_0plusNorm3

#Calculates T_0 (time it takes for water lvl to rise/fall to 37% of H0)
T_0 = (-1-b_all)/a_all
T_0plus = (-1-b_Hvorslev)/a_Hvorslev
T_0plus2 = (-1-b_BR)/a_BR
T_0plus3 = (-1-b_first)/a_first

#Calculates Horizontal conductivity with Hvorslev, aniso = 1
k_h_H = calculateKhHvorslev(rc, rw, Le, Re, T_0, T_0plus, T_0plus3, aniso)
Kh_cms_full = k_h_H['Kh_cms_full'] #Kh entire head data set
Kh_md_full = k_h_H['Kh_md_full']
Kh_cms_part = k_h_H['Kh_cms_part'] #Kh interval 0.15-0.25
Kh_md_part = k_h_H['Kh_md_part']
Kh_cms_first = k_h_H['Kh_cms_first'] #Kh interval first part curve
Kh_md_first = k_h_H['Kh_md_first']
diff = 100-(Kh_cms_full/Kh_cms_part)*100

#Calculates Horizontal conductivity with Hvorslev, aniso = 0.1
k_h_H2 = calculateKhHvorslev(rc, rw, Le, Re, T_0, T_0plus, T_0plus3, aniso2)
Kh_cms_full2 = k_h_H2['Kh_cms_full'] #Kh entire head data set
Kh_md_full2 = k_h_H2['Kh_md_full']
Kh_cms_part2 = k_h_H2['Kh_cms_part'] #Kh interval 0.15-0.25
Kh_md_part2 = k_h_H2['Kh_md_part']
Kh_cms_first2 = k_h_H2['Kh_cms_first'] #Kh interval first part curve
Kh_md_first2 = k_h_H2['Kh_md_first']
diff2 = 100-(Kh_cms_full2/Kh_cms_part2)*100


print 'HVORSLEV'
print '*********'
print 'For anisotropy ratio (Kz/Kh) = 1:'
print 'Horizontal hydraulic conductivity based on all data: Kh = ', "%.8f"%Kh_cms_full, 'cm/s or ', "%.5f"%Kh_md_full, 'm/d'
print 'Horizontal hydraulic conductivity based on data in interval [0.15-0.25]: Kh = ', "%.8f"%Kh_cms_part, 'cm/s or ', "%.5f"%Kh_md_part, 'm/d'
print 'Horizontal hydraulic conductivity based on data in interval [0.3-1]: Kh = ', "%.8f"%Kh_cms_first, 'cm/s or ', "%.5f"%Kh_md_first, 'm/d'
print ' '
print 'For anisotropy ratio (Kz/Kh) = 0.1:'
print 'Horizontal hydraulic conductivity based on all data: Kh = ', "%.8f"%Kh_cms_full2, 'cm/s or ', "%.5f"%Kh_md_full2, 'm/d'
print 'Horizontal hydraulic conductivity based on data in interval [0.15-0.25]: Kh = ', "%.8f"%Kh_cms_part2, 'cm/s or ', "%.5f"%Kh_md_part2, 'm/d'
print 'Horizontal hydraulic conductivity based on data in interval [0.3-1]: Kh = ', "%.8f"%Kh_cms_first2, 'cm/s or ', "%.5f"%Kh_md_first2, 'm/d'
print ' '

#calculates Kh with Bouwer & Rice for aniso ratio 1 and 0.1
k_h_BR1 = calculateKhBouwerRice(rc, rw, Le, T_0, T_0plus2, T_0plus3, d, Bt, aniso)
Kh_BR1_cms = k_h_BR1['Kh_cms']
Kh_BR1_md = k_h_BR1['Kh_md']
Kh_BR1_cms_plus = k_h_BR1['Kh_cms_plus'] #for interval 0.2-0.3
Kh_BR1_md_plus = k_h_BR1['Kh_md_plus']
Kh_BR1_cms_first = k_h_BR1['Kh_cms_first'] #for first part of curve
Kh_BR1_md_first = k_h_BR1['Kh_md_first']


r_asterisk = k_h_BR1['r_asterisk']
A = k_h_BR1['A']
B = k_h_BR1['B']
coeff = k_h_BR1['coeff']
coeff2 = k_h_BR1['coeff2']


k_h_BR2 = calculateKhBouwerRice(rc, rw, Le, T_0, T_0plus2, T_0plus3, d, Bt, aniso2)
Kh_BR2_cms = k_h_BR2['Kh_cms']
Kh_BR2_md = k_h_BR2['Kh_md']
Kh_BR2_cms_plus = k_h_BR2['Kh_cms_plus'] #for interval 0.2-0.3
Kh_BR2_md_plus = k_h_BR2['Kh_md_plus']
Kh_BR2_cms_first = k_h_BR2['Kh_cms_first'] #for first part of curve
Kh_BR2_md_first = k_h_BR2['Kh_md_first']

r_asterisk2 = k_h_BR2['r_asterisk']
A2 = k_h_BR2['A']
B2 = k_h_BR2['B']
coeff_2 = k_h_BR2['coeff']
coeff2_2 = k_h_BR2['coeff2']


print 'BOUWER & RICE'
print '**************'
print 'Horizontal hydraulic conductivity based on all data and Kz/Kh = 1: Kh = ', "%.8f"%Kh_BR1_cms, 'cm/s or ', "%.5f"%Kh_BR1_md, 'm/d'
print 'Horizontal hydraulic conductivity based on all data and Kz/Kh = 0.1: Kh = ', "%.8f"%Kh_BR2_cms, 'cm/s or ', "%.5f"%Kh_BR2_md, 'm/d'
print ' '
print 'Horizontal hydraulic conductivity based on data in interval [0.20-0.30] and Kz/Kh = 1: Kh = ', "%.8f"%Kh_BR1_cms_plus, 'cm/s or ', "%.5f"%Kh_BR1_md_plus, 'm/d'
print 'Horizontal hydraulic conductivity based on data in interval [0.20-0.30] and Kz/Kh = 0.1: Kh = ', "%.8f"%Kh_BR2_cms_plus, 'cm/s or ', "%.5f"%Kh_BR2_md_plus, 'm/d'
print ' '
print 'Horizontal hydraulic conductivity based on data in interval [0.3-1] and Kz/Kh = 1: Kh = ', "%.8f"%Kh_BR1_cms_first, 'cm/s or ', "%.5f"%Kh_BR1_md_first, 'm/d'
print 'Horizontal hydraulic conductivity based on data in interval [0.3-1] and Kz/Kh = 0.1: Kh = ', "%.8f"%Kh_BR2_cms_first, 'cm/s or ', "%.5f"%Kh_BR2_md_first, 'm/d'


#writing output to file

def writeLogfile():
    logText='\n'.join(['##############################################################', 'Logfile Slug Test Analysis Script' , '##############################################################', ' ', ' '.join(['Time:',time.strftime("%a, %d %b %Y %H:%M:%S ")]), 'Working directory: '+ os.getcwd()+'\n', 'Filename: '+filename+'\n'])
    logText = logText + '\n***Identifying Start of Slug Test***\n'+'Start at t = ' + '%.1f'%t_0 +' s, H = ' + '%.2f'%minHead + ' cm\n\n'
    logText = logText + 'Equilibrium Head = '+ '%.2f'%baseLvl + ' cm' + '\n' + 'Initial Displacement H_0 = '+ '%.2f'%H_0 + ' cm' + '\nDuration Slugtest: ' + '%.1f'%slugTestSeconds[endTest-1] + ' s\n\n'
    logText = logText + 'Figure 0 is a H vs t plot of the entire head data set.\n' + '***Saved as ' + filename[:-4] + '_fig0.pdf***\n\n'
    logText = logText + 'Figure 1 is a normalized head vs t plot of the slugtest head data.\n' + '***Saved as ' + filename[:-4] + '_fig1.pdf***\n\n'
    logText = logText + 'Figure 2 is a normalized head vs log t plot of the slugtest head data.\n' + '***Saved as ' + filename[:-4] + '_fig2.pdf***\n\n'
    logText = logText + 'Figure 3 is a log normalized head vs t plot.\n' + '***Saved as ' + filename[:-4] + '_fig3.pdf***\n\n'
    logText = logText + 'Four regression lines are plotted:\n\n' + '1) Regression based on all head data.\n' + 'Regression Equation: y =' + ' %.13f'%a_all + 'x +' + ' %.11f'%b_all + '\nR² = ' + '%.4f'%float(r2_all) + '\n\n'
    logText = logText + '2) Regression based on normalized head data in interval 0.15-0.25.\n' + 'Regression Equation: y =' + ' %.13f'%a_Hvorslev + 'x +' + ' %.11f'%b_Hvorslev + '\nR² = ' + '%.4f'%float(r2_Hvorslev) + '\n\n'
    logText = logText + '3) Regression based on normalized head data in interval 0.20-0.30.\n' + 'Regression Equation: y =' + ' %.13f'%a_BR + 'x +' + ' %.11f'%b_BR + '\nR² = ' + '%.4f'%float(r2_BR) + '\n\n'
    logText = logText + '4) Regression based on normalized head data in interval 0.3-1.\n' + 'Regression Equation: y =' + ' %.13f'%a_first + 'x +' + ' %.11f'%b_first + '\nR² = ' + '%.4f'%float(r2_first) + '\n\n'
    logText = logText + '\n################'+'\n### HVORSLEV ###\n'+'################\n\n' + 'Parameters used:\n\n' + 'rc = ' + '%.3f'%rc + '\nrw = ' + '%.3f'%rw  + 'cm\nLe = ' + str(Le)  + 'cm\nRe = ' + '%.3f'%Re + ' cm\n\n'
    logText = logText + '1) Calculation for anisotropy ratio (Kz/Kh) of 1:\n'+'-------------------------------------------------\n\n'
    logText = logText + 'a) Calculation based on all head data.\n\n' + 'H_0 = '+ '%.2f'%H_0 + ' cm\n' + 'T_0 = ' + '%.1f'%T_0 + ' s\n\n' + 'Horizontal hydraulic conductivity = ' + '%.5f'%Kh_cms_full + ' cm/s' + ' or ' + '%.3f'%Kh_md_full + ' m/d\n\n'
    logText = logText + 'b) Calculation based on normalized head data in interval 0.15-0.25.\n\n' + 'H_0+ = '+ '%.2f'%H_0plus + ' cm\n' + 'T_0+ = ' + '%.1f'%T_0plus + ' s\n\n' + 'Horizontal hydraulic conductivity = ' + '%.5f'%Kh_cms_part + ' cm/s' + ' or ' + '%.3f'%Kh_md_part + ' m/d\n\n'
    logText = logText + 'c) Calculation based on normalized head data in interval 0.3-1.\n\n' + 'H_0+ = '+ '%.2f'%H_0plus3 + ' cm\n' + 'T_0+ = ' + '%.1f'%T_0plus3 + ' s\n\n' + 'Horizontal hydraulic conductivity = ' + '%.5f'%Kh_cms_first + ' cm/s' + ' or ' + '%.3f'%Kh_md_first + ' m/d\n\n\n'
    logText = logText + '2) Calculation for anisotropy ratio (Kz/Kh) of 0.1:\n'+'---------------------------------------------------\n\n'
    logText = logText + 'a) Calculation based on all head data.\n\n' + 'H_0 = '+ '%.2f'%H_0 + ' cm\n' + 'T_0 = ' + '%.1f'%T_0 + ' s\n\n' + 'Horizontal hydraulic conductivity = ' + '%.5f'%Kh_cms_full2 + ' cm/s' + ' or ' + '%.3f'%Kh_md_full2 + ' m/d\n\n'
    logText = logText + 'b) Calculation based on normalized head data in interval 0.15-0.25.\n\n' + 'H_0+ = '+ '%.2f'%H_0plus + ' cm\n' + 'T_0+ = ' + '%.1f'%T_0plus + ' s\n\n' + 'Horizontal hydraulic conductivity = ' + '%.5f'%Kh_cms_part2 + ' cm/s' + ' or ' + '%.3f'%Kh_md_part2 + ' m/d\n\n'
    logText = logText + 'c) Calculation based on normalized head data in interval 0.3-1.\n\n' + 'H_0+ = '+ '%.2f'%H_0plus3 + ' cm\n' + 'T_0+ = ' + '%.1f'%T_0plus3 + ' s\n\n' + 'Horizontal hydraulic conductivity = ' + '%.5f'%Kh_cms_first2 + ' cm/s' + ' or ' + '%.3f'%Kh_md_first2 + ' m/d\n\n\n\n'   
    logText = logText + '\n#####################\n'+'### BOUWER & RICE ###\n'+'#####################\n\n' + 'Parameters used:\n\n' + 'rc = ' + '%.3f'%rc + ' cm\nLe = ' + str(Le) + ' cm\nRe = ' + '%.3f'%Re + ' cm\n' + 'rw = ' + '%.3f'%rw + ' cm\n' + 'd = ' + str(d) + ' cm\nBt = ' + str(Bt) + ' cm\n\n'
    logText = logText + '1) Calculation for anisotropy ratio (Kz/Kh) of 1:\n'+'-------------------------------------------------\n\n' + 'rw* = ' + '%.3f'%r_asterisk + ' cm\nA = ' + '%.11f'%A + '\nB = ' + '%.12f'%B + '\nln[(B-(d+b))/rw*] = ' + '%.2f'%coeff + '\nln[Re/rw*] = ' + '%.11f'%coeff2 + '\n\n' + 'a) Calculation based on all head data.\n\n'+ 'H_0 = '+ '%.2f'%H_0 + ' cm\n' + 'T_0 = ' + '%.1f'%T_0 +' s\n\n'  
    logText = logText + 'Horizontal hydraulic conductivity = ' + '%.5f'%Kh_BR1_cms + ' cm/s' + ' or ' + '%.3f'%Kh_BR1_md + ' m/d\n\n'
    logText = logText + 'b) Calculation based on normalized head data in interval 0.20-0.30.\n\n'+ 'H_0 = '+ '%.2f'%H_0plus2 + ' cm\n' + 'T_0 = ' + '%.1f'%T_0plus2 + ' s\n\n' + 'Horizontal hydraulic conductivity = ' + '%.5f'%Kh_BR1_cms_plus + ' cm/s' + ' or ' + '%.3f'%Kh_BR1_md_plus + ' m/d\n\n'
    logText = logText + 'c) Calculation based on normalized head data in interval 0.3-1.\n\n'+ 'H_0 = '+ '%.2f'%H_0plus3 + ' cm\n' + 'T_0 = ' + '%.1f'%T_0plus3 + ' s\n\n' + 'Horizontal hydraulic conductivity = ' + '%.5f'%Kh_BR1_cms_first + ' cm/s' + ' or ' + '%.3f'%Kh_BR1_md_first + ' m/d\n\n\n'
    logText = logText + '2) Calculation for anisotropy ratio (Kz/Kh) of 0.1:\n'+'---------------------------------------------------\n\n' + 'rw* = ' + '%.3f'%r_asterisk2 + ' cm\nA = ' + '%.11f'%A2 + '\nB = ' + '%.12f'%B2 + '\nln[(B-(d+b))/rw*] = ' + '%.2f'%coeff_2 + '\nln[Re/rw*] = ' + '%.11f'%coeff2_2 + '\n\n'
    logText = logText + 'a) Calculation based on all head data.\n\n'+ 'H_0 = '+ '%.2f'%H_0 + ' cm\n' + 'T_0 = ' + '%.1f'%T_0 +' s\n\n'
    logText = logText + 'Horizontal hydraulic conductivity = ' + '%.5f'%Kh_BR2_cms + ' cm/s' + ' or ' + '%.3f'%Kh_BR2_md + ' m/d\n\n'
    logText = logText + 'b) Calculation based on normalized head data in interval 0.20-0.30.\n\n'+ 'H_0 = '+ '%.2f'%H_0plus2 + ' cm\n' + 'T_0 = ' + '%.1f'%T_0plus2 + ' s\n\n' + 'Horizontal hydraulic conductivity = ' + '%.5f'%Kh_BR2_cms_plus + ' cm/s' + ' or ' + '%.3f'%Kh_BR2_md_plus + ' m/d\n\n'
    logText = logText + 'c) Calculation based on normalized head data in interval 0.3-1.\n\n'+ 'H_0 = '+ '%.2f'%H_0plus3 + ' cm\n' + 'T_0 = ' + '%.1f'%T_0plus3 + ' s\n\n' + 'Horizontal hydraulic conductivity = ' + '%.5f'%Kh_BR2_cms_first + ' cm/s' + ' or ' + '%.3f'%Kh_BR2_md_first + ' m/d\n\n\n'
    logText = logText + 'List of conductivity values:\n----------------------------\n' + '%.2f'%Kh_md_full+'\n' + '%.2f'%Kh_md_part+'\n' + '%.2f'%Kh_md_first+'\n' + '%.2f'%Kh_md_full2+'\n' + '%.2f'%Kh_md_part2+'\n' + '%.2f'%Kh_md_first2+'\n' + '%.2f'%Kh_BR1_md+'\n' + '%.2f'%Kh_BR1_md_plus+'\n' + '%.2f'%Kh_BR1_md_first+'\n' + '%.2f'%Kh_BR2_md+'\n' + '%.2f'%Kh_BR2_md_plus+'\n' + '%.2f'%Kh_BR2_md_first+'\n\n\n'    
    logText = logText + '                _____ _          _____       _ \n'+'               |_   _| |_ ___   |   __|___ _| |\n' + '                 | | |   | -_|  |   __|   | . |\n' + '                 |_| |_|_|___|  |_____|_|_|___|\n'
    
    return logText

logText = writeLogfile()

with open(filename[:-4] + '_logfile' + '.dat', 'w') as output:    
    
    output.write(logText)
    