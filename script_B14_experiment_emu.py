#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 10:28:13 2021

@author: dla


This contains the calcs of the B.1.4. experiment with the spectras. 

*Crystal: BGO
*Sample elevatipn: 1 cylinder and 0 for Na and Co
*Similar coniguraiotn for all the measurements (obviously)

"""

#reset to manually clear all the variables
#clear               #to clear the command windows
#%reset -f          #to clear all the variables without confirmation
#magic('reset -sf')

#######General packages useful##

import matplotlib.pyplot as plt  #for simplicity, to not write matplotlib.pyplot
        #everytime we want to plot something

#from scipy.stats import norm               ##norm.fit() fit to gaussian
import numpy as np
    #np contain linspaces as np.linspace(a,b,N)
import pandas as pd      
#from plotly.graph_objs import Bar, Layout
#from plotly import offline

import sys                   #to import functions from other folders!!
sys.path.insert(0, '/home/dla/Python/Functions_homemade')   #path where I have the functions

import Fits
######3



#%% ###############################################################
#########################0), Data loading #####################3
################################################################

#The files to load are in txt. The best way to read is:

#All with same histogram, gain, bias, threshold, but different gate to optimize 
#each signal

#Firstly create the data by setting the led driver amplitude to 0, obtaining 
#extrange results, so will now unplug it, and closing the SiPM to see if now 
#the results are fine.

#indexes
Cs = 0
Co =1
Na = 2
Ba = 3

time = np.zeros(4)                          #[s] duration of the measurements
time[Cs], time[Co], time[Na], time[Ba] = 120,120,120,120

counts_stored = np.array([])                          #storing variable of the counts
ADC_channel = np.array([])                          #to store the channels

######## Cs137 #########
Cs137 = pd.read_csv('Cs137_120s_gain_0_3.csv', sep = ',', header = None, 
                   names = ['Counts','None!'])
                #, as sepparator, no header. I rename the columns. THere is an
                #empty column!
                            
#Storing of the values
counts_stored = np.append(counts_stored,Cs137['Counts'])
ADC_channel = np.append(ADC_channel, np.arange(0, len(counts_stored)) ) 
        #all have the same channels!



########### Co60 ############     
Co60 = pd.read_csv('Co_60_120s_gain_0_3.csv', sep = ',', header = None, 
                   names = ['Counts','None!'])   

#Storing of the values (the 2nd storing and more has to be columns stack!!)
counts_stored = np.column_stack( (counts_stored,Co60['Counts']) )

###  


############### Na22 ###################
Na22 = pd.read_csv('Na22_120s_gain_0_3.csv', sep = ',', header = None, 
                   names = ['Counts','None!'])      
counts_stored = np.column_stack( (counts_stored,Na22['Counts']) )



########### Ba 133 ####################
Ba133 = pd.read_csv('Ba133_120s_gain_0_3.csv', sep = ',', header = None, 
                   names = ['Counts','None!'])         
counts_stored = np.column_stack( (counts_stored,Ba133['Counts']) )
 



###Plot combined, as CAEN's
#plotting the count rate, since the measure time is differen!!
plt.figure(figsize=(8,5))  #width, heigh 6.4*4.8 inches by default
plt.plot(ADC_channel, counts_stored[:,Cs] )    
plt.plot(ADC_channel, counts_stored[:,Co])
plt.plot(ADC_channel, counts_stored[:,Na])  
plt.plot(ADC_channel, counts_stored[:,Ba])       
plt.legend(['Cs', 'Co', 'Na', 'Ba'], fontsize=10) 
plt.title("Spectras using BGO (emulator)", fontsize=20)           #title
plt.xlabel("ADC channels", fontsize=10)                        #xlabel
plt.ylabel("Counts", fontsize=10)              #ylabel
# Set size of tick labels.
plt.tick_params(axis='both', labelsize=10)              #size of axis
plt.grid(True) 
plt.xlim(0,3000)                   #limits of x axis
plt.savefig('counts_emu.png', format='png')


#plt.ylim(0,11000)                            #limits of y axis



#%% #################################################################
############################1) FIT #################################
###########################################################

#Lets do the gaussian fit to the gamma () of Cs137, to see the FWHM as a function
#of the scintillation crystal

def gaussian(x, Heigh, Mean, Std_dev):
    return Heigh * np.exp(- (x-Mean)**2 / (2 * Std_dev**2)) 
            #this is a gaussian function (more general than normal distribution)

#$$$$$$$$$$$$$$$$$$$ Cs $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#To find the intervals, I could easily find the index of the maximum, and from
#there move by looking at the .txt, or maybe say, peak me 100 hundred values above
#and below the max, etc.


#I will have to serach it 
#by hand by looking at the .txt file #and the plot:
    
index_min = 1200          #index that starts the peak (by looking the graph)
index_max = 1500          #index that ends the peak (by looking the graph)

counts_peak = counts_stored[index_min-1:index_max-1,Cs] 
ch_peak = ADC_channel[index_min-1:index_max-1]

counts_max = max(counts_peak)
counts_max_index = np.where( counts_peak == counts_max )
ch_counts_max = ch_peak[counts_max_index[0][0] ]
#counts_peak[ counts_max_index[0][0] ] == counts_max            #debug


###Plot of the interval of the peak (indexes)
plt.figure(figsize=(10,6))                  #width, heigh 6.4*4.8 inches by default
plt.plot( ADC_channel, counts_stored[:,Cs], 'b.-')    
plt.plot( ch_peak , counts_peak, 'r.-')   
plt.title("Cs137", fontsize=22)           #title
plt.xlabel("ADC channels", fontsize=14)                        #xlabel
plt.ylabel("Counts", fontsize=14)              #ylabel
# Set size of tick labels.
plt.tick_params(axis='both', labelsize=14)              #size of axis
plt.grid(True) 
plt.xlim(0,3000)                   #limits of x axis
plt.savefig('Signal_Cs.png', format='png')
###

fit = Fits.Gaussian_fit(ch_peak,counts_peak)


#Storing of the relevant data, sigma and its error
sigma_stored = np.array([])
mean_stored = np.array([])
delta_mean_stored = np.array([])
delta_sigma_stored = np.array([])
FWHM_stored = np.array([])
delta_FWHM_stored = np.array([])
heigh_stored = np.array([])
delta_heigh_stored = np.array([])
ch_peak_stored = np.array([])


sigma_stored = np.append(sigma_stored, fit['sigma'])
mean_stored = np.append(mean_stored, fit['mean'])
delta_mean_stored = np.append(delta_mean_stored, fit['\Delta(mean)'])
delta_sigma_stored = np.append(delta_sigma_stored, fit['\Delta(sigma)'])
FWHM_stored = np.append(FWHM_stored, fit['FWHM'])
delta_FWHM_stored = np.append(delta_FWHM_stored, fit['\Delta(FWHM)'])
heigh_stored = np.append(heigh_stored, fit['heigh'])
delta_heigh_stored = np.append(delta_heigh_stored, fit['\Delta(heigh)'])
ch_peak_stored = np.append(ch_peak_stored, ch_counts_max)

############

#Plot of FWHM and Centroid for a peak. i will plot lines showing the channels for the CsI, the best.

point_left_mean = mean_stored[Cs] - FWHM_stored[Cs] / 2        #Point to the left of the mean value with max ampl/2
point_right_mean = mean_stored[Cs] + FWHM_stored[Cs] / 2       #Point to the right of the mean value with max ampl/2

plt.figure(figsize=(10,6))  #width, heigh 6.4*4.8 inches by default
#plt.plot( ADC_channel[:,0], counts_stored[:,0], 'b.-')    
#plt.plot( ch_peak , counts_peak, 'b.-')   
plt.bar(ch_peak, counts_peak, width = ADC_channel[1]-ADC_channel[0])        #peak
plt.plot(ch_peak, gaussian(ch_peak, heigh_stored[Cs], mean_stored[Cs], sigma_stored[Cs]), 'ro')    #fit

plt.plot(mean_stored[Cs] * np.array([1, 1]), [min(counts_stored[:,Cs]), max(counts_stored[:,Cs]) ], 'm--', linewidth=3 )
                            #Vertical line for the centroid

plt.plot(point_left_mean * np.array([1, 1]) , [min(counts_stored[:,Cs]), max(counts_stored[:,Cs]) ], 'k--', linewidth=3 )
plt.plot(point_right_mean * np.array([1, 1]) , [min(counts_stored[:,Cs]), max(counts_stored[:,Cs]) ], 'k--', linewidth=3 )

plt.title("Cs137 photopeak", fontsize=22)           #title
plt.xlabel("ADC channels", fontsize=14)                        #xlabel
plt.ylabel("Counts", fontsize=14)              #ylabel
# Set size of tick labels.
plt.legend(['fit', 'Centroid', 'FWHM'], fontsize=10) 
plt.tick_params(axis='both', labelsize=14)              #size of axis
plt.grid(True) 
plt.savefig('Signal_Cs_R_visual.png', format='png')
###



#$$$$$$$$$$$$$$$$$$$ Co $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#Here we have 2 peaks, so finding the channels will be more tough


######Peak  1#####
index_min = 2050#          #index that starts the peak (by looking the graph)
index_max = 2400          #index that ends the peak (by looking the graph)
                #830, 1150 for the 26_5 measurements (9000 ch)
                #1100, 1400 for the 31/5 measurements
                
counts_peak = counts_stored[index_min-1:index_max-1,Co] 
ch_peak = ADC_channel[index_min-1:index_max-1]

counts_max = max(counts_peak)
counts_max_index = np.where( counts_peak == counts_max )
ch_counts_max = ch_peak[counts_max_index[0][0] ]





###Plot of the interval of the peak (indexes)
plt.figure(figsize=(10,6))                  #width, heigh 6.4*4.8 inches by default
plt.plot( ADC_channel, counts_stored[:,Co], 'b.-')      
plt.plot( ch_peak , counts_peak, 'r.-')   
plt.title("Co60", fontsize=22)           #title
plt.xlabel("ADC channels", fontsize=14)                        #xlabel
plt.ylabel("Counts", fontsize=14)              #ylabel
# Set size of tick labels.
plt.tick_params(axis='both', labelsize=14)              #size of axis
plt.grid(True) 
plt.xlim(0,3000)                   #limits of x axis
plt.savefig('Signal_Co60_pico1_sum.png', format='png')
###


fit = Fits.Gaussian_fit(ch_peak,counts_peak)


sigma_stored = np.append(sigma_stored, fit['sigma'])
mean_stored = np.append(mean_stored, fit['mean'])
delta_mean_stored = np.append(delta_mean_stored, fit['\Delta(mean)'])
delta_sigma_stored = np.append(delta_sigma_stored, fit['\Delta(sigma)'])
FWHM_stored = np.append(FWHM_stored, fit['FWHM'])
delta_FWHM_stored = np.append(delta_FWHM_stored, fit['\Delta(FWHM)'])
heigh_stored = np.append(heigh_stored, fit['heigh'])
delta_heigh_stored = np.append(delta_heigh_stored, fit['\Delta(heigh)'])
ch_peak_stored = np.append(ch_peak_stored, ch_counts_max)

#Plot of FWHM and Centroid for a peak. i will plot lines showing the channels for the CsI, the best.

point_left_mean = mean_stored[Co] - FWHM_stored[Co] / 2        #Point to the left of the mean value with max ampl/2
point_right_mean = mean_stored[Co] + FWHM_stored[Co] / 2       #Point to the right of the mean value with max ampl/2

#
plt.figure(figsize=(10,6))  #width, heigh 6.4*4.8 inches by default
#plt.plot( ADC_channel[:,0], counts_stored[:,0], 'b.-')    
#plt.plot( ch_peak , counts_peak, 'b.-')   
plt.bar(ch_peak, counts_peak, width = ADC_channel[1]-ADC_channel[0])        #peak
plt.plot(ch_peak, gaussian(ch_peak, heigh_stored[Co], mean_stored[Co], sigma_stored[Co]), 'ro')    #fit
plt.plot(mean_stored[Co] * np.array([1, 1]), [min(counts_stored[:,Co]), max(counts_stored[:,Co]) ], 'm--', linewidth=3 )
                            #Vertical line for the centroid
plt.plot(point_left_mean * np.array([1, 1]) , [min(counts_stored[:,Co]), max(counts_stored[:,Co]) ], 'k--', linewidth=3 )
plt.plot(point_right_mean * np.array([1, 1]) , [min(counts_stored[:,Co]), max(counts_stored[:,Co]) ], 'k--', linewidth=3 )
plt.title("Co photo", fontsize=22)           #title
plt.xlabel("ADC channels", fontsize=14)                        #xlabel
plt.ylabel("Counts", fontsize=14)              #ylabel
# Set size of tick labels.
plt.legend(['fit', 'Centroid', 'FWHM'], fontsize=10) 
plt.tick_params(axis='both', labelsize=14)              #size of axis
plt.grid(True) 
plt.savefig('Signal_Co_R_pico1_visual.png', format='png')


######Peak  2#####
index_min = 2400#          #index that starts the peak (by looking the graph)
index_max = 2700          #index that ends the peak (by looking the graph)
                #830, 1150 for the 26_5 measurements (9000 ch)
                #1100, 1400 for the 31/5 measurements
        
counts_peak = counts_stored[index_min-1:index_max-1,Co] 
ch_peak = ADC_channel[index_min-1:index_max-1]

counts_max = max(counts_peak)
counts_max_index = np.where( counts_peak == counts_max )
ch_counts_max = ch_peak[counts_max_index[0][0] ]


#Plot
plt.figure(figsize=(10,6))                  #width, heigh 6.4*4.8 inches by default
plt.plot( ADC_channel, counts_stored[:,Co], 'b.-')      
plt.plot( ch_peak , counts_peak, 'r.-')   
plt.title("Co60", fontsize=22)           #title
plt.xlabel("ADC channels", fontsize=14)                        #xlabel
plt.ylabel("Counts", fontsize=14)              #ylabel
# Set size of tick labels.
plt.tick_params(axis='both', labelsize=14)              #size of axis
plt.grid(True) 
plt.xlim(0,3000)                   #limits of x axis
plt.savefig('Signal_Co60_pico2_sum.png', format='png')
###


fit = Fits.Gaussian_fit(ch_peak,counts_peak)


sigma_stored = np.append(sigma_stored, fit['sigma'])
mean_stored = np.append(mean_stored, fit['mean'])
delta_mean_stored = np.append(delta_mean_stored, fit['\Delta(mean)'])
delta_sigma_stored = np.append(delta_sigma_stored, fit['\Delta(sigma)'])
FWHM_stored = np.append(FWHM_stored, fit['FWHM'])
delta_FWHM_stored = np.append(delta_FWHM_stored, fit['\Delta(FWHM)'])
heigh_stored = np.append(heigh_stored, fit['heigh'])
delta_heigh_stored = np.append(delta_heigh_stored, fit['\Delta(heigh)'])
ch_peak_stored = np.append(ch_peak_stored, ch_counts_max)

#Plot of FWHM and Centroid for a peak. i will plot lines showing the channels for the CsI, the best.

point_left_mean = mean_stored[Co+1] - FWHM_stored[Co+1] / 2        #Point to the left of the mean value with max ampl/2
point_right_mean = mean_stored[Co+1] + FWHM_stored[Co+1] / 2       #Point to the right of the mean value with max ampl/2
            #+1 because Na have 2 peaks, the 1st have index Na, and the second one, if wanting to
            #keep the notation, Na+1 (so this one will have to be added always!)

#
plt.figure(figsize=(10,6))  #width, heigh 6.4*4.8 inches by default
#plt.plot( ADC_channel[:,0], counts_stored[:,0], 'b.-')    
#plt.plot( ch_peak , counts_peak, 'b.-')   
plt.bar(ch_peak, counts_peak, width = ADC_channel[1]-ADC_channel[0])        #peak
plt.plot(ch_peak, gaussian(ch_peak, heigh_stored[Co+1], mean_stored[Co+1], sigma_stored[Co+1]), 'ro')    #fit
plt.plot(mean_stored[Co+1] * np.array([1, 1]), [min(counts_stored[:,Co]), max(counts_stored[:,Co]) ], 'm--', linewidth=3 )
                            #Vertical line for the centroid
plt.plot(point_left_mean * np.array([1, 1]) , [min(counts_stored[:,Co]), max(counts_stored[:,Co]) ], 'k--', linewidth=3 )
plt.plot(point_right_mean * np.array([1, 1]) , [min(counts_stored[:,Co]), max(counts_stored[:,Co]) ], 'k--', linewidth=3 )
plt.title("Co photo", fontsize=22)           #title
plt.xlabel("ADC channels", fontsize=14)                        #xlabel
plt.ylabel("Counts", fontsize=14)              #ylabel
# Set size of tick labels.
plt.legend(['fit', 'Centroid', 'FWHM'], fontsize=10) 
plt.tick_params(axis='both', labelsize=14)              #size of axis
plt.grid(True) 
plt.savefig('Signal_Co_R_pico2_visual.png', format='png')


#$$$$$$$$$$$$$$$$$$$ Na $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

#here, again, the maximum peak is not the one of the gamma peak, so again have
#to find the peak by hand :)

###Peak 1###

index_min = 850         #index that starts the peak (by looking the graph)
index_max = 1100          #index that ends the peak (by looking the graph)
                #2531, 3181 for the 26/5 measurements (9000ch)
                #2800 3500 for the 31/5 measurements (600s)

counts_peak = counts_stored[index_min-1:index_max-1,Na] 
ch_peak = ADC_channel[index_min-1:index_max-1]

counts_max = max(counts_peak)
counts_max_index = np.where( counts_peak == counts_max )
ch_counts_max = ch_peak[counts_max_index[0][0] ]


#Plot
plt.figure(figsize=(10,6))                  #width, heigh 6.4*4.8 inches by default
plt.plot( ADC_channel, counts_stored[:,Na], 'b.-')     
plt.plot( ch_peak , counts_peak, 'r.-')   
plt.title("Na22", fontsize=22)           #title
plt.xlabel("ADC channels", fontsize=14)                        #xlabel
plt.ylabel("Counts", fontsize=14)              #ylabel
# Set size of tick labels.
plt.tick_params(axis='both', labelsize=14)              #size of axis
plt.grid(True) 
plt.xlim(0,3000)                   #limits of x axis
plt.savefig('Signal_Na_pico1_sum.png', format='png')
###

fit = Fits.Gaussian_fit(ch_peak,counts_peak)

sigma_stored = np.append(sigma_stored, fit['sigma'])
mean_stored = np.append(mean_stored, fit['mean'])
delta_mean_stored = np.append(delta_mean_stored, fit['\Delta(mean)'])
delta_sigma_stored = np.append(delta_sigma_stored, fit['\Delta(sigma)'])
FWHM_stored = np.append(FWHM_stored, fit['FWHM'])
delta_FWHM_stored = np.append(delta_FWHM_stored, fit['\Delta(FWHM)'])
heigh_stored = np.append(heigh_stored, fit['heigh'])
delta_heigh_stored = np.append(delta_heigh_stored, fit['\Delta(heigh)'])
ch_peak_stored = np.append(ch_peak_stored, ch_counts_max)


#Plot of FWHM and Centroid for a peak. i will plot lines showing the channels for the CsI, the best.

point_left_mean = mean_stored[Na+1] - FWHM_stored[Na+1] / 2        #Point to the left of the mean value with max ampl/2
point_right_mean = mean_stored[Na+1] + FWHM_stored[Na+1] / 2       #Point to the right of the mean value with max ampl/2

#
plt.figure(figsize=(10,6))  #width, heigh 6.4*4.8 inches by default
#plt.plot( ADC_channel[:,0], counts_stored[:,0], 'b.-')    
#plt.plot( ch_peak , counts_peak, 'b.-')   
plt.bar(ch_peak, counts_peak, width = ADC_channel[1]-ADC_channel[0])        #peak
plt.plot(ch_peak, gaussian(ch_peak, heigh_stored[Na+1], mean_stored[Na+1], sigma_stored[Na+1]), 'ro')    #fit
plt.plot(mean_stored[Na+1] * np.array([1, 1]), [min(counts_stored[:,Na]), max(counts_stored[:,Na]) ], 'm--', linewidth=3 )
                            #Vertical line for the centroid
plt.plot(point_left_mean * np.array([1, 1]) , [min(counts_stored[:,Na]), max(counts_stored[:,Na]) ], 'k--', linewidth=3 )
plt.plot(point_right_mean * np.array([1, 1]) , [min(counts_stored[:,Na]), max(counts_stored[:,Na]) ], 'k--', linewidth=3 )
plt.title("Na22", fontsize=22)           #title
plt.xlabel("ADC channels", fontsize=14)                        #xlabel
plt.ylabel("Counts", fontsize=14)              #ylabel
# Set size of tick labels.
plt.legend(['fit', 'Centroid', 'FWHM'], fontsize=10) 
plt.tick_params(axis='both', labelsize=14)              #size of axis
plt.grid(True) 
plt.savefig('Signal_Na_R_pico1_visual.png', format='png')


###Peak 2###

index_min = 2320         #index that starts the peak (by looking the graph)
index_max = 2550          #index that ends the peak (by looking the graph)
                #2531, 3181 for the 26/5 measurements (9000ch)
                #2800 3500 for the 31/5 measurements (600s)

counts_peak = counts_stored[index_min-1:index_max-1,Na] 
ch_peak = ADC_channel[index_min-1:index_max-1]

counts_max = max(counts_peak)
counts_max_index = np.where( counts_peak == counts_max )
ch_counts_max = ch_peak[counts_max_index[0][0] ]


#Plot
plt.figure(figsize=(10,6))                  #width, heigh 6.4*4.8 inches by default
plt.plot( ADC_channel, counts_stored[:,2], 'b.-')     
plt.plot( ch_peak , counts_peak, 'r.-')   
plt.title("Na22", fontsize=22)           #title
plt.xlabel("ADC channels", fontsize=14)                        #xlabel
plt.ylabel("Counts", fontsize=14)              #ylabel
# Set size of tick labels.
plt.tick_params(axis='both', labelsize=14)              #size of axis
plt.grid(True) 
plt.xlim(0,3000)                   #limits of x axis
plt.savefig('Signal_Na_pico2_sum.png', format='png')
###

fit = Fits.Gaussian_fit(ch_peak,counts_peak)

sigma_stored = np.append(sigma_stored, fit['sigma'])
mean_stored = np.append(mean_stored, fit['mean'])
delta_mean_stored = np.append(delta_mean_stored, fit['\Delta(mean)'])
delta_sigma_stored = np.append(delta_sigma_stored, fit['\Delta(sigma)'])
FWHM_stored = np.append(FWHM_stored, fit['FWHM'])
delta_FWHM_stored = np.append(delta_FWHM_stored, fit['\Delta(FWHM)'])
heigh_stored = np.append(heigh_stored, fit['heigh'])
delta_heigh_stored = np.append(delta_heigh_stored, fit['\Delta(heigh)'])
ch_peak_stored = np.append(ch_peak_stored, ch_counts_max)

#Plot of FWHM and Centroid for a peak. i will plot lines showing the channels for the CsI, the best.

point_left_mean = mean_stored[Na+2] - FWHM_stored[Na+2] / 2        #Point to the left of the mean value with max ampl/2
point_right_mean = mean_stored[Na+2] + FWHM_stored[Na+2] / 2       #Point to the right of the mean value with max ampl/2

#
plt.figure(figsize=(10,6))  #width, heigh 6.4*4.8 inches by default
#plt.plot( ADC_channel[:,0], counts_stored[:,0], 'b.-')    
#plt.plot( ch_peak , counts_peak, 'b.-')   
plt.bar(ch_peak, counts_peak, width = ADC_channel[1]-ADC_channel[0])        #peak
plt.plot(ch_peak, gaussian(ch_peak, heigh_stored[Na+2], mean_stored[Na+2], sigma_stored[Na+2]), 'ro')    #fit
plt.plot(mean_stored[Na+2] * np.array([1, 1]), [min(counts_stored[:,Na]), max(counts_stored[:,Na]) ], 'm--', linewidth=3 )
                            #Vertical line for the centroid
plt.plot(point_left_mean * np.array([1, 1]) , [min(counts_stored[:,Na]), max(counts_stored[:,Na]) ], 'k--', linewidth=3 )
plt.plot(point_right_mean * np.array([1, 1]) , [min(counts_stored[:,Na]), max(counts_stored[:,Na]) ], 'k--', linewidth=3 )
plt.title("Na22", fontsize=22)           #title
plt.xlabel("ADC channels", fontsize=14)                        #xlabel
plt.ylabel("Counts", fontsize=14)              #ylabel
# Set size of tick labels.
plt.legend(['fit', 'Centroid', 'FWHM'], fontsize=10) 
plt.tick_params(axis='both', labelsize=14)              #size of axis
plt.grid(True) 
plt.savefig('Signal_Na_R_pico2_visual.png', format='png')



#$$$$$$$$$$$$$$$$$$$ Ba $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

#here, again, the maximum peak is not the one of the gamma peak, so again have
#to find the peak by hand :)

index_min = 630         #index that starts the peak (by looking the graph)
index_max = 750          #index that ends the peak (by looking the graph)
                #2531, 3181 for the 26/5 measurements (9000ch)
                #2800 3500 for the 31/5 measurements (600s)

counts_peak = counts_stored[index_min-1:index_max-1,Ba] 
ch_peak = ADC_channel[index_min-1:index_max-1]

counts_max = max(counts_peak)
counts_max_index = np.where( counts_peak == counts_max )
ch_counts_max = ch_peak[counts_max_index[0][0] ]


#Plot
plt.figure(figsize=(10,6))                  #width, heigh 6.4*4.8 inches by default
plt.plot( ADC_channel, counts_stored[:,Ba], 'b.-')     
plt.plot( ch_peak , counts_peak, 'r.-')   
plt.title("Ba133", fontsize=22)           #title
plt.xlabel("ADC channels", fontsize=14)                        #xlabel
plt.ylabel("Counts", fontsize=14)              #ylabel
# Set size of tick labels.
plt.tick_params(axis='both', labelsize=14)              #size of axis
plt.grid(True) 
plt.xlim(0,3000)                   #limits of x axis
plt.savefig('Signal_LYSO_sum.png', format='png')
###

fit = Fits.Gaussian_fit(ch_peak,counts_peak)

sigma_stored = np.append(sigma_stored, fit['sigma'])
mean_stored = np.append(mean_stored, fit['mean'])
delta_mean_stored = np.append(delta_mean_stored, fit['\Delta(mean)'])
delta_sigma_stored = np.append(delta_sigma_stored, fit['\Delta(sigma)'])
FWHM_stored = np.append(FWHM_stored, fit['FWHM'])
delta_FWHM_stored = np.append(delta_FWHM_stored, fit['\Delta(FWHM)'])
heigh_stored = np.append(heigh_stored, fit['heigh'])
delta_heigh_stored = np.append(delta_heigh_stored, fit['\Delta(heigh)'])
ch_peak_stored = np.append(ch_peak_stored, ch_counts_max)

#Plot of FWHM and Centroid for a peak. i will plot lines showing the channels for the CsI, the best.

point_left_mean = mean_stored[Ba+2] - FWHM_stored[Ba+2] / 2        #Point to the left of the mean value with max ampl/2
point_right_mean = mean_stored[Ba+2] + FWHM_stored[Ba+2] / 2       #Point to the right of the mean value with max ampl/2

#
plt.figure(figsize=(10,6))  #width, heigh 6.4*4.8 inches by default
#plt.plot( ADC_channel[:,0], counts_stored[:,0], 'b.-')    
#plt.plot( ch_peak , counts_peak, 'b.-')   
plt.bar(ch_peak, counts_peak, width = ADC_channel[1]-ADC_channel[0])        #peak
plt.plot(ch_peak, gaussian(ch_peak, heigh_stored[Ba+2], mean_stored[Ba+2], sigma_stored[Ba+2]), 'ro')    #fit
plt.plot(mean_stored[Ba+2] * np.array([1, 1]), [min(counts_stored[:,Ba]), max(counts_stored[:,Ba]) ], 'm--', linewidth=3 )
                            #Vertical line for the centroid
plt.plot(point_left_mean * np.array([1, 1]) , [min(counts_stored[:,Ba]), max(counts_stored[:,Ba]) ], 'k--', linewidth=3 )
plt.plot(point_right_mean * np.array([1, 1]) , [min(counts_stored[:,Ba]), max(counts_stored[:,Ba]) ], 'k--', linewidth=3 )
plt.title("Ba22", fontsize=22)           #title
plt.xlabel("ADC channels", fontsize=14)                        #xlabel
plt.ylabel("Counts", fontsize=14)              #ylabel
# Set size of tick labels.
plt.legend(['fit', 'Centroid', 'FWHM'], fontsize=10) 
plt.tick_params(axis='both', labelsize=14)              #size of axis
plt.grid(True) 
plt.savefig('Signal_Ba_R_visual.png', format='png')




#%% #########################################3############################
#################### 2) Resolution #######################
##########################################################################

#resolution = FWHM/<E> , <E> the centroid of the peak.

R_stored = 100 * FWHM_stored / mean_stored           #channel Resolution [%] 

#calc of delta R:
    
sqrt_sum = np.sqrt( (delta_FWHM_stored / FWHM_stored)**2 + (delta_mean_stored / mean_stored)**2 ) 
                                         #sqrt of the sum of relative errors

delta_R_stored = R_stored * sqrt_sum                                #delta(R[%])



#Plot

#I will rearrange them o that the order is from lower to higher! I must not alter the load order
#since it will affect everything.
            #load order: 'CsI', 'BGO', 'LYSO'
            #order from low to high: 'CsI',  'LYSO', 'BGO'

plt.figure(figsize=(8,5))  #width, heigh 6.4*4.8 inches by default
plt.bar(['Cs','Co1', 'Co2','Na1', 'Na2', 'Ba'], np.absolute( R_stored), 
        yerr = delta_R_stored, edgecolor="black")
plt.title("Resolution of the Cs137 peak", fontsize=22, wrap=True)           #title
plt.ylabel("R (%)", fontsize=14)              #ylabel
plt.tick_params(axis='both', labelsize=14)              #size of axis
plt.grid(True) 
plt.savefig('Resolutions.png', format='png')



#%% #################################################################
################# 3) energy##########################################

#The energy of the peaks used are:
    
Energy = np.array( [661.7, 1173.238, 1332.502, 511, 1274, 383.8])   #[keV]
            #source: nudat, my slides from experimental courses

plt.figure(figsize=(8,5))  #width, heigh 6.4*4.8 inches by default
plt.errorbar(Energy, R_stored, delta_R_stored, fmt='.r', capsize = 5)
plt.title("Resolution vs Energy", fontsize=22, wrap=True)           #title
plt.ylabel("R (%)", fontsize=14)              #ylabel
plt.xlabel('Energy (keV)', fontsize = 14)
plt.tick_params(axis='both', labelsize=14)              #size of axis
plt.grid(True) 

#R = FWHM/<E> ; FWHM = 2 \sqrt(2log2 )sigma,  sigma^2 = <E> (poisson) ==>
  #R = 2 \sqrt(2log2 ) / \sqrt(E) ==> do fit!

#
def linear(x, m, n):       #Definition of the function to use to fit the data
    return m * x + n 

fit_lin_exp = Fits.LinearRegression(np.log(Energy), np.log(R_stored) )

plt.figure(figsize=(8,5))  #width, heigh 6.4*4.8 inches by default
plt.errorbar(Energy, R_stored, delta_R_stored, fmt='.r', capsize = 5)
plt.plot(np.linspace(min(Energy), max(Energy)),np.exp(fit_lin_exp['Slope'] * np.log(np.linspace(min(Energy), max(Energy)))) * np.exp(fit_lin_exp['Intercept']) )      #fit
plt.legend(['fit', 'data'], fontsize=10) 
plt.title("Resolution vs Energy (emu)", fontsize=22, wrap=True)           #title
plt.ylabel("R (%)", fontsize=14)              #ylabel
plt.xlabel('Energy (keV)', fontsize = 14)
plt.tick_params(axis='both', labelsize=14)              #size of axis
plt.grid(True)           
plt.text(600,15, 'ln(R(E)) = {0:1.3f} ln(E) + {1:1.3f} ; r = {2:1.3f}'
         .format(fit_lin_exp['Slope'],fit_lin_exp['Intercept'],fit_lin_exp['r']), fontsize=14)    #first 2 arguments are x,y position.
plt.savefig('R_vs_E_emu.png', format='png')





#%% #################################################
################4)Calibration plot

#for the energy channel resolution, will choose the maximum from the hist,isntead of
#form the fit, in order that the error would be greater (\sqrt counts)


lin_fit = Fits.LinearRegression(ch_peak_stored, Energy)            

#Plot with fit

plt.figure(figsize=(8,5))  #width, heigh 6.4*4.8 inches by default
plt.errorbar(ch_peak_stored, Energy, np.sqrt(ch_peak_stored), fmt='.r', capsize = 5)
plt.plot(ch_peak_stored, [linear(a, lin_fit['Slope'], lin_fit['Intercept']) for a in ch_peak_stored])      #fit
plt.title("Energy calibration (emu)", fontsize=22, wrap=True)           #title
plt.xlabel("ADC channels", fontsize=14)              #ylabel
plt.ylabel('Energy (keV)', fontsize = 14)
plt.legend(['linear fit','data',], fontsize=14)             #legend
plt.tick_params(axis='both', labelsize=14)              #size of axis
plt.grid(True) 
plt.text(1.2e3,400, 'y(x) = {0:1.3f}x + {1:1.3f} ; r = {2:1.3f}'
         .format(lin_fit['Slope'],lin_fit['Intercept'],lin_fit['r']), fontsize=14)    #first 2 arguments are x,y position.
    #0:1.3f: 0 is 1st argument in format, 1.3f means float on 3 decimals
#plt.xlim(5.35,5.55) 
plt.savefig('Energy_calibration_emu.png', format='png')
