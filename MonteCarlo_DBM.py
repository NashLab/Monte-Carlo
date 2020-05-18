# This script simulates the SMFS forced pulling process of a CohDoc compelx
# with dual binding behavior, under a constant pulling speed.
# 
# The energy profiles of the complex and the X module of Doc to estimate the 
# dissiocation process are given from AFM SMFS measurements, fitted using the 
# Bell-Evans model.
# 
# - Bending correction that takes into account the deflection of the AFM cantilever
# 	is conducted along the extension axis.
# 
# - More information could be found on https://github.com/NashLab/Monte-Carlo
# 
# 					Last update: May.15th 2020
# 					by Haipei Liu, Nash Lab

from matplotlib import pyplot as plt
from scipy import stats
import scipy
import numpy as np
import seaborn as sns
import matplotlib as mpl
import time as timer
import logging as log
import os, os.path

np.random.seed(None)

##	HERE ARE ALL THE INPUT PARAMETERS
#	Units : ms nm pN

kT = 4.14 		# Boltzmann constant x T 		in	[pN*nm]
k = 91. 		# Spring constant of AFM probe 	in	[pN/nm]

linker = 174	# spacer length 				in	[nm]
# if dispersity is True, a random number drawn from a gaussian dist. around zero with this sigma will be added to spacer
dispers = 0
persist = 0.365 # persist length in	[nm]	

#syntax params= [contL,  dx, 	k0]
xdoc_w = 	[204, 0.366, 2.79e-5]	# complex rupture under weak mode
xdoc_s = 	[204, 0.178, 4.7e-8 ]	# complex rupture under strong mode
xdoc_s_xMU = 	[204, 0.277, 7.43e-4] # strong mode complex rupture after Xmod unfolding
xdoc = ("xdoc",xdoc_w,xdoc_w,xdoc_s,xdoc_s_xMU)

# xMod[1] xMod unfolding under strong mode
xMod = ("xMod", [ 38, 0.116, 4.53e-5],[38, 0.000000016, 4.53e-8])

load = 1000 #in pN/s
speed = 800 #in nm/s

# ratio defined for dual binding behavior
ratio_of_weak_bound = 0.2

# formatting for graphics
fontsize = 24
sns.set(style="ticks",font="Arial")
colors = sns.color_palette()

def format_figure(axis):
	axis.tick_params(axis = "both", direction = "in", width = 1.2)
	for item in [axis.title, axis.xaxis.label, axis.yaxis.label]:
		item.set_fontsize(fontsize)
	for item in (axis.get_xticklabels() + axis.get_yticklabels()):
		item.set_fontsize(fontsize-2)

#probability density function according to the BE-model
"""
Variable 
	dx	intermolecular separation
	k0	unperturbed transition rate
	r	loading rate in pN/s
"""
def pdfBell(force,r,dx,k0):
	"""probability density function in the BE model"""
	pdf = k0/r * np.exp(dx*force/kT - k0*kT*(np.exp(dx*force/kT)-1)/(dx*r))
	return pdf

def intbell(force,r,dx,k0):
	"""Integrated Bell pdf, cumulative distribtion function
	the probability that the unbinding of pair(dx,k0) happens at FORCE"""
	ipdf = 1.0-np.exp(-(k0*kT/(r*dx)*(np.exp((dx*force)/kT)-1.0)))
	return ipdf

#off rate in the BE model:
def rate(force,dx,k0):
	"""off rate in the BE model"""
	rate = k0 * np.exp( np.minimum(dx*force/kT,50.))
	return rate

# WLC model
"""
Variables
	x	extention
	lp  persist length
"""
def f_WLC(x,contL,lp):
	"""force extension behavior in the WLC model"""
	wlc = kT/lp * (1/(4.0*(1.0-np.minimum(x,0.9999*contL)/contL)**2.0) + x/contL - 1.0/4.0)
	return wlc

# generation of time axis with extension and force
def generate_axis(time,step,speed,noiseval,lc,time_from):
	timeaxis = np.linspace(time_from,time,np.int(1+(time-time_from)/step))
	ext = timeaxis*speed
	#setup a baseline with some gaussian noise
	fnoise = np.zeros(len(timeaxis))
	if noiseval == 0:
		fnoise_base = np.zeros(len(timeaxis))
	else:
		fnoise_base = np.random.normal(0,noiseval,len(timeaxis))
	# in constant speed, force is calculated using th WLC model
	fnoise = fnoise_base + f_WLC(ext,lc,persist)
	# bending correction
	dist = ext+fnoise/k
	timeaxis = dist/speed

	i = 0
	#Produce a finer axis when the force extension resolution is not enough
	resolution_F = 1.0
	while ext[i] < lc*0.95:
		if (fnoise[i+1]-fnoise[i]) > resolution_F :
			n_step= 2* (lc*0.97-ext[i])/(ext[i+1]-ext[i])
			ext= np.concatenate((ext[:i],np.linspace(ext[i],lc*0.97,np.int(n_step))),axis= None)
			timeaxis= ext/speed
			fnoise_c = np.zeros(len(timeaxis))
			if noiseval == 0:
				fnoise_base = np.zeros(len(timeaxis))
			else:
				fnoise_base = np.random.normal(0,noiseval,len(timeaxis))
			fnoise_c = fnoise_base + f_WLC(ext,lc,persist)
			fnoise= np.concatenate((fnoise[:i],fnoise_c[i:]),axis= None)
			dist = ext+fnoise/k
			timeaxis = dist/speed
		i += 1
	return timeaxis,fnoise,ext,dist

# Monte Carlo simulation of a constant pulling speed test
def simulate_constspeed(time,step,noiseval,spacer,cpl,cpl_2,fingerprint,speed,k,dispers=0):
	"""function to simulate constant speed experiments"""

	#polydispersity if needed
	if dispers != 0:
		spacer += np.random.normal(0,dispers)

	#setup a time and extension axis
	timeaxis, fnoise, ext, dist = generate_axis(time,step,speed,noiseval,spacer, time_from=0)

	i_rup = 0
	i_unf = 0

	if fingerprint != None:
		p_rup = np.zeros(len(timeaxis))
		p_unf = np.zeros(len(timeaxis))
		lr_rup = 0.
		f_rupture = 0.
		f_unfold= 0
		lr_unf = 0. 
		i = 10
		i_print = 0
		while fnoise[i] < 50. :
			i+=1
		while i < len(timeaxis):
			step_bc = (timeaxis[i+1]-timeaxis[i])
			p_unf[i] = 1-np.exp(-rate(fnoise[i],*fingerprint[1:])*step_bc)
			p_rup[i] = 1-np.exp(-rate(fnoise[i],*cpl[1:])*step_bc)
			pick_rup = np.random.random_sample()
			pick_unf = np.random.random_sample()

			"""random test for complex rupture or unfolding of Xmod"""
			if pick_unf > p_unf[i] and pick_rup > p_rup[i]: #both survive, nothing happens								
				i+=1
			elif pick_unf < p_unf[i] and pick_rup > p_rup[i]: #fingerprint unfolds, experiment continues
				f_unfold = fnoise[i]
				lr_unf = (fnoise[i] - fnoise[i-1]) / step_bc
				i+=1
				#update the axis for the Lc increment
				timeaxis_c, fnoise_c, ext_c, dist_c = generate_axis(time,step,speed,noiseval,spacer+fingerprint[0],timeaxis[i])
				timeaxis= np.concatenate((timeaxis[:(i-1)],timeaxis_c),axis= None)
				fnoise= np.concatenate((fnoise[:(i-1)],fnoise_c),axis= None)
				ext= np.concatenate((ext[:(i-1)],ext_c),axis= None)
				dist= np.concatenate((dist[:(i-1)],dist_c),axis= None)
				p_rup_fpunf = 1-np.exp(-rate(fnoise,*cpl_2[1:])*step)
				f_rupture, lr_rup = compare_nofp(timeaxis,fnoise,p_rup_fpunf,i_nofp=i)

				while fnoise[i]<f_rupture:
					i+=1
				fnoise[i+1:] = fnoise[i+1:]*0.

				return f_rupture, f_unfold, ext, fnoise, lr_rup, lr_unf
				break
			else: #complex ruptures with fingerprint intact, experiment ends
				lr_rup = (fnoise[i] - fnoise[i-1]) / (timeaxis[i] - timeaxis[i-1])
				f_rupture = fnoise[i]
				f_unfold = 0
				fnoise[i+1:] = fnoise[i+1:]*0.
				return f_rupture, f_unfold, ext, fnoise, lr_rup, lr_unf
				break
	else:
		f_rupture = compare_nofp(timeaxis,fnoise,p_rup)
		f_unfold = 0
		return f_rupture, f_unfold, ext, fnoise, lr_rup, lr_unf	#, i_rup,i_unf

def compare_nofp(timeaxis,fnoise,p_rup,i_nofp=0):
	"""function to simulate the rupture without any further protein unfolding"""
	f_rupture_nofp = 0.
	rup_lr_nofp = 0.
	while i_nofp < len(timeaxis):
		pick = np.random.random_sample() 
		if pick < p_rup[i_nofp]:
			f_rupture_nofp = fnoise[i_nofp]
			rup_lr_nofp = (fnoise[i_nofp] - fnoise[i_nofp-1]) / (timeaxis[i_nofp] - timeaxis[i_nofp-1])
			break
		else:
			i_nofp+=1
	return f_rupture_nofp, rup_lr_nofp

def run_simulation(p_w,num,time,step,noiseval,spacer,cpl,ldr,fingerprint,forceramp=True,verbose=False,speed=speed,k=k,savecurves=False):
	rup_forces = np.zeros(num)
	unf_forces = np.zeros(num)
	rup_lr = np.zeros(num)
	unf_lr = np.zeros(num)

	i_a, i_b, i_c =0,0,0
	i_groups= np.zeros(num)

	if forceramp == True:
		method = "ramp"
	else:
		method = "speed"
	
	#setup directories for saving files
	current_dir = os.getcwd()
	pathname = "DBmode_"+method+"_"+cpl[0]+"_ratio-"+str(p_w)
	data_dir = os.path.join(current_dir, pathname)
	if not (os.path.exists(data_dir)):
		os.mkdir(pathname)
	os.chdir(data_dir)

	#log all parameters of the simulation
	print ("Saving simulation parameters")
	print ("Simulating {0} curves".format(num))
	log.basicConfig(filename="speed_"+str(speed)+"logfile.log",filemode='w',level=log.DEBUG)


	log.info("num: {0:0}".format(num))
	log.info("time: {0} sec".format(time))
	log.info("step: {0} sec".format(step))
	log.info("noise: {0} pN".format(noiseval))
	log.info("complex: {0}".format(cpl))
	log.info("fingerprint: {0}".format(fingerprint))
	log.info("forceramp: {0}".format(forceramp))
	if forceramp == True:
		log.info("ldr: {0} pN/s".format(ldr))
	else:
		log.info("speed: {0} nm/s".format(speed))
		log.info("spacer: {0} nm".format(spacer))
	
	#setup directory to save force-ext traces
	if savecurves == True:
		#print "Saving force-ext-traces"
		curve_dir = os.path.join(data_dir, "curves")
		if not (os.path.exists(curve_dir)):
			os.mkdir("curves")
		os.chdir(curve_dir)

	print ("Starting constant speed with {0} nm/s".format(speed))
	print ("Sampling {0:.1f} points/nm".format(1/(speed*step)))
	starttime2 = timer.perf_counter()

	for i in range(num):
		pick_boundstate = np.random.random_sample()
		#	test bound state
		if pick_boundstate < p_w:
			rup_forces[i], unf_forces[i], extension, force, rup_lr[i], unf_lr[i] = simulate_constspeed(time,step,noiseval,spacer,cpl[1],cpl[2],fingerprint[2],speed,k)
			i_groups[i] = 3
			i_a +=1
		else:
			rup_forces[i], unf_forces[i], extension, force, rup_lr[i], unf_lr[i] = simulate_constspeed(time,step,noiseval,spacer,cpl[3],cpl[4],fingerprint[1],speed,k)
			
			if unf_forces[i] != 0:
				i_groups[i] = 2
				i_c +=1
			else:
				i_groups[i] = 1
				i_b +=1
		if savecurves == True:
			np.savetxt("ext-force-{0}.txt".format(i),np.transpose([extension,force]))
			distance = extension+force/k
			plot_fed_trace(extension,distance, force)
			plt.savefig("ext-force-trace-{0}.pdf".format(i))
			plt.show()

	endtime2 = timer.perf_counter()
	print ("Computation took {0:.2f} sec".format(endtime2-starttime2))
	unf_forcesclean = unf_forces[unf_forces!=0]
	unf_lrclean = unf_lr[unf_forces!=0]

#Group 1 High force rupture (Xmod intact)
	rup_forcesclean_b = rup_forces[i_groups ==1]
#Group 2 High force rupture (Xmod unfolded)
	rup_forcesclean_c = rup_forces[i_groups ==2]
#Group 3 Low force rupture
	rup_forcesclean_a = rup_forces[i_groups ==3]

	rup_forcesclean = rup_forces[rup_forces!=0]
	average_lr = np.sum(rup_lr)/num
	average_rf = np.sum(rup_forces)/num

#Output
	print ("-------{0} nm/s-------".format(speed))
	print ("{0}/{1} curves showed low force ({2:.1f}%)".format(i_a,num,100*i_a/num))
	print ("{0}/{1} curves showed high force without Xmod unfolding ({2:.1f}%)".format(i_b,num,100*i_b/num))
	print ("{0}/{1} curves showed high force with Xmod unfolded ({2:.1f}%)".format(i_c,num,100*i_c/num))
	print ("In total: {0} , LR:{1}, RF:{2}".format(num, average_lr, average_rf))

	log.info("{0}/{1} curves showed low force ({2:.1f}%)".format(i_a,num,100*i_a/num))
	log.info("{0}/{1} curves showed high force without Xmod unfolding ({2:.1f}%)".format(i_b,num,100*i_b/num))
	log.info("{0}/{1} curves showed high force with Xmod unfolded ({2:.1f}%)".format(i_c,num,100*i_c/num))
	log.info("In total(from MonteCarlo): {0} , LR:{1}, RF:{2}".format(num, average_lr, average_rf))
	os.chdir(data_dir)

#Plot resulting distributions
	maxis = 50*round(np.amax(rup_forces)/50)+50
	forceax = np.linspace(0,1000,200)
	plotbins = [20,60,100,140,180,220,260,300,340,380,420,460,500,540,580,620,660,700,740,780,820]
	
	fig, ax = plt.subplots()
	if fingerprint == None:
		sns.distplot(rup_forces,ax=ax,label= "sim ruptures")
	else:	
		cp_dist = pdfBell(forceax,average_lr,*cpl[1][1:]) *p_w
		cp_dist_s = pdfBell(forceax,average_lr,*cpl[3][1:]) *(1-p_w)

		sns.distplot(rup_forcesclean_a,ax=ax,bins=plotbins,label= "low force",kde=False,norm_hist=False,color=colors[3])
		sns.distplot(rup_forcesclean_b,ax=ax,bins=plotbins,label= "high force (Xmod intact)",kde=False,norm_hist=False,color=colors[4])
		sns.distplot(rup_forcesclean_c,ax=ax,bins=plotbins,label= "high force (Xmod unfolded)",kde=False,norm_hist=False,color=colors[2])
	ax.set(xlabel="Force [pN]",ylabel="Probabilty density [1/pN]",xlim=0)
	ax.legend(loc ="upper right", prop={"size": fontsize-8})
	for item in fig.get_axes():
		format_figure(item)
	fig.tight_layout()
	plt.savefig("analysis-hist-ratio_"+str(p_w)+"speed_"+str(speed)+".pdf")

	fig, ax = plt.subplots()
	sns.distplot(unf_forcesclean,ax=ax,bins=plotbins,label= "Xmod unfolding",kde=False,norm_hist=False,color=colors[3])
	ax.legend(loc ="upper right", prop={"size": fontsize-8})
	fig.tight_layout()
	plt.savefig("Xmod_unfolding-hist-ratio_"+str(p_w)+"speed_"+str(speed)+".pdf")

	os.chdir(data_dir)
	filename = "LR-RF_v"+str(speed)+".txt"
	rup_lr_si= rup_lr * 1e-12
	rup_forces_si= rup_forces * 1e-12
	np.savetxt(filename,np.transpose([rup_lr_si, rup_forces_si, i_groups]),header="LR rup_force class")
#	save xmod unfolding result
	filename = "Xmod_v"+str(speed)+".txt"
	unf_forces_si = unf_forcesclean * 1e-12
	unf_lr_si = unf_lrclean * 1e-12
	np.savetxt(filename,np.transpose([unf_lr_si, unf_forces_si]),header="LR rup_force")
	os.chdir(current_dir)
	return 0

def plot_fed_trace(x,d,f):
	fig, ax = plt.subplots()
	ax.plot(x,f,label="Tip-sample sepration",marker=",")
	ax.plot(d,f,label="AFM head height",marker=",")
	ax.set_xlim(0, 400)
	ax.set(xlabel="Distance [nm]",ylabel="Force [pN]")
	ax.legend(loc ="upper right", prop={"size": fontsize-6})
	for item in fig.get_axes():
		format_figure(item)
	fig.tight_layout()
	plt.show()
	return 0

#RUN YOUR SIMULATION HERE
"""
Variable 
	num		num of forced pulling tests wanted
	time	time for pulling period
	p_w	the ratio of weakly_bound/all_bound
	speed 	applied pulling speed
	k		probe spring constant 

Choose Fingerprint
	xMod is used as fingerprint
"""

run_simulation(p_w=0.2,num=1000,time=4.,step=0.0008,noiseval=0.,spacer=linker,cpl=xdoc,ldr=load,fingerprint=xMod,forceramp=False,speed=100.,k=91.,savecurves=False)

