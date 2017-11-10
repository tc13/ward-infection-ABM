#Creating simple individual based model - aquision of Klebsiella among neonates
#Tom Crellen, MORU Postdoc

import itertools
import random
import copy
import numpy
import sys
import math
from lifelines import KaplanMeierFitter
import argparse
import matplotlib.pyplot as plt
from matplotlib import interactive

#Argparse
parser = argparse.ArgumentParser(description="Individal Based Model of AMR introduction and spread in a hospital ward")
parser.add_argument("ward-height", metavar="H",type=int,help="height, ward size is height*width, integer")
parser.add_argument("ward-width", metavar="W",type=int,help="width, ward size is height*width, integer")
parser.add_argument("iterations", metavar="I",type=int,help="number of model iterations (days), integer")
parser.add_argument("entry-rate", metavar="E",type=int,help="lamda param for poisson entry process, integer")
parser.add_argument("risk-transmission", metavar="R",type=float,help="risk of infection within ward per-person per-day, float")
parser.add_argument("graphic", metavar="G",type=bool,help="show graphic of ward composition, boolean")
args = parser.parse_args()

#Boolean Parser
def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

#Takes 5 command line arguements- height and width of ward, number of iterations, entry rate param and risk of infection
height = int(sys.argv[1])
width = int(sys.argv[2])
n_iterations = int(sys.argv[3])
entry_rate = int(sys.argv[4])
risk_transmission = float(sys.argv[5])
image = str2bool(sys.argv[6])

#Create Klebsiella class - 5 arguments are inherited from the command line
class Klebsiella:
	def __init__(self, height, width, n_iterations, entry_rate, risk_transmission, image):
		self.height = height
		self.width = width
		self.n_iterations = n_iterations
		self.patients = {}
		self.patient_bed = {}
		self.entry_rate = entry_rate
		self.risk_transmission = risk_transmission
		self.image = image

	#Function to populate the ward with beds of given coordinates
	def populate(self):
		#list of length height*width, each element is unique (height, width)
		self.all_beds = list(itertools.product(range(self.width),range(self.height)))
		#As all beds are empty at start, empty beds list equal to all beds
		self.empty_beds = self.all_beds

	#Admit patients to wards
	def admit(self):
		#Iterate through time intervals (days)
		for i in range(self.n_iterations):
			#Remove outgoing patients (if discharge day == i)
			remove = [bed for bed, date in self.patient_bed.items() if date[1] == i]
			#Replace beds in empty list
			self.empty_beds.extend(remove)
			for k in remove:
				#Remove empty bed from patient bed dict
				self.patient_bed.pop(k, None)
			
			#In each day admit n new patients, where n is sampled from poisson
			new_patients = numpy.random.poisson(entry_rate)
			for n in range(new_patients):
				#Check number of new patients is greater than zero and does not exceed bed capacity
				if i>0 and n>0 and len(self.patient_bed.keys())+new_patients <= len(self.all_beds) :
					#Give unique ID to each patient
					ID = str(i)+'.'+str(n)
					#Discharge day - sample from lognormal distribution (and logged) +1 so no zeros
					discharge = int(math.log(numpy.random.lognormal(5,2)))+1
					#Infection Status at entry (sampled from binomial)
					entry_status= numpy.random.binomial(1,0.6)
					#Patient bed dict, values are ID, discharge day, infection status at, day of entry
					self.patient_bed[self.empty_beds[0]] = [ID, i+discharge, entry_status, i]
					#patient dictionary for survival analysis
					self.patients[ID] = [i, discharge, entry_status, entry_status, discharge]
					#Remove now occupied bed from empty bed list
					self.empty_beds.remove(self.empty_beds[0])
			
			#Create Risk of infection from other patients
			#If at least one infected patient in the ward
			if 1 in [j[2] for j in self.patient_bed.values()]:
				#Iterate through dictionary of occuped patient beds
				for key, value in self.patient_bed.iteritems():
					#Checks patients are uninfected and have been in the ward for at least one iteration
					if value[2] == 0 and i-value[3] > 0:
						#Risk of infection per day sampled from binomial (eg. 5%)
						new_infect_status = numpy.random.binomial(1,self.risk_transmission)
						#If patient newly infected...
						if new_infect_status == 1:
							#Replace patient infection status in bed dict
							value[2] = new_infect_status
							#Replace patient infection status in patient dict
					        	self.patients[value[0]][3] = new_infect_status
							#Replace day of discharge with day of infection (for survival analysis)
							self.patients[value[0]][4] = i-self.patients[value[0]][0] 
						
			#Check infection numbers at the end of each day	
			bed_infected = []
                        bed_uninfected = []
			infected_dict = {}
			for key, value in self.patient_bed.iteritems():
				#Infection status of beds
				if value[2] == 1:
				#If positive append to infected bed list
					bed_infected.append(key)
				#Else append to uninfected bed list
				else:
					bed_uninfected.append(key)
			
			#dict of staus - for visualisation
			infected_dict.update({el:1 for el in self.empty_beds}) 
			infected_dict.update({el:2 for el in bed_infected})
			infected_dict.update({el:3 for el in bed_uninfected})
			
			if self.image ==True:
				plot(infected_dict, i, self.height, self.width)
		
			#print 'time = {}, empty beds = {}, infected beds = {}, uninfected beds = {}'.format(i, len(self.empty_beds), len(bed_infected), len(bed_uninfected)) 
	#Survival Analysis 
	def survival(self):
		#T is length of time before the 'event'
		T= []
		#E is the 'event'
		E= []
		#print 'id, time, infected'
		for key, value in self.patients.iteritems():
			#If infected at baseline - remove from analysis
			if value[2] == 1:
				pass
			else:
				T.append(value[4])
				E.append(value[3])	
		#		print '{}, {}, {}'.format(key, value[1], value[3])
		#Kaplan Meier Estimator from lifelines package
		kmf= KaplanMeierFitter()
		E_ar = numpy.array(E)
		T_ar = numpy.array(T)
		
		#Fit Kaplan Meier model and plot
		kmf.fit(T_ar, event_observed=E_ar)
		#print kmf.survival_function_
		kmf.plot()
		plt.show()

#Plot ward grid each iteration
def plot(infection_dict, i, height, width):
		interactive(True)
                title = 'Time = {}'.format(i)
                fig, ax = plt.subplots()
                patient_colors = {1:'tab:gray', 2:'tab:red', 3:'tab:green'}
                #print self.infection_dict
                for patient in infection_dict:
                        ax.scatter(patient[0]+0.5, patient[1]+0.5, color=patient_colors[infection_dict[patient]])
                ax.set_title(title, fontsize=10, fontweight='bold')
                ax.set_xlim([0, width])
                ax.set_ylim([0, height])
                ax.set_xticks([])
                ax.set_yticks([])
                #plt.show()
		plt.plot()
		raw_input('press return to continue')


#Run model - params inherited from command line arguments
run_1 = Klebsiella(height, width, n_iterations, entry_rate, risk_transmission, image)
run_1.populate()
run_1.admit()
run_1.survival()
