#Creating simple individual based model - aquision of Klebsiella among neonates
#Tom Crellen, MORU Postdoc

import itertools
import random
import numpy
import sys
from lifelines import KaplanMeierFitter
import argparse
import matplotlib.pyplot as plt
from matplotlib import interactive

#Argparse
parser = argparse.ArgumentParser(description="Individal Based Model of AMR introduction and spread in a hospital ward")
parser.add_argument('-H', '--height', default=10, required=False, dest="height", metavar="<ward height>", type=int,help="height, ward size is height*width, integer, default=10")
parser.add_argument('-W', '--width', default=10, required=False, dest="width", metavar="<ward width>", type=int,help="width, ward size is height*width, integer, default=10")
parser.add_argument('-I', '--iterations', default=300, required=False, dest="iters", metavar="<number of interations>",type=int,help="number of model iterations (days), integer, default=300")
parser.add_argument('-E', '--entry', default=3, required=False, dest="entry", metavar="<entry rate>",type=int,help="lamda param for poisson entry process, integer, default=3")
parser.add_argument('-M', '--median-stay', default=5, required=False, dest="stay", metavar="<median stay length>",type=int,help="median patient stay (days), integer, default=5")
parser.add_argument('-R', '--risk', default=0.025, required=False, dest="risk", metavar="<risk of transmission>",type=float,help="risk of infection per-person per-day, float, default=0.025")
parser.add_argument('-G', '--graphic', default='f', required=False, dest="graph", metavar="<show ward graphic>",type=str,help="show graphic of ward composition, boolean, default=False")
parser.add_argument('-KM', '--kaplan', default='f', required=False, dest="KM", metavar="<Kaplan-Meier survival output>",type=str,help="Kaplan-Meier table, plot, median or false (default)")
parser.add_argument('-D', '--distribution', default='log-normal', required=False, dest="dist", metavar="<distribution stay length>",type=str,help="Distribution of length stay, log-normal (default), weibull, exponential, gamma")
parser.add_argument('--data', default=None, required=False, dest="data", metavar="<file of stay lengths>",help="file where each line is length of stay (days)")
parser.add_argument('-P', '--prevalence', default="f", required=False, dest="prev", metavar="<show prevalence each iteration>",help="show prevalence each iteration, boolean, default=False")
args = parser.parse_args()

#Boolean Parser
def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

#9 optional command line arguements- assumes default value if not specified
height = args.height
width = args.width
n_iterations = args.iters
entry_rate = args.entry
median_stay = args.stay
risk_transmission = args.risk
image = str2bool(args.graph)
KM = args.KM.lower()
distribution = args.dist.lower()
data= args.data
prevalence = str2bool(args.prev)

#Check list of input stay lengths in data
data_list = []
if data != None:
	with open(data, 'r') as input_data:
		for line in input_data:
			num = float(line.split()[0].strip())
			data_list.append(num)
#Create Klebsiella class
class Klebsiella:
	def __init__(self, height, width, n_iterations, entry_rate, median_stay, risk_transmission, image, KM, distribution, data, prevalence):
		self.height = height
		self.width = width
		self.n_iterations = n_iterations
		self.patients = {}
		self.beds = {}
		self.entry_rate = entry_rate
		self.risk_transmission = risk_transmission
		self.image = image
		self.bed_infected = []
		self.median_stay = median_stay
		self.KM = KM
		self.distribution = distribution
		self.data= data
		self.prevalence = prevalence

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
			remove =[bed for bed, date in self.beds.items() if date[1] == i]
			#Replace beds at the front of empty list
			self.empty_beds = remove + self.empty_beds
			for k in remove:
				#Remove empty bed from patient bed dict
				self.beds.pop(k, None)
			
			#In each day admit n new patients, where n is sampled from poisson
			new_patients = numpy.random.poisson(entry_rate)
			#Check number of new patients is greater than zero and does not exceed bed capacity
			if i>0 and new_patients>0 and len(self.beds.keys())+new_patients <= len(self.all_beds) :
				for n in range(1, new_patients+1):
					#Give unique ID to each patient
					ID = str(i)+'.'+str(n)
					#Discharge day - call dist function - sample from either log-normal, wiebull, gamma or exponential dists
					discharge = dist(self.distribution, self.median_stay, self.data)
					#Infection Status at entry (sampled from binomial)
					entry_status= numpy.random.binomial(1,0.6)
					#Patient bed dict, values are ID, discharge day, infection status at, day of entry
					self.beds[self.empty_beds[0]] = [ID, i+discharge, entry_status, i]
					#patient dictionary for survival analysis
					#values[0:2] remain unchanged
					#value[3] modified if patient becomes infected
					#value[4] gives day of infection, otherwise remains length of stay
					self.patients[ID] = [i, discharge, entry_status, entry_status, discharge]
					#Remove now occupied bed from empty bed list
					self.empty_beds.remove(self.empty_beds[0])
			
			#Create Risk of infection from other patients
			#If at least one infected patient in the ward
			if 1 in [j[2] for j in self.beds.values()]:
				#Check number of infected patients (have been in ward for >0 days)
				n_infec= len(self.bed_infected)
				#Iterate through dictionary of occuped patient beds
				for key, value in self.beds.iteritems():
					#Checks patients are uninfected and have been in the ward for at least one iteration
					if value[2] == 0 and i-value[3] > 0:
						#Risk of infection per day sampled from binomial (eg. 2.5%) - to the power of n infected patients
						new_infect_status = numpy.random.binomial(1,1-(1-self.risk_transmission)**n_infec)
						#If patient newly infected...
						if new_infect_status == 1:
							#Replace patient infection status in bed dict
							value[2] = new_infect_status
							#Replace patient infection status in patient dict
					        	self.patients[value[0]][3] = new_infect_status
							#Replace day of discharge with day of infection (for survival analysis)
							self.patients[value[0]][4] = i-self.patients[value[0]][0] 
						
			#Check infection numbers at the end of each day	
			self.bed_infected = []
                        bed_uninfected = []
			infected_dict = {}
			for key, value in self.beds.iteritems():
				#Infection status of beds
				if value[2] == 1:
				#If positive append to infected bed list
					self.bed_infected.append(key)
				#Else append to uninfected bed list
				else:
					bed_uninfected.append(key)
			
			#dict of staus - for visualisation
			infected_dict.update({el:1 for el in self.empty_beds}) 
			infected_dict.update({el:2 for el in self.bed_infected})
			infected_dict.update({el:3 for el in bed_uninfected})
			#If image flag turned on - use plot function
			if self.image ==True:
				plot(infected_dict, i, self.height, self.width)
			#If prevalence turned on - print values
			if self.prevalence == True:
				print 'time = {}, empty beds = {}, infected beds = {}, uninfected beds = {}'.format(i, len(self.empty_beds), len(self.bed_infected), len(bed_uninfected))

	#Survival Analysis 
	def survival(self):
		#T= length of time before the 'event' (either time to infection or time to discharge), E indicates if infection 'event' took place 
		T= []
		E= []
		for key, value in self.patients.iteritems():
			#If infected at baseline - remove from analysis
			if value[2] == 1:
				pass
			else:
				T.append(value[4])
				E.append(value[3])	
		
		#Kaplan Meier Estimator from lifelines package
		kmf= KaplanMeierFitter()
		E_ar = numpy.array(E)
		T_ar = numpy.array(T)
		
		#Fit Kaplan Meier model and plot
		kmf.fit(T_ar, event_observed=E_ar)
		if self.KM == "table":
			print kmf.survival_function_
		elif self.KM == "plot":
			kmf.plot()
                	plt.show()	
		elif self.KM == "median":
			print kmf.median_
		else:
			pass

#Set distribution of length of stay
def dist(d, median, data):
	if d == "log-normal":
		return int(numpy.clip(numpy.log(numpy.random.lognormal(median,2)), 1, None))
	elif d == "exponential":
		return int(numpy.clip(numpy.random.exponential(median), 1, None))
	elif d == "gamma":
		return int(numpy.clip(numpy.random.gamma(median), 1, None))
	elif d== "weibull":
		return int(numpy.clip(median*numpy.random.weibull(1), 1, None))
	elif d== "data":
		return int(numpy.clip(random.choice(data), 1, None))
	else:
		raise ValueError("Distribution must be log-normal, gamma, exponential, weibull or data")

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
run_1 = Klebsiella(height, width, n_iterations, entry_rate, median_stay, risk_transmission, image, KM, distribution, data_list, prevalence)
run_1.populate()
run_1.admit()
run_1.survival()
