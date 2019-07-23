#Creating simple individual based model - aquision of Klebsiella among neonates
#Tom Crellen, MORU Postdoc

import itertools
import random
import numpy
import sys
import argparse
from lifelines import KaplanMeierFitter
#import scipy.stats as stats

#Argparse
parser = argparse.ArgumentParser(description="Individal Based Model of AMR introduction and spread in a hospital ward. \n \n Author Tom Crellen (tomcrellen@gmail.com) MORU Postdoc")
parser.add_argument('-H', '--height', default=10, required=False, dest="height",metavar="", type=int,help="height, ward size is height*width, integer, default=10")
parser.add_argument('-W', '--width', default=10, required=False, dest="width",metavar="", type=int,help="width, ward size is height*width, integer, default=10")
parser.add_argument('-I', '--iterations', default=300, required=False, dest="iters",metavar="", type=int,help="number of model iterations (days), integer, default=300")
parser.add_argument('-E', '--entry', default=3, required=False, dest="entry", metavar="", type=int,help="lamda param for poisson entry process, integer, default=3")
parser.add_argument('-B', '--beta', default=0.025, required=False, dest="risk", metavar="", type=float,help="risk of infection per-person per-day, float, default=0.025")
parser.add_argument('-ER', '--entryRisk', default=0.6, required=False, dest="entry_risk", metavar="", type=float,help="risk of infection at entry, float, default=0.6")
parser.add_argument('-P', '--prevalence', default="f", required=False, dest="prev", metavar="", help="show prevalence each iteration, boolean, default=False")
parser.add_argument('-KM', '--kaplan', default='f', required=False, dest="KM", metavar="",type=str,help="Kaplan-Meier table, plot, median or false (default)")
parser.add_argument('-A', '--average-stay', default=5, required=False, dest="stay", metavar="", type=int,help="depending on distribution, median (log-normal, uniform), mean (exponential) or scale (gamma, weibull),  integer, default=5")
parser.add_argument('-D', '--distribution', default='log-normal', required=False, dest="dist",metavar="", type=str,help="Distribution of stay lengths, log-normal (default), weibull, exponential, gamma, uniform or data (set file path with --data flag)")
parser.add_argument('-P2', '--parameter2', default=2, required=False, dest="param2",metavar="",type=int,help="Second parameter, variance for log-normal and shape parameter for weibull and gamma, integer, default=2")
parser.add_argument('--data', default=None, required=False, dest="data",metavar="",help="If data specified with -D, path to file where each line is length of stay (days)")
args = parser.parse_args()

#Boolean Parser
def str2bool(v):
        if v.lower() in ('yes', 'true', 't', 'y', '1'):
                return True
        elif v.lower() in ('no', 'false', 'f', 'n', '0'):
                return False
        else:
                raise argparse.ArgumentTypeError('Boolean value expected.')

#Optional command line arguments- assumes default value if not specified
height = args.height
width = args.width
n_iterations = args.iters
entry_rate = args.entry
median_stay = args.stay
risk_transmission = args.risk
KM = args.KM.lower()
distribution = args.dist.lower()
data= args.data
prevalence = str2bool(args.prev)
param2 = args.param2
entry_risk = args.entry_risk

#Check list of input stay lengths in data
data_list = []
if data != None:
	with open(data, 'r') as input_data:
		for line in input_data:
			num = float(line.split()[0].strip())
			data_list.append(num)

#Create Klebsiella class
class Klebsiella:
	def __init__(self, height, width, n_iterations, entry_rate, median_stay, risk_transmission, KM, distribution, data, prevalence, param2, entry_risk):
		self.height = height
		self.width = width
		self.n_iterations = n_iterations
		self.patients = {}
		self.beds = {}
		self.entry_rate = entry_rate
		self.risk_transmission = risk_transmission
		self.contacts = []
		self.median_stay = median_stay
		self.KM = KM
		self.distribution = distribution
		self.data= data
		self.prevalence = prevalence
                self.param2 = param2
                self.entry_risk = entry_risk
                self.bed_infected = []

	#Function to populate the ward with beds of given coordinates
	def populate(self):
		#list of length height*width, each element is unique (height, width) - all beds empty at start
		self.empty_beds = list(itertools.product(range(self.width),range(self.height)))

	#Admit patients to wards
	def admit(self):
                #distribution of length of stay- sampled from log-normal, wiebull, gamma, exponential, uniform or user entered distributions
                length_stay_dist = dist(self.distribution, self.median_stay, self.data, self.param2)
                #Distribution of infection risk at entry
                infect_entry_risk_dist = list(numpy.random.binomial(1, self.entry_risk, size=10000)) 
		#Iterate through time intervals (days)
		for i in range(self.n_iterations):
			#Remove outgoing patients (if discharge day == i)
                        if i > 0:
			        remove = [bed for bed, date in self.beds.items() if date[1] == i]
			        #Replace beds at the front of empty list
			        self.empty_beds = remove + self.empty_beds
			        for k in remove:
				        #Remove empty bed from patient bed dict
				        self.beds.pop(k, None)
                                        try:
                                                self.bed_uninfected.remove(k)
                                        except ValueError:
                                                self.bed_infected.remove(k)



			#In each day admit n new patients, where n is sampled from poisson
			new_patients = int(numpy.random.poisson(self.entry_rate))
			#Check number of empty beds and cap maximum entrants
                        if new_patients > len(self.empty_beds):
                                new_patients = len(self.empty_beds)
                        #Check number of new patients is greater than zero
			if i>0 and new_patients>0 :
			        for n in range(1, new_patients+1):
					#Give unique ID to each patient
					ID = str(i)+'.'+str(n)
					#Sample discharge list from distribution
					discharge = int(numpy.random.choice(length_stay_dist))
					#Infection Status at entry (sampled from binomial distribution)
					entry_status= numpy.random.choice(infect_entry_risk_dist)
					#Patient bed dict, values are ID, discharge day, infection status at, day of entry
					self.beds[self.empty_beds[0]] = [ID, i+discharge, entry_status, i]
					#patient dictionary for survival analysis
					#[0:2] remain unchanged,[3] modified if patient becomes infected
					#value[4] gives day of infection, otherwise remains length of stay
					self.patients[ID] = [i, discharge, entry_status, entry_status, discharge]
					#Remove now occupied bed from empty bed list
					self.empty_beds.remove(self.empty_beds[0])
			
			#Risk of infection from other patients, if at least one infected patient in the ward
                        if len(self.bed_infected) > 0:
                                #iterate through uninfected patients
                                for p in self.bed_uninfected:
                                        #geometric distribution determines number of trials before infection
                                        geom= numpy.random.geometric(1-(1-self.risk_transmission))
                                        #If more trials than infectious individuals
                                        if geom > len(self.bed_infected):
                                                #No transmission
                                                pass
                                        else:
                                                #q is bed of infector; allow access to 0th element
                                                q = self.bed_infected[geom-1]
                                                #check ID of infected patient
                                                patient_id = self.beds[p][0]
                                                #check ID of infector
                                                infector_id = self.beds[q][0]
                                                #Replace patient infection status in bed dict
                                                self.beds[p][2] = 1
                                                #Replace patient infection status in patient dict
                                                self.patients[patient_id][3] = 1
                                                #Replace day of discharge with day of infection (for survival analysis)
                                                self.patients[patient_id][4] = i-self.patients[patient_id][0]
                                                #Add transmission event to contact list
                                                self.contacts.append([infector_id,patient_id])
                                                                    						
			#Check infection numbers at the end of each day	
			self.bed_infected = []
                        self.bed_uninfected = []
			for key, value in self.beds.iteritems():
				#Infection status of beds
				if value[2] == 1:
				#If positive append to infected bed list
					self.bed_infected.append(key)
				#Else append to uninfected bed list
				else:
					self.bed_uninfected.append(key)

			#If prevalence turned on - print values
			if self.prevalence == True:
                                St = len(self.bed_uninfected)
                                It = len(self.bed_infected)
				print '{},{},{},{}'.format(i, St, It, St+It)

	#Survival Analysis 
        def survival(self):
		#T= length of time before the 'event' (either time to infection or time to discharge), 
                #E indicates if infection 'event' took place 
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
def dist(d, median, data, param2):
	if d == "log-normal":
		return list(numpy.clip(numpy.log(numpy.random.lognormal(median, param2, size=10000)), 1, None))
	elif d == "exponential":
		return list(numpy.clip(numpy.random.exponential(median, size=10000), 1, None))
	elif d == "gamma":
		return list(numpy.clip(numpy.random.gamma(param2, median, size=10000), 1, None))
	elif d== "weibull":
		return list(numpy.clip(median*numpy.random.weibull(param2, size=10000), 1, None))
        elif d== "uniform":
                return list(numpy.random.uniform(1, (median-0.5)*2, size=10000))
	elif d== "data":
		return list(data)
	else:
		raise ValueError("Distribution must be log-normal, gamma, exponential, weibull, uniform or data")

#Run model - params inherited from command line arguments
run_1 = Klebsiella(height, width, n_iterations, entry_rate, median_stay, risk_transmission, KM, distribution, data_list, prevalence, param2, entry_risk)
run_1.populate()
print 'time,S,I,N'
run_1.admit()
if KM in ["plot", "median", "table"]:
        run_1.survival()
