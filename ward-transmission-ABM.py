#Creating simple individual based model - aquision of Klebsiella among neonates
#Tom Crellen, MORU Postdoc

import itertools
import random
import numpy
import sys
from lifelines import KaplanMeierFitter
import argparse
import matplotlib.pyplot as plt
#from matplotlib import interactive
import networkx as nx
from scipy.optimize import minimize
import scipy.stats as stats

#Argparse
parser = argparse.ArgumentParser(description="Individal Based Model of AMR introduction and spread in a hospital ward. \n \n Author Tom Crellen (tomcrellen@gmail.com) MORU Postdoc")
parser.add_argument('-H', '--height', default=10, required=False, dest="height", metavar="<ward height>", type=int,help="height, ward size is height*width, integer, default=10")
parser.add_argument('-W', '--width', default=10, required=False, dest="width", metavar="<ward width>", type=int,help="width, ward size is height*width, integer, default=10")
parser.add_argument('-I', '--iterations', default=300, required=False, dest="iters", metavar="<number of interations>",type=int,help="number of model iterations (days), integer, default=300")
parser.add_argument('-E', '--entry', default=3, required=False, dest="entry", metavar="<entry rate>",type=int,help="lamda param for poisson entry process, integer, default=3")
parser.add_argument('-R', '--risk', default=0.025, required=False, dest="risk", metavar="<risk of transmission>",type=float,help="risk of infection per-person per-day, float, default=0.025")
parser.add_argument('-ER', '--entryRisk', default=0.6, required=False, dest="entry_risk", metavar="<entry infection risk>",type=float,help="risk of infection at entry, float, default=0.6")
parser.add_argument('-P', '--prevalence', default="f", required=False, dest="prev", metavar="<show prevalence each iteration>",help="show prevalence each iteration, boolean, default=False")
parser.add_argument('-G', '--graphic', default='f', required=False, dest="graph", metavar="<show ward graphic>",type=str,help="show graphic of ward composition, boolean, default=False")
parser.add_argument('-KM', '--kaplan', default='f', required=False, dest="KM", metavar="<Kaplan-Meier survival output>",type=str,help="Kaplan-Meier table, plot, median or false (default)")
parser.add_argument('-A', '--average-stay', default=5, required=False, dest="stay", metavar="<average stay length>",type=int,help="depending on distribution, median (log-normal, uniform), mean (exponential) or scale (gamma, weibull),  integer, default=5")
parser.add_argument('-D', '--distribution', default='log-normal', required=False, dest="dist", metavar="<distribution-stay>",type=str,help="Distribution of stay lengths, log-normal (default), weibull, exponential, gamma, uniform or data (set file path with --data flag)")
parser.add_argument('-P2', '--parameter2', default=2, required=False, dest="param2", metavar="<second parameter>",type=int,help="Second parameter, variance for log-normal and shape parameter for weibull and gamma, integer, default=2")
parser.add_argument('--data', default=None, required=False, dest="data", metavar="<file of stay lengths>",help="If data specified with -D, path to file where each line is length of stay (days)")
parser.add_argument('-N', '--network', default='f', type=str, required=False, dest="network", metavar="<network analysis>",help="perform network analysis (false, default), graph (of directed transmission network), distribution (histogram of out-degrees), parameters (mu and k from neg-binom)")
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
image = str2bool(args.graph)
KM = args.KM.lower()
distribution = args.dist.lower()
data= args.data
prevalence = str2bool(args.prev)
param2 = args.param2
entry_risk = args.entry_risk
network_cmd = args.network.lower()

#Check list of input stay lengths in data
data_list = []
if data != None:
	with open(data, 'r') as input_data:
		for line in input_data:
			num = float(line.split()[0].strip())
			data_list.append(num)

#Create Klebsiella class
class Klebsiella:
	def __init__(self, height, width, n_iterations, entry_rate, median_stay, risk_transmission, image, KM, distribution, data, prevalence, param2, entry_risk, network_cmd):
		self.height = height
		self.width = width
		self.n_iterations = n_iterations
		self.patients = {}
		self.beds = {}
		self.entry_rate = entry_rate
		self.risk_transmission = risk_transmission
		self.image = image
		self.contacts = []
		self.median_stay = median_stay
		self.KM = KM
		self.distribution = distribution
		self.data= data
		self.prevalence = prevalence
                self.param2 = param2
                self.entry_risk = entry_risk
                self.network = network_cmd
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
			
			#If image flag turned on - use plot function
			if self.image ==True:
                                #dict of staus - for visualisation
                                infected_dict = {}
                                infected_dict.update({el:1 for el in self.empty_beds})
                                infected_dict.update({el:2 for el in self.bed_infected})
                                infected_dict.update({el:3 for el in self.bed_uninfected})
				plot(infected_dict, i, self.height, self.width)
			#If prevalence turned on - print values
			if self.prevalence == True:
				print 'time = {}, empty beds = {}, infected beds = {}, uninfected beds = {}'.format(i, len(self.empty_beds), len(self.bed_infected), len(self.bed_uninfected))

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
       
        #Analysis of directed contact network
        def network_run(self):
                DG= nx.DiGraph()
                #patient dictionary is nodes
                DG.add_nodes_from(self.patients)
                #contact list is edges
                DG.add_edges_from(self.contacts)
                #Produce network graph if option selected
                if self.network == "graph":
                        #Get list of nodes with edges (tricky syntax!)
                        node_edge = list(set([j for i, k in DG.edges for j in i, k]))
                        #and make dictionary of labels
                        node_labels = {name: name for name in node_edge}
                        #plot network
                        pos = nx.spring_layout(DG)
                        nx.draw(DG, pos, with_labels=False, nodelist=node_edge)
                        #nx.draw_networkx_nodes(DG, pos, with_labels=False, nodelist=node_edge, alpha=0.8, size=20)
                        nx.draw_networkx_labels(DG, pos, node_labels, font_size=12)
                        plt.show()
                #Produce edge distribution (offspring distribution, or secondary cases)
                else: 
                        out_degree= DG.out_degree()
                        out_degree_array = numpy.array([el[1] for el in out_degree])
                        #print 'mean = {}, max = {}, min = {}'.format(int(numpy.mean(out_degree_array)), max(out_degree_array), min(out_degree_array))
                        
                        if self.network == "distribution":
                                plt.hist(out_degree_array)
                                plt.show()
                        #Maximum likelihood estimator for negative binomial
                        elif self.network == "parameters":
                                #define MLE function
                                def MLE(params):
                                        mu1 = params[0]
                                        k1 = params[1]
                                        p1 = k1/(k1+mu1)
                                        #Calculate log-likelihood as the negative sum of the log probablity-mass-function
                                        logLik = -numpy.sum( stats.nbinom.logpmf(out_degree_array, k1, p1))
                                        return(logLik)

                                #estimate mu and k parameters from negative binomial
                                init_mu = numpy.mean(out_degree_array)
                                #moment estimator for k
                                init_k = numpy.mean(out_degree_array)**2/(numpy.var(out_degree_array)+numpy.mean(out_degree_array))
                                #Run MLE function
                                results = minimize(MLE, [init_mu,init_k], method='L-BFGS-B', bounds=((0, None),(0, None)))
                                #print esimates
                                print 'mu = {}, k = {}'.format(results.x[0], results.x[1])


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

#Plot ward grid each iteration
def plot(infection_dict, i, height, width):
        #plt.interactive(False)
        title = 'Time = {}'.format(i)
        fig, ax = plt.subplots()
        patient_colors = {1:'tab:gray', 2:'tab:red', 3:'tab:green'}
        for patient in infection_dict:
                ax.scatter(patient[0]+0.5, patient[1]+0.5, color=patient_colors[infection_dict[patient]])
        ax.set_title(title, fontsize=10, fontweight='bold')
        ax.set_xlim([0, width])
        ax.set_ylim([0, height])
        ax.set_xticks([])
        ax.set_yticks([])
	#plt.plot()
        plt.show()
	#raw_input('press return to continue')

#Run model - params inherited from command line arguments
run_1 = Klebsiella(height, width, n_iterations, entry_rate, median_stay, risk_transmission, image, KM, distribution, data_list, prevalence, param2, entry_risk, network_cmd)
run_1.populate()
run_1.admit()
if KM in ["plot", "median", "table"]:
        run_1.survival()
if network_cmd in ["graph", "distribution", "parameters"]:
        run_1.network_run()
