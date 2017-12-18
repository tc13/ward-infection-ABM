#Model to simulate ward outbreaks where 1 infectious individual is introduced
#Tom Crellen, MORU Postdoc, tomcrellen@gmail.com

import numpy
import sys
import argparse
import itertools
import random

#Argparse
parser=argparse.ArgumentParser(description="Simulation of single infection in ward \n \n Author Tom Crellen (tomcrellen@gmail.com) MORU Postdoc")
parser.add_argument('-H', '--height', default=10, required=False, dest="height", metavar="<ward height>", type=int,help="height, ward size is height*width, integer, default=10")
parser.add_argument('-W', '--width', default=10, required=False, dest="width", metavar="<ward width>", type=int,help="width, ward size is height*width, integer, default=10")
parser.add_argument('-T', '--time', default=500, required=False, dest="time", metavar="<maximum days>",type=int,help="number of model days (iterations), integer, default=500")
parser.add_argument('-R', '--replicates', default=1000, required=False, dest="replicates", metavar="<model replicates>",type=int,help="number of model replications, integer, default=1000")
parser.add_argument('-TR', '--transmission-risk', default=0.025, required=False, dest="risk", metavar="<risk of transmission>",type=float,help="risk of infection per-person per-day, float, default=0.025")
parser.add_argument('-D', '--distribution', default='log-normal', required=False, dest="dist", metavar="<distribution-stay>",type=str,help="Distribution of stay lengths, log-normal (default), weibull, exponential, gamma, uniform or data (set file path with --data flag)")
parser.add_argument('-A', '--average-stay', default=5, required=False, dest="stay", metavar="<average stay length>",type=int,help="depending on distribution, median (log-normal, uniform), mean (exponential) or scale (gamma, weibull),  integer, default=5")
parser.add_argument('-P2', '--parameter2', default=2, required=False, dest="param2", metavar="<second parameter>",type=int,help="Second parameter, variance for log-normal and shape parameter for weibull and gamma, integer, default=2")
parser.add_argument('--data', default=None, required=False, dest="data", metavar="<file of stay lengths>",help="If data specified with -D, path to file where each line is length of stay (days)")
args = parser.parse_args()

#Command line args- assumes default value if not specified
height = args.height
width = args.width
n_days = args.time
risk = args.risk
distribution = args.dist.lower()
average_stay = args.stay
param2 = args.param2
data = args.data
replicates = args.replicates

#Check list of input stay lengths in data
data_list = []
if data != None:
        with open(data, 'r') as input_data:
                for line in input_data:
                        num = float(line.split()[0].strip())
                        data_list.append(num)

#Ward Class
class R0:
        def __init__(self, height, width, n_days, risk, distribution, average_stay, param2, data_list, replicate):
                self.height = height
                self.width = width
                self.n_days = n_days
                self.risk = risk
                self.patients = {}
                self.beds = {}
                self.stay_distribution = dist(distribution, average_stay, data_list, param2, 100000)
                self.contact_list = []
                self.empty_beds = []
                self.replicate = replicate
                self.bed_infected = []

    #Function to populate the ward with beds of given coordinates (n= width*height)
        def populate(self):
                self.ward = list(itertools.product(range(self.width), range(self.height)))
                #Fill the ward with patients on day 0
                for q in range(1, len(self.ward)):
                        #Give unique ID to patients
                        ID = "0."+str(q)
                        #Sample discharge date from distribution
                        discharge = int(numpy.random.choice(self.stay_distribution))
                        #Patient dictionary- day of entry, length of stay, infect status
                        self.patients[ID] = [0, discharge, 0]
                        #Bed dictionary
                        self.beds[self.ward[q]] = [ID, discharge, 0]

                #Introduce the index case
                discharge = int(numpy.random.choice(self.stay_distribution))
                self.patients["index"] = [0, discharge, 1]
                self.beds[self.ward[0]] = ["index", discharge, 1]
                self.bed_infected.append(self.ward[0])
                self.bed_uninfected = self.ward[1:] 
        
        #Run simulation
        def simulate(self):
                for day in range(self.n_days):
                        if day != 0:
                                n_infectors = len(self.bed_infected)
                                for infector in range(n_infectors):
                                        identity = self.beds[self.bed_infected[infector]][0]
                                        #array of transmission outcomes
                                        transmission_array = numpy.random.binomial(1, self.risk, len(self.bed_uninfected))
                                        transmission_events = [j for j, event in enumerate(transmission_array) if event == 1]
                                        transmission_bed_coords = [self.bed_uninfected[c] for c in transmission_events] 
                                        #For each successful transmission event
                                        for case in transmission_bed_coords:
                                                #Update bed dictionary with new status
                                                self.beds[case][2]=1
                                                id_victim = self.beds[case][0]
                                                #Update patient dictionary with new status
                                                self.patients[id_victim][2]=1
                                                #Add to infected bed list
                                                self.bed_infected.append(case)
                                                #Remove from uninfected bed list
                                                self.bed_uninfected.remove(case)
                                                #Add to contact list
                                                self.contact_list.append([identity, id_victim])

                        #Print output
                        prop_infected =float(len(self.bed_infected))/float(len(self.ward))
                        print '{} {} {}'.format(self.replicate, day, prop_infected)
                        
                        #terminate loop if no more infected patients
                        if len(self.bed_infected) == 0:
                                break

                        #Remove patients with discharge date == date
                        remove = [bed for bed, date in self.beds.items() if date[1] <= day]
                        for k in remove:
                                #remove empty bed from bed dictionary
                                self.beds.pop(k, None)
                                #and other lists
                                try:
                                        self.bed_uninfected.remove(k)
                                except ValueError:
                                        self.bed_infected.remove(k)
                                #add to empty bed list
                                self.empty_beds.append(k)

                        #Admit new patients
                        for spare in range(len(self.empty_beds)):
                                #Give unique ID to patients
                                ID = str(day)+"."+str(spare)
                                #Sample discharge date from distribution
                                discharge = int(numpy.random.choice(self.stay_distribution))
                                #Patient dictionary- day of entry, length of stay, infect status
                                self.patients[ID] = [day, discharge, 0]
                                #Bed dictionary
                                self.beds[self.empty_beds[0]] = [ID, discharge, 0]
                                self.bed_uninfected.append(self.empty_beds[0])
                                self.empty_beds.remove(self.empty_beds[0])

#Set distribution of length of stay
def dist(d, average, data, param2, size):
        if d == "log-normal":
                return list(numpy.clip(numpy.log(numpy.random.lognormal(average, param2, size=size)), 1, None))
        elif d == "exponential":
                return list(numpy.clip(numpy.random.exponential(average, size=size), 1, None))
        elif d == "gamma":
                return list(numpy.clip(numpy.random.gamma(param2, average, size=size), 1, None))
        elif d== "weibull":
                return list(numpy.clip(average*numpy.random.weibull(param2, size=size), 1, None))
        elif d== "uniform":
                return list(numpy.random.uniform(1, (average-0.5)*2, size=size))
        elif d== "data":
                return list(data)
        else:
                raise ValueError("Distribution must be log-normal, gamma, exponential, weibull, uniform or data")

#Run simulation
for rep in range(1, replicates+1):
        name = "run."+str(rep)
        name = R0(height, width, n_days, risk, distribution, average_stay, param2, data_list, rep)
        name.populate()
        name.simulate()
