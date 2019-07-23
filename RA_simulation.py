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
parser.add_argument('-T', '--time', default=300, required=False, dest="time", metavar="<maximum days>",type=int,help="number of model days (iterations), integer, default=300")
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
                self.replicate = replicate
                self.bed_infected = []
                self.transmission = []

    #Function to populate the ward with beds of given coordinates (n= width*height)
        def populate(self):
                self.ward = list(itertools.product(range(self.width), range(self.height)))               
                #Give unique ID to patients
                ID = ["0."+str(q) for q in range(len(self.ward))]
                #Sample discharge date from distribution
                discharge = numpy.random.choice(self.stay_distribution, size=len(self.ward))       
                #Bed dictionary
                for q in range(1, len(self.ward)):
                        self.beds[self.ward[q]] = [ID[q], int(discharge[q])]
                #Introduce the index case             
                self.beds[self.ward[0]] = ["index", int(discharge[0])]
                self.bed_infected.append(self.ward[0])
                self.bed_uninfected = self.ward[1:] 
        
        #Run simulation
        def simulate(self):
                for day in range(self.n_days):
                        if day != 0:
                                n_infectors = len(self.bed_infected)
                                #For each uninfected, geometric probability of infection
                                transmission_array = numpy.random.geometric(self.risk, size=len(self.bed_uninfected))
                                #keep if success is <= number of infectors
                                transmission_events = [j for j, event in enumerate(transmission_array) if event <= n_infectors] 
                                #get bed coordinates of successful transmission events
                                transmission_bed_coords = [self.bed_uninfected[c] for c in transmission_events]
                                #get IDs of infected patients from the bed dict
                                transmission_bed_IDs = [self.beds[d][0] for d in transmission_bed_coords]
                                #get bed_coords of infectors
                                infector_bed_coords = [self.bed_infected[f-1] for f in transmission_array if f <= n_infectors]
                                #get ID of bed-infectors
                                infector_ID = [self.beds[g][0] for g in infector_bed_coords]
                                #Zip lists - to give transmission/ contact pairs
                                self.transmission.append(zip(infector_ID, transmission_bed_IDs))
                                #Add to bed infected list
                                self.bed_infected = self.bed_infected +list(transmission_bed_coords)
                                #Remove from uninfected bed list
                                self.bed_uninfected = [x for x in self.bed_uninfected if x not in transmission_bed_coords]
                                               
                        #Print output
                        prop_infected =float(len(self.bed_infected))/float(len(self.ward))
                        print '{} {} {}'.format(self.replicate, day, prop_infected)

                        #terminate loop if no more infected patients
                        if len(self.bed_infected) == 0:
                                break

                        #Remove patients with discharge date == date
                        remove = [bed for bed, date in self.beds.items() if date[1] <= day]
                        self.bed_infected = [bed for bed in self.bed_infected if bed not in remove]
                        self.bed_uninfected = [bed for bed in self.bed_uninfected if bed not in remove]
                       
                        #Admit new patients (equal to number of spare beds)                     
                        spares = len(remove)
                        #Give unique ID to patients
                        ID = [str(day)+"."+str(spare) for spare in range(spares)]
                        #Sample discharge date from distribution
                        discharge = numpy.random.choice(self.stay_distribution, size= spares)
                        #Update bed dictionary with new patients
                        for bed in range(spares):
                                self.beds[remove[bed]] = [ID[bed], int(discharge[bed])+day]
                        #add to uninfected list and remove empty beds
                        self.bed_uninfected = self.bed_uninfected + remove

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
