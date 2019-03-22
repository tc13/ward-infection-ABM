#Agent based model to simulate transmission of ESBL Klebsiella pneumoniae and E. coli on a hospital ward
#Tom Crellen, MORU Postdoc. Started January 2019. tomcrellen@gmail.com

#! /usr/bin/env python

import itertools
import random
import numpy
import sys
import argparse

#Argparse
parser = argparse.ArgumentParser(description="Agent Based Models of AMR introduction and spread in a hospital ward. \n\nAuthor Tom Crellen (tomcrellen@gmail.com) MORU Postdoc. \n \n Model permits interventions (non-time varying)")
#import patient length of stay distribution from the command line
parser.add_argument('-l', '--lengthofstay', default=None, required=False, dest="los",metavar="",help="Path to file where each line is length of stay in days (empirical distribution)")
parser.add_argument('-i', '--iterations', default=365, required=False, dest="iter", metavar="", help="Number of model iterations / days (365)")
parser.add_argument('-b', '--beds', default =8, required=False, dest="beds", metavar="", help="Number of beds in ward (8)") 
parser.add_argument('-e', '--entryrate', default =3, required=False, dest="entry", metavar="", help="Entry rate of patients per day, Poisson rate parameter (3)")
parser.add_argument('-t0','--trans0', default=0.02, required=False, dest="trans0", metavar="", help="Probability of person-to-person transmission of K. pneumoniae in group 0 (0.02)")
parser.add_argument('-t1','--trans1', default=0.07, required=False, dest="trans1", metavar="", help="Probability of person-to-person transmission of K. pneumoniae in group 1 (0.07)")
parser.add_argument('-p','--prob', default=0.5, required=False, dest="prob_intervention", metavar="", help="Probability that patient is assigned to group 1")
parser.add_argument('-x', '--importkleb', default=0.4, required=False, dest="import_kleb", metavar="", help="Probability that patient is colonized with K. pneumoniae on admission (imported case) (0.4)")
parser.add_argument('-r', '--replicates', default=1, required=False, dest="replicates", metavar="", help="number of model runs")

args = parser.parse_args()

#Collect arguments passed from the command line
los = args.los
n_iterations = int(args.iter)+1
beds= int(args.beds)
entry_rate = int(args.entry)
trans0 = float(args.trans0)
trans1 = float(args.trans1)
prob_intervention = float(args.prob_intervention)
import_klebs = float(args.import_kleb)
model_runs = int(args.replicates)

#Process input lengths of stay
los_dist = []
if los != None:
        with open(los, 'r') as input_los:
	        for line in input_los:
			num = int(line.split()[0].strip())
			los_dist.append(num)
else:
        #if los not specified, use default (333 infants in Cambodian neonatal unit study)
        los_dist = [3, 4, 3, 5, 12, 29, 12, 4, 6, 5, 22, 4, 5, 16, 11, 9, 4, 5, 5, 5, 6, 4, 10, 66, 4, 6, 4, 8, 4, 12, 14, 3, 5, 5, 8, 10, 9, 8, 16, 38, 3, 5, 47, 15, 9, 3, 3, 5, 7, 7, 9, 4, 7, 4, 5, 3, 2, 3, 3, 9, 11, 28, 21, 7, 4, 17, 8, 5, 6, 5, 4, 4, 1, 6, 20, 13, 11, 7, 8, 19, 5, 22, 8, 18, 6, 9, 5, 4, 6, 6, 19, 17, 5, 3, 11, 26, 3, 12, 7, 7, 11, 8, 21, 6, 8, 4, 4, 31, 11, 3, 6, 14, 10, 3, 11, 6, 12, 5, 14, 6, 5, 5, 7, 3, 6, 3, 3, 6, 8, 2, 4, 10, 6, 11, 51, 11, 2, 11, 3, 15, 4, 56, 8, 3, 4, 27, 3, 8, 18, 3, 10, 7, 19, 6, 3, 3, 5, 16, 8, 4, 16, 5, 58, 3, 3, 2, 34, 13, 4, 3, 8, 2, 5, 9, 10, 3, 4, 4, 19, 6, 8, 8, 7, 8, 10, 3, 8, 1, 14, 2, 5, 8, 7, 3, 7, 9, 5, 3, 3, 3, 2, 2, 43, 8, 4, 40, 7, 4, 3, 60, 7, 9, 3, 3, 10, 6, 2, 9, 4, 8, 4, 4, 2, 2, 3, 4, 5, 5, 5, 32, 11, 3, 8, 4, 3, 2, 3, 5, 9, 3, 6, 4, 5, 25, 7, 6, 5, 20, 4, 5, 3, 54, 6, 32, 20, 6, 4, 6, 3, 7, 3, 6, 4, 4, 20, 17, 16, 3, 12, 27, 31, 5, 48, 5, 3, 3, 10, 6, 6, 5, 4, 8, 37, 8, 3, 8, 7, 4, 4, 3, 10, 20, 3, 3, 10, 4, 5, 20, 3, 29, 5, 3, 2, 15, 7, 25, 3, 30, 42, 21, 57, 41, 3, 3, 5, 13, 5, 5, 20, 5, 34, 4, 4, 4, 6, 8, 27, 14, 5, 5, 54, 34, 22]

class ward:
        def __init__(self, n_iterations=n_iterations, entry_rate=entry_rate, beds=beds, los_dist=los_dist, trans0=trans0, trans1=trans1, prob_intervention=prob_intervention, import_klebs=import_klebs):
                self.n_iterations = n_iterations
                self.los = numpy.array(los_dist)
                self.patients = {}
                self.entry_rate = entry_rate
                self.entry_risk_klebs = import_klebs
                #number of new patients admitted each day, average is rate parameter of poisson
                self.new_patients = numpy.random.poisson(entry_rate, n_iterations)
                self.empty_beds = beds
                self.occupied_beds = []
                self.trans0 = trans0
                self.trans1 = trans1
                self.p_group = prob_intervention
                #simulation outcome variables
                self.uncolon_entry_0 = 0
                self.colon_exit_0 = 0
                self.uncolon_entry_1 = 0
                self.colon_exit_1 = 0

        def admit(self):
                #For each day (iteration)
                for day in range(self.n_iterations):
                        #after day zero
                        if day >= 1:
                                #remove patients where discharge day == current day
                                remove = [bed for bed, date in self.patients.items() if date[1]==day]
                                self.occupied_beds = [bed for bed in self.occupied_beds if bed not in remove]
                                self.empty_beds += len(remove)
                                #admin n new patients
                                new_patients = self.new_patients[day]
                                #check there are enough empty beds for new patients
                                if new_patients > self.empty_beds:
                                        new_patients = self.empty_beds
                                #admit patients, if at least one spare bed
                                if new_patients > 0:
                                        for n in range(1, new_patients+1):
                                                #give characteristics to new patients
                                                name = str(day)+"."+str(n)
                                                los = numpy.random.choice(self.los)
                                                discharge_day = day+los
                                                #patient dict. ID is key, values are day of entry, discharge, klebs dict, intervention group
                                                self.patients[name] = [day, discharge_day, {}, 0]
                                                #colonised with Klebsiella on entry and sequence type
                                                klebs_entry = numpy.random.binomial(n=1, p=self.entry_risk_klebs)
                                                if klebs_entry==1:
                                                        klebs_entry_ST = numpy.random.random_integers(1,300)
                                                        self.patients[name][2][klebs_entry_ST] = ["entry", day]
                                                #intervention group
                                                group = numpy.random.binomial(n=1, p=self.p_group)
                                                self.patients[name][3] = group
                                                #add patient to relevant variable
                                                if group==0 and klebs_entry==0:
                                                        self.uncolon_entry_0 += 1
                                                elif group==1 and klebs_entry==0:
                                                        self.uncolon_entry_1 += 1
                                                #add patient name to occupied bed list
                                                self.occupied_beds.append(name)
                                        #reduce number of empty beds by number of new patients
                                        self.empty_beds -= new_patients
                       
                                ## TRANSMISSION ##
                                klebs_colonised = []
                                klebs_colonised_ST = []
                                klebs_uncolon_0 = []
                                klebs_uncolon_1 = []
                        
                                for key, value in self.patients.iteritems():
                                        #if colonised with any klebs ST and not discharged
                                        if bool(value[2])==True and value[1]>day:
                                                klebs_colonised.append(key)
                                                klebs_colonised_ST.append(numpy.random.choice(value[2].keys()))
                                        #else if uncolonised and present in ward (group 0)
                                        elif bool(value[2])==False and value[1]>day and value[3]==0:
                                                klebs_uncolon_0.append(key)
                                        #else if uncolonised and present in ward (group 1)
                                        elif bool(value[2])==False and value[1]>day and value[3]==1:
                                                klebs_uncolon_1.append(key)

                                #check if any elements in klebs colonised
                                if klebs_colonised:
                                        #BETWEEN HOST TRANSMISSION PROCESS (PSEUDO MASS ACTION PRINCIPAL - PMA)
                                        #check for susceptible patients
                                        if len(klebs_uncolon_0) > 0:
                                                #calculate force of infection (group 0)
                                                klebs_foi_0 = (1-(1-self.trans0)**len(klebs_colonised)) 
                                                #binomial random outcome
                                                klebs_PMA_outcome_0 = numpy.random.binomial(n=1, p=klebs_foi_0, size=len(klebs_uncolon_0))
                                                klebs_PMA_index_0 = [o for o,x in enumerate(klebs_PMA_outcome_0) if x==1]
                                                #update patient dictionary with transmission events (klebs -> klebs)
                                                if klebs_PMA_index_0:
                                                        for j in klebs_PMA_index_0:
                                                                self.patients[klebs_uncolon_0[j]][2][numpy.random.choice(klebs_colonised_ST)] = ["PMA", day]
                                                        #update outcome variable
                                                        self.colon_exit_0 += len(klebs_PMA_index_0)
                                        #group 1
                                        if len(klebs_uncolon_1) > 0:
                                               #foi for group 1
                                                klebs_foi_1 = (1-(1-self.trans1)**len(klebs_colonised)) 
                                                klebs_PMA_outcome_1 = numpy.random.binomial(n=1, p=klebs_foi_1, size=len(klebs_uncolon_1))
                                                klebs_PMA_index_1 = [o for o,x in enumerate(klebs_PMA_outcome_1) if x==1]
                                                #update patient dictionary with transmission events (klebs -> klebs)
                                                if klebs_PMA_index_1:
                                                        for j in klebs_PMA_index_1:
                                                                self.patients[klebs_uncolon_1[j]][2][numpy.random.choice(klebs_colonised_ST)] = ["PMA", day]
                                                        #update outcome variable
                                                        self.colon_exit_1 += len(klebs_PMA_index_1)

                #print output to command line
                print(str(len(self.patients)) + "\t" + str(self.uncolon_entry_0) + "\t" + str(self.colon_exit_0)  + "\t" + str(self.uncolon_entry_1) + "\t" + str(self.colon_exit_1) + "\t" + str(self.uncolon_entry_0+self.uncolon_entry_1) + "\t" + str(self.colon_exit_0+self.colon_exit_1) + "\t" + str(float(self.colon_exit_0+self.colon_exit_1)/float(self.uncolon_entry_0+self.uncolon_entry_1)))
                        

#print column headers
print("total_patients" + "\t" + "uncolon_entry_0" + "\t" + "acquired_exit_0" + "\t" + "uncolon_entry_1" + "\t" + "acquired_exit_1" + "\t" + "uncolon_entry_total" + "\t" + "acquired_exit_total" + "\t" + "proportion_acquired_total")
#run model
for i in range(model_runs):
        run=ward()
        run.admit()
