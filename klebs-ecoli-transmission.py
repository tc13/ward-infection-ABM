#Agent based model to simulate transmission of ESBL Klebsiella pneumoniae and E. coli on a hospital ward
#Tom Crellen, MORU Postdoc. Started January 2019. tomcrellen@gmail.com

#! /usr/bin/env python

import itertools
import random
import numpy
import sys
import argparse

#Argparse
parser = argparse.ArgumentParser(description="Agent Based Model of AMR introduction and spread in a hospital ward. \n \n Author Tom Crellen (tomcrellen@gmail.com) MORU Postdoc")
#import patient length of stay distribution from the command line
parser.add_argument('-l', '--lengthofstay', default=None, required=False, dest="los",metavar="",help="Path to file where each line is length of stay (days)")
parser.add_argument('-i', '--iterations', default=350, required=False, dest="iter", metavar="", help="Number of model iterations (days)")
parser.add_argument('-b', '--beds', default =8, required=False, dest="beds", metavar="", help="Number of beds in ward") 
parser.add_argument('-r', '--entryrate', default =3, required=False, dest="entry", metavar="", help="Entry rate of patients per day (Poisson rate parameter")
parser.add_argument('-k', '--hgtkleb', default=0.05, required=False, dest="HGT_kleb", metavar="", help="Probability of horizontal gene tranfer from K. pneumoniae to E. coli")
parser.add_argument('-p','--masskleb', default=0.05, required=False, dest="PMA_kleb", metavar="", help="Probability of person-to-person transmission of K. pneumoniae")
parser.add_argument('-e', '--hgtecoli', default= 0.0005, required=False, dest="HGT_ecoli", metavar="", help="Probability of horizontal gene transfer from E. coli to K. pneumoniae")
parser.add_argument('-c', '--massecoli', default=0.01, required=False, dest="PMA_ecoli", metavar="", help="Probability of person-to-person transmission of E. coli")
parser.add_argument('-x', '--importkleb', default=0.4, required=False, dest="import_kleb", metavar="", help="Probability that patient is colonized with K. pneumoniae on admission (imported case)")
parser.add_argument('-y', '--importecoli', default=0.3, required=False, dest="import_ecoli", metavar="", help="Probability that patient is colonized with E. coli on admission (imported case)")
args = parser.parse_args()

#Collect arguments passed from the command line
los = args.los
n_iterations = int(args.iter)+1
beds= int(args.beds)
entry_rate = int(args.entry)
HGT_klebs = float(args.HGT_kleb)
PMA_klebs = float(args.PMA_kleb)
HGT_ecoli = float(args.HGT_ecoli)
PMA_ecoli = float(args.PMA_ecoli)

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
        def __init__(self, n_iterations=n_iterations, entry_rate=entry_rate, beds=beds, los_dist=los_dist, HGT_klebs=HGT_klebs, PMA_klebs=PMA_klebs, HGT_ecoli=HGT_ecoli, PMA_ecoli=PMA_ecoli):
                self.n_iterations = n_iterations
                self.los = numpy.array(los_dist)
                self.patients = {}
                self.entry_rate = entry_rate
                self.entry_risk_klebs = 0.4
                self.entry_risk_ecoli = 0.3
                #number of new patients admitted each day, average is rate parameter of poisson
                self.new_patients = numpy.random.poisson(entry_rate, n_iterations)
                self.empty_beds = beds
                self.occupied_beds = []
                self.p_klebs_HGT = HGT_klebs
                self.p_klebs_PMA = PMA_klebs
                self.p_ecoli_HGT = HGT_ecoli
                self.p_ecoli_PMA = PMA_ecoli

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
                                                #patient dict. ID is key, values are day of entry, discharge, klebs dict, ecoli dict
                                                self.patients[name] = [day, discharge_day, {}, {}]
                                                #colonised with Klebsiella on entry and sequence type
                                                klebs_entry = numpy.random.binomial(n=1, p=self.entry_risk_klebs)
                                                if klebs_entry==1:
                                                        klebs_entry_ST = numpy.random.random_integers(1,300)
                                                        self.patients[name][2][klebs_entry_ST] = ["entry", day]
                                                #colonised with E. coli on entry and sequence type
                                                ecoli_entry = numpy.random.binomial(n=1, p=self.entry_risk_ecoli)
                                                if ecoli_entry==1:
                                                        ecoli_entry_ST = numpy.random.random_integers(1,300)
                                                        self.patients[name][3][ecoli_entry_ST] = ["entry", day]
                                                #add patient name to occupied bed list
                                                self.occupied_beds.append(name)
                                        #reduce number of empty beds by number of new patients
                                        self.empty_beds -= new_patients
                       
                                ## TRANSMISSION ##
                                klebs_colonised = []
                                klebs_colonised_ST = []
                                klebs_uncolonised = []
                        
                                ecoli_colonised = []
                                ecoli_colonised_ST = []
                                ecoli_uncolonised = []
                                for key, value in self.patients.iteritems():
                                        #if colonised with any klebs ST and not discharged
                                        if bool(value[2])==True and value[1]>day:
                                                klebs_colonised.append(key)
                                                klebs_colonised_ST.append(numpy.random.choice(value[2].keys()))
                                        #else if uncolonised and present in ward
                                        elif bool(value[2])==False and value[1]>day:
                                                klebs_uncolonised.append(key)
                                        
                                        #if colonised with any E. coli ST and not discharged
                                        if bool(value[3])==True and value[1]>day:
                                                ecoli_colonised.append(key)
                                                ecoli_colonised_ST.append(numpy.random.choice(value[3].keys()))
                                        #else if uncolonised and present in ward
                                        elif bool(value[3])==False and value[1]>day:
                                                ecoli_uncolonised.append(key)

                                #check if any elements in klebs colonised
                                if klebs_colonised:
                                        #WITHIN HOST TRANSMISSION PROCESS (HORIZONTAL GENE TRANSFER - HGT)
                                        klebs_HGT_outcome = numpy.random.binomial(n=1, p=self.p_klebs_HGT, size=len(klebs_colonised))
                                        #indexes of successful transmission events (1)
                                        klebs_HGT_index = [o for o,x in enumerate(klebs_HGT_outcome) if x==1]
                                        #update patient dictionary (klebs -> ecoli)
                                        for j in klebs_HGT_index:
                                                self.patients[klebs_colonised[j]][3][klebs_colonised_ST[j]] = ["HGT", day] 
                        
                                        #BETWEEN HOST TRANSMISSION PROCESS (PSEUDO MASS ACTION PRINCIPAL - PMA)
                                        #check for susceptible patients
                                        if len(klebs_uncolonised) > 0:
                                                #calculate force of infection
                                                klebs_foi = (1-(1-self.p_klebs_PMA))**len(klebs_colonised)
                                                klebs_PMA_outcome = numpy.random.binomial(n=1, p=klebs_foi, size=len(klebs_uncolonised))
                                                klebs_PMA_index = [o for o,x in enumerate(klebs_PMA_outcome) if x==1]
                                                #update patient dictionary with transmission events (klebs -> klebs)
                                                if klebs_PMA_index:
                                                        for j in klebs_PMA_index:
                                                                self.patients[klebs_uncolonised[j]][2][numpy.random.choice(klebs_colonised_ST)] = ["PMA", day]
                                #check if any elements in ecoli colonised
                                if ecoli_colonised:
                                        #WITHIN HOST TRANSMISSION PROCESS (HORIZONTAL GENE TRANSFER - HGT)
                                        ecoli_HGT_outcome = numpy.random.binomial(n=1, p=self.p_ecoli_HGT, size=len(ecoli_colonised))
                                        #indexes of successful transmission events (1)                     
                                        ecoli_HGT_index = [o for o,x in enumerate(ecoli_HGT_outcome) if x==1]
                                        #update patient dictionary (ecoli -> klebs)
                                        if ecoli_HGT_index:
                                                for j in ecoli_HGT_index:
                                                        self.patients[ecoli_colonised[j]][2][ecoli_colonised_ST[j]] = ["HGT", day]

                                        #BETWEEN HOST TRANSMISSION PROCESS (PSEUDO MASS ACTION PRINCIPAL - PMA)
                                        #check for susceptible patients
                                        if ecoli_uncolonised:
                                                #calculate force of infection
                                                ecoli_foi = (1-(1-self.p_ecoli_PMA))**len(ecoli_colonised)
                                                ecoli_PMA_outcome = numpy.random.binomial(n=1, p=ecoli_foi, size=len(ecoli_uncolonised))
                                                ecoli_PMA_index = [o for o,x in enumerate(ecoli_PMA_outcome) if x==1]
                                                #update patient dictionary with transmission events (ecoli -> ecoli)
                                                if ecoli_PMA_index:
                                                        for j in ecoli_PMA_index:
                                                                self.patients[ecoli_uncolonised[j]][3][numpy.random.choice(ecoli_colonised_ST)] = ["PMA", day]


                                #print output to command line
                                print(str(day) + "\t" + str(len(klebs_uncolonised)) + "\t" + str(len(klebs_colonised)) + "\t" + str(len(ecoli_uncolonised)) + "\t" + str(len(ecoli_colonised))) 


#Run model
run=ward()
#print column headers
print("day" + "\t" + "klebsiella_uncolonised" + "\t" + "klebsiella_colonised" + "\t" + "ecoli_uncolonised" + "\t" + "ecoli_colonised")
run.admit()
