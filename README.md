# ward-infection-ABM
Reproducible code and data from "Transmission dynamics and control of multidrug-resistant Klebsiella pneumoniae in neonates in a developing country" by Crellen, T., et al.

Python scripts to run agent based models to estimate the basic ward reproduction number (RA) and to simulate the impact of interventions. Check scripts to ensure required libraries are installed prior to running. Options for each model are given by running `python file_name.py -h`.

To use values for colonisation pressure estimated by Bayesian inference in the paper, see the `parameters` folder. 

To replicate the RA simulation in the paper, call from bash (or terminal in MacOS):

`while read A; do python RA_simulation.py -H 4 -W 2 -R 100 -TR ${A} -D data --data parameters/neonates.los.NU.txt | awk '{ sum += $8 } END { if (NR > 0) print sum / NR }'; done < parameters/FOI.posterior.txt > results.txt`

The intervention_simulation.py script can read in two sets of values for colonisation pressure (with options -t0 and -t1), in the form of a tab seperated file. The probability of an individual in the simuations being assisgned colonisation pressure values from -t1 is given by -p. 

For instance, to simulate the impact of breast feeding rates on the number of individuals remaining uncolonised, where 25% of infants in the simulation are breast fed: 

`python intervention_simulation.py -b 9 -e 3 -t0 0.15 -t1 0.10 -p 0.25 -r 100 -x 0.05 `

In the above line, the user has specified the daily probability of acquisition per colonised patient as 0.15 for non-breast fed individuals and 0.10 for breast fed individuals.

To use the values estimated from the model (n=2000) in a loop, where each of the 100 replicates per parameter values are averaged, the script is called by bash (or using Terminal in MacOS) as:

`while read A B; do python intervention_simulation.py -b 9 -e 3 -t0 ${A} -t1 ${B} -p 0.25 -r 100 -x 0.05 | awk '{ sum += $8 } END { if (NR > 0) print sum / NR }'; done < parameters/breast.milk.intervention.txt > results.txt`

Options such as -x (the proportion of infants who are colonised on first admission) can be altered to replicate the analysis for Figure 4 in the text.

Note that the default in intervention_simulation.py is to use the empirical length of stay distribution observed in the study, however the user can specify a different distribution in the form of a file where each LOS values is an integer on a seperate line with the -l option.

For any comments on this code, please contact me on thomas.crellen@ndm.ox.ac.uk or tomcrellen@gmail.com. The code is my own, the original dataset is the property of Prof Ben Cooper, Prof Paul Turner and Dr Claudia Turner.
