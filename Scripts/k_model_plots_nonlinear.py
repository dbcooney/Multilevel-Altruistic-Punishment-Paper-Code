"""
Script used to generate Figure 5.4, which illustrates the average payoff and average
fraction of altruistic punishers achieved after many time-steps in our model of multilevel
selection with pairwise group-level competition featuring per-interaction punishment costs
and the Fermi group-level victory probability. We plot the average payoffs and and
average fraction of altruistic punishers as a function of the punishment strength p and
compare cases of different values of the per-interaction punishment cost k.


This figure takes files from the Simulation_Outputs folder to produce this figure, and
the simulation outputs used for this figure were generated by running the script
"loop_altruistic_fvpairwise.py" in the Scripts folder.
"""


import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')
import os


from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

#eta = 1.0

plt.figure(1)

script_folder = os.getcwd()
altruistic_folder = os.path.dirname(script_folder)

file1 = altruistic_folder + "/Simulation_Outputs/GplotDPnonlineark0p6model.txt"
file2 = altruistic_folder + "/Simulation_Outputs/GplotDPnonlineark1p6model.txt"
file3 = altruistic_folder + "/Simulation_Outputs/GplotDPnonlineark2p6model.txt"
file4 = altruistic_folder + "/Simulation_Outputs/GplotDPnonlineark3p6model.txt"
file5 = altruistic_folder + "/Simulation_Outputs/GplotDPnonlineark4p6model.txt"
file6 = altruistic_folder + "/Simulation_Outputs/GplotDPnonlineark5p6model.txt"
file7 = altruistic_folder + "/Simulation_Outputs/GplotDPnonlinearkp1model.txt"
file8 = altruistic_folder + "/Simulation_Outputs/GplotDPnonlineark3p6model.txt"


"""
Functions used to read in files from the Simulation_Outputs folder. 
"""

def process(list):
	list = list.split(',')
	list = [float(a) for a in list]
	return list
	
def avg_G_list(f):
	f.readline()
	f.readline()
	f.readline()
	group_payoff_list = f.readline()
	group_payoff_list = process(group_payoff_list)
	return group_payoff_list
	
def avg_punisher_list(f):
	f.readline()
	f.readline()
	f.readline()
	avg_punisher_list = f.readline()
	avg_punisher_list = process(avg_punisher_list)
	return avg_punisher_list

"""
Reading in the simulation outputs we will use to compare average payoff for different
per-interaction costs of punishment k.
"""


f1 = open(file1, 'r+')
f1.readline()
p_list1 = f1.readline()
p_list1 = process(p_list1)


f1.readline()
group_average_1_list = f1.readline()
group_average_1_list = process(group_average_1_list)


f1.close()



f2 = open(file2, 'r+')
f2.readline()
p_list = f2.readline()
p_list = process(p_list)
f2.readline()
group_average_2_list = f2.readline()
group_average_2_list = process(group_average_2_list)
f2.close()



f3 = open(file3, 'r+')
group_average_3_list = avg_G_list(f3)
f3.close()



f4 = open(file4, 'r+')
group_average_4_list = avg_G_list(f4)
f4.close()



f5 = open(file5, 'r+')
group_average_5_list = avg_G_list(f5)
f5.close()

f6 = open(file6, 'r+')
group_average_6_list = avg_G_list(f6)
f6.close()

#f7 = open(file7, 'r+')
#group_average_7_list = avg_G_list(f7)
#f7.close()

#f8 = open(file6, 'r+')
#group_average_6_list = avg_G_list(f6)
#f6.close()




"""
Plotting average payoff achieved after 9,600 time-steps as a function of the punishment 
strength p, providing a comparison across different values of the per-interaction
punishment cost k. 
"""

plt.figure(1)



cmap_range = 5.6 - 0.6
cmap_min = 0.6

plt.plot(p_list,group_average_6_list, lw = 5., color = plt.cm.viridis((5.6 - cmap_min )/cmap_range ), label = r"$k = 5.6$")
plt.plot(p_list,group_average_5_list, lw = 4., color = plt.cm.viridis((4.6 - cmap_min )/cmap_range ), label = r"$k = 4.1$")
plt.plot(p_list,group_average_4_list, lw = 4., color = plt.cm.viridis((3.6 - cmap_min )/cmap_range ), label = r"$k = 3.6$")
plt.plot(p_list,group_average_3_list, lw = 4., color = plt.cm.viridis((2.6 - cmap_min )/cmap_range ), label = r"$k = 2.1$")
plt.plot(p_list,group_average_2_list, lw = 4., color = plt.cm.viridis((1.6 - cmap_min )/cmap_range ), label = r"$k = 1.6$")
plt.plot(p_list1,group_average_1_list, lw = 4., color = plt.cm.viridis((0.6 - cmap_min )/cmap_range ), label = r"$k = 0.6$")



plt.axis([0.,1.5,-0.005,1.1])
plt.legend(loc = "lower left", fontsize = 14., ncol = 3)

plt.xlabel(r"Strength of Punishment $p$", fontsize = 20., labelpad = 10.)
plt.ylabel(r"Average Group Payoff $\langle G(\cdot) \rangle_f$", fontsize = 18.)


plt.tight_layout()

plt.savefig(altruistic_folder + "/Figures/avgGplotkmodelDP.png" )






file1 = altruistic_folder + "/Simulation_Outputs/CoopplotDPnonlineark0p6model.txt"
file2 = altruistic_folder + "/Simulation_Outputs/CoopplotDPnonlineark1p6model.txt"
file3 = altruistic_folder + "/Simulation_Outputs/CoopplotDPnonlineark2p6model.txt"
file4 = altruistic_folder + "/Simulation_Outputs/CoopplotDPnonlineark3p6model.txt"
file5 = altruistic_folder + "/Simulation_Outputs/CoopplotDPnonlineark4p6model.txt"
file6 = altruistic_folder + "/Simulation_Outputs/CoopplotDPnonlineark5p6model.txt"




f1 = open(file1, 'r+')
avg_punisher_1_list = avg_punisher_list(f1)
f1.close()


f2 = open(file2, 'r+')
avg_punisher_2_list = avg_punisher_list(f2)
f2.close()


f3 = open(file3, 'r+')
avg_punisher_3_list = avg_punisher_list(f3)
f3.close()


f4 = open(file4, 'r+')
avg_punisher_4_list = avg_punisher_list(f4)
f4.close()


f5 = open(file5, 'r+')
avg_punisher_5_list = avg_punisher_list(f5)
f5.close()


f6 = open(file6, 'r+')
avg_punisher_6_list = avg_punisher_list(f6)
f6.close()





"""
Plotting the average fraction of altruistic punishment achieved in the population after
9,600 time-steps as a function of the punishment strength p, providing a comparison across 
different values of the per-interaction punishment cost k. 
"""

plt.figure(2)

cmap_range = 5.6 - 0.6
cmap_min = 0.6

plt.plot(p_list,avg_punisher_6_list, lw = 4., color = plt.cm.viridis((5.6 - cmap_min )/cmap_range ), label = r"$k = 5.6$")
plt.plot(p_list,avg_punisher_5_list, lw = 4., color = plt.cm.viridis((4.6 - cmap_min )/cmap_range ), label = r"$k = 4.6$")
plt.plot(p_list,avg_punisher_4_list, lw = 4., color = plt.cm.viridis((3.6 - cmap_min )/cmap_range ), label = r"$k = 3.6$")
plt.plot(p_list,avg_punisher_3_list, lw = 4., color = plt.cm.viridis((2.6 - cmap_min )/cmap_range ), label = r"$k = 2.6$")
plt.plot(p_list,avg_punisher_2_list, lw = 4., color = plt.cm.viridis((1.6 - cmap_min )/cmap_range ), label = r"$k = 1.6$")
plt.plot(p_list1,avg_punisher_1_list, lw = 4., color = plt.cm.viridis((0.6 - cmap_min )/cmap_range ), label = r"$k = 0.6$")



plt.axis([0.,1.5,-0.005,1.1])
plt.legend(loc = "lower left", fontsize = 14., ncol = 3)

plt.xlabel(r"Strength of Punishment $p$", fontsize = 20., labelpad = 10.)
plt.ylabel(r"Fraction of Punishers $\langle z \rangle_f$", fontsize = 18.)


plt.tight_layout()

plt.savefig(altruistic_folder + "/Figures/avgpunisherplotkmodelDP.png" )


plt.show()