"""
Script used to generate Figure 5.4, which compares the average payoff and average fraction
of altruistic punishers achieved after 9,600 time-steps for the case of multilevel
competition following the locally normalized group-level pairwise comparison rule and 
fixed cost of punishment (q > 0). This figure makes use of simulation outputs generated 
running the script "loop_altruistic_fvpairwise.py" in the Scripts folder.
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

file1 = altruistic_folder + "/Simulation_Outputs/GplotDPlocalnonlinearq0p5model.txt"
file2 = altruistic_folder + "/Simulation_Outputs/GplotDPlocalnonlinearq0p6model.txt"
file3 = altruistic_folder + "/Simulation_Outputs/GplotDPlocalnonlinearq0p7model.txt"
file4 = altruistic_folder + "/Simulation_Outputs/GplotDPlocalnonlinearq0p8model.txt"
file5 = altruistic_folder + "/Simulation_Outputs/GplotDPlocalnonlinearq0p86model.txt"
file6 = altruistic_folder + "/Simulation_Outputs/GplotDPlocalnonlinearq0p925model.txt"


"""
Functions for processing data from simulation outputs. 
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
Reading in average payoffs achieved after 9,600 from simulation outputs for different
fixed costs q of punishing defectors. 
"""

f1 = open(file1, 'r+')
f1.readline()
p_list = f1.readline()
p_list = process(p_list)


f1.readline()
group_average_1_list = f1.readline()
group_average_1_list = process(group_average_1_list)


f1.close()



f2 = open(file2, 'r+')
group_average_2_list = avg_G_list(f2)
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




"""
Plotting comparison of average payoff for different fixed costs of punishment q, displayed
as a function of strength of punishment p.
"""




plt.figure(1)


cmap_range = 0.86 - 0.5
cmap_min = 0.5


plt.plot(p_list,group_average_5_list, lw = 5., color = plt.cm.viridis((0.86 - cmap_min )/cmap_range ), label = r"$q = 0.86$")
plt.plot(p_list,group_average_4_list, lw = 5., color = plt.cm.viridis((0.8 - cmap_min )/cmap_range ), label = r"$q = 0.8$")
plt.plot(p_list,group_average_3_list, lw = 5., color = plt.cm.viridis((0.7 - cmap_min )/cmap_range ), label = r"$q = 0.7$")
plt.plot(p_list,group_average_2_list, lw = 5., color = plt.cm.viridis((0.6 - cmap_min )/cmap_range ), label = r"$q = 0.6$")
plt.plot(p_list,group_average_1_list, lw = 5., color = plt.cm.viridis((0.5 - cmap_min )/cmap_range ), label = r"$q = 0.5$")



plt.axis([0.,3.,-0.005,0.8])
plt.legend(loc = "upper left", fontsize = 14., ncol = 3)

plt.xlabel(r"Strength of Punishment $p$", fontsize = 20., labelpad = 10.)
plt.ylabel(r"Average Group Payoff $\langle G(\cdot) \rangle_f$", fontsize = 18.)


plt.tight_layout()

plt.savefig(altruistic_folder + "/Figures/avgGplotlocalqmodelDP.png" )





"""
Reading in average fraction of altruistic punishers achieved after 9,600 from simulation 
outputs for different fixed costs q of punishing defectors. 
"""


file1 = altruistic_folder + "/Simulation_Outputs/CoopplotDPlocalnonlinearq0p5model.txt"
file2 = altruistic_folder + "/Simulation_Outputs/CoopplotDPlocalnonlinearq0p6model.txt"
file3 = altruistic_folder + "/Simulation_Outputs/CoopplotDPlocalnonlinearq0p7model.txt"
file4 = altruistic_folder + "/Simulation_Outputs/CoopplotDPlocalnonlinearq0p8model.txt"
file5 = altruistic_folder + "/Simulation_Outputs/CoopplotDPlocalnonlinearq0p86model.txt"



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
Plotting comparison of average fraction of altruistic punishment achieved for different 
fixed costs of punishment q, displayed as a function of strength of punishment p.
"""

plt.figure(2)

cmap_range = 0.86 - 0.5
cmap_min = 0.5


plt.plot(p_list,avg_punisher_5_list, lw = 3., color = plt.cm.viridis((0.86 - cmap_min )/cmap_range ), label = r"$q = 0.86$")
plt.plot(p_list,avg_punisher_4_list, lw = 3., color = plt.cm.viridis((0.8 - cmap_min )/cmap_range ), label = r"$q = 0.8$")
plt.plot(p_list,avg_punisher_3_list, lw = 3., color = plt.cm.viridis((0.7 - cmap_min )/cmap_range ), label = r"$q = 0.7$")
plt.plot(p_list,avg_punisher_2_list, lw = 3., color = plt.cm.viridis((0.6 - cmap_min )/cmap_range ), label = r"$q = 0.6$")
plt.plot(p_list,avg_punisher_1_list, lw = 3., color = plt.cm.viridis((0.5 - cmap_min )/cmap_range ), label = r"$q = 0.5$")



plt.axis([0.,3.,0.4,1.05])
plt.legend(loc = "lower right", fontsize = 14., ncol = 3)

plt.xlabel(r"Strength of Punishment $p$", fontsize = 20., labelpad = 10.)
plt.ylabel(r"Fraction of Punishers $\langle z \rangle_f$", fontsize = 18.)


plt.tight_layout()

plt.savefig(altruistic_folder + "/Figures/avgpunisherplotlocalqmodelDP.png" )


plt.show()