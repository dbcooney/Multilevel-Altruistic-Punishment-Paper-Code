"""
Script used to run the numerical simulations used to produce Figures 5.3, 5.4, 5.5, and
5.6, which provide comparisons of the average payoff and average fraction of altruistic
punishers achieved under the multilevel dynamics with pairwise group-level competition.
These figures provide comparisons of these quantities for various fixed and per-interaction
costs of punishment q and p for both the Fermi and local group update rules, and each 
figure will consider the average payoff and fraction of altruistic punishers plotted as 
a function of the strength of punishment p. 
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cmap
import os


from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)



"""
Defining the average payoff of group members in terms of the game-theoretic parameters.
"""

def G(y,b,c,p,q,k):
	return (b - (c+p+q+k)) * y + (p + k) * (y ** 2.0)
	
def G_j(j,N,b,c,p,q,k):
	gamma = b - (c+p+q+k)
	alpha = p+k
	return N * (0.5 * gamma * (2.0 * j + 1.0) * (N ** -2.0) + alpha * (1.0 / (3.0 * (N** 3.0))) * (3.0 * (j ** 2.0)  + 3.0 * j + 1.0) )
	
G_vec = np.vectorize(G)
Gj_vec = np.vectorize(G_j)
	
def G_j_max(N,b,c,p,q,k):
	j_range = np.arange(0.,float(N) + 1.,1.)
	return np.max(Gj_vec(j_range,N,b,c,p,q,k))
	
def G_j_min(N,b,c,p,q,k):
	#j_range = np.arange(0.,float(N) + 1.,1.)
	#return np.min(Gj_vec(j_range,N,b,c,p,q,k))
	if b - c - p - q - k < 0:
		num = -((b - (c + p + q + k)) ** 2.)
		denom = 4 * (p+k)
		return  num / denom
	else:
		return G(0.,b,c,p,q,k)
	

	
	
N = 200
time_step = 0.006
time_length = 9600
#time_length = 20000

#type = "coop"
type = "payoff"
#type = "group local"

quantity = "payoff"

switch_prob = "Fermi"
#switch_prob = "local"
#switch_prob = "Tullock"

b = 2.
#p = 0.6
#k = 5.6
k = 0.
c = 1.
q = 0.5


#alpha = p - k
#beta = -(c + k)

s = 1.
inv_a = 2.
lamb = 2.
lamb = 5.

def theta_init(j,N,theta):
	return N ** (1.0 - theta) * (((N - j) ** theta) - ((N - j - 1.0) ** theta) )
	
theta_vec = np.vectorize(theta_init)


index_holder = np.zeros(N)
for j in range(N):
	index_holder[j] = j
	
f_j = np.ones(N)
f_j = theta_vec(index_holder,N,1.0)

index_vec = np.zeros(N)
for j in range(N):
	index_vec[j] = float(j) / N


"""
Defining fluxes across volume boundaries (corresponding to the effects of individual-
level replication events) and characterizing how these fluxes impact the within-group
replicator dynamics. 
"""


def flux_right(j,N,b,c,p,q,k):
	beta = -(c+k+q)
	alpha = k + p
	return ((j+1.0) / N) * (1.0 - (j+1.0) / N) * (beta + alpha * ((j+1.0)/N))
	
def flux_left(j,N,b,c,p,q,k):
	beta = -(c+k+q)
	alpha = k + p
	return ((np.float(j)) / N) * (1.0 - (np.float(j)) / N) * (beta + alpha * ((np.float(j))/N))
	



flux_right_vec = np.vectorize(flux_right)
flux_left_vec = np.vectorize(flux_left)

def within_group(f,N,b,c,p,q,k,index_holder):

	left_roll = np.roll(f,-1)
	left_roll[-1] = 0.
	right_roll = np.roll(f,1)
	right_roll[0] = 0.
	
	upper_flux = flux_right_vec(index_holder,N,b,c,p,q,k)
	lower_flux = flux_left_vec(index_holder,N,b,c,p,q,k)
	
	upper_flux_up = np.where(upper_flux < 0.0,1.0,0.0)
	upper_flux_down = np.where(upper_flux > 0.0,1.0,0.0)
	
	lower_flux_up = np.where(lower_flux < 0.0,1.0,0.0)
	lower_flux_down = np.where(lower_flux > 0.0,1.0,0.0)
	
	
	upper_half = upper_flux_up * upper_flux * left_roll + upper_flux_down * upper_flux * f 
	lower_half = lower_flux_up * lower_flux * f + lower_flux_down * lower_flux * right_roll
	return N*(-upper_half + lower_half)
	

	
"""
Defining terms used to describe between-group competition. 
"""		

def group_function(y,type,b,c,p,q,k,N):
	
	if type == "coop":
		return y
	elif type == "payoff":
		return G(y,b,c,p,q,k) 
	elif type == "group local":
		G_max = G_j_max(N,b,c,p,q,k)
		G_min = G_j_min(N,b,c,p,q,k)
		return G(y,b,c,p,q,k) / (G_max - G_min)


"""
Defining the victory probability of a focal group featuring a fraction y of altruistic 
punishers when paired against a u-punisher group in a group-level conflict.
""" 				
		
def group_switch_prob(y,u,s,type,b,c,p,q,k,N,G_min):
	
	focal_group = group_function(y,type,b,c,p,q,k,N)
	role_group = group_function(u,type,b,c,p,q,k,N)
	if switch_prob == "Fermi":
		return 0.5 + 0.5 * np.tanh(s * (focal_group - role_group))
	elif switch_prob == "local":
		if focal_group == 0. and role_group == 0.:
			return 0.5
		else:
			return 0.5 + 0.5 * ((focal_group - role_group) / (np.abs(focal_group) + np.abs(role_group)))
	elif switch_prob == "Tullock":
		if focal_group == G_min and role_group == G_min:
			return 0.5
		else:
			num = (focal_group - G_min)**(inv_a)
			denom = (focal_group - G_min)**(inv_a) + (role_group - G_min)**(inv_a)
			if focal_group - G_min < 0 or role_group - G_min < 0:
				print("Negative")
				print(role_group)
				print(G_min)
			return num / denom
	#elif type == "group local":
	#	return 0.5 + 0.5 * (focal_group - role_group)
	#return 0.5 + 0.5 * (focal_group - role_group)



"""
Calculating average group-level victory probability for z-punisher groups over u-punisher
groups for (z,u) \in [i/N,(i+1)/N] \times [j/N,(j+1)/N] using the trapezoidal rule and
our finite volume assumption that the density is a piecewise-constant function taking
constant values on each grid volume.
"""

	
def group_switch_terms(i,j,N,s,type,b,c,p,q,k):
	
	ll = group_switch_prob(float(i)/N,float(j)/N,s,type,b,c,p,q,k,N,G_min)
	lr = group_switch_prob((i+1.)/N,float(j)/N,s,type,b,c,p,q,k,N,G_min)
	ul = group_switch_prob(float(i)/N,(j+1.)/N,s,type,b,c,p,q,k,N,G_min)
	ur = group_switch_prob((i+1.)/N,(j+1.)/N,s,type,b,c,p,q,k,N,G_min)
	
	return 0.25 * (ll + lr + ul + ur) 
	
	
"""
Further characterizing the group-level victory probabilities for each grid volume, and
using these calculations to describe the effect of pairwise group-level competition on
the dynamics of multilevel selection.
"""	
	
def group_switch_matrix(N,s,type,b,c,p,q,k,G_min):
	
	matrix = np.zeros((N,N))
	for i in range(N):
		for j in range(N):
			matrix[i,j] = group_switch_terms(i,j,N,s,type,b,c,p,q,k)
	return matrix
	

#group_matrix = group_switch_matrix(N,s,type,b,c,p,q,k)


def between_group_term(f,N,s,type,group_matrix,b,c,p,q,k,G_min):
	return (1. / N) * f * ( np.dot(group_matrix,f) - np.dot(np.transpose(group_matrix),f))

	


peak_holder = [float(np.argmax(f_j))/N]





Z = [[0,0],[0,0]]
levels = np.arange(0.,time_step * time_length+ time_step,time_step)
CS3 = plt.contourf(Z, levels, cmap=cmap.get_cmap('viridis_r'))
plt.clf()


"""
Defining range of strengths p of altruistic punishment for which we run our finite
volume simulations. We will loop through these values of p, allowing us to describe the
collective outcomes achieved across this range of punishment strengths. 
"""

p_min = 0.
p_max = 3.
p_step = 0.05
p_range = np.arange(p_min,p_max + p_step, p_step)

p_holder_list = []
group_average_holder = []
mean_holder = []


"""
Running the finite volume simulations for our model of multilevel selection with
pairwise group-level competition for a range of strengths p of altruistic punishment.
We use these simulations to generate plots of the densities achieved for each strength p
after 9,600 time-steps.
"""

for p in p_range:

	f_j = theta_vec(index_holder,N,1.0)
	G_min = G_j_min(N,b,c,p,q,k)
	group_matrix = group_switch_matrix(N,s,type,b,c,p,q,k,G_min)
	


	for time in range(time_length):

		#if time % 200 == 0:
			#plt.plot(np.arange(0.5/N,1.0+0.5/N,1.0/N),f_j, color = cmap.viridis_r((np.float(time) / time_length)**1), lw = 3.)
			#plt.plot(np.arange(0.5/N,1.0+0.5/N,1.0/N),f_j, color = cmap.YlOrRd((np.float(time) / time_length)**.25), lw = 3.)
		#elif time in [25,50,100]:
		#plt.plot(np.arange(0.5/N,1.0+0.5/N,1.0/N),f_j, color = cmap.jet((np.float(time) / time_length)**1), lw = 3.)
		#plt.plot(np.arange(0.5/N,1.0+0.5/N,1.0/N),f_j, color = cmap.YlOrRd((np.float(time) / time_length)**.25), lw = 3.)
		between_group_effect = between_group_term(f_j,N,s,type,group_matrix,b,c,p,q,k,G_min)
		within_group_effect = within_group(f_j,N,b,c,p,q,k,index_holder)
		righthandside = lamb * between_group_effect + within_group_effect
		#print(np.sum(righthandside))
		f_j = f_j + time_step * righthandside
		#f_j = f_j * (N / np.sum(f_j))
	
		#print((1.0 / N) * np.sum(f_j))
		
		
	print(p)	
		
		
	if p == 1. or p ==2.:
			#print(Gj_vec(range(N),))	
			plt.plot(np.arange(0.5/N,1.0+0.5/N,1.0/N),f_j, color = cmap.viridis_r((np.float(time) / time_length)**1), lw = 3.)
			print((1.0 / float(N)) * np.dot(f_j,Gj_vec(range(N),N,b,c,p,q,k)))
			
			print(f_j)
			print(Gj_vec(range(N),N,b,c,p,q,k))
	p_holder_list.append(p)
	group_average_holder.append((1.0 / float(N)) * np.dot(f_j,Gj_vec(range(N),N,b,c,p,q,k)))
	mean_holder.append((1.0 / float(N)) * np.dot(f_j,np.arange(0.,1.,1./float(N))))

		


plt.xlabel(r"Fraction of Cooperators ($x$)", fontsize = 20.)
plt.ylabel(r"Probability Density ($f(t,x)$)", fontsize = 20.)
plt.colorbar(CS3) 

plt.tight_layout()





print(np.where(c < 0.0,1.0,0.0))


print(np.dot(f_j,index_vec) / np.sum(f_j))

print(p_holder_list)
print(group_average_holder)
print(Gj_vec(0,N,b,c,p_holder_list,q,k))

plt.figure(2)

plt.plot(p_holder_list,group_average_holder, lw = 5.)
plt.plot(p_holder_list,mean_holder, lw = 5., color = 'Purple')

def G_steady(b,c,p,q,k,theta,lamb):
	if p > c + q:
		return b - c -q
	else:
		return max(b - c - q - (theta / lamb) * (c + q - p),0.)

G_steady_vec = np.vectorize(G_steady)



plt.axis([0.,p_max,-0.1,1.2])


script_folder = os.getcwd()
altruistic_folder = os.path.dirname(script_folder)
file_path = altruistic_folder



if k == 0. and switch_prob == "Fermi" and lamb == 2.:
	if q == 0.5:
		file1 = file_path + "/Simulation_Outputs/GplotDPnonlinearq0p5model.txt"
		file2 = file_path + "/Simulation_Outputs/CoopplotDPnonlinearq0p5model.txt"
	elif q == 0.6:
		file1 = file_path + "/Simulation_Outputs/GplotDPnonlinearq0p6model.txt"
		file2 = file_path + "/Simulation_Outputs/CoopplotDPnonlinearq0p6model.txt"
	elif q == 0.7:
		file1 = file_path + "/Simulation_Outputs/GplotDPnonlinearq0p7model.txt"
		file2 = file_path + "/Simulation_Outputs/CoopplotDPnonlinearq0p7model.txt"
	elif q == 0.8:
		file1 = file_path + "/Simulation_Outputs/GplotDPnonlinearq0p8model.txt"
		file2 = file_path + "/Simulation_Outputs/CoopplotDPnonlinearq0p8model.txt"
	elif q == 0.86:
		file1 = file_path + "/Simulation_Outputs/GplotDPnonlinearq0p86model.txt"
		file2 = file_path + "/Simulation_Outputs/CoopplotDPnonlinearq0p86model.txt"
	elif q == 0.925:
		file1 = file_path + "/Simulation_Outputs/GplotDPnonlinearq0p925model.txt"
		file2 = file_path + "/Simulation_Outputs/CoopplotDPnonlinearq0p925model.txt"

elif k == 0. and switch_prob == "Fermi" and lamb == 5.:	
	if q == 0.5:
		file1 = file_path + "/Simulation_Outputs/GplotDPnonlinearq0p5modellamb5.txt"
		file2 = file_path + "/Simulation_Outputs/CoopplotDPnonlinearq0p5modellamb5.txt"
	elif q == 0.6:
		file1 = file_path + "/Simulation_Outputs/GplotDPnonlinearq0p6modellamb5.txt"
		file2 = file_path + "/Simulation_Outputs/CoopplotDPnonlinearq0p6modellamb5.txt"
	elif q == 0.7:
		file1 = file_path + "/Simulation_Outputs/GplotDPnonlinearq0p7modellamb5.txt"
		file2 = file_path + "/Simulation_Outputs/CoopplotDPnonlinearq0p7modellamb5.txt"
	elif q == 0.8:
		file1 = file_path + "/Simulation_Outputs/GplotDPnonlinearq0p8modellamb5.txt"
		file2 = file_path + "/Simulation_Outputs/CoopplotDPnonlinearq0p8modellamb5.txt"
	elif q == 0.86:
		file1 = file_path + "/Simulation_Outputs/GplotDPnonlinearq0p86modellamb5.txt"
		file2 = file_path + "/Simulation_Outputs/CoopplotDPnonlinearq0p86modellamb5.txt"
	elif q == 0.925:
		file1 = file_path + "/Simulation_Outputs/GplotDPnonlinearq0p925modellamb5.txt"
		file2 = file_path + "/Simulation_Outputs/CoopplotDPnonlinearq0p925modellamb5.txt"
		
elif k == 0. and switch_prob == "Fermi" and lamb == 10.:	
	if q == 0.5:
		file1 = file_path + "/Simulation_Outputs/GplotDPnonlinearq0p5modellamb10.txt"
		file2 = file_path + "/Simulation_Outputs/CoopplotDPnonlinearq0p5modellamb10.txt"
	elif q == 0.6:
		file1 = file_path + "/Simulation_Outputs/GplotDPnonlinearq0p6modellamb10.txt"
		file2 = file_path + "/Simulation_Outputs/CoopplotDPnonlinearq0p6modellamb10.txt"
	elif q == 0.7:
		file1 = file_path + "/Simulation_Outputs/GplotDPnonlinearq0p7modellamb10.txt"
		file2 = file_path + "/Simulation_Outputs/CoopplotDPnonlinearq0p7modellamb10.txt"
	elif q == 0.8:
		file1 = file_path + "/Simulation_Outputs/GplotDPnonlinearq0p8modellamb10.txt"
		file2 = file_path + "/Simulation_Outputs/CoopplotDPnonlinearq0p8modellamb10.txt"
	elif q == 0.86:
		file1 = file_path + "/Simulation_Outputs/GplotDPnonlinearq0p86modellamb10.txt"
		file2 = file_path + "/Simulation_Outputs/CoopplotDPnonlinearq0p86modellamb10.txt"
	elif q == 0.925:
		file1 = file_path + "/Simulation_Outputs/GplotDPnonlinearq0p925modellamb10.txt"
		file2 = file_path + "/Simulation_Outputs/CoopplotDPnonlinearq0p925modellamb10.txt"
				
		
elif q == 0. and switch_prob == "Fermi":
	if k == 0.1:
		file1 = file_path + "/Simulation_Outputs/GplotDPnonlineark0p1model.txt"
		file2 = file_path + "/Simulation_Outputs/CoopplotDPnonlineark0p1model.txt"
	elif k == 0.6:
		file1 = file_path + "/Simulation_Outputs/GplotDPnonlineark0p6model.txt"
		file2 = file_path + "/Simulation_Outputs/CoopplotDPnonlineark0p6model.txt"
	elif k == 1.1:
		file1 = file_path + "/Simulation_Outputs/GplotDPnonlineark1p1model.txt"
		file2 = file_path + "/Simulation_Outputs/CoopplotDPnonlineark1p1model.txt"
	elif k == 1.6:
		file1 = file_path + "/Simulation_Outputs/GplotDPnonlineark1p6model.txt"
		file2 = file_path + "/Simulation_Outputs/CoopplotDPnonlineark1p6model.txt"
	elif k == 2.1:
		file1 = file_path + "/Simulation_Outputs/GplotDPnonlineark2p1model.txt"
		file2 = file_path + "/Simulation_Outputs/CoopplotDPnonlineark2p1model.txt"
	elif k == 2.6:
		file1 = file_path + "/Simulation_Outputs/GplotDPnonlineark2p6model.txt"
		file2 = file_path + "/Simulation_Outputs/CoopplotDPnonlineark2p6model.txt"
	elif k == 3.6:
		file1 = file_path + "/Simulation_Outputs/GplotDPnonlineark3p6model.txt"
		file2 = file_path + "/Simulation_Outputs/CoopplotDPnonlineark3p6model.txt"
	elif k == 4.6:
		file1 = file_path + "/Simulation_Outputs/GplotDPnonlineark4p6model.txt"
		file2 = file_path + "/Simulation_Outputs/CoopplotDPnonlineark4p6model.txt"
	elif k == 5.6:
		file1 = file_path + "/Simulation_Outputs/GplotDPnonlineark5p6model.txt"
		file2 = file_path + "/Simulation_Outputs/CoopplotDPnonlineark5p6model.txt"
	#elif k == 0.925:
	#	file1 = file_path + "/Simulation_Outputs/GplotDPnonlineark0p925model.txt"
	#	file2 = file_path + "/Simulation_Outputs/CoopplotDPnonlineark0p925model.txt"
	
elif k == 0. and switch_prob == "local":
	if q == 0.5:
		file1 = file_path + "/Simulation_Outputs/GplotDPlocalnonlinearq0p5model.txt"
		file2 = file_path + "/Simulation_Outputs/CoopplotDPlocalnonlinearq0p5model.txt"
	elif q == 0.6:
		file1 = file_path + "/Simulation_Outputs/GplotDPlocalnonlinearq0p6model.txt"
		file2 = file_path + "/Simulation_Outputs/CoopplotDPlocalnonlinearq0p6model.txt"
	elif q == 0.7:
		file1 = file_path + "/Simulation_Outputs/GplotDPlocalnonlinearq0p7model.txt"
		file2 = file_path + "/Simulation_Outputs/CoopplotDPlocalnonlinearq0p7model.txt"
	elif q == 0.8:
		file1 = file_path + "/Simulation_Outputs/GplotDPlocalnonlinearq0p8model.txt"
		file2 = file_path + "/Simulation_Outputs/CoopplotDPlocalnonlinearq0p8model.txt"
	elif q == 0.86:
		file1 = file_path + "/Simulation_Outputs/GplotDPlocalnonlinearq0p86model.txt"
		file2 = file_path + "/Simulation_Outputs/CoopplotDPlocalnonlinearq0p86model.txt"
	elif q == 0.925:
		file1 = file_path + "/Simulation_Outputs/GplotDPlocalnonlinearq0p925model.txt"
		file2 = file_path + "/Simulation_Outputs/CoopplotDPlocalnonlinearq0p925model.txt"
		

elif q == 0. and switch_prob == "local":
	if k == 0.1:
		file1 = file_path + "/Simulation_Outputs/GplotDPlocalnonlineark0p1model.txt"
		file2 = file_path + "/Simulation_Outputs/CoopplotDPlocalnonlineark0p1model.txt"
	elif k == 0.6:
		file1 = file_path + "/Simulation_Outputs/GplotDPlocalnonlineark0p6model.txt"
		file2 = file_path + "/Simulation_Outputs/CoopplotDPlocalnonlineark0p6model.txt"
	elif k == 1.1:
		file1 = file_path + "/Simulation_Outputs/GplotDPlocalnonlineark1p1model.txt"
		file2 = file_path + "/Simulation_Outputs/CoopplotDPlocalnonlineark1p1model.txt"
	elif k == 1.6:
		file1 = file_path + "/Simulation_Outputs/GplotDPlocalnonlineark1p6model.txt"
		file2 = file_path + "/Simulation_Outputs/CoopplotDPlocalnonlineark1p6model.txt"
	elif k == 2.1:
		file1 = file_path + "/Simulation_Outputs/GplotDPlocalnonlineark2p1model.txt"
		file2 = file_path + "/Simulation_Outputs/CoopplotDPlocalnonlineark2p1model.txt"
	elif k == 2.6:
		file1 = file_path + "/Simulation_Outputs/GplotDPlocalnonlineark2p6model.txt"
		file2 = file_path + "/Simulation_Outputs/CoopplotDPlocalnonlineark2p6model.txt"
			
				
elif k == 0. and switch_prob == "Tullock" and lamb == 5. and inv_a == 2.:
	if q == 0.5:
		file1 = file_path + "/Simulation_Outputs/GplotDPTullocknonlinearq0p5modellamb5a0p5.txt"
		file2 = file_path + "/Simulation_Outputs/CoopplotDPTullocknonlinearq0p5modellamb5a0p5.txt"
	elif q == 0.6:
		file1 = file_path + "/Simulation_Outputs/GplotDPTullocknonlinearq0p6modellamb55a0p5.txt"
		file2 = file_path + "/Simulation_Outputs/CoopplotDPTullocknonlinearq0p6modellamb55a0p5.txt"
	elif q == 0.7:
		file1 = file_path + "/Simulation_Outputs/GplotDPTullocknonlinearq0p7modellamb55a0p5.txt"
		file2 = file_path + "/Simulation_Outputs/CoopplotDPTullocknonlinearq0p7modellamb55a0p5.txt"
	elif q == 0.8:
		file1 = file_path + "/Simulation_Outputs/GplotDPTullocknonlinearq0p8modellamb55a0p5.txt"
		file2 = file_path + "/Simulation_Outputs/CoopplotDPTullocknonlinearq0p8modellamb55a0p5.txt"
	elif q == 0.86:
		file1 = file_path + "/Simulation_Outputs/GplotDPTullocknonlinearq0p86modellamb55a0p5.txt"
		file2 = file_path + "/Simulation_Outputs/CoopplotDPTullocknonlinearq0p86modellamb55a0p5.txt"
	elif q == 0.925:
		file1 = file_path + "/Simulation_Outputs/GplotDPTullocknonlinearq0p925modellamb55a0p5.txt"
		file2 = file_path + "/Simulation_Outputs/CoopplotDPTullocknonlinearq0p925modellamb55a0p5.txt"
				

else:		
	file1 = file_path + "/Simulation_Outputs/holder_file_1.txt"
	file2 = file_path + "/Simulation_Outputs/holder_file_2.txt"		

		
f1 = open(file1,'w+')
f1.write("strength of punishment p")
f1.write('\n')
f1.write(str(p_holder_list)[1:-1])
f1.write('\n')
f1.write("average payoff G at steady state")
f1.write('\n')
f1.write(str(group_average_holder)[1:-1])
f1.write('\n')
f1.close()
	
f2 = open(file2,'w+')
f2.write("strength of punishment p")
f2.write('\n')
f2.write(str(p_holder_list)[1:-1])
f2.write('\n')
f2.write("average fraction of punishers y at steady state")
f2.write('\n')
f2.write(str(mean_holder)[1:-1])
f2.write('\n')
f2.close()
	

plt.show()


