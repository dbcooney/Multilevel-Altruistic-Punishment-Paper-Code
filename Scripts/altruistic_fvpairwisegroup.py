"""
Script used to generate Figure 5.1, which provides snapshots of the solutions to the
PDE model for multilevel selection with pairwise group-level competition. The figure
provides a comparison of the trajectories for the multilevel model for different relative
strengths of group-level competition \lambda, highlighting cases in which defectors take
over the population and a case in which the population appears to converge to a steady 
state supporting groups with various levels of altruistic punishment. 
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
	return N * (0.5 * gamma * (2.0 * j + 1.0) * (N ** -2.0) + alpha * (1.0 / (3.0 * (N** 3.0))) * (3.0 * (j ** 2.0)  + 3.0 * j + 1.0))
	
G_vec = np.vectorize(G)
Gj_vec = np.vectorize(G_j)
	
def G_j_max(N,b,c,p,q,k):
	j_range = np.arange(0.,float(N) + 1.,1.)
	return np.max(Gj_vec(j_range,N,b,c,p,q,k))
	
def G_j_min(N,b,c,p,q,k):
	j_range = np.arange(0.,float(N) + 1.,1.)
	return np.min(Gj_vec(j_range,N,b,c,p,q,k))
	
	
N = 200
time_step = 0.003



"""
There are several options for setting the group-level victory probability. The variable
"type" allows for a choice of whether group-level victory probability depends on the
average payoff of group members. The choice "coop" allows the probability of victory of 
a z-punisher group over a u-punisher group to depend on the difference in the fraction of
altruistic punisher z - u between the two groups. The choice "payoff" results in 
group-level victory probabilities \rho(z,u) that depends on the difference of average
payoffs G(z) - G(u) between the two groups.

The variable "switch-prob" allows a choice for the form of the group-level victory
probability. The choice "Fermi" allows for the use of the Fermi update rule (defined in
Equation 3.12), while the choice "local" allows for the use of the local group-level 
pairwise comparison with local normalization of payoff differences (defined in Equation
3.13).
"""


#type = "coop"
type = "payoff"

#switch_prob = "Fermi"
switch_prob = "local"


"""
Defining payoff parameters for the model and the parameters determining the dynamics of
multilevel selection for the simulation. 
"""

b = 2.
p = 0.5
k = 0.1
c = 1.
#q = 0.925
q = 0.



s = 2.
lamb = 2.


if lamb == 0.1:
	time_length = 800
	time_snapshot = 150
elif lamb == 2.:
	time_length = 2400
	time_snapshot = 250
else:
	time_length = 9600
	time_snapshot = 200

def theta_init(j,N,theta):
	return N ** (1.0 - theta) * (((N - j) ** theta) - ((N - j - 1.0) ** theta) )
	
theta_vec = np.vectorize(theta_init)


index_holder = np.zeros(N)
x_vec = np.zeros(N)
for j in range(N):
	index_holder[j] = j
	x_vec[j] = np.float64(j) / N
	
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
		
"""
Defining the victory probability of a focal group featuring a fraction y of altruistic 
punishers when paired against a u-punisher group in a group-level conflict.
"""
		
def group_switch_prob(y,u,s,type,b,c,p,q,k,N):
	
	focal_group = group_function(y,type,b,c,p,q,k,N)
	role_group = group_function(u,type,b,c,p,q,k,N)
	
	if switch_prob == "Fermi":
		return 0.5 + 0.5 * np.tanh(s * (focal_group - role_group))
	
	elif switch_prob == "local":
		if focal_group == 0. and role_group == 0.:
			return 0.5
		else:
			return 0.5 + 0.5 * ((focal_group - role_group) / (np.abs(focal_group) + np.abs(role_group)))


"""
Calculating average group-level victory probability for z-punisher groups over u-punisher
groups for (z,u) \in [i/N,(i+1)/N] \times [j/N,(j+1)/N] using the trapezoidal rule and
our finite volume assumption that the density is a piecewise-constant function taking
constant values on each grid volume.
"""	
def group_switch_terms(i,j,N,s,type,b,c,p,q,k):
	
	ll = group_switch_prob(float(i)/N,float(j)/N,s,type,b,c,p,q,k,N)
	lr = group_switch_prob((i+1.)/N,float(j)/N,s,type,b,c,p,q,k,N)
	ul = group_switch_prob(float(i)/N,(j+1.)/N,s,type,b,c,p,q,k,N)
	ur = group_switch_prob((i+1.)/N,(j+1.)/N,s,type,b,c,p,q,k,N)
	
	return 0.25 * (ll + lr + ul + ur) 



"""
Further characterizing the group-level victory probabilities for each grid volume, and
using these calculations to describe the effect of pairwise group-level competition on
the dynamics of multilevel selection.
"""
	
def group_switch_matrix(N,s,type,b,c,p,q,k):
	
	matrix = np.zeros((N,N))
	for i in range(N):
		for j in range(N):
			matrix[i,j] = group_switch_terms(i,j,N,s,type,b,c,p,q,k)
	return matrix
	

#print(group_function(1.,"coop",b,c,p,q,k,N))
#print(group_switch_prob(0.,1.,s,"coop",b,c,p,q,k,N))
#print(group_switch_terms(N,N,N,s,"coop",b,c,p,q,k))

group_matrix = group_switch_matrix(N,s,type,b,c,p,q,k)


def between_group_term(f,N,s,type,group_matrix,b,c,p,q,k):
	return (1. / N) * f * ( np.dot(group_matrix,f) - np.dot(np.transpose(group_matrix),f))


def group_reproduction_rate(f,N,s,type,group_matrix,b,c,p,q,k):
	return (1. / N) * np.dot(group_matrix,f) 	
	




Z = [[0,0],[0,0]]
levels = np.arange(0.,time_step * time_length+ time_step,time_step)
CS3 = plt.contourf(Z, levels, cmap=cmap.get_cmap('viridis_r'))
plt.clf()


"""
Running the finite volume simulations for our model of multilevel selection with
pairwise group-level competition, and plotting sample snapshots obtained from the numerical
solution obtained at different points of time.  
"""

for time in range(time_length):

	if time % time_snapshot == 0:
		plt.plot(np.arange(0.5/N,1.0+0.5/N,1.0/N),f_j, color = cmap.viridis_r((np.float(time) / time_length)**1), lw = 3.)
		
	between_group_effect = between_group_term(f_j,N,s,type,group_matrix,b,c,p,q,k)
	within_group_effect = within_group(f_j,N,b,c,p,q,k,index_holder)
	righthandside = lamb * between_group_effect + within_group_effect
	print(np.sum(righthandside))
	f_j = f_j + time_step * righthandside
	#f_j = f_j * (N / np.sum(f_j))
	
	print((1.0 / N) * np.sum(f_j))

	


plt.xlabel(r"Fraction of Altruistic Punishers ($z$)", fontsize = 20.)
plt.ylabel(r"Probability Density ($f(t,z)$)", fontsize = 20.)

plt.colorbar(CS3) 

plt.tight_layout()

script_folder = os.getcwd()
altruistic_folder = os.path.dirname(script_folder)
file_path = altruistic_folder

if lamb == 0.1:
	plt.savefig(altruistic_folder + "/Figures/altruistic_nonlinear_trajectory_delta.png")
elif lamb == 2.0:
	plt.savefig(altruistic_folder + "/Figures/altruistic_nonlinear_trajectory_steady_state.png")




plt.show()


