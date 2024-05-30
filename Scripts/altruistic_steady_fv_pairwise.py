"""
Script used to generate Figure 5.2, which provides comparisons of the long-time outcome
achieved under numerical simulations of multilevel dynamics with pairwise group-level
competition. The figure provides a comparison of long-time outcomes for different cases
of the strength p of altruistic punishment for multilevel dynamics under the Fermi
group-level victory probability. 
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
time_step = 0.006
time_length = 9600


#type = "coop"
type = "payoff"

switch_prob = "Fermi"
#switch_prob = "local"



#type = "coop"
type = "payoff"

p = 0.6
k = 0.1
c = 1.
b = 2.
q = 0.

gamma = b - c - p - k - q
alpha = p + k
beta = -(c + k)


lamb = 2.
s = 2.
lamb = 1.55


"""
Setting the strength of altruistic punishment p for each of the simulations we will
run. 
"""

p1 = 0.
p2 = 0.2
p3 = 0.4
p4 = 0.6
p5 = 0.8

alpha1 = p1 - k
alpha2 = p2 - k
alpha3 = p3 - k
alpha4 = p4 - k
alpha5 = p5 - k




def theta_init(j,N,theta):
	return N ** (1.0 - theta) * (((N - j) ** theta) - ((N - j - 1.0) ** theta) )
	
theta_vec = np.vectorize(theta_init)

G_vec = np.vectorize(G)
Gj_vec = np.vectorize(G_j)

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

def flux_right(j,N,beta,alpha):
	return ((j+1.0) / N) * (1.0 - (j+1.0) / N) * (beta + alpha * ((j+1.0)/N))
	
def flux_left(j,N,beta,alpha):
	return ((np.float(j)) / N) * (1.0 - (np.float(j)) / N) * (beta + alpha * ((np.float(j))/N))
	



flux_right_vec = np.vectorize(flux_right)
flux_left_vec = np.vectorize(flux_left)

def within_group(f,N,alpha,beta,index_holder):
	left_roll = np.roll(f,-1)
	left_roll[-1] = 0.
	right_roll = np.roll(f,1)
	right_roll[0] = 0.
	
	upper_flux = flux_right_vec(index_holder,N,beta,alpha)
	lower_flux = flux_left_vec(index_holder,N,beta,alpha)
	
	upper_flux_up = np.where(upper_flux < 0.0,1.0,0.0)
	upper_flux_down = np.where(upper_flux > 0.0,1.0,0.0)
	
	lower_flux_up = np.where(lower_flux < 0.0,1.0,0.0)
	lower_flux_down = np.where(lower_flux > 0.0,1.0,0.0)
	
	
	upper_half = upper_flux_up * upper_flux * left_roll + upper_flux_down * lower_flux * f 
	lower_half = lower_flux_up * lower_flux * f + lower_flux_down * lower_flux * right_roll
	return N*(-upper_half + lower_half)
	
	
	
"""
Defining terms used to describe between-group competition. 
"""	
	
def group_function(x,type,b,c,p,q,k,N):
	
	if type == "coop":
		return x
	elif type == "payoff":
		return G(x,b,c,p,q,k) 
		
		
"""
Defining the victory probability of a focal group featuring a fraction y of altruistic 
punishers when paired against a u-punisher group in a group-level conflict.
""" 		
		
def group_switch_prob(x,u,s,type,alpha,gamma):
	
	
	focal_group = group_function(x,type,b,c,p,q,k,N)
	role_group = group_function(u,type,b,c,p,q,k,N)
	
	if switch_prob == "Fermi":
		return 0.5 + 0.5 * np.tanh(s * (focal_group - role_group))
	#return 0.5 + 0.5 * (focal_group - role_group)
	elif switch_prob == "local":
		#print(0.5 + 0.5 * ((focal_group - role_group) / (np.abs(focal_group) + np.abs(role_group))))
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

	
def group_switch_terms(j,k,N,s,type,alpha,gamma):
	
	ll = group_switch_prob(float(j)/N,float(k)/N,s,type,alpha,gamma)
	lr = group_switch_prob((j+1.)/N,float(k)/N,s,type,alpha,gamma)
	ul = group_switch_prob(float(j)/N,(k+1.)/N,s,type,alpha,gamma)
	ur = group_switch_prob((j+1.)/N,(k+1.)/N,s,type,alpha,gamma)
	
	return 0.25 * (ll + lr + ul + ur) 
	
	
"""
Further characterizing the group-level victory probabilities for each grid volume, and
using these calculations to describe the effect of pairwise group-level competition on
the dynamics of multilevel selection.
"""
		
	
def group_switch_matrix(N,s,type,alpha,gamma):
	
	matrix = np.zeros((N,N))
	for j in range(N):
		for k in range(N):
			matrix[j,k] = group_switch_terms(j,k,N,s,type,alpha,gamma)
	return matrix


group_matrix = group_switch_matrix(N,s,type,alpha,gamma)


def between_group_term(f,N,s,type,group_matrix,alpha,gamma):
	return (1. / N) * f * ( np.dot(group_matrix,f) - np.dot(np.transpose(group_matrix),f))
	
	


peak_holder = [float(np.argmax(f_j))/N]


Z = [[0,0],[0,0]]
levels = np.arange(0.,time_step * time_length+ time_step,time_step)
CS3 = plt.contourf(Z, levels, cmap=cmap.get_cmap('viridis'))
plt.clf()

f_j1 = theta_vec(index_holder,N,1.0)
f_j2 = theta_vec(index_holder,N,1.0)
f_j3 = theta_vec(index_holder,N,1.0)
f_j4 = theta_vec(index_holder,N,1.0)
f_j5 = theta_vec(index_holder,N,1.0)



"""
Running the finite volume simulations for our model of multilevel selection with
pairwise group-level competition for each strength p of altruistic punishment considered.
We use these simulations to generate plots of the densities achieved for each strength p
after 9,600 time-steps.
"""


for time in range(time_length):
	
	between_group_effect1 = between_group_term(f_j1,N,s,type,group_matrix,alpha,gamma)
	within_group_effect1 = within_group(f_j1,N,alpha1,beta,index_holder)
	righthandside1 = lamb * between_group_effect1 + within_group_effect1
	f_j1 = f_j1 + time_step * righthandside1
	
	between_group_effect2 = between_group_term(f_j2,N,s,type,group_matrix,alpha,gamma)
	within_group_effect2 = within_group(f_j2,N,alpha2,beta,index_holder)
	righthandside2 = lamb * between_group_effect2 + within_group_effect2
	f_j2 = f_j2 + time_step * righthandside2
	
	between_group_effect3 = between_group_term(f_j3,N,s,type,group_matrix,alpha,gamma) 
	within_group_effect3 = within_group(f_j3,N,alpha3,beta,index_holder)
	righthandside3 = lamb * between_group_effect3 + within_group_effect3
	f_j3 = f_j3 + time_step * righthandside3
	
	between_group_effect4 = between_group_term(f_j4,N,s,type,group_matrix,alpha,gamma)
	within_group_effect4 = within_group(f_j4,N,alpha4,beta,index_holder)
	righthandside4 = lamb * between_group_effect4 + within_group_effect4
	f_j4 = f_j4 + time_step * righthandside4
	
	between_group_effect5 = between_group_term(f_j5,N,s,type,group_matrix,alpha,gamma)
	within_group_effect5 = within_group(f_j5,N,alpha5,beta,index_holder)
	righthandside5 = lamb * between_group_effect5 + within_group_effect5
	f_j5 = f_j5 + time_step * righthandside5
	
	#print (1.0 / N) * np.sum(f_j)
	
plt.plot(np.arange(0.5/N,1.0+0.5/N,1.0/N),f_j1, color = plt.cm.YlOrRd(0.2), lw = 6., label = r"$p = 0$")
plt.plot(np.arange(0.5/N,1.0+0.5/N,1.0/N),f_j2, color = plt.cm.YlOrRd(0.4), lw = 6., label = r"$p = 0.2$")
plt.plot(np.arange(0.5/N,1.0+0.5/N,1.0/N),f_j3, color = plt.cm.YlOrRd(0.6), lw = 6., label = r"$p = 0.4$")
plt.plot(np.arange(0.5/N,1.0+0.5/N,1.0/N),f_j4, color = plt.cm.YlOrRd(0.8), lw = 6., label = r"$p = 0.6$")
plt.plot(np.arange(0.5/N,1.0+0.5/N,1.0/N),f_j5, color = plt.cm.YlOrRd(1.0), lw = 6., label = r"$p = 0.8$")






plt.xlabel(r"Fraction of Altruistic Punishers ($z$)", fontsize = 20., labelpad = 10.)
plt.ylabel(r"Probability Density ($f(z)$)", fontsize = 20.)
plt.legend(loc = "upper center", fontsize = 16.)
plt.tight_layout()

script_folder = os.getcwd()
altruistic_folder = os.path.dirname(script_folder)

plt.savefig(altruistic_folder + "/Figures/altruistic_punishment_steady_compare.png")

plt.show()




