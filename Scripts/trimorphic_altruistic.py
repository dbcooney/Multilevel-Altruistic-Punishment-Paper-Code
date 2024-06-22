"""
This is our baseline script for running finite volume simulations of the trimorphic
multilevel model for competition in groups featuring defectors, cooperators, and
altruistic punishers. In this script, we consider examples of group-level competition
based on additively separable group-level victory probabilities, exploring group-level
competition based on either the difference in the fraction of non-defecting individuals
or based on normalized differences in the average payoff of group members.

This script is used to produce Figure 6.1 (illustrating the change in the density of
strategic compositions of groups over time) and Figure 6.3 (illustrating differences in
the long-time outcome for different parameters governing costs of punishment and the
relative strength of between-group competition). 


Note: In this script, we use the notation kG to represent the fixed cost of punishment k
and kappa to denote the per-interaction cost of punishment q, and we use j and k to 
describe the indices of the grid volumes used for the finite volume simulation (in
comparison with the indices i and j described in the paper). These discrepancies between
the notation in the code and paper will be fixed in an updated version of the script.  
"""



import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as spi
import os

import matplotlib
matplotlib.use('TkAgg')


from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)


b = 3.
c = 1.
p = 0.5
kappa = 0.1
kG = 0.0

eta = .33



"""
Different options for saving images or text data from the simulations: "trajectory"
saves snapshots of the density after a specified time step, "grouptime" saves a list of 
collective replication rates in the population over the timesteps of the simulation,
and "steady" saves images at the end of long runs (10,000 timesteps) of the simulations.
"""
#quantity = "trajectory"
#quantity = "grouptime"
#quantity = "none"
quantity = "steady"



"""
Choosing the net group-level reproduction function for our numerical simulation, either
consider group-level victory based on normalized differences in average payoff of group
members (corresponding to group_rate_type = "average payoff") or in the difference in
non-defecting individuals (corresponding to group_rate_type = "fraction cooperating").
"""

group_rate_type = "average payoff"
#group_rate_type = "fraction cooperating"


"""
Setting parameters for multilevel competition model and size of numerical grid.
"""

#eta = 0.
N = 100
lamb = 3.
0.


"""
Specifying the time-steppping parameters for the simulations.
"""
timesteps = 40000
time_increments = 0.0015



"""
Specifying possible initial distributions.
"""

uniform_dist = np.zeros((N,N))
partial_uniform_dist = np.zeros((N,N))


index_sum = np.zeros((N,N))

for j in range(N):
	for k in range(N):
		index_sum[j,k] = j + k


for j in range(N):
	for k in range(N-j):
		uniform_dist[j,k] = 2.0
		
		if j <= 10 and k <= 10:
			partial_uniform_dist[j,k] = 10.


"""
Discretized group-level reproduction function (as derived in Section E.2.2).
"""

def Gjk(j,k,b,c,p,kappa,kG,N):

	N = float(N)
	x = float(j) / N
	y = float(k) / N
	
	if j + k < N-1:
		if group_rate_type == "average payoff":
			total = b - c  - (b - c + kappa + p) * (y + 0.5 / N)  \
			+ (kappa + p) * (x * y + (0.5 / N) * (x + y) + 0.25 / (N**2.))  \
			+ (kappa + p) * (y**2. + y / (N) + 1. / (3. * N**2.)) - kG * (1. - x - y - 1./N)
			
		elif group_rate_type == "fraction cooperating":
			total = 1.0 - (y + 0.5 / N)
		
		return total
	elif j + k == N-1:
		"""total = x - eta * x * x + (0.5 / N) - ((7.0 * eta) / (24.0 * N * N)) - ((x * eta) / N )
		"""
		if group_rate_type == "average payoff":
			total = b - c  - (b-c + kappa + p) * (1. - x - 2. / (3. * N))  \
			+ (kappa + p) * (x + 1. / (3. * N) - x**2. - x / N - 1. / (4. * N**2.)) \
			+ (kappa + p) * ((1. - x)**2. - (4. * (1. - x)) / (3. * N) + 1. / (4. * N**2.)) \
			- kG * (1.0/N - 1./(3. * (N**2.)))
			
		
		elif group_rate_type == "fraction cooperating":
			total = x - 1. / (3. * (N**2.)) + 1. / N
		
		return total
	else:
		return None
		
		
Gjk_vec = np.vectorize(Gjk)


def group_birth(j,k,b,c,p,kappa,kG,N):
	if j + k < N - 1:
		return Gjk(j,k,b,c,p,kappa,kG,N)
	elif j + k == N-1:
		return Gjk(j,k,b,c,p,kappa,kG,N)
	else:
		return 0.



group_birth_values = np.zeros((N,N))	
cooperator_values = np.zeros((N,N))
punisher_values = np.zeros((N,N))
	
for j in range(N):
	for k in range(N):
		group_birth_values[j,k] = group_birth(j,k,b,c,p,kappa,kG,N)	
		jx = np.float(j)
		ky = np.float(k)
		if j + k < N - 1:
			cooperator_values[j,k] = jx / N + 1.0 / (2.0 * N)
			punisher_values[j,k] = 1. - (jx / N + 1.0 / (2.0 * N)) - (ky / N + 1.0 / (2.0 * N))
		elif j+k == N - 1:
			cooperator_values[j,k] = jx / N + 1.0 / (3.0 * N)
			punisher_values[j,k] = 1. / N - (2. / (3. * (N ** 2.)))
			


		


spatial_grid = np.zeros((N+1,N+1))
step = 1./float(N)
x_range = np.arange(0.0,1.0 + step, step)
y_range = np.arange(0.0,1.0 + step, step)

def cell_weights(j,k,N):
	if j + k < N - 1:
		return 1.
	elif j + k == N-1:
		return 0.5
	else:
		return 0.
		

		
		
def has_top_right(j,k,N):
	if j+k < N-1:
		return 1.0
	else:
		return 0.0	
		
		
"""
Defining fluxes across the four possible boundaries for a given volume. These fluxes 
correspond to the effects of individual-level birth-death events.
"""		

def left_flux(j,k,N,b,c,p,kappa,kG):

	if j + k >= N:
		return 0.

	N = float(N)
	xj = float(j) / N
	yk1 = float(k+1.0) / N
	yk = float(k)/ N
	
	
	linear_coeff = kG * xj - kG * ((xj) ** 2.)
	quadratic_coeff = 0.5 * (kappa + p - c - kG) * xj - 0.5 * (kappa + p) * (xj ** 2.)
	cubic_coeff = -(1. / 3.) * (kappa + p) * xj
	
	

	return linear_coeff * (yk1 - yk) + quadratic_coeff * (yk1**2. - yk**2.) \
	+ cubic_coeff * (yk1**3. - yk**3.)
	
	
	
def right_flux(j,k,N,b,c,p,kappa,kG):

	if j + k >= N-1:
		return 0.

	N = float(N)
	xj1 = float(j+1.0) / N
	yk1 = float(k+1.0) / N
	yk = float(k)/ N
	
	
	linear_coeff = kG * xj1 - kG * ((xj1) ** 2.)
	quadratic_coeff = 0.5 * (kappa + p - c - kG) * xj1 - 0.5 * (kappa + p) * (xj1 ** 2.)
	cubic_coeff = -(1. / 3.) * (kappa + p) * xj1
	
	
	return -(linear_coeff * (yk1 - yk) + quadratic_coeff * (yk1**2. - yk**2.)) \
	- cubic_coeff * (yk1**3. - yk**3.)
	
	
	
def bottom_flux(j,k,N,b,c,p,kappa,kG):

	if j + k >= N:
		return 0.
	
	N = float(N)
	xj = float(j) / N
	xj1 = float(j + 1.0) / N
	yk = float(k)/ N
	
	linear_coeff = (c + kG -p) * yk + (kappa-c + 2. * p - kG) * (yk ** 2.) - (p+kappa) * (yk ** 3.)
	quadratic_coeff = 0.5 * ((p - kG) * yk - (p + kappa) * (yk **2.))
	
	
	return linear_coeff * (xj1 - xj) + quadratic_coeff * (xj1**2. - xj**2.)	
	
def top_flux(j,k,N,b,c,p,kappa,kG):

	if j + k >= N-1:
		return 0.
	
	N = float(N)
	xj = float(j) / N
	xj1 = float(j + 1.0) / N
	yk1 = float(k+1.0)/ N
	
	linear_coeff = (c + kG -p) * yk1 + (kappa-c + 2. * p - kG) * (yk1 ** 2.) - (p+kappa) * (yk1 ** 3.)
	quadratic_coeff = 0.5 * ((p - kG) * yk1 - (p + kappa) * (yk1 **2.))
	
	return -(linear_coeff * (xj1 - xj) + quadratic_coeff * (xj1**2. - xj**2.))	





cell_weights = np.zeros((N,N))
inv_cell_weights = np.zeros((N,N))

for j in range(N):
	for k in range(N):
		if j + k < N-1:
			cell_weights[j,k] = 1.
			inv_cell_weights[j,k] = 1.
		elif j + k == N-1:
			cell_weights[j,k] = 0.5
			inv_cell_weights[j,k] = 2.
		


def group_righthand(state,group_birth_values,b,c,p,kappa,kG,N,cell_weights):
	weighted = np.multiply(state,group_birth_values)
	increase = weighted
	decrease = (1.0 / (float(N) ** 2.)) * np.sum(np.multiply(weighted,cell_weights)) * state
	return increase-decrease


"""
Calculating fluxes for the boundaries of each volume in our grid.
"""	
left_flux_array = np.zeros((N,N))
right_flux_array = np.zeros((N,N))
bottom_flux_array = np.zeros((N,N))
top_flux_array = np.zeros((N,N))

for j in range(N):
	for k in range(N):
		left_flux_array[j,k] = left_flux(j,k,N,b,c,p,kappa,kG)
		right_flux_array[j,k] = right_flux(j,k,N,b,c,p,kappa,kG)
		bottom_flux_array[j,k] = bottom_flux(j,k,N,b,c,p,kappa,kG)
		top_flux_array[j,k] = top_flux(j,k,N,b,c,p,kappa,kG)
	



"""
Choosing initial condition.
"""		

	
state = uniform_dist
#state = partial_uniform_dist


cell_indices = np.zeros((N,N,2))
for j in range(N):
	for k in range(N):
		cell_indices[j,k,0] = j
		cell_indices[j,k,1] = k



"""
Calculating impact of indvidual-level birth-death competition dynamics on finite volume 
representation of the density f(t,x,y). 
"""	


def within_group(state,left_flux_array,right_flux_array,bottom_flux_array,top_flux_array,cell_weights):
	

	
	up_roll = np.roll(state,-1,axis = 0)
	up_roll[N-1,:] = 0.
	
	down_roll = np.roll(state,1,axis = 0)
	down_roll[0,:] = 0.
	
	left_roll = np.roll(state,-1,axis = 1)
	left_roll[:,N-1] = 0.
	
	right_roll = np.roll(state,1,axis = 1)
	right_roll[:,0] = 0.
	
	
	
	top_flux_from_top = np.where(top_flux_array > 0.,1.0,0.0)
	top_flux_from_bottom = np.where(top_flux_array < 0.,1.0,0.0)
	
	bottom_flux_from_top = np.where(bottom_flux_array < 0.,1.0,0.0)
	bottom_flux_from_bottom = np.where(bottom_flux_array > 0.,1.0,0.0)

	
	right_flux_from_right = np.where(right_flux_array > 0.,1.0,0.0)
	right_flux_from_left = np.where(right_flux_array < 0.,1.0,0.0)
	
	left_flux_from_right = np.where(left_flux_array < 0.,1.0,0.0)
	left_flux_from_left = np.where(left_flux_array > 0.,1.0,0.0)
	
	
	top_contribution = top_flux_array * top_flux_from_top * left_roll + top_flux_array * top_flux_from_bottom * state
	bottom_contribution = bottom_flux_array * bottom_flux_from_top * state + bottom_flux_array * bottom_flux_from_bottom * right_roll

	right_contribution = right_flux_array * right_flux_from_right * up_roll + right_flux_array * right_flux_from_left * state
	left_contribution = left_flux_array * left_flux_from_right * state + left_flux_array * left_flux_from_left * down_roll
	
	sum = bottom_contribution + top_contribution + left_contribution + right_contribution
	
	
	
	weighted_sum = np.multiply(sum,inv_cell_weights)
	#weighted_sum = sum
	
	return (N**2.) * weighted_sum 
	#return weighted_sum
	


"""
Heatmap plotting group-level reproduction function $G(x,y,z)$ on the simplex.
"""

	
G_values = np.zeros((N+1,N+1))

for j in range(N+1):
	for k in range(N+1):
		G_values[j,k] = Gjk(j,k,b,c,p,kappa,kG,N)
		


cmap = plt.get_cmap('YlOrRd')
#cmap = plt.jet()
cmap.set_bad('w',None)

G_values = np.fliplr(G_values)
G_values = np.transpose(G_values)

plt.imshow(G_values)
plt.colorbar(pad = 0.02)

x_ticks = [0.0,0.2,0.4,0.6,0.8,1.0]
y_ticks = [1.0,0.8,0.6,0.4,0.2,0.0]
plt.xticks(range(0,N+1,20),x_ticks)
plt.yticks(range(0,N+1,20),y_ticks)

plt.xticks(fontsize = 14., rotation = 0)
plt.yticks(fontsize = 14.)





plt.tick_params(top = False, right = False)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['top'].set_visible(False)


script_folder = os.getcwd()
protocell_folder = os.path.dirname(script_folder)



def cooperator_count(state,N,cell_weights):
	weighted = np.multiply(state,cooperator_values)
	return (1.0 / (float(N) ** 2.)) * np.sum(np.multiply(weighted,cell_weights))



def punisher_count(state,N,cell_weights):
	weighted = np.multiply(state,punisher_values)
	return (1.0 / (float(N) ** 2.)) * np.sum(np.multiply(weighted,cell_weights))




plt.figure(2)


avgGholder = []

"""
Running the finite volume simulation for the trimorphic multilevel competition.
"""

for time in range(timesteps):
	old_state = state 
	 
	within = within_group(state,left_flux_array,right_flux_array,bottom_flux_array,top_flux_array,inv_cell_weights)
	between = group_righthand(state,group_birth_values,b,c,p,kappa,kG,N,cell_weights) 
	if group_rate_type == "average payoff":
		between = between * (1. / (np.nanmax(G_values) - np.nanmin(G_values)))
	righthand =  within + lamb * between 
	state = state + time_increments * righthand
	
	if time % 40 == 0:
		print(np.sum(np.multiply(state,cell_weights)))
		print(state)
	
	
	print(np.sum(np.multiply(state,cell_weights)))
	weighted = np.multiply(state,group_birth_values)
	

	avgGholder.append((1.0 / (float(N) ** 2.)) * np.sum(np.multiply(weighted,cell_weights)))
	
	print((1.0 / (float(N) ** 2.)) * np.sum(np.multiply(weighted,cell_weights))) 

print("Cooperator Count")	
average_coop = cooperator_count(state,N,cell_weights)
print(average_coop)
print(punisher_count(state,N,cell_weights))
print("Hello")
 
 
weighted = np.multiply(state,group_birth_values)
print((1.0 / (float(N) ** 2.)) * np.sum(np.multiply(weighted,cell_weights)))



"""
Making heat map of the final state achieved as a numerical solution for the multilevel
competition.
"""

state = np.where(index_sum >= N,np.nan,state)
state = np.fliplr(state)
state = np.transpose(state)




#cmap = plt.get_cmap('viridis_r')\
cmap = plt.get_cmap('YlOrRd')
#cmap = plt.jet()
cmap.set_bad('w',None)



plt.imshow(state,cmap =cmap)
plt.colorbar()

x_ticks = [0.0,0.2,0.4,0.6,0.8,1.0]
y_ticks = [1.0,0.8,0.6,0.4,0.2,0.0]
plt.xticks(range(0,N+1,20),x_ticks)
plt.yticks(range(0,N+1,20),y_ticks)
print(np.arange(0.,1.0 + .1,.1))

plt.xticks(fontsize = 14., rotation = 0)
plt.yticks(fontsize = 14.)

#plt.tight_layout()

plt.tick_params(top = False, right = False)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['top'].set_visible(False)





	
print(G_values)



def lamb_thresh(b,c,p,q,theta,Gvalues):
	#Gnum = G(1.,b,c,p,q) - Gmin_vec(b,c,p,q)
	Gnum = (np.nanmax(G_values) - np.nanmin(G_values))
	#Gnum = 1.
	#Gnum = G_range
	pi_diff_1 = c + q - p
	G_diff_1 = b - c - q
	
	if pi_diff_1 > 0 and G_diff_1 > 0:
		return (theta * Gnum * pi_diff_1) / (G_diff_1)
	else:
		return 0.
		print("yes")
		
lamb_thresh_vec = np.vectorize(lamb_thresh)

def G_steady(b,c,p,q,theta, lamb,G_values):
	if lamb < lamb_thresh(b,c,p,q,theta,G_values):
		return 0.
	elif c + q < p:
		return b - c - q
	else:
		return (b - c - q) * (1.0 - (lamb_thresh(b,c,p,q,theta,G_values)) / lamb)
	
G_steady_vec = np.vectorize(G_steady)


print(G_steady(b,c,p,kG,1.,lamb,G_values))
print(G_steady(b,c,0.,0.,1.,lamb,G_values))

def cooperation_CD(b,c,theta,lamb):
	if lamb < theta * c:
		return 0.
	else:
		return 1. - c * (theta / lamb)
		
def cooperation_PD(b,c,q,p,theta,lamb):
	if c + q - p < 0:
		return 1.
	elif lamb < (c + q - p) * theta:
		return 0.
	else:
		return 1. - (c + q  - p) * (theta / lamb)
		
cooperation_CD_vec = np.vectorize(cooperation_CD)
cooperation_PD_vec = np.vectorize(cooperation_PD)
	
print(cooperation_CD(b,c,1.,lamb))
print(cooperation_PD(b,c,kG,p,1.,lamb))


script_folder = os.getcwd()
altruistic_folder = os.path.dirname(script_folder)


plt.xlabel(r"Fraction of Cooperators $x$", fontsize = 24.,labelpad = 10.)
plt.ylabel(r"Fraction of Defectors $y$", fontsize = 24.)



"""
Plotting numerical solution achieved at given time-steps, which will be used to produce 
trajectories for trimorphic multilevel dynamics (as shown in Figure 6.1).
"""

if quantity == "trajectory" and group_rate_type == "fraction cooperation":
	
	if timesteps == 10:
		plt.title(r"$10$ timesteps", fontsize = 24.)
		plt.tight_layout()
		plt.savefig(altruistic_folder + "/Figures/fvtrajectorycooptime10.png",transparent = True,bbox_inches='tight',pad = 0)
	
	elif timesteps == 500:
		plt.title(r"$500$ timesteps", fontsize = 24.)
		plt.tight_layout()
		plt.savefig(altruistic_folder + "/Figures/fvtrajectorycooptime500.png",transparent = True,bbox_inches='tight',pad = 0)


	elif timesteps == 1000:
		plt.title(r"$1{,}000$ timesteps", fontsize = 24.)
		plt.tight_layout()
		plt.savefig(altruistic_folder + "/Figures/fvtrajectorycooptime1000.png",transparent = True,bbox_inches='tight',pad = 0)


	elif timesteps == 2500:
		plt.title(r"$2{,}500$ timesteps", fontsize = 24.)
		plt.tight_layout()
		plt.savefig(altruistic_folder + "/Figures/fvtrajectorycooptime2500.png",transparent = True,bbox_inches='tight',pad = 0)


	elif timesteps == 5000:
		plt.title(r"$5{,}000$ timesteps", fontsize = 24.)
		plt.tight_layout()
		plt.savefig(altruistic_folder + "/Figures/fvtrajectorycooptime5000.png",transparent = True,bbox_inches='tight',pad = 0)


	elif timesteps == 10000:
		plt.title(r"$10{,}000$ timesteps", fontsize = 24.)
		plt.tight_layout()
		plt.savefig(altruistic_folder + "/Figures/fvtrajectorycooptime10000.png",transparent = True,bbox_inches='tight',pad = 0)


	elif timesteps == 20000:
		plt.title(r"$20{,}000$ timesteps", fontsize = 24.)
		plt.tight_layout()
		plt.savefig(altruistic_folder + "/Figures/fvtrajectorycooptime20000.png",transparent = True,bbox_inches='tight',pad = 0)


"""
Plotting the long-time steady state and saving the resulting image (to be displayed in
Figure 6.2).
"""


if quantity == "steady" and group_rate_type == "average payoff" and timesteps == 40000:

	if kG == 0.0 and kappa == 0.1 and lamb == 2.:
		plt.title(r"$q = 0$, $k = 0.1$, $\lambda = 2$", fontsize = 24.)
		plt.tight_layout()
		plt.savefig(altruistic_folder + "/Figures/fvtrimorphicsteadyk0p1lamb2t40000.png",transparent = True,bbox_inches='tight',pad = 0)
	
	if kG == 0.0 and kappa == 0.1 and lamb == 3.:
		plt.title(r"$q = 0$, $k = 0.1$, $\lambda = 3$", fontsize = 24.)
		plt.tight_layout()
		plt.savefig(altruistic_folder + "/Figures/fvtrimorphicsteadyk0p1lamb3t40000.png",transparent = True,bbox_inches='tight',pad = 0)

	
	elif kG == 0.0 and kappa == 0.1 and lamb == 10.:
		plt.title(r"$q = 0, k = 0.1, \lambda = 10$", fontsize = 24.)
		plt.tight_layout()
		plt.savefig(altruistic_folder + "/Figures/fvtrimorphicsteadyk0p1lamb10t40000.png",transparent = True,bbox_inches='tight',pad = 0)
		
	elif kG == 0.1 and kappa == 0.0 and lamb == 2.:
		plt.title(r"$q = 0.1, k = 0, \lambda = 2$", fontsize = 24.)
		plt.tight_layout()
		plt.savefig(altruistic_folder + "/Figures/fvtrimorphicsteadyq0p1lamb2t40000.png",transparent = True,bbox_inches='tight',pad = 0)

	elif kG == 0.1 and kappa == 0.0 and lamb == 3.:
		plt.title(r"$q = 0.1, k = 0, \lambda = 3$", fontsize = 24.)
		plt.tight_layout()
		plt.savefig(altruistic_folder + "/Figures/fvtrimorphicsteadyq0p1lamb3t40000.png",transparent = True,bbox_inches='tight',pad = 0)


	elif kG == 0.1 and kappa == 0.0 and lamb == 10.:
		plt.title(r"$q = 0.1, k = 0, \lambda = 10$", fontsize = 24.)
		plt.tight_layout()
		plt.savefig(altruistic_folder + "/Figures/fvtrimorphicsteadyq0p1lamb10t40000.png",transparent = True,bbox_inches='tight',pad = 0)

plt.show()
		
		