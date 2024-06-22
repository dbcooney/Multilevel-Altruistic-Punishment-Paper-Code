"""
This is the script used to run simulations for trimorphic multilevel competition on the 
defector-cooperator-punisher simplex with replicator group-level dynamics, looping through 
the multilevel dynamics achieved for various strengths \lambda of between-group competition. 
These simulations are used to generate the simulation data required to produce Figure B.1
and the left panels of Figure B.2 (which provides comparison with the trimorphic dynamics
with globally normalized pairwise group-level conflict). 


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
kappa = 2.0
kG = 0.0




"""
Different options for saving images or text data from the simulations: "trajectory"
saves snapshots of the density after a specified time step, "grouptime" saves a list of 
collective replication rates in the population over the timesteps of the simulation,
and "steady" saves images at the end of long runs (10,000 timesteps) of the simulations.
"""
#quantity = "trajectory"
#quantity = "grouptime"
quantity = "none"
#quantity = "steady"


            

group_rate_type = "average payoff"
#group_rate_type = "fraction cooperating"


"""
Setting parameters for multilevel competition model and size of numerical grid.
"""

#eta = 0.
N = 100
lamb = 6.
0.


"""
Specifying the time-steppping parameters for the simulations.
"""
timesteps = 40000
time_increments = 0.015



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
			+ (kappa + p) * (x * y + (0.5 / N) * (x + y) + 0.25 / (N**4.))  \
			+ (kappa + p) * (y**2. + y / (N) + 1. / (3. * N**2.)) - kG * (1. - x - y - 1./N)
			
		elif group_rate_type == "fraction cooperating":
			total = 1.0 - (y - 0.5 / N)
		
		
		return total
	elif j + k == N-1:
		
		if group_rate_type == "average payoff":
			total = b - c  - (b-c + kappa + p) * (1. - x - 2. / (3. * N))  \
			+ (kappa + p) * (x + 1. / (3. * N) - x**2. - x / N - 1. / (4. * N**2.)) \
			+ (kappa + p) * ((1. - x)**2. - (4. * (1. - x)) / (3. * N) + 1. / (4. * N**2.)) \
			- kG * (1.0/N - 1./(3. * (N**2.)))
		
		elif group_rate_type == "fraction cooperating":
			#total = x + 1. / (3. * (N**2.))
			#1. - (y + 1. / (3. * (N**2.)))
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
			#punisher_values[j,k] = (N - jx) / N - 1.0 / (3.0 * N)
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
correspond to indvidual-level birth-death dynamics (as derived in Section E.2.1).
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
	cubic_coeff = -(1. / 3.) * (kappa + p) *xj1
	
	return -(linear_coeff * (yk1 - yk) + quadratic_coeff * (yk1**2. - yk**2.)) \
	- cubic_coeff * (yk1**3. - yk**3.)
	
	
	
def bottom_flux(j,k,N,b,c,p,kappa,kG):

	if j + k >= N:
		return 0.
	
	N = float(N)
	xj = float(j) / N
	xj1 = float(j + 1.0) / N
	yk = float(k)/ N
	


	linear_coeff = (c+kG-p) * yk + (kappa-c + 2. * p - kG) * (yk ** 2.) - (p+kappa) * (yk ** 3.)
	quadratic_coeff = 0.5 * ((p-kG) * yk - (p + kappa) * (yk **2.))
	
	return linear_coeff * (xj1 - xj) + quadratic_coeff * (xj1**2. - xj**2.)	
	
def top_flux(j,k,N,b,c,p,kappa,kG):

	if j + k >= N-1:
		return 0.
	
	N = float(N)
	xj = float(j) / N
	xj1 = float(j + 1.0) / N
	yk1 = float(k+1.0)/ N
	

	
	linear_coeff = (c+kG-p) * yk1 + (kappa-c + 2. * p - kG) * (yk1 ** 2.) - (p+kappa) * (yk1 ** 3.)
	quadratic_coeff = 0.5 * ((p-kG) * yk1 - (p + kappa) * (yk1 **2.))
	
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
	
print(left_flux_array)	
print(bottom_flux_array)	

print(left_flux_array[85,15])
print(right_flux_array[84,15])
"""
Choosing initial condition.
"""		

print(top_flux_array[24,12])
print(bottom_flux_array[24,13])


"""
Choosing initial condition.
"""		

print(left_flux_array[24,14])
print(right_flux_array[23,14])
	
state = uniform_dist
#state = partial_uniform_dist


cell_indices = np.zeros((N,N,2))
for j in range(N):
	for k in range(N):
		cell_indices[j,k,0] = j
		cell_indices[j,k,1] = k



"""
Calculating impact of individual-level birth-death dynamics on finite volume representation
of the densities. 
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

	return (N**2.) * weighted_sum 



"""
Heatmap plotting group-level reproduction function $G(x,y)$ on the simplex.
"""

	
G_values = np.zeros((N+1,N+1))

for j in range(N+1):
	for k in range(N+1):
		G_values[j,k] = Gjk(j,k,b,c,p,kappa,kG,N)
		


cmap = plt.get_cmap('YlOrRd')
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

lamb_min = 0.
lamb_max = 30.
lamb_step = 0.25
lamb_range = np.arange(lamb_min,lamb_max,lamb_step)

lamb_list = []
avg_G_final = []
average_cooperator_holder = []
average_punisher_holder = []

for lamb in lamb_range:
	
	avgGholder = []
	state = uniform_dist
	
	for time in range(timesteps):
		old_state = state 
	 
		within = within_group(state,left_flux_array,right_flux_array,bottom_flux_array,top_flux_array,inv_cell_weights)
		between = group_righthand(state,group_birth_values,b,c,p,kappa,kG,N,cell_weights)
		if group_rate_type == "average payoff":
			between = between 
		righthand =  within + lamb * between 
		state = state + time_increments * righthand
	
		"""
			if time % 40 == 0:
		print(np.sum(np.multiply(state,cell_weights)))
		print(state)
		"""
	
	
		weighted = np.multiply(state,group_birth_values)
	

		avgGholder.append((1.0 / (float(N) ** 2.)) * np.sum(np.multiply(weighted,cell_weights)))
	
		
	lamb_list.append(lamb)
	avg_G_final.append(avgGholder[-1])
	average_cooperator_holder.append(cooperator_count(state,N,cell_weights))
	average_punisher_holder.append(punisher_count(state,N,cell_weights))
	
	print("Iteration")
	print(lamb)
	print(avgGholder[-1])
	
print("Cooperator Count")	
average_coop = cooperator_count(state,N,cell_weights)
print(average_coop)
print(punisher_count(state,N,cell_weights))
print("Hello")
 
 
weighted = np.multiply(state,group_birth_values)
print((1.0 / (float(N) ** 2.)) * np.sum(np.multiply(weighted,cell_weights)))
#print 0.5 * (1.0 - 0.5 * eta)
#print 1 - eta

"""
Making heat map of the final state achieved as a numerical solution for the multilevel
competition.
"""

state = np.where(index_sum >= N,np.nan,state)
state = np.fliplr(state)
state = np.transpose(state)



cmap = plt.get_cmap('YlOrRd')
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

G_range = np.nanmax(G_values) - np.nanmin(G_values)

def lamb_thresh(b,c,p,q,theta,G_range):
	#Gnum = G(1.,b,c,p,q) - Gmin_vec(b,c,p,q)
	Gnum = (np.nanmax(G_values) - np.nanmin(G_values))
	Gnum = 1.
	#Gnum = G_range
	pi_diff_1 = c + q - p
	G_diff_1 = b - c - q
	
	if pi_diff_1 > 0 and G_diff_1 > 0:
		return (theta * Gnum * pi_diff_1) / (G_diff_1)
	else:
		return 0.
		print("yes")
		
lamb_thresh_vec = np.vectorize(lamb_thresh)

def G_steady(b,c,p,q,theta, lamb,G_range):
	if lamb < lamb_thresh(b,c,p,q,theta,G_range):
		return 0.
	elif c + q < p:
		return b - c - q
		#return 1.
	else:
		return (b - c - q) * (1.0 - (lamb_thresh(b,c,p,q,theta,G_range)) / lamb)
		#return 1. 
	
G_steady_vec = np.vectorize(G_steady)

print("Hello")
print(G_steady(b,c,p,kG,1.,lamb,G_range))
print(G_steady(b,c,0.,0.,1.,lamb,G_range))

plt.figure(3)

plt.plot(lamb_range, avg_G_final, lw = 5.)
plt.plot(lamb_range,G_steady_vec(b,c,p,kG,1.,lamb_range,G_range), lw = 5., ls = "--")
plt.plot(lamb_range,G_steady_vec(b,c,0.,0.,1.,lamb_range,G_range), lw = 5., ls = "--", color = 'k')


script_folder = os.getcwd()
altruistic_folder = os.path.dirname(script_folder)
file_path = altruistic_folder

if kappa == 0. and group_rate_type == "average payoff":
	if kG == 0.1 and p == 0.5:
		file1 = file_path + "/Simulation_Outputs/LMGplottrimorphicgroupq0p1p0p5.txt"
		file2 = file_path + "/Simulation_Outputs/LMCoopplottrimorphicgroupq0p1p0p5.txt"
	elif kG == 0.4 and p == 0.5:
		file1 = file_path + "/Simulation_Outputs/LMGplottrimorphicgroupq0p4p0p5.txt"
		file2 = file_path + "/Simulation_Outputs/LMCoopplottrimorphicgroupq0p4p0p5.txt"
	elif kG == 0.1 and p == 0.7:
		file1 = file_path + "/Simulation_Outputs/LMGplottrimorphicgroupq0p1p0p7.txt"
		file2 = file_path + "/Simulation_Outputs/LMCoopplottrimorphicgroupq0p1p0p7.txt"
	elif kG == 0.1 and p == 0.9:
		file1 = file_path + "/Simulation_Outputs/LMGplottrimorphicgroupq0p1p0p9.txt"
		file2 = file_path + "/Simulation_Outputs/LMCoopplottrimorphicgroupq0p1p0p9.txt"
	elif kG == 0.1 and p == 0.3:
		file1 = file_path + "/Simulation_Outputs/LMGplottrimorphicgroupq0p1p0p3.txt"
		file2 = file_path + "/Simulation_Outputs/LMCoopplottrimorphicgroupq0p1p0p3.txt"
	elif kG == 0.1 and p == 0.1:
		file1 = file_path + "/Simulation_Outputs/LMGplottrimorphicgroupq0p1p0p1.txt"
		file2 = file_path + "/Simulation_Outputs/LMCoopplottrimorphicgroupq0p1p0p1.txt"
	elif kG == 0.05 and p == 0.5:
		file1 = file_path + "/Simulation_Outputs/LMGplottrimorphicgroupq0p05p0p5.txt"
		file2 = file_path + "/Simulation_Outputs/LMCoopplottrimorphicgroupq0p05p0p5.txt"
	elif kG == 0.2 and p == 0.5:
		file1 = file_path + "/Simulation_Outputs/LMGplottrimorphicgroupq0p2p0p5.txt"
		file2 = file_path + "/Simulation_Outputs/LMCoopplottrimorphicgroupq0p2p0p5.txt"
	elif kG == 0.3 and p == 0.5:
		file1 = file_path + "/Simulation_Outputs/LMGplottrimorphicgroupq0p3p0p5.txt"
		file2 = file_path + "/Simulation_Outputs/LMCoopplottrimorphicgroupq0p3p0p5.txt"

		
elif kG == 0. and group_rate_type == "average payoff":

	if kappa == 0.1 and p == 0.5:
		file1 = file_path + "/Simulation_Outputs/LMGplottrimorphicgroupk0p1p0p5.txt"
		file2 = file_path + "/Simulation_Outputs/LMCoopplottrimorphicgroupk0p1p0p5.txt"
		
	elif kappa == 2.0 and p == 0.5:
		file1 = file_path + "/Simulation_Outputs/LMGplottrimorphicgroupk2p0p0p5.txt"
		file2 = file_path + "/Simulation_Outputs/LMCoopplottrimorphicgroupk2p0p0p5.txt"




f1 = open(file1,'w+')
f1.write("between-group selection strength lambda")
f1.write('\n')
f1.write(str(lamb_list)[1:-1])
f1.write('\n')
f1.write("average payoff G at steady state")
f1.write('\n')
f1.write(str(avg_G_final)[1:-1])
f1.write('\n')
f1.close()
	
f2 = open(file2,'w+')
f2.write("between-group selection strength lambda")
f2.write('\n')
f2.write(str(lamb_list)[1:-1])
f2.write('\n')
f2.write("average fraction of cooperators x at steady state")
f2.write('\n')
f2.write(str(average_cooperator_holder)[1:-1])
f2.write('\n')
f2.write("average fraction of punishers y at steady state")
f2.write('\n')
f2.write(str(average_punisher_holder)[1:-1])
f2.write('\n')
f2.close()


print(G_values)

plt.show()
		
		