"""
Script used to generate Figure 4.1, which provides illustrations of steady state
densities for our model of multilevel selection with pairwise group-level competition
following the globally normalized group-level update rule (defined in Equation 3.8).
The figure plots steady state densities for various strengths of punishment p for two
example values of the per-interaction cost of punishment k.  
"""


import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as spi
import os

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)


theta = 2.0
s = 0.5

x_step  = 0.00125
x_range = np.arange(x_step, 1.0 +x_step, x_step)


"""Defining the formula for the steady state densities in terms of the parameters
characterizing the game-theoretic interactions, the costs and strength of punishment,
and the dynamics of the two-level selection process.
"""


def G(z,b,c,p,q,k):
	return (b - (c + p + q + k)) * z + (p+k) * (z ** 2.)
	
def Gmin(b ,c ,p, q, k):
	
	if b - c - p - q - k < 0:
		num = -((b - (c + p + q + k)) ** 2.)
		denom = 4. * (p+k)
		return  (num / denom)
	else:
		return G(0.,b,c,p,q,k)
		
Gmin_vec = np.vectorize(Gmin)
	
def steady_state_density(y,lamb,c,q,k,p,theta):

	G_diff_1 = b - c - q
	Gdenom = G(1.,b,c,p,q,k) - Gmin_vec(b,c,p,q,k)
	#print(G(1.,b,c,p,q,k))
	#print(Gmin_vec(b,c,p,q,k))
	#print(G_diff_1/Gdenom)
	#Gdenom = 1.

	y_power = (1. / (c + k + q)) * (lamb * (G_diff_1 / Gdenom) - theta * (c + q - p))
	y_eq_power = -(lamb * (b + k) + theta * (k+p) * Gdenom) / (Gdenom * (c + q + k))
	
	return (y**(y_power - 1.)) * ((1. - y) ** (theta - 1.)) * ((c + k + q - (p+k)*y)**(y_eq_power - 1.))
	
steady_vec = np.vectorize(steady_state_density)

plt.figure(1)

"""
Plotting steady-state densities for various punishment strengths p  for the case of
per-interaction cost k = 0.1 (plotted in Figure 4.1, left).
"""

	
b = 2.
c = 1.
q = 0.1
theta = 2.
lamb = 2.
k = 0.1

x_step  = 0.005

q = 0.0
density_plot_1 = steady_vec(x_range,lamb,c,q,k,0.0,theta)
density_plot_1 = density_plot_1 / spi.simps(density_plot_1,x_range)
plt.plot(x_range,density_plot_1, lw = 4., color = plt.cm.YlOrRd(0.15), label = r"$p = 0$")

density_plot_2 = steady_vec(x_range,lamb,c,q,k,0.15,theta)
density_plot_2 = density_plot_2/ spi.simps(density_plot_2,x_range)
plt.plot(x_range,density_plot_2, lw = 4., color = plt.cm.YlOrRd(0.4),label = r"$p = 0.15$")

density_plot_3 = steady_vec(x_range,lamb,c,q,k,0.3,theta)
density_plot_3 = density_plot_3 / spi.simps(density_plot_3,x_range)
plt.plot(x_range,density_plot_3, lw = 4., color = plt.cm.YlOrRd(0.6),label = r"$p = 0.3$")

density_plot_4 = steady_vec(x_range,lamb,c,q,k,0.45,theta)
density_plot_4 = density_plot_4 / spi.simps(density_plot_4,x_range)
plt.plot(x_range,density_plot_4, lw = 4., color = plt.cm.YlOrRd(0.8), label = r"$p =0.45$")

density_plot_5 =steady_vec(x_range,lamb,c,q,k,0.6,theta)
density_plot_5 = density_plot_5 / spi.simps(density_plot_5,x_range)
plt.plot(x_range,density_plot_5, lw = 4., color = plt.cm.YlOrRd(1.), label = r"$p =0.6$")

plt.axis([0.0,1.0,0.0,10.])
plt.legend(loc = "upper center", fontsize = 14.)

plt.xlabel(r"Fraction of Altruistic Punishers ($z$)", fontsize = 20., labelpad = 20.)
plt.ylabel(r"Steady State Density $f^{\lambda}_{\theta}(z)$", fontsize = 20.)


plt.xticks(fontsize = 14.)
plt.yticks(fontsize = 14.)

plt.tight_layout()



script_folder = os.getcwd()
altruistic_punishment_folder = os.path.dirname(script_folder)
plt.savefig(altruistic_punishment_folder + "/Figures/group_local_density_k_model_k_0p1_pvary.png")



print((1. / len(x_range)) * np.dot(G(x_range,b,c,0.,q,k),density_plot_1))
print((1. / len(x_range)) * np.dot(G(x_range,b,c,0.15,q,k),density_plot_2))
print((1. / len(x_range)) * np.dot(G(x_range,b,c,0.3,q,k),density_plot_3))
print((1. / len(x_range)) * np.dot(G(x_range,b,c,0.45,q,k),density_plot_4))
print((1. / len(x_range)) * np.dot(G(x_range,b,c,0.6,q,k),density_plot_5))


"""
Plotting steady-state densities for various punishment strengths p  for the case of
per-interaction cost k = 0.4 (plotted in Figure 4.1, right).
"""

plt.figure(2)


x_step  = 0.005

q = 0.0
k = 4.
p = 0.5
density_plot_1 = steady_vec(x_range,lamb,c,q,k,0.0,theta)
density_plot_1 = density_plot_1 / spi.simps(density_plot_1,x_range)
plt.plot(x_range,density_plot_1, lw = 4., color = plt.cm.YlOrRd(0.15), label = r"$p = 0$")

density_plot_2 = steady_vec(x_range,lamb,c,q,k,0.15,theta)
density_plot_2 = density_plot_2/ spi.simps(density_plot_2,x_range)
plt.plot(x_range,density_plot_2, lw = 4., color = plt.cm.YlOrRd(0.4),label = r"$p = 0.15$")

density_plot_3 = steady_vec(x_range,lamb,c,q,k,0.3,theta)
density_plot_3 = density_plot_3 / spi.simps(density_plot_3,x_range)
plt.plot(x_range,density_plot_3, lw = 4., color = plt.cm.YlOrRd(0.6),label = r"$p = 0.3$")

density_plot_4 = steady_vec(x_range,lamb,c,q,k,0.45,theta)
density_plot_4 = density_plot_4 / spi.simps(density_plot_4,x_range)
plt.plot(x_range,density_plot_4, lw = 4., color = plt.cm.YlOrRd(0.8), label = r"$p =0.45$")

density_plot_5 =steady_vec(x_range,lamb,c,q,k,0.6,theta)
density_plot_5 = density_plot_5 / spi.simps(density_plot_5,x_range)
plt.plot(x_range,density_plot_5, lw = 4., color = plt.cm.YlOrRd(1.), label = r"$p =0.6$")

plt.axis([0.0,1.0,0.0,15.])
plt.legend(loc = "upper center", fontsize = 14.)

plt.xlabel(r"Fraction of Altruistic Punishers ($z$)", fontsize = 20., labelpad = 20.)
plt.ylabel(r"Steady State Density $f^{\lambda}_{\theta}(z)$", fontsize = 20.)


plt.xticks(fontsize = 14.)
plt.yticks(fontsize = 14.)

plt.tight_layout()



plt.savefig(altruistic_punishment_folder + "/Figures/group_local_density_k_model_k4_pvary.png")



print((1. / len(x_range)) * np.dot(G(x_range,b,c,0.,q,k),density_plot_1))
print((1. / len(x_range)) * np.dot(G(x_range,b,c,0.15,q,k),density_plot_2))
print((1. / len(x_range)) * np.dot(G(x_range,b,c,0.3,q,k),density_plot_3))
print((1. / len(x_range)) * np.dot(G(x_range,b,c,0.45,q,k),density_plot_4))
print((1. / len(x_range)) * np.dot(G(x_range,b,c,0.6,q,k),density_plot_5))



"""
Plotting steady-state densities for various punishment strengths p for the case of
fixed punishment cost q = 0.1  (figure not displayed in paper).
"""

plt.figure(3)


x_step  = 0.005

q = 0.1
k = 0.
p = 0.5
density_plot_1 = steady_vec(x_range,lamb,c,q,k,0.0,theta)
density_plot_1 = density_plot_1 / spi.simps(density_plot_1,x_range)
plt.plot(x_range,density_plot_1, lw = 4., color = plt.cm.YlOrRd(0.15), label = r"$p = 0$")

density_plot_2 = steady_vec(x_range,lamb,c,q,k,0.15,theta)
density_plot_2 = density_plot_2/ spi.simps(density_plot_2,x_range)
plt.plot(x_range,density_plot_2, lw = 4., color = plt.cm.YlOrRd(0.4),label = r"$p = 0.15$")

density_plot_3 = steady_vec(x_range,lamb,c,q,k,0.3,theta)
density_plot_3 = density_plot_3 / spi.simps(density_plot_3,x_range)
plt.plot(x_range,density_plot_3, lw = 4., color = plt.cm.YlOrRd(0.6),label = r"$p = 0.3$")

density_plot_4 = steady_vec(x_range,lamb,c,q,k,0.45,theta)
density_plot_4 = density_plot_4 / spi.simps(density_plot_4,x_range)
plt.plot(x_range,density_plot_4, lw = 4., color = plt.cm.YlOrRd(0.8), label = r"$p =0.45$")

density_plot_5 =steady_vec(x_range,lamb,c,q,k,0.6,theta)
density_plot_5 = density_plot_5 / spi.simps(density_plot_5,x_range)
plt.plot(x_range,density_plot_5, lw = 4., color = plt.cm.YlOrRd(1.), label = r"$p =0.6$")

plt.axis([0.0,1.0,0.0,15.])
plt.legend(loc = "upper center", fontsize = 14.)

plt.xlabel(r"Fraction of Altruistic Punishers ($z$)", fontsize = 20., labelpad = 20.)
plt.ylabel(r"Steady State Density $f^{\lambda}_{\theta}(z)$", fontsize = 20.)


plt.xticks(fontsize = 14.)
plt.yticks(fontsize = 14.)

plt.tight_layout()



plt.savefig(altruistic_punishment_folder + "/Figures/group_local_density_q_model_q0p1_pvary.png")




"""
Plotting steady-state densities for various punishment strengths p for the case of
fixed punishment cost q = 0.9  (figure not displayed in paper).
"""


plt.figure(4)


x_step  = 0.005

q = 0.4
k = 0.
p = 0.5
density_plot_1 = steady_vec(x_range,lamb,c,q,k,0.0,theta)
density_plot_1 = density_plot_1 / spi.simps(density_plot_1,x_range)
plt.plot(x_range,density_plot_1, lw = 4., color = plt.cm.YlOrRd(0.15), label = r"$p = 0$")

density_plot_2 = steady_vec(x_range,lamb,c,q,k,0.15,theta)
density_plot_2 = density_plot_2/ spi.simps(density_plot_2,x_range)
plt.plot(x_range,density_plot_2, lw = 4., color = plt.cm.YlOrRd(0.4),label = r"$p = 0.15$")

density_plot_3 = steady_vec(x_range,lamb,c,q,k,0.3,theta)
density_plot_3 = density_plot_3 / spi.simps(density_plot_3,x_range)
plt.plot(x_range,density_plot_3, lw = 4., color = plt.cm.YlOrRd(0.6),label = r"$p = 0.3$")

density_plot_4 = steady_vec(x_range,lamb,c,q,k,0.45,theta)
density_plot_4 = density_plot_4 / spi.simps(density_plot_4,x_range)
plt.plot(x_range,density_plot_4, lw = 4., color = plt.cm.YlOrRd(0.8), label = r"$p =0.45$")

density_plot_5 =steady_vec(x_range,lamb,c,q,k,0.6,theta)
density_plot_5 = density_plot_5 / spi.simps(density_plot_5,x_range)
plt.plot(x_range,density_plot_5, lw = 4., color = plt.cm.YlOrRd(1.), label = r"$p =0.6$")

plt.axis([0.0,1.0,0.0,15.])
plt.legend(loc = "upper center", fontsize = 14.)

plt.xlabel(r"Fraction of Altruistic Punishers ($z$)", fontsize = 20., labelpad = 20.)
plt.ylabel(r"Steady State Density $f^{\lambda}_{\theta}(z)$", fontsize = 20.)


plt.xticks(fontsize = 14.)
plt.yticks(fontsize = 14.)

plt.tight_layout()



plt.savefig(altruistic_punishment_folder + "/Figures/group_local_density_q_model_q0p9_pvary.png")




plt.show() 