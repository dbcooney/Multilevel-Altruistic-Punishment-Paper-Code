"""
Script used to generate Figure 4.4, which provides comparisons between the average 
payoffs achieved by the dynamics of multilevel selection on the defector-cooperator
and defector-punisher edges of the simplex. This figure compares the average payoff
achieved as a function of the strength \lambda of between-group competition, consider
both the cases of fixed and per-interaction punishment costs. 
"""


import matplotlib.pyplot as plt
import numpy as np
import os


from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)


"""
Characterizing the average payoff function G(z) and the average payoff achieved at steady
state in terms of the payoff parameters and the parameters governing the dynamics
of multilevel selection. 
"""

def G(z,b,c,p,q,k):
	return (b - (c + p + q + k)) * z + (p + k) * (z ** 2.)
	
def Gmin(b,c,p,q,k):
	
	if b - c - p - q - k < 0:
		num = -((b - (c + p + q + k)) ** 2.)
		denom = 4 * (k+p)
		return  (num / denom)
	else:
		return G(0.,b,c,p,q,k)
		
Gmin_vec = np.vectorize(Gmin)
	
def lamb_thresh(b,c,p,q,k,theta):
	Gnum = G(1.,b,c,p,q,k) - Gmin_vec(b,c,p,q,k)
	pi_diff_1 = c + q - p
	G_diff_1 = b - c - q
	
	if pi_diff_1 > 0 and G_diff_1 > 0:
		return (theta * Gnum * pi_diff_1) / (G_diff_1)
	else:
		return 0.
		print("yes")
		
lamb_thresh_vec = np.vectorize(lamb_thresh)

def G_steady(b,c,p,q,k,theta, lamb):
	if lamb < lamb_thresh(b,c,p,q,k,theta):
		return 0.
	elif c + q < p:
		return b - c - q
	else:
		return (b - c - q) * (1.0 - (lamb_thresh(b,c,p,q,k,theta)) / lamb)
	
G_steady_vec = np.vectorize(G_steady)


def G_steady_DC(b,c,theta,lamb):
	if lamb < theta * c:
		return 0.
	else:
		return (b-c) * (1. - (theta / lamb) * c)
		
G_steady_DC_vec = np.vectorize(G_steady_DC)


#def steady_state_payoff()

	
b = 2.
c = 1.
theta = 2.
p = 0.9

lamb_min = 0.
lamb_max = 20.
lamb_step = 0.1
lamb_range = np.arange(lamb_min,lamb_max + lamb_step,lamb_step)


"""
Plotting average payoff achieved at steady state for the two edges of the simplex in the
case of fixed punishment costs with q = 0.25.
"""


plt.figure(1)

q = 0.25
k = 0.0

plt.plot(lamb_range,G_steady_DC_vec(b,c,theta,lamb_range), lw = 5., color = 'r', label = r'$\langle G_{DC}(\cdot)\rangle_{f^{\lambda}_{\theta}}$')
plt.plot(lamb_range,G_steady_vec(b,c,p,q,k,theta,lamb_range), lw = 5., color = 'b', label = r"$\langle G_{DP}(\cdot)\rangle_{f^{\lambda}_{\theta}}$")

plt.axhline(y = b - c, lw = 5., ls = '--', color = 'r', label = r"$G_{DC}(1)$")
plt.axhline(y = b - c - q, lw = 5., ls = '--', color = 'b', label = r"$G_{DP}(1)$")


plt.xlabel(r"Strength of Between-Group Selection $\lambda$", fontsize = 20., labelpad = 10.)
plt.ylabel(r"Average Steady State Payoff", fontsize = 20.)

plt.legend(loc = 'upper center', fontsize =14., ncol = 2, framealpha = 0.)

plt.axis([lamb_min,lamb_max,0.,1.6])

plt.tight_layout()


script_folder = os.getcwd()
altruistic_punishment_folder = os.path.dirname(script_folder)

plt.savefig(altruistic_punishment_folder + "/Figures/DPvsDClowq.png")




"""
Plotting average payoff achieved at steady state for the two edges of the simplex in the
case of per-interaction punishment costs with k = 0.15.
"""


plt.figure(2)

b = 2.
c = 1.

lamb_min = 0.
lamb_max = 25.
lamb_step = 0.1
lamb_range = np.arange(lamb_min,lamb_max + lamb_step,lamb_step)

q = 0.0
k = 0.15

plt.plot(lamb_range,G_steady_DC_vec(b,c,theta,lamb_range), lw = 5., color = 'r', label = r'$\langle G_{DC}(\cdot)\rangle_{f^{\lambda}_{\theta}}$')
plt.plot(lamb_range,G_steady_vec(b,c,p,q,k,theta,lamb_range), lw = 5., color = 'b', label = r"$\langle G_{DP}(\cdot)\rangle_{f^{\lambda}_{\theta}}$")

plt.axhline(y = b - c, lw = 5., ls = '--', color = 'k', label = r"$G_{DC}(1) = G_{DP}(1)$")



plt.xlabel(r"Strength of Between-Group Selection $\lambda$", fontsize = 20., labelpad = 10.)
plt.ylabel(r"Average Steady State Payoff", fontsize = 20.)

plt.legend(loc = 'upper center', fontsize =14., ncol = 2, framealpha = 0.)

plt.axis([lamb_min,lamb_max,0.,1.6])

plt.tight_layout()

plt.savefig(altruistic_punishment_folder + "/Figures/DPvsDClowk.png")





"""
Plotting average payoff achieved at steady state for the two edges of the simplex in the
case of fixed punishment costs with q = 0.75.
"""


b = 2.
c = 1.
theta = 2.
p = 0.9

lamb_min = 0.
lamb_max = 20.
lamb_step = 0.1
lamb_range = np.arange(lamb_min,lamb_max + lamb_step,lamb_step)





plt.figure(3)

q = 0.75
k = 0.0

plt.plot(lamb_range,G_steady_DC_vec(b,c,theta,lamb_range), lw = 5., color = 'r', label = r'$\langle G_{DC}(\cdot)\rangle_{f^{\lambda}_{\theta}}$')
plt.plot(lamb_range,G_steady_vec(b,c,p,q,k,theta,lamb_range), lw = 5., color = 'b', label = r"$\langle G_{DP}(\cdot)\rangle_{f^{\lambda}_{\theta}}$")


plt.axhline(y = b - c, lw = 5., ls = '--', color = 'r', label = r"$G_{DC}(1)$")
plt.axhline(y = b - c - q, lw = 5., ls = '--', color = 'b', label = r"$G_{DP}(1)$")



plt.xlabel(r"Strength of Between-Group Selection $\lambda$", fontsize = 20., labelpad = 10.)
plt.ylabel(r"Average Steady State Payoff", fontsize = 20.)

plt.legend(loc = 'upper center', fontsize =14., ncol = 2, framealpha = 0.)

plt.axis([lamb_min,lamb_max,0.,1.6])

plt.tight_layout()
plt.savefig(altruistic_punishment_folder + "/Figures/DPvsDChigq.png")



"""
Plotting average payoff achieved at steady state for the two edges of the simplex in the
case of per-interaction punishment costs with k = 50.25.
"""


b = 2.
c = 1.
theta = 2.
p = 0.9

lamb_min = 0.
lamb_max = 25.
lamb_step = 0.1
lamb_range = np.arange(lamb_min,lamb_max + lamb_step,lamb_step)


plt.figure(4)

q = 0.0
k = 50.25

plt.plot(lamb_range,G_steady_DC_vec(b,c,theta,lamb_range), lw = 5., color = 'r', label = r'$\langle G_{DC}(\cdot)\rangle_{f^{\lambda}_{\theta}}$')
plt.plot(lamb_range,G_steady_vec(b,c,p,q,k,theta,lamb_range), lw = 5., color = 'b', label = r"$\langle G_{DP}(\cdot)\rangle_{f^{\lambda}_{\theta}}$")

plt.axhline(y = b - c, lw = 5., ls = '--', color = 'k', label = r"$G_{DC}(1) = G_{DP}(1)$")



plt.xlabel(r"Strength of Between-Group Selection $\lambda$", fontsize = 20., labelpad = 10.)
plt.ylabel(r"Average Steady State Payoff", fontsize = 20.)

plt.legend(loc = 'upper center', fontsize =14., ncol = 2, framealpha = 0.)

plt.axis([lamb_min,lamb_max,0.,1.6])

plt.tight_layout()

plt.savefig(altruistic_punishment_folder + "/Figures/DPvsDChighk.png")



plt.show()