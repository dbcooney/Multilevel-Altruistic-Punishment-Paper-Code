"""
Script used to generate Figures 4.2 and 4.3, providing an illustration of the 
threshold selection strength \lambda^* and the average payoff <G()>_f achieved at 
steady state for the model of multilevel selection with the globally normalized
local pairwise group-level victory probability. These figures illustrate how these two
quantities vary with the strength of punishment p for different possible values of the
per-interaction punishment cost k (Figure 4.2) or the fixed punishment cost q 
(Figure 4.3).
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
Characterizing the average payoff at steady state and the threshold selection strength
required to sustain altruistic punishment at steady state as a function of the 
payoff parameters and the parameters describing multilevel selection. 
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

#def steady_state_payoff()

	
b = 2.
c = 1.
q = 0.85
theta = 2.
k = 0.

lamb = 10.

p_range = np.arange(0.1,6,0.01)






script_folder = os.getcwd()
altruistic_punishment_folder = os.path.dirname(script_folder)



"""
Plotting the threshold selection strength for the survival of altruistic punishment as
a function of the strength of punishment p for the case of fixed punishment 
costs. 
"""

plt.figure(1)



cmap_range = 0.86 - 0.5
cmap_min = 0.5
#plt.plot(p_range, lamb_thresh_vec(b,c,p_range, 0.3,theta), lw = 5., color = plt.cm.viridis((0.3 - cmap_min )/cmap_range ), label = r'$q = 0.2$')
plt.plot(p_range, lamb_thresh_vec(b,c,p_range, 0.5,k,theta), lw = 6., color = plt.cm.viridis((0.5 - cmap_min )/cmap_range ), label = r'$q = 0.5$')
plt.plot(p_range, lamb_thresh_vec(b,c,p_range, 0.7,k,theta), lw = 6., color = plt.cm.viridis((0.65 - cmap_min )/cmap_range ), label = r'$q = 0.7$')
plt.plot(p_range, lamb_thresh_vec(b,c,p_range, 0.8,k,theta), lw = 6., color = plt.cm.viridis((0.77 - cmap_min )/cmap_range ), label = r'$q = 0.8$')
plt.plot(p_range, lamb_thresh_vec(b,c,p_range, 0.86,k,theta), lw = 6., color = plt.cm.viridis((0.86 - cmap_min )/cmap_range ), label = r'$q = 0.86$')


print(lamb_thresh_vec(b,c,p_range,0.95,k,theta))

plt.legend(loc = "upper right",fontsize = 14.)

plt.xlabel(r"Strength of Punishment $p$", fontsize = 20., labelpad = 10.)
plt.ylabel(r"Threshold Selection Strength $\lambda^*$", fontsize = 16., labelpad = 10.)



plt.tight_layout()



plt.savefig(altruistic_punishment_folder + "/Figures/lambda_thresh_q_model.png")


"""
Plotting the average payoff at steady state as a function of the strength of punishment
p for the case of fixed punishment costs. 
"""

plt.figure(2)

cmap_range = 0.925 - 0.5
cmap_min = 0.5

plt.plot(p_range, G_steady_vec(b,c,p_range, 0.5,k,theta,lamb), lw = 6., color = plt.cm.viridis((0.5 - cmap_min )/cmap_range ), label = r'$q = 0.5$')
plt.plot(p_range, G_steady_vec(b,c,p_range, 0.65,k,theta,lamb), lw = 6., color = plt.cm.viridis((0.65 - cmap_min )/cmap_range ), label = r'$q = 0.65$')
plt.plot(p_range, G_steady_vec(b,c,p_range, 0.85,k,theta,lamb), lw = 6., color = plt.cm.viridis((0.85 - cmap_min )/cmap_range ), label = r'$q = 0.85$')
plt.plot(p_range, G_steady_vec(b,c,p_range, 0.925,k,theta,lamb), lw = 6., color = plt.cm.viridis((0.925 - cmap_min )/cmap_range ), label = r'$q = 0.925$')

plt.axis([0.,6.,0.,0.7])

plt.legend(loc = "upper left", fontsize = 14., ncol = 2)

plt.xlabel(r"Strength of Punishment $p$", fontsize = 20., labelpad = 10.)
plt.ylabel(r"Steady State Payoff $\langle G(\cdot) \rangle_{f^{\lambda}_{\theta}}$", fontsize = 16., labelpad = 10.)

plt.tight_layout()

plt.savefig(altruistic_punishment_folder + "/Figures/average_G_q_model.png")




q = 0.



"""
Plotting the threshold selection strength for the survival of altruistic punishment as
a function of the strength of punishment p for the case of per-interaction punishment 
costs. 
"""


plt.figure(3)



cmap_range = 40. - 10.
cmap_min = 10.
#plt.plot(p_range, lamb_thresh_vec(b,c,p_range, 0.3,theta), lw = 5., color = plt.cm.viridis((0.3 - cmap_min )/cmap_range ), label = r'$q = 0.2$')
plt.plot(p_range, lamb_thresh_vec(b,c,p_range, q,10.,theta), lw = 6., color = plt.cm.viridis((10. - cmap_min )/cmap_range ), label = r'$k = 10$')
plt.plot(p_range, lamb_thresh_vec(b,c,p_range, q,20.,theta), lw = 6., color = plt.cm.viridis((20. - cmap_min )/cmap_range ), label = r'$k = 20$')
plt.plot(p_range, lamb_thresh_vec(b,c,p_range, q,30.,theta), lw = 6., color = plt.cm.viridis((30. - cmap_min )/cmap_range ), label = r'$k = 30$')
plt.plot(p_range, lamb_thresh_vec(b,c,p_range, q,40.,theta), lw = 6., color = plt.cm.viridis((40. - cmap_min )/cmap_range ), label = r'$k = 40$')





print(lamb_thresh_vec(b,c,p_range,0.95,k,theta))

plt.legend(loc = "upper right",fontsize = 14.)

plt.xlabel(r"Strength of Punishment $p$", fontsize = 20., labelpad = 10.)
plt.ylabel(r"Threshold Selection Strength $\lambda^*$", fontsize = 16., labelpad = 10.)

plt.axis([0.,2.5,-0.5,20.])

plt.tight_layout()


plt.savefig(altruistic_punishment_folder + "/Figures/lambda_thresh_k_model.png")





"""
Plotting the average payoff at steady state as a function of the strength of punishment
p for the case of per-interaction punishment costs. 
"""


plt.figure(4)

cmap_range = 40. - 10.
cmap_min = 10.

plt.plot(p_range, G_steady_vec(b,c,p_range, q,10.,theta,lamb), lw = 6., color = plt.cm.viridis((10. - cmap_min )/cmap_range ), label = r'$k = 10$')
plt.plot(p_range, G_steady_vec(b,c,p_range, q,20.,theta,lamb), lw = 6., color = plt.cm.viridis((20. - cmap_min )/cmap_range ), label = r'$k = 20$')
plt.plot(p_range, G_steady_vec(b,c,p_range,q, 30.,theta,lamb), lw = 6., color = plt.cm.viridis((30. - cmap_min )/cmap_range ), label = r'$k = 30$')
plt.plot(p_range, G_steady_vec(b,c,p_range,q, 40.,theta,lamb), lw = 6., color = plt.cm.viridis((40. - cmap_min )/cmap_range ), label = r'$k = 40$')

plt.axis([0.,2.5,-0.025,1.05])

plt.legend(loc = "center right", fontsize = 14., ncol = 1)

plt.xlabel(r"Strength of Punishment $p$", fontsize = 20., labelpad = 10.)
plt.ylabel(r"Steady State Payoff $\langle G(\cdot) \rangle_{f^{\lambda}_{\theta}}$", fontsize = 16., labelpad = 10.)

plt.tight_layout()

plt.savefig(altruistic_punishment_folder + "/Figures/average_G_k_model.png")




plt.show()