"""
Script used to generate Figure 2.1, illustrating the phase portrait and sample
trajectories for the the dynamics of individual-level selection described by a replicator
equation for the fractions of cooperators, defectors, and altruistic punishers.
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

"""
s = 1.
bS = 1.
bF = 1. + s
bD = (1. + s) / (2. + s)
"""

p = 0.5
q = 0.0
k = 0.2
b = 2. 
c = 1.

timesteps = 8000
time_diff = 0.01

x_step = 0.001

N = 20

"""
Righthand sides of ODEs for within-group competition.
"""

def piC(x,y,b,c,q,k,p):
	return b * (1. - y) - c
	
def piD(x,y,b,c,q,k,p):
	return b * (1. - y) - p * (1. - x - y)

def piP(x,y,b,c,q,k,p):
	return b * (1. - y) - c - q - k * y
	

def xrighthand(x,y,b,c,q,k,p):
	PiC = piC(x,y,b,c,q,k,p)
	PiD = piD(x,y,b,c,q,k,p)
	PiP = piP(x,y,b,c,q,k,p)
	return x * ((1. - x) * (PiC - PiP) - y * (PiD - PiP))
	
def yrighthand(x,y,b,c,q,k,p):
	PiC = piC(x,y,b,c,q,k,p)
	PiD = piD(x,y,b,c,q,k,p)
	PiP = piP(x,y,b,c,q,k,p)
	return y * ((1.-y) * (PiD - PiP) - x * (PiC - PiP))
	
def simplex_edge(x):
	return 1.0 - x
	
state1_init = [0.5,0.1]
state2_init = [0.31,0.5]
state3_init = [0.51,0.04]

if p == 2.5 and q == 0.2:	
	state1_init = [0.1,0.4]
	state2_init = [0.1,0.35]
	state3_init = [0.1,0.2]
	
elif p == 0.5 and q == 0.2:
	state1_init = [0.1,0.1]
	state2_init = [0.31,0.01]
	state3_init = [0.51,0.004]
	
elif p == 2.5 and k == 0.2:
	state1_init = [0.075,0.45]
	state2_init = [0.15,0.4]
	state3_init = [0.51,0.1]
	
elif p == 0.5 and k == 0.2:
	state1_init = [0.1,0.1]
	state2_init = [0.31,0.01]
	state3_init = [0.51,0.004]


state1 = state1_init 



statex = [state1_init[0]]
statey = [state1_init[1]]

statex2 = [state2_init[0]]
statey2 = [state2_init[1]]

statex3 = [state3_init[0]]
statey3 = [state3_init[1]]




"""
Integrating sample within-group trajectories in time.
"""

for time in range(timesteps):
	
	
	xold = statex[-1]
	yold = statey[-1]
	
	xnew = xold + time_diff * xrighthand(xold,yold,b,c,q,k,p)
	ynew = yold + time_diff * yrighthand(xold,yold,b,c,q,k,p)
	
	statex.append(xnew)
	statey.append(ynew)
	
	
	x2old = statex2[-1]
	y2old = statey2[-1]
	
	x2new = x2old + time_diff * xrighthand(x2old,y2old,b,c,q,k,p)
	y2new = y2old + time_diff * yrighthand(x2old,y2old,b,c,q,k,p)
	
	statex2.append(x2new)
	statey2.append(y2new)
	
	x3old = statex3[-1]
	y3old = statey3[-1]
	
	x3new = x3old + time_diff * xrighthand(x3old,y3old,b,c,q,k,p)
	y3new = y3old + time_diff * yrighthand(x3old,y3old,b,c,q,k,p)
	
	statex3.append(x3new)
	statey3.append(y3new)
	


x_vec = np.arange(0.,1. + x_step,x_step)



"""
Calculating and plotting arrows for within-group phase portrait.
"""
for i in range(N):
	for j in range(N):
		x = np.float(i)/N
		y = np.float(j)/N
		
		dx = xrighthand(x,y,b,c,q,k,p)
		dy = yrighthand(x,y,b,c,q,k,p)
		
		if i + j <= N:
			plt.quiver(x,y,dx,dy, color = 'b', width = 0.005, alpha = 0.8)
	
"""
Plotting sample within-group trajectories.
"""		
plt.plot(statex,statey, lw = 6., color = 'r', alpha = 0.7)	
plt.plot(statex2,statey2, lw = 6., color = 'r', alpha = 0.7, ls = '--')	
plt.plot(statex3,statey3, lw = 6., color = 'r', alpha = 0.7, ls = '-.')	
plt.plot(x_vec,simplex_edge(x_vec), color = 'k', lw = 4., alpha = 0.7)		

plt.tick_params(top = False, right = False)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['top'].set_visible(False)

plt.axis([-0.005,1.,-0.005,1.])

plt.xticks(fontsize = 14.)
plt.yticks(fontsize = 14.)

plt.xlabel(r"Fraction of Cooperators $x$", fontsize = 20.,labelpad = 10.)
plt.ylabel(r"Fraction of Defectors $y$", fontsize = 20.)

plt.tight_layout()

script_folder = os.getcwd()
protocell_folder = os.path.dirname(script_folder)

if q == 0.2 and p == 2.5:
	plt.savefig(protocell_folder + "/Figures/withintrimorphicq0p2p2p5.png",transparent = True,bbox_inches='tight',pad = 0)
elif q == 0.2 and p == 0.5:
	plt.savefig(protocell_folder + "/Figures/withintrimorphicq0p2p0p5.png",transparent = True,bbox_inches='tight',pad = 0)
elif k == 0.2 and p == 0.5:
	plt.savefig(protocell_folder + "/Figures/withintrimorphick0p2p0p5.png",transparent = True,bbox_inches='tight',pad = 0)
elif k == 0.2 and p == 2.5:
	plt.savefig(protocell_folder + "/Figures/withintrimorphick0p2p2p5.png",transparent = True,bbox_inches='tight',pad = 0)


plt.show()
	
	
