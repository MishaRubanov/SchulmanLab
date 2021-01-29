import numpy as np 
import matplotlib.pyplot as plt
import GeneralGeneletModel as GGM

''' 
###############################################################################
This simulates a IFFL12 with 150 nM added dB for IFFL2 (Figure 3f)

###############################################################################
'''

# Defining network
#          G1 G2 G3 G4 G5 G6 G7 G8 G9     
act_vec =  [1, 1, 2, 2, 3, 4, 4, 5, 6]
blk_vec =  [0, 0, 2, 2, 3, 4, 4, 5, 6]
prod_vec = [3, 2,-3, 4, 0, 6, 5,-6, 0]
indc_vec = [0, 0, 0, 0, 0, 0, 0, 0, 0]

''' initializing topology '''
IFFL12 = GGM.GeneletNetwork(act_vec,prod_vec,indc_vec,blk_vec)
IFFL12.plot_topology()
IFFL12.plot_topology(show_rnas=0)

# Define initial conditions
#                  dA1  dA2  dA3  dA4  dA5  dA6
dA_tot = np.array([250, 125, 250, 250, 250, 250]) # total activator added

#                 G1  G2 G3  G4  G5  G6 G7  G8  G9
G_tot = np.array([25, 5, 25, 25, 25, 25, 5, 25, 25]) # total genelet added

# initial genelet states (0 = OFF, 1 = ON, -1 = BLK)
#            G1 G2 G3 G4 G5 G6 G7 G8 G9
G_int_vec = [0, 0,-1,-1,-1,-1,-1,-1,-1]

''' initializing initial conditions '''
IFFL12.initial_conditions(dA_tot,G_tot,G_int_vec,dB_added=[0,0,0,150,150,150]) # default of 0 for all other initial conditions

t_vec1 = np.linspace(0,5,1001)*3600 # seconds

''' simulating the IFFL'''
IFFL12.simulate(t_vec1,1)

# pulling out the desired concentrations for plotting
G1 = IFFL12.output_concentration['GdA1']
G3 = IFFL12.output_concentration['GdA3']
G5 = IFFL12.output_concentration['GdA5']
G9 = IFFL12.output_concentration['GdA9']
sim_t = IFFL12.sol.t

fs = 13
plt.figure()
plt.subplot(3,4,1)
plt.plot(sim_t/3600,G1/G_tot[0],color=[0,0.8,0],linewidth=4)  
plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
plt.xticks(fontsize=fs,weight='bold')
plt.yticks(fontsize=fs,weight='bold')
ax1 = plt.gca()
ax1.set_xlim(0,5)
ax1.set_ylim(-0.1,1.1)
ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')

plt.subplot(3,4,5)
plt.plot(sim_t/3600,G3/G_tot[2],color=[1,0,0],linewidth=4)  
plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
plt.xticks(fontsize=fs,weight='bold')
plt.yticks(fontsize=fs,weight='bold')
ax1 = plt.gca()
ax1.set_xlim(0,5)
ax1.set_ylim(-0.1,1.1)
ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')

plt.subplot(3,4,9)
plt.plot(sim_t/3600,G5/G_tot[3],color=[0,0,1],linewidth=4)  
plt.plot(sim_t/3600,G9/G_tot[8],color=[0.4,0,0.4],linewidth=4)  
plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
plt.xlabel('time (hr)',fontsize=fs,weight='bold')
plt.xticks(fontsize=fs,weight='bold')
plt.yticks(fontsize=fs,weight='bold')
ax1 = plt.gca()
ax1.set_xlim(0,5)
ax1.set_ylim(-0.1,1.1)
ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')

''' Exporting symbolic equations '''
IFFL12.export_sym_eqs(file_name='IFFL12_equations')

''' Compiling genelet sequences to build the network '''
input_file_loc = 'C:\\Users\\sscha\\OneDrive - Johns Hopkins University\\Desktop\\Programming Club\\Python\Python 2020\\General genelet model\\git test\\PythonVersion\\all_genelet_sequences.xlsx'
file_name='IFFL12_sequences'

IFFL12.compile_network_sequences(input_file_loc,desired_nodes=['G2','G3','G1','G5','G6','G4'],save_file_name=file_name)