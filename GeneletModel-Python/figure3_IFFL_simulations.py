
import numpy as np 
import matplotlib.pyplot as plt
import GeneralGeneletModel as GGM

''' 
###############################################################################
This simulates a IFFL with a titration of dA3 concentrations (Figure 3b, simulation)
###############################################################################
'''

# Defining network
#          G1 G2 G3 G4      
act_vec =  [1, 1, 2, 3]
blk_vec =  [0, 0, 2, 3]
prod_vec = [3, 2,-3, 0]
indc_vec = [0, 0, 0, 0]

''' initializing topology '''
IFFL = GGM.GeneletNetwork(act_vec,prod_vec,indc_vec,blk_vec)
IFFL.plot_topology()
IFFL.plot_topology(show_rnas=0)

# Define initial conditions
dA_tot =[0]*3 # initializing dA list
#                     dA1  dA2  dA3
dA_tot[0] = np.array([250, 250, 125]) # total activator added
dA_tot[1] = np.array([250, 250, 250]) # total activator added
dA_tot[2] = np.array([250, 250, 500]) # total activator added

#                 G1  G2  G3  G4
G_tot = np.array([25, 5, 25, 25]) # total genelet added

# initial genelet states (0 = OFF, 1 = ON, -1 = BLK)
#            G1 G2 G3 G4
G_int_vec = [0, 0,-1,-1]

G4 = [0]*3 # initializing G4 for storage
cv = [[0,0,0.5],[0,0,1],[0,0.5,1]] # colors for plotting
fs = 13 # font size for plotting
plt.figure()

for i in range(len(dA_tot)):

    ''' initializing initial conditions '''
    IFFL.initial_conditions(dA_tot[i],G_tot,G_int_vec) # default of 0 for all other initial conditions
    
    t_vec1 = np.linspace(0,3,1001)*3600 # seconds
    
    ''' simulating the IFFL'''
    IFFL.simulate(t_vec1,1)
    
    # pulling out the desired concentrations for plotting
    G4[i] = IFFL.output_concentration['GdA4']
    sim_t = IFFL.sol.t
    
    plt.subplot(3,4,1)
    plt.plot(sim_t/3600,G4[i]/G_tot[3],color=cv[i],linewidth=4)
    
plt.suptitle('IFFL with dA3 titration',fontsize=fs+2,weight='bold')    
plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
plt.xlabel('time (hr)',fontsize=fs,weight='bold')
plt.xticks(fontsize=fs,weight='bold')
plt.yticks(fontsize=fs,weight='bold')
ax1 = plt.gca()
ax1.set_xlim(0,3)
ax1.set_ylim(-0.1,1.1)
ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')

''' Exporting symbolic equations '''
IFFL.export_sym_eqs(file_name='IFFL1_equations')

''' Compiling genelet sequences to build the network '''
input_file_loc = 'C:\\Users\\sscha\\OneDrive - Johns Hopkins University\\Desktop\\Programming Club\\Python\Python 2020\\General genelet model\\git test\\PythonVersion\\all_genelet_sequences.xlsx'
file_name='IFFL1_sequences'

IFFL.compile_network_sequences(input_file_loc,desired_nodes=['G2','G3','G1'],save_file_name=file_name)
