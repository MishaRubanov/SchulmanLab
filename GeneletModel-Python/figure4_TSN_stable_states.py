
import numpy as np 
import matplotlib.pyplot as plt
import GeneralGeneletModel as GGM

'''
###############################################################################
This simulates the TSN initialized in State 2
-------------------------------------------------------------------------------
   -This example sets up the initial condtions like I did experimentally where
    initially only one activator is present and the others are added after 30 m
   -This example compiles sequences without BTH domains
###############################################################################
'''

#          G1 G2 G3 G4 G5 G6      
act_vec =  [1, 1, 2, 2, 3, 3]
blk_vec =  [0, 0, 0, 0, 0, 0]
prod_vec = [-3,-2,-1,-3,-1,-2]
indc_vec = [0, 0, 0, 0, 0, 0]

''' initializing topology '''
TSN = GGM.GeneletNetwork(act_vec,prod_vec,indc_vec,blk_vec)
TSN.plot_topology()

# Define initial conditions
#                  dA1  dA2  dA3
dA_tot = np.array([0, 250, 0]) # total activator added (only dA2 is present to set that state)
#                 G1  G2  G3  G4  G5  G6
G_tot = np.array([50, 50, 50, 50, 50, 50]) # total genelet added

# initial genelet states (0 = OFF, 1 = ON, -1 = BLK)
#            G1 G2 G3 G4 G5 G6
G_int_vec = [0, 0, 1, 1, 0, 0]

''' initializing initial conditions '''
TSN.initial_conditions(dA_tot,G_tot,G_int_vec) # default of 0 for all other initial conditions
TSN.plot_topology()
TSN.plot_topology(show_rnas=0)

t_vec1 = np.linspace(0,0.5,101)*3600 # seconds

''' simulating the TSN model for initial RNA production to set State 2'''
TSN.simulate(t_vec1,1)

t_vec2 = np.linspace(t_vec1[-1]/3600,6.5,1001)*3600 # seconds

''' simulating the TSN model after added the other activators to start experiment '''
TSN.simulate(t_vec2,2,dA=[250,'NA',250])

''' Plotting solution '''

# pulling out solutions for plotting from output_concentration attribute
G1 = TSN.output_concentration['GdA1']
G3 = TSN.output_concentration['GdA3']
G5 = TSN.output_concentration['GdA5']

sim_t = TSN.sol.t

fs = 13
plt.figure()
plt.suptitle('TSN: State 2',fontsize=fs+2,weight='bold')
plt.subplot(3,3,1)
# 0.5 hr shift as that is when the experiment actually starts here
plt.plot(sim_t/3600-0.5,G1/G_tot[0],color='b',linewidth=4)
plt.plot(sim_t/3600-0.5,G3/G_tot[2],color='r',linewidth=4)
plt.plot(sim_t/3600-0.5,G5/G_tot[4],color=[0.4,0,0.4],linewidth=2)
plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
plt.xlabel('time (hr)',fontsize=fs,weight='bold')
plt.xticks(fontsize=fs,weight='bold')
plt.yticks(fontsize=fs,weight='bold')
ax1 = plt.gca()
ax1.set_xlim(0,6)
ax1.set_ylim(-0.1,1.1)
ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')

''' Exporting symbolic equations '''
TSN.export_sym_eqs(file_name='TSN_equations')

''' Compiling genelet sequences to build the network '''
input_file_loc = 'C:\\Users\\sscha\\OneDrive - Johns Hopkins University\\Desktop\\Programming Club\\Python\Python 2020\\General genelet model\\git test\\PythonVersion\\all_genelet_sequences.xlsx'
file_name='TSN_sequences'

# compile sequences without the BTH domains
TSN.compile_network_sequences(input_file_loc,desired_nodes=['G1','G3','G4'],bth=0,save_file_name=file_name)