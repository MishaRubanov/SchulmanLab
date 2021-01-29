
import numpy as np 
import matplotlib.pyplot as plt
import GeneralGeneletModel as GGM

'''
###############################################################################
Simulates the I_BS_IFFL1|2_FB1 network (Figure 5g)
-------------------------------------------------------------------------------

###############################################################################
'''

#          G1 G2 G3 G4 G5 G6 G7 G8 G9 G10 G11 G12 G13 G14   
act_vec =  [1, 1, 2, 2, 3, 3, 4, 5, 6, 6,  7,  8,  9,  10]
blk_vec =  [0, 0, 0, 0, 3, 3, 4, 5, 6, 6,  7,  8,  9,  0]
prod_vec = [-2,3,-1, 6, 4, 5,-5, 9, 7, 8, -8,  0,  0,  0]
indc_vec = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0,  0, -2, -1]

''' initializing topology '''
I_BS_IFFL_FB = GGM.GeneletNetwork(act_vec,prod_vec,indc_vec,blk_vec)

# Define initial conditions
#                  dA1  dA2  dA3  dA4  dA5  dA6  dA7  dA8  dA9  dA10
dA_tot = np.array([150, 250, 250, 250, 250, 250, 250, 250, 750,  0]) # total activator added 
#                 G1  G2  G3  G4  G5  G6  G7  G8  G9  G10 G11 G12 G13  G14
G_tot = np.array([25, 15, 50, 25,  5, 15, 35, 50,  5, 25, 35, 25, 175, 175]) # total genelet added

# initial genelet states (0 = OFF, 1 = ON, -1 = BLK)
#           G1 G2 G3 G4 G5 G6 G7 G8 G9 G10 G11 G12 G13 G14
G_int_vec = [1, 1, 0, 0, -1,-1, -1,-1,-1,-1, -1, -1, -1, 0]
rR_int = [0,1000,0,0,0,0,0,0,0,0]
#dA_add = ['NA','NA','NA','NA','NA','NA','NA','NA','NA',250]
 
''' initializing initial conditions '''
I_BS_IFFL_FB.initial_conditions(dA_tot,G_tot,G_int_vec,dB_added=[0,0,250,50,150,150,150,150,0,0],rRin=rR_int) # default of 0 for all other initial conditions
I_BS_IFFL_FB.plot_topology(show_rnas=1,layout='shell',plot_title='State 1')

# re-defining the rate constants such that there is no enzymes in the beginning to set state 
# this mimics the anneal of rR2 and dA2
kpr = 0.0 # repressor production rates
kpc = 0.0 # coactivator production rates
kpi = 0.0 # inducer production rates
kd_H = 0.0 # RNaseH degradation rates
kd_A = 0.0 # RNaseA degradation rates
kga = 1e4 # activation rates
kgar = 5e3 # repression rates
kar = 1e4 # activator inhibition rates
kgb = 1e4 # free blocking rates
kgbc = 5e3 # coactivation rates
kbc = 1e4 # blocker inhibition rates
kgab = 5e3 # active blocking rates
kir = 1e4 # inducer binding rates

rates_1 = [kpr,kpc,kpi,kd_H,kd_A,kga,kgar,kar,kgb,kgbc,kbc,kgab,kir]

t_vec1 = np.linspace(0,10,1001)*3600 # seconds

''' simulating the model before turning on induce nodes  '''
I_BS_IFFL_FB.simulate(t_vec1,1,rate_constants=rates_1)

t_vec2 = np.linspace(t_vec1[-1]/3600,16,1001)*3600 # seconds

''' simulating the model  '''
I_BS_IFFL_FB.simulate(t_vec2,2) 
  
#    
# pulling out solutions for plotting from output_concentration attribute
G1 = I_BS_IFFL_FB.output_concentration['GdA1']
G3 = I_BS_IFFL_FB.output_concentration['GdA3']
G8 = I_BS_IFFL_FB.output_concentration['GdA8']
G12 = I_BS_IFFL_FB.output_concentration['GdA12']
G13 = I_BS_IFFL_FB.output_concentration['GdA13']
G14 = I_BS_IFFL_FB.output_concentration['GdA14']

rR1 = I_BS_IFFL_FB.output_concentration['rR1']
rR2 = I_BS_IFFL_FB.output_concentration['rR2']

''' Plotting solution '''
sim_t = I_BS_IFFL_FB.sol.t

fs = 13
plt.figure()
plt.suptitle('I_BS_IFFL1|2_FB1',fontsize=fs+2,weight='bold')
   
plt.subplot(3,3,1)   
plt.title('State 1 -> State 2')
plt.plot(sim_t/3600-10,G13/G_tot[12],color=[0.5,0.5,0],linewidth=4,linestyle='-')
plt.plot(sim_t/3600-10,G14/G_tot[13],color=[1,0.75,0],linewidth=4,linestyle='-')
plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
#plt.xlabel('time (hr)',fontsize=fs,weight='bold')
plt.legend(['TN9','TN10'],frameon=0)
plt.xticks(fontsize=fs,weight='bold')
plt.yticks(fontsize=fs,weight='bold')
ax1 = plt.gca()
ax1.set_xlim(0,6)
ax1.set_ylim(-0.1,1.1)
ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')

plt.subplot(3,3,4)    
plt.plot(sim_t/3600-10,G3/G_tot[2],color=[0,0.5,1],linewidth=4,linestyle='-')
plt.plot(sim_t/3600-10,G1/G_tot[0],color=[0.75,0.9,0],linewidth=4,linestyle='-')
plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
#plt.xlabel('time (hr)',fontsize=fs,weight='bold')
plt.legend(['TN1','TN2'],frameon=0)
plt.xticks(fontsize=fs,weight='bold')
plt.yticks(fontsize=fs,weight='bold')
ax1 = plt.gca()
ax1.set_xlim(0,6)
ax1.set_ylim(-0.1,1.1)
ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')

plt.subplot(3,3,7)
plt.plot(sim_t/3600-10,G12/G_tot[11],color=[0.4,0,0.4],linewidth=4)
plt.plot(sim_t/3600-10,G8/G_tot[7],color='b',linewidth=4)
plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
plt.xlabel('time (hr)',fontsize=fs,weight='bold')
plt.xticks(fontsize=fs,weight='bold')
plt.yticks(fontsize=fs,weight='bold')
plt.legend(['TN5','TN8'],frameon=0,loc=9)
ax1 = plt.gca()
ax1.set_xlim(0,6)
ax1.set_ylim(-0.1,1.1)
ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')

''' Exporting symbolic equations '''
I_BS_IFFL_FB.export_sym_eqs(file_name='I_BS_IFFL_FB_equations')

''' Compiling genelet sequences to build the network '''
input_file_loc = 'C:\\Users\\sscha\\OneDrive - Johns Hopkins University\\Desktop\\Programming Club\\Python\Python 2020\\General genelet model\\git test\\PythonVersion\\all_genelet_sequences.xlsx'
file_name='I_BS_IFFL_FB_sequences'

# compile sequences without the BTH domains where possible (in the BS module)
I_BS_IFFL_FB.compile_network_sequences(input_file_loc,desired_nodes=['G7','G8','G2','G3','G1','G5','G6','G4','G9','G10'],bth=0,save_file_name=file_name)
