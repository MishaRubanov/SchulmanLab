import numpy as np 
import matplotlib.pyplot as plt
import GeneralGeneletModel as GGM

''' 
###############################################################################
Supplementary Figure 36
This simulates an extended IFFL where an upstream genelet produces the coactivator
of the first two nodes of the IFFL with and without 5% leak
This simulates the network both with and without 150 nM dB2,dB3,dB4 which suppress the leak
###############################################################################
'''

#          G1 G2 G3 G4 G5     
act_vec =  [1, 2, 2, 3, 4]
blk_vec =  [0, 2, 2, 3, 4]
prod_vec = [2, 3, 4,-4, 0]
indc_vec = [0, 0, 0, 0, 0]

lf = 0.05 # leak fraction

''' initializing topology '''
eIFFL = GGM.GeneletNetwork(act_vec,prod_vec,indc_vec,blk_vec)
eIFFL.plot_topology(layout='circular')
eIFFL.plot_topology(layout='circular',show_rnas=0)

# Define initial conditions
dA_tot = [0]*2
#                     dA1  dA2  dA3  dA4
dA_tot[0] = np.array([250, 250, 250, 250]) # total activator added
dA_tot[1] = np.array([  0, 250, 250, 250]) # total activator added
#                 G1  G2  G3 G4  G5
G_tot = np.array([25, 5, 25, 25, 25]) # total genelet added

# initial genelet states (0 = OFF, 1 = ON, -1 = BLK)
#            G1 G2 G3 G4 G5
G_int_vec = [0, -1,-1,-1,-1]

G2 = [0]*2
G4 = [0]*2
G5 = [0]*2
ls = ['-',':']
'''
0 nM ADDED BLOCKERS
'''
plt.figure()
for i in range(len(dA_tot)):
 
    ''' initializing initial conditions '''
    eIFFL.initial_conditions(dA_tot[i],G_tot,G_int_vec) # default of 0 for all other initial conditions
    
    t_vec1 = np.linspace(0,3,1001)*3600 # seconds
    
    ''' simulating the IFFL'''
    eIFFL.simulate(t_vec1,1,leak=lf)
    
    # pulling out the desired concentrations for plotting
    
    G2[i] = eIFFL.output_concentration['GdA2']
    G4[i] = eIFFL.output_concentration['GdA4']
    G5[i] = eIFFL.output_concentration['GdA5']
    
    sim_t = eIFFL.sol.t
    
    fs = 13
    
    plt.suptitle('Extended IFFL with leak',fontsize=fs+2,weight='bold')
    plt.subplot(3,4,2)
    plt.title('No addded blockers')
    plt.plot(sim_t/3600,G2[i]/G_tot[1],color=[0,0.4,0],linewidth=4,linestyle=ls[i])
    plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
    plt.xticks(fontsize=fs,weight='bold')
    plt.yticks(fontsize=fs,weight='bold')
    ax1 = plt.gca()
    ax1.set_xlim(0,3)
    ax1.set_ylim(-0.1,1.1)
    ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
    ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')
    
    plt.subplot(3,4,6)
    plt.plot(sim_t/3600,G4[i]/G_tot[3],color=[0.6,0,0],linewidth=4,linestyle=ls[i])
    plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
    plt.xticks(fontsize=fs,weight='bold')
    plt.yticks(fontsize=fs,weight='bold')
    ax1 = plt.gca()
    ax1.set_xlim(0,3)
    ax1.set_ylim(-0.1,1.1)
    ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
    ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')
    
    plt.subplot(3,4,10)
    plt.plot(sim_t/3600,G5[i]/G_tot[4],color=[0.4,0,0.4],linewidth=4,linestyle=ls[i])
    plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
    plt.xlabel('time (hr)',fontsize=fs,weight='bold')
    plt.xticks(fontsize=fs,weight='bold')
    plt.yticks(fontsize=fs,weight='bold')
    ax1 = plt.gca()
    ax1.set_xlim(0,3)
    ax1.set_ylim(-0.1,1.1)
    ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
    ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')
    
    
G2 = [0]*2
G4 = [0]*2
G5 = [0]*2
'''
150 nM ADDED BLOCKERS
'''
for i in range(len(dA_tot)):
    ''' initializing initial conditions '''
    eIFFL.initial_conditions(dA_tot[i],G_tot,G_int_vec,dB_added=[0,150,150,150]) # default of 0 for all other initial conditions
    
    t_vec1 = np.linspace(0,3,1001)*3600 # seconds
    
    ''' simulating the IFFL'''
    eIFFL.simulate(t_vec1,1,leak=lf)
    
    # pulling out the desired concentrations for plotting
    G2[i] = eIFFL.output_concentration['GdA2']
    G4[i] = eIFFL.output_concentration['GdA4']
    G5[i] = eIFFL.output_concentration['GdA5']
    
    sim_t = eIFFL.sol.t
    
    fs = 13
    plt.suptitle('Extended IFFL with leak',fontsize=fs+2,weight='bold')
    plt.subplot(3,4,4)
    plt.title('150 nM blockers')
    plt.plot(sim_t/3600,G2[i]/G_tot[1],color=[0,0.4,0],linewidth=4,linestyle=ls[i])
    plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
    plt.xticks(fontsize=fs,weight='bold')
    plt.yticks(fontsize=fs,weight='bold')
    ax1 = plt.gca()
    ax1.set_xlim(0,3)
    ax1.set_ylim(-0.1,1.1)
    ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
    ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')
    
    plt.subplot(3,4,8)
    plt.plot(sim_t/3600,G4[i]/G_tot[3],color=[0.6,0,0],linewidth=4,linestyle=ls[i])
    plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
    plt.xticks(fontsize=fs,weight='bold')
    plt.yticks(fontsize=fs,weight='bold')
    ax1 = plt.gca()
    ax1.set_xlim(0,3)
    ax1.set_ylim(-0.1,1.1)
    ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
    ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')
    
    plt.subplot(3,4,12)
    plt.plot(sim_t/3600,G5[i]/G_tot[4],color=[0.4,0,0.4],linewidth=4,linestyle=ls[i])
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
eIFFL.export_sym_eqs(file_name='eIFFL_leak_equations')

''' Compiling genelet sequences to build the network '''
input_file_loc = 'C:\\Users\\sscha\\OneDrive - Johns Hopkins University\\Desktop\\Programming Club\\Python\Python 2020\\General genelet model\\git test\\PythonVersion\\all_genelet_sequences.xlsx'
file_name='eIFFL_leak_sequences'

eIFFL.compile_network_sequences(input_file_loc,desired_nodes=['G3','G5','G6','G4'],save_file_name=file_name)