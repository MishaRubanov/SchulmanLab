
''' 
###############################################################################
GeneralGeneletModel
###############################################################################
General genelet model functions are helper functions for the GeneletNetwork class 

-------------------------------------------------------------------------------------------------------------------------
CLASS: GeneletNetwork(act_vec=req1,prod_vec=req1,indc_vec=0,blk_vec=inferred,top_mat= req2,Itop_mat=0,genelet_type='KWG')
-------------------------------------------------------------------------------------------------------------------------
# this defines the genelet network to be simulated
  --------------------
   Input definitions:
  --------------------
     For act_vec and blk_vec: 
         List the length of ind_nodes (total nodes)
         Numbers represent which orthogonal activators/blockers correspond to which nodes
         0 indicates no activator/blocker for a given node
         blk_vec is optional and if not provided will be inferred from nodes that have coactivators coming to them
     For prod_vec:
         List the length of ind_nodes (total nodes)
         -ve values are repressors
         +ve values are coactivators
           0 = no production of an inducer RNA
     For indc_vec:
         List the length of ind_nodes (total nodes)
         Indicates which nodes produce inducer RNAs in the network
         -ve values are inducers that bind repressors (thereby activating) 
         +ve values are inducers that bind coactivators (thereby BLKing)
           0 = no production of an inducer RNA
          indc_vec is optional is default to zeros if not provided
     Numbers represent which orthogonal nodes the other nodes connect to
     For genelet_type:
         The type of genelet chemistry used
         'KWG' = Kim and Winfree genelet design with free activators and blockers
         'STG' = Spatiotemporal genelets with activators directly attached to genelets
                 This model also allows free activator to be added as a threshold
                 the free activator binds RNA but does not activate the genelets
    
    For top_mat:
        A trinary topology matrix that defines how topological nodes connect to one another
        A list of lists (ortho_nodes x ortho_nodes) where entries are defined as for prod_vec
        EX:  IFFL network     EX:   TSN network
             TN1  TN2  TN3         TN1  TN2  TN3
        TN1 [ 0    0    0 ]   TN1 [ 0   -1   -1 ]
        TN2 [ 1    0    0 ]   TN2 [-1    0   -1 ]
        TN3 [ 1   -1    0 ]   TN3 [-1   -1    0 ]
   
    For Itop_mat:
        A trinary topology matrix that defines how topological nodes connect to RNAs in the network
        Defines which nodes produce which inducer RNAs
        Values are defined as for indc_vec
        This is an optional input and will default to 0s if not implemented
        
    There are two ways to initialize the GeneletNetwork
    1) input BOTH act_vec and prod_vec (this is the way Sam's initial code was set up)
       
        EX: model = GeneletNetwork(act_vec,prod_vec,indc_vec,blk_vec)
        EX: model = GeneletNetwork(act_vec,prod_vec)
        EX: model = GeneletNetwork(act_vec,prod_vec,blk_vec=[0,2,2,3])
    
    2) input top_mat (this is the way Misha's sets up a topology matrix)
    
        EX: model = GeneletNetwork(top_mat=topology_mat)
        EX: model = GeneletNetwork(top_mat=topology_mat,Itop_mat=I_mat)
        
  --------------------
   Output definitions:
  --------------------
        After initializing the GeneletNetwork class the user has access to the following attributes
            model.ortho_nodes     - number of orthogonal genelet nodes
            model.ind_nodes       - number of total (individual) genelet nodes
            model.act_vec         - activator vector mapping ind_nodes to orthogonal activators
            model.prod_vec        - RNA production vector mapping ind_nodes to which RNAs they produce (repressors/coactivators)
            model.blk_vec         - blocker vector mapping ind_nodes to orthogonal blockers
            model.indc_vec        - RNA inducer production vector mapping ind_nodes to which inducers they produce
            model.topology_mat    - topology matrix (an alternative representation of act_vec/prod_vec notation)
            model.I_topology_mat  - inducer topology matrix (an alternative representation of act_vec/indc_vec notation)
            model.act_mat         - activator connection matrix used in ODEs to map activators to ind_nodes
            model.rep_mat         - repressor connection matrix used in ODEs to map repressors to ind_nodes
            model.blk_mat         - blocker connection matrix used in ODEs to map blockers to ind_nodes         
            model.ca_mat          - coactivator connection matrix used in ODEs to map coactivators to ind_nodes
            model.Rprod_mat       - repressor production matrix used in ODEs to map which ind_nodes produce which repressors
            model.Cprod_mat       - coactivator production matrix used in ODEs to map which ind_nodes produce which coactivators
            model.Rindc_mat       - repressor inducer production matrix used in ODEs to map which ind_nodes produce which repressor inducers
            model.Cindc_mat       - coactivator inducer production matrix used in ODEs to map which ind_nodes produce which coactivator inducers
            
###############################################################################
- FUNCTIONS within GenletNetwork CLASS
###############################################################################
-----------------------------------------------------------------------------------------
 -initial_conditions(dA_tot,G_tot,G_int_vec=[],dB_added=0,rRin=0,rCin=0,rIrin=0,rIcin=0)
-----------------------------------------------------------------------------------------
     # Sets the initial conditions for the network
    -------------------
     Input definitions
    -------------------
       For G_int_vec:
           List the length of ind_nodes (total nodes)
           Represents the intial state of a node
           1 for ON
          -1 for BLK
           0 for OFF
           G_int_vec is optional and if not provided it will be assumed from the network definition:
               Nodes with no inputs are ON
               Nodes with coactivator inputs are BLK
               Nodes that are only repressed are ON
       For dA_tot and G_tot:
           numpy arrays the length of ortho_nodes and ind_nodes, respectively
           Represent total concentrations of activators and genelets, respectively
       For all other species:
           Default values are 0 
           Just input a list of concentrations for any other species by name
           
     EX: model.initial_conditions(dA_tot,G_tot,G_int_vec,dB_added=[0,150,150])
     EX: model.initial_conditions(dA_tot,G_tot,G_int_vec,rRin=[1000,0,0])
    
    -------------------
     Output definitions
    -------------------
        After calling initial_conditions the user will have access to the following attributes:
            model.G_int_vec  -the list of initial genelet states (these may be determined as a default)
            model.ind_cond   -the vector holding the initial conditions for each species
            model.(species)  -all of the individual species concentrations that make up int_cond are also accessible as attributes by name
            
---------------------------------------------------------------------------------------------------------
 -simulate(t_vec,iteration,rnase='RnH',leak=0,rate_constants=[],rR=[],dA=[],rC=[],dB=[],rIr=[],rIc=[],dR=[]) 
--------------------------------------------------------------------------------------------------------- 
     # Simulates the network with the above initial conditions
     # Simulate can be called numerous times and specific conditions can be updated by name
       Use 'NA' for values that should not be updated for a given species
       For updating a species conditions, the list has to be the length of the total number of that species
       dR is a special option for KWG models where DNA repressors can be added to the reaction to permanently remove activators
          dR species are currently not included in the symbolic equations that can be output as they are not really part of the networks
       t_vec2 should start at the last timepoint from t_vec1 and go to a later 
    -------------------
     Input definitions
    -------------------
          For rnase:
              selects the rnase to use in the simulations
              'RnH' = RNase H (degrade RNA in RNA:DNA duplex)
              'RnA' = RNase A (degrade ssRNA and RNA in RNA:DNA duplex (likely with slower rate than RNase H))
              'Both' = RNase H and A
              This input is optional and 'RnH' is the default
          For leak:
              represents the leak transcription rate of BLK (KWG) or GrR (STG) 
              genelets as a fraction of the ON transcription rate
              between 0 and 1
              can be entered as a list for different leak rates for each RNA production
                 leak = [0.05,0.1,0.075,0.06]
	      There should be ind_node # of leak inputs as in this model each individual genelet
                 has its own production rate even if two genelets produce the same RNA
                       
          For rate_constants:
              import different rate constants than the default values
              If genelet_type = 'KWG' then import as a list of single values (1/M-s)
                 [kpr,kpc,kpi,kd_H,kd_A,kga,kgar,kar,kgb,kgab,kgbc,kbc,kir]
              If genelet_type = 'STG' then import as a list of single values (1/M-s)
                 [kpr,kpi,kd_H,kd_A,kd_Ag,kd_Hg,kgar,kar,kir]
              Each rate can be a list of length ind_nodes for production rates
                 Or a list of length ortho_nodes for other rates
                 Or a single value which will be assumed to be the same for all nodes
    
     EX: model.simulate(t_vec1,1) 
     EX: model.simulate(t_vec2,2,rIr=['NA','NA',10000]) 
     EX: model.simulate(t_vec3,3,rR=[1000,'NA','NA']) 
     
     EX: model.simulate(t_vec1,1,rnase='RnA')  # use RNase A 
     
     EX: model.simulate(t_vec1,1,rnase='both',leak=0.1) # use both RNase A and H 
                                                        # 10% leak transcription
                                                        
     EX: model.simulate(t_vec1,1,rnase='both',leak=0.1,rate_constants=rate_list) 
                                                        # use both RNase A and H 
                                                        # 10% leak transcription for all RNAs
                                                        # user input rate constants
                                                        
     EX: model.simulate(t_vec1,1,leak=[0.05,0.1,0.075]) # different leak rates for 3 RNAs
                                                        
     ------------------
     Output definitions
     ------------------
         After running simulate the user has access to the following useful attributes:
             model.sol    -sol.y contains all of the concentration date and sol.t contains the time in seconds
             model.(rate_constants)   -all of the individual rate constants can be accessed by name
             model.output_concentration    -the dictionary containing all of the concentrations with specific name keys
                               
      EX: model.output_concentration['GdA1'] # for KWG 
      EX: model.output_concentration['Gon1'] # for STG
      EX: model.output_concentration['rR4']                          

-------------------------------------------------------------------------------
 -export_sym_eqs(file_name='sym_eqs',path_name='')
-------------------------------------------------------------------------------
      # exports the symbolic ODE equations of the model to a text file
      # both of the inputs are optional and have default values    
    -------------------
     Input definitions
    -------------------
          For file_name:
              user defined name of the .txt file (do not include .txt in the name)
          For path_name:
              user defined path to where the file should be saved (include \ at the end)
              if user does not supply a path the file is saved in current folder
              
      EX: model.export_sym_eqs(file_name='TSN_eqs',path_name='C:\\Desktop\\')
    
    ------------------
    Output definitions
    ------------------
        The output from this is the text file with the symbolic equations in it
        with the specified or default filename in the specified or current (default) folder
    
-------------------------------------------------------------------------------
 -plot_topology(pos=[],layout='spring',plot_title='',show_rnas=1)
-------------------------------------------------------------------------------
      # plots the network topology with initial node states
      # this considers both genelet nodes and RNA coactivators/repressors as topological nodes
        which allows inducer RNA connections to be represented
      # If this is called right after initializing the network the initial states will be inferred
      # If this is called after the initial_conditions function has been called the user defined initial states will be displayed
    -------------------
     Input definitions
    -------------------
          For pos:
              This is an optional input of each nodes position in the plot
              It has to be a dictionary with keys that are the names of the topological nodes
              and the desired [x,y] position values
          For layout:
              This is an optional input that selects the algorithm to organize the nodes in the plot
              Default is 'spring'
              Other options:
                  'spectral'
                  'circular'
                  'shell'
              See the networkX package layout options for further details
          For plot_title:
              This is an optional title for the topology plot
	  For show_rnas:
	      This is an option to show the RNA repressors/coactivators in the topology as their own nodes
              This is necessary if you want to plot topologies where there are nodes that produce inducers
		if this is set to 0 then the RNAs will not be plotted and any inducer nodes will appear free standing
              
      EX: model.plot_topology(layout='shell',plot_title='IFFL network')

    ------------------
    Output definitions
    ------------------
        The output from this is the topology plot
        model.net_edges gives the user access to all the network edges used in the plot

-------------------------------------------------------------------------------
 -compile_network_sequences(input_file_loc,desired_nodes=[],min_design=0,bth=1,save_file_name='genelet_seqs',save_path_name='')
-------------------------------------------------------------------------------
      # exports genelet sequences to a text file   
      # by default this function first populates nodes that need to both be coactivated and repressed 
        it then populates other nodes that only need to be coactivated or repressed - selecting the first sequences it finds that meet the desired criteria
            i.e. if a node only needs to be coactivated but there are still some sequences that can be coact and repressed these can be selected
    -------------------
     Input definitions
    -------------------
          For input_file_loc:
	      the file path to the Excel file holding all of the viable genelet node sequences
	  For desired_nodes:
	      A list of length ortho_nodes with the names of nodes to be used to compile the sequences (['G1','G5','G8'])
              By default this is empty so the algorithm will select its own nodes
	  For min_design:
	      If min_design = 1 the algorithm will select nodes with the minimal possible function in the network
                 for example, if a node is only coactivated in the network then the algorithm will use the first node that can only be coactivated at this position
 		 Currently if there are not enough minimal nodes to fill the network this will not work
              By default this is 0 so the first nodes that meet the needs of the network will be used 
          For bth:
	      This is 1 by default which keeps the 8 base 5' blocker toehold (BTH) on all the genelets in the network
	      This can be set to 0 and then the 8 base 5' blocker toehold will be removed from all genelets that are only repressed in a given network as these nodes do not need a blocker
                 this makes it possible to order shorter -nt sequences without having to truncate their 3' end (as for the TSN)
	  For save_file_name:
              user defined name of the .txt file to be saved (do not include .txt in the name)
          For save_path_name:
              user defined path to where the file should be saved (include \ at the end)
              if user does not supply a path the file is saved in current folder
              
      EX: model.compile_network_sequences(input_file_loc,save_file_name='TSN_sequences',save_path_name='C:\\Desktop\\')
    
    ------------------
    Output definitions
    ------------------
        The output from this is the text file with the genelet sequences in it
        Each node is comprised of G-nt / G-t / activator / blocker sequences as needed for the simulated network
        


#########################################################################################################################
TROUBLE SHOOTING - common bugs I have found before that may be the source of incorrect functionality
#########################################################################################################################
-Check model.blk_vec if you are having issues and you have not supplied the blk_vec yourself
	I have previously fixed two bugs related to populating a blk_vec based off topology alone so there may still be issues here for networks I havent tested

-When doing iterative model.simulate calls make sure the t_vec units are correct, ie if t_vec1 = linspace(0,5,1000)*3600 then t_vec2 = linspace(t_vecl[-1]/3600,10,1000)*3600

-The model currently cannot handle dummy nodes that have no unique connections
                G1 G2 G3
  EX: act_vec = [1, 2, 2]
     prod_vec = [2,-2, 0]

 Here G3 does not possess a new activator and does not connect to anything - delete it
 Such a node would never be present in the simplified toplogy mat format
  so things like the blk_vec will be inferred incorrectly from the topology_mat

-When defining a network from the topology matrix, the production vector will always be defined in increasing numeric order. 
 This could cause confusion when going back and forth between act_vec/prod_vec and top_mat definitions
 For example you might define the IFFL as having prod_vec = [3,2,-3,0] which will produce a top_mat = [[0,0,0],[1,0,0],[1,-1,0]]
 but coverting that top_mat into an act_vec and prod_vec will produce prod_vec = [2,3,-3,0]. 
 Both prod_vec definitions are correct since the first 2 nodes in the IFFL share the same activator but the values of G_tot need to be changed accordingly
