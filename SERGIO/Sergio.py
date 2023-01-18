#Input parameters#############################################
n_sims = 100
n_genes = 100
n_cells = 60
outlier_pars = False #{ 'p':1, 'mean':2, 'scale':1.5 }
library_pars = False #{ 'mean':13.8, 'scale':0.3 }
dropout_pars = False
##############################################################


#Dependencies
from multiprocessing import Pool
import sys
import os; os.environ["OMP_NUM_THREADS"] = "8" #fix oversubscribing problem
import numpy as np
import pandas as pd
from SERGIO.sergio import sergio


#Code
if outlier_pars or library_pars or dropout_pars:
  sim_style = 'dirty'
else:
  sim_style = 'clean'
output_dir = '{}g_{}'.format(n_genes, sim_style)
output_name = '{}{}g{}'.format(sim_style[0], int(n_genes / 100), int(n_cells/100))
interaction_file = "Interaction{}g.txt".format(n_genes)
regulator_file = 'Regs{}g_1bin.txt'.format(n_genes)

nCores = int(sys.argv[1]) #number of cores is loaded as a string but should be an integer
print(nCores)

def SERGIO_fun(curr_index): # SERGIO_fun has to be defined before pool is instanced, as it takes functions that exist from before it was initialized and not after. (i think)

  np.random.seed()
  print(curr_index)
  
  # # Simulate Clean Data _ Steady-state simulation

  sim = sergio(number_genes=n_genes, number_bins = 1, number_sc = n_cells, noise_params = 1, decays=0.8, sampling_state=15, noise_type='dpd')
  
  sim.build_graph(input_file_taregts=interaction_file, input_file_regs=regulator_file, shared_coop_state=2)
  
  sim.simulate()
  expr = sim.getExpressions()
    
  # # Add Technical Noise _ Steady-State Simulations

  if outlier_pars:
    """
    Add outlier genes
    """
    expr = sim.outlier_effect(expr, outlier_prob = outlier_pars['p'], mean = outlier_pars['mean'], scale = outlier_pars['scale'])
  
  if library_pars:
    """
    Add library size effect
    """
    libFactor, expr = sim.lib_size_effect(expr, mean = library_pars['mean'], scale = library_pars['scale'])

  """
  Convert to UMI count
  """
  count_matrix = sim.convert_to_UMIcounts(expr)

  """
  Make a 2d gene expression matrix
  """
  count_matrix = np.concatenate(count_matrix, axis = 1)

  np.savetxt('sergio_output/{}/{}_{}.csv'.format(output_dir,output_name, curr_index), count_matrix, delimiter=",") 
  print('ok')

pool = Pool(nCores)

inputs = range(1, n_sims+1)
pool.map(SERGIO_fun, inputs)
#results = [pool.apply(SERGIO_fun, args = (interaction_file, regulator_file)) for i in range(n_sims)]

pool.close()

#for index, sim in enumerate(results):
#np.savetxt("sergio_output/expr_clean_400g_{}.csv".format(index + 1), sim, delimiter=",")





