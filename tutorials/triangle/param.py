#import pyalps
import numpy as np
#import matplotlib.pyplot as plt
#import pyalps.plot
import h5py

#from h5dump import *
def create_dataset_from_dict(f, path, dict):
    dset = []
    for k,v in dict.items():
        #print k, v, type(v)
        if isinstance(v,int):
            dset.append(f.create_dataset(path+'/'+k, data=v, dtype='i4'))
        else:
            dset.append(f.create_dataset(path+'/'+k, data=v))
    return dset



#prepare the input parameters
parms = (
        {                         
              'total_steps'              : 15000000,
              'thermalization_steps'     : 15000,
              'measurement_period'       : 10,
              'model.spins'            : 2,
              'model.sites'              : 3,
              'model.beta'               : 10.0,
              'model.G0_tau_file'        : 'G0_TAU.txt',
              'G1.n_legendre'            : 50,
              'model.U'                  : 4.0,
              'G1.n_matsubara'           : 1000,
          #'SEED'                    : 100,
          #'FLAVORS'                 : 2,
          #'SITES'                   : 1,
          #'NMATSUBARA'              : 1000, 
          #'NMATSUBARA_MEASUREMENTS' : 1, 
          #'N_LEGENDRE'              : 100, 
          #'N_TAU'                   : 1000, 
          #'SWEEPS'                  : 100000000,
          #'THERMALIZATION'          : 300,
          #'MAX_TIME'                : 300,
          #'MEASUREMENT_PERIOD'      : 20,
          #'BETA'                    : 50.0,
          #'GENERAL_U_MATRIX_FILE'   : 'Uijkl.txt',
          #'G0_OMEGA'                : 'G0_OMEGA.txt',
          #'G0_TAU'                  : 'G0_TAU.txt',
          #'RECALC_PERIOD'           : 20,
          #'N_TAU_UPDATE_STATISTICS' : 100,
          #'N_MULTI_VERTEX_UPDATE'   : 1,
          #'PREFIX_OUTPUT_TIME_SERIES'   : 'input'
        }
    )


#write the input file and run the simulation
#for p in parms:
#input_file = pyalps.writeParameterFile('input', parms)
#input_file = 
with open('input.ini', 'w') as f:
    for key,value in parms.items():
        print>>f, key, ' = ', value

#pyalps.writeParameterFile('input', parms)
#print input_file
#write_parms(input_f, parms)
#res = pyalps.runDMFT(input_file)
