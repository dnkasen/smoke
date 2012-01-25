
# input/output files

model_file  =  "../models/lucy_test.mod"    # 3D model file to read in
output_file =  "my_lightcurve.dat"          # output light curve file

# time stepping (all times in days)

tstep_min   =   1.0          # smallest allowed time step
tstep_max   =   1.0          # largest allowed time step
tstep_del   =   0.2          # largest allowed step fraction of current time
n_times     =  1000          # maximum nuber of steps taken

t_stop      =  60.0          # time to stop calculation, in days
t_delta     =   1.0          # time spacing in ouput light curve


n_photons   =  1e4           # total number of photon packets to use
n_mu        =    1           # number of cos(theta) bins in output spectrum             
n_phi       =    1           # number of phi bins in output spectrum

opacity     =  0.1           # optical grey opacity in cm^2/g
