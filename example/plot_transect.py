from section import Section

# fake file just for now
fname = fname = 'fake/path/to/netcdf'

# set the timestep of intereste
nstep = np.arange(0,2,1)
# variable to plot
var   = 'salt'
# instantiate the Section class, given the i,j-indexes
can_clim = Section(fname,19,72)
# load data for the given transect (in sigma layers)
can_clim.load_section(nstep,var,mean=True)
# import other variables, as lon,depth,sigma,h1
can_clim.load_auxiliar_data()

# create new grid to interpolate
can_clim.set_newgrid(1000,80,80,150000)
# and finally, interpolate from sigma to z
can_clim.interp_sig2z()

# plot section
can_clim.plot()
# save
can_clim.savefig(os.path.join(outputdir,'temp_transect','clim_cananeia.png'))
