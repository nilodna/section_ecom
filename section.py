"""
    TODO:

        insert crossSection_optimized into this class, so you don't have to
        load from another file.

        improve the visualization methods, adjusting the z and x axis
"""
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.axes_grid1 import ImageGrid
from collections import OrderedDict
import pickle
import numpy as np
import pandas as pd
import xarray as xr
import cmocean as cmo
import os

import sys
sys.path.insert(1,'../')
from Thermohaline.interp_teste import crossSection_optimized

class Section():
    def __init__(self,fname,iy,jx):
        self.fname = fname
        self.ncin  = xr.open_dataset(fname)
        self.iy    = iy
        self.jx    = jx

    def load_section(self,nstep,var='temp',mean=False):
        self.nstep = nstep
        self.var   = var

        if mean:
            self.raw_section = np.nanmean(self.ncin[var][nstep,:,self.iy,:],axis=0)
        else:
            self.raw_section = self.ncin[var][nstep,:,self.iy,:]

    def load_auxiliar_data(self):
        # getting auxiliar variables from netcdf
        self.depth= self.ncin.depth.values
        self.lon  = self.ncin.lon.values
        self.lat  = self.ncin.lat.values
        self.sigma= self.ncin.sigma.values
        self.h1   = self.ncin.h1.values
        self.angle= self.ncin.ang.values

        self.lon[self.lon == 0.] = np.nan
        self.lat[self.lat == 0.] = np.nan

    def set_newgrid(self,horizRes,vertRes,depRef,xlim):
        # setting configurations for interpolate from sigma to z-level at depth
        self.horizResolution = horizRes# this configuration give us a horizontal resolution nearest to 10 km
        self.vertResolution  = vertRes #isso significa que teremos uma resolucao vertical de 1m
        self.depRef          = depRef  #profundidade de referencia para interpolacao
        self.limiteEixoX     = xlim    # limite, em metros, do eixo X para a secao vertical

    def interp_sig2z(self):
        # self.interp_section,self.ndist,self.ndepth,self.dist2,self.sig,self.depth = self.crossSection_optimized()
        lon = self.lon[self.iy,:self.jx]
        dep = self.depth[self.iy,:self.jx]
        h1 = self.h1[self.iy,:self.jx]
        data = self.raw_section[self.iy,:self.jx]

        self.interp_section,self.ndist,self.ndepth,self.dist2,self.sig,self.depth = crossSection_optimized(lon,dep,self.sigma,h1,data,horizResolution=1000,vertResolution=80,depRef=80,ind=self.iy)
        self.xplot,self.zplot = np.meshgrid(self.ndist,self.ndepth*5000)

    def plot(self,contourf_lines=None,contour_lines=None,**kwargs):
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot()

        time_title = pd.DatetimeIndex(self.ncin.time[self.nstep].values)[0]
        self.ax.set_title(time_title.strftime("%Y-%m-%d"))

        if contourf_lines is not None:
            self.ax.contourf(self.xplot,-self.zplot,self.interp_section,contourf_lines,**kwargs)
        else:
            self.ax.contourf(self.xplot,-self.zplot,self.interp_section,**kwargs)

        if contour_lines is not None:
            self.ax.contour(self.xplot,-self.zplot,self.interp_section,contour_lines,**kwargs)

    def savefig(self,output_dir):
        self.output_dir = output_dir
        plt.savefig(self.output_dir,dpi=150)


"""
    #### Example:

    fname = fname = os.path.join('../data','model_outputs','anomalous.cdf')
    nstep = np.arange(0,2,1)
    var   = 'salt'
    can_clim = Section(fname,19,72)
    can_clim.load_section(nstep,var,mean=True)
    can_clim.load_auxiliar_data()
    can_clim.set_newgrid(1000,80,80,150000)
    can_clim.interp_sig2z()

    can_clim.plot()
    can_clim.savefig(os.path.join(outputdir, 'animation','temp_transect','clim_cananeia.png'))

    nstep = np.arange(250,255,1)
    can_anom = Section(fname,19)
    can_anom.load_section(nstep,var,mean=True)
    can_anom.load_auxiliar_data()
    can_anom.set_newgrid(1000,80,80,150000)
    can_anom.interp_sig2z()

    can_anom.plot()


"""
