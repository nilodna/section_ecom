from __future__ import unicode_literals
import xarray as xr
import numpy as np
import pandas as pd
import sys
import os

# interpolation packages
from scipy import interpolate


# ---- INTERPOLATION FROM SIGMA TO Z LEVEL ---- #
def interpDistance(variable,h1,ndepth,ind=99,nlines=10000):
    """Função para interpolar as colunas de uma matrix 2D, refinando a seção
    vertical.

    variable is a 2D matrix, with x axis given by distance and y axis by depth

    Parameters
    ----------
    variable : np.ndarray
        2D array data, with section information.
    h1 : np.ndarray
        2D array with distance between each grid cell.
    cutLon : int
        Offshore limit for the section.
    ndepth : np.ndarray
        New depth vector.
    nlines : int
        New distance lenght.
    maxDist : int
        Max distance of the transect.

    Returns
    -------
    type
        Description of returned object.

    """
    # calculate distance
    dist = np.cumsum(h1)
    maxDist = int(round(np.nanmax(dist)))

    nstdl,columns = variable.shape

    # simple hack because I don't know what to yet
    # variable[:,:16] = np.nan

    variableI = np.zeros((nstdl,nlines))*np.nan # array to store interpolated data

    newDist = np.linspace(0,maxDist,nlines) # new distance axis with 1km resolution

    # considering a 2D array, we don't need to vectorize data again, like in
    # sigma2stdl, because only one for loop is enough to do this.
    # so ...
    for z in np.arange(0,nstdl,1): # reading each z level
        line = variable[z,:]

        # preciso criar uma condicao dupla, onde a distancia que quero interpolar
        # seja a partir de do primeiro indice onde nao eh nan
        nans = np.where(~np.isnan(line)) # indices de onde eh nan
        dist2interp = dist[nans[0]] # distancias sem NaN, onde a posicao [0] eh o limite inicial

        onlywater = newDist >= dist2interp[0]
        # filtrando o vetor da nova distancia, para somente aonde temos dado
        axis2interp = np.array(newDist)[onlywater]

        # criando funcao de interpolacao 1D
        finterp = interpolate.interp1d(dist[:], line,fill_value='extrapolate')
        lineI   = finterp(axis2interp)

		# stores at vectorized variable
        variableI[z, onlywater] = lineI
        variableI[z, ~onlywater] = np.NAN

    grid_x,grid_z = np.meshgrid(newDist/1000,ndepth)

    return newDist,grid_x,grid_z,variableI

def interpDepth(depth,h1,ind,ndist):
    """Interpolate bathymetry to a new distance.

    Parameters
    ----------
    depth : np.ndarray
        2D array with depth data.
    h1 : np.ndarray
        2D array with distance between each grid cell.
    ind : int
        Index for the latitude section.
    ndist : np.ndarray
        New distance array.
    cutLon : int
        Offshore limit for the section.

    Returns
    -------
    nDep : np.ndarray
        New bathymetry interpolated to ndist lenght.

    """
    dist = np.cumsum(h1) # in km

    dist2interp = dist
    dep2interp  = depth

    # inds = ndist <= np.nanmax(dist)
    axis2interp = ndist# np.array(ndist)[inds]

    fDep = interpolate.interp1d(dist2interp,dep2interp,fill_value='extrapolate')
    nDep = fDep(axis2interp)

    return nDep

def interpSigma(variable,sigma,localdep,nstdl,depRef=100):
    """Function to interpolate from sigma vertical coordinates to
    z depth coordinates, based on standard levels given by user.

    Parameters
    ----------
    variable : np.ndarray
        2D array, with distanceXsigma
    sigma : np.ndarray
        1D array with sigma levels
    localdep : np.ndarray
        1D array with depth in the section
    nstdl : int
        Number of standard levels to interpolate.

    Returns
    -------
    stdl : np.ndarray
        Standard levels created.
    variableI : np.ndarray
        2D array with interpolated data
    """
    stdl = np.linspace(0,depRef,nstdl)

    sigmalevels,columns = variable.shape
    variableI = np.zeros((nstdl,columns)) * np.nan

    for i in np.arange(0,columns,1):
        if ~np.isnan(localdep[i]):
            depsigma = -localdep[i] * sigma

            D = list(depsigma)
            D.insert(0,0)

            profile = np.zeros(sigmalevels+1)
            profile[1:] = variable[:,i]
            profile[0] = profile[1]

            watercolumn = stdl <= localdep[i]
            stdl2interp = np.array(stdl)[watercolumn]

            # interpolate to the same standard levels
            fsigma2stdl = interpolate.interp1d(D, profile)
            profileI = fsigma2stdl(stdl2interp)

            variableI[watercolumn,i] = profileI
            variableI[~watercolumn,i] = np.nan

    return stdl,variableI

def crossSection_optimized(lon,depth,sigma,h1,variable,horizResolution=10000,vertResolution=100,depRef=100,ind=99):
    """
        >> Uplot,ndist,ndepth,dist2,sig,depth = oceano.crossSection_optimized(lon,depth,sigma,h1,U,horizResolution=horizResolution,vertResolution=vertResolution,depRef=depRef,ind=ind)
    """

    # setting the new grid
    ndepth = np.linspace(0,depRef,vertResolution) # new vertical resolution
    # ndist  = np.linspace(0,100000,horizResolution # new horizontal resolution

    # interpolating horizontal/distance
    ndist,grid_x,grid_z,Tinterp = interpDistance(variable,h1,ndepth,ind=ind,nlines=horizResolution)

    # interpolatin bathymetry, to have same len of ndist
    ndep = interpDepth(depth,h1,ind,ndist)

    # interpolating vertically
    stdl,Tplot = interpSigma(Tinterp,sigma,ndep,vertResolution,depRef)

    # create grid to plot transect
    xgrid,zgrid = np.meshgrid(ndist/1000,ndepth) # with distance in km

    # creating grid for bathymetry
    s = np.tile(sigma,(len(ndist),1))
    s = np.transpose(s)
    d = np.tile(ndep,(len(sigma),1))
    sig = s*d

    dist2 = np.tile(ndist,(21,1))

    return Tplot,ndist,ndepth,dist2,sig,depth
#
# #################################################################################
# iy = 19
# jx = 72
#
# can = data['clim']['iy_can'][:,:jx]
#
# dataplot,ndist,ndepth,dist2,sig,depth = crossSection_optimized(lon,depth,sigma,h1,can,
# horizResolution=horizResolution,vertResolution=vertResolution,depRef=depRef,ind=iy)
#
# # create new grid to plot
# xgrid,zgrid = np.meshgrid(ndist,ndepth)
#
# # plot
# fig,ax = plt.subplots()
# cf,cr,cref = plot_section(ax,dataplot,xgrid,zgrid,np.arange(10,30,0.1),cm.binary,np.arange(10,30,2),[18])
