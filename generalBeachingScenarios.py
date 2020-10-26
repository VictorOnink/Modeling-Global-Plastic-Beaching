# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 09:47:40 2019

@author: Victor Onink

Here is one general file that will contain all the different beaching scenarios
"""
#Loading in all the relevant packages
from parcels import FieldSet, ParticleSet, JITParticle, AdvectionRK4,ErrorCode, plotTrajectoriesFile
from parcels import Variable,Field,GeographicPolarSquare, GeographicSquare
from parcels import GeographicPolar,Geographic
from parcels import rng as random
from datetime import timedelta, datetime
import numpy as np
from operator import attrgetter
import math
import time
from netCDF4 import Dataset
import os

###############################################################################
# General run setup parameters                                                #
###############################################################################
os.system('echo "Loading the general run parameters"')
#0=First Order, 1 = Coastal proximity, 2 = simple beaching/resuspension,
#3 = coasttype dependence
scenario=int(os.environ['SCENARIO'])
os.system('echo "scenario="'+str(scenario))
#for scenario 1, the time a particle needs to be near the coast to be deleted
vicinity=int(os.environ['VICINITY']) #days
os.system('echo "vicinity="'+str(vicinity))
#for scenario 2, the beaching and resuspension timescales
shoreTime,resusTime=int(os.environ['SHORETIME']),int(os.environ['RESUSTIME']) #days, days
#For scenario 3, how does the coastline dependence work? 
#0 = more sand is less likely beaching, 1 = more land is more likely beaching
shoreDepen = int(os.environ['shoreDepen'])
#for multiple sub runs, which one is being run
run=int(os.environ['run'])
os.system('echo "run="'+str(run))
#restart stage, where 0 = a completely new set up runs, and progressively higher
#values indicate the year of the simulation at which new run starts
#e.g. 1 means we start after one year of simulation, 2 after two years, etc.
restart=int(os.environ['restartnum'])
os.system('echo "restart="'+str(restart))
#starting year of the simulation
startyear= int(os.environ['STARTTIME'])
os.system('echo "starting year="'+str(startyear))
#Which input file we use. if input=0, then we use Jambeck
Input=int(os.environ['INPUT'])
os.system('echo "input="'+str(Input))
#Inclusion of Stokes Drift. 0 = included, 1 = not included
stokes =int(os.environ['STOKES'])
#The option to run multiple ensembles, which we of course don't want saved in the
#same folder since they they would overwrite each other...
ensemble =int(os.environ['ENSEMBLE'])

#To save myself a lot of hassle, we will have a server parameter:
# 0 = kup cluster, 1 = ubelix
server = 1


#%%
os.system('echo "Create the fieldset"')
"""
Loading in all the relevant data into the fieldset
"""
if server==0:
    datadirec="/alphadata04/onink/lagrangian_sim/"
    dataInputdirec="/alphadata04/onink/lagrangian_sim/BeachingSim/Input/"
elif server==1:
    datadirec="/home/ubelix/climate/shared/onink/"
    dataInputdirec="/home/ubelix/climate/shared/onink/Input/"
filenames = {'U': datadirec+"HYCOM/HYCOM*20*.nc",
             'V': datadirec+"HYCOM/HYCOM*20*.nc",
             'Ust': datadirec+"WaveWatchIIIstokes/ww3.20*_uss.nc",
             'Vst': datadirec+"WaveWatchIIIstokes/ww3.20*_uss.nc",
                }
variables = {'U': 'water_u',
             'V': 'water_v',
             'Ust': 'uuss',
             'Vst': 'vuss',
                }
dimensions = {'U':{'lat': 'lat','lon': 'lon','time': 'time','depth':'depth'},
              'V':{'lat': 'lat','lon': 'lon','time': 'time','depth':'depth'},
              'Ust':{'lat': 'latitude','lon': 'longitude','time': 'time'},
              'Vst':{'lat': 'latitude','lon': 'longitude','time': 'time'},
              }

fieldset = FieldSet.from_netcdf(filenames, variables, dimensions,allow_time_extrapolation=True)

###############################################################################
# Adding the Stokes drift to the HYCOM currents                               #
###############################################################################
if stokes==0:
    fieldset.Ust.units = GeographicPolar()
    fieldset.Vst.units = Geographic()
    fieldset=FieldSet(U=fieldset.U+fieldset.Ust,
                     V=fieldset.V+fieldset.Vst,
                     )

###############################################################################
# Adding the border current, which applies for all scenarios except for 0     #
###############################################################################
datasetBor=Dataset(dataInputdirec+'boundary_velocities_HYCOM.nc')
borU=datasetBor.variables['MaskUvel'][:]
borU[borU!=0]=borU[borU!=0]/abs(borU[borU!=0])
borV=datasetBor.variables['MaskVvel'][:]
borV[borV!=0]=borV[borV!=0]/abs(borV[borV!=0])
borMag=np.sqrt(np.square(borU)+np.square(borV))
nonzeroMag=borMag>0
borMag[borMag==0]=1
borU=np.divide(borU,borMag)
borV=np.divide(borV,borMag)
lonBor,latBor=datasetBor.variables['lon'][:],datasetBor.variables['lat'][:]
fieldset.add_field(Field('borU', borU,lon=lonBor,lat=latBor,mesh='spherical'))
fieldset.add_field(Field('borV', borV,lon=lonBor,lat=latBor,mesh='spherical'))
fieldset.borU.units = GeographicPolar()
fieldset.borV.units = Geographic()
    

###############################################################################
# Adding the horizontal diffusion                                             #
###############################################################################
kh=10 #m^2 s^-1, following Lacerda et al. (2019) and Liubertseva et al. (2018)
dataset=Dataset(datadirec+'HYCOM/HYCOM_Surface_3h_2000-01-01.nc')
uo=dataset.variables['water_u'][0,0,:,:]
lat_kh=dataset.variables['lat'][:]
lon_kh=dataset.variables['lon'][:]
kh_f=kh*np.ones(uo.shape)
kh_f[uo.mask==True]=0
fieldset.add_field(Field('Kh_zonal', kh_f,lon=lon_kh,lat=lat_kh,mesh='spherical'))
fieldset.add_field(Field('Kh_meridional', kh_f,lon=lon_kh,lat=lat_kh,mesh='spherical'))

###############################################################################
# Adding in the  land cell identifiers                                        #
###############################################################################
landID=np.load(dataInputdirec+'land_cell_identifier.npy')
fieldset.add_field(Field('landID', landID,lon=lonBor,lat=latBor,mesh='spherical'))

###############################################################################
# Distance to the shore                                                       #
###############################################################################
datasetCoast=Dataset(dataInputdirec+'distance2coast.nc')
distance=datasetCoast.variables['distance'][0,:,:]
lonD,latD=datasetCoast.variables['lon'][:],datasetCoast.variables['lat'][:]
fieldset.add_field(Field('distance2shore', distance,lon=lonD,lat=latD,mesh='spherical'))

    
###############################################################################
# and finally (for now), the coastline type                                   #
###############################################################################
if scenario==3:
    coasttype=np.load(dataInputdirec+'coastline_sand_vs_not_sand.npy')
    fieldset.add_field(Field('coasttype', coasttype,lon=lon_kh,lat=lat_kh,mesh='spherical'))

###############################################################################
# Now the periodic halo for when we go across the 180/-180 degree line        #
###############################################################################
fieldset.add_periodic_halo(zonal=True)

#%%
os.system('echo "Setting up the particleset"')
###############################################################################
# First we set the name of the outputfile, and the restart file               #
###############################################################################
if server==0:
    rootodirec="/alphadata04/onink/lagrangian_sim/BeachingSim/Output/"
elif server==1:
    rootodirec="/home/ubelix/climate/shared/onink/Output/"
if scenario==0:
    if ensemble==1:
        odirec=rootodirec+"FirstOrder/"
    else:
        odirec=rootodirec+"FirstOrder_ensemble_"+str(ensemble)+"/"
    if not os.path.exists(odirec):
        os.makedirs(odirec)
    if stokes==0:
        ofile=odirec+"FirstOrder_y="+str(startyear)+"_I="+str(Input)+"_r="+str(restart)+"_run="+str(run)+".nc"
        rfile=odirec+"FirstOrder_y="+str(startyear)+"_I="+str(Input)+"_r="+str(restart-1)+"_run="+str(run)+".nc"
    elif stokes==1:
        ofile=odirec+"FirstOrder_NoStokes_y="+str(startyear)+"_I="+str(Input)+"_r="+str(restart)+"_run="+str(run)+".nc"
        rfile=odirec+"FirstOrder_NoStokes_y="+str(startyear)+"_I="+str(Input)+"_r="+str(restart-1)+"_run="+str(run)+".nc"
elif scenario==1:
    if ensemble==1:
        odirec=rootodirec+"coastal_"+str(vicinity)+"/"
    else:
        odirec=rootodirec+"coastal_"+str(vicinity)+"_ensemble_"+str(ensemble)+"/"
    if not os.path.exists(odirec):
        os.makedirs(odirec)
    ofile=odirec+"coastalProx_v="+str(vicinity)+"_y="+str(startyear)+"_I="+str(Input)+"_r="+str(restart)+"_run="+str(run)+".nc"
    rfile=odirec+"coastalProx_v="+str(vicinity)+"_y="+str(startyear)+"_I="+str(Input)+"_r="+str(restart-1)+"_run="+str(run)+".nc"
elif scenario==2:
    if ensemble==1:
        odirec=rootodirec+"BeachTest/shore_"+str(shoreTime)+"_resus_"+str(resusTime)+"/"
    else:
        odirec=rootodirec+"BeachTest/shore_"+str(shoreTime)+"_resus_"+str(resusTime)+"_ensemble_"+str(ensemble)+"/"
    if not os.path.exists(odirec):
        os.makedirs(odirec)
    if stokes==0:
        ofile=odirec+"simpleStochastic_st="+str(shoreTime)+"_rt="+str(resusTime)+"_y="+str(startyear)+"_I="+str(Input)+"_r="+str(restart)+"_run="+str(run)+".nc"
        rfile=odirec+"simpleStochastic_st="+str(shoreTime)+"_rt="+str(resusTime)+"_y="+str(startyear)+"_I="+str(Input)+"_r="+str(restart-1)+"_run="+str(run)+".nc"
    elif stokes==1:
        ofile=odirec+"simpleStochastic_NoStokes_st="+str(shoreTime)+"_rt="+str(resusTime)+"_y="+str(startyear)+"_I="+str(Input)+"_r="+str(restart)+"_run="+str(run)+".nc"
        rfile=odirec+"simpleStochastic_NoStokes_st="+str(shoreTime)+"_rt="+str(resusTime)+"_y="+str(startyear)+"_I="+str(Input)+"_r="+str(restart-1)+"_run="+str(run)+".nc"
elif scenario==3:
    if ensemble==1:
        odirec=rootodirec+"/ShoreType_"+str(shoreDepen)+"/shore_"+str(shoreTime)+"_resus_"+str(resusTime)+"/"
    else:
        odirec=rootodirec+"/ShoreType_"+str(shoreDepen)+"/shore_"+str(shoreTime)+"_resus_"+str(resusTime)+"_ensemble_"+str(ensemble)+"/"
    if not os.path.exists(odirec):
        os.makedirs(odirec)
    ofile=odirec+"CoastTypeStochastic_dep="+str(shoreDepen)+"_st="+str(shoreTime)+"_rt="+str(resusTime)+"_y="+str(startyear)+"_I="+str(Input)+"_r="+str(restart)+"_run="+str(run)+".nc"
    rfile=odirec+"CoastTypeStochastic_dep="+str(shoreDepen)+"_st="+str(shoreTime)+"_rt="+str(resusTime)+"_y="+str(startyear)+"_I="+str(Input)+"_r="+str(restart-1)+"_run="+str(run)+".nc"

###############################################################################
# Now the particle class, which depends on the scenario                       #
###############################################################################
if scenario==0:
    class FirstOrderParticle(JITParticle):
        #First we keep track of whether a particle is on a land cell or not
        on_land = Variable('on_land', dtype=np.int32, to_write=False,
                            initial=0)
        #Now the beaching variables
        #0=open ocean, 1=beached
        beach    = Variable('beach', dtype=np.int32,
                            initial=attrgetter('beach'))
        #Finally, I want to keep track of the age of the particle
        age = Variable('age',dtype=np.int32,initial=attrgetter('age'))
        #Weight of the particle in tons
        weights= Variable('weights',dtype=np.float32,initial=attrgetter('weights'))
        #Distance of the particle to the coast
        distance= Variable('distance',dtype=np.float32,initial=0)
elif scenario==1:
    class CoastalProximityParticle(JITParticle):
        #First we keep track of how long a particle has been close to the shore
        prox = Variable('prox', dtype=np.int32, initial=attrgetter('prox'))
        #And we need to define the cutoff point for how long the particle can 
        #be close to shore
        vic = Variable('vic', dtype=np.float32,initial=vicinity, to_write=False)
        #Now the beaching variables
        #0=open ocean, 1=beached
        beach    = Variable('beach', dtype=np.int32,
                            initial=attrgetter('beach'))
        #Finally, I want to keep track of the age of the particle
        age = Variable('age',dtype=np.int32,initial=attrgetter('age'))
        #Weight of the particle in tons
        weights = Variable('weights',dtype=np.float32,initial=attrgetter('weights'))
        #Distance of the particle to the coast
        distance= Variable('distance',dtype=np.float32,initial=0)
elif scenario==2:
    class SimpleBeachingResuspensionParticle(JITParticle):
        #Now the beaching variables
        #0=open ocean, 1=beached
        beach    = Variable('beach', dtype=np.int32,
                            initial=attrgetter('beach'))
        #Now the setting of the resuspension time and beaching time
        resus_t = Variable('resus_t',dtype=np.float32,
                            initial=resusTime,to_write=False)
        coastPar = Variable('coastPar',dtype=np.float32,
                            initial=shoreTime,to_write=False)
        #Finally, I want to keep track of the age of the particle
        age = Variable('age',dtype=np.int32,initial=attrgetter('age'))
        #Weight of the particle in tons
        weights = Variable('weights',dtype=np.float32,initial=attrgetter('weights'))
        #Distance of the particle to the coast
        distance= Variable('distance',dtype=np.float32,initial=0)
elif scenario==3:
    class shoreTypeBeachingResuspensionParticle(JITParticle):
        #Now the beaching variables
        #0=open ocean, 1=beached
        beach    = Variable('beach', dtype=np.float32,
                            initial=attrgetter('beach'))
        #Now the setting of the resuspension time and beaching time
        resus_t = Variable('resus_t',dtype=np.float32,
                            initial=resusTime,to_write=False)
        coastPar = Variable('coastPar',dtype=np.float32,
                            initial=shoreTime,to_write=False)
        #Finally, I want to keep track of the age of the particle
        age = Variable('age',dtype=np.int32,initial=attrgetter('age'))
        #weight of the particle in tons
        weights = Variable('weights',dtype=np.float32,initial=attrgetter('weights'))
        #Distance of the particle to the coast
        distance= Variable('distance',dtype=np.float32,initial=0)

###############################################################################
# Setting the time parameters of the runs, which will depend on restarts, and #
# the initial particle positions                                              #
###############################################################################
starttime=datetime(startyear+restart,1,1,0,0)
endtime=datetime(startyear+restart+1,1,1,0,0)
simulationlength=endtime-starttime

###############################################################################
# Loading in the starting coordinates. In case of a restart file we make use  #
# of the restart function we use, which grabs the particle positions from the #
# previous file                                                               #
###############################################################################
def removeNan(dataset,variable,lastSelec,finalTime,lastTimeSelec):
    var=np.array(dataset.variables[variable][:])
    varSelec=var[lastSelec[0],lastSelec[1]]
    #removed particles will be indicated by the 2, so we just need to filter the particles that don't have a beach variable of 0 or 1 in the analysis. Also, by setting this to 2 the particles won't be advected in the run
    varSelec[lastTimeSelec!=finalTime]=2
    return varSelec

def restartInitialPosition(rfile,restart,scenario):
    dataset=Dataset(rfile)
    time=dataset.variables['time'][:]
    finalTime=time[0,-1]
    lastSelec=np.ma.notmasked_edges(time,axis=1)[1]
    lastTimeSelec=time[lastSelec[0],lastSelec[1]]

    lon_reset=removeNan(dataset,'lon',lastSelec,finalTime,lastTimeSelec)
    lat_reset=removeNan(dataset,'lat',lastSelec,finalTime,lastTimeSelec)
    beach_reset=removeNan(dataset,'beach',lastSelec,finalTime,lastTimeSelec)
    age_reset=removeNan(dataset,'age',lastSelec,finalTime,lastTimeSelec)
    weight_reset=removeNan(dataset,'weights',lastSelec,finalTime,lastTimeSelec)
    if scenario==1:
        prox_reset=removeNan(dataset,'prox',lastSelec,finalTime,lastTimeSelec)
    #Returning the relevant arrays    
    if scenario==1:
        return lon_reset,lat_reset,beach_reset,age_reset,prox_reset,weight_reset
    else:
        return lon_reset,lat_reset,beach_reset,age_reset,weight_reset

if restart==0:
    if Input==0:
        lons=np.load(dataInputdirec+'Jambeck2010/Jam'+str(2010)+'Lons'+str(run)+'.npy')
        lats=np.load(dataInputdirec+'Jambeck2010/Jam'+str(2010)+'Lats'+str(run)+'.npy')
        weights=np.load(dataInputdirec+'Jambeck2010/Jam'+str(2010)+'Weight'+str(run)+'.npy')
    if Input==1:
        lons=np.load(dataInputdirec+'Lebreton2010/Leb'+str(2010)+'Lons'+str(run)+'.npy')
        lats=np.load(dataInputdirec+'Lebreton2010/Leb'+str(2010)+'Lats'+str(run)+'.npy')
        weights=np.load(dataInputdirec+'Lebreton2010/Leb'+str(2010)+'Weight'+str(run)+'.npy')
    beached=np.zeros(len(lons),dtype=np.int32)
    age_par=np.zeros(len(lons),dtype=np.int32)
    if scenario==1:
        prox_par=np.zeros(len(lons),dtype=np.int32)
    repeatStep=timedelta(days=31)
else:
    if scenario==1:
        lons,lats,beached,age_par,prox_par,weights=restartInitialPosition(rfile,restart,scenario)
    else:
        lons,lats,beached,age_par,weights=restartInitialPosition(rfile,restart,scenario)
    repeatStep=None
lons[lons>180]-=360

###############################################################################
# Creating the particle set                                                   #
###############################################################################

#Creating the particle set, with an initial re-release time of 4 weeks
if scenario==0:
    pset = ParticleSet(fieldset=fieldset, pclass=FirstOrderParticle, 
                       lon=lons, lat=lats,beach=beached,age=age_par,
                       weights=weights,
                       time=starttime, repeatdt=repeatStep)
elif scenario==1:
    pset = ParticleSet(fieldset=fieldset, pclass=CoastalProximityParticle, 
                       lon=lons, lat=lats,beach=beached,age=age_par,
                       prox=prox_par,weights=weights,
                       time=starttime, repeatdt=repeatStep)
elif scenario==2:
    pset = ParticleSet(fieldset=fieldset, pclass=SimpleBeachingResuspensionParticle, 
                       lon=lons, lat=lats,beach=beached,age=age_par,weights=weights,
                       time=starttime, repeatdt=repeatStep)
elif scenario==3:
    pset = ParticleSet(fieldset=fieldset, pclass=shoreTypeBeachingResuspensionParticle, 
                       lon=lons, lat=lats,beach=beached,age=age_par,weights=weights,
                       time=starttime, repeatdt=repeatStep)

random.seed(int(time.time())*100000)

#%% 
os.system('echo "Loading the relevant kernels"')
###############################################################################
# The delete particle Kernel                                                  #
###############################################################################
def DeleteParticle(particle, fieldset, time):
    particle.delete()

###############################################################################
# The beaching kernel                                                         #
###############################################################################
if scenario==0:
    def beach(particle,fieldset,time):
        #A particle is considered beached if it is within a land cell
        if math.floor(fieldset.landID[time,particle.depth,particle.lat,particle.lon])==1:
            particle.beach=1
        #Update the age of the particle
        particle.age+=particle.dt
        
elif scenario==1:
    def beach(particle,fieldset,time):
        if particle.beach==0:        
            dist=fieldset.distance2shore[time, particle.depth, particle.lat, particle.lon]
            #If a particle is within 10 km of the shore
            if dist<10:
                particle.prox+=particle.dt
            else:
                particle.prox=0.
            if particle.prox>86400*particle.vic:
                particle.beach=1
        #Update the age of the particle
        particle.age+=particle.dt

elif scenario==2:
    def beach(particle,fieldset,time):
        if particle.beach==0:
            dist=fieldset.distance2shore[time, particle.depth, particle.lat, particle.lon]
            if dist<10: 
                beach_prob=math.exp(-particle.dt/(particle.coastPar*86400.))
                if random.random(0.,1.)>beach_prob:
                    particle.beach=1
        #Now the part where we build in the resuspension
        elif particle.beach==1:
            resus_prob=math.exp(-particle.dt/(particle.resus_t*86400.))
            if random.random(0.,1.)>resus_prob:
                particle.beach=0
        #Update the age of the particle
        particle.age+=particle.dt

elif scenario==3:
    if shoreDepen==0:
        #Beaching and resuspension is considered less likely if the coastline
        #is more sandy. Correction factor is (1-0.25*(1-s))=0.75+0.25*s. Thus,
        #if the coastline is more sandy (s~1), the beaching and resuspension
        #timescales become larger
        def beach(particle,fieldset,time):
            if particle.beach==0:
                dist=fieldset.distance2shore[time, particle.depth, particle.lat, particle.lon]
                #s=sandiness
                s=fieldset.coasttype[time, particle.depth, particle.lat, particle.lon]
                if dist<10:
                    # beach_prob=math.exp(-particle.dt/(particle.coastPar*(0.75+0.25*s)*86400.))
                    beach_prob=math.exp(-particle.dt/(particle.coastPar*86400.))
                    if random.random(0,1.)>beach_prob:
                        particle.beach=1
            #Next the resuspension part
            elif particle.beach==1:
                s=fieldset.coasttype[time, particle.depth, particle.lat, particle.lon]
                resus_prob=math.exp(-particle.dt/(particle.resus_t*(0.75+0.25*s)*86400.))
                if random.random(0,1.)>resus_prob:
                    particle.beach=0
            #Update the age of the particle
            particle.age+=particle.dt
    elif shoreDepen==1:
        #Beaching and resuspension is considered more likely if the coastline
        #is more sandy. Correction factor is 1-0.25*s. Thus, if the coastline is
        #more sandy, then the beaching and resuspension timescales become smaller
        def beach(particle,fieldset,time):
            if particle.beach==0:
                dist=fieldset.distance2shore[time, particle.depth, particle.lat, particle.lon]
                #s=sandiness
                s=fieldset.coasttype[time, particle.depth, particle.lat, particle.lon]
                if dist<10:
                    # beach_prob=math.exp(-particle.dt/(particle.coastPar*(0.25+0.75*s)*86400.))
                    beach_prob=math.exp(-particle.dt/(particle.coastPar*86400.))
                    if random.random(0,1.)>beach_prob:
                        particle.beach=1
            #Next the resuspension part
            elif particle.beach==1:
                s=fieldset.coasttype[time, particle.depth, particle.lat, particle.lon]
                resus_prob=math.exp(-particle.dt/(particle.resus_t*(0.25+0.75*s)*86400.))
                if random.random(0,1.)>resus_prob:
                    particle.beach=0
            #Update the age of the particle
            particle.age+=particle.dt
        

###############################################################################
# The advection kernel, which is just the default for scenario=0              #
###############################################################################
def AntiBeachNudging(particle,fieldset,time):
    """    
    The nudging current is 1 m s^-1, which ought to be sufficient to overpower
    any coastal current (I hope) and push our particle back out to sea so as to
    not get stuck
    
    update 11/03/2020: Following tests and discussions with Cleo, the nudging 
    current will now kick in starting at 500m from the coast, since otherwise 
    the particles tended to get stuck if we used the velocity treshhold. 
    """
    d1=particle.depth
    if fieldset.distance2shore[time,d1,particle.lat,particle.lon]<0.5:
        borUab,borVab=fieldset.borU[time, d1, particle.lat, particle.lon],fieldset.borV[time, d1, particle.lat, particle.lon]
        particle.lon-=borUab*particle.dt
        particle.lat-=borVab*particle.dt
        
def AdvectionRK4_floating(particle, fieldset, time):
    """Advection of particles using fourth-order Runge-Kutta integration.

    Function needs to be converted to Kernel object before execution
        
    A particle only moves if it has not beached (rather obviously)
    """
    if particle.beach==0:
        particle.distance=fieldset.distance2shore[time,particle.depth,particle.lat,particle.lon]
        d2=particle.depth
        if particle.lon>180:
            particle.lon-=360
        if particle.lon<-180:
            particle.lon+=360
        (u1, v1) = fieldset.UV[time, d2, particle.lat, particle.lon]
        lon1, lat1 = (particle.lon + u1*.5*particle.dt, particle.lat + v1*.5*particle.dt)
        if lon1>180:
            lon1-=360
        if lon1<-180:
            lon1+=360
        (u2, v2) = fieldset.UV[time + .5 * particle.dt, d2, lat1, lon1]
        lon2, lat2 = (particle.lon + u2*.5*particle.dt, particle.lat + v2*.5*particle.dt)
        if lon2>180:
            lon2-=360
        if lon2<-180:
            lon2+=360
        (u3, v3) = fieldset.UV[time + .5 * particle.dt, d2, lat2, lon2]
        lon3, lat3 = (particle.lon + u3*particle.dt, particle.lat + v3*particle.dt)
        
        if lon3>180:
            lon3-=360
        if lon3<-180:
            lon3+=360
        (u4, v4) = fieldset.UV[time + particle.dt, d2, lat3, lon3]
        
        particle.lon += (u1 + 2*u2 + 2*u3 + u4) / 6. * particle.dt
        if particle.lon>180:
            particle.lon-=360
        if particle.lon<-180:
            particle.lon+=360
        particle.lat += (v1 + 2*v2 + 2*v3 + v4) / 6. * particle.dt

###############################################################################
# The diffusion kernel                                                        #
###############################################################################
def BrownianMotion2D(particle, fieldset, time):
    """Kernel for simple Brownian particle diffusion in zonal and meridional direction.
    Assumes that fieldset has fields Kh_zonal and Kh_meridional
    we don't want particles to jump on land and thereby beach"""
    if particle.beach==0:
        r = 1/3.
        kh_meridional = fieldset.Kh_meridional[time, particle.depth, particle.lat, particle.lon]
        lat_p = particle.lat + random.uniform(-1., 1.) * math.sqrt(2*math.fabs(particle.dt)*kh_meridional/r)
        kh_zonal = fieldset.Kh_zonal[time, particle.depth, particle.lat, particle.lon]
        lon_p = particle.lon + random.uniform(-1., 1.) * math.sqrt(2*math.fabs(particle.dt)*kh_zonal/r)
        particle.lon=lon_p
        particle.lat=lat_p

###############################################################################
# The random initial input kernel                                             #
###############################################################################
def initialInput(particle,fieldset,time):
    """
    Since we have many instances that particles start at the very same position,
    when a particle is first added to the simulation it will get a random kick
    that moves it slightly away from the initial position, and so with multiple
    particles we will see a sort of star pattern around the central position.
    However, the particle shouldn't be put on land though, so a particle will
    have 100000 attempts to be placed in the simulation within a cell that is not
    land. If it doesn't work after 100000 attempts, then the particle just ends up
    starting from the unchanged position. We also set it so that the initial
    position of the particle is always closer to land than the original position
    
    Note: Tests show that at most we get at most around 7 or 8 particles per 
    release getting placed immediately on land. this varies a bit
    
    """
    if particle.age==0:
        #The low latitudes/equatorial regions have larger grid sizes
        if math.fabs(particle.lat)<40.:
            check=0
            distCur=fieldset.distance2shore[time,particle.depth,particle.lat,particle.lon]
            while check<100000:
                potLat=particle.lat+random.uniform(-0.08, 0.08)
                potLon=particle.lon+random.uniform(-0.08, 0.08)
                potLand=math.floor(fieldset.landID[time,particle.depth,particle.lat,particle.lon])
                distPot=fieldset.distance2shore[time,particle.depth,potLat,potLon]
                if potLand==0 and distPot<=distCur:
                    check+=100001
                    particle.lat=potLat
                    particle.lon=potLon
                check+=1
        #Higher latitudes above 40 degrees
        else:
            check=0
            distCur=fieldset.distance2shore[time,particle.depth,particle.lat,particle.lon]
            while check<100000:
                potLat=particle.lat+random.uniform(-0.04, 0.04)
                potLon=particle.lon+random.uniform(-0.04, 0.04)
                potLand=math.floor(fieldset.landID[time,particle.depth,particle.lat,particle.lon])
                distPot=fieldset.distance2shore[time,particle.depth,potLat,potLon]
                if potLand==0 and distPot<=distCur:
                    check+=100001
                    particle.lat=potLat
                    particle.lon=potLon
                check+=1

###############################################################################
# And now the overall kernel                                                  #
###############################################################################
if scenario==0:
    totalKernel=pset.Kernel(initialInput)+pset.Kernel(AdvectionRK4_floating)+pset.Kernel(BrownianMotion2D)+pset.Kernel(beach)
else:
    totalKernel=pset.Kernel(initialInput)+pset.Kernel(AdvectionRK4_floating)+pset.Kernel(BrownianMotion2D)+pset.Kernel(AntiBeachNudging)+pset.Kernel(beach)
    
#%%
os.system('echo "The actual calculation of the particle trajectories"')
###############################################################################
# finally the execution                                                       #
###############################################################################
pfile = pset.ParticleFile(name=ofile,
                          outputdt=timedelta(hours=24))

pset.execute(totalKernel,
             runtime=timedelta(days=simulationlength.days),
             dt=timedelta(minutes=10),
             recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle},
             output_file=pfile
             )
