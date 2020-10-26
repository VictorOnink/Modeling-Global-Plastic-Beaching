#This will be a general code that will allow easy submission of large beaching simulation jobs. It is supposed to work in unison with generalBeachingScenarios.py, which is the python code that contains the actually parcels to run the different simulations

#####################################################################################
# First we define the general parameters of the run                                 #
#####################################################################################
#0=first order, 1=coastal, 2=stochastic beaching/resuspension, 3=coast type dependent
SCENARIO=2
export SCENARIO
#for scenario 1, the time a particle must be near the coast to beach (in days)
VICINITY=2
export VICINITY
#for scenario 2, the beaching and resuspension timescales (in days)
SHORETIME=10
export SHORETIME
RESUSTIME=69
export RESUSTIME
#For scenario 3, we need to indicate if beaching is more likely with sand or not sand. 0 = more sand is less likely beaching and resuspension, 1 = more sand is more likely beaching and resuspension
shoreDepen=0
export shoreDepen
#the starting year of the simulation, and how many years the simulation will take
STARTTIME=2010
export STARTTIME
#Which input distribution do we want to use? 0=Jambeck, 1=another, maybe lebreton
INPUT=1
export INPUT

START=0 #if 0, then it is a new simulation set. But this allows a crashed run to be started not entirely from scratch
SIMLEN=5 #Number of years the simulation runs

STOKES=0 #0 = include stokes, 1 = do not include stokes
export STOKES

ENSEMBLE=1 #For when we want to run multiple iterations of the same parameter setup
export ENSEMBLE

#Now, we can set the job name prefix
if [ "$SCENARIO" -eq "0" ]; then
	RUNNAMEPREFIX="FirstOrder_y="${STARTTIME}"_"
	if [ "$STOKES" -eq "1" ]; then
	    RUNNAMEPREFIX=${RUNNAMEPREFIX}"_NoStokes_"
        fi
elif [ "$SCENARIO" -eq "1" ]; then
	RUNNAMEPREFIX="Coastal_vic="${VICINITY}"_y="${STARTTIME}"_"
elif [ "$SCENARIO" -eq "2" ]; then
	RUNNAMEPREFIX="stochastic_ST="${SHORETIME}"_RT="${RESUSTIME}"_y"${STARTTIME}"_"
	if [ "$STOKES" -eq "1" ]; then
	    RUNNAMEPREFIX=${RUNNAMEPREFIX}"_NoStokes_"
	fi
elif [ "$SCENARIO" -eq "1" ]; then
	RUNNAMEPREFIX="Coastal_vic="${VICINITY}"_y="${STARTTIME}"_"
elif [ "$SCENARIO" -eq "2" ]; then
        RUNNAMEPREFIX="stochastic_ST="${SHORETIME}"_RT="${RESUSTIME}"_y"${STARTTIME}"_"
        if [ "$STOKES" -eq "1" ]; then
		RUNNAMEPREFIX=${RUNNAMEPREFIX}"NoStokes_"
        fi
elif [ "$SCENARIO" -eq "3" ]; then
        RUNNAMEPREFIX="stochasticShore_sD="${shoreDepen}"_ST="${SHORETIME}"_RT="${RESUSTIME}"_y"${STARTTIME}"_"
fi
RUNNAMEPREFIX=${RUNNAMEPREFIX}"_ENSEMBLE="${ENSEMBLE}
echo $RUNNAMEPREFIX
#####################################################################################
# Now the part where we create the submission file                                  #
#####################################################################################
queu=() #a list to which all the job id's are added for queuing
if [ "$INPUT" -eq "0" ]; then
    runlength=8
elif [ "$INPUT" -eq "1" ]; then
    runlength=3
fi
for ((run=0; run<=$runlength;run++))
do
	export run           
	for ((restartnum=$START; restartnum<$SIMLEN; restartnum++))
	do
	   export restartnum
	   runname=$RUNNAMEPREFIX"run="$run"_restart="$restartnum
	   part1="#!/bin/sh"
	   part2="#SBATCH --mail-type=begin,end,fail"
	   part3="#SBATCH --mail-user=victor.onink@climate.unibe.ch"
	   part4="#SBATCH --job-name="$runname
	   part5="#SBATCH --output="runOutput/$runname".o%j"
	   part6="#SBATCH --mem-per-cpu=6G"
           part7="#SBATCH --time=38:00:00"
	   #loading the bash and setting the environment
	   part8="source /home/ubelix/climate/vo18e689/.bash_profile"
	   part9="source /home/ubelix/climate/vo18e689/anaconda3/bin/activate py3_parcels"
	   part10='cd "/home/ubelix/climate/vo18e689/codes/Plastic-Beaching-Estimates/Beach Tests"'
	   #And now the actual running of the code
	   part11="python generalBeachingScenarios.py -p 10 -v"
	   #and now the creation of the submission file
	   for i in {1..11} 
	   do
		partGrab="part"$i
		echo ${!partGrab} >> jobsubmissionFile_${run}_${restartnum}.sh
	   done
	   if [ "$restartnum" -eq "$START" ]; then
   	      jobid=$(sbatch --parsable jobsubmissionFile_${run}_${restartnum}.sh)
#	      echo $jobid
           else
	      jobid=$(sbatch --parsable --dependency=afterok:${jobid} jobsubmissionFile_${run}_${restartnum}.sh)
	   fi
	   rm jobsubmissionFile_${run}_${restartnum}.sh
	done
done
