#!/bin/bash
# -------------------------------------------
# --          use bash                     --
#$ -S /bin/bash
# -------------------------------------------
# --             batch name                --
#$ -N jpsito3pi
# -------------------------------------------
# --     What to redirect to where         --
#$ -cwd
#$ -o $JOB_NAME.o$JOB_ID
#$ -j y
# -------------------------------------------
# --             Enviroment                --
#$ -v PATH=$PATH:$HOME/release/KdRunFastMon,LD_LIBRARY_PATH=/usr/local/root/lib/root:/home/alexbarn/release/lib,KDBHOST=bison-2
# -------------------------------------------
# --             Queue list                --
#$ -soft
#$ -l time=24:00:00
#$ -q remote
##$ -q extralong
#$ -m beas
#$ -M ovtin.ivan@gmail.com

##$ -t 24-29
#$ -t 51-60

i=${SGE_TASK_ID}
myrand=$[1000+$i]

inruns="sim0000"$i".dat"
outfile="jpsito3pi_sim_"$i".root"

if [ $i == 24 ]; then
numcal=21897
echo "numcal=" "$numcal"
fi

if [ $i == 25 ]; then
numcal=22248
echo "numcal=" "$numcal"
fi

if [ $i == 26 ]; then
numcal=22390
echo "numcal=" "$numcal"
fi

if [ $i == 27 ]; then
numcal=22996
echo "numcal=" "$numcal"
fi

if [ $i == 28 ]; then
numcal=23114
echo "numcal=" "$numcal"
fi

if [ $i == 29 ]; then
numcal=23114
echo "numcal=" "$numcal"
fi

if [ $i == 30 ]; then
numcal=22390
echo "numcal=" "$numcal"
fi

if [ $i == 31 ]; then
numcal=22390
echo "numcal=" "$numcal"
fi

if [ $i == 32 ]; then
numcal=23114
echo "numcal=" "$numcal"
fi

if [ $i == 33 ]; then
numcal=23114
echo "numcal=" "$numcal"
fi

if [ $i == 34 ]; then
numcal=22390
echo "numcal=" "$numcal"
fi

if [ $i == 35 ]; then
numcal=22390
echo "numcal=" "$numcal"
fi

if [ $i == 36 ]; then
numcal=22390
echo "numcal=" "$numcal"
fi

if [ $i == 37 ]; then
numcal=22390
echo "numcal=" "$numcal"
fi

if [ $i == 38 ]; then
numcal=22390
echo "numcal=" "$numcal"
fi

if [ $i == 39 ]; then
numcal=21897
echo "numcal=" "$numcal"
fi

if [ $i == 40 ]; then
numcal=22248
echo "numcal=" "$numcal"
fi

if [ $i == 41 ]; then
numcal=22390
echo "numcal=" "$numcal"
fi

if [ $i == 42 ]; then
numcal=22390
echo "numcal=" "$numcal"
fi

if [ $i == 43 ]; then
numcal=22390
echo "numcal=" "$numcal"
fi

if [ $i == 44 ]; then
numcal=22390
echo "numcal=" "$numcal"
fi

if [ $i == 45 ]; then
numcal=22390
echo "numcal=" "$numcal"
fi

if [ $i == 46 ]; then
numcal=22390
echo "numcal=" "$numcal"
fi

if [ $i == 47 ]; then
numcal=22390
echo "numcal=" "$numcal"
fi

if [ $i == 48 ]; then
numcal=22390
echo "numcal=" "$numcal"
fi

if [ $i == 49 ]; then
numcal=22390
echo "numcal=" "$numcal"
fi

if [ $i == 50 ]; then
numcal=22390
echo "numcal=" "$numcal"
fi


if [ $i == 51 ]; then
numcal=22390
echo "numcal=" "$numcal"
fi
if [ $i == 52 ]; then
numcal=22390
echo "numcal=" "$numcal"
fi
if [ $i == 53 ]; then
numcal=22390
echo "numcal=" "$numcal"
fi
if [ $i == 54 ]; then
numcal=22390
echo "numcal=" "$numcal"
fi
if [ $i == 55 ]; then
numcal=22390
echo "numcal=" "$numcal"
fi
if [ $i == 56 ]; then
numcal=22390
echo "numcal=" "$numcal"
fi
if [ $i == 57 ]; then
numcal=22390
echo "numcal=" "$numcal"
fi
if [ $i == 58 ]; then
numcal=22390
echo "numcal=" "$numcal"
fi
if [ $i == 59 ]; then
numcal=22390
echo "numcal=" "$numcal"
fi
if [ $i == 60 ]; then
numcal=22390
echo "numcal=" "$numcal"
fi


##numcal="numcal_"$i

echo "inruns=" "$inruns"
echo "outfile=" "$outfile"
echo "numcal=" "$numcal"

#inpatt='/spool/users/skononov/runs/psi2s/tt*.nat'
#inruns="sim000013.dat"
#inruns="sim000014.dat"
#inruns="sim000015.dat"
#inruns="sim000016.dat"
#inruns="sim000017.dat"
#inruns="sim000018.dat"
#outfile="/spool/users/skononov/runs/psi2s/tt_50runs.root"

#start the job
##$HOME/development/KEDR/bin/mktree2t -b -s "$inpatt" -o $outfile $inruns
#
$HOME/development/Jpsi_to_pipipi0/analysis_j-psi_pi_pi_pi0  -v $numcal -o $outfile $inruns
##$HOME/development/Jpsi_to_pipipi0/analysis_j-psi_pi_pi_pi0  -v $numcal -n 10000 -o $outfile $inruns
##$HOME/development/Jpsi_to_pipipi0/analysis_j-psi_pi_pi_pi0  -v 21897 -n 20000 -o jpsi_rho0pi0_events_sim_test.root $inruns
##$HOME/development/Jpsi_to_pipipi0/analysis_j-psi_pi_pi_pi0  -v 21897 -o jpsi_rho0pi0_events_sim.root $inruns
##$HOME/development/Jpsi_to_pipipi0/analysis_j-psi_pi_pi_pi0  -v 22248 -o jpsi_rho-pi+rho+pi-_events_sim.root $inruns
##$HOME/development/Jpsi_to_pipipi0/analysis_j-psi_pi_pi_pi0  -v 22370 -o jpsi_rho+pi-rho-pi+_events_sim.root $inruns
##$HOME/development/Jpsi_to_pipipi0/analysis_j-psi_pi_pi_pi0  -v 22803 -o jpsi_rho0pi0_events_sim_2.root $inruns
##$HOME/development/Jpsi_to_pipipi0/analysis_j-psi_pi_pi_pi0  -v 22803 -o jpsi_rho0pi0_events_sim_2.root $inruns
##$HOME/development/Jpsi_to_pipipi0/analysis_j-psi_pi_pi_pi0  -v 22996 -o jpsi_rho+pi-rho-pi+_events_sim_2.root $inruns

status=$?
if [ $status != 0 ]; then
  echo "Program exited with status $status"
  exit
fi

echo "Job finished\n"

