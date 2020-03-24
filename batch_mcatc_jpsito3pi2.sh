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

#$ -t 24-29

i=${SGE_TASK_ID}
myrand=$[1000+$i]

inruns="sim0000"$i".dat"
outfile="jpsito3pi_sim_"$(($i+36))".root"

##numcal=21897        #+6
##numcal=22248      #+12
##numcal=22390      #+18
##numcal=22996      #+24
##numcal=23114      #+30
numcal=23114      #+36

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

