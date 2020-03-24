# The code for analysis J/psi->pi+pi-pi0 process
Get the code: <br />
```
git clone https://github.com/ivovtin/Jpsi_to_pipipi0
```
Example for run: <br />
```
analysis_j-psi_pi_pi_pi0 -x /space/runs/daq022928.nat.bz2
analysis_j-psi_pi_pi_pi0  -v 23110 -n 50 /spool/users/ovtin/sim000018.dat - simulation 50 events
```
For run on batch need use next line:<br />
```
bsub batch_dataatc_jpsito3pi.sh
bsub batch_mcatc_jpsito3pi.sh
```

Launch KDisplay for view event with reconstruction: <br />
```
bzcat /space/runs/daq021913.nat.bz2 | KDisplay -r -e3197
```
Launch simulation in KDisplay: <br />
```
/home/ovtin/development/bin/KDisplay < /spool/users/ovtin/sim000004.dat -r -R22996
```
