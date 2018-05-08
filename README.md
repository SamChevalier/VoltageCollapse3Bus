# VoltageCollapse3Bus
Simulation and analysis of voltage collapse in a 3 bus power system. The simulation is run via Power System Analysis Toolbox (version 2.1.10 at http://faraday1.ucd.ie/psat.html). In order to stop the time domain simulation appropraitely, the PSAT file fm_int.m has been hacked slightly. The hacked version is given in the repository, and it must replace the original version in PSAT's main folder.

There are 6 files:

1. Main_TEST.m - This is the main test file which runs the voltage collapse simualtion. It loops over all three controllers.
2. Calc_CritVars.m - This file calculates the critical variance at the load bus using a statistical solver.
3. perturb.m - This is the PSAT perturbation file which adds load noise (fast and slow) and implements control
4. **dyn_4.m** - This is the PSAT data file which defines the system: line, generator, load, control, and shunt data. The controllable shunt is used to replicate the SVC.
5. Noie_Data.mat - This file contains the fast and slow noise data which is applied in each simulation. It is saved in memory such that each controller sees identical simulation dynamics.
6. fm_int.m - This is a PSAT data file which should *replace* the fm_int.m file which is in the psat folder.
