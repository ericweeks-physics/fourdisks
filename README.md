# fourdisks
Simulation code from Weeks and Criddle, Phys Rev E 102, 062153 (2020):  four disks moving in a small 2D arena

The simulation considers four disks, radius = 1, in a circular arena, radius 3 + epsilon.  There are two main programs, one that does Brownian Dynamics and the other that does Molecular Dynamics (ballistic motion with collisions).  Everything is in 2D.


BROWNIAN MOTION CODE
====================
This program initializes the four disks.  It then randomly selects a disk, attempts to move it a fixed distance in a random direction, and accepts that move if it does not overlap any other disk or the edge of the arena.  The program can save the trajectories.  It also calculates all of the free energy landscapes described in our paper, so there's a lot of code associated with those calculations.

evolvefour07.pro -- main IDL program implementing the Brownian dynamics code.  Calling procedure:

a = evolvefour07(epsilon,circs,steps,DIFFN=DIFFN,SAVESTEPS=SAVESTEPS,TEMP=TEMP,NOINIT=NOINIT,SEED=SEED,FILENAME=FILENAME,MORESTEPS=MORESTEPS,TRAJSAVE=TRAJSAVE,NOPLOT=NOPLOT,NOPIC=NOPIC)

* epsilon:  size of arena.  Required.  Typical values are 0.1 to 1.0.
* circs:  could be an empty variable.  If set, this should be a 2x4 array with initial positions of circles (disks)
* steps:  how many steps to simulate.  In general the program wants this to be a long integer, and should be at least 200000L and probably better to be 1000000L (one million) as some free energy calculations are only done intermittently.
* DIFFN:  size of diffusion constant.  1e-4 by default.  I often used 1e-5.
* SAVESTEPS:  how often to save the trajectory.  I usually used 1000L, but for the smaller epsilons for which disk rearrangements occurred infrequently, I would make this 10000L or 20000L.
* TEMP:  undocumented/untested feature, the program can use this as a temperature.  This treats particles as soft disks that have an (overlap^6) potential energy.  Mostly this is leftover from Du and Weeks, PRE 2016 code that dealt with soft particles.  I have not tested this with four disks.  I usually set temperature = -1.0 for safety's sake, to make sure the soft disk code is completely ignored.
* NOINIT:  Use /NOINIT if you want to use the user-defined initial positions in the 'circ' variable.  Otherwise, even if you provide initial positions in 'circ', they will be overruled by the program.
* SEED:  used for the random number generator.
* FILENAME:  set this as a suffix, if you want to save the data in a variety of files (see below)
* MORESTEPS:  Usually I told the program to do 1000000000L steps (1e9) but if I wanted a larger number, I used moresteps.  MORESTEPS=2 would do twice as many.  For the smallest values of epsilon, I might do steps=1e9 and moresteps=20L, which might take several weeks to complete.
* TRAJSAVE:  Use /trajsave to write out the trajectory data.  Otherwise, just saves the free energy data.
* NOPLOT:  /noplot to not plot various quantities as the program runs
* NOPIC:  /nopic to not display a picture related to a free energy landscape as the program runs

Files that this program saves:
* ab.filename:  
* abu.filename:
* uv.filename:
* fab.filename:
* fabu.filename:
* tau.filename:  the intervals between rearrangements
* energy.filename:  
* energyu.filename:  
* energyq.filename:  
* abdiff.filename:
* abudiff.filename:
* theta.filename:  A phase space we didn't discuss in the paper.  HELP
* sum.filename:  summary file, contains (epsilon, temperature, total steps simulated, savesteps, diffn, mean lifetime, dsh0, deh,dsh1,dehu, version)
    - dsh0 = barrier height in theta phase space
    - dsh1 = barrier height in theta phase space based on unit vectors
    - deh = potential energy barrier for soft disks, may or may not be working
    - dehu = potential energy barrier for soft disks based on unit vector phase space, may or may not be working
    - version = 7.1
