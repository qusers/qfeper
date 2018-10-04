# qfeper  

### Introduction:  
qfeper is a program for automated writing of so-called fep files used by the **Q** package. 
For each state, one topology file is needed and all information is read from topology files.
The program needs an input file containing the number of states, the name of
each topology file, the number of q atoms, and a one to one mapping (see tests).
The command "./qfeper h" will print out the format of the instruction file.

### How to use:  

    ./qfeper   h                    For instruction file format 
    ./qfeper  "input"               For creating FEP file 
    ./qfeper  "input" p             For printing details information on std output 
    ./qfeper  "input" s             For splitting the FEP file to 2 states FEP files 
    ./qfeper  "input" sp            It is also accepted to add both functions (sp or ps) 

### How it works:  

The reference state is the first topology file. All atom numbers are translated
to the first state numbering. Only the bond, angle, torsion, and improper will be
considered that their atoms are q atoms. This implies that the user has to choose
the q atoms so that their range covers all possible perturbations. Subsequently,
the parameters will be compared in different states and only the ones that
are changing will be printed out to the fep file. The procedure is different for
atom types and charges. All of the q atoms types and charges will be added to the fep
file (it makes it easier to further manipulate the fep/EVB if needed and has
no effect on final results). The coupling will be suggested based on the presence
of a breaking/forming bond in an angle, torsion and/or improper. The program can
handle more than two states; however, the format written out might become messy
in case of more than six states.


### How to compile:   
qfeper was compiled and tested using Intel Fortran's ifort and GNU's GCC gfortran in linux(ubuntu). The
compilation commands follow:

ifort:  

    ifort qfeper_pars.f90 qfeper_analyz.f90 qfeper.f90 -o qfeper

gfortran:  

    gfortran qfeper_pars.f90 qfeper_analyz.f90 qfeper.f90 -o qfeper

debugging:   

    ifort -check all -debug all qfeper_pars.f90 qfeper_analyz.f90 qfeper.f90 -o qfeper


As an alternative one can use the make file included in the src folder, for usage just type:

    make

In your terminal.


### No Support Note:

Please **NOTE** that **qfeper** at the moment cannot give any support
but, if you're really interested and have the means to support
**qfeper** please contact its author.


Masoud Kazemi   
<kazemimsoud@gmail.com>
