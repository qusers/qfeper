# qfeper #

### Introduction: ###
qfeper is a program to write the fep files that is used by Q package. For each
state, one topology file is needed and all information is read from topology files.
The program need an input file containing the number of states, the name of
each topology file, the number of q atoms, and an one to one map (see tests).
The command ”./qfeper h” will print out the format of instruction file.

### The program usage: ###
    ./qfeper   h                    For instruction file format 
    ./qfeper  "input"               For creating FEP file 
    ./qfeper  "input" p             For printing details information on std output 
    ./qfeper  "input" s             For splitting the FEP file to 2 states FEP files 
    ./qfeper  "input" sp            It is also accepted to add both functions (sp or ps) 

### Contribution guidelines ###
How it works:
The reference state is the first topology file. All atom numbers is translated
to first state numbering. Only the bond, angle, torsion, and improper will be
considered that their atoms are q atoms. This implies that user have to choose
the q atoms so that its range cover all possible perturbation. Subsequently,
the parameters will be compared in different states and only the ones that
are changing will be printed out to fep file. The procedure is different for
atom types and charges. All of q atoms type and charge will be added to fep
file (it make it easier to further manipulate the fep/EVB if needed and have
no effect on final results). The coupling will be suggested based on presence
of breaking/forming bond in a angle, torsion and improper. The program can
handle more than two states; however, the write out format might become messy
in case of more than three states.


###How to compile:###
qfeper was compiled and tested by ifort and gfortran in linux(ubuntu). The
compile command is as follow:

ifort:

* "ifort qfeper_pars.f90 qfeper_analyz.f90 qfeper.f90 -o qfeper"

gfortran:

* "gfortran qfeper_pars.f90 qfeper_analyz.f90 qfeper.f90 -o qfeper"

debugging: 

* "ifort -check all -debug all qfeper_pars.f90 qfeper_analyz.f90 qfeper.f90 -o qfeper"




Masoud Kazemi 
kazemimsoud@gmail.com