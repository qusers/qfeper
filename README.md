# qfeper #

### Introduction: ###
qfeper is a program to write the fep files that is used by Q package. For each
state, one topology file is needed and all information is read from topology files.
The program need an input file containing the number of states, the name of
each topology file, the number of q atoms, and an one to one map (see tests).
The command ”./qfeper h” will print out the format of instruction file.

### The program usage: ###
>    ./qfeper   h                    For instruction file format 
>    ./qfeper  "input"               For creating fep file 
>    ./qfeper  "input" p             For printing details information on std output 
>    ./qfeper  "input" s             For spliting the fep file to 2 states fep files 
>    ./qfeper  "input" sp            It is also accepted to add both functions (sp or ps) 