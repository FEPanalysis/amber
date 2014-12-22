Attached please find a .tar ball with the scripts and example inputs.
There are five subdirectories in `/data` named `ti00[1-5]` each containing `ti00[1-5].out` (which are the copies of the ti001.out
you sent me, with `clambda` edited).

The script is called as follows:

`python alchemical_analysis.py -a AMBER -d input_data -p ti*/ti -q out -u kcal -r 8`


The argument that follows the `-a` flag should not necesarily be in all capitals; any combination of lower- and upper-case letters are OK.
`data` is the path to the directory with the data files

the `-p` flag seeks for the prefix of the data file. If, like in your case, the data files are in multiple subdirectories,
the name of those subdirectories (in a form of the glob pattern) should preface the file prefix (like `ti*/ti` above).
 
Luckily, my concerns about memory consumption were proved to be groundless.
The Amber data file parser is the most efficient among the three in terms of data storage (`parser_sire.py` should be rewritten that way, too; in the case of Gromacs' FEP-related data, though, it's not trivial). The dvdl timeseries, after having been put through the
autocorrelation analysis, is dumped and only two float numbers, the average and standard error of the mean, are stored for each lambda. The size of the numpy entry being 8 bytes,
any of the timeseries should be of not less than 100M entries to start causing problems.

Let me know what you think of it. (I do remember about the gradients for various force field terms).

