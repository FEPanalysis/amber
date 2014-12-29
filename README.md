There are five subdirectories in `/input_data` named `ti00[1-5]`, each containing `ti00[1-5].out` (which are copies of the `ti001.out`, with `clambda` edited)
and `ti00[1-5].en` (copies of `ti001.en`).

The script is called as follows:

`python alchemical_analysis.py -a AMBER -d input_data -p ti*/ti -q out -u kcal -r 8 -v`

The argument that follows the `-a` flag should not necessarily be in all capitals; any combination of lower- and upper-case letters is OK.
`input_data` is the path to the directory with the data files.

The `-p` flag seeks for the prefix of the data file. If the data files are in multiple subdirectories,
the name of those subdirectories (in a form of the glob pattern) should preface the file prefix (like `ti*/ti` above).

The `-v` flag (verbose option) prints out the averages and RMS fluctuations computed for each quantity in the MDEN file.
