Command line options:

- `-t` : derived type to compare (tracers, dynamics, ice)
- `-d` : two directories with results to compare
- `-r` : tolerance (relative difference). The default is `1.0E-15`.

At the moment one file is generated from each MPI process and the number of MPI processes has to be the same of the generated results. The output is trying to mimic the one from `cdo diffn`.

Example:

`mpirun -n 8 bin/fdiff -t tracers -d /relative/path/directory1 relative/path/directory2 -r 1.0E-14`

