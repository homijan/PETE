# PETE
Plasma Euler and Transport Equations - new generation MHD code.

The hydrodynamic solver of PETE inherits from [Laghos](https://github.com/CEED/Laghos), 
which solves the time-dependent Euler equations of compressible gas dynamics 
in a moving Lagrangian frame using unstructured high-order finite element 
spatial discretization and explicit high-order time-stepping.

The Laghos miniapp is part of the [CEED software suite](http://ceed.exascaleproject.org/software) and is based on the discretization method described in 
the following article:

> V. Dobrev, Tz. Kolev and R. Rieben <br>
> [High-order curvilinear finite element methods for Lagrangian hydrodynamics](https://doi.org/10.1137/120864672) <br>
> *SIAM Journal on Scientific Computing*, (34) 2012, pp. B606–B641.

## Building

PETE has the following external dependencies:

- *hypre*, used for parallel linear algebra, we recommend version 2.10.0b<br>
   https://computation.llnl.gov/casc/hypre/software.html

-  METIS, used for parallel domain decomposition (optional), we recommend [version 4.0.3](http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/OLD/metis-4.0.3.tar.gz) <br>
   http://glaros.dtc.umn.edu/gkhome/metis/metis/download

- MFEM, used for (high-order) finite element discretization, its GitHub master branch <br>
  https://github.com/mfem/mfem

- Laghos (LAGrangian High-Order Solver), used as the Euler's equation solver, its GitHub master branch <br>
  https://github.com/CEED/Laghos

- HerEOS, providing a set of plasma physics equation of states, its BitBucket master branch <br>
  https://bitbucket.org/Zemaster/hereos

To build PETE, download/clone necessary libraries.
Clone PETE:
```sh
~> git clone https://github.com/homijan/PETE.git
```
Clone Laghos:
```sh
~> git clone https://github.com/CEED/Laghos.git
```
Clone MFEM:
```sh
~> git clone git@github.com:mfem/mfem.git ./mfem
```
Clone HEREOS:
```sh
~> git clone https://bitbucket.org/Zemaster/hereos.git
```
Download *hypre* and METIS from the links above
and put everything on the same level as the `PETE` directory:
```sh
~> ls
PETE/ Laghos/  mfem/ hereos/ hypre-2.10.0b.tar.gz  metis-4.0.tar.gz
```

Build *hypre*:
```sh
~> tar -zxvf hypre-2.10.0b.tar.gz
~> cd hypre-2.10.0b/src/
~/hypre-2.10.0b/src> ./configure --disable-fortran
~/hypre-2.10.0b/src> make -j
~/hypre-2.10.0b/src> cd ../..
```
For large runs (problem size above 2 billion unknowns), add the
`--enable-bigint` option to the above `configure` line.

Build METIS:
```sh
~> tar -zxvf metis-4.0.3.tar.gz
~> cd metis-4.0.3
~/metis-4.0.3> make
~/metis-4.0.3> cd ..
~> ln -s metis-4.0.3 metis-4.0
```
This build is optional, as MFEM can be build without METIS by specifying
`MFEM_USE_METIS = NO` below.

Build the parallel version of MFEM:
```sh
~> cd mfem/
~/mfem> make parallel -j
~/mfem> cd ..
```

Build HerEOS
```sh
~> cd hereos/buid
~/build> cmake ..
~/build> make
```

Build PETE
```sh
~> cd PETE/
~/PETE> make
```
This can be followed by `make test` and `make install` to check and install the
build respectively. See `make help` for additional options.

## Running

#### 1D

###### IG

Sod shock tube
```sh
mpirun -np 8 pete -p 2 -m data/segment01.mesh -rs 4 -tf 0.6 -vis -fa -ot 3 -ok 4
```

Sedov blast
```sh
mpirun -np 8 pete -p 1 -m data/segment01.mesh -rs 4 -tf 0.8 -vis -fa -ot 3 -ok 4
```

Warm wall Sedov blast
```sh
mpirun -np 8 pete -p 4 -m data/segment01.mesh -rs 4 -tf 0.8 -vis -fa -ot 3 -ok 4
```


###### HerEOS

Hydrogen and Aluminum - Sod shock tube
```sh
mpirun -np 8 pete -p 2 -m data/segment01.mesh -rs 4 -tf 1e-7 -vis -fa -heos
mpirun -np 8 pete -p 2 -m data/segment01.mesh -rs 4 -tf 1e-6 -vis -fa -heos -mn Al
```

Hydrogen - warm wall Sedov blast
```sh
mpirun -np 8 pete -p 4 -m data/segment01.mesh -rs 4 -tf 1e-6 -vis -fa -ot 3 -ok 4 -heos -mn H
```


#### 2D

###### IG

Sedov blast
```sh
mpirun -np 8 pete -p 1 -m data/square01_quad.mesh -rs 3 -tf 0.8 -vis -pa -ot 2 -ok 3
```

Warm wall Sedov blast
```sh
mpirun -np 8 pete -p 4 -m data/square01_quad.mesh -rs 3 -tf 0.8 -vis -pa -ot 1 -ok 2
```

Triple-point problem
```sh
mpirun -np 8 pete -p 3 -m data/rectangle01_quad.mesh -rs 2 -tf 2.5 -cfl 0.025 -vis -pa
```


###### HerEOS

Hydrogen - warm wall Sedov blast
```sh
mpirun -np 8 pete -p 4 -m data/square01_quad.mesh -rs 3 -tf 1.5e-7 -vis -pa -ot 1 -ok 2 -heos -mn H
```

Hydrogen - triple-point problem
```sh
mpirun -np 8 pete -p 5 -m data/rectangle01_quad.mesh -rs 2 -tf 2.5 -cfl 0.025 -vis -pa -heos -mn H
```


#### 3D

###### IG

Sedov blast
```sh
mpirun -np 8 pete -p 1 -m data/cube01_hex.mesh -rs 2 -tf 0.6 -vis -pa
```

Warm wall Sedov blast
```sh
mpirun -np 8 pete -p 4 -m data/cube01_hex.mesh -rs 2 -tf 0.6 -vis -pa
```


###### HerEOS

Hydrogen - warm wall Sedov blast
```sh
mpirun -np 8 pete -p 4 -m data/cube01_hex.mesh -rs 2 -tf 1e-7 -vis -pa -heos -mn H
```

