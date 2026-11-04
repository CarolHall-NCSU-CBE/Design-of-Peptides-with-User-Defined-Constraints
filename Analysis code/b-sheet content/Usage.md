#### 1. Beta-Sheet Content Analysis

##### a. Compile the code

On an HPC system, run the following commands:

```bash
module load PrgEnv-intel
ifort -o bsheet_content bsheet_v3.f90
```

> **Note:** Check whether the Intel Fortran compiler is available on your HPC system. The module name may differ between systems.

This compiles `bsheet_v3.f90` into the executable `bsheet_content`.

##### b. Run the analysis

Navigate to a DMD/PRIME20 simulation folder and run the executable:

```bash
bsheet_path="/path/to/location/of/bsheet_content"
cd /path/to/DMDPRIME20_run
"$bsheet_path"
```

##### c. Output

The analysis generates a file named `bsheet_content.txt` in the simulation folder.

The file contains the following columns:

- `Collisions(Billion)`
- `Time(microsecond)`
- `Number of antiparallel residues`
- `Number of parallel residues`
- `Number of U-turns`
- `beta-sheet content`

Example output (without column headers):

```text
...
 200.98649    19.63727        94       297         0        0.71429
 201.05720    19.64418        87       308         0        0.74702
 201.12801    19.65108       105       338         0        0.75298
 201.19857    19.65799        79       332         0        0.73810
 201.25000    19.66302        83       306         0        0.73214
...
```