pyTPCI is an improved version of The PLUTO-CLOUDY Interface, coupling together hydrodynamics code [PLUTO](https://plutocode.ph.unito.it/) and gas microphysics code [CLOUDY](https://gitlab.nublado.org/cloudy/cloudy).

## Setup
1) Make sure you have libhdf5-dev installed.  On Ubuntu, you can do "sudo apt install libhdf5-dev"
2) After cloning the pyTPCI directory to your computer, copy the `wrapper` folder for **each instance of pyTPCI you wish to run**. 
3) Define the environment variable `PLUTO_DIR` and run `$PLUTO_DIR/setup.py` in your folder. Further customization of PLUTO takes place in `pluto_template.ini`.
4) Run `make` and `make clean` in your folder to produce the PLUTO executable.
5) Run `make` in `cloudy/source` to produce `cloudy.exe`.

pyTPCI runs PLUTO and then CLOUDY in alternating steps, transferring heating information between them with `write_heating_file`. In order to save computational power and improve speed, 
by default pyTPCI only calls CLOUDY if there is at least a 10% maximum fractional difference (`max_rel_diff`) in density or pressure between PLUTO files. 

## Running pyTPCI
1) Create a stellar spectrum at the planet's surface (semimajor axis), with columns of log10(frequency (Hz)) and log10(F_nu), spectral irradiance (ergs/(s cm^2 Hz)). At the end this also contains the band-integrated flux at the planet's surface, and the energy range in Rydbergs. See `spectra_example.ini`.
2) Edit `tpci.py` to set the system parameters: the stellar spectrum file, the semimajor axis, planet radius and mass, stellar mass, initial temperature, and metallicity. Upon running, this will generate a CLOUDY `params.h` file and an input script ending in `.in`.
3) Open `screen` and run `python tpci.py` to start at t=0 and use the default timestep dt=0.01. If restarting from a previous pyTPCI run, instead do `python tpci.py global_ind time timestep`, where `global_ind` is the number of the CLOUDY file you wish to start from.
4) Monitor pyTPCI's progress by plotting PLUTO and CLOUDY files using `plot.py` and `cloudy_plot.py`.

## Plotting information
- To plot one file, run `python plot.py global_ind` where `global_ind` is the number of the file, such as "54".
- To plot multiple files, run `python plot.py start stop step`, where `start` and `stop` are the initial and final file numbers to plot, with a given `step`. Plotting files will automatically skip file names that do not exist. 

## Logging and Docs information
- pyTPCI run information is in `tpci_log.txt`.
- PLUTO output logs for each PLUTO run are appended in `pluto_log.txt`
- CLOUDY docs are in `cloudy/docs` and PLUTO docs are in `PLUTO/Doc`

## Troubleshooting and Common Errors
- Is the path to the `spectra.ini` file correct?
- Does `definitions.h` contain `#include params.h` at the top?
- If CLOUDY crashes quickly, try using a smaller dt=0.001
