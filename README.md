pyTPCI is an improved version of The PLUTO-CLOUDY Interface, coupling together hydrodynamics code [PLUTO](https://plutocode.ph.unito.it/) and gas microphysics code [CLOUDY](https://gitlab.nublado.org/cloudy/cloudy).

## Setup
1) Make sure you have libhdf5-dev installed.  On Ubuntu, you can do "sudo apt install libhdf5-dev"
2) After cloning the pyTPCI directory to your computer, copy the `wrapper` folder for **each instance of pyTPCI you wish to run**. 
3) Define the environment variable `PLUTO_DIR` and run `$PLUTO_DIR/setup.py` in your folder. Further customization of PLUTO takes place in `pluto_template.ini`.
4) Run `make` and `make clean` in your folder to produce the PLUTO executable.
5) Run `make` in `cloudy/source` to produce `cloudy.exe`.

pyTPCI runs PLUTO and then CLOUDY in alternating steps, transferring heating information between them with `write_heating_file`. In order to save computational power and improve speed, 
by default pyTPCI only calls CLOUDY if there is at least a 10% maximum fractional difference (`max_rel_diff`) in density or pressure between PLUTO files. 

Three refinements of this coupling:
- **Refresh.** CLOUDY is also called whenever the heating table is older than `MAX_AGE` (0.5 code time units). Otherwise a quiet flow never triggers a call and PLUTO settles against a stale table. The 10% trigger only looks at r < 5 Rp.
- **Lyman-alpha thermostat.** `curr_heat.dat` holds heating and cooling per gram separately, the density and P/rho of the state CLOUDY was given, and local exponents. `radiat.c` rescales the rates to each cell's current state, so the steep temperature dependence of Lyman-alpha cooling acts inside PLUTO between CLOUDY calls. The only non-zero exponent is dln(cooling)/dlnT = f_Lya * 118348 K / T, where f_Lya is the fraction of the cooling CLOUDY attributes to collisionally excited hydrogen lines.
- **Filtered CLOUDY input.** The density, pressure and velocity passed to CLOUDY (and used by the trigger) are Gaussian-smoothed in radius (width 0.02 r, only above 1.3-1.4 Rp) to remove small acoustic pulses, which would otherwise be imprinted on the heating table and re-excited. PLUTO's own state is not modified.

## Running pyTPCI
1) Create a stellar spectrum at the planet's surface (semimajor axis), with columns of log10(frequency (Hz)) and log10(F_nu), spectral irradiance (ergs/(s cm^2 Hz)). At the end this also contains the band-integrated flux at the planet's surface, and the energy range in Rydbergs. See `spectra_example.ini`.
2) Edit `tpci.py` to set the system parameters: the stellar spectrum file, the semimajor axis, planet radius and mass, stellar mass, initial temperature, and metallicity. Upon running, this will generate a CLOUDY `params.h` file and an input script ending in `.in`.
3) Open `screen` and run `python tpci.py` to start at t=0 and use the default timestep dt=0.01. If restarting from a previous pyTPCI run, instead do `python tpci.py global_ind time timestep`, where `global_ind` is the number of the CLOUDY file you wish to start from.
4) Monitor pyTPCI's progress by plotting PLUTO and CLOUDY files using `plot.py` and `cloudy_plot.py`.

## Initial condition
`IC_MODE` in `tpci.py` selects the state at t=0:
- `"parker"` (default): `parker_ic.py` builds a transonic wind on top of the cold hydrostatic base and writes it to `ic.dat`, which `init.c` reads. The first guess is an isothermal 10^4 K wind with an energy-limited mass-loss rate (10% efficiency, XUV flux taken from the stellar spectrum file). The temperature, mean molecular weight and mass-loss rate are then iterated to CLOUDY's thermal equilibrium for the wind's own density profile. This takes a few tens of CLOUDY runs before the first PLUTO step, but the simulation then starts close to the steady state instead of having to lift the atmosphere off the planet.
- `"isothermal"`: the hydrostatic atmosphere at `init_temp` (no `ic.dat`).

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
