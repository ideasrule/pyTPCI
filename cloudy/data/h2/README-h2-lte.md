## Notes

The small H2 model computes the H2 line cooling according to the prescription
of Glover & Abel 2008, MNRAS, 388, 1627.  The total cooling at intermediate
densities is computed according to their eqn (39), which combines the cooling
at low densities and at LTE.  The low-density fits for cooling due to collisions
between para-H2 and ortho-H2 with various colliders are implemented in the
source code.  The LTE cooling is evaluated by interpolation over the tabulated
cooling computed with their eqn (38).  To produce the tabulated data, run the
following command in the present directory:

$ ../../source/cloudy.exe -r get-H2-LTE-cooling

The LTE H2 cooling data will be contained in the file lte_cooling.dat.
