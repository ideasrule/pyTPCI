## Notes

These test various abort modes.   The tests are:

* abort.in - throw a cloudy_abort exception
* assert.in - the assert macro

* avx_domain.in - invalid arguments passed to a vectorized math function

* bounds_array.in - test array bounds
* bounds_auto.in
* bounds_avx_ptr.in
* bounds_multi.in
* bounds_multi_iter.in
* bounds_static.in
* bounds_vector.in

* exception.in - throws a C++ exception

* fpe_isnan_float.in - floating point error
* fpe_isnan.in
* fpe_longoverflow.in
* fpe_NaN.in
* fpe_overflow.in
* fpe_setnan_float.in
* fpe_setnan.in
* fpe_zero.in

* grid_corners.in - corners of grid exceed temperature limits of code, test that
  control passes smoothly back to main with all punch produced after abort is
  declared

* insanity.in - tests TotalInsanity

* limit_lte_he1_nomole_ste_nocoll2.in - a fairly large (but quick) grid where
  random models fail due to various errors. this tests if the grid is resilient
  to a variety of errors; also check if SAVE GRID output remains complete
  despite these errors

* monitor.in - test various monitor commands

* segfault.in - test segmentation fault

* undef.in - variable is undefined

Scripts:

* run_parallel.pl - perl script to run all tests
* clean_tsuite.pl - clean up after run
* run_sequence.pl - helper script
* checkall.pl - empty stub
