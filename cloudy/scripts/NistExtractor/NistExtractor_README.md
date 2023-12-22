NistExtractor
=============
NistExtractor.py retrieves the energy levels and transition probabilities from NIST for a single species. 
The desired species is the first parameter.
The second and optional parameter allows the user to specify an early cutoff to limit the number of levels to be used by Cloudy



Example: Fe IX  with 30 levels

```
NistExtractor.py Fe_IX 30
```

NEWrapper.py allows for NistExtractor to be run on multiple species.
The required parameter is a file listing the desired species and optionally the level limit.

See "all_species.txt" and "cloudy_species.txt" as examples of how to format the input file.


