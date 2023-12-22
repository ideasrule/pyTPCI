# Problem with Duplicate chemical reactions

## Gargi Shaw,  September 2020


A calculation using the development version of the code will
generate several warnings
about duplicate reactions.  These are silenced on release versions.
Gargi Shaw spent some time looking into this and created the following report.

The following points discuss why we are getting duplicate reaction
rates and how to rectify it.
The files involved are `mole.cpp`, `mole_reactions.cpp`, `mole_co_base.dat`,
`mole_co_federman.dat` and `mole_co_umisthack.dat`.

The reactions listed in  `mole_co_federman.dat` and `mole_co_umisthack.dat`
are also listed in `mole_co_base.dat`. 
Hence we are getting duplicate reaction comment.
The relevant blocks of code are the following:

## 1. mole.cpp (line no. 24 to 28)

````
lgLeidenHack = false;
/* option to use diffuse cloud chemical rates from Table 8 of Federman, S. R. & Zsargo, J. 
2003, ApJ, 589, 319. By default, this is false - changed with set chemistry command */
lgFederman = true; 
````


## 2. mole\_reactions.cpp (line 1830 to 1840)

````
source = base;
read_data("mole_co_base.dat",parse_base);
if (mole_global.lgFederman)
{
        source = federman;
        read_data("mole_co_federman.dat",parse_base);
}
if (!mole_global.lgLeidenHack) 
{
        source = umisthack;
        read_data("mole_co_umisthack.dat",parse_base);
}
````


From points 1. and 2. we can see that reactions listed in
`mole_co_federman.dat` and `mole_co_umisthack.dat` override
the rates listed in `mole_co_base.dat` for those reactions.
Hence, rates from  `mole_co_federman.dat` and `mole_co_umisthack.dat`
are finally used for calculations.

## 3. The rates from `mole_co_federman.dat` and `mole_co_umisthack.dat`
are different from those listed in `mole_co_base.dat`.
    
Now the question is, which rates do we want to set as default values?
Was it intended to use rates listed in `mole_co_federman.dat` and
`mole_co_umisthack.dat` 
as default values? Do we want to keep all default rates in one file,
namely `mole_co_base.dat` for ease of future updates?

## 4. If the answers are yes,

Either, update the rates in `mole_co_base.dat` and make the following
changes in 

````
lgFederman = false (mole.cpp);   
````

and

````
if (mole_global.lgLeidenHack)   ( mole_reactions.cpp)        
````

OR
                                                   
just comment out those reactions from `mole_co_base.dat`.   
      

## 5. If we don’t want to use rates listed in  `mole_co_federman.dat`
and `mole_co_umisthack.dat`
as default values, then make the following changes,

````
   lgFederman = false (mole.cpp);   
````

and
    
````
     if (mole_global.lgLeidenHack)   ( mole_reactions.cpp)       
````


Please note, this will result in huge changes since the rate values
listed in `mole_co_base.dat` are different.

