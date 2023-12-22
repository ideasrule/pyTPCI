## Notes

Trapezium cluster massive stars
stellar types listed in Table 1 of Ferland+12
 http://adsabs.harvard.edu/abs/2012ApJ...757...79F

That paper used TLusty stellar atmospheres.
Trapezium_Tlusty shows this

Wagle+ in preparation claim that WMbasic gives better fit to [Ne III] lines.
Trapezium_WMbasic shows this

Veusz file and pdf compare these two.

* * *

The original form came from Will Henney in email of
2012 Dec 5.  His notes follow:

>    Stellar parameters for C, A, D are from Simón-Díaz \cite{2006A&A…448..351S}
>    I originally had th1B as a B2V star, which should be cooler than what Gary had (28kK and 10**4.4 Lsun), but I can't find the provenance of this
>        actually th1B (BM Ori) is a quintuple system: see SIMBAD record
>        SIMBAD says primary is B1V citing \cite{Close:2003}, but that paper does not mention spectral types at all. I suspect this is a confusion with B1, B2, etc being the designations of the quintuple components
>        \cite{Weigelt:1999} find a B3 spectral type for the primary (B1) with log L = 3.25 and T = 18000
>        B2, B3, and B4 are all lower mass, but B5 has mass similar to B1 accoding to \cite{Weigelt:1999}. However, this must be a mistake: according to \cite{Abt:1991}, who detect B1 and B5 as an eclipsing spectroscopic binary, the mass ratio is 0.3, implying that B5 has a mass of only 2 Msun..
>    th1C is a triple system \cite{Lehmann:2010a}
>        \cite{Schertl:2003} have broad constraints on the parameters of C2 (dark gray line in their Fig 4). If we take the middle of this range, we get a giantish B3IV star: log L = 3.8, T = 17000, M = 9 Msun, age much less than 1 Myr.
>        We can work out the surface gravity:
>            g = G M / R2 = G M 4 pi sigma T4 / L = 6.673e-8 9 1.989e33 4 pi 5.6703e-5 17000**4 / 10**3.8 3.82e33 = 2950 cm/s2 => log g = 3.5
>        [X] We should check these parameters against \cite{Lehmann:2010a}, who say 12 Msun in their abstract. The article itself is not freely available -grab it when I get back to work
>            If we use the same evo tracks as in \cite{Schertl:2003}, then this would imply B2IV with log L = 4.0, T = 21000
>            So, after checking out the Lehman article, it seems the 12 Msun is robust, since it is based on orbital solution. The translation to luminosity and Teff is rather more speculative.
>            In fact, according to \cite{Kraus:2007}, the companion is O9.5V with Teff = 32000, log L = 4.68. This is also consistent with the Schertl tracks, but it is assuming 15.5 Msun.
>            If we instead use 12 Msun as the best mass estimate \cite{Lehmann:2010a}, then we would get log L = 4.2, Teff = 25,000
>                Recalculate surface gravity for new parameters
>                    g = G M / R2 = G M 4 pi sigma T4 / L = 6.673e-8 12 1.989e33 4 pi 5.6703e-5 25000**4 / 10**4.2 3.82e33 = 7322 cm/s2 => log g = 3.86
> 
> 2014 May 06, changed some tlusty stars to WMBasic - had used tlusty for all stars originally.  WMBasic produces better agreement with H II region optical emission, especially [Ne III]
>             Using this new solution for th1C2 changes the broadband ratios by about 5%.
>     th1A is a double system
> 
> stars = dict(
>     th1C = dict(T=39000., g=4.1, L=5.31),
>     th1A = dict(T=30000., g=4.0, L=4.45),
>     th1D = dict(T=32000., g=4.2, L=4.47),
>     th1B = dict(T=18000., g=4.1, L=3.25),
>     th1C2 = dict(T=25000., g=3.86, L=4.2),
>     )
