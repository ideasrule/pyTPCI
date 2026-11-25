"""Physical constants in cgs, shared by all pyTPCI Python code.

The values are those of PLUTO (CONST_* in PLUTO/Src/pluto.h), so that Python and PLUTO agree exactly.
"""

# --- same values as pluto.h ---
AU = 1.49597892e13        # CONST_au      astronomical unit
G = 6.6726e-8             # CONST_G       gravitational constant
k_B = 1.3806505e-16       # CONST_kB      Boltzmann constant
m_H = 1.6733e-24          # CONST_mH      hydrogen atom mass
M_sun = 2.e33             # CONST_Msun    solar mass
M_earth = 5.9736e27       # CONST_Mearth  Earth mass
R_earth = 6.378136e8      # CONST_Rearth  Earth radius

# --- not defined in pluto.h ---
R_jup = 11.2 * R_earth    # Jupiter radius
M_jup = 1.898e30          # Jupiter mass
nu_Ryd = 3.2898e15        # frequency of 1 Rydberg [Hz]
T_Lya = 118348.0          # Lyman-alpha excitation energy / k_B [K]
