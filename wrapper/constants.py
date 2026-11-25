"""Physical constants in cgs, shared by all pyTPCI Python code.

The values are those of PLUTO (CONST_* in PLUTO/Src/pluto.h), so that Python and PLUTO agree exactly.
"""

# --- same values as pluto.h ---
amu = 1.66053886e-24      # CONST_amu     atomic mass unit
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

# Standard atomic weights [amu] of the 30 elements Cloudy follows, keyed by the names in its abundance files
atomic_weight = {
    "HYDROGEN": 1.00794, "HELIUM": 4.002602, "LITHIUM": 6.941, "BERYLLIUM": 9.012182, "BORON": 10.811,
    "CARBON": 12.0107, "NITROGEN": 14.0067, "OXYGEN": 15.9994, "FLUORINE": 18.9984032, "NEON": 20.1797,
    "SODIUM": 22.98976928, "MAGNESIUM": 24.305, "ALUMINIUM": 26.9815386, "SILICON": 28.0855,
    "PHOSPHORUS": 30.973762, "SULPHUR": 32.065, "CHLORINE": 35.453, "ARGON": 39.948, "POTASSIUM": 39.0983,
    "CALCIUM": 40.078, "SCANDIUM": 44.955912, "TITANIUM": 47.867, "VANADIUM": 50.9415, "CHROMIUM": 51.9961,
    "MANGANESE": 54.938045, "IRON": 55.845, "COBALT": 58.933195, "NICKEL": 58.6934, "COPPER": 63.546,
    "ZINC": 65.38,
}
