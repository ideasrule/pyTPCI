#define  PLANET_RADIUS      (2.9 * CONST_Rearth)
#define  PLANET_MASS        (11 * CONST_Mearth)
#define  MU                 1.23
#define  H_FRAC             0.92
#define  INIT_TEMP          740
#define  STAR_MASS          (0.75 * CONST_Msun)
#define  SEMIMAJOR_AXIS     (0.0596 * CONST_au)
#define  n_SURF             1e14
#define  P_SURF             (n_SURF * CONST_kB * INIT_TEMP)
#define  RHO_SURF           (P_SURF * (MU * CONST_mH) / (CONST_kB * INIT_TEMP))
#define  BETA               (CONST_G * PLANET_MASS * (MU * CONST_mH) / (PLANET_RADIUS * CONST_kB * INIT_TEMP))

