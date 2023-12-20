#include "params.h"
#define  PHYSICS                        HD
#define  DIMENSIONS                     1
#define  GEOMETRY                       SPHERICAL
#define  BODY_FORCE                     VECTOR
#define  COOLING                        TABULATED
#define  RECONSTRUCTION                 WENO3
#define  TIME_STEPPING                  RK3
#define  NTRACER                        0
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            0

/* -- physics dependent declarations -- */

#define  DUST_FLUID                     NO
#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  THERMAL_CONDUCTION             EXPLICIT
#define  VISCOSITY                      NO
#define  ROTATING_FRAME                 NO

/* -- user-defined parameters (labels) -- */


/* [Beg] user-defined constants (do not change this line) */

#define  UNIT_DENSITY                   1e-10
#define  UNIT_LENGTH                    PLANET_RADIUS
#define  UNIT_VELOCITY                  1.0e5

/* [End] user-defined constants (do not change this line) */
