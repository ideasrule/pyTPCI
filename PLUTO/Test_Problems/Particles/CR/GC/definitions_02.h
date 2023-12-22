#define  PHYSICS                        ResRMHD
#define  DIMENSIONS                     2
#define  GEOMETRY                       CARTESIAN
#define  BODY_FORCE                     NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  RK2
#define  NTRACER                        0
#define  PARTICLES                      PARTICLES_CR
#define  USER_DEF_PARAMETERS            0

/* -- physics dependent declarations -- */

#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  DIVB_CONTROL                   CONSTRAINED_TRANSPORT
#define  DIVE_CONTROL                   NO

/* -- user-defined parameters (labels) -- */


/* [Beg] user-defined constants (do not change this line) */

#define  PARTICLES_CR_C                 1.
#define  PARTICLES_CR_E_MC              1.0
#define  PARTICLES_CR_FEEDBACK          NO
#define  PARTICLES_CR_GC                YES
#define  PARTICLES_CR_4VEL              YES
#define  SETUP                          2
#define  UNIT_LENGTH                    100.

/* [End] user-defined constants (do not change this line) */
