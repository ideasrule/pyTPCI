#define  PHYSICS                        RMHD
#define  DIMENSIONS                     2
#define  GEOMETRY                       CARTESIAN
#define  BODY_FORCE                     NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  RK2
#define  NTRACER                        0
#define  PARTICLES                      PARTICLES_CR
#define  USER_DEF_PARAMETERS            1

/* -- physics dependent declarations -- */

#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  RADIATION                      NO
#define  DIVB_CONTROL                   CONSTRAINED_TRANSPORT

/* -- user-defined parameters (labels) -- */

#define   BZ_GUIDE                      0

/* [Beg] user-defined constants (do not change this line) */

#define  PARTICLES_CR_C                 1.
#define  PARTICLES_CR_FEEDBACK          NO
#define  PARTICLES_CR_E_MC              (1.0/PARTICLES_CR_C)
#define  PARTICLES_CR_GC                YES
#define  PARTICLES_CR_GC_4VEL           YES
#define  PARTICLES_CR_NCELL_MAX         1
#define  SETUP                          4

/* [End] user-defined constants (do not change this line) */
