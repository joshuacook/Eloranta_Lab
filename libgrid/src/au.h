/*
 * Constants related to atomic units.
 *
 */

#ifndef __ATOMIC_UNITS__
#define __ATOMIC_UNITS__
#define GRID_AUTOANG   0.52917725           /* Bohr to Angstrom  */
#define GRID_AUTOK     3.15773213e5         /* Hartree to Kelvin */
#define GRID_AUTOCM1   (3.15773213e5/1.439) /* Hartree to wavenumbers */
#define GRID_HZTOCM1   3.335641E-11         /* Hz to wavenumbers */
#define GRID_AUTOAMU   (1.0/1822.88853006)  /* Me (mass of electron) to atomic mass unit */
#define GRID_AUTOFS   (0.02418884)          /* Atomic time unit to femtosecond */
#define GRID_AUTOBAR  (2.9421912E8)         /* Atomic pressure unit (Hartree/bohr**2) to bar */
/* TODO: Is this correct? */
//#define GRID_AUTOMPS  (1.0/4.57255954E-7)   /* Atomic unit to m/s (velocity) */
#endif

