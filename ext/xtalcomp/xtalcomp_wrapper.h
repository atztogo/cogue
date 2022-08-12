#ifndef __XTALCOMP_WRAPPER_H__
#define __XTALCOMP_WRAPPER_H__

#ifdef __cplusplus
extern "C" {
#endif
  int xtalcomp(const int num_atom,
	       double lattice1[3][3],
	       const int types1[],
	       double positions1[][3],
	       double lattice2[3][3],
	       const int types2[],
	       double positions2[][3],
	       const double tolerance,
	       const double angle_tolerance);
#ifdef __cplusplus
}
#endif

#endif
