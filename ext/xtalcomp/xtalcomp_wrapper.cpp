#include "xtalcomp.h"
#include "xtalcomp_wrapper.h"

#include <algorithm>

extern "C" {

  int xtalcomp(const int num_atom,
	       double lattice1[3][3],
	       const int types1[],
	       double positions1[][3],
	       double lattice2[3][3],
	       const int types2[],
	       double positions2[][3],
	       const double tolerance,
	       const double angle_tolerance)
  {
    int i;
    bool match;

    XcMatrix xclattice1 (lattice1);
    XcMatrix xclattice2 (lattice2);
    std::vector<XcVector> xcpositions1, xcpositions2;
    std::vector<unsigned int> xctypes1, xctypes2;

    xcpositions1.reserve(num_atom);
    xcpositions2.reserve(num_atom);
    xctypes1.reserve(num_atom);
    xctypes2.reserve(num_atom);
  

    for (i = 0; i < num_atom; i++) {
      xcpositions1.push_back(XcVector(positions1[i][0],
				      positions1[i][1],
				      positions1[i][2]));
      xcpositions2.push_back(XcVector(positions2[i][0],
				      positions2[i][1],
				      positions2[i][2]));
      xctypes1.push_back(types1[i]);
      xctypes2.push_back(types2[i]);
    }

    match = XtalComp::compare(xclattice1, xctypes1, xcpositions1,
			      xclattice2, xctypes2, xcpositions2,
			      NULL, tolerance, angle_tolerance);

    if (match) {
      return 1;
    } else {
      return 0;
    }
  }

}
