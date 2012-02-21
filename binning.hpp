/*
  The initial version of this file was written by Ben Reynwar.
  Ben Reynwar has assigned copyright to Alan Grossfield.
  Copyright (c) 2012, Alan Grossfield
  Department of Biochemistry and Biophysics
  School of Medicine & Dentistry, University of Rochester

  This package (LOOS) is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation under version 3 of the License.

  This package is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#if !defined(LOOS_BINNING_HPP)
#define LOOS_BINNING_HPP

#include <vector>

#include <loos_defs.hpp>
#include <AtomicGroup.hpp>

namespace loos {

  typedef int dim;
  typedef int bin;

  /*
	Bins the atoms into a rectangular grid.
	Pairs of atoms closer than a specificed cutoff distace are
	guaranteed to be in neighbouring bins.
	Useful for when you want to get all pairs of atoms within
	a certain cutoff distance efficiently.
   */
  class Binning {
  public:
	Binning(double cutoff, const AtomicGroup &group);
	~Binning();
	/* Returns an AtomicGroup containing all atoms that are potentially
	   within the cutoff distance of a given atom 'p'.
	 */
	AtomicGroup* get_nbs(pAtom p) const;
	GCoord box;
  private:
	const AtomicGroup group;
	const double cutoff;
	// The minimum and maximum coordinates (to work out system size).
	double mins[3];
	double maxs[3];
	// The number of bins in each dimensions.
	bin nbins[3];
	// The total number of bins.
	bin total_bins;
	// Whether the system is periodic in a given dimension.
	bool periodic[3];
	// An array of vectors.  Each bin has a vector containing the indices of 
	// the bins which neighbour it.
	std::vector<bin> *bin_nbs;
	// An AtomicGroup of atoms which reside in each bin.
	AtomicGroup *bin_atoms;
	// An AtomicGroup of atoms which either reside in this bin or a neighbouring bin.
	AtomicGroup *bin_nb_atoms;
	// Work out system size, decide how many bins to use...
	void populate_bin_info();
	// Work out which bins are neighbours to which bins.
	void populate_bin_nbs();
	// Return the bin which atom 'p' is in.
	bin get_bin(pAtom p) const;
	// Work out which atoms are in each bin.
	void populate_bin_atoms();
	// Work out which atoms are potentially within the cutoff distance of each bin.
	void populate_bin_nb_atoms();
  };

}

#endif
