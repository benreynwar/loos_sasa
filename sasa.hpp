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

#if !defined(LOOS_SASA_HPP)
#define LOOS_SASA_HPP

#include <loos.hpp>

namespace loos {

  /* A class to hold the name and VdW radius of an atom. */
  class name_and_radius {
  public:
	std::string name;
	double radius;
	name_and_radius(std::string name_, double radius_) :
	  name(name_), radius(radius_) {};
	// Define this constructor for SWIG vector. 
	name_and_radius() :
	  name(" "), radius(-1) {};
  };

  /*
	Calculate the Solvent Accesible Surface Area (sasa) of an
	AtomicGroup.
	`selected` is a group of the atoms for which we want to find the
	surface area.
	`everything` is a group of the atoms which we want to consider.
	(e.g. `selected` might be atoms in the hydrophobic residues and
	`everything` could be all atoms.  In this way we would get the
	hydrophobic surface area).
	`cutoff` is the probe radius.
	`M` is the number of test insertions to place around every atoms.
	The higher this number is the smaller the error will be.
  */
  double
  calculate_sasa(const AtomicGroup& selected, const AtomicGroup& everything,
				 const double cutoff, const int M);

  /*
	Also calculates the SASA but allows for explicitly specifying the
	`binning` (which atoms are assigned to which rectangular bin) and
	`nar` (a mapping from atom names to VdW radii).
	Useful if you want to change the VdW radii or are using atom names
	that do not have default VdW.
	Explicitly specifying the `binning` might be useful if you're
	doing lots of SASA calculations using the same structure and you
	only want to calculate the `binning` once.
   */
  double
  calculate_sasa_with_binning_nar(const AtomicGroup& selected,
								  const double cutoff,
								  const Binning& binning,
								  const std::vector<name_and_radius>& nar,
								  const int M);

  /*
	Randomly generate a coordinate on the surface of the unit sphere.
   */
  GCoord random_sphere_point();

}

#endif
