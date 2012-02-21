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

#include <math.h>
#include <ctype.h>
#include <boost/math/constants/constants.hpp>
#include <MersenneTwister.hpp>
#include <iostream>

#include <loos.hpp>
#include <sasa.hpp>

#define DEBUG_SASA 0

namespace loos {

  const double pi = boost::math::constants::pi<double>();

  /*
	Find the VdW radius of an atom given it's name.

	`nar` maps names to radii.  The radius will be returned for the
	first name in `nar` where that name matches the start of the given
	atom (excluding digits).

	e.g. If we have the mapping:
       H -> 1.0
	   N -> 1.4
	   Hg -> 42
	Then a name of '1H' would return 1.0.
	               'Na' would return 1.4.
                   'Hg' would return 1.0. (not 42 since 'H' comes before 'Hg' vector)
                   'B' would raise an error.
    There's clearly the potential for error here so if's you've got
	unusual atom types it's necessary to double check that 'nar' is
	set up such that they're given the correct radii.
   */
  double
  get_radius_from_name(const std::string& name,
					   const std::vector<name_and_radius>& nar) {
	double radius = -1;
	for (std::vector<name_and_radius>::const_iterator it=nar.begin();
		 it!=nar.end(); it++) {
	  std::string start_name = it->name;
	  bool startswith = true;
	  if (start_name.length() > name.length())
		startswith = false;
	  else {
	    std::string::size_type j = 0;
		for (std::string::size_type i=0; i<start_name.length(); i++) {
		  // Get rid of the digits.
		  while ((isdigit(name[j])) && (j<name.length())) {
		    j++;
		  }
		  // If there is not name left then we don't have a complete match.
		  if (j >= name.length()) {
		    startswith = false;
		    break;
		  }
		  // If this char doesn't match then no match.
		  if (start_name[i] != name[j]) {
			startswith = false;
			break;
		  }
		}
	  }
	  // Found a match.
	  if (startswith) {
		radius = it->radius;
		break;
	  }
	}
	if (radius == -1) {
	  throw(std::runtime_error("Did not find a radius to match this atom name (" + name + ")."));
	}
	return radius;
  }

  /*
	Returns a vector of atom names with VdW radii that is compatible
	with those used by VMD.
   */
  std::vector<name_and_radius> get_default_names_and_radii() {
	// Use the same radii as vmd does by default.
	// See BaseMolecule::default_radius in vmd.
	// http://www.ks.uiuc.edu/Research/vmd/doxygen/BaseMolecule_8C-source.html
	std::vector<name_and_radius> nar;
	nar.push_back(name_and_radius("H", 1.0));
	nar.push_back(name_and_radius("C", 1.5));
	nar.push_back(name_and_radius("N", 1.4));
	nar.push_back(name_and_radius("O", 1.3));
	nar.push_back(name_and_radius("F", 1.2));
	nar.push_back(name_and_radius("S", 1.9));
	// P and Cl aren't from vmd but were suggested in a vmd mailing list post.
	// Not very reliable.  Should check.
	// http://www.ks.uiuc.edu/Research/vmd/mailing_list/vmd-l/4610.html
	nar.push_back(name_and_radius("P", 2.0));
	nar.push_back(name_and_radius("Cl", 2.5));
	return nar;
  };

  /* See header file for description. */
  GCoord random_sphere_point() {
	MTRand mtrand1;
	double z = mtrand1.rand()*2-1;
	double r = sqrt(1 - z*z);
	double phi = mtrand1.rand()*2*pi;
	return GCoord(r*sin(phi), r*cos(phi), z);
  }

  /* See header file for description. */
  double calculate_sasa_with_binning_nar(const AtomicGroup &selected,
										 const double radius,
										 const Binning &binning,
										 const std::vector<name_and_radius> &nar,
										 const int M) {
	// Get random points on sphere to use for all spheres.
	if (DEBUG_SASA) {
	  std::cout << "calculate_sasa_with_binning_nar: Starting" << std::endl;
	}
	GCoord box = binning.box;
	GCoord *random_points = new GCoord[M];
	// Note that we use the same random points for every atom.
	// This is because the generation of random points usually is the
	// most computationally expensive step.
	for (int i=0; i<M; i++)
	  random_points[i] = random_sphere_point();
	// Loop throuh all the atoms.
	double area = 0;
	pAtom p;
	AtomicGroup::Iterator iter(selected);
	while (p = iter()) {
	  int num_on_surface = 0;
	  GCoord pc = p->coords();
	  double pr = get_radius_from_name(p->name(), nar);
	  AtomicGroup *nbs = binning.get_nbs(p);
	  // Get nbs that are within (2*radius+pr+qr).
	  std::vector<GCoord> nb_coords;
	  std::vector<double> nb_radii;
	  AtomicGroup::Iterator iter_nbs(*nbs);
	  pAtom q;
	  while(q = iter_nbs()) {
		GCoord qc = q->coords();
		double qr = get_radius_from_name(q->name(), nar);
		if ((qc-pc).length() < 2*radius+pr+qr) {
		  if (q != p){
			nb_coords.push_back(qc);
			nb_radii.push_back(qr);
		  }
		}
	  }
	  // Place test particles on surface.
	  std::vector<GCoord>::iterator iter_qc;
	  for (int i=0; i<M; i++) {
		GCoord c = random_points[i] * (pr+radius) + pc;
		bool on_surface = true;
		for (unsigned int j=0; j<nb_coords.size(); j++) {
		  double distance = nb_coords[j].distance(c, box);
		  if (distance < nb_radii[j] + radius) {
			on_surface = false;
			break;
		  }
		}
		if (on_surface)
		  num_on_surface++;
	  }
	  if (DEBUG_SASA) {
		std::cout << num_on_surface << " of " << M << " test particles lie on the surface." << std::endl;
	  }
	  area += 1.0*num_on_surface/M*4*pi*(pr+radius)*(pr+radius);
	  if (DEBUG_SASA) {
		std::cout << "Radius is " << pr+radius << " and total area is " << area << std::endl;
	  }
	}
	delete [] random_points;
	return area;
  }

  /*
	Make a Binning of atoms in `everything`.  Uses the maximum
	VdW radius + the probe radius to get atom neighbours.
	Uses `nar` to get the VdW radii for the atoms.
   */
  Binning make_binning(const AtomicGroup &everything, const double radius,
					   const std::vector<name_and_radius> &nar) {
	// Get maximum radius of an atom.
	double max_radius = 0;
	pAtom a;
	AtomicGroup::Iterator itera(everything);
	while (a = itera()) {
	  double rad = get_radius_from_name(a->name(), nar);
	  if (rad > max_radius)
		max_radius = rad;
	}
	double total_rad = max_radius + radius;
	return Binning(2*total_rad, everything);
  }
  
  /* See header file for description. */
  double calculate_sasa(const AtomicGroup& selected,
						const AtomicGroup& everything,
						const double radius, const int M) {
	std::vector<name_and_radius> nar = get_default_names_and_radii();
	Binning binning = make_binning(everything, radius, nar);
	return calculate_sasa_with_binning_nar(selected, radius, binning, nar, M);
  }

}
