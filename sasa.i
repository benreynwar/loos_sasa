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

%include <std_string.i>
%include <std_vector.i>

%{
#include <loos_defs.hpp>
#include <sasa.hpp>
%}

%template(narVector) std::vector<loos::name_and_radius>;

namespace loos {

  class name_and_radius {
  public:
	std::string name;
	double radius;
	name_and_radius(std::string name_, double radius_) :
	  name(name_), radius(radius_) {};
	// Define this constructor for SWIG vector fun. 
	name_and_radius() :
	  name(" "), radius(-1) {};
  };

  double
  calculate_sasa(const AtomicGroup& selected, const AtomicGroup& everything,
				 const double cutoff, const int M);

  double
  calculate_sasa_with_binning_nar(const AtomicGroup& selected,
								  const double cutoff,
								  const Binning& binning,
								  const std::vector<name_and_radius>& nar,
								  const int M);

  GCoord random_sphere_point();

}


