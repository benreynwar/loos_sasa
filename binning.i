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

%include <std_vector.i>

 // Commented out intVector because it was getting defined somewhere else.
 //%template(intVector) std::vector<int>;

%{
#include <loos_defs.hpp>
#include <binning.hpp>

%}

namespace loos {

  class Binning {
  public:
	Binning(double cutoff, const AtomicGroup &group);
	~Binning();
	AtomicGroup* get_nbs(pAtom p) const;
  };

}
