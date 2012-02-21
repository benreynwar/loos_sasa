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
#include <vector>
#include <iostream>
#include <algorithm>

#include <binning.hpp>
#include <loos_defs.hpp>
#include <AtomicGroup.hpp>

#define DEBUG_BINNING 1

namespace loos {
  Binning::Binning(double cutoff_, const AtomicGroup& group_) :
	group(group_), cutoff(cutoff_) {
	if (DEBUG_BINNING) {
	  std::cout << "Starting Binning::Binning" << std::endl;
	}
	box = group.periodicBox();
	populate_bin_info();
	bin_nbs = new std::vector<bin>[total_bins];
	populate_bin_nbs();
	bin_atoms = new AtomicGroup[total_bins];
	populate_bin_atoms();
	bin_nb_atoms = new AtomicGroup[total_bins];
	populate_bin_nb_atoms();
  }

  Binning::~Binning() {
	delete [] bin_nbs;
	delete [] bin_atoms;
	delete [] bin_nb_atoms;
  }
	
  void Binning::populate_bin_info() {
	if (DEBUG_BINNING) {
	  std::cout << "Starting Binning::populate_bin_info" << std::endl;
	}
	// Work out the system size.
	if (group.isPeriodic()) {
	  // If it's periodic just use the the periodic box.
	  for (dim d=0; d<3; d++) {
		if (box[d] < cutoff) {
		  throw(std::runtime_error("Periodic box size is smaller than cutoff.  It is possible for multiple images to be within cutoff distance.  This code cannot handle that."));
		}
		mins[d] = -box[d]/2;
		maxs[d] = box[d]/2;
		periodic[d] = true;
	  }
	} else {
	  // If it's not periodic loop through the atoms to find
	  // the minimum and maximum coordinates.
	  for (dim d=0; d<3; d++) {
		periodic[d] = false;
	  }
	  AtomicGroup::Iterator iter(group);
	  pAtom p;
	  bool first[3];
	  for (dim d=0; d<3; d++) {
		first[d] = true;
	  }
	  while (p = iter()) {
		const GCoord c = p->coords();
		if (DEBUG_BINNING) {
		  std::cout << "populate_bin_info: c = " << c << std::endl;
		}
		for (dim d=0; d<3; d++) {
		  if (first[d]) {
			mins[d] = c[d];
			maxs[d] = c[d];
			first[d] = false;
		  } else {
			if (c[d] < mins[d]) {
			  mins[d] = c[d];
			}
			if (c[d] > maxs[d]) {
			  maxs[d] = c[d];
			}
		  }
		}
	  }
	}
	// Make the system a tiny bit bigger to make sure everything
	// lies comfortably inside.
	double widths[3];
	double magic_ratio = 1000;
	for (dim d=0; d<3; d++) {
	  widths[d] = maxs[d]-mins[d];
	  // Nicer to avoid zero widths.  Might create bugs.
	  if (widths[d] == 0) {
		widths[d] = cutoff/magic_ratio;
	  }
	  maxs[d] += widths[d]/magic_ratio;
	  mins[d] -= widths[d]/magic_ratio;
	}
	// Work out the bin sizes
	for (dim d=0; d<3; d++) {
	  nbins[d] = ceil((maxs[d] - mins[d])/cutoff);
	  if (DEBUG_BINNING) {
		std::cout << "populate_bin_info: d = " << d << " maxs[d] = " << maxs[d] << " mins[d] " << mins[d] << " cutoff = " << cutoff << " and nbins[d] = " << nbins[d] << std::endl;
	  }
	}
	total_bins = nbins[0]*nbins[1]*nbins[2];
	if (DEBUG_BINNING) {
	  std::cout << "populate_bin_info: Total bins = " << total_bins << std::endl;
	}
  }
  
  bin Binning::get_bin(pAtom p) const {
	GCoord c = p->coords();
	c.reimage(box);
	bin b[3];
	for (dim d=0; d<3; d++) {
	  b[d] = floor((c[d]-mins[d])/(maxs[d]-mins[d])*nbins[d]);
	  if (b[d] < 0) {
		throw(std::runtime_error("Bin should not be negative."));
	  } else if (b[d] >= nbins[d]) {
		throw(std::runtime_error("Bin out of range."));
	  }
	}
	return b[0]*nbins[1]*nbins[2]+b[1]*nbins[2]+b[2];
  }

  void Binning::populate_bin_atoms() {
	AtomicGroup::Iterator iter(group);
	pAtom p;
	while (p = iter()) {
	  bin_atoms[get_bin(p)] += p;
	}
  }

  AtomicGroup* Binning::get_nbs(pAtom p) const {
	return &bin_nb_atoms[get_bin(p)];
  }

  void Binning::populate_bin_nbs() {
	std::vector<bin> *x_nbs;
	// x_nbs[3*x+d] contains the indexes of the neighbours of bin x in dimension d.
	x_nbs = new std::vector<bin>[3*std::max(nbins[0], std::max(nbins[1], nbins[2]))];
	for (dim d=0; d<3; d++) {
	  for (bin x=0; x<nbins[d]; x++) {
		// This bin is always a nb of itself.
		x_nbs[3*x+d].push_back(x);
		// The bin one up is a nb unless we're on the edge.
		if (x<nbins[d]-1) {
		  x_nbs[3*x+d].push_back(x+1);
		} else if ((periodic[d]) && (nbins[d]>2)) {
		  x_nbs[3*x+d].push_back(0);
		}
		// The bin one down is a nb unless we're on the edge.
		if (x>0) {
		  x_nbs[3*x+d].push_back(x-1);
		} else if ((periodic[d]) && (nbins[d]>2)) {
		  x_nbs[3*x+d].push_back(nbins[d]-1);
		}
	  }
	}
	if (DEBUG_BINNING) {
	  for (dim d=0; d<3; d++) {
		for (bin x=0; x<nbins[d]; x++) {
		  std::cout << "d = " << d << " bin is " << x << " and num of nbs is " << x_nbs[3*x+d].size() << std::endl;
		}
	  }
	}	
	for (bin i=0; i<total_bins; i++) {
	  bin x = i/(nbins[1]*nbins[2]);
	  bin y = (i-x*nbins[1]*nbins[2])/nbins[2];
	  bin z = i-x*nbins[1]*nbins[2]-y*nbins[2];
	  for (std::vector<bin>::iterator itx = x_nbs[3*x].begin();
		   itx != x_nbs[3*x].end(); ++itx) {
		for (std::vector<bin>::iterator ity = x_nbs[3*y+1].begin();
			 ity != x_nbs[3*y+1].end(); ++ity) {
		  for (std::vector<bin>::iterator itz = x_nbs[3*z+2].begin();
			   itz != x_nbs[3*z+2].end(); ++itz) {
			bin inb = (*itx)*nbins[1]*nbins[2] + (*ity)*nbins[2] + (*itz);
			bin_nbs[i].push_back(inb);
		  }
		}
	  }
	  if (DEBUG_BINNING) {
	  	std::cout << "Length of bin_nbs["<<i<<"] is " << bin_nbs[i].size() << std::endl;
	  }
	}
	delete [] x_nbs;
  }	
  
  void Binning::populate_bin_nb_atoms() {
	for (bin b=0; b<total_bins; b++) {
	  for (std::vector<bin>::iterator nb = bin_nbs[b].begin();
		   nb != bin_nbs[b].end(); nb++) {
		bin_nb_atoms[b] += bin_atoms[*nb];
	  }
	}
  }

}

