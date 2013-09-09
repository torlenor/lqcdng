/*
 * simulation.cc - main simulation functions
 *
 * Copyright © 2013 H.-P. Schadler  <hanspeter.schadler@uni-graz.at>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 */

#include <iostream>

#include "su3.h"
#include "globalsettings.h"
#include "helper.h"
#include "simulation.h"

MCSimulation::MCSimulation(GlobalSettings &settings) {
  settings_ = settings;
}

MCSimulation::~MCSimulation() {

}

void MCSimulation::PrepareNeib() {
  // Allocate memory for the neib array
  neib_.resize(settings_.nsites);
  for (int ns=0; ns<settings_.nsites; ns++){
    neib_[ns].resize(2*settings_.dim);
  }

  // Fills the neib array
  int i1p,i2p,i3p,i4p,i1m,i2m,i3m,i4m,is,isp1,isp2,isp3,isp4,ism1,ism2,ism3,ism4;
  for (int i1 = 0; i1<settings_.ns; i1++) {
    i1p = i1 + 1;
    i1m = i1 - 1;
    if (i1p == settings_.ns) i1p = 0;
    if (i1m == -1) i1m = settings_.ns-1;

    for (int i2 = 0; i2<settings_.ns; i2++) {
      i2p = i2 + 1;
      i2m = i2 - 1;
      if (i2p == settings_.ns) i2p = 0;
      if (i2m == -1) i2m = settings_.ns-1;

      for (int i3 = 0; i3<settings_.ns; i3++) {
      i3p = i3 + 1;
      i3m = i3 - 1;
      if (i3p == settings_.ns) i3p = 0;
      if (i3m == -1) i3m = settings_.ns-1;

      for (int i4 = 0; i4<settings_.nt; i4++) {
        i4p = i4 + 1;
        i4m = i4 - 1;
        if (i4p == settings_.nt) i4p = 0;
        if (i4m == -1) i4m = settings_.nt-1;

        // Compute the site address and the addresses of the sites shifted
        // by one unit in each direction

        is = i1 + i2*settings_.ns + i3*settings_.ns*settings_.ns + i4*settings_.ns*settings_.ns*settings_.ns;

        isp1 = i1p + i2*settings_.ns + i3*settings_.ns*settings_.ns + i4*settings_.ns*settings_.ns*settings_.ns;
        isp2 = i1 + i2p*settings_.ns + i3*settings_.ns*settings_.ns + i4*settings_.ns*settings_.ns*settings_.ns;
        isp3 = i1 + i2*settings_.ns + i3p*settings_.ns*settings_.ns + i4*settings_.ns*settings_.ns*settings_.ns;
        isp4 = i1 + i2*settings_.ns + i3*settings_.ns*settings_.ns + i4p*settings_.ns*settings_.ns*settings_.ns;

        ism1 = i1m + i2*settings_.ns + i3*settings_.ns*settings_.ns + i4*settings_.ns*settings_.ns*settings_.ns;
        ism2 = i1 + i2m*settings_.ns + i3*settings_.ns*settings_.ns + i4*settings_.ns*settings_.ns*settings_.ns;
        ism3 = i1 + i2*settings_.ns + i3m*settings_.ns*settings_.ns + i4*settings_.ns*settings_.ns*settings_.ns;
        ism4 = i1 + i2*settings_.ns + i3*settings_.ns*settings_.ns + i4m*settings_.ns*settings_.ns*settings_.ns;

        // Fill the neib array

        neib_[is][0] = isp1;
        neib_[is][1] = isp2;
        neib_[is][2] = isp3;
        neib_[is][3] = isp4;

        neib_[is][4] = ism1;
        neib_[is][5] = ism2;
        neib_[is][6] = ism3;
        neib_[is][7] = ism4;
        }
      }
    }
  }
}

void MCSimulation::PrepareStorage() {
  for (int i=0; i<settings_.dim*settings_.nsites; i++) {
    lattice_.push_back(new Su3Matrix());
  }
 
  // Set all links to unit matrix 
  for (std::vector<Su3Matrix*>::iterator iter=lattice_.begin(); iter != lattice_.end(); ++iter) {
    for (unsigned int i=0; i<3; i++) {
      (*iter)->set(i,i,1.0);
    }   
  }
}

void MCSimulation::DeleteStorage() {
  for (std::vector<Su3Matrix*>::iterator iter=lattice_.begin(); iter != lattice_.end(); ++iter) {
    delete *iter;
  }
}

void MCSimulation::Update(const int nskip) {
  for (int skip=0; skip<nskip; skip++) {  
    // TODO: Update procedures
    std::cout << "Skip " << skip << std::endl;
  }
}

int MCSimulation::StartSimulation() {
  PrepareNeib();
  PrepareStorage();

  // Equilibration
  std::cout << "Equi..." << std::endl;
  Update(settings_.nequi);;

  for (int n=0; n<settings_.nmeas; n++) {
    Update(settings_.nskip);

    std::cout << "Meas " << n << std::endl;
  }

  DeleteStorage();

  return 0;
}

