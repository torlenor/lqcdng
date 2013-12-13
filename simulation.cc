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

#include "simulation.h"

#include <cstdlib>
#include <iostream>

#include "globalsettings.h"
#include "helper.h"
#include "su3.h"

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
  lattice_.resize(settings_.nsites);
  for (int i=0; i<settings_.nsites; i++) {
    for(int d=0; d<2*settings_.dim; d++) {
      lattice_[i].push_back(new Su3Matrix());
    }
  }
 
  // Set all links to unit matrix 
  std::vector<std::vector<Su3Matrix*> >::iterator site_iter;
  std::vector<Su3Matrix*>::iterator direction_iter;

  for (site_iter=lattice_.begin(); site_iter != lattice_.end(); ++site_iter) {
    for (direction_iter=site_iter->begin(); direction_iter != site_iter->end(); ++direction_iter) {
      for (unsigned int i=0; i<3; i++) {
        (*direction_iter)->set(i,i,1.0);
      }   
    }
  }
}

void MCSimulation::DeleteStorage() {
  std::vector<std::vector<Su3Matrix*> >::iterator site_iter;
  std::vector<Su3Matrix*>::iterator direction_iter;

  for (site_iter=lattice_.begin(); site_iter != lattice_.end(); ++site_iter) {
    for (direction_iter=site_iter->begin(); direction_iter != site_iter->end(); ++direction_iter) {
        delete *direction_iter;
    }
  }
}

double MCSimulation::Uni() {
  return ((double)rand() / ((double)RAND_MAX + 1.0));
}

void MCSimulation::StapleSum(Su3Matrix &S, int mu,int x) {
  int xpmu,xpnu,xmnupmu,xmnu;
  Su3Matrix U12, U123;

  xpmu = neib_[x][mu];

  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      S.set(i, j, 0);
    }   
  }
  
  for (int nu=0; nu<settings_.dim; nu++) {
    if (mu != nu) {
      xpnu = neib_[x][nu];
      
      MultMatrixabdagc(*lattice_[xpmu][nu], *lattice_[xpnu][mu], U12);
      MultMatrixabdagc(U12, *lattice_[x][nu], U123);
    
      AddMatrix(S, U123, S);

      xmnu = neib_[x][4+nu];
      xmnupmu = neib_[xmnu][mu];
    
      MultMatrixadagbdagc(*lattice_[xmnupmu][nu], *lattice_[xmnu][mu], U12);
      MultMatrixabc(U12, *lattice_[xmnu][nu], U123);

      AddMatrix(S, U123, S);
    }   
  }
}

void MCSimulation::MetroOffer(Su3Matrix &Unew, Su3Matrix &Uold) {
  // Generates a trial matrix for the metropolis update
  Su3Matrix change;
  double phi, rho, rx, ry; 
  const double bias = 3.3;
  
  for (int j=0; j<3-1; j++) {
    for (int kk=0; kk<3; kk++) {
      phi = 2*M_PI * Uni();
      rho = sqrt(-log(1.0-Uni()));
      rx = rho*cos(phi);
      ry = rho*sin(phi);
      change.set(j, kk, std::complex<double>(rx,ry));
    }   
    change.set(j, j, change.get(j,j) + bias);
  }

  change.Norm();
  
  if (Uni() < 0.5) { 
    MultMatrixabc(change, Uold, Unew);
  } else {
    MultMatrixadagbc(change, Uold, Unew);
  }    
}

void MCSimulation::OverOffer(Su3Matrix &Unew, Su3Matrix &Uold, Su3Matrix &stot) {
  // Generates a overrelaxation trial matrix
  Su3Matrix gtem, g0;

  g0 = stot;

  g0.Norm();

  MultMatrixadagbdagc(Uold, g0, gtem);
  MultMatrixadagbc(g0, gtem, Unew);
}

void MCSimulation::Mixed() {
  Su3Matrix uold, ustap0;
  Su3Matrix utrial, udiff, betasum;
  double b03 = -settings_.beta/(double)3.0; // MODIFIED FOR SU(N=md), md=2,3
  
  std::complex<double> trace;

  for (int n=0; n<settings_.nsites; n++) {
    for (int mu=0; mu<settings_.dim; mu++) {
      uold = *lattice_[n][mu];

      StapleSum(ustap0, mu, n);

      for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
          betasum.set(i, j, ustap0.get(i, j)*b03);
        }
      }

      for (int ihit=0; ihit<settings_.nhit; ihit++) {
        MetroOffer(utrial, uold);
        SubstractMatrix(uold, utrial, udiff);
        trace = MultTraceMatrix(udiff,betasum);
        if ( Uni() < exp(std::real(trace))){
          *lattice_[n][mu] = utrial;
          uold = utrial;
        }
      }

      // Perform one additional overrelaxation step
      if (settings_.beta>0.0001) {
        OverOffer(utrial, uold, betasum);
        SubstractMatrix(uold, utrial, udiff);
        trace = MultTraceMatrix(udiff,betasum);

        if (Uni() < exp(std::real(trace))) {
          *lattice_[n][mu] = utrial;
        }   
      }

    }   
  }
}

std::complex<double> MCSimulation::MeasPoll() {
  // Calculates the Polyakov loop spatial average
  std::complex<double> poll, trace;
  poll=std::complex<double> (0,0);
  trace=std::complex<double> (0,0);

  Su3Matrix up, uu, upaux;

  int is0=0;
  int i4=0;

  for(int i1=0;i1<settings_.ns;i1++){
    for(int i2=0;i2<settings_.ns;i2++){
      for(int i3=0;i3<settings_.ns;i3++){
        i4 = 0;
        is0 = i1 + i2*settings_.ns + i3*settings_.ns*settings_.ns + i4*settings_.ns*settings_.ns*settings_.ns;

        up = *lattice_[is0][3];

        for(i4=1;i4<settings_.nt-1;i4++){
          is0 = i1 + i2*settings_.ns + i3*settings_.ns*settings_.ns + i4*settings_.ns*settings_.ns*settings_.ns;
          uu = *lattice_[is0][3];
          MultMatrixabc(up, uu, upaux);
          up = upaux;
        }

        i4 = settings_.nt-1;
        is0 = i1 + i2*settings_.ns + i3*settings_.ns*settings_.ns + i4*settings_.ns*settings_.ns*settings_.ns;

        uu = *lattice_[is0][3];
        trace = MultTraceMatrix(up,uu);

        poll = poll + trace;
      }
    }
  }

  poll = poll/((double)3.0*settings_.ns*settings_.ns*settings_.ns);

  return poll;
}

void MCSimulation::Update(const int nskip) {
  for (int skip=0; skip<nskip; skip++) {  
    // TODO: Update procedures
    std::cout << "\r Sweep " << skip+1 << " of " << nskip << std::flush;
    Mixed();
  }
  std::cout << std::endl;
}

int MCSimulation::StartSimulation() {
  PrepareNeib();
  PrepareStorage();

  // Equilibration
  std::cout << "Equi..." << std::endl;
  Update(settings_.nequi);;

  for (int n=0; n<settings_.nmeas; n++) {
    // std::cout << "Meas " << n << std::endl;
    Update(settings_.nskip);
    if (settings_.meas == true) {
      std::cout << "Poll = " << MeasPoll() << std::endl;
    }
  }

  DeleteStorage();

  return 0;
}

