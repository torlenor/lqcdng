/*
 * puregaugesim.cc - class for pure SU(3) gauge sim inherited from
 * GenericSimClass
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

#include "puresu2gauge.h"

#include <cstdlib>
#include <ctime>
#include <iostream>
#include <limits>
#include <sstream>

#include "globalsettings.h"
#include "su2.h"

void PureSU2GaugeSim::StapleSum(Su2Matrix &S, int mu,int x) {
  int xpmu,xpnu,xmnupmu,xmnu;
  Su2Matrix U12, U123;

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

void PureSU2GaugeSim::MetroOffer(Su2Matrix &Unew, Su2Matrix &Uold) {
  // Generates a trial matrix for the metropolis update
  Su2Matrix change;
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

void PureSU2GaugeSim::OverOffer(Su2Matrix &Unew, Su2Matrix &Uold, Su2Matrix &stot) {
  // Generates a overrelaxation trial matrix
  Su2Matrix gtem, g0;

  g0 = stot;

  g0.Norm();

  MultMatrixadagbdagc(Uold, g0, gtem);
  MultMatrixadagbc(g0, gtem, Unew);
}

void PureSU2GaugeSim::Mixed() {
  Su2Matrix uold, ustap0;
  Su2Matrix utrial, udiff, betasum;
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

std::complex<double> PureSU2GaugeSim::CalcPoll() {
  // Calculates the Polyakov loop spatial average
  std::complex<double> poll, trace;
  poll=std::complex<double> (0,0);
  trace=std::complex<double> (0,0);

  Su2Matrix up, uu, upaux;

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

  poll = poll/((double)2.0*settings_.ns*settings_.ns*settings_.ns);

  return poll;
}

std::complex<double> PureSU2GaugeSim::CalcPlaq() {
  // Calculates the plaquette spatial average over all 6N plaquettes
  int ispmu, ispnu;
  std::complex<double> sumplaqs = std::complex<double>(0,0), trace;
  Su2Matrix *u1, *u2, *u3, *u4, u23, u234;

  for (int is = 0; is<settings_.nsites; is++) {
    for (int imu = 0; imu<settings_.dim; imu++) {
      for (int inu = imu+1; inu<settings_.dim; inu++) {
        ispmu = neib_[is][imu];
        ispnu = neib_[is][inu];

        u1 = lattice_[is][imu];
        u2 = lattice_[ispmu][inu];
        u3 = lattice_[ispnu][imu];
        u4 = lattice_[is][inu];

        MultMatrixabdagc(*u2, *u3, u23);
        MultMatrixabdagc(u23, *u4, u234);
        trace = MultTraceMatrix(*u1,u234);

        sumplaqs = sumplaqs + trace;
      }
    }
  }

  sumplaqs=sumplaqs/((double)6*2*(double)settings_.nsites); // factor because sum over N lattice points and md from trace

  return sumplaqs;
}

void PureSU2GaugeSim::Update(const int nskip) {
  for (int skip=0; skip<nskip; skip++) {  
    std::cout << "\r Sweep " << skip+1 << " of " << nskip << std::flush;
    Mixed();
  }
  std::cout << std::endl;
}

void PureSU2GaugeSim::Measurement() {
  //
}

void PureSU2GaugeSim::PrepareStorage() {
  lattice_.resize(settings_.nsites);
  for (int i=0; i<settings_.nsites; i++) {
    for(int d=0; d<2*settings_.dim; d++) {
      lattice_[i].push_back(new Su2Matrix());
    }
  }
 
  // Set all links to unit matrix 
  std::vector<std::vector<Su2Matrix*> >::iterator site_iter;
  std::vector<Su2Matrix*>::iterator direction_iter;

  for (site_iter=lattice_.begin(); site_iter != lattice_.end(); ++site_iter) {
    for (direction_iter=site_iter->begin(); direction_iter != site_iter->end(); ++direction_iter) {
      for (unsigned int i=0; i<3; i++) {
        (*direction_iter)->set(i,i,1.0);
      }   
    }
  }
}

void PureSU2GaugeSim::DeleteStorage() {
  std::vector<std::vector<Su2Matrix*> >::iterator site_iter;
  std::vector<Su2Matrix*>::iterator direction_iter;

  for (site_iter=lattice_.begin(); site_iter != lattice_.end(); ++site_iter) {
    for (direction_iter=site_iter->begin(); direction_iter != site_iter->end(); ++direction_iter) {
        delete *direction_iter;
    }
  }
}

void PureSU2GaugeSim::InitIndividual() {
  std::stringstream randletters;
  if (settings_.writeconf == true || settings_.meas == true) {
    // We are generating 3 random letters to mark the stream
    srand(time(NULL));
    char a;
    for (int i=0; i<3; i++) {
      a = char('a' + (rand() % 26));
      randletters << a;
    }
  }

  if (settings_.meas == true) {
    // Generate filename for measurements writeout
    std::stringstream filename;
    filename << "meas_" << settings_.ns << "x" << settings_.nt << "_b" << settings_.beta << "_";
    filename << randletters.str();
    filename << ".data";
    filemeasname_ = filename.str();

    filemeas_.open(filemeasname_.c_str());
    filemeas_.flags (std::ios::scientific);
    filemeas_.precision(std::numeric_limits<double>::digits10 + 1); 
    filemeas_ << "# Totalplaq Re(pl) Im(pl) Re(S/Beta) Spacelikeplaq Timelikeplaq" << std::endl;

    std::cout << "Observables will be written to: " << filemeasname_ << std::endl;
  }
  
  if (settings_.writeconf == true) {
    // Generate filename for config writeout
    std::stringstream filename;
    filename << "conf_" << settings_.ns << "x" << settings_.nt << "_b" << settings_.beta << "_";
    filename << randletters.str();
    filename << "_";
    fileconfnamebase_ = filename.str();

    std::cout << "Configurations will be written to: " << fileconfnamebase_ << "#ofconfig" << std::endl;
  }
}

void PureSU2GaugeSim::CleanupIndividual() {
  if (settings_.meas == true) {
    filemeas_.close();
  }
}

int PureSU2GaugeSim::WriteConfig(const int &m) {
  int elems=0;
  std::complex<double> plaq, poll;
  FILE* pFile;

  // int nindex=nsite*matrixdim*matrixdim*DIM;
  int nindex=settings_.nsites*2*2*4;

  std::stringstream fconfigname;
  fconfigname << fileconfnamebase_ << m << ".bin";

  pFile = fopen(fconfigname.str().c_str(), "wb");
  if (pFile == NULL) perror ("Error opening file");
  else{
    elems += fwrite(&settings_.ns, sizeof(int), 1, pFile);
    elems += fwrite(&settings_.ns, sizeof(int), 1, pFile);
    elems += fwrite(&settings_.ns, sizeof(int), 1, pFile);
    elems += fwrite(&settings_.nt, sizeof(int), 1, pFile);
    
    poll=CalcPoll();
    plaq=CalcPlaq();

    elems += fwrite(&poll, sizeof(std::complex<double>),1, pFile);
    elems += fwrite(&plaq, sizeof(std::complex<double>),1, pFile);

    for (int is=0; is<settings_.nsites; is++) {
      for (int mu=0; mu<settings_.dim; mu++) {
        elems += fwrite((&lattice_[is][mu]->at(0,0)), sizeof(std::complex<double>), 2*2 , pFile);
      }
    }

    fclose(pFile);
  }
  return nindex + 6 - elems;
}

void PureSU2GaugeSim::WriteMeas(const int &m) {
  std::complex<double> poll=CalcPoll();
  std::complex<double> plaq=CalcPlaq();
  filemeas_ << std::real(plaq) << " " << std::real(poll) << " " << std::imag(poll) << " " << 0 << " " << 0 << " " << 0 << std::endl; 
}