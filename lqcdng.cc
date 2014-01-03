/*
 * lqcdnq.cc - Lattice QCD Monte Carlo main functions
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

#include <cstdlib>
#include <ctime>
#include <iostream>
#include <vector>

#include "globalsettings.h"
#include "helper.h"
#include "init.h"
#include "simulationkernels/puregauge.h"
#include "storage.hpp"
#include "su3.h"

int main(int argc, char **argv) {
  // Print help and read settings from command line
  GlobalSettings settings;
  Init(argc, argv, settings);
  PrintSettings(settings);
  std::cout << std::endl;

  LatticeStorage<Su3Matrix, double> storage1(4,4,4,4,true,true);

  // LatticeStorage<std::vector<double>, double> storage1(4,4,4,4,true,true);
  /*for(unsigned int i=0; i<4*4*4*4; i++) {
    for(unsigned int mu=0; mu<4; mu++) {
      storage1.at(i, mu).resize(3*3);
    }
  }*/

  storage1.at(1,1).set(1,1,4.0);
  // std::cout << "Site at 1: " << storage1.at(1) << std::endl;
  std::cout << "Link at 1,2: " << storage1.at(1,1).get(1,1) <<  std::endl;
  std::cout << std::endl;

  // Initialize random number generator
  // TODO: Change this to something else!
  srand (time(NULL));

  // Create a simulation instance
  PureGaugeSim *sim1 = new PureGaugeSim(settings);
  sim1->StartSimulation();

  return 0;  
}
