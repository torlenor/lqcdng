/*
 * puresu3gauge.cc - Pure SU(3) Monte Carlo simulation
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
#include "init.h"
#include "measuretime.h"
#include "simulationkernels/puresu3gauge.h"

int main(int argc, char **argv) {
  double tstart, tend;
  tstart = gettime();
  // Print help and read settings from command line
  GlobalSettings settings;
  Init(argc, argv, settings);
  settings.PrintSettings();
  std::cout << std::endl;

  // Initialize random number generator
  // TODO: Change this to something else!
  srand (time(NULL));

  // Create a simulation instance
  PureSU3GaugeSim *sim1 = new PureSU3GaugeSim(settings);
  sim1->StartSimulation();

  tend = gettime();
  
  std::cout << std::endl << "Total run time: " << (tend - tstart)/(double)60 << " minutes." << std::endl;
  std::cout << std::endl << "Have a nice day!" << std::endl;

  return 0;  
}
