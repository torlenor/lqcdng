/*
 * helper.cc - helper functions
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

void PrintMatrix(Su3Matrix &in){
  for (unsigned int i=0;i<3;i++) {
    for (unsigned int j=0;j<3;j++) {
      std::cout << in.get(i,j) << " ";
    }   
    std::cout << std::endl;
  }
}

void PrintSettings(GlobalSettings &settings) {
  std::cout << "Settings:" << std::endl
  << "Ns = " << settings.ns << std::endl
  << "Nt = " << settings.nt << std::endl
  << "Equilibrations = " << settings.nequi << std::endl
  << "Measurements = " << settings.nmeas << std::endl
  << "Skip = " << settings.nskip << std::endl
  << "Beta = " << settings.beta << std::endl;
}

