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

#include <iostream>

#include "init.h"
#include "su3.h"

using std::cout;
using std::endl;

unsigned int Ns=4, Nt=4;

void PrintMatrix(Su3Matrix &in){
  for (unsigned int i=0;i<3;i++){
    for (unsigned int j=0;j<3;j++){
      cout << in.get(i,j) << " ";
    }
    cout << endl;
  }
}

int main(int argc, char **argv) {
  std::vector<Su3Matrix*> lattice;

  for (unsigned int i=0;i<Ns*Ns*Ns*Nt;i++) {
    lattice.push_back(new Su3Matrix());
  }
  
  for(std::vector<Su3Matrix*>::iterator iter=lattice.begin(); iter != lattice.end(); ++iter) {
    for (unsigned int i=0;i<3;i++) {
      (*iter)->set(i,i,1.0);
    }
  }

  PrintMatrix( *lattice[0] );

  Su3Matrix total;

  for(std::vector<Su3Matrix*>::iterator iter=lattice.begin(); iter != lattice.end(); ++iter) {
    add_m( *(*iter),total,total);
  }

  PrintMatrix(total);

  for(std::vector<Su3Matrix*>::iterator iter=lattice.begin(); iter != lattice.end(); ++iter) {
    // Su3Matrix* m = *iter;
    delete *iter;
  }

  return 0;  
}
