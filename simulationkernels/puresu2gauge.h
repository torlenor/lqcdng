/*
 * puresu3gaugesim.h - class for pure SU(3) gauge sim inherited from
 * GenericSimClass - header
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

#ifndef PURESU2GAUGESIM_H
#define PURESU2GAUGESIM_H

#include <vector>

#include "genericsimclass.h"

#include "globalsettings.h"
#include "storage.hpp"
#include "su2.h"

class PureSU2GaugeSim : public GenericSimClass {
  public:
    PureSU2GaugeSim(GlobalSettings &settings) : GenericSimClass(settings){
      std::cout << "Pure SU(3) Simulation class version 0.1" << std::endl;
    }

  private:
    void Update(const int nskip);

    void StapleSum(Su3Matrix &S, int mu,int x);
    void OverOffer(Su3Matrix &Unew, Su3Matrix &Uold, Su3Matrix &stot);
    void MetroOffer(Su3Matrix &Unew, Su3Matrix &Uold);

    void Measurement();
    std::complex<double> MeasPoll();

    void Mixed();
};

#endif // PURESU2GAUGESIM_H
