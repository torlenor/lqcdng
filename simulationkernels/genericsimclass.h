/*
 * genericsimclass.cc - generic class for MC simulations - header
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

#ifndef GENERICSIMCLASS_H
#define GENERICSIMCLASS_H

#include <vector>

#include "globalsettings.h"
#include "storage.hpp"

class GenericSimClass {
  public:
    GenericSimClass(GlobalSettings &settings);
    ~GenericSimClass();

    int StartSimulation();

  protected:
    void PrepareNeib();
    
    virtual void PrepareStorage(){};
    virtual void DeleteStorage(){};

    virtual void InitIndividual();
    virtual void CleanupIndividual();

    double Uni();

    virtual void Update(const int nskip);
    
    virtual void Measurement();
    virtual void WriteMeas(const int &m);
    virtual int WriteConfig(const int &m);

    GlobalSettings settings_;
    std::vector<std::vector<int> > neib_;
};

#endif // GENERICSIMCLASS_H
