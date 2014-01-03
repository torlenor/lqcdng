/*
 * storage.h - storage handling class using templates
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

#ifndef STORAGE_HPP
#define STORAGE_HPP

#include <stdlib.h>
#include <vector>
#ifdef DEBUG
#include <stdio.h> 
#include <iostream>
#endif

template <class TL, class TS>
class LatticeStorage{
  public:
    // Constructor/Deconstructor
    LatticeStorage(int len1, int len2, int len3, int len4, bool hasLinks, bool hasSites);
    ~LatticeStorage();

    TS& at(int is);
    TL& at(int is, int mu);

  private:
    std::vector<TL> linkstorage_;
    std::vector<TS> sitestorage_;

    int len1_;
    int len2_;
    int len3_;
    int len4_;
    int nsites_;

    bool hasLinks_;
    bool hasSites_;
};


template <class TL, class TS>
LatticeStorage<TL, TS>::LatticeStorage(int len1, int len2, int len3, int len4, bool hasLinks, bool hasSites){
  len1_=len1;
  len2_=len2;
  len3_=len3;
  len4_=len4;
  nsites_=len1_*len2_*len3_*len4_;

  hasLinks_=hasLinks;
  hasSites_=hasSites;

  if (hasLinks_==true) {
    linkstorage_.resize(4*nsites_);
  }
  if (hasSites_==true) {
    sitestorage_.resize(nsites_);
  }

#ifdef DEBUG
  std::cout << "INFO: Storage class instantiated!" << std::endl;
#endif
}

template <class TL, class TS>
LatticeStorage<TL, TS>::~LatticeStorage() {

}

template <class TL, class TS>
TL& LatticeStorage<TL, TS>::at(int is, int mu) {
  #ifdef DEBUG
    if (hasLinks_ == false) {
      fputs ("ERROR: Wanted to access links in a nonlink storage!\n",stderr);
      abort();
    }
  #endif
  return linkstorage_.at(is*mu);
}

template <class TL, class TS>
TS& LatticeStorage<TL, TS>::at(int is) {
  #ifdef DEBUG
    if (hasSites_ == false) {
      fputs ("ERROR: Wanted to access sites in a nonsite storage!\n",stderr);
      abort();
    }
  #endif
  return sitestorage_.at(is);
}

#endif // STORAGE_HPP
