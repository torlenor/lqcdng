/*
 * su3.h - su3 matrix handling
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

Su3Matrix::Su3Matrix() {
  matrix_.resize(kMatrixDim);
  for (unsigned int d=0;d<kMatrixDim;d++) {
    matrix_[d].resize(kMatrixDim);
  }
}

Su3Matrix::~Su3Matrix() {

}

void Su3Matrix::set(const unsigned int x, const unsigned int y, const std::complex<double> in) {
  matrix_[x][y] = in;
}

std::complex<double> Su3Matrix::get(const unsigned int x, const unsigned int y) {
  return matrix_[x][y];
}

void add_m(Su3Matrix &m_in1, Su3Matrix &m_in2, Su3Matrix &m_out) {
	for (unsigned int x=0; x<3; x++) {
		for (unsigned int y=0; y<3; y++) {
			m_out.set(x, y, m_in1.get(x, y) + m_in2.get(x, y));
		}
	}
}

void mult_m(Su3Matrix &m_in1, Su3Matrix &m_in2, Su3Matrix &m_out) {
  std::complex<double> multvar = 0;
	for (unsigned int x=0; x<3; x++) {
		for (unsigned int y=0; y<3; y++) {
      multvar = 0;
		  for (unsigned int k=0;k<3; k++) {
        multvar = m_in1.get(x,k)*m_in2.get(k,y);
      }
			m_out.set(x, y, multvar);
		}
	}
}

