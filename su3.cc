/*
 * su3.cc - su3 matrix handling
 *
 * Copyright © 2013 H.-P. Schadler  <hanspeter.schadler@uni-graz.at>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNmatrix_ General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOmatrix_T ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICmatrix_LAR Pmatrix_RPOSE.  See the
 * GNmatrix_ General Public License for more details.
 *
 * You should have received a copy of the GNmatrix_ General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 matrix_SA.
 *
 */

#include "su3.h"

#include <iostream>

Su3Matrix::Su3Matrix() {
  matrix_.resize(kMatrixDim);
  for (int d=0; d<kMatrixDim; d++) {
    matrix_[d].resize(kMatrixDim);
  }
}

Su3Matrix::~Su3Matrix() {

}

void Su3Matrix::Norm() {
  // Normalizes the first row of matrix matrix_ for SU(3)
  double normu, normv;
  std::complex<double> v[3];
  std::complex<double> vdotu;

  // Normalize the first row 
  normu = sqrt(real(matrix_[0][0]*conj(matrix_[0][0]) + matrix_[0][1]*conj(matrix_[0][1]) + matrix_[0][2]*conj(matrix_[0][2])));
  
  matrix_[0][0] = matrix_[0][0]/normu;
  matrix_[0][1] = matrix_[0][1]/normu;
  matrix_[0][2] = matrix_[0][2]/normu;

  // Build a vector v orthonormal to row 1 (Gramm-Schmidt method) and use
  // it as row 2
  vdotu = matrix_[1][0]*conj(matrix_[0][0])+matrix_[1][1]*conj(matrix_[0][1])+matrix_[1][2]*conj(matrix_[0][2]);

  v[0] = matrix_[1][0]-matrix_[0][0]*vdotu;
  v[1] = matrix_[1][1]-matrix_[0][1]*vdotu;
  v[2] = matrix_[1][2]-matrix_[0][2]*vdotu;

  normv = sqrt(real(v[0]*conj(v[0]) + v[1]*conj(v[1]) + v[2]*conj(v[2])));

  matrix_[1][0] = v[0]/normv;
  matrix_[1][1] = v[1]/normv;
  matrix_[1][2] = v[2]/normv;

  // The third row is cross product of row1* x row2* 
  matrix_[2][0] = conj(matrix_[0][1])*conj(matrix_[1][2])-conj(matrix_[0][2])*conj(matrix_[1][1]);
  matrix_[2][1] = conj(matrix_[0][2])*conj(matrix_[1][0])-conj(matrix_[0][0])*conj(matrix_[1][2]);
  matrix_[2][2] = conj(matrix_[0][0])*conj(matrix_[1][1])-conj(matrix_[0][1])*conj(matrix_[1][0]);
}

void Su3Matrix::print(Su3Matrix &in){
  for (unsigned int i=0;i<3;i++) {
    for (unsigned int j=0;j<3;j++) {
      std::cout << in.get(i,j) << " ";
    }   
    std::cout << std::endl;
  }
}