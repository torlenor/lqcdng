/*
 * su3.h - su3 matrix handling header
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

#ifndef SU3_H
#define SU3_H

#include<complex>
#include<vector>

class Su3Matrix {
	public:
		Su3Matrix();
		~Su3Matrix();

    std::complex<double> get(const int x, const int y);
    void set(const int x, const int y, const std::complex<double> in);

    void Norm(); // Normalizes the SU3 matrix

  private:
    static const int kMatrixDim = 3;
    std::vector<std::vector<std::complex<double> > > matrix_;
};

void AddMatrix(Su3Matrix &m_in1, Su3Matrix &m_in2, Su3Matrix &m_out);
void SubstactMatrix(Su3Matrix &m_in1, Su3Matrix &m_in2, Su3Matrix &m_out);
void MultMatrixabc(Su3Matrix &m_in1, Su3Matrix &m_in2, Su3Matrix &m_out);
void MultMatrixadagbc(Su3Matrix &m_in1, Su3Matrix &m_in2, Su3Matrix &m_out);
void MultMatrixabdagc(Su3Matrix &m_in1, Su3Matrix &m_in2, Su3Matrix &m_out);
void MultMatrixadagbdagc(Su3Matrix &m_in1, Su3Matrix &m_in2, Su3Matrix &m_out);
std::complex<double> MultTraceMatrix(Su3Matrix &m_in1, Su3Matrix &m_in2);

#endif // SU3_H
