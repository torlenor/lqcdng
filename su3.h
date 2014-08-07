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

    inline void set(const int i, const int j, const std::complex<double> in) {
      #ifdef DEBUG
        matrix_.at(i*3+j) = in; 
      #else
        matrix_[i*3+j] = in; 
      #endif
    }
    inline std::complex<double> get(const int i, const int j) {
      #ifdef DEBUG
        return matrix_.at(i*3+j);
      #else
        return matrix_[i*3+j];
      #endif
    }
    inline std::complex<double>& at(const int i, const int j) {
      #ifdef DEBUG
        return matrix_.at(i*3+j);
      #else
        return matrix_[i*3+j];
      #endif
    }
    

    void Norm(); // Normalizes the SU3 matrix

    void print(Su3Matrix &in); // Prints the matrix to STDOUT

  private:
    static const int kMatrixDim = 3;
    std::vector<std::complex<double> > matrix_;
};

inline void AddMatrix(Su3Matrix &m_in1, Su3Matrix &m_in2, Su3Matrix &m_out) {
	for (int x=0; x<3; x++) {
		for (int y=0; y<3; y++) {
			m_out.set(x, y, m_in1.get(x, y) + m_in2.get(x, y));
		}
	}
}

inline void SubstractMatrix(Su3Matrix &m_in1, Su3Matrix &m_in2, Su3Matrix &m_out) {
	for (int x=0; x<3; ++x) {
		for (int y=0; y<3; ++y) {
			m_out.at(x, y) = m_in1.get(x, y) - m_in2.get(x, y);
		}
	}
}

inline void MultMatrixabc(Su3Matrix &m_in1, Su3Matrix &m_in2, Su3Matrix &m_out) {
	for (int x=0; x<3; ++x) {
		for (int y=0; y<3; ++y) {
      m_out.at(x,y)=m_in1.at(x,0)*m_in2.at(0,y);
		  for (int k=1;k<3; ++k) {
        m_out.at(x,y) += m_in1.at(x,k)*m_in2.at(k,y);
      }
		}
	}
}

inline void MultMatrixadagbc(Su3Matrix &m_in1, Su3Matrix &m_in2, Su3Matrix &m_out) {
	for (int x=0; x<3; ++x) {
		for (int y=0; y<3; ++y) {
      m_out.at(x,y) = conj(m_in1.at(0,x))*m_in2.at(0,y);
		  for (int k=1;k<3; ++k) {
        m_out.at(x,y) += conj(m_in1.at(k,x))*m_in2.at(k,y);
      }
		}
	}
}

inline void MultMatrixabdagc(Su3Matrix &m_in1, Su3Matrix &m_in2, Su3Matrix &m_out) {
	for (int x=0; x<3; ++x) {
		for (int y=0; y<3; ++y) {
      m_out.at(x,y) = m_in1.at(x,0)*conj(m_in2.at(y,0)) + m_in1.at(x,1)*conj(m_in2.at(y,1))  + m_in1.at(x,2)*conj(m_in2.at(y,2));
		}
	}
}

inline void MultMatrixadagbdagc(Su3Matrix &m_in1, Su3Matrix &m_in2, Su3Matrix &m_out) {
	for (int x=0; x<3; ++x) {
		for (int y=0; y<3; ++y) {
      m_out.at(x,y) = conj(m_in1.at(0,x))*conj(m_in2.at(y,0));
		  for (int k=1;k<3; ++k) {
        m_out.at(x,y) += conj(m_in1.at(k,x))*conj(m_in2.at(y,k));
      }
		}
	}
}

inline std::complex<double> MultTraceMatrix(Su3Matrix &m_in1, Su3Matrix &m_in2) {
  std::complex<double> tr=0;
  for (int k=0; k<3; k++) {
		for (int i=0; i<3; i++) {
      tr += m_in1.at(i,k)*m_in2.at(k,i);
    }
  }

  return tr;
}

#endif // SU3_H
