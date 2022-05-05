/**
 * Copyright (C) 2015 by Liangliang Nan (liangliang.nan@gmail.com)
 * https://3d.bk.tudelft.nl/liangliang/
 *
 * This file is part of Easy3D. If it is useful in your research/work,
 * I would be grateful if you show your appreciation by citing it:
 * ------------------------------------------------------------------
 *      Liangliang Nan.
 *      Easy3D: a lightweight, easy-to-use, and efficient C++
 *      library for processing and rendering 3D data. 2018.
 * ------------------------------------------------------------------
 * Easy3D is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License Version 3
 * as published by the Free Software Foundation.
 *
 * Easy3D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef EASY3D_MATRIX_ALGOTHMS_H
#define EASY3D_MATRIX_ALGOTHMS_H


#include "matrix.h"

namespace easy3d {

    /**
     * Compute the determinant of a square matrix.
     * @param A The input matrix.
     * @return The determinant of A.
     */
    double determinant(const Matrix &A);


    /**
     * Compute the inverse of a square matrix. This is a wrapper around Eigen's inverse function.
     * @param A The input matrix.
     * @param invA The inverse of A.
     * @return false if failed (failure is reported only if the input matrix is not square).
     *      Upon return, invA carries the inverse of A.
     */
    bool inverse(const Matrix &A, Matrix &invA);


    /**
     * Compute the inverse of a square matrix. This is a wrapper around Eigen's inverse function.
     * @param A The input matrix.
     * @return The inverse of A.
     */
    Matrix inverse(const Matrix &A);


    /**
     * compute the singular value decomposition (svd) of an m by n matrix. this is a wrapper around eigen's jacobisvd.
     *
     * for an m-by-n matrix a, the singular value decomposition is an m-by-m orthogonal matrix u, an m-by-n diagonal
     * matrix s, and an n-by-n orthogonal matrix v so that a = u*s*v^t.
     *
     * the singular values, s[k] = s[k][k], are sorted in decreasing order, i.e., sigma[i] >= sigma[i+1] for any i.
     *
     * the singular value decomposition always exists, so the decomposition will never fail.
     *
     * @param a the input matrix to be decomposed, which can have an arbitrary size.
     * @param u the left side m by m orthogonal matrix.
     * @param s the middle m by n diagonal matrix, with zero elements outside of its main diagonal.
     * @param v the right side n by n orthogonal matrix v.
     *
     * @return upon return, u, s, and v carry the result of the svd decomposition.
     *
     * @attention v is returned (instead of v^t).
     */
    void svd_decompose(const Matrix &A, Matrix &U, Matrix &S, Matrix &V);


    /**
     * Solve a linear system (Ax=b) in the least squares sense.
     *
     * @param A The m-by-n (m >= n) coefficient matrix.
     * @param b The right-hand constant vector (m dimensional).
     * @param x The result of the system was successfully solved (m dimensional).
     * @return false if failed. If true, x carries the least-squares solution to the linear system.
     */
    bool solve_least_squares(const Matrix &A, const std::vector<double> &b, std::vector<double> &x);
}

#endif // EASY3D_MATRIX_ALGOTHMS_H
