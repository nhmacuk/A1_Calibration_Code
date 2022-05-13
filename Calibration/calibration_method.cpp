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

#include "calibration.h"
#include "matrix_algo.h"


using namespace easy3d;



/**
 * TODO: Finish]ss, otherwise false. On success, the camera parameters are returned by
 */
bool Calibration::calibration(
        const std::vector<Vector3D>& points_3d, /// input: An array of 3D points.
        const std::vector<Vector2D>& points_2d, /// input: An array of 2D image points.
        double& fx, double& fy,    /// output: the focal length (in our slides, we use 'alpha' and 'beta'),
        double& cx, double& cy,    /// output: the principal point (in our slides, we use 'u0' and 'v0'),
        double& skew,              /// output: the skew factor ('-alpha * cot_theta')
        Matrix33& R,               /// output: the 3x3 rotation matrix encoding camera orientation.
        Vector3D& t)               /// output：a 3D vector encoding camera translation.
{
    /* std::cout << "\n[Liangliang]:\n"
                  "\tThe input parameters of this function are:\n"
                 "\t\t- points_3d: An array of 3D points (input to this function)\n"
                 "\t\t- points_2d: An array of 2D image points (input to this function)\n"
                 "\tThis function must return either 'true' on success or 'false' otherwise. On success, the camera\n"
                 "\tparameters are returned by the following variables:\n"
                 "\t\t- fx and fy: the focal lengths (in our slides, we use 'alpha' and 'beta')\n"
                 "\t\t- cx and cy: the principal point (in our slides, we use 'u0' and 'v0')\n"
                 "\t\t- skew:      the skew factor ('-alpha * cot_theta')\n"
                 "\t\t- R:         the 3x3 rotation matrix encoding camera orientation\n"
                 "\t\t- t:         a 3D vector encoding camera location.\n"
                 "\tIMPORTANT: don't forget to write your recovered parameters to the above variables." << std::endl;
    */

    // check if input is valid (e.g., number of correspondences >= 6, sizes of 2D/3D points must match)
    int num_points = points_3d.size();
    if (num_points < 6) {
        std::cout << "ERROR: Incufficient number of points";
        return false;
    }

    if (points_2d.size() != points_3d.size()) {
        std:: cout << "ERROR: 2D and 3D points do not have same size\n";
        return false;
    }

    // creat P-matrix:
    int n = points_3d.size();
    Matrix P(2*n, 12, 0.0);
    for (int i = 0; i < P.rows(); i++) {
        if (i % 2==0) {
            int m = i/2;
            P.set_row(i, {points_3d[m][0], points_3d[m][1], points_3d[m][2], 1,
                          0, 0, 0, 0,
                          -points_2d[m][0]*points_3d[m][0], -points_2d[m][0]*points_3d[m][1],
                          -points_2d[m][0]*points_3d[m][2], -points_2d[m][0]*1});
        }
        if (i % 2==1) {
            int m = i/2;
            P.set_row(i, {0, 0, 0, 0,
                          points_3d[m][0], points_3d[m][1], points_3d[m][2], 1,
                          -points_2d[m][1]*points_3d[m][0], -points_2d[m][1]*points_3d[m][1],
                          -points_2d[m][1]*points_3d[m][2], -points_2d[m][1]*1});
        }
    }

    // std::cout<<"We're printing matrix: \n";
    //std::cout<<P<<std::endl;

    // solve for M (the whole projection matrix, i.e., M = K * [R, t]) using SVD decomposition.

    /// matrix-vector product
    int mm = P.rows(); int nn = P.cols();
    Matrix U(mm, mm, 0.0);   // initialized with 0s
    Matrix S(mm, nn, 0.0);   // initialized with 0s
    Matrix V(nn, nn, 0.0);   // initialized with 0s

    // Compute the SVD decomposition of A
    svd_decompose(P, U, S, V);

    //last col of V is m
    Vector m = V.get_column(V.cols()-1);

    // fill the M matrix
    Matrix M(3,4, 0.0);
    Vector3D b(3,1, 0.0);
    for (int i = 0; i < m.size(); i = i+4) {
        M.set_row(i/4, {m[i], m[i+1],m[i+2],m[i+3]});
        // last col of M is the vector b
        b[i/4] =  m[i+3];
        //rest  of M is matrix A
    }

    Vector3D a1(m[0], m[1], m[2]);
    Vector3D a2(m[4], m[5], m[6]);
    Vector3D a3(m[8], m[9], m[10]);

    double rho = 1 / length(a3);

    // Optional: you can check if your M is correct by applying M on the 3D points. If correct, the projected point
    // should be very close to your input images points.
    Matrix P_tick = (M/rho)*P_w;
    //normalization of P:
    for (int i = 0; i < P.cols(); i++) {
        P_tick[0][i] = P_tick[0][i]/P_tick[2][i];
        P_tick[1][i] = P_tick[1][i]/P_tick[2][i];
        P_tick[2][i] = P_tick[2][i]/P_tick[2][i];
    }
    for (int i = 0; i < P.cols(); i++) {
        P_tick[1][i] = P_tick[1][i]/P_tick[2][i];
        P_tick[2][i] = 1;
    }
//    std::cout << "P': " << P_tick << "\n";
    Matrix P_compare(3, 6,{points_2d[0][0], points_2d[1][0], points_2d[2][0], points_2d[3][0], points_2d[4][0], points_2d[5][0],
                     points_2d[0][1], points_2d[1][1], points_2d[2][1], points_2d[3][1], points_2d[4][1], points_2d[5][1],
                    1,1,1,1,1,1});
//    std::cout << "Pcomp: " << P_compare << "\n";
//    std::cout << "P comparison: " << P_compare-P_tick<< "\n";

    // extract intrinsic parameters from M.

    double uo = rho * rho * dot(a1, a3);
    double vo = rho * rho * dot(a2, a3);

    Vector a1a3 = cross(a1, a3);
    Vector a2a3 = cross(a2, a3);

    double upper = dot(a1a3,a2a3);
    double lower = length(a1a3)*length(a2a3);
    double theta = acos(-upper/lower);

    double alpha = rho*rho* length(a1a3)*sin(theta);
    double beta =  rho*rho* length(a2a3)*sin(theta);

    Matrix33 K(alpha, -alpha*(cos(theta)/sin(theta)), uo,
               0, beta/(sin(theta)), vo,
               0, 0, 1);

    // TODO: extract extrinsic parameters from M.
    Vector r1 = a2a3/length(a2a3);//G
    Vector r3 = rho*a3;//G
    Vector r2 = cross(r3,r1);//G
//    Matrix R(3, 3, 0.0);
    Matrix33 R_fake(r1[0],r1[1],r1[2], //G
                    r2[0],r2[1],r2[2],
                    r3[0],r3[1],r3[2]);
    //check:
    std::cout<< "R.T*R "<< R_fake*transpose(R_fake) << "\n";

    fx = alpha; fy = beta;
    cx = uo; cy = vo;
    skew = -alpha * (cos(theta) / sin(theta));
    R = R_fake;

    Vector3D t_fake = rho*inverse(K)*b;  // translation matrix
    t = t_fake ;

    return true;
}