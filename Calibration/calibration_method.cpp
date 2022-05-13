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
    /*
    std::cout << "\nTODO: I am going to implement the calibration() function in the following file:\n"
                 "\t    - calibration_method.cpp\n\n";
    */
    /*
    std::cout << "[Liangliang]:\n"
                 "\tCamera calibration requires computing the SVD and inverse of matrices.\n"
                 "\tIn this assignment, I provide you with a 'Matrix' and a 'Vector' data structures for storing and\n"
                 "\tmanipulating matrices and vectors of arbitrary sizes. I also wrote some code to show you how to:\n"
                 "\t    - compute the SVD of a matrix;\n"
                 "\t    - compute the inverse of a matrix;\n"
                 "\t    - compute the transpose of a matrix.\n\n"
                 "\tFeel free to use any of the provided data structures and functions. The commonly used linear algebra\n"
                 "\tfunctions are provided in the following files:\n"
                 "\t    - Calibration/matrix.h  Matrices of arbitrary dimensions and related functions.\n"
                 "\t    - Calibration/vector.h  Vectors of arbitrary dimensions and related functions.\n"
                 "\t    - Calibration/matrix_algo.h  Determinant, inverse, SVD, linear least-squares...\n"
                 "\tPlease refer to the above files for a complete list of useful functions and their usage.\n\n"
                 "\tIn your final submission, please\n"
                 "\t    - delete ALL unrelated test or debug code and avoid unnecessary output.\n"
                 "\t    - include all the source code (and please do NOT modify the structure of the directories).\n"
                 "\t    - do NOT include the 'build' directory (which contains the intermediate files in a build step).\n"
                 "\t    - make sure your code compiles and can reproduce your results without ANY modification.\n\n" << std::flush;
    */
    /*
    /// Below are a few examples showing some useful data structures and functions.

    // This is a 1D array of 'double' values. Alternatively, you can use 'double mat[25]' but you cannot change it
    // length. With 'std::vector', you can append/delete/insert elements, and much more. The 'std::vector' can store
    // not only 'double', but also any other types of objects. In case you may want to learn more about 'std::vector'
    // check here: https://en.cppreference.com/w/cpp/container/vector
    std::vector<double> array = {1, 3, 3, 4, 7, 6, 2, 8, 2, 8, 3, 2, 4, 9, 1, 7, 3, 23, 2, 3, 5, 2, 1, 5, 8, 9, 22};
    array.push_back(5); // append 5 to the array (so the size will increase by 1).
    array.insert(array.end(), 10, 3);  // append ten 3 (so the size will grow by 10).

    /// To access the value of an element.
    double a = array[2];

    /// define a 2D vector/point
    Vector2D b(1.1, 2.2);

    /// define a 3D vector/point
    Vector3D c(1.1, 2.2, 3.3);

    /// get the Cartesian coordinates of a (a is treated as Homogeneous coordinates)
    Vector2D p = c.cartesian();

    /// get the Homogeneous coordinates of p
    Vector3D q = p.homogeneous();

    /// the length of a vector
    double len = p.length();
    /// the squared length of a vector
    double sqr_len = p.length2();

    /// the dot product of two vectors
    double dot_prod = dot(p, q);

    /// the cross product of two vectors
    Vector cross_prod = cross(c, q);

    /// normalize this vector
    cross_prod.normalize();

    // Define an m-by-n double valued matrix.
    // Here I use the above array to initialize it. You can also use A(i, j) to initialize/modify/access its elements.
    const int m = 6, n = 5;
    Matrix A(m, n, array.data());    // 'array.data()' returns a pointer to the array.
//    std::cout << "M: \n" << A << std::endl;

    /// define a 3 by 4 matrix (and all elements initialized to 0.0)
    Matrix M(3, 4, 0.0);

    /// set first row by a vector
    M.set_row(0, Vector4D(1.1, 2.2, 3.3, 4.4));

    /// set second column by a vector
    M.set_column(1, Vector3D(5.5, 5.5, 5.5));

    /// define a 3 by 3 matrix (and all elements initialized to 0.0)
    Matrix33 B;

    /// define and initialize a 3 by 3 matrix
    Matrix33 T(1.1, 2.2, 3.3,
               0, 2.2, 3.3,
               0, 0, 1);

    /// define and initialize a 3 by 4 matrix
    Matrix34 P(1.1, 2.2, 3.3, 0,
               0, 2.2, 3.3, 1,
               0, 0, 1, 1);

    /// define a 15 by 9 matrix (and all elements initialized to 0.0)
    Matrix W(15, 9, 0.0);
    /// set the first row by a 9-dimensional vector
    W.set_row(0, {0, 1, 2, 3, 4, 5, 6, 7, 8}); // {....} is equivalent to a std::vector<double>

    /// get the number of rows.
    int num_rows = W.rows();

    /// get the number of columns.
    int num_cols = W.cols();

    /// get the the element at row 1 and column 2
    double value = W(1, 2);

    /// get the last column of a matrix
    Vector last_column = W.get_column(W.cols() - 1);

    /// define a 3 by 3 identity matrix
    Matrix33 I = Matrix::identity(3, 3, 1.0);

    /// matrix-vector product
    Vector3D v = M * Vector4D(1, 2, 3, 4); // M is 3 by 4

    Matrix U(m, m, 0.0);   // initialized with 0s
    Matrix S(m, n, 0.0);   // initialized with 0s
    Matrix V(n, n, 0.0);   // initialized with 0s

    // Compute the SVD decomposition of A
    svd_decompose(A, U, S, V);

    // Now let's check if the SVD result is correct

    // Check 1: U is orthogonal, so U * U^T must be identity
//    std::cout << "U*U^T: \n" << U * transpose(U) << std::endl;

    // Check 2: V is orthogonal, so V * V^T must be identity
//    std::cout << "V*V^T: \n" << V * transpose(V) << std::endl;

    // Check 3: S must be a diagonal matrix
//    std::cout << "S: \n" << S << std::endl;

    // Check 4: according to the definition, A = U * S * V^T
//    std::cout << "M - U * S * V^T: \n" << A - U * S * transpose(V) << std::endl;

    // Compute the inverse of a matrix
    Matrix invT;
    inverse(T, invT);
    // Let's check if the inverse is correct
//    std::cout << "B * invB: \n" << B * invB << std::endl;

    // TODO: the above code just demonstrates some useful data structures and APIs. Please remove all above code in your
    //       final submission.
    */
    //--------------------------------------------------------------------------------------------------------------
    // implementation starts ...

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

    // TODO: check if input is valid (e.g., number of correspondences >= 6, sizes of 2D/3D points must match)
    std::cout << '\n' << "OUR IMPLEMENTATION STARTS HERE :)" << '\n';

    int num_points = points_3d.size();
//    std::cout << "\n total amount of points: " << num_points << '\n';
//    std::cout << "\n amount of points on the first row: " << points_3d[0].size() << '\n';

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

    // TODO: solve for M (the whole projection matrix, i.e., M = K * [R, t]) using SVD decomposition.
    //   Optional: you can check if your M is correct by applying M on the 3D points. If correct, the projected point
    //             should be very close to your input images points.

    /// matrix-vector product
    int mm = P.rows();
    int nn = P.cols();
    Matrix U(mm, mm, 0.0);   // initialized with 0s
    Matrix S(mm, nn, 0.0);   // initialized with 0s
    Matrix V(nn, nn, 0.0);   // initialized with 0s

    // Compute the SVD decomposition of A

/*
* @param a the input matrix to be decomposed, which can have an arbitrary size.
* @param u the left side m by m orthogonal matrix.
* @param s the middle m by n diagonal matrix, with zero elements outside of its main diagonal.
* @param v the right side n by n orthogonal matrix v.
*/
    svd_decompose(P, U, S, V);
//    std::cout<<V<<'\n';

    // std::cout<<"We're printing matrix V: \n";
    // std::cout<<V.get_column(V.cols()-1)<<'\n';

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

    Vector3D a1(m[0],m[1],m[2]);//G
    Vector3D a2(m[4],m[5],m[6]);//G
    Vector3D a3(m[8],m[9],m[10]);//G


    double rho = 1 / length(a3);//G // TODO `pos or neg decide

    double uo = rho*rho* dot(a1,a3);//G
    double vo = rho*rho* dot(a2,a3);//G

    // TODO: extract intrinsic parameters from M.
    Vector a1a3 = cross(a1,a3);
    Vector a2a3 = cross(a2,a3);

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


    // output variables:
    //    double& fx, double& fy,    /// output: the focal length (in our slides, we use 'alpha' and 'beta'),
    //    double& cx, double& cy,    /// output: the principal point (in our slides, we use 'u0' and 'v0'),
    //    double& skew,              /// output: the skew factor ('-alpha * cot_theta')
    //    Matrix33& R,               /// output: the 3x3 rotation matrix encoding camera orientation.
    //    Vector3D& t)               /// output：a 3D vector encoding camera translation.

    fx = alpha;
    fy = beta;
    cx = uo;
    cy = vo;
    skew = -alpha * (cos(theta) / sin(theta));
    //skew = 0; //TODO discuss
    R = R_fake;

    Vector3D t_fake = rho*inverse(K)*b;  // translation matrix
    t = t_fake ;

    /// add calculate the coordinates to check:
    Matrix P_tick = M * P_w;
    std::cout << "P': " << P_tick << "\n";
    for (int i = 0; i < P.cols(); i++) {
        P_tick[0][i] = P_tick[0][i]/P_tick[2][i];
//        P_tick[1][i] = P_tick[1][i]/P_tick[2][i];
    }
    std::cout << "rho: " << rho << "\n";
    std::cout << "fx: " << fx << "\n";
    std::cout << "fy: " << fy << "\n";
    std::cout << "cx: " << cx << "\n";
    std::cout << "cy: " << cy
              << "\n";                //    std::cout << "\n\tTODO: After you implement this function, please return 'true' - this will trigger the viewer to\n";;
    std::cout << "sk: " << skew
              << "\n";        //                 "\t\tupdate the rendering using your recovered camera parameters. This can help you to visually check\n"
    std::cout << "R : " << R
              << "\n";   //                 "\t\tif your calibration is successful or not.\n\n" << std::flush;
    std::cout << "t : " << t << "\n";
    std::cout << "P': " << P_tick << "\n";
    std::cout << "Pw: " << P_w << "\n";
    std::cout << "K: " << K << "\n";
    std::cout << "M: " << M << "\n";
    return true;
}