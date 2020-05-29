#include <iostream>
#include <vector>
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;

MatrixXd CalculateJacobian(const VectorXd &x_state);

int main()
{
    /**
   * Compute the Jacobian Matrix
   */

    // predicted state example
    // px = 1, py = 2, vx = 0.2, vy = 0.4
    VectorXd x_predicted(4);
    x_predicted << 1, 2, 0.2, 0.4;

    MatrixXd Hj = CalculateJacobian(x_predicted);

    cout << "Hj:" << endl
         << Hj << endl;

    return 0;
}

MatrixXd CalculateJacobian(const VectorXd &x_state)
{

    MatrixXd Hj(3, 4);
    // recover state parameters
    float px = x_state(0);
    float py = x_state(1);
    float vx = x_state(2);
    float vy = x_state(3);

    // check division by zero

    if (px == 0 || py == 0)
    {
        cout << "error: px or py cannot be zero" << endl;
        exit(0);
    }

    // compute the Jacobian matrix

    float px_sq_plus_py_sq = pow(px, 2) + pow(py, 2);

    float a11 = px / sqrt(px_sq_plus_py_sq);
    float a12 = py / sqrt(px_sq_plus_py_sq);
    float a13 = 0;
    float a14 = 0;

    float a21 = -py / px_sq_plus_py_sq;
    float a22 = px / px_sq_plus_py_sq;
    float a23 = 0;
    float a24 = 0;

    float a31 = py * (vx * py - vy * px) / pow(px_sq_plus_py_sq, 1.5);
    float a32 = px * (vy * px - vx * py) / pow(px_sq_plus_py_sq, 1.5);
    float a33 = a11;
    float a34 = a12;

    Hj << a11, a12, a13, a14,
        a21, a22, a23, a24,
        a31, a32, a33, a34;

    return Hj;
}