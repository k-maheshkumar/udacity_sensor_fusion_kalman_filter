#include <iostream>
#include "Eigen/Dense"
#include "ukf.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

int main()
{

    // Create a UKF instance
    UKF ukf;

    /**
   * Programming assignment calls
   */
    MatrixXd Xsig = MatrixXd(5, 11);
    // ukf.GenerateSigmaPoints(&Xsig);
    // ukf.AugmentedSigmaPoints(&Xsig);

    // MatrixXd Xsig_pred = MatrixXd(15, 5);
    // ukf.SigmaPointPrediction(&Xsig_pred);

    // // print result
    // std::cout << "Xsig_pred = " << std::endl
    //           << Xsig_pred << std::endl;

    // VectorXd x_pred = VectorXd(5);
    // MatrixXd P_pred = MatrixXd(5, 5);
    // ukf.PredictMeanAndCovariance(&x_pred, &P_pred);

    // VectorXd z_out = VectorXd(3);
    // MatrixXd S_out = MatrixXd(3, 3);
    // ukf.PredictRadarMeasurement(&z_out, &S_out);

    VectorXd x_out = VectorXd(5);
    MatrixXd P_out = MatrixXd(5, 5);
    ukf.UpdateState(&x_out, &P_out);

    return 0;
}