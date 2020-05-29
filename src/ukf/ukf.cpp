#include "ukf.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

UKF::UKF()
{
    Init();
}

UKF::~UKF()
{
}

void UKF::Init()
{
}

/**
 * Programming assignment functions: 
 */
void UKF::GenerateSigmaPoints(MatrixXd *Xsig_out)
{

    // set state dimension
    int n_x = 5;

    // define spreading parameter
    double lambda = 3 - n_x;

    // set example state
    VectorXd x = VectorXd(n_x);
    x << 5.7441,
        1.3800,
        2.2049,
        0.5015,
        0.3528;

    // set example covariance matrix
    MatrixXd P = MatrixXd(n_x, n_x);
    P << 0.0043, -0.0013, 0.0030, -0.0022, -0.0020,
        -0.0013, 0.0077, 0.0011, 0.0071, 0.0060,
        0.0030, 0.0011, 0.0054, 0.0007, 0.0008,
        -0.0022, 0.0071, 0.0007, 0.0098, 0.0100,
        -0.0020, 0.0060, 0.0008, 0.0100, 0.0123;

    // create sigma point matrix
    MatrixXd Xsig = MatrixXd(n_x, 2 * n_x + 1);

    // calculate square root of P
    MatrixXd A = P.llt().matrixL();

    /**
   * Student part begin
   */

    //// set first column of sigma point matrix
    // Xsig.col(0) = x;

    // // set remaining sigma points
    // for (int i = 0; i < n_x; ++i)
    // {
    //     Xsig.col(i + 1) = x + sqrt(lambda + n_x) * A.col(i);
    //     Xsig.col(i + 1 + n_x) = x - sqrt(lambda + n_x) * A.col(i);
    // }

    float sigma_factor = sqrt(lambda + n_x);

    MatrixXd variance_offset_left = x.rowwise().replicate(5) + sigma_factor * A;
    MatrixXd variance_offset_right = x.rowwise().replicate(5) - sigma_factor * A;

    Xsig.block<5, 1>(0, 0) = x;
    Xsig.block<5, 5>(0, 1) = variance_offset_left;
    Xsig.block<5, 5>(0, 6) = variance_offset_right;

    /**
   * Student part end
   */

    // print result
    // std::cout << "Xsig = " << std::endl << Xsig << std::endl;

    // write result
    *Xsig_out = Xsig;
}

/**
 * Programming assignment functions: 
 */

void UKF::AugmentedSigmaPoints(MatrixXd *Xsig_out)
{

    // set state dimension
    int n_x = 5;

    // set augmented dimension
    int n_aug = 7;

    // Process noise standard deviation longitudinal acceleration in m/s^2
    double std_a = 0.2;

    // Process noise standard deviation yaw acceleration in rad/s^2
    double std_yawdd = 0.2;

    // define spreading parameter
    double lambda = 3 - n_aug;

    // set example state
    VectorXd x = VectorXd(n_x);
    x << 5.7441,
        1.3800,
        2.2049,
        0.5015,
        0.3528;

    // create example covariance matrix
    MatrixXd P = MatrixXd(n_x, n_x);
    P << 0.0043, -0.0013, 0.0030, -0.0022, -0.0020,
        -0.0013, 0.0077, 0.0011, 0.0071, 0.0060,
        0.0030, 0.0011, 0.0054, 0.0007, 0.0008,
        -0.0022, 0.0071, 0.0007, 0.0098, 0.0100,
        -0.0020, 0.0060, 0.0008, 0.0100, 0.0123;

    // create augmented mean vector
    VectorXd x_aug = VectorXd(7);

    // create augmented state covariance
    MatrixXd P_aug = MatrixXd(7, 7);

    // create sigma point matrix
    MatrixXd Xsig_aug = MatrixXd(n_aug, 2 * n_aug + 1);

    /**
   * Student part begin
   */

    // // create augmented mean state
    // x_aug.head(5) = x;
    // x_aug(5) = 0;
    // x_aug(6) = 0;

    // // create augmented covariance matrix
    // P_aug.fill(0.0);
    // P_aug.topLeftCorner(5, 5) = P;
    // P_aug(5, 5) = std_a * std_a;
    // P_aug(6, 6) = std_yawdd * std_yawdd;

    // // create square root matrix
    // MatrixXd L = P_aug.llt().matrixL();

    // // create augmented sigma points
    // Xsig_aug.col(0) = x_aug;
    // for (int i = 0; i < n_aug; ++i)
    // {
    //     Xsig_aug.col(i + 1) = x_aug + sqrt(lambda + n_aug) * L.col(i);
    //     Xsig_aug.col(i + 1 + n_aug) = x_aug - sqrt(lambda + n_aug) * L.col(i);
    // }

    // create augmented mean state
    x_aug.head(n_x) = x;

    // create augmented covariance matrix
    P_aug.topLeftCorner(5, 5) = P;

    MatrixXd Q = MatrixXd(2, 2);
    Q << std_a, 0,
        0, std_yawdd;

    P_aug.bottomRightCorner(2, 2) = Q * Q;

    // create square root matrix
    MatrixXd A = P_aug.llt().matrixL();

    // create augmented sigma points
    float sigma_factor = sqrt(lambda + n_aug);

    MatrixXd variance_offset_left = x_aug.rowwise().replicate(7) + sigma_factor * A;
    MatrixXd variance_offset_right = x_aug.rowwise().replicate(7) - sigma_factor * A;

    Xsig_aug.block<7, 1>(0, 0) = x_aug;
    Xsig_aug.block<7, 7>(0, 1) = variance_offset_left;
    Xsig_aug.block<7, 7>(0, 8) = variance_offset_right;

    /**
   * Student part end
   */

    // print result
    // std::cout << "Xsig_aug = " << std::endl
    //           << Xsig_aug << std::endl;

    // write result
    *Xsig_out = Xsig_aug;
}

/**
 * Programming assignment functions: 
 */

void UKF::SigmaPointPrediction(MatrixXd *Xsig_out)
{

    // set state dimension
    int n_x = 5;

    // set augmented dimension
    int n_aug = 7;

    // create example sigma point matrix

    // clang-format off

    MatrixXd Xsig_aug = MatrixXd(n_aug, 2 * n_aug + 1);
    Xsig_aug << 5.7441, 5.85768, 5.7441, 5.7441, 5.7441, 5.7441, 5.7441, 5.7441, 5.63052, 5.7441, 5.7441, 5.7441, 5.7441, 5.7441, 5.7441,
                1.38,   1.34566, 1.52806, 1.38, 1.38, 1.38, 1.38, 1.38, 1.41434, 1.23194, 1.38, 1.38, 1.38, 1.38, 1.38,
                2.2049, 2.28414, 2.24557, 2.29582, 2.2049, 2.2049, 2.2049, 2.2049, 2.12566, 2.16423, 2.11398, 2.2049, 2.2049, 2.2049, 2.2049,
                0.5015, 0.44339, 0.631886, 0.516923, 0.595227, 0.5015, 0.5015, 0.5015, 0.55961, 0.371114, 0.486077, 0.407773, 0.5015, 0.5015, 0.5015,
                0.3528, 0.299973, 0.462123, 0.376339, 0.48417, 0.418721, 0.3528, 0.3528, 0.405627, 0.243477, 0.329261, 0.22143, 0.286879, 0.3528, 0.3528,
                0, 0, 0, 0, 0, 0, 0.34641, 0, 0, 0, 0, 0, 0, -0.34641, 0,
                0, 0, 0, 0, 0, 0, 0, 0.34641, 0, 0, 0, 0, 0, 0, -0.34641;

    // clang-format on

    // create matrix with predicted sigma points as columns
    MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);

    double delta_t = 0.1; // time diff in sec

    /**
   * Student part begin
   */

    // predict sigma points
    for (int i = 0; i < Xsig_aug.cols(); i++)
    {
        VectorXd x_k = Xsig_aug.block<5, 1>(0, i);
        float v_k = x_k[2];
        float phi_k = x_k[3];
        float phi_k_dot = x_k[4];

        float nu_a = Xsig_aug(5, i);
        float nu_phi_dot_dot = Xsig_aug(6, i);

        VectorXd velocity_term(n_x);
        VectorXd acceleration_term(n_x);

        // avoid division by zero
        if (fabs(phi_k_dot) > 1e-3)
        {
            velocity_term[0] = (v_k / phi_k_dot) * (sin(phi_k + phi_k_dot * delta_t) - sin(phi_k));
            velocity_term[1] = (v_k / phi_k_dot) * (-cos(phi_k + phi_k_dot * delta_t) + cos(phi_k));
        }
        else
        {
            velocity_term[0] = v_k * cos(phi_k) * delta_t;
            velocity_term[1] = v_k * sin(phi_k) * delta_t;
        }

        velocity_term[2] = 0;
        velocity_term[3] = phi_k_dot * delta_t;
        velocity_term[4] = 0;

        acceleration_term[0] = 0.5 * pow(delta_t, 2) * cos(phi_k) * nu_a;
        acceleration_term[1] = 0.5 * pow(delta_t, 2) * sin(phi_k) * nu_a;
        acceleration_term[2] = delta_t * nu_a;
        acceleration_term[3] = 0.5 * pow(delta_t, 2) * nu_phi_dot_dot;
        acceleration_term[4] = delta_t * nu_phi_dot_dot;

        // write predicted sigma points into right column
        Xsig_pred.block<5, 1>(0, i) = x_k + (velocity_term + acceleration_term);
    }

    /**
   * Student part end
   */

    // print result
    std::cout
        << "Xsig_pred = " << std::endl
        << Xsig_pred << std::endl;

    // write result
    *Xsig_out = Xsig_pred;
}

/**
 * Programming assignment functions: 
 */

void UKF::PredictMeanAndCovariance(VectorXd *x_out, MatrixXd *P_out)
{

    // set state dimension
    int n_x = 5;

    // set augmented dimension
    int n_aug = 7;

    // define spreading parameter
    double lambda = 3 - n_aug;

    // create example matrix with predicted sigma points
    MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);
    Xsig_pred << 5.9374, 6.0640, 5.925, 5.9436, 5.9266, 5.9374, 5.9389, 5.9374, 5.8106, 5.9457, 5.9310, 5.9465, 5.9374, 5.9359, 5.93744,
        1.48, 1.4436, 1.660, 1.4934, 1.5036, 1.48, 1.4868, 1.48, 1.5271, 1.3104, 1.4787, 1.4674, 1.48, 1.4851, 1.486,
        2.204, 2.2841, 2.2455, 2.2958, 2.204, 2.204, 2.2395, 2.204, 2.1256, 2.1642, 2.1139, 2.204, 2.204, 2.1702, 2.2049,
        0.5367, 0.47338, 0.67809, 0.55455, 0.64364, 0.54337, 0.5367, 0.53851, 0.60017, 0.39546, 0.51900, 0.42991, 0.530188, 0.5367, 0.535048,
        0.352, 0.29997, 0.46212, 0.37633, 0.4841, 0.41872, 0.352, 0.38744, 0.40562, 0.24347, 0.32926, 0.2214, 0.28687, 0.352, 0.318159;

    // create vector for weights
    VectorXd weights = VectorXd(2 * n_aug + 1);

    // create vector for predicted state
    VectorXd x = VectorXd(n_x);

    // create covariance matrix for prediction
    MatrixXd P = MatrixXd(n_x, n_x);

    /**
   * Student part begin
   */

    for (int i = 0; i < Xsig_pred.cols(); i++)
    {
        // set weights
        if (i == 0)
            weights(i) = lambda / (lambda + n_aug);
        else
            weights(i) = 0.5 * 1 / (lambda + n_aug);

        // predict state mean
        x += weights(i) * Xsig_pred.col(i);
    }

    // predict state covariance matrix

    for (int i = 0; i < Xsig_pred.cols(); i++)
    {

        VectorXd delta_x = Xsig_pred.col(i) - x;

        // angle normalization
        while (delta_x(3) > M_PI)
            delta_x(3) -= 2. * M_PI;
        while (delta_x(3) < -M_PI)
            delta_x(3) += 2. * M_PI;

        P += weights(i) * delta_x * delta_x.transpose();
    }

    /**
   * Student part end
   */

    // print result
    std::cout << "Predicted state" << std::endl;
    std::cout << x << std::endl;
    std::cout << "Predicted covariance matrix" << std::endl;
    std::cout << P << std::endl;

    // write result
    *x_out = x;
    *P_out = P;
}

/**
 * Programming assignment functions: 
 */

void UKF::PredictRadarMeasurement(VectorXd *z_out, MatrixXd *S_out)
{

    // set state dimension
    int n_x = 5;

    // set augmented dimension
    int n_aug = 7;

    // set measurement dimension, radar can measure r, phi, and r_dot
    int n_z = 3;

    // define spreading parameter
    double lambda = 3 - n_aug;

    // set vector for weights
    VectorXd weights = VectorXd(2 * n_aug + 1);
    double weight_0 = lambda / (lambda + n_aug);
    double weight = 0.5 / (lambda + n_aug);
    weights(0) = weight_0;

    for (int i = 1; i < 2 * n_aug + 1; ++i)
    {
        weights(i) = weight;
    }

    // radar measurement noise standard deviation radius in m
    double std_radr = 0.3;

    // radar measurement noise standard deviation angle in rad
    double std_radphi = 0.0175;

    // radar measurement noise standard deviation radius change in m/s
    double std_radrd = 0.1;

    // create example matrix with predicted sigma points
    MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);
    Xsig_pred << 5.9374, 6.0640, 5.925, 5.9436, 5.9266, 5.9374, 5.9389, 5.9374, 5.8106, 5.9457, 5.9310, 5.9465, 5.9374, 5.9359, 5.93744,
        1.48, 1.4436, 1.660, 1.4934, 1.5036, 1.48, 1.4868, 1.48, 1.5271, 1.3104, 1.4787, 1.4674, 1.48, 1.4851, 1.486,
        2.204, 2.2841, 2.2455, 2.2958, 2.204, 2.204, 2.2395, 2.204, 2.1256, 2.1642, 2.1139, 2.204, 2.204, 2.1702, 2.2049,
        0.5367, 0.47338, 0.67809, 0.55455, 0.64364, 0.54337, 0.5367, 0.53851, 0.60017, 0.39546, 0.51900, 0.42991, 0.530188, 0.5367, 0.535048,
        0.352, 0.29997, 0.46212, 0.37633, 0.4841, 0.41872, 0.352, 0.38744, 0.40562, 0.24347, 0.32926, 0.2214, 0.28687, 0.352, 0.318159;

    // create matrix for sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug + 1);

    // mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);

    // measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z, n_z);

    /**
   * Student part begin
   */

    // transform sigma points into measurement space
    for (int i = 0; i < Xsig_pred.cols(); i++)
    {

        VectorXd x = Xsig_pred.col(i);
        float px = x[0];
        float py = x[1];
        float v = x[2];
        float phi = x[3];
        float phi_dot = x[4];

        Zsig(0, i) = sqrt(pow(px, 2) + pow(py, 2));
        Zsig(1, i) = atan2(py, px);
        Zsig(2, i) = ((px * cos(phi) * v) + (py * sin(phi) * v)) / sqrt(pow(px, 2) + pow(py, 2));
    }

    // calculate mean predicted measurement
    for (int i = 0; i < Zsig.cols(); i++)
    {

        z_pred += weights(i) * Zsig.col(i);
    }

    // calculate innovation covariance matrix S
    MatrixXd R = MatrixXd(n_z, n_z);
    R << pow(std_radr, 2), 0, 0,
        0, pow(std_radphi, 2), 0,
        0, 0, pow(std_radrd, 2);

    for (int i = 0; i < Zsig.cols(); i++)
    {
        VectorXd delta_z = Zsig.col(i) - z_pred;

        // angle normalization
        while (delta_z(1) > M_PI)
            delta_z(1) -= 2. * M_PI;
        while (delta_z(1) < -M_PI)
            delta_z(1) += 2. * M_PI;

        S += weights(i) * delta_z * delta_z.transpose();
    }

    S += R;

    /**
   * Student part end
   */

    // print result
    std::cout
        << "z_pred: " << std::endl
        << z_pred << std::endl;
    std::cout << "S: " << std::endl
              << S << std::endl;

    // write result
    *z_out = z_pred;
    *S_out = S;
}

/**
 * Programming assignment functions: 
 */

void UKF::UpdateState(VectorXd *x_out, MatrixXd *P_out)
{

    // set state dimension
    int n_x = 5;

    // set augmented dimension
    int n_aug = 7;

    // set measurement dimension, radar can measure r, phi, and r_dot
    int n_z = 3;

    // define spreading parameter
    double lambda = 3 - n_aug;

    // set vector for weights
    VectorXd weights = VectorXd(2 * n_aug + 1);
    double weight_0 = lambda / (lambda + n_aug);
    double weight = 0.5 / (lambda + n_aug);
    weights(0) = weight_0;

    for (int i = 1; i < 2 * n_aug + 1; ++i)
    {
        weights(i) = weight;
    }

    // create example matrix with predicted sigma points in state space
    MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);
    Xsig_pred << 5.9374, 6.0640, 5.925, 5.9436, 5.9266, 5.9374, 5.9389, 5.9374, 5.8106, 5.9457, 5.9310, 5.9465, 5.9374, 5.9359, 5.93744,
        1.48, 1.4436, 1.660, 1.4934, 1.5036, 1.48, 1.4868, 1.48, 1.5271, 1.3104, 1.4787, 1.4674, 1.48, 1.4851, 1.486,
        2.204, 2.2841, 2.2455, 2.2958, 2.204, 2.204, 2.2395, 2.204, 2.1256, 2.1642, 2.1139, 2.204, 2.204, 2.1702, 2.2049,
        0.5367, 0.47338, 0.67809, 0.55455, 0.64364, 0.54337, 0.5367, 0.53851, 0.60017, 0.39546, 0.51900, 0.42991, 0.530188, 0.5367, 0.535048,
        0.352, 0.29997, 0.46212, 0.37633, 0.4841, 0.41872, 0.352, 0.38744, 0.40562, 0.24347, 0.32926, 0.2214, 0.28687, 0.352, 0.318159;

    // create example vector for predicted state mean
    VectorXd x = VectorXd(n_x);
    x << 5.93637,
        1.49035,
        2.20528,
        0.536853,
        0.353577;

    // create example matrix for predicted state covariance
    MatrixXd P = MatrixXd(n_x, n_x);
    P << 0.0054342, -0.002405, 0.0034157, -0.0034819, -0.00299378,
        -0.002405, 0.01084, 0.001492, 0.0098018, 0.00791091,
        0.0034157, 0.001492, 0.0058012, 0.00077863, 0.000792973,
        -0.0034819, 0.0098018, 0.00077863, 0.011923, 0.0112491,
        -0.0029937, 0.0079109, 0.00079297, 0.011249, 0.0126972;

    // create example matrix with sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug + 1);
    Zsig << 6.1190, 6.2334, 6.1531, 6.1283, 6.1143, 6.1190, 6.1221, 6.1190, 6.0079, 6.0883, 6.1125, 6.1248, 6.1190, 6.1188, 6.12057,
        0.24428, 0.2337, 0.27316, 0.24616, 0.24846, 0.24428, 0.24530, 0.24428, 0.25700, 0.21692, 0.24433, 0.24193, 0.24428, 0.24515, 0.245239,
        2.1104, 2.2188, 2.0639, 2.187, 2.0341, 2.1061, 2.1450, 2.1092, 2.0016, 2.129, 2.0346, 2.1651, 2.1145, 2.0786, 2.11295;

    // create example vector for mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);
    z_pred << 6.12155,
        0.245993,
        2.10313;

    // create example matrix for predicted measurement covariance
    MatrixXd S = MatrixXd(n_z, n_z);
    S << 0.0946171, -0.000139448, 0.00407016,
        -0.000139448, 0.000617548, -0.000770652,
        0.00407016, -0.000770652, 0.0180917;

    // create example vector for incoming radar measurement
    VectorXd z = VectorXd(n_z);
    z << 5.9214, // rho in m
        0.2187,  // phi in rad
        2.0062;  // rho_dot in m/s

    // create matrix for cross correlation Tc
    MatrixXd Tc = MatrixXd(n_x, n_z);

    /**
   * Student part begin
   */

    // calculate cross correlation matrix

    for (int i = 0; i < Zsig.cols(); i++)
    {
        VectorXd delta_z = Zsig.col(i) - z_pred;

        // angle normalization
        while (delta_z(1) > M_PI)
            delta_z(1) -= 2. * M_PI;
        while (delta_z(1) < -M_PI)
            delta_z(1) += 2. * M_PI;

        VectorXd delta_x = Xsig_pred.col(i) - x;

        // angle normalization
        while (delta_x(3) > M_PI)
            delta_x(3) -= 2. * M_PI;
        while (delta_x(3) < -M_PI)
            delta_x(3) += 2. * M_PI;

        Tc += weights(i) * delta_x * delta_z.transpose();
    }

    // calculate Kalman gain K;

    MatrixXd K = Tc * S.inverse();

    // update state mean and covariance matrix

    VectorXd delta_z = z - z_pred;

    // angle normalization
    while (delta_z(1) > M_PI)
        delta_z(1) -= 2. * M_PI;
    while (delta_z(1) < -M_PI)
        delta_z(1) += 2. * M_PI;

    x += K * delta_z;

    P -= K * S * K.transpose();

    /**
   * Student part end
   */

    // print result
    std::cout
        << "Updated state x: " << std::endl
        << x << std::endl;
    std::cout << "Updated state covariance P: " << std::endl
              << P << std::endl;

    // write result
    *x_out = x;
    *P_out = P;
}