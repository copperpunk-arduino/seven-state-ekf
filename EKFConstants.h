//const float kInitState[] = { 0, 0, -1, 0, 0, 0, 0 };
//const float kInitStdDevs[] = { .1, .1, .3, .1, .1, .3, .05 };
//
////Process noise model
////note that the process covariance matrix is diag(pow(QStd, 2))*dtIMU
//
//const float kQPosXYStd = .1;
//const float kQPosZStd = .05;
//const float kQVelXYStd = .25;
//const float kQVelZStd = .1;
//const float kQYawStd = .08;
//
//// GPS measurement std deviations
//const float kGpsPosXYStd = .715;
//const float kGpsPosZStd = 2.05;
//const float kGpsVelXYStd = .088;
//const float kGpsVelZStd = .31;
//
//// Magnetometer
//const float kMagYawStd = 0.2;
//
//const float dtIMU = 0.005;
//
//const int EKF_NUM_STATES = 7;
//const float TWO_PI = 2.0*M_PI;