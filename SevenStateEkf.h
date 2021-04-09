#ifndef __SEVENSTATEEKF_H_
#define __SEVENSTATEEKF_H_

#define _USE_MATH_DEFINES
#include <math.h>
#include "MatrixSpec.h"

#define GRAVITY 9.80665F
#define F_PI
// Initial conditions
	const float kInitState[7] = { 0, 0, 0, 0, 0, 0, 0 };
	const float kInitStdDevs[7] = { .1, .1, .3, .1, .1, .3, .05 };

	//Process noise model
	//note that the process covariance matrix is diag(pow(QStd, 2))*dtIMU

	const float kQPosXYStd = .1;
	const float kQPosZStd = .05;
	const float kQVelXYStd = .25;
	const float kQVelZStd = .1;
	const float kQYawStd = .08;

	// GPS measurement std deviations
	const float kGpsPosXYStd = .715;
	const float kGpsPosZStd = 2.05;
	const float kGpsVelXYStd = .088;
	const float kGpsVelZStd = .31;

	// Heading Measurement
	const float kHeadingYawStd = 0.2;

	const float kDtImu = 0.05; // Expected IMU update interval

	const int EKF_NUM_STATES = 7;


class SevenStateEkf
{
public:
	SevenStateEkf();
	void GetRbgPrimeAccelInertial(float attitude[], float rbg_prime[3][3], float accel[3], float accel_inertial[3]);
	void Predict(float attitude[], float accel[], float dt);
	void UpdateFromGps(float position[], float velocity[]);
	float UpdateFromHeading(float heading);
	void GetPosition(float position[]);
	void GetVelocity(float velocity[]);
	float SetHeading(float heading);
	float GetHeading();
	void SetAltitude(float altitude);
private:
	float ekf_state_[EKF_NUM_STATES];
	float ekf_cov_[EKF_NUM_STATES][EKF_NUM_STATES];
	float Q[EKF_NUM_STATES];
	float R_GPS[6];
	float R_Heading[1];
};

#endif