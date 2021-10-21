#include "SevenStateEkf.h"

SevenStateEkf::SevenStateEkf()
{
	for (int i = 0; i < EKF_NUM_STATES; i++)
	{
		for (int j = 0; j < EKF_NUM_STATES; j++) {
			if (i == j) {
				ekf_cov_[i][i] = kInitStdDevs[i] * kInitStdDevs[i];
				Q[i] = 0.F;
				if (i < 6) {
					R_GPS[i] = 0.F;
				}
			}
			else {
				ekf_cov_[i][j] = 0.F;
			}

		}
		ekf_state_[i] = kInitState[i];
	}

	// GPS measurement model covariance
	R_GPS[0] = kGpsPosXYStd * kGpsPosXYStd;
	R_GPS[1] = kGpsPosXYStd * kGpsPosXYStd;
	R_GPS[2] = kGpsPosZStd * kGpsPosZStd;
	R_GPS[3] = kGpsVelXYStd * kGpsVelXYStd;
	R_GPS[4] = kGpsVelXYStd * kGpsVelXYStd;
	R_GPS[5] = kGpsVelZStd * kGpsVelZStd;

	// magnetometer measurement model covariance
	R_Heading[0] = kHeadingYawStd * kHeadingYawStd;

	// load the transition model covariance
	Q[0] = kDtImu * kQPosXYStd*kQPosXYStd;
	Q[1] = kDtImu * kQPosXYStd*kQPosXYStd;
	Q[2] = kDtImu * kQPosZStd*kQPosZStd;
	Q[3] = kDtImu * kQVelXYStd*kQVelXYStd;
	Q[4] = kDtImu * kQVelXYStd*kQVelXYStd;
	Q[5] = kDtImu * kQVelZStd*kQVelZStd;
	Q[6] = kDtImu * kQYawStd*kQYawStd;
}

void SevenStateEkf::GetRbgPrimeAccelInertial(float attitude[], float rbg_prime[3][3], float accel[3], float accel_inertial[3])
{
	float cosphi = cos(attitude[0]);
	float sinphi = sin(attitude[0]);
	float costheta = cos(attitude[1]);
	float sintheta = sin(attitude[1]);
	float cospsi = cos(attitude[2]);
	float sinpsi = sin(attitude[2]);
	rbg_prime[0][0] = -costheta * sinpsi;
	rbg_prime[1][0] = costheta * cospsi;
	rbg_prime[2][0] = 0;
	rbg_prime[0][1] = -sinphi * sintheta*sinpsi - cosphi * cospsi;
	rbg_prime[1][1] = sinphi * sintheta*cospsi - cosphi * sinpsi;
	rbg_prime[2][1] = 0;
	rbg_prime[0][2] = -cosphi * sintheta*sinpsi + sinphi * cospsi;
	rbg_prime[1][2] = cosphi * sintheta*cospsi + sinphi * sinpsi;
	rbg_prime[2][2] = 0;
	accel_inertial[0] = accel[2] * (sinphi*sinpsi + cosphi * cospsi*sintheta) - accel[1] * (cosphi*sinpsi - cospsi * sinphi*sintheta) + accel[0] * cospsi*costheta;
	accel_inertial[1] = accel[1] * (cosphi*cospsi + sinphi * sinpsi*sintheta) - accel[2] * (cospsi*sinphi - cosphi * sinpsi*sintheta) + accel[0] * costheta*sinpsi;
	accel_inertial[2] = accel[2] * cosphi*costheta - accel[0] * sintheta + accel[1] * costheta*sinphi;
}
void SevenStateEkf::Predict(float attitude[], float accel[], float dt)
{
	float accel_inertial[3];
	float rbg_prime[3][3];
	GetRbgPrimeAccelInertial(attitude, rbg_prime, accel, accel_inertial);
	// Predict State
	ekf_state_[0] += ekf_state_[3] * dt;
	ekf_state_[1] += ekf_state_[4] * dt;
	ekf_state_[2] += ekf_state_[5] * dt;
	ekf_state_[3] += accel_inertial[0] * dt;
	ekf_state_[4] += accel_inertial[1] * dt;
	ekf_state_[5] += (accel_inertial[2] + GRAVITY) * dt;
	ekf_state_[6] = attitude[2];
	// Update Covariance matrix
	float g_prime[EKF_NUM_STATES][EKF_NUM_STATES];
	for (int i = 0; i < EKF_NUM_STATES; i++) {
		for (int j = 0; j < EKF_NUM_STATES; j++) {
			if (i == j) {
				g_prime[i][i] = 1.0F;
			}
			else {
				g_prime[i][j] = 0.F;
			}
		}
	}
	g_prime[0][3] = dt;
	g_prime[1][4] = dt;
	g_prime[2][5] = dt;
	g_prime[3][6] = (rbg_prime[0][0] * accel[0] + rbg_prime[0][1] * accel[1] + rbg_prime[0][2] * accel[2])*dt;
	g_prime[4][6] = (rbg_prime[1][0] * accel[0] + rbg_prime[1][1] * accel[1] + rbg_prime[1][2] * accel[2])*dt;
	g_prime[5][6] = (rbg_prime[2][0] * accel[0] + rbg_prime[2][1] * accel[1] + rbg_prime[2][2] * accel[2])*dt;

	//// 
	//ekfCov = g_prime * ekfCov* gPrime.transpose() + Q;
	// Update non-diagonal values first, as they are not used anywhere else
	ekf_cov_[3][0] = ekf_cov_[3][3] * g_prime[0][3];
	ekf_cov_[4][1] = ekf_cov_[4][4] * g_prime[1][4];
	ekf_cov_[5][2] = ekf_cov_[5][5] * g_prime[2][5];
	ekf_cov_[0][3] = ekf_cov_[3][3] * g_prime[0][3];
	ekf_cov_[4][3] = ekf_cov_[6][6] * g_prime[3][6] * g_prime[4][6];
	ekf_cov_[5][3] = ekf_cov_[6][6] * g_prime[3][6] * g_prime[5][6];
	ekf_cov_[6][3] = ekf_cov_[6][6] * g_prime[3][6];
	ekf_cov_[1][4] = ekf_cov_[4][4] * g_prime[1][4];
	ekf_cov_[3][4] = ekf_cov_[6][6] * g_prime[3][6] * g_prime[4][6];
	ekf_cov_[5][4] = ekf_cov_[6][6] * g_prime[4][6] * g_prime[5][6];
	ekf_cov_[6][4] = ekf_cov_[6][6] * g_prime[4][6];
	ekf_cov_[2][5] = ekf_cov_[5][5] * g_prime[2][5];
	ekf_cov_[3][5] = ekf_cov_[6][6] * g_prime[3][6] * g_prime[5][6];
	ekf_cov_[4][5] = ekf_cov_[6][6] * g_prime[4][6] * g_prime[5][6];
	ekf_cov_[6][5] = ekf_cov_[6][6] * g_prime[5][6];
	ekf_cov_[3][6] = ekf_cov_[6][6] * g_prime[3][6];
	ekf_cov_[4][6] = ekf_cov_[6][6] * g_prime[4][6];
	ekf_cov_[5][6] = ekf_cov_[6][6] * g_prime[5][6];
	// Diagonals
	ekf_cov_[0][0] = ekf_cov_[0][0] + ekf_cov_[3][3] * g_prime[0][3] * g_prime[0][3] + Q[0];
	ekf_cov_[1][1] = ekf_cov_[1][1] + ekf_cov_[4][4] * g_prime[1][4] * g_prime[1][4] + Q[1];
	ekf_cov_[2][2] = ekf_cov_[2][2] + ekf_cov_[5][5] * g_prime[2][5] * g_prime[2][5] + Q[2];
	ekf_cov_[3][3] = ekf_cov_[3][3] + ekf_cov_[6][6] * g_prime[3][6] * g_prime[3][6] + Q[3];
	ekf_cov_[4][4] = ekf_cov_[4][4] + ekf_cov_[6][6] * g_prime[4][6] * g_prime[4][6] + Q[4];
	ekf_cov_[5][5] = ekf_cov_[5][5] + ekf_cov_[6][6] * g_prime[5][6] * g_prime[5][6] + Q[5];
	ekf_cov_[6][6] = ekf_cov_[6][6] + Q[6];
}

void SevenStateEkf::UpdateFromGps(float position[], float velocity[])
{
	float delta_z[6];
	delta_z[0] = position[0] - ekf_state_[0];
	delta_z[1] = position[1] - ekf_state_[1];
	delta_z[2] = -position[2] - ekf_state_[2];
	delta_z[3] = velocity[0] - ekf_state_[3];
	delta_z[4] = velocity[1] - ekf_state_[4];
	delta_z[5] = velocity[2] - ekf_state_[5];



	M66 mti;
	mti.m00 = ekf_cov_[0][0] + R_GPS[0];
	mti.m11 = ekf_cov_[1][1] + R_GPS[1];
	mti.m22 = ekf_cov_[2][2] + R_GPS[2];
	mti.m33 = ekf_cov_[3][3] + R_GPS[3];
	mti.m44 = ekf_cov_[4][4] + R_GPS[4];
	mti.m55 = ekf_cov_[5][5] + R_GPS[5];
	mti.m03 = ekf_cov_[0][3];
	mti.m14 = ekf_cov_[1][4];
	mti.m25 = ekf_cov_[2][5];
	mti.m30 = ekf_cov_[3][0];
	mti.m34 = ekf_cov_[3][4];
	mti.m35 = ekf_cov_[3][5];
	mti.m41 = ekf_cov_[4][1];
	mti.m43 = ekf_cov_[4][3];
	mti.m45 = ekf_cov_[4][5];
	mti.m52 = ekf_cov_[5][2];
	mti.m53 = ekf_cov_[5][3];
	mti.m54 = ekf_cov_[5][4];
	M66 inv_mat = matrixSpec::invertMatrix(mti);


	float K[7][6] = {
	{ekf_cov_[0][0] * inv_mat.m00 + ekf_cov_[0][3] * inv_mat.m30, ekf_cov_[0][0] * inv_mat.m01 + ekf_cov_[0][3] * inv_mat.m31, ekf_cov_[0][0] * inv_mat.m02 + ekf_cov_[0][3] * inv_mat.m32, ekf_cov_[0][0] * inv_mat.m03 + ekf_cov_[0][3] * inv_mat.m33, ekf_cov_[0][0] * inv_mat.m04 + ekf_cov_[0][3] * inv_mat.m34, ekf_cov_[0][0] * inv_mat.m05 + ekf_cov_[0][3] * inv_mat.m35},
	{ekf_cov_[1][1] * inv_mat.m10 + ekf_cov_[1][4] * inv_mat.m40, ekf_cov_[1][1] * inv_mat.m11 + ekf_cov_[1][4] * inv_mat.m41, ekf_cov_[1][1] * inv_mat.m12 + ekf_cov_[1][4] * inv_mat.m42, ekf_cov_[1][1] * inv_mat.m13 + ekf_cov_[1][4] * inv_mat.m54, ekf_cov_[1][1] * inv_mat.m14 + ekf_cov_[1][4] * inv_mat.m44, ekf_cov_[1][1] * inv_mat.m15 + ekf_cov_[1][4] * inv_mat.m45},
	{ekf_cov_[2][2] * inv_mat.m20 + ekf_cov_[2][5] * inv_mat.m50, ekf_cov_[2][2] * inv_mat.m21 + ekf_cov_[2][5] * inv_mat.m51, ekf_cov_[2][2] * inv_mat.m22 + ekf_cov_[2][5] * inv_mat.m52, ekf_cov_[2][2] * inv_mat.m23 + ekf_cov_[2][5] * inv_mat.m53, ekf_cov_[2][2] * inv_mat.m24 + ekf_cov_[2][5] * inv_mat.m54, ekf_cov_[2][2] * inv_mat.m25 + ekf_cov_[2][5] * inv_mat.m55},
	{ekf_cov_[3][0] * inv_mat.m00 + ekf_cov_[3][3] * inv_mat.m30 + ekf_cov_[3][4] * inv_mat.m40 + ekf_cov_[3][5] * inv_mat.m50, ekf_cov_[3][0] * inv_mat.m01 + ekf_cov_[3][3] * inv_mat.m31 + ekf_cov_[3][4] * inv_mat.m41 + ekf_cov_[3][5] * inv_mat.m51, ekf_cov_[3][0] * inv_mat.m02 + ekf_cov_[3][3] * inv_mat.m32 + ekf_cov_[3][4] * inv_mat.m42 + ekf_cov_[3][5] * inv_mat.m52, ekf_cov_[3][0] * inv_mat.m03 + ekf_cov_[3][3] * inv_mat.m33 + ekf_cov_[3][4] * inv_mat.m54 + ekf_cov_[3][5] * inv_mat.m53, ekf_cov_[3][0] * inv_mat.m04 + ekf_cov_[3][3] * inv_mat.m34 + ekf_cov_[3][4] * inv_mat.m44 + ekf_cov_[3][5] * inv_mat.m54, ekf_cov_[3][0] * inv_mat.m05 + ekf_cov_[3][3] * inv_mat.m35 + ekf_cov_[3][4] * inv_mat.m45 + ekf_cov_[3][5] * inv_mat.m55},
	{ekf_cov_[4][1] * inv_mat.m10 + ekf_cov_[4][3] * inv_mat.m30 + ekf_cov_[4][4] * inv_mat.m40 + ekf_cov_[4][5] * inv_mat.m50, ekf_cov_[4][1] * inv_mat.m11 + ekf_cov_[4][3] * inv_mat.m31 + ekf_cov_[4][4] * inv_mat.m41 + ekf_cov_[4][5] * inv_mat.m51, ekf_cov_[4][1] * inv_mat.m12 + ekf_cov_[4][3] * inv_mat.m32 + ekf_cov_[4][4] * inv_mat.m42 + ekf_cov_[4][5] * inv_mat.m52, ekf_cov_[4][1] * inv_mat.m13 + ekf_cov_[4][3] * inv_mat.m33 + ekf_cov_[4][4] * inv_mat.m54 + ekf_cov_[4][5] * inv_mat.m53, ekf_cov_[4][1] * inv_mat.m14 + ekf_cov_[4][3] * inv_mat.m34 + ekf_cov_[4][4] * inv_mat.m44 + ekf_cov_[4][5] * inv_mat.m54, ekf_cov_[4][1] * inv_mat.m15 + ekf_cov_[4][3] * inv_mat.m35 + ekf_cov_[4][4] * inv_mat.m45 + ekf_cov_[4][5] * inv_mat.m55},
	{ekf_cov_[5][2] * inv_mat.m20 + ekf_cov_[5][3] * inv_mat.m30 + ekf_cov_[5][4] * inv_mat.m40 + ekf_cov_[5][5] * inv_mat.m50, ekf_cov_[5][2] * inv_mat.m21 + ekf_cov_[5][3] * inv_mat.m31 + ekf_cov_[5][4] * inv_mat.m41 + ekf_cov_[5][5] * inv_mat.m51, ekf_cov_[5][2] * inv_mat.m22 + ekf_cov_[5][3] * inv_mat.m32 + ekf_cov_[5][4] * inv_mat.m42 + ekf_cov_[5][5] * inv_mat.m52, ekf_cov_[5][2] * inv_mat.m23 + ekf_cov_[5][3] * inv_mat.m33 + ekf_cov_[5][4] * inv_mat.m54 + ekf_cov_[5][5] * inv_mat.m53, ekf_cov_[5][2] * inv_mat.m24 + ekf_cov_[5][3] * inv_mat.m34 + ekf_cov_[5][4] * inv_mat.m44 + ekf_cov_[5][5] * inv_mat.m54, ekf_cov_[5][2] * inv_mat.m25 + ekf_cov_[5][3] * inv_mat.m35 + ekf_cov_[5][4] * inv_mat.m45 + ekf_cov_[5][5] * inv_mat.m55},
	{ekf_cov_[6][3] * inv_mat.m30 + ekf_cov_[6][4] * inv_mat.m40 + ekf_cov_[6][5] * inv_mat.m50, ekf_cov_[6][3] * inv_mat.m31 + ekf_cov_[6][4] * inv_mat.m41 + ekf_cov_[6][5] * inv_mat.m51, ekf_cov_[6][3] * inv_mat.m32 + ekf_cov_[6][4] * inv_mat.m42 + ekf_cov_[6][5] * inv_mat.m52, ekf_cov_[6][3] * inv_mat.m33 + ekf_cov_[6][4] * inv_mat.m54 + ekf_cov_[6][5] * inv_mat.m53, ekf_cov_[6][3] * inv_mat.m34 + ekf_cov_[6][4] * inv_mat.m44 + ekf_cov_[6][5] * inv_mat.m54, ekf_cov_[6][3] * inv_mat.m35 + ekf_cov_[6][4] * inv_mat.m45 + ekf_cov_[6][5] * inv_mat.m55}
	};
	//for (int i = 0; i < 7; i++){
	//	printf("K: %.3f/%.3f/%.3f/%.3f/%.3f/%.3f\n", K[i][0], K[i][1], K[i][2], K[i][3], K[i][4], K[i][5]);
	//}
	ekf_state_[0] += K[0][0] * delta_z[0] + K[0][3] * delta_z[3] + K[0][4] * delta_z[4] + K[0][5] * delta_z[5];
	ekf_state_[1] += K[1][1] * delta_z[1] + K[1][3] * delta_z[3] + K[1][4] * delta_z[4] + K[1][5] * delta_z[5];
	ekf_state_[2] += K[2][2] * delta_z[2] + K[2][3] * delta_z[3] + K[2][4] * delta_z[4] + K[2][5] * delta_z[5];
	ekf_state_[3] += K[3][0] * delta_z[0] + K[3][1] * delta_z[1] + K[3][2] * delta_z[2] + K[3][3] * delta_z[3] + K[3][4] * delta_z[4] + K[3][5] * delta_z[5];
	ekf_state_[4] += K[4][0] * delta_z[0] + K[4][1] * delta_z[1] + K[4][2] * delta_z[2] + K[4][3] * delta_z[3] + K[4][4] * delta_z[4] + K[4][5] * delta_z[5];
	ekf_state_[5] += K[5][0] * delta_z[0] + K[5][1] * delta_z[1] + K[5][2] * delta_z[2] + K[5][3] * delta_z[3] + K[5][4] * delta_z[4] + K[5][5] * delta_z[5];
	ekf_state_[6] += K[6][0] * delta_z[0] + K[6][1] * delta_z[1] + K[6][2] * delta_z[2] + K[6][3] * delta_z[3] + K[6][4] * delta_z[4] + K[6][5] * delta_z[5];


	ekf_cov_[0][0] = -K[0][3] * ekf_cov_[3][0] + ekf_cov_[0][0] * (1.0F - K[0][0]);
	ekf_cov_[0][1] = -K[0][4] * ekf_cov_[4][1];
	ekf_cov_[0][2] = -K[0][5] * ekf_cov_[5][2];
	ekf_cov_[0][3] = -K[0][3] * ekf_cov_[3][3] - K[0][4] * ekf_cov_[4][3] - K[0][5] * ekf_cov_[5][3] + ekf_cov_[0][3] * (1.0F - K[0][0]);
	ekf_cov_[0][4] = -K[0][3] * ekf_cov_[3][4] - K[0][4] * ekf_cov_[4][4] - K[0][5] * ekf_cov_[5][4];
	ekf_cov_[0][5] = -K[0][3] * ekf_cov_[3][5] - K[0][4] * ekf_cov_[4][5] - K[0][5] * ekf_cov_[5][5];
	ekf_cov_[0][6] = -K[0][3] * ekf_cov_[3][6] - K[0][4] * ekf_cov_[4][6] - K[0][5] * ekf_cov_[5][6];

	ekf_cov_[1][0] = -K[1][3] * ekf_cov_[3][0];
	ekf_cov_[1][1] = -K[1][4] * ekf_cov_[4][1] + ekf_cov_[1][1] * (1.0F - K[1][1]);
	ekf_cov_[1][2] = -K[1][5] * ekf_cov_[5][2];
	ekf_cov_[1][3] = -K[1][3] * ekf_cov_[3][3] - K[1][4] * ekf_cov_[4][3] - K[1][5] * ekf_cov_[5][3];
	ekf_cov_[1][4] = -K[1][3] * ekf_cov_[3][4] - K[1][4] * ekf_cov_[4][4] - K[1][5] * ekf_cov_[5][4] + ekf_cov_[1][4] * (1.0F - K[1][1]);
	ekf_cov_[1][5] = -K[1][3] * ekf_cov_[3][5] - K[1][4] * ekf_cov_[4][5] - K[1][5] * ekf_cov_[5][5];
	ekf_cov_[1][6] = -K[1][3] * ekf_cov_[3][6] - K[1][4] * ekf_cov_[4][6] - K[1][5] * ekf_cov_[5][6];

	ekf_cov_[2][0] = -K[2][3] * ekf_cov_[3][0];
	ekf_cov_[2][1] = -K[2][4] * ekf_cov_[4][1];
	ekf_cov_[2][2] = -K[2][5] * ekf_cov_[5][2] + ekf_cov_[2][2] * (1.0F - K[2][2]);
	ekf_cov_[2][3] = -K[2][3] * ekf_cov_[3][3] - K[2][4] * ekf_cov_[4][3] - K[2][5] * ekf_cov_[5][3];
	ekf_cov_[2][4] = -K[2][3] * ekf_cov_[3][4] - K[2][4] * ekf_cov_[4][4] - K[2][5] * ekf_cov_[5][4];
	ekf_cov_[2][5] = -K[2][3] * ekf_cov_[3][5] - K[2][4] * ekf_cov_[4][5] - K[2][5] * ekf_cov_[5][5] + ekf_cov_[2][5] * (1.0F - K[2][2]);
	ekf_cov_[2][6] = -K[2][3] * ekf_cov_[3][6] - K[2][4] * ekf_cov_[4][6] - K[2][5] * ekf_cov_[5][6];


	ekf_cov_[3][0] = -K[3][0] * ekf_cov_[0][0] + ekf_cov_[3][0] * (1.0F - K[3][3]);
	ekf_cov_[3][1] = -K[3][1] * ekf_cov_[1][1] - K[3][4] * ekf_cov_[4][1];
	ekf_cov_[3][2] = -K[3][2] * ekf_cov_[2][2] - K[3][5] * ekf_cov_[5][2];
	ekf_cov_[3][3] = -K[3][0] * ekf_cov_[0][3] - K[3][4] * ekf_cov_[4][3] - K[3][5] * ekf_cov_[5][3] + ekf_cov_[3][3] * (1.0F - K[3][3]);
	ekf_cov_[3][4] = -K[3][1] * ekf_cov_[1][4] - K[3][4] * ekf_cov_[4][4] - K[3][5] * ekf_cov_[5][4] + ekf_cov_[3][4] * (1.0F - K[3][3]);
	ekf_cov_[3][5] = -K[3][2] * ekf_cov_[2][5] - K[3][4] * ekf_cov_[4][5] - K[3][5] * ekf_cov_[5][5] + ekf_cov_[3][5] * (1.0F - K[3][3]);
	ekf_cov_[3][6] = -K[3][4] * ekf_cov_[4][6] - K[3][5] * ekf_cov_[5][6] + ekf_cov_[3][6] * (1.0F - K[3][3]);

	ekf_cov_[4][0] = -K[4][0] * ekf_cov_[0][0] - K[4][3] * ekf_cov_[3][0];
	ekf_cov_[4][1] = -K[4][1] * ekf_cov_[1][1] + ekf_cov_[4][1] * (1.0F - K[4][4]);
	ekf_cov_[4][2] = -K[4][2] * ekf_cov_[2][2] - K[4][5] * ekf_cov_[5][2];
	ekf_cov_[4][3] = -K[4][0] * ekf_cov_[0][3] - K[4][3] * ekf_cov_[3][3] - K[4][5] * ekf_cov_[5][3] + ekf_cov_[4][3] * (1.0F - K[4][4]);
	ekf_cov_[4][4] = -K[4][1] * ekf_cov_[1][4] - K[4][3] * ekf_cov_[3][4] - K[4][5] * ekf_cov_[5][4] + ekf_cov_[4][4] * (1.0F - K[4][4]);
	ekf_cov_[4][5] = -K[4][2] * ekf_cov_[2][5] - K[4][3] * ekf_cov_[3][5] - K[4][5] * ekf_cov_[5][5] + ekf_cov_[4][5] * (1.0F - K[4][4]);
	ekf_cov_[4][6] = -K[4][3] * ekf_cov_[3][6] - K[4][5] * ekf_cov_[5][6] + ekf_cov_[4][6] * (1.0F - K[4][4]);

	ekf_cov_[5][0] = -K[5][0] * ekf_cov_[0][0] - K[5][3] * ekf_cov_[3][0];
	ekf_cov_[5][1] = -K[5][1] * ekf_cov_[1][1] - K[5][4] * ekf_cov_[4][1];
	ekf_cov_[5][2] = -K[5][2] * ekf_cov_[2][2] + ekf_cov_[5][2] * (1.0F - K[5][5]);
	ekf_cov_[5][3] = -K[5][0] * ekf_cov_[0][3] - K[5][3] * ekf_cov_[3][3] - K[5][4] * ekf_cov_[4][3] + ekf_cov_[5][3] * (1.0F - K[5][5]);
	ekf_cov_[5][4] = -K[5][1] * ekf_cov_[1][4] - K[5][3] * ekf_cov_[3][4] - K[5][4] * ekf_cov_[4][4] + ekf_cov_[5][4] * (1.0F - K[5][5]);
	ekf_cov_[5][5] = -K[5][2] * ekf_cov_[2][5] - K[5][3] * ekf_cov_[3][5] - K[5][4] * ekf_cov_[4][5] + ekf_cov_[5][5] * (1.0F - K[5][5]);
	ekf_cov_[5][6] = -K[5][3] * ekf_cov_[3][6] - K[5][4] * ekf_cov_[4][6] + ekf_cov_[5][6] * (1.0F - K[5][5]);

	ekf_cov_[6][0] = -K[6][0] * ekf_cov_[0][0] - K[6][3] * ekf_cov_[3][0];
	ekf_cov_[6][1] = -K[6][1] * ekf_cov_[1][1] - K[6][4] * ekf_cov_[4][1];
	ekf_cov_[6][2] = -K[6][2] * ekf_cov_[2][2] - K[6][5] * ekf_cov_[5][2];
	ekf_cov_[6][3] = ekf_cov_[6][3] - K[6][0] * ekf_cov_[0][3] - K[6][3] * ekf_cov_[3][3] - K[6][4] * ekf_cov_[4][3] - K[6][5] * ekf_cov_[5][3];
	ekf_cov_[6][4] = ekf_cov_[6][4] - K[6][1] * ekf_cov_[1][4] - K[6][3] * ekf_cov_[3][4] - K[6][4] * ekf_cov_[4][4] - K[6][5] * ekf_cov_[5][4];
	ekf_cov_[6][5] = ekf_cov_[6][5] - K[6][2] * ekf_cov_[2][5] - K[6][3] * ekf_cov_[3][5] - K[6][4] * ekf_cov_[4][5] - K[6][5] * ekf_cov_[5][5];
	ekf_cov_[6][6] = ekf_cov_[6][6] - K[6][3] * ekf_cov_[3][6] - K[6][4] * ekf_cov_[4][6] - K[6][5] * ekf_cov_[5][6];
}

float SevenStateEkf::UpdateFromHeading(float heading)
{
	float delta_z = heading - ekf_state_[6];
	if (delta_z > M_PI) delta_z -= 6.2831853071;
	if (delta_z < -M_PI) delta_z += 6.2831853071;
	//printf("state/meas/z_pred:%.2f/%.2f/%.2f\n", ekf_state_[6) * 180 / F_PI, mag_yaw * 180 / F_PI, zFromX(0) * 180 / F_PI);
	float inv_mat = 0;
	if ((ekf_cov_[6][6] + R_Heading[0]) != 0) {
		inv_mat = 1 / (ekf_cov_[6][6] + R_Heading[0]);
	}
	ekf_state_[3] += delta_z * ekf_cov_[3][6] * inv_mat;
	ekf_state_[4] += delta_z * ekf_cov_[4][6] * inv_mat;
	ekf_state_[5] += delta_z * ekf_cov_[5][6] * inv_mat;
	float delta_yaw = delta_z * ekf_cov_[6][6] * inv_mat;
	ekf_state_[6] += delta_yaw;

	ekf_cov_[0][0] = ekf_cov_[0][0];
	ekf_cov_[0][1] = 0;
	ekf_cov_[0][2] = 0;
	ekf_cov_[0][3] = ekf_cov_[0][3];
	ekf_cov_[0][4] = 0;
	ekf_cov_[0][5] = 0;
	ekf_cov_[0][6] = 0;



	ekf_cov_[1][0] = 0;
	ekf_cov_[1][1] = ekf_cov_[1][1];
	ekf_cov_[1][2] = 0;
	ekf_cov_[1][3] = 0;
	ekf_cov_[1][4] = ekf_cov_[1][4];
	ekf_cov_[1][5] = 0;
	ekf_cov_[1][6] = 0;


	ekf_cov_[2][0] = 0;
	ekf_cov_[2][1] = 0;
	ekf_cov_[2][2] = ekf_cov_[2][2];
	ekf_cov_[2][3] = 0;
	ekf_cov_[2][4] = 0;
	ekf_cov_[2][5] = ekf_cov_[2][5];
	ekf_cov_[2][6] = 0;

	ekf_cov_[3][0] = ekf_cov_[3][0];
	ekf_cov_[3][1] = 0;
	ekf_cov_[3][2] = 0;
	ekf_cov_[3][3] = ekf_cov_[3][3] - ekf_cov_[3][6] * ekf_cov_[6][3] * inv_mat;
	ekf_cov_[3][4] = ekf_cov_[3][4] - ekf_cov_[3][6] * ekf_cov_[6][4] * inv_mat;
	ekf_cov_[3][5] = ekf_cov_[3][5] - ekf_cov_[3][6] * ekf_cov_[6][5] * inv_mat;
	ekf_cov_[3][6] = ekf_cov_[3][6] - ekf_cov_[3][6] * ekf_cov_[6][6] * inv_mat;

	ekf_cov_[4][0] = 0;
	ekf_cov_[4][1] = ekf_cov_[4][1];
	ekf_cov_[4][2] = 0;
	ekf_cov_[4][3] = ekf_cov_[4][3] - ekf_cov_[6][3] * ekf_cov_[4][6] * inv_mat;
	ekf_cov_[4][4] = ekf_cov_[4][4] - ekf_cov_[4][6] * ekf_cov_[6][4] * inv_mat;
	ekf_cov_[4][5] = ekf_cov_[4][5] - ekf_cov_[4][6] * ekf_cov_[6][5] * inv_mat;
	ekf_cov_[4][6] = ekf_cov_[4][6] - ekf_cov_[4][6] * ekf_cov_[6][6] * inv_mat;

	ekf_cov_[5][0] = 0;
	ekf_cov_[5][1] = 0;
	ekf_cov_[5][2] = ekf_cov_[5][2];
	ekf_cov_[5][3] = ekf_cov_[5][3] - ekf_cov_[6][3] * ekf_cov_[5][6] * inv_mat;
	ekf_cov_[5][4] = ekf_cov_[5][4] - ekf_cov_[6][4] * ekf_cov_[5][6] * inv_mat;
	ekf_cov_[5][5] = ekf_cov_[5][5] - ekf_cov_[5][6] * ekf_cov_[6][5] * inv_mat;
	ekf_cov_[5][6] = ekf_cov_[5][6] - ekf_cov_[5][6] * ekf_cov_[6][6] * inv_mat;

	ekf_cov_[6][0] = 0;
	ekf_cov_[6][1] = 0;
	ekf_cov_[6][2] = 0;
	ekf_cov_[6][3] = ekf_cov_[6][3] * (1.0F - ekf_cov_[6][6] * inv_mat);
	ekf_cov_[6][4] = ekf_cov_[6][4] * (1.0F - ekf_cov_[6][6] * inv_mat);
	ekf_cov_[6][5] = ekf_cov_[6][5] * (1.0F - ekf_cov_[6][6] * inv_mat);
	ekf_cov_[6][6] = ekf_cov_[6][6] * (1.0F - ekf_cov_[6][6] * inv_mat);
	return delta_yaw;
}

void SevenStateEkf::GetPosition(float position[])
{
	position[0] = ekf_state_[0];
	position[1] = ekf_state_[1];
	position[2] = ekf_state_[2];
}

void SevenStateEkf::GetVelocity(float velocity[])
{
	velocity[0] = ekf_state_[3];
	velocity[1] = ekf_state_[4];
	velocity[2] = ekf_state_[5];
}

float SevenStateEkf::SetHeading(float heading)
{
	float delta_z = heading - ekf_state_[6];
	if (delta_z > M_PI) delta_z -= 6.2831853071;
	if (delta_z < -M_PI) delta_z += 6.2831853071;
	ekf_state_[6] = heading;
	return delta_z;
}

float SevenStateEkf::GetHeading()
{
	return ekf_state_[6];
}

void SevenStateEkf::SetAltitude(float altitude)
{
	ekf_state_[2] = -altitude;
}
