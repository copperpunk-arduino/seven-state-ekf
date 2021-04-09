# SevenStateEkf
Position (3), Velocity (3), and Heading (1)

Inspired by the EKF for Udacity's Flying Car Nanodegree (https://github.com/udacity/FCND-Estimation-CPP). Modified to suit the purposes of an external IMU connected to a dual-GPS setup, providing differential-based heading.<br><br>
This is a pure C++ library, so we stuck with the C++ style guide, instead of converting it the Arduino library style guide. Hopefully that's not too confusing for you Arduino users.<br><br>
Key functions:
* `Predict(float attitude[], float accel[], float dt)` - Rotates acceleration into inertial axes and applies to solution; Takes heading directly
* `UpdateFromGps(float position[], float velocity[])` - Applies GPS position and velocity corrections to solution
* `UpdateFromHeading(float heading)` - Updates heading measurement based on external (non-IMU) reading

For an example of this library being used, please see the [IcmInsDualGps](https://github.com/copperpunk-arduino/icm-ins-dual-gps) repository.
