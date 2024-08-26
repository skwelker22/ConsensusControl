#include <Wire.h>
#include "KalmanFilter1D.h"

const float ACC_X_BIAS=-0.06;
const float ACC_Y_BIAS=-0.01;
const float ACC_Z_BIAS=0.09;
const float RAD2DEG=(180.0/3.14159);
const int N_CAL=2e3;

class IMU
{
  private:
	  float omegaX;
	  float omegaY;
	  float omegaZ;
	  float calOmegaX;
	  float calOmegaY;
	  float calOmegaZ;
	  //int iCal;
	  float accX;
	  float accY;
	  float accZ;
	  float thetaX;
	  float thetaY;
	  bool calibrationComplete;
	
	  //kalman filter
	  float kfOut[N_STATES];
	  KalmanFilter1D* kalmanFilter;
	
  public:
    IMU();
	void GetSensorDataAndFilter(TwoWire& wire);
  	void SetGyroSignals(TwoWire& wire);
	void SetAccelSignals(TwoWire& wire);
	void CalibrateGyro(TwoWire& wire);
	float GetFilterStates(int filterIx);

	//getters for TM
	void GetOmegas(float* OMEGA){
		OMEGA[0]=omegaX;
		OMEGA[1]=omegaY;
		OMEGA[2]=omegaZ;
	}

	void GetAccels(float* ACCEL){
		ACCEL[0]=accX;
		ACCEL[1]=accY;
		ACCEL[2]=accZ;
	}

	void GetThetas(float* THETAS){
		THETAS[0]=thetaX;
		THETAS[1]=thetaY;
	}

	void GetAccelBias(float* ACCELBIAS){
		ACCELBIAS[0]=ACC_X_BIAS;
		ACCELBIAS[1]=ACC_Y_BIAS;
		ACCELBIAS[2]=ACC_Z_BIAS;
	}

	void GetGyroBias(float* GYROBIAS){
		GYROBIAS[0]=calOmegaX;
		GYROBIAS[1]=calOmegaY;
		GYROBIAS[2]=calOmegaZ;
	}
 
};