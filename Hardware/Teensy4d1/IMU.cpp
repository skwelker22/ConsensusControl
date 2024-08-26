#include "IMU.h"

// Constructors ////////////////////////////////////////////////////////////////
IMU::IMU() {
  //initialize members
  omegaX = 0.0;
  omegaY = 0.0;
  omegaZ = 0.0;
  calOmegaX = 0.0;
  calOmegaY = 0.0;
  calOmegaZ = 0.0;
  accX = 0.0;
  accY = 0.0;
  accZ = 0.0;
  thetaX = 0.0;
  thetaY = 0.0;
  memset(kfOut, 0.0, N_STATES * sizeof(float));
  calibrationComplete = false;

  //initialize kalman filter
  kalmanFilter = new KalmanFilter1D();
}

//Class methods
void IMU::GetSensorDataAndFilter(TwoWire& wire) {
  //start I2C communication with the gyro
  wire.beginTransmission(0x68);

  //Switch on the low pass filter
  wire.write(0x1A);
  wire.write(0x05);
  wire.endTransmission();

  //set accel
  this->SetAccelSignals(wire);

  //set gyro
  this->SetGyroSignals(wire);

  //calculate the "measurements"
  thetaX = atan(accY / sqrt(accX * accX + accZ * accZ)) * RAD2DEG;
  thetaY = atan(accX / sqrt(accY * accY + accZ * accZ)) * RAD2DEG;

  //filter angles using accels to remove the errors
  float U[N_STATES] = { omegaX, omegaY };
  float Z[N_STATES] = { thetaX, thetaY };

  //only update the filter after the calibration routine is over
  //update the filter at a slower rate
  if (calibrationComplete) {
    kalmanFilter->Filter(kfOut, U, Z);
  }
}

void IMU::SetAccelSignals(TwoWire& wire) {
  //set sensitivity for accels
  wire.beginTransmission(0x68);
  wire.write(0x1C);
  wire.write(0x10);
  wire.endTransmission();

  //get data
  wire.beginTransmission(0x68);
  wire.write(0x3B);
  wire.endTransmission();
  wire.requestFrom(0x68, 6);

  //accels
  int16_t tmpAccXLSB = wire.read() << 8 | wire.read();
  int16_t tmpAccYLSB = wire.read() << 8 | wire.read();
  int16_t tmpAccZLSB = wire.read() << 8 | wire.read();

  //set accels
  accX = (float)tmpAccXLSB / 4096 - ACC_X_BIAS;
  accY = (float)tmpAccYLSB / 4096 - ACC_Y_BIAS;
  accZ = (float)tmpAccZLSB / 4096 - ACC_Z_BIAS;
}

void IMU::SetGyroSignals(TwoWire& wire) {
  //Set the sensistivity and scale factor
  wire.beginTransmission(0x68);
  wire.write(0x1B);
  wire.write(0x8);
  wire.endTransmission();

  //Access registers storing gyro measurements
  wire.beginTransmission(0x68);
  wire.write(0x43);
  wire.endTransmission();

  //request data
  wire.requestFrom(0x68, 6);
  int16_t tmpGyroX = wire.read() << 8 | wire.read();
  int16_t tmpGyroY = wire.read() << 8 | wire.read();
  int16_t tmpGyroZ = wire.read() << 8 | wire.read();

  //convert from bits to degrees per second
  omegaX = (float)tmpGyroX / 65.5;
  omegaY = (float)tmpGyroY / 65.5;
  omegaZ = (float)tmpGyroZ / 65.5;

  //remove bias
  if (calibrationComplete) {
    omegaX -= calOmegaX;
    omegaY -= calOmegaY;
    omegaZ -= calOmegaZ;
  }
}

void IMU::CalibrateGyro(TwoWire& wire)
{
  //calibration routine
  for (int iCal = 0; iCal < N_CAL; iCal++) {
    GetSensorDataAndFilter(wire);
    calOmegaX += omegaX;
    calOmegaY += omegaY;
    calOmegaZ += omegaZ;
    delay(1);
  }

  //average
  calOmegaX /= N_CAL;
  calOmegaY /= N_CAL;
  calOmegaZ /= N_CAL;

  calibrationComplete = true;
}

float IMU::GetFilterStates(int filterIx) {
  return kfOut[filterIx];
}