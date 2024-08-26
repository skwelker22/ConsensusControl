#include <Wire.h>
#include "IMU.h"

IMU* imu = new IMU();
uint32_t LoopTimer,TmTimer;
const int THETA_X_IX=0, THETA_Y_IX=1;
bool printHeader=true;
uint32_t LOG_CNT=5e5;

//function helper prototypes
void PrintTM(void);

void setup() {
  Serial.begin(57600);
  pinMode(13,OUTPUT);
  digitalWrite(13,HIGH);

  //set clock speed of the I2C
  Wire.setClock(400000);
  Wire.begin();
  delay(250);

  Wire.beginTransmission(0x68);
  Wire.write(0x6B);
  Wire.write(0x00);
  //terminate connection with gyro
  Wire.endTransmission();

  //run calibration scheme
  imu->CalibrateGyro(Wire);

  //initialize loop timer
  LoopTimer=micros();
  TmTimer=micros();
}

void loop() {
  //capture signals from gyro and filter to get angles
  imu->GetSensorDataAndFilter(Wire);

  //print tm
  PrintTM();

  //update timer
  while(micros()-LoopTimer<4000);
  LoopTimer=micros();

}

void PrintTM(void)
{
  //update serial every second
  if(printHeader)
  {
    float GYROBIAS[3]={0.0};
    float ACCELBIAS[3]={0.0};
    imu->GetGyroBias(GYROBIAS);
    imu->GetAccelBias(ACCELBIAS);

    //print out bias only the first time
    Serial.print("Gyro_Bias_X:");
    Serial.print(GYROBIAS[0]);
    Serial.print(",");
    Serial.print(" Gyro_Bias_Y:");
    Serial.print(GYROBIAS[1]);
    Serial.print(",");
    Serial.print(" Gyro_Bias_Z:");
    Serial.print(GYROBIAS[2]);
    Serial.print("\n");
    Serial.print("Accel_Bias_X:");
    Serial.print(ACCELBIAS[0]);
    Serial.print(",");
    Serial.print(" Accel_Bias_Y:");
    Serial.print(ACCELBIAS[1]);
    Serial.print(",");
    Serial.print(" Accel_Bias_Z:");
    Serial.print(ACCELBIAS[2]);
    Serial.print("\n");
    printHeader=false;
  }
  if(micros()-TmTimer>=LOG_CNT)
  {

    float OMEGAS[3]={0.0};
    float ACCELS[3]={0.0};
    float THETAS[2]={0.0};
    imu->GetOmegas(OMEGAS);
    imu->GetAccels(ACCELS);
    //get states from filter
    THETAS[0]=imu->GetFilterStates(THETA_X_IX);
    THETAS[1]=imu->GetFilterStates(THETA_Y_IX);

    float printTime=micros()*1e-6;
    Serial.print("t=");
    Serial.print(printTime);
    Serial.print("\n");
    //raw gyro
    Serial.print("omegaX[deg/sec]:");
    Serial.print(OMEGAS[0]);
    Serial.print(",");
    Serial.print(" omegaY[deg/sec]:");
    Serial.print(OMEGAS[1]);
    Serial.print(",");
    Serial.print(" omegaZ[deg/sec]:");
    Serial.print(OMEGAS[2]);
    Serial.print("\n");

    //accelerometer
    Serial.print("Acc_X[g]:");
    Serial.print(ACCELS[0]);
    Serial.print(",");
    Serial.print(" Acc_Y[g]:");
    Serial.print(ACCELS[1]);
    Serial.print(",");
    Serial.print(" Acc_Z[g]:");
    Serial.print(ACCELS[2]);
    Serial.print("\n");
    //Filtered Angles
    Serial.print("RollAngle_[deg]:");
    Serial.print(THETAS[0]);
    Serial.print(",");
    Serial.print(" PitchAngle_[deg]:");
    Serial.print(THETAS[1]);
    Serial.print("\n");
    Serial.println("");
    TmTimer=micros();
  }
}
