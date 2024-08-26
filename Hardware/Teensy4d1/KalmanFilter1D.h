const int N_STATES = 2;
const float UNC_INIT = 2.0 * 2.0;
const float MEAS_VAR = 3.0 * 3.0;
const float DT = 4e-3;

class KalmanFilter1D {
private:
  float X[N_STATES];  //state vector
  float P[N_STATES][N_STATES];

public:
  KalmanFilter1D() {
    for (int iState = 0; iState < N_STATES; ++iState) {
      //zero out states and state covariance
      X[iState] = 0.0;
      P[iState][iState] = UNC_INIT;
    }
  };

  void Filter(float *filterOut, float *U, float *Z) {
    for (int iState = 0; iState < N_STATES; ++iState) {
      //update states and state covariance with current input
      X[iState] += DT * U[iState];
      P[iState][iState] += (DT * DT) * (4 * 4);

      //calculate kalman gain
      float tmpVar = P[iState][iState];
      float K = tmpVar * 1 / (tmpVar + MEAS_VAR);

      //perform kalman update
      X[iState] += K * (Z[iState] - X[iState]);
      P[iState][iState] *= (1 - K);

      //set output states
      filterOut[iState] = X[iState];
    }
  }
};
