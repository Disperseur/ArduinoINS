/*
Filtre de kalman (sensor fusion) pour détermination d'attitude

Pistes d'ameliorations:
- Mieux estimer Q et R: probablement que la meilleure methode serait de calculer les matrices Qr et Rr (pour les angles d'euler fournis par l'accelero) 
puis de faire la transformation en matrice 4x4 ...

- utiliser des incertitudes differentes pour les differents angles d'euler (roll moins precis que le pitch)

- Ajouter le magnetometre pour le yaw

*/


#define RAD2DEG 180.0f / PI
#define DELTA_ACCEL_TRESHOLD 0.15
#include <BasicLinearAlgebra.h>
#include "Arduino_BMI270_BMM150.h"
using namespace BLA;
using namespace std;

using Mat = Matrix<4, 4, double>;
using Vec = Matrix<4, 1, double>;

double deltaT = 10.0; //ms

const Matrix<4, 4, double> I4 = {1.0, 0, 0, 0,
                                 0, 1.0, 0, 0,
                                 0, 0, 1.0, 0,
                                 0, 0, 0, 1.0};

Matrix<4, 1, double> EulerToQuaternion(const Matrix<3, 1, double> &EulerAngles);
Matrix<3, 1, double> QuaternionToEuler(const Matrix<4, 1, double> &q);


//Kalman filter vectors and matricess
Matrix<3, 1, double> EulerAngles = {0, 0, 0};
Vec Xestimate   = {1, 0, 0, 0};
Mat Pestimate   = I4;
Vec Xpredicted;
Mat Ppredicted;
Mat K;
Mat A;
Mat B;
Mat Q           = 0.001*I4;
Mat H           = I4;
Mat R           = 1.55e-10*I4;



int get_pitch_roll_from_accelerometer(float &pitch, float &roll);



void setup() {
  Serial.begin(460800);
  while (!Serial);
  Serial.println("Started");

  if (!IMU.begin()) {
    Serial.println("Failed to initialize IMU!");
    while (1);
  }

}

void loop() {
  //prediction state and error covariance
  //  compute the new state matrix with gyro measurements
  float w1, w2, w3; //gyro measurements on x, y and z axis
  IMU.readGyroscope(w1, w2, w3);

  B = {   0, -w1/2, -w2/2, -w3/2,
       w1/2,     0,  w3/2, -w2/2,
       w2/2, -w3/2,     0,  w1/2,
       w3/2,  w2/2, -w1/2,     0};

  A = I4 + deltaT*B;
  
  Xpredicted = A*Xestimate;
  Ppredicted = A*Pestimate*(~A) + Q; //~ stands for transposition

  //kalman gain
  K = Ppredicted*(~H)* Inverse(H*Ppredicted*(~H) + R);

  //estimate
  // get measurement (from accelero and magneto for sensor fusion)
  float the, phi, psi;
  psi = EulerAngles(2);
  get_pitch_roll_from_accelerometer(the, phi);
  Matrix<3, 1, double> euler = {phi, the, psi};
  Vec z = EulerToQuaternion(euler);

  Xestimate = Xpredicted + K*(z - H*Xpredicted); //this is the output
  EulerAngles = QuaternionToEuler(Xestimate);

  //error covariance
  Pestimate = Ppredicted - K*H*Ppredicted;

  Matrix<3, 1, double> EulerAnglesDeg = RAD2DEG * EulerAngles;
  Serial.print(EulerAnglesDeg(0)); //roll
  Serial.print('\t');
  Serial.println(EulerAnglesDeg(1)); //pitch

  delay(deltaT);
}



Matrix<4, 1, double> EulerToQuaternion(const Matrix<3, 1, double> &EulerAngles) {
    Matrix<4, 1, double> q;

    double phi = EulerAngles(0); // roll
    double the = EulerAngles(1); // pitch
    double psi = EulerAngles(2); // yaw

    double c1 = cos(phi * 0.5);
    double c2 = cos(the * 0.5);
    double c3 = cos(psi * 0.5);

    double s1 = sin(phi * 0.5);
    double s2 = sin(the * 0.5);
    double s3 = sin(psi * 0.5);

    q(0) =  c1*c2*c3 + s1*s2*s3; // scalar part
    q(1) =  s1*c2*c3 - c1*s2*s3;
    q(2) =  c1*s2*c3 + s1*c2*s3;
    q(3) =  c1*c2*s3 - s1*s2*c3;

    return q;
}


Matrix<3, 1, double> QuaternionToEuler(const Matrix<4, 1, double> &q) {
    Matrix<3, 1, double> e;

    double w = q(0);
    double x = q(1);
    double y = q(2);
    double z = q(3);

    // --- Roll φ ---
    double sinr = 2.0 * (w * x + y * z);
    double cosr = 1.0 - 2.0 * (x * x + y * y);
    double phi  = atan2(sinr, cosr);

    // --- Pitch θ ---
    double sinp = 2.0 * (w * y - z * x);
    // Clamp to handle numerical noise
    if (abs(sinp) >= 1.0)
        sinp = (sinp > 0 ? 1.0 : -1.0);

    double the = asin(sinp);

    // --- Yaw ψ ---
    double siny = 2.0 * (w * z + x * y);
    double cosy = 1.0 - 2.0 * (y * y + z * z);
    double psi  = atan2(siny, cosy);

    e(0) = phi; // roll
    e(1) = the; // pitch
    e(2) = psi; // yaw

    return e;
}



//matrices and vectors functions
int get_pitch_roll_from_accelerometer(float &pitch, float &roll) {
  float xyz_norm, delta_acceleration;
  float x, y, z;
  float gravity_projection;
  bool is_moving = true;

  if (IMU.accelerationAvailable()) {
    IMU.readAcceleration(x, y, z);

    // //Debug
    // Serial.print(x);
    // Serial.print('\t');
    // Serial.print(y);
    // Serial.print('\t');
    // Serial.print(z);
    // Serial.print('\t');

    //if there is no other acceleration than g (norm(x, y, z) = 1) then we can compute attitude solution
    xyz_norm            = sqrt(x*x + y*y + z*z);
    delta_acceleration  = abs(1 - xyz_norm);
    is_moving           = delta_acceleration > DELTA_ACCEL_TRESHOLD;

    if(!is_moving) {
      //pitch computation
      if(abs(x) < 1) pitch = asin(x);
      else pitch = x/abs(x) * PI/2;

      //roll computation
      gravity_projection = -y/cos(pitch);
      if(abs(gravity_projection) < 1) roll = asin(gravity_projection);
      else roll = gravity_projection/abs(gravity_projection) * PI/2;

      return(0);
    }
    else return(-1);
  }
}

