#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include "eigen-3.4.0/Eigen/Dense"
#include <cstdlib>

using namespace std;
using namespace Eigen;

// Function to calculate the integral of a vector
double Integral(const Eigen::VectorXd& x, const Eigen::VectorXd& y) {
    double integral = 0.0;
    for (int i = 1; i < x.size(); i++) {
        double dx = x(i) - x(i - 1);
        double avg_y = (y(i) + y(i - 1)) / 2.0;
        integral += dx * avg_y;
    }
    return integral;
}

int main() {
    // Data Acquisition and Re-arranging
    ifstream dataFile("naca2412.txt");
    vector<double> XB, YB;
    system("cls");
    double x, y;
    while (dataFile >> x >> y) {
        XB.push_back(x);
        YB.push_back(y);
    }
    reverse(YB.begin(), YB.end());  // Reverse the order of the elements in the vector
    reverse(XB.begin(), XB.end());  // to get the geometry in clockwise direction

    int nPan = XB.size() - 1;       // Number of panels

    // Calculation of control points and other geometric parameters
    double alpha_deg = 0,alpha = 0, pi = M_PI;
    cout << "Enter the value of alpha(deg): ";
    //cin >> alpha_deg;
    alpha_deg = 8;
    alpha = alpha_deg * (pi / 180);

    double Y[nPan], S[nPan], phi[nPan];
    VectorXd X(nPan), beta(nPan+1), V(nPan), CP(nPan);

    for (int i = 0; i < nPan ; i++) {
        X[i] = 0.5 * (XB[i] + XB[i + 1]);
        Y[i] = 0.5 * (YB[i] + YB[i + 1]);
        S[i] = sqrt(pow(XB[i + 1] - XB[i], 2) + pow(YB[i + 1] - YB[i], 2));
        phi[i] = atan2(YB[i + 1] - YB[i], XB[i + 1] - XB[i]);
        beta[i] = sin(phi[i] - alpha);
    }

    // Calculation of Coefficients
    MatrixXd C_n1(nPan, nPan ), C_n2(nPan, nPan), C_t1(nPan, nPan), C_t2(nPan, nPan);
    double A,B,C,D,E,F,G,P,Q;

    for (int i = 0; i < nPan; i++) {
        C_n1(i, i) = -1;
        C_n2(i, i) = 1;
        C_t1(i, i) = 0.5 * pi;
        C_t2(i, i) = 0.5 * pi;
        for (int j = 0; j < nPan; j++) {
             if (j != i) {
                A = - (X[i] - XB[j]) * cos(phi[j]) - (Y[i] - YB[j]) * sin(phi[j]);
                B = pow(X[i] - XB[j], 2) + pow(Y[i] - YB[j], 2);
                
                C = sin(phi[i] - phi[j]);
                D = cos(phi[i] - phi[j]);
                E = (X[i] - XB[j]) * sin(phi[j]) - (Y[i] - YB[j]) * cos(phi[j]);
                F = log(1 + (pow(S[j], 2) + (2 * A * S[j])) / B);
                G = atan2((E * S[j]), (B + A * S[j]));
                P = ((X[i] - XB[j]) * sin(phi[i] - 2 * phi[j])) + ((Y[i] - YB[j]) * cos(phi[i] - 2 * phi[j]));
                Q = ((X[i] - XB[j]) * cos(phi[i] - 2 * phi[j])) - ((Y[i] - YB[j]) * sin(phi[i] - 2 * phi[j]));

                C_n2(i, j) = D + ((0.5 * Q * F) / S[j]) - ((A * C + D * E) * (G / S[j]));
                C_n1(i, j) = 0.5 * D * F + C * G - C_n2(i, j);
                C_t2(i, j) = C + ((0.5 * P * F) / S[j]) + ((A * D - C * E) * (G / S[j]));
                C_t1(i, j) = 0.5 * C * F - D * G - C_t2(i, j);
            }
        }
    }
    // Computation of Influence Coefficients
    MatrixXd An(nPan+1, nPan+1), At(nPan, nPan+1);
    for (int i = 0; i < nPan; i++) {
        An(i, 0) = C_n1(i, 0);
        An(i, nPan ) = C_n2(0, nPan - 1);
        At(i, 0) = C_t1(i, 0);
        At(i, nPan ) = C_t2(i, nPan - 1);
        for (int j = 1; j < nPan; j++) {
            An(i, j) = C_n1(i, j) + C_n2(i, j - 1);
            At(i, j) = C_t1(i, j) + C_t2(i, j - 1);
        }
    }

    // Satisfying Kutta condition - Making the flow leave smoothly at the trailing edge
    An(nPan, 0) = 1;
    An(nPan, nPan) = 1;
    for (int j = 1; j < nPan; j++) {
        An(nPan, j) = 0;
    }
    beta[nPan] = 0;

    // Solve for Gamma and velocity/pressure
    VectorXd Gamma = An.fullPivLu().solve(beta);
    
    for (int i = 0; i < nPan; i++) {
        V[i] = cos(phi[i] - alpha);
        for (int j = 0; j < nPan+1; j++) {
            V[i] += At(i, j) * Gamma[j];
            CP[i] = 1 - pow(V[i], 2);
        }
    }
    
    ofstream vortex("VortexPanel.csv");
    vortex << "x/c,CP" << endl;

    for (int i = 0; i < nPan; i++) {
        vortex << X[i] << "," << CP[i] << endl;
    }
    vortex.close();

    ofstream Lift("Lift.csv");
    Lift << "alpha,Cl" << endl;

    // Calculation of Lift Coefficient
    for (int k=-30;k<=30;k+=5){ // Loop for different angles of attack from -30 to 30 degrees
        alpha = k * (pi / 180);

        for (int i = 0; i < nPan ; i++) {
        beta[i] = sin(phi[i] - alpha);
        }
        beta[nPan] = 0;

        // calculating for the different values of alpha
        VectorXd Gamma = An.fullPivLu().solve(beta);

        for (int i = 0; i < nPan; i++) {
        V[i] = cos(phi[i] - alpha);
        for (int j = 0; j < nPan+1; j++) {
            V[i] += At(i, j) * Gamma[j];
            CP[i] = 1 - pow(V[i], 2);
            }
        }

        // Calculation of Lift Coefficient
        // splitting the CP vector into upper and lower surface
        VectorXd CPl(CP.head((nPan) / 2)), CPu(CP.tail((nPan) / 2)); 
        VectorXd dx = X.segment((nPan) / 2, CPu.size());
        reverse(CPl.data(), CPl.data() + CPl.size()); // reversing the lower surface CP vector

        VectorXd dCP = CPl - CPu;

        // Calculate the integral of dCP/dx 
        double Cl = Integral(dx,dCP);

        Lift <<alpha*180/pi<<","<<Cl<< endl;
    }
    cout<<"Finished";
    
    return 0;
}
