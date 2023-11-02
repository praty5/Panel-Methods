#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include "eigen-3.4.0/Eigen/Dense"
#include <cstdlib>
using namespace std;
using namespace Eigen;

int main(){
    system("cls");

    //Variable initialisation
    double V_inf, pi = M_PI, alpha_deg,alpha_rad;
    // Read airfoil data from file
    ifstream dataFile("naca0012.txt");
    vector<double> XB, YB;
    system("cls");
    double a, b;
    while (dataFile >> a >> b) {
        XB.push_back(a);
        YB.push_back(b);
    }
    dataFile.close(); 
    
    reverse(YB.begin(), YB.end());  // Reverse the order of the elements in the vector
    reverse(XB.begin(), XB.end());  // to get the geometry in clockwise direction
    
    int nPan = XB.size() - 1; // number of panels

    alpha_deg = 0;  // angle of attack in degrees
    alpha_rad = (pi/180)*alpha_deg;

    // Define arrays and variables
    MatrixXd I(nPan, nPan), J(nPan, nPan); // dense arrays from Eigen 
    VectorXd V_n(nPan),x(nPan), y(nPan), Cp(nPan);  // which allow for matrix operations
    double A, B, C_n, D_n, C_t, D_t, E, F, G, S[nPan];
    double phi[nPan], beta[nPan], V_t[nPan];

    
    // Calculate the midpoint, length, and orientation of each panel
    for (int i = 0; i < nPan; i++) {
        // Calculate the midpoint of the i-th panel
        x[i] = 0.5 * (XB[i + 1] + XB[i]);
        y[i] = 0.5 * (YB[i + 1] + YB[i]);
        // Calculate the length of the i-th panel
        S[i] = sqrt(pow(XB[i + 1] - XB[i], 2) + pow(YB[i + 1] - YB[i], 2));
        // Calculate the orientation of the i-th panel
        phi[i] = atan2(YB[i + 1] - YB[i], XB[i + 1] - XB[i]);
        // Calculate the angle between the panel and the freestream velocity
        beta[i] = phi[i] + pi/2 - alpha_rad;
    }

    //Calculating the Integrals, I(for normal component) and J(for tangential component)
    for (int i = 0; i < nPan; i++){ // iterate over each panel of the cylinder
            V_n[i] = -2*pi*V_inf*cos(beta[i]);
            I(i, i) = pi; // set the diagonal of the I matrix to pi
            J(i, i) = 0;  // set the diagonal of the J matrix to 0
            for (int j = 0; j < nPan; j++){ // iterate over each panel of the cylinder again
                if (j != i){ // skip the panel if it is the same as the current panel
                    // calculate intermediate variables for the source panel method
                    double dx = x[i] - XB[j];
                    double dy = y[i] - YB[j];

                    A = -(dx)*cos(phi[j]) - (dy)*sin(phi[j]);
                    B = pow(dx,2) + pow(dy,2);
                    
                    C_n = sin(phi[i] - phi[j]), C_t = -cos(phi[i] - phi[j]);

                    D_n = -(dx)*sin(phi[i]) + (dy)*cos(phi[i]);
                    D_t = (dx)*cos(phi[i]) + (dy)*sin(phi[i]);

                    E = sqrt(B - pow(A,2));
                    G = (atan2((S[j] + A),E) - atan2(A,E));
                    F = log(1+(pow(S[j],2) + 2*A*S[j])/B);
                    
                    // calculate the integral for the source panel method
                    I(i, j) = C_n*F*0.5 + ((D_n - A*C_n)/E)*G;  
                    J(i, j) = C_t*F*0.5 + ((D_t - A*C_t)/E)*G;
                }
            }
        }
    //Solving for lambda/source panel strengths using LU decomposition
    VectorXd Lambda = I.fullPivLu().solve(V_n); 

    // double LambdaSum = Lambda.sum();
    // cout << "Sum of Lambda: " << LambdaSum << endl;   //Sum of Lambda should be zero   
    
    //Evaluating total surface velocity at panel and pressure coefficient
    for (int i = 0; i < nPan; i++){ 
            V_t[i] = V_inf*sin(beta[i]); // calculate the total surface velocity at the midpoint of the panel
            for (int j = 0; j < nPan; j++){ 
                V_t[i] += (Lambda[j]/(2*pi))*J(i, j); 
                }
            Cp[i] = 1 - pow((V_t[i]/V_inf),2);
    }

    // saving Cp data to file
    ofstream lower("lower.csv");
    ofstream upper("upper.csv");
    int middle = floor(nPan/2); // Find the midpoint of the panels
    for (int i = 0; i < middle; i++) {
        lower << x[i] << "," << Cp[i] << "\n";      //Cp of lower surface of airfoil
        upper << x[i+middle] << "," << Cp[i+middle] << "\n"; //Cp of upper surface of airfoil
    }
    lower.close();
    upper.close();
    cout<< "Finished";
   
    /* Plotting functionality included with the help gnuplot(needs to be installed)*/
    return 0;
}
