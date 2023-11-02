#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include "eigen-3.4.0/Eigen/Dense"
#include <cstdlib>
using namespace std;
using namespace Eigen;
int main(){
    //Variable initialisation
    double V_inf, pi = M_PI, alpha_deg,alpha_rad;

    system("cls");
    //system("clear");
    int nPan;
    cout << "Source panel for a unit cylinder.\nEnter the number of panels: ";
    // cin >> nPan;
    nPan=8;
    cout << "Enter the Magnitude of Velocity (m/s): ";
    // cin >> V_inf;
    V_inf = 1;
    cout << "Enter the angle of attack (degrees): ";
    // cin >> alpha_deg;
    alpha_deg = 0;
    alpha_rad = (pi/180)*alpha_deg;

    // Define arrays and variables
    MatrixXd I(nPan, nPan), J(nPan, nPan); // dense arrays from Eigen 
    VectorXd V_n(nPan);  // which allow for matrix operations
    int Max = 500; // number of points to plot the ideal Cp curve
    double XB[nPan + 1], YB[nPan + 1], x[nPan], y[nPan];
    double A, B, C_n, D_n, C_t, D_t, E, F, G, S[nPan];
    double phi[nPan], beta[nPan], Cp[nPan], Cp_ideal[Max], V_t[nPan];

    // creating the geometry of the cylinder
    // Initialize arrays XB and YB with nPan+1 points on the surface of a cylinder
    for (int i = 0; i <= nPan; i++){
        // Calculate the x and y coordinates of the boundary points
        double angle = pi + pi/nPan - 2*i*pi/nPan; // Calculate the angle in radians
        XB[i] = cos(angle); // x coordinate 
        YB[i] = sin(angle); // y coordinate 
        
    }
    // Calculate the midpoint, length, and orientation of each panel
    for (int i = 0; i < nPan; i++) {
        // Calculate the midpoint of the i-th panel (control point)
        x[i] = 0.5 * (XB[i + 1] + XB[i]);
        y[i] = 0.5 * (YB[i + 1] + YB[i]);

        // Calculate the orientation of the i-th panel
        phi[i] = atan2(YB[i + 1] - YB[i], XB[i + 1] - XB[i]);

        // Calculate the length of the i-th panel
        S[i] = sqrt(pow(XB[i + 1] - XB[i], 2) + pow(YB[i + 1] - YB[i], 2));
         
        // // Calculate the angle between the panel and the freestream velocity
        beta[i] = phi[i] + pi/2 - alpha_rad;

    }
    cout<<endl;
    //Calculating the Integrals, I(for normal component) and J(for tangential component)
    for (int i = 0; i < nPan; i++){ // iterate over each panel of the cylinder
        V_n[i] = -2*pi*V_inf*cos(beta[i]);
        I(i, i) = pi; // set the diagonal of the I matrix to pi
        J(i, i) = 0;
        for (int j = 0; j < nPan; j++){ // iterate over each panel of the cylinder again
            if (j != i){ // skip the panel if it is the same as the current panel
                // calculate intermediate variables for the source panel method
                double dx = x[i] - XB[j];
                double dy = y[i] - YB[j];

                A = -(dx)*cos(phi[j]) - (dy)*sin(phi[j]);
                B = pow(dx,2) + pow(dy,2);
                
                C_n = sin(phi[i] - phi[j]);
                C_t = -cos(phi[i] - phi[j]);

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
    cout << "Values of Lambda:\n" << Lambda << endl; 

    // double LambdaSum = Lambda.sum();
    // cout << "Sum of Lambda: " << LambdaSum << endl;  //Sum of Lambda should be zero
    
    ofstream final("SourcePanel.csv"); 
    cout << "Panel\tTheta(deg)\tCp\n";
    //Evaluating total surface velocity at panel and pressure coefficient
    for (int i = 0; i < nPan; i++){ 
            double sum = 0, angle; 
            V_t[i] = V_inf*sin(beta[i]); // calculate the total surface velocity at the control point
            for (int j = 0; j < nPan; j++){ 
                V_t[i] += (Lambda[j]/(2*pi))*J(i, j); 
                }

            Cp[i] = 1 - pow((V_t[i]/V_inf),2);

            angle = 180 - beta[i]*180/pi;   
            if (angle < 0) {
                angle += 360;
            }
            cout << i+1 << "\t"<< angle << "\t\t" << Cp[i] << "\n"; 
            final << angle*pi/180 << "," << Cp[i] << "\n";
    }
    //final << 2*pi << "," << Cp[0] << "\n"; // close the loop
    final.close();

    // Calculating the ideal Cp values (analytical solution)
    ofstream IdealCp("Cp_ideal_values.csv");
    for (int i = 0; i < Max; ++i) {
        double theta = 2.0 * pi * i / Max;
        Cp_ideal[i] = 1 - 4*pow(sin(theta),2);
        IdealCp << theta << "," << Cp_ideal[i] << endl;
    }
    IdealCp.close();
    /* Plotting functionality included with the help gnuplot(needs to be installed)*/
    return 0;
}
