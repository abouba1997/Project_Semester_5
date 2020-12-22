// Coursach_exemple.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <fstream>
#include <vector>

#define dvector vector<vector<double>>

using std::cin; using std::cout; using std::endl; using std::vector; using std::ofstream;

//Two dimensions wave
void setID(dvector& u, dvector& um, double c, double dt, double x, double y);
void wave(dvector& up, dvector& u, dvector& um, double c, double dt, double h);
void plotSolu(dvector& up, double t);

double f(double x, double y);

int main()
{
    int intervals = 20;  
    int n = intervals + 1;
    dvector up(n, vector<double>(n, 0)), u(n, vector<double>(n, 0)), um(n, vector<double>(n, 0));
    double C = 1;
    double tStop = 2;
    double t = 0.;

    //dx = dy = h
    double h = 1.0 / (n - 1.);
    double dt = h;
    double x = 1., y = 2.;

    setID(u, um, C, dt, x, y);
    wave(um, u, um, C, dt, h);
    plotSolu(um, 0);

    while (t < tStop) {
        t += dt;
        cout << "t = " << t << endl;
        wave(up, u, um, C, dt, h);
        um = u; u = up;
        plotSolu(u, t);
    }

    return 0;
}

// 2D

void setID(dvector& u, dvector& um, double c, double dt, double x, double y)
{
    int n = u.size();
    int m = u[0].size();
    vector<double> X, Y; double h = 1.0 / (n - 1.);

    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            u[i][j] = exp((-0.5) * dt*dt* (i - x)/2.0 - (-0.5) * dt * dt * (j - y) / 2.0);
        }
    }

    //setting um
    for (size_t i = 1; i < n - 1; i++) {
        for (size_t j = 1; j < n - 1; j++) {
            um[i][j] = u[i][j] + (1/2)*(((c * c) * (dt * dt)) / (h * h)) * (u[i + 1.][j] + u[i - 1.][j] + u[i][j + 1.] + u[i][j - 1.]);
        }
    }
    //update the boundary
    int a, b, nx = n;
    a = 0;
    for (size_t b = 1; b < nx - 1.; b++) {
        um[a][b] = u[a][b] + (((c * c) * (dt * dt)) / (h * h) / 2) * (u[a + 1.][b] + u[a + 1.][b] + u[a][b + 1.] + u[a][b - 1.]);
    }

    a = nx - 1;
    for (size_t b = 1; b < nx - 2; b++) {
        um[a][b] = u[a][b] + (((c * c) * (dt * dt)) / (h * h) / 2) * (u[a - 1.][b] + u[a - 1.][b] + u[a][b + 1.] + u[a][b - 1.]);
    }

    b = 0;
    for (size_t a = 1; a < nx - 1; a++) {
        um[a][b] = u[a][b] + (((c * c) * (dt * dt)) / (h * h) / 2) * (u[a + 1.][b] + u[a - 1.][b] + u[a][b + 1.] + u[a][b + 1.]);
    }

    b = nx - 2;
    for (size_t a = 1; a < nx - 2; a++) {
        um[a][b] = u[a][b] + (((c * c) * (dt * dt)) / (h * h) / 2) * (u[a - 1.][b] + u[a + 1.][b] + u[a][b + 1.] + u[a][b + 1.]);
    }

    a = 0; b = 0;
    um[a][b] = u[a][b] + (((c * c) * (dt * dt)) / (h * h) / 2) * (u[a + 1.][b] + u[a + 1.][b] + u[a][b + 1.] + u[a][b + 1.]);

    a = nx - 1; b = 0;
    um[a][b] = u[a][b] + (((c * c) * (dt * dt)) / (h * h) / 2) * (u[a - 1.][b] + u[a - 1.][b] + u[a][b + 1.] + u[a][b + 1.]);

    a = 0; b = nx - 1;
    um[a][b] = u[a][b] + (((c * c) * (dt * dt)) / (h * h) / 2) * (u[a + 1.][b] + u[a + 1.][b] + u[a][b - 1.] + u[a][b - 1.]);

    a = nx - 1; b = nx - 1;
    um[a][b] = u[a][b] + (((c * c) * (dt * dt)) / (h * h) / 2) * (u[a - 1.][b] + u[a - 1.][b] + u[a][b - 1.] + u[a][b - 1.]);
}

void wave(dvector& up, dvector& u, dvector& um, double c, double dt, double h)
{
    int nx = up.size();
    int a, b;

    for (size_t i = 1; i < nx - 1.; i++) {
        for (size_t j = 1; j < nx - 1.; j++) {
            up[i][j] = 2 * u[i][j] - um[i][j] + (((c * c) * (dt*dt))/(h*h)) * (u[i + 1.][j] + u[i - 1.][j] + u[i][j + 1.] + u[i][j - 1.]);
        }
    }

    //update the boundary
    a = 0;
    for (size_t b = 1; b < nx - 1.; b++) {
        up[a][b] = 2 * u[a][b] - um[a][b] + (((c * c) * (dt * dt)) / (h * h)) * (u[a + 1.][b] + u[a + 1.][b] + u[a][b + 1.] + u[a][b - 1.]);
    }

    a = nx - 1;
    for (size_t b = 1; b < nx - 2; b++) {
        up[a][b] = 2 * u[a][b] - um[a][b] + (((c * c) * (dt * dt)) / (h * h)) * (u[a - 1.][b] + u[a - 1.][b] + u[a][b + 1.] + u[a][b - 1.]);
    }

    b = 0;
    for (size_t a = 1; a < nx - 1; a++) {
        up[a][b] = 2 * u[a][b] - um[a][b] + (((c * c) * (dt * dt)) / (h * h)) * (u[a + 1.][b] + u[a - 1.][b] + u[a][b + 1.] + u[a][b + 1.]);
    }

    b = nx - 2;
    for (size_t a = 1; a < nx - 2; a++) {
        up[a][b] = 2 * u[a][b] - um[a][b] + (((c * c) * (dt * dt)) / (h * h)) * (u[a - 1.][b] + u[a + 1.][b] + u[a][b + 1.] + u[a][b + 1.]);
    }

    a = 0; b = 0;
    up[a][b] = 2 * u[a][b] - um[a][b] + (((c * c) * (dt * dt)) / (h * h)) * (u[a + 1.][b] + u[a + 1.][b] + u[a][b + 1.] + u[a][b + 1.]);

    a = nx - 1; b = 0;
    up[a][b] = 2 * u[a][b] - um[a][b] + (((c * c) * (dt * dt)) / (h * h)) * (u[a - 1.][b] + u[a - 1.][b] + u[a][b + 1.] + u[a][b + 1.]);

    a = 0; b = nx - 1;
    up[a][b] = 2 * u[a][b] - um[a][b] + (((c * c) * (dt * dt)) / (h * h)) * (u[a + 1.][b] + u[a + 1.][b] + u[a][b - 1.] + u[a][b - 1.]);

    a = nx - 1; b = nx - 1;
    up[a][b] = 2 * u[a][b] - um[a][b] + (((c * c) * (dt * dt)) / (h * h)) * (u[a - 1.][b] + u[a - 1.][b] + u[a][b - 1.] + u[a][b - 1.]);
}

void plotSolu(dvector &up, double t)
{
    int n = up.size();
    double h = 1.0 / (n - 1.);

    char fn[30];
    static int i = -1;

    i++; sprintf_s(fn, ".u.dat.%03d", i);
    ofstream outfile(fn);

    for (int i = 1; i < n; i++) {
        for (int j = 1; j < n; j++) {
            outfile << h * i << "\t\t\t" << h * j << "\t\t\t" << up[i][j] << endl;
        }
        outfile << endl;
    }
}

double f(double x, double y)
{
    double Lx = 3, Ly = 3;
    return exp(-0.5 * pow((x - Lx / 2.0), 2) - 0.5 * pow((y - Ly / 2.0), 2));
}
