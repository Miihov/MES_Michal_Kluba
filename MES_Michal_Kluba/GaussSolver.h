#pragma once
#include <iostream>
#include <vector>

using namespace std;

class GaussSolver {
public:
	GaussSolver();

	void display(vector<vector<double>> A);
	vector<double> gauss(vector<vector<double>> A);

	static double* solve(double** tab, double* params, int params_c);
};