#include "GaussSolver.h"
#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

GaussSolver::GaussSolver(){
}

void GaussSolver::display(vector<vector<double>> A) {
	int n = A.size();
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n + 1; j++) {
			cout << A[i][j] << "\t";
			if (j == n - 1) {
				cout << "| ";
			}
		}
		cout << "\n";
	}
	cout << endl;
}

vector<double> GaussSolver::gauss(vector< vector<double> > A) {
	int n = A.size();

	for (int i = 0; i < n; i++) {
		// Search for maximum in this column
		double maxEl = abs(A[i][i]);
		int maxRow = i;
		for (int k = i + 1; k < n; k++) {
			if (abs(A[k][i]) > maxEl) {
				maxEl = abs(A[k][i]);
				maxRow = k;
			}
		}

		// Swap maximum row with current row (column by column)
		for (int k = i; k < n + 1; k++) {
			double tmp = A[maxRow][k];
			A[maxRow][k] = A[i][k];
			A[i][k] = tmp;
		}

		// Make all rows below this one 0 in current column
		for (int k = i + 1; k < n; k++) {
			double c = -A[k][i] / A[i][i];
			for (int j = i; j < n + 1; j++) {
				if (i == j) {
					A[k][j] = 0;
				}
				else {
					A[k][j] += c * A[i][j];
				}
			}
		}
	}

	// Solve equation Ax=b for an upper triangular matrix A
	vector<double> x(n);
	for (int i = n - 1; i >= 0; i--) {
		x[i] = A[i][n] / A[i][i];
		for (int k = i - 1; k >= 0; k--) {
			A[k][n] -= A[k][i] * x[i];
		}
	}
	return x;
}



double * GaussSolver::solve(double ** tab, double * params, int params_c) {

	GaussSolver gauss = GaussSolver();

	int n = params_c;

	vector<double> line(n + 1, 0);
	vector< vector<double> > A(n, line);

	// Read input data
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			A[i][j] = tab[i][j];
		}
	}

	for (int i = 0; i < n; i++) {
		A[i][n] = params[i];
	}

	// Print input
	gauss.display(A);

	// Calculate solution
	vector<double> x(n);
	x = gauss.gauss(A);

	double * res = new double[n];
	for (int i = 0; i < n; i++) {
		res[i] = x[i];
	}

	return res;
}
