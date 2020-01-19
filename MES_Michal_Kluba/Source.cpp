#include "GaussSolver.h"
#include <iostream>

using namespace std;

enum Stan {
	LEWY_WARUNEK_BRZEGOWY = 1,
	BRAK_WARUNKU_BRZEGOWEGO = 0,
	PRAWY_WARUNEK_BRZEGOWY = 2
};

struct DaneGlobalne {
	int liczbaWezlow;
	int liczbaElementow;
	double dlugoscCalegoElementu;
	double dlugoscPojedynczegoElementu;
	double strumien;					//q
	double przekroj;					//s
	double k;							//stala k
	double alfa;
	double temperaturaOtoczenia;		//tx

	DaneGlobalne() {
		//Konstruktor domyslny
	}

	DaneGlobalne(int n, double l, double q, double s, double k, double a, double tx) {
		this->liczbaWezlow = n + 1;
		this->liczbaElementow = n;
		this->dlugoscCalegoElementu = l;
		this->dlugoscPojedynczegoElementu = l / n;
		this->strumien = q;
		this->przekroj = s;
		this->k = k;
		this->alfa = a;
		this->temperaturaOtoczenia = tx;
	}
};

struct Wezel {
	int id;					//id wezla
	double x;				//wspolrzedna x
	Stan stan;				//stan (warunek brzegowy)

	Wezel() {
		//Konstruktor domyslny
	}

	Wezel(int id, double x, Stan stan) {
		this->id = id;
		this->x = x;
		this->stan = stan;
	}
};

struct Element {
	double C;					//(s*k)/l
	Wezel wezly[2];				//kazdy element ma 2 wezly
	double H[2][2];				//macierz lokalna
	double P[2];				//wektor lokalny

	Element() {
		//Konstruktor domyslny
	}

	Element(Wezel w1, Wezel w2) {
		wezly[0] = w1;
		wezly[1] = w2;
	}

	void budujMacierz(double a, double s) {
		H[0][0] = C;
		H[0][1] = -C;
		H[1][1] = C;
		H[1][0] = -C;
		
		if (wezly[0].stan == PRAWY_WARUNEK_BRZEGOWY) H[0][0] += a * s;
		if (wezly[1].stan == PRAWY_WARUNEK_BRZEGOWY) H[1][1] += a * s;
	}

	void budujWektor(double q, double s, double a, double tx) {
		P[0] = 0;
		P[1] = 0;

		if (wezly[0].stan == LEWY_WARUNEK_BRZEGOWY) P[0] = q * s;
		if (wezly[1].stan == PRAWY_WARUNEK_BRZEGOWY) P[1] = -a * s * tx;
	}
};

struct Siatka {
	Wezel* wezly;				//wezly
	Element* elementy;			//elementy
	double** H;					//macierz globalna
	double* P;					//wektor globalny
	double* tempWezla;			//wynik temperatury na danym wezle;
	DaneGlobalne dane;			//dane

	//Budowanie siatki MES
	Siatka(DaneGlobalne dane) {
		this->dane = dane;
		Stan stan;

		wezly = new Wezel[dane.liczbaWezlow];
		for (int i = 0; i < dane.liczbaWezlow; i++) {
			stan = BRAK_WARUNKU_BRZEGOWEGO;
			if (i == 0) stan == LEWY_WARUNEK_BRZEGOWY;
			if (i == dane.liczbaWezlow - 1) stan = PRAWY_WARUNEK_BRZEGOWY;
			wezly[i] = Wezel(i, i, stan);
		}

		elementy = new Element[dane.liczbaElementow];
		for (int i = 0; i < dane.liczbaElementow; i++) {
			elementy[i] = Element(wezly[i], wezly[i + 1]);
		}
	}

	void obliczanieMacierzyLokalnej() {
		for (int i = 0; i < dane.liczbaElementow; i++) {
			elementy[i].C = (dane.przekroj*dane.k) / dane.dlugoscPojedynczegoElementu;
			elementy[i].budujMacierz(dane.alfa, dane.przekroj);
		}
	}

	void obliczanieWektoraLokalnego() {
		for (int i = 0; i < dane.liczbaElementow; i++) {
			elementy[i].budujWektor(dane.strumien, dane.przekroj, dane.alfa, dane.temperaturaOtoczenia);
		}
	}

	void budujMacierzGlobalna() {
		int wymiar = dane.liczbaWezlow;

		//Tworzenie globalnej macierzy kwadratowej o wielkosci wymiar x wymiar
		H = new double*[wymiar];
		for (int i = 0; i < wymiar; i++) {
			H[i] = new double[wymiar];
		}

		//Wypelnianie macierzy zerami
		for (int i = 0; i < wymiar; i++) {
			for (int j = 0; j < wymiar; j++) {
				H[i][j] = 0;
			}
		}

		//Wypelnianie macierzy odpowiednimi wartosciami
		for (int i = 0; i < dane.liczbaElementow; i++) {
			H[i][i] += elementy[i].H[0][0];
			H[i][i + 1] += elementy[i].H[0][1];
			H[i + 1][i] += elementy[i].H[1][0];
			H[i + 1][i + 1] += elementy[i].H[1][1];
		}
	}

	void budujWektorGlobalny() {
		int wymiar = dane.liczbaWezlow;

		//Tworzenie wektora globalnego i wypenianie go zerami
		P = new double[wymiar];
		for (int i = 0; i < wymiar; i++) {
			P[i] = 0;
		}

		//Wypelnianie wektora odpowiednimi wartosciami
		for (int i = 0; i < wymiar; i++) {
			P[i] += elementy[i].P[0];
			P[i + 1] += elementy[i].P[1];
		}
	}

	double** macierzGlobalna() {
		return H;
	}

	double* wektorGlobalny() {
		return P;
	}

	int liczbaWezlow() {
		return dane.liczbaWezlow;
	}
};

double* obliczMES(double n, double l, double q, double s, double k, double a, double tx) {
	cout << "------------------ Dane pobrane z pliku ------------------" << endl;
	cout << "K = " << k << endl;
	cout << "Alfa = " << a << endl;
	cout << "Przekroj (s) = " << s << endl;
	cout << "Strumien (q) = " << q << endl;
	cout << "Temperatura otoczenia = " << tx << endl;
	cout << "Ilosc elementow = " << n << endl;
	cout << "Ilosc wezlow = " << n + 1 << endl;
	cout << "Dlugosc calego elementu (l) = " << l << endl;
	cout << "Dlugosc pojedynczego elementu (l) = " << l/n << endl;

	cout << "-------------------- Rozwiazanie MES --------------------" << endl;
	DaneGlobalne daneGlobalne(n, l, q, s, k, a, tx);
	Siatka siatka(daneGlobalne);
	siatka.obliczanieMacierzyLokalnej();
	siatka.obliczanieWektoraLokalnego();
	siatka.budujMacierzGlobalna();
	siatka.budujWektorGlobalny();

	return GaussSolver::solve(siatka.macierzGlobalna(), siatka.wektorGlobalny(), siatka.liczbaWezlow());
}


int main() {
	cout << "Hello world";

	return 0;
}