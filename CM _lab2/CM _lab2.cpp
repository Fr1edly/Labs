#include <iostream>
#include <fstream>
#include <iomanip>
#include "slave.h"
double scalar(double* fobj, double* sobj, size_t len) {
    double ans = 0;
    for (size_t j = 0; j < len; j++) {
        ans += fobj[j] * sobj[j];
    }
    return ans;
}
using namespace std;

int main()
{
    double scalar(double*, double*, size_t);
    size_t len, type;
    ifstream fin("input.txt");
    double **slay;
    double *B_mass;
    
    #pragma region ввод

    if (!fin.is_open()) {
        cout << "Out of File!" << endl;
        return -1;
    }
    else {
        fin >> type >> len;
        slay = new double* [len];
        B_mass = new double[len];
        for (int i = 0; i < len; i++)
            slay[i] = new double[len]();

        for (int i = 0; i < len ; i++)
            for (int j = 0; j < len; j++) {
                fin >> slay[i][j];
                if (j == len - 1)  fin >> B_mass[i];
            }
    }
    fin.close();
#pragma endregion
    for (int i = 0; i < len; i++) {
        for (int j = 0; j < len; j++)
            cout << setw(5)<< slay[i][j];
        cout << setw(5) << B_mass[i] << endl;
    }

    ofstream fout("output.txt");
    double* ans = slave::orto(slay, B_mass, len);
    double* s_ans = slave::iter(slay, B_mass, len, 0.001);
    double* wtf = slave::miter(slay, B_mass, len);
    cout << "-------------------------------------------" << endl;
    for (int i = 0; i < len; i++)
        cout << "\t" << s_ans[i] ;
      cout <<endl << "-------------------------------------------" << endl;
      for (int i = 0; i < len; i++)
          cout << "\t" << wtf[i];
      cout << endl << "-------------------------------------------" << endl;
    fout << "Вектор решения: (";
    for (int i = 0; i < len; i++) {
        fout << setw(3) << ans[i];
    }
    fout << ")";
    double** rev_ans = new double* [len];
    cout << endl;
    for (int i = 0; i < len; i++) {
        for (int j = 0; j < len; j++) {
            B_mass[j] = slay[j][i];
        }
        rev_ans[i] = slave::orto(slay, B_mass, len);
        for (int j = 0; j < len; j++)
            cout << "\t" << rev_ans[i][j];
        cout << endl;
    }

    for (int i=0; i < len; i++)
        cout << setw(5) << ans[i];
}
