//===========================================================================
//SOFTWARE:     3D geometric model for CNT and GNP networks
//OBJECTIVE:    Definition of Matrix and matrix operations
//AUTHOR:       Angel Mora
//E-MAIL:       angel.mora@cicy.mx
//===========================================================================

#ifndef MATHMATRIX_H
#define MATHMATRIX_H

#include<cmath>
#include<iomanip>
#include<iostream>
#include<vector>
#include<cstdlib>
using namespace std;

const double Zero = 1.0E-10;

//---------------------------------------------------------------------------
class MathMatrix
{
public:
	//-----------------------------------------------
	//Data Member
	vector<vector<double> > element;

	//-----------------------------------------------
	//Constructors
	MathMatrix() {};
	//Constructor for identity matrix
	MathMatrix(const int& M);
	//Constructor for MxN zero matrix
	MathMatrix(const int& M, const int& N);
	MathMatrix(const vector<double>& vec1D);
	MathMatrix(const vector<vector<double> >& vec2D);
	MathMatrix(const double *vec, const int& N);
	MathMatrix(const double *vec, const int& M, const int& N);

	//-----------------------------------------------
	//Member Functions
	void Evalu(const double *vec, const int& N);	
	void Evalu(const double *vec, const int& M, const int& N);
	void operator=(const vector<double>& vec);
	void operator=(const vector<vector<double> >& vec);	
	MathMatrix operator*(const MathMatrix& matrix);
	MathMatrix operator*(const double& R);
	//Addition operator for two matrices
	MathMatrix operator+(const MathMatrix& matrix);
	MathMatrix operator-(const MathMatrix& matrix);
	MathMatrix Inverse();
	MathMatrix Transpose();
	MathMatrix GetCol(int cn);
	int Rows();
	int Cols();
	int Symmetry();	
    double determinant();

	//overload the operator of output stream
	friend ostream& operator<<(ostream& o, const MathMatrix& matrix);
};
//overload the operator of output stream
inline ostream& operator<<(ostream& o,const MathMatrix& matrix)
{
	for (int i=0; i<int(matrix.element.size()); i++)
	{
		for(int j=0; j<int(matrix.element[i].size()); j++)
		{
			if (fabs(matrix.element[i][j])<1.0E-8)
			{
				o <<setw(8) << setfill(' ') << 0 << " ";
			}
			else
			{
				o <<setw(8) << setfill(' ') << matrix.element[i][j] << " ";
			}
		}
		o << endl;
	}
	return o;
}
//---------------------------------------------------------------------------
#endif
//===========================================================================
