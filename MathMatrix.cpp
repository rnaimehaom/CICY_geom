//===========================================================================
//SOFTWARE:     3D geometric model for CNT and GNP networks
//OBJECTIVE:    Definition of Matrix and matrix operations
//AUTHOR:       Angel Mora
//E-MAIL:       angel.mora@cicy.mx
//===========================================================================

#include "MathMatrix.h"

//---------------------------------------------------------------------------
//Constructor

//Constructor for an identity matrix of dimensions MxM
MathMatrix::MathMatrix(const int& M)
{
	vector<double> temp(M, 0.0);
	element.assign(M, temp);

	for (int i = 0; i < M; i++)
	{
		element[i][i] = 1.0;
	}
}
//Constructor for a matrix of dimensions MxN initialized with all components zero
MathMatrix::MathMatrix(const int& M, const int& N)
{
	vector<double> temp(N,0.0);
	element.assign(M,temp);
}
//Constructor for a row matrix from a vector
MathMatrix::MathMatrix(const vector<double>& vec1D)
{
	int N = int(vec1D.size());
	vector<double> temp(N);
	element.assign(1,temp);
	for (int i=0; i<N; i++)
	{
		element[0][i] = vec1D[i];
	}
}
//Constructor from a 2D vector
MathMatrix::MathMatrix(const vector<vector<double> >& vec2D)
{
	element = vec2D;
}
//Constructor for a matrix from a pointer with N elements
MathMatrix::MathMatrix(const double *vec, const int& N)
{
	if (N<=0)
	{
		cout << "Error. Vector to initialize matrix has size less than 1: size N=" << N << endl;
		exit (0);
	}
	vector<double> temp(N);
	element.assign(1,temp);
	for (int i=0; i<N; i++)
	{
		element[0][i] = vec[i];
	}
}
//Constructor for a matrix from a pointer with MxN elements
MathMatrix::MathMatrix(const double *vec, const int& M, const int& N)
{
	if (M<=0 || N<=0)
	{
		cout << "Error. One of dimensions of a 2D vector used to initialize a matrix is less than 1. M=" << M << "N=" << N << endl;
		exit (0);
	}
	vector<double> temp(N);
	element.assign(M,temp);
	if (N!=1)
	{
		for (int i=0; i<M; i++)
		{
			for (int j=0; j<N; j++)
			{
				element[i][j] = vec[M*i+j];
			}
		}
	}
	else 
	{
		for (int i=0; i<M; i++)
		{
			for (int j=0; j<N; j++)
			{
				element[i][j] = vec[i];
			}
		}
	}
}
//---------------------------------------------------------------------------
//Set matrix elements from a pointer with N elements, previous values and size of matrix are discarded
void MathMatrix::Evalu(const double *vec, const int& N)
{
	if (N<=0)
	{
		cout << "Error. Vector to set matrix element has size less than 1: size N=" << N << endl;
		exit (0);
	}
	vector<double> temp(N);
	element.assign(1,temp);
	for (int i=0; i<N; i++)
	{
		element[0][i] = vec[i];
	}
}
//---------------------------------------------------------------------------
//Set matrix elements from a pointer with MxN elements, previous values and size of matrix are discarded
void MathMatrix::Evalu(const double *vec, const int& M, const int& N)
{
	if (M<=0 || N<=0)
	{
		cout << "Error. One of dimensions of a 2D vector used to set a matrix is less than 1. M=" << M << "N=" << N << endl;
		exit (0);
	}
	vector<double> temp(N);
	element.assign(M,temp);
	if (N!=1)
	{
		for (int i=0; i<M; i++)
		{
			for (int j=0; j<N; j++)
			{
				element[i][j] = vec[M*i+j];
			}
		}
	}
	else 
	{
		for (int i=0; i<M; i++)
		{
			for (int j=0; j<N; j++)
			{
				element[i][j] = vec[i];
			}
		}
	}
}
//---------------------------------------------------------------------------
//Set matrix elements from a 1D vector, previous values and size of matrix are discarded
void MathMatrix::operator=(const vector<double>& vec)
{
	element.assign(1,vec);
}
//Set matrix elements from a 2D vector, previous values and size of matrix are discarded
void MathMatrix::operator=(const vector<vector<double> >& vec)
{
	element = vec;
}
//---------------------------------------------------------------------------
//Multiplication operator when multiplying by another matrix
MathMatrix MathMatrix::operator*(const MathMatrix& matrix)
{
	if (element[0].size() != matrix.element.size())
	{
		cout << "Error, it fails to implement multiplication between two matrices, the number of cloumn is not equal to the number of row!" << endl;
		exit(0);
	}

	//This object has size MxL
	//matrix has size LxN
	//Thus, resulting matrix will have size MxN

	//Calculate dimension L
	int L = int(matrix.element.size());
	
	//Size of matrix taht results from multiplication: M rows and N columns
	int M = int(element.size());
	int N = int(matrix.element[0].size());

	//Matrix to store result of multiplication
	MathMatrix Tem_mat(M,N);

	for(int i=0; i<M; i++)
	{
		for(int j=0; j<N; j++)
		{
			//Variable to accumulate the result of element ij
			double acc = 0.0;
			for(int k=0; k<L; k++)
			{
				acc = acc + element[i][k] * matrix.element[k][j];
			}
			Tem_mat.element[i][j] = acc;
		}
	}

	return Tem_mat;
}
//---------------------------------------------------------------------------
//Multiplication operator when multiplying by a double (real number)
MathMatrix MathMatrix::operator*(const double& R)
{
	//Get matrix dimensions: M rows and N columns
	int M = int(element.size());
	int N = int(element[0].size());
	//Matrix to store result of multiplication
	MathMatrix Tem_mat(M,N);

	for(int i=0; i<M; i++)
	{
		for(int j=0; j<N; j++)
		{
			Tem_mat.element[i][j] = element[i][j] * R;
		}
	}

	return Tem_mat;
}
//---------------------------------------------------------------------------
//Addition operator for two matrices
MathMatrix MathMatrix::operator+(const MathMatrix& matrix)
{
	if ((element.size() != matrix.element.size())||(element[0].size() != matrix.element[0].size()))
	{
		cout << "Error, it fails to implement addition between two matrices, the number of cloumn is not equal to the number of row!" << endl;
		exit(0);
	}
	//Get matrix dimensions: M rows and N columns
	int M = int(element.size());
	int N = int(element[0].size());
	//Matrix to store result of addition
	MathMatrix Tem_mat(M,N);

	for(int i=0; i<M; i++)
	{
		for(int j=0; j<N; j++)
		{
			Tem_mat.element[i][j] = element[i][j] + matrix.element[i][j];
		}
	}

	return Tem_mat;
}
//---------------------------------------------------------------------------
//Matrix subtraction
MathMatrix MathMatrix::operator-(const MathMatrix& matrix)
{
	if ((element.size() != matrix.element.size())||(element[0].size() != matrix.element[0].size()))
	{
		cout << "Error. Subtraction between two matrices is not possible as they have different dimensions." << endl;
		cout << "Matrix1 dimensions:" << element.size() << "x" << element[0].size() << endl;
		cout << "Matrix2 dimensions:" << matrix.element.size() << "x" << matrix.element[0].size() << endl;
		exit(0);
	}

	//Get matrix dimensions: M rows and N columns
	int M = int(element.size());
	int N = int(element[0].size());

	//MAtrix to store result
	MathMatrix Tem_mat(M,N);

	for(int i=0; i<M; i++)
	{
		for(int j=0; j<N; j++)
		{
			Tem_mat.element[i][j] = element[i][j] - matrix.element[i][j];
		}
	}

	return Tem_mat;
}
//---------------------------------------------------------------------------
//matrix inversion
MathMatrix MathMatrix::Inverse()
{
	vector<vector<double> > vec = element;
	//--------------------------------------------------
	//to judge if this matrix is a positive definite matrix

	//--------------------------------------------------
	//performs matrix inversion
	int stRank = int(vec.size());
	vector<double> b(stRank, 0.0);
	for (int k = 0; k < stRank; k++)
	{
		double w = vec[0][0];
		int m = stRank - k - 1;
		for (int i = 1; i < stRank; i++)
		{
			double g = vec[i][0];
			b[i] = g / w;
			if (i <= m)
			{
				b[i] = -b[i];
			}
			for (int j = 1; j <= i; j++)
			{
				vec[i - 1][j - 1] = vec[i][j] + g * b[j];
			}
		}
		vec[stRank - 1][stRank - 1] = 1.0 / w;
		for (int i = 1; i < stRank; i++)
		{
			vec[stRank - 1][i - 1] = b[i];
		}
	}
	for (int i = 0; i < stRank - 1; i++)
	{
		for (int j = i + 1; j < stRank; j++)
		{
			vec[i][j] = vec[j][i];
		}
	}

	MathMatrix matrix = vec;
	return matrix;
}
//---------------------------------------------------------------------------
//matrix transposition
MathMatrix MathMatrix::Transpose()
{
	int M = int(element.size());
	int N = int(element[0].size());

	vector<vector<double> > vec;
	vector<double> temp(M);
	vec.assign(N, temp);

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < M; j++)
		{
			vec[i][j] = element[j][i];
		}
	}

	MathMatrix matrix = vec;
	return matrix;
}
//---------------------------------------------------------------------------
//To covert a column of matrix into a new matrix
MathMatrix MathMatrix::GetCol(int cn) {

	//Create a column matrix to store the result
	MathMatrix matrix(int(element.size()), 1);

	//Copy column cn
	for (int i = 0; i<int(element.size()); i++)
		matrix.element[i][0] = element[i][cn];

	return matrix;
}         
//---------------------------------------------------------------------------
//This function returns the number of rows in a matrix
int MathMatrix::Rows()
{
	return int(element.size());
}
//---------------------------------------------------------------------------
//This function returns the number of columns in a matrix
int MathMatrix::Cols()
{
	return int(element[0].size());
}
//---------------------------------------------------------------------------
//This function determines if a matrix is symmetric 
//-1: is not square matrix
// 0: is a square matrix but nonsymmetric
// 1: is a square matrix but symmetric
int MathMatrix::Symmetry()
{
	int M = int(element.size());
	int N = int(element[0].size());

	if (M == N)
	{
		for (int i = 0; i < M; i++)
		{
			for (int j = i + 1; j < M; j++)
			{
				if (fabs(element[i][j] - element[j][i]) > Zero)
				{
					return 0;
				}
			}
		}
	}
	else
	{
		return -1;
	}
	return 1;
}
//---------------------------------------------------------------------------
//It is assumed
double MathMatrix::determinant()
{
    if (element.size() != element.front().size()) {
        return NAN;
    }
    //Base case
    if (element.size() == 2) {
        //Determinant of a 2x2 matrix
        return element[0][0]*element[1][1] - element[1][0]*element[0][1];
    } else {
        //Initial size of squared matrix
        int n = (int)element.size();
        //Create a reduced (i.e. smaller) matrix for the recursion
        MathMatrix reduced(n-1,n-1);
        //Initialize value of determinant
        double det=0;
        //Variable to change the signs
        double sign = 1;
        //For each element construct the reduced matrix
        for(int p = 0; p < n; p++) {
            //Start at the second row
            for(int i = 1; i < n; i++) {
                //Scan all columns
                int k = 0;
                for(int j = 0; j < n; j++) {
                    //if the column number j is the same as the column for the top row element i, then skip this value
                    if(j!=p) {
                        //Row numbers in the reduced matrix are the same as in the original matrix minus one
                        reduced.element[i-1][k] = element[i][j];
                        k++;
                    }
                }
            }
            //Add the product of the current top row element with the reduced matrix
            det=det+element[0][p]*sign*reduced.determinant();
            //Change sign for the next iteration
            sign = -sign;
        }
        return det;
    }
    
}
//============================================================================
