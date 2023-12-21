#include <solve.hpp>        // everything from ngsolve
#include <cmath>

#include "intrules.hpp"
#include "ngbem.hpp"
#include "hmat.hpp"
#include "test_compression.hpp"



namespace ngbem
{
  using namespace ngcomp;


  Matrix<> CreateMatrix(int nx, int ny, double eta)
  {
    Array<Vec<3>> ax(nx), ay(ny);

    for (auto & xi : ax)
      for (double & val : xi)
        val = double (rand()) / RAND_MAX;

    for (auto & yi : ay)
      for (double & val : yi)
        val = double (rand()) / RAND_MAX;
    for (auto & yi : ay)
      yi(0) += eta;
    
    
    // cout << "ax = " << ax << endl;
    // cout << "ay = " << ay << endl;
    Matrix mat(nx, ny);
    for (int i = 0; i < nx; i++)
      for (int j = 0; j < ny; j++)
        mat(i,j) = 1.0 / L2Norm(ax[i]-ay[j]);

    return mat;
  }


  
  tuple<int, double, double> TestCompressionSVD (int nx, int ny, double eta, double eps)
  {
    int num = 0;

    auto mat = CreateMatrix(nx, ny, eta);
    Matrix<> savemat = mat;
    
    Matrix<> U(nx, nx), V(ny,ny);
    Vector<> S(std::min(nx,ny));
    
    Timer t("compression");
    t.Start();
    LapackSVD (mat, U, V, S, false);
    t.Stop();

    for (auto si : S)
      if (fabs(si) > eps)
        num++;


    for (int i = 0; i < num; i++)
      {
        U.Col(i) *= sqrt(S(i));
        V.Row(i) *= sqrt(S(i));
      }
    double err = L2Norm(savemat - U.Cols(num)*V.Rows(num));
    
    return { num, err, t.GetTime() };
  }




  template <typename T>
  void MyStochasticTSVD1 (MatrixView<T> A, MatrixView<T> U, MatrixView<T> V, VectorView<> S)
  {
    for (int i = 0; i < V.Height(); i++)
      for (int j = 0; j < V.Width(); j++)
        V(i,j) = double (rand()) / RAND_MAX;

    Matrix<T> AVt = A * Trans(V);
    Matrix<T> tmp(V.Height(), V.Height());
    LapackSVD (AVt, U, tmp, S, false);    
    Matrix<T> UtA = Trans(U) * A;
    LapackSVD (UtA, tmp, V, S, false);
    U = Matrix<T> (U * tmp);
  }

  template <typename T>
  size_t MyStochasticTSVD (MatrixView<T> A, MatrixView<T> U, MatrixView<T> V, VectorView<> S, double eps)
  {
    static Timer tsvd("ngbem - StochasticTSVD"); 
    RegionTimer reg(tsvd);
    
    int rank = 5;
    int p = min(A.Height(), A.Width());
    while (rank < p)
      {
        MyStochasticTSVD1 (A, U.Cols(rank), V.Rows(rank), S.Range(rank));
        if (S[rank-1] < eps)
          {
            for (int j = 1; j < p; j++)
              if (S[j] < eps)
                return j-1;
          }

        rank = int (rank*1.5);
      }

    LapackSVD(A, U, V, S, false);
    for (int j = 1; j < p; j++)
      if (S[j] < eps)
        return j-1;

    return p;
  }
  

  
  tuple<int, double, double> TestCompressionTSVD (int nx, int ny, double eta, double eps)
  {
    auto mat = CreateMatrix(nx, ny, eta);

    Matrix<> savemat = mat;
    Matrix<> U(nx, nx), V(ny,ny);
    Vector<> S(std::min(nx,ny));
    
    Timer t("compression");
    t.Start();
    int num = MyStochasticTSVD<double> (mat, U, V, S, eps);
    t.Stop();

    for (int i = 0; i < num; i++)
      {
        U.Col(i) *= sqrt(S(i));
        V.Row(i) *= sqrt(S(i));
      }
    double err = L2Norm(savemat - U.Cols(num)*V.Rows(num));
    return { num, err, t.GetTime() };
  }

  
  void TestFunc (FlatMatrix<double> U, FlatMatrix<double> V, int k, int jmax)
  {
    U.Col(k) -= U.Cols(0,k) * V.Col(jmax).Range(0,k);    
  }


  template <typename MATRIX, typename T, ORDERING O, typename ...TU, typename TV>
  int CalcACA (const MATRIX & mat, MatrixView<T,O,TU...> U, FlatMatrix<TV> V, double eps)
  {
    for (int k = 0; k < U.Width(); k++)
      {
        int ik = k;  // what else ?
        
        V.Row(k) = mat.Row(ik);
        // for (int l = 0; l < k; l++)
        // V.Row(k) -= U(ik,l) * V.Row(l);

        V.Row(k) -= Trans(V.Rows(0,k)) * U.Row(ik).Range(0,k);
         
        double err = L2Norm(V.Row(k));
        // cout << "Norm vk = " << err << endl;
        if (err < eps) return k;
        
        int jmax = 0;
        for (int j = 0; j < V.Width(); j++)
          if (fabs (V(k,j)) > fabs(V(k,jmax)))
            jmax = j;
        V.Row(k) *= 1.0 / V(k,jmax);
        U.Col(k) = mat.Col(jmax);

        // for (int l = 0; l < k; l++)
        // U.Col(k) -= V(l,jmax) * U.Col(l);

        U.Col(k) -= U.Cols(0,k) * V.Col(jmax).Range(0,k);
      }
    
    return U.Width();
  }
  

  int CalcACA1 (FlatMatrix<> mat, FlatMatrix<double> U, FlatMatrix<double> V, double eps)

  {
    return CalcACA (mat, U, V, eps);
  }
  
  tuple<int, double, double> TestCompressionACA (int nx, int ny, double eta, double eps)
  {
    auto mat = CreateMatrix(nx, ny, eta);

    Matrix<> savemat = mat;
    Matrix<> U(nx, nx), V(ny,ny);
    Vector<> S(std::min(nx,ny));
    
    Timer t("compression");
    t.Start();

    int num = CalcACA (mat, U, V, eps);

    t.Stop();

    double err = L2Norm(savemat - U.Cols(num)*V.Rows(num));
    return { num, err, t.GetTime() };
  }
  


}
