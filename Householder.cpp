#include "Householder.h"
#include "GradConj.h"
#include <iostream>
#include <math.h>
#include <iostream>

using namespace std;

Householder::Householder(vector<vector<double>> A, vector<double> b) : A_(A), b_(b)
{
}
void Householder::Decomposition_QR()
{
    // initialisation des paramètres

    double beta = 0;
    int n = A_.rows();
    int m = A_.cols();
    vector<vector<double, n>> P(m * n), P1(n); // matrices  qui permettent de stocker le Pk
    R_.resize(n, m);
    Q_.resize(n, n);
    vector<vector<double, n>> I(n);

    //I = MatrixXd::Identity(n, n);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            
        }
        
    }
    
    Q_ = I;
    R_ = A_;
    VectorXd z(n), w(n);

    // algorithme

    for (int k = 0; k < m; k++)
    {
        if (k > 0)
        {
            for (int i = 0; i < k; i++)
            {

                P1 = P.block(0, i * n, n, n);
                R_.col(k) = P1 * R_.col(k);
            }
        }
        beta = 0;
        if (A_(k, k) != 0)
        {
            for (int i = k; i < n; i++)
            {
                beta = sqrt(beta * beta + A_(i, k) * A_(i, k));
            }

            beta = (A_(k, k) / abs(A_(k, k))) * beta;
        }

        // construction du vecteur Z

        for (int j = 0; j < n; j++)
        {
            if (j < k)
            {
                z(j) = 0;
            }
            if (j == k)
            {
                z(j) = beta + A_(j, j);
            }
            if (j > k)
            {
                z(j) = A_(j, k);
            }
        }

        w = (1 / sqrt(z.dot(z))) * z;   // construction de Pk
        P1 = I - 2 * w * w.transpose(); // construction de Pk
        P.block(0, k * n, n, n) = P1;   // stocker la matrice Pk
        R_.col(k) = P1 * R_.col(k);     // mise à jour de rk
        // mise à jour de qk
        for (int i = 0; i < k + 1; i++)
        {
            P1 = P.block(0, (k - i) * n, n, n);
            Q_.col(k) = P1 * Q_.col(k);
        }
    }
}
VectorXd Householder::Solve()
{
    systeme_triangulaire_sup u(R_, b_);
    VectorXd x = u.Solve();
    x = Q_.transpose() * x;
    return x;
}
MatrixXd Householder::get_Q()
{
    return Q_;
}
MatrixXd Householder::get_R()
{
    return R_;
}