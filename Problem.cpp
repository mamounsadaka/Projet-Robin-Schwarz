#ifndef _PB_CPP
#include <string>
#include <fstream>
#include <iostream>
#include "GradConj.h"
#include "Problem.h"

/*-----------------------------------------------------
Classe de construction de la matrice du problème, du second terme complet avec les conditions de bords ainsi que la résolution du système linéaire par gradient conjugué.
Seul la fonction Solve serai utilisée en main
*/

// Macros de représentation et de debuggage

#define debug std::cout << "step here" << std::endl;

#define siz(a) std::cout << a.size() << std::endl;

// Fonctions de représentation de vecteur
void print_vector2(std::vector<double> x)
{
	int n = x.size();
	std::cout << "le vecteur de taille " << n << std::endl;

	for (int i = 0; i < n; i++)
	{
		std::cout << x[i] << " ";
	}
	std::cout << std::endl;
	std::cout << "----------------------------------" << std::endl;
}

// Constructeur
Problem::Problem(BC *f, int Nx, int Ny, int Nt, double Lx, double Ly, double deltat, int procID, int Np, int n, double r, double Lx0, double alpha, double beta, std::vector<double> stencil1, std::vector<double> stencil2) : functions_(f), procID_(procID), Np_(Np), r_(r), n_(n), stencil1_(stencil1), stencil2_(stencil2)
{
	alpha_ = alpha;
	beta_ = beta;
	Lx0_ = Lx0;
	Lx_ = Lx;
	Ly_ = Ly;
	Nx_ = Nx;
	Ny_ = Ny;
	Nt_ = Nt;
	deltat_ = deltat;
	std::cout << "Construction de la classe Matrice" << std::endl;
	std::cout << "-------------------------------------------------" << std::endl;
}

// Construit la matrice des flux
std::vector<std::vector<double>> Problem::Construct_Matrix()
{
	// la matrice serai 3 vecteurs diag,sdiag et ssdiag. Elle sera de cinq vecteurs en version parallèle mais plus sur ça dans la partie //
	std::vector<std::vector<double>> A(5);

	// les pas de l'espace
	deltax_ = Lx_ / (Nx_ + 1.);
	deltay_ = Ly_ / (Ny_ + 1.);

	double phix = -D_ * deltat_ / (deltax_ * deltax_);
	double phiy = -D_ * deltat_ / (deltay_ * deltay_);
	double alpha = 1 - 2 * (phix + phiy);

	for (int i = 0; i < Nx_ * Ny_; i++)
	{

		if (alpha_ != 0)
		{
			if (i % Nx_ == 0 && procID_ > 0)
			{
				A[0].push_back(alpha + 2 * (beta_ * deltax_ / alpha_) * phix);
			}
			else
			{
				if ((i + 1) % Nx_ == 0 && procID_ < Np_ - 1)
				{
					A[0].push_back(alpha - 2 * (beta_ * deltax_ / alpha_) * phix);
				}
				else
				{
					A[0].push_back(alpha);
				}
			}
		}
		else
		{
			A[0].push_back(alpha);
		}
	}

	for (int i = 0; i < Nx_ * Ny_ - 1; i++)
	{
		if ((i + 1) % Nx_ == 0)
		{
			A[1].push_back(0.);
		}
		else
		{
			if (alpha_ != 0)
			{
				if (i % Nx_ == 0 && procID_ > 0)
				{
					A[1].push_back(2 * phix);
				}
				else
				{
					A[1].push_back(phix);
				}
			}
			else
			{
				A[1].push_back(phix);
			}
		}
	}
	A[1].push_back(0.0);
	A[3].push_back(0.0);
	for (int i = 1; i < Nx_ * Ny_; i++)
	{
		if (i % Nx_ == 0)
		{
			A[3].push_back(0.);
		}
		else
		{
			if (alpha_ != 0)
			{
				if ((i + 1) % Nx_ == 0 && procID_ < Np_ - 1)
				{
					A[3].push_back(2 * phix);
				}
				else
				{
					A[3].push_back(phix);
				}
			}
			else
			{
				A[3].push_back(phix);
			}
		}
	}

	for (int i = 0; i < Nx_ * Ny_; i++)
	{
		if (i < Nx_ * (Ny_ - 1))
		{
			A[2].push_back(phiy);
		}
		else
		{
			A[2].push_back(0.0);
		}
	}
	for (int i = 0; i < Nx_ * Ny_; i++)
	{
		if (i >= Nx_)
		{
			A[4].push_back(phiy);
		}
		else
		{
			A[4].push_back(0.0);
		}
	}
	A_ = A;
	return A;
}

// Construction du terme source
void Problem::Construct_F(int cas, double t, std::vector<double> &test)
{
	int n = Nx_ * Ny_;
	F_.resize(n);
	std::vector<double> z(n);
	for (int i = 0; i < Ny_; i++)
	{
		for (int j = 0; j < Nx_; j++)
		{
			F_[i * Nx_ + j] = (sol_[i * Nx_ + j] + deltat_ * functions_->Source_term(procID_ * deltax_ * (Nx_ + 1. - n_) + (j + 1) * deltax_, (i + 1.0) * deltay_, t + deltat_, cas));
		}
	}

	test.resize(n);
	test = F_;
}

// Ajout des conditions de bords
void Problem::Construct_Bd(int cas, double t)
{
	int n = Nx_ * Ny_;
	std::vector<double> M(n, 0.);
	Bd_.resize(n);
	Bd_ = M;
	for (int i = 0; i < Ny_; i++)
	{
		for (int j = 0; j < Nx_; j++)
		{
			if (i == 0 && j == 0)
			{
				if (procID_ == 0)
				{
					Bd_[i * Nx_ + j] = D_ * functions_->Dirichlet_Function0(procID_ * deltax_ * (Nx_ + 1 - n_) + (j + 1) * deltax_, i * deltay_, t + deltat_, cas) / (deltax_ * deltax_) + D_ * functions_->Dirichlet_Function1(procID_ * deltax_ * (Nx_ + 1 - n_) + (j + 1) * deltax_, i * deltay_, t + deltat_, cas) / (deltay_ * deltay_);
				}
				else
				{
					Bd_[i * Nx_ + j] = D_ * stencil1_[i] / (deltax_ * deltax_) + D_ * functions_->Dirichlet_Function1(procID_ * deltax_ * (Nx_ + 1 - n_) + (j + 1) * deltax_, i * deltay_, t + deltat_, cas) / (deltay_ * deltay_);
				}
			}
			if (i == 0 && j != 0 && j != Nx_ - 1)
			{
				Bd_[i * Nx_ + j] = D_ * functions_->Dirichlet_Function1(procID_ * deltax_ * (Nx_ + 1 - n_) + (j + 1) * deltax_, i * deltay_, t + deltat_, cas) / (deltay_ * deltay_); // D_*functions_->Dirichlet_Function0(procID_ *deltax_ * (Nx_ + 1 - n_) + (j + 1) * deltax_, i * deltay_, t, cas) / (deltax_ * deltax_) +
			}
			if (i == 0 && j == Nx_ - 1)
			{
				if (procID_ != Np_ - 1)
				{
					Bd_[i * Nx_ + j] = D_ * stencil2_[i] / (deltax_ * deltax_) + D_ * functions_->Dirichlet_Function1(procID_ * deltax_ * (Nx_ + 1 - n_) + (j + 1) * deltax_, i * deltay_, t + deltat_, cas) / (deltay_ * deltay_);
				}
				else
				{
					Bd_[i * Nx_ + j] = D_ * functions_->Dirichlet_Function0(procID_ * deltax_ * (Nx_ + 1 - n_) + (j + 1) * deltax_, i * deltay_, t + deltat_, cas) / (deltax_ * deltax_) + D_ * functions_->Dirichlet_Function1(procID_ * deltax_ * (Nx_ + 1 - n_) + (j + 1) * deltax_, i * deltay_, t + deltat_, cas) / (deltay_ * deltay_);
				}
			}
			if (i == Ny_ - 1 && j == 0)
			{
				if (procID_ != 0)
				{
					Bd_[i * Nx_ + j] = D_ * stencil1_[i] / (deltax_ * deltax_) + D_ * functions_->Dirichlet_Function1(procID_ * deltax_ * (Nx_ + 1 - n_) + (j + 1) * deltax_, i * deltay_, t + deltat_, cas) / (deltay_ * deltay_);
				}
				else
				{
					Bd_[i * Nx_ + j] = D_ * functions_->Dirichlet_Function0(procID_ * deltax_ * (Nx_ + 1 - n_) + (j + 1) * deltax_, i * deltay_, t + deltat_, cas) / (deltax_ * deltax_) + D_ * functions_->Dirichlet_Function1(procID_ * deltax_ * (Nx_ + 1 - n_) + (j + 1) * deltax_, i * deltay_, t + deltat_, cas) / (deltay_ * deltay_);
				}
			}
			if (i == Ny_ - 1 && j == Nx_ - 1)
			{
				if (procID_ != Np_ - 1)
				{
					Bd_[i * Nx_ + j] = D_ * stencil2_[i] / (deltax_ * deltax_) + D_ * functions_->Dirichlet_Function1(procID_ * deltax_ * (Nx_ + 1 - n_) + (j + 1) * deltax_, i * deltay_, t + deltat_, cas) / (deltay_ * deltay_);
				}
				else
				{
					Bd_[i * Nx_ + j] = D_ * functions_->Dirichlet_Function0(procID_ * deltax_ * (Nx_ + 1 - n_) + (j + 1) * deltax_, i * deltay_, t + deltat_, cas) / (deltax_ * deltax_) + D_ * functions_->Dirichlet_Function1(procID_ * deltax_ * (Nx_ + 1 - n_) + (j + 1) * deltax_, i * deltay_, t + deltat_, cas) / (deltay_ * deltay_);
				}
			}
			if (i == Ny_ - 1 && j != 0 && j != Nx_ - 1)
			{
				Bd_[i * Nx_ + j] = D_ * functions_->Dirichlet_Function1(procID_ * deltax_ * (Nx_ + 1 - n_) + (j + 1) * deltax_, i * deltay_, t + deltat_, cas) / (deltay_ * deltay_);
			}
			if (j == 0 && i != 0 && i != Ny_ - 1)
			{
				if (procID_ == 0)
				{
					Bd_[i * Nx_ + j] = D_ * functions_->Dirichlet_Function0(procID_ * deltax_ * (Nx_ + 1 - n_) + (j + 1) * deltax_, i * deltay_, t + deltat_, cas) / (deltax_ * deltax_); //+ D_ * functions_->Dirichlet_Function1(procID_ *deltax_ * (Nx_ + 1 - n_) + (j + 1) * deltax_, i * deltay_, t, cas) / (deltay_ * deltay_);
				}
				else
				{
					Bd_[i * Nx_ + j] = D_ * stencil1_[i] / (deltax_ * deltax_);
				}
			}
			if (j == Nx_ - 1 && i != 0 && i != Ny_ - 1)
			{
				if (procID_ == Np_ - 1)
				{
					Bd_[i * Nx_ + j] = D_ * functions_->Dirichlet_Function0(procID_ * deltax_ * (Nx_ + 1 - n_) + (j + 1) * deltax_, i * deltay_, t + deltat_, cas) / (deltax_ * deltax_); //+ D_ * functions_->Dirichlet_Function1(procID_ *deltax_ * (Nx_ + 1 - n_) + (j + 1) * deltax_, i * deltay_, t, cas) / (deltay_ * deltay_);
				}
				else
				{
					Bd_[i * Nx_ + j] = D_ * stencil2_[i] / (deltax_ * deltax_);
				}
			}
		}
	}
}

// Résolution du problème, schéma du mémoire du projet
void Problem::Solve_problem(int cas, double tf)
{
	int n = Nx_ * Ny_;
	sol_.resize(n);
	for (int i = 0; i < n; i++)
	{
		sol_[i] = 0.0;
	}

	// solve Au_=f with gradconj::solve after initialsie
	int nb_iter = 10000;
	double t(0.);
	std::vector<double> test;
	std::vector<double> S(Nx_ * Ny_, 0.), S1(Nx_ * Ny_, 0.);
	while (t < tf)
	{
		A_ = Problem::Construct_Matrix(); // writes twice to change
		Problem::Construct_F(cas, t, test);
		Problem::Construct_Bd(cas, t);
		S1 = GradConj::prod_scal(Bd_, deltat_);
		S = GradConj::sum(F_, S1, 1);
		GradConj gc = GradConj(A_, S, Nx_, Ny_);
		if(alpha_==0)
		{
			gc.Solve(nb_iter, sol_);
		}
		else {
			gc.TSolve(nb_iter, sol_);
		}
		t += deltat_;
	}
}
// Récupération de la solution
std::vector<double> Problem::get_sol()
{
	std::vector<double> x(sol_.size());
	x = sol_;
	return sol_;
}
std::vector<double> Problem::getLeftStencil()
{
	std::vector<double> x(Ny_);
	for (int i = 1; i < Ny_; i++)

	{
		if (alpha_ != 0)
		{
			x[i] = sol_[i * Nx_ + n_ - 1] - sol_[i * Nx_ + n_ - 3] + (2 * beta_ * deltax_ / alpha_) * sol_[i * Nx_ + n_ - 2];
		}
		else
		{
			// x[i - 1] = sol_[i * Nx_ + n_];
			x[i] = sol_[i * Nx_ + n_ - 1];
		}
	}
	return x;
}
std::vector<double> Problem::getRightStencil()
{
	std::vector<double> x(Ny_);
	for (int i = 0; i < Ny_; i++)
	{
		if (alpha_ != 0)
		{
			// x[i - 1] = sol_[i * Nx_+Nx_-n_-1] - sol_[i * Nx_+Nx_+1-n_] + (2 * beta_ * deltax_ / alpha_) * sol_[i * Nx_ +Nx_-n_];
			x[i] = sol_[i * Nx_ + Nx_ - n_] - sol_[i * Nx_ + Nx_ + 2 - n_] - (2 * beta_ * deltax_ / alpha_) * sol_[i * Nx_ + Nx_ + 1 - n_];
		}
		else
		{
			x[i] = sol_[i * Nx_ + Nx_ - n_];
			// x[i - 1] = sol_[i * Nx_ + Nx_ - 1 - n_];
		}
	}
	return x;
}

// double Problem::transformX ()
// {
// 	double del= x+procID*(Lx_-r_);
// }

#define _PB_CPP
#endif
