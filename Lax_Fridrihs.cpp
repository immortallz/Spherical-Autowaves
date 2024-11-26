#include <iostream>
#include <vector>
#include <cmath>

#define PI 3.141592653589793
#define gamma 1.4
#define k_m 4.15723130907662
// #define k //TODO
// #define m 
using namespace std;

double minmod(double x, double y)
{
	return (2*double(x>0) - 1) * max(0.0, min(abs(x), y*(2*double(x>0) - 1)));
}

double Qn(double rho, double T)
{
	if(rho < 0 || T < 0)
	{
		printf("pow(rho or T):");
		return -1;
	}
	// return rho * (0.685*(2.46*pow(rho, 0.8) - rho*T) * pow(T, 0.5) - pow(T, 0.4)*pow(rho, 0.17));
	return 0;
}

double SR(double E1, double E2, double E3)
{
	if(abs(E1) < 1e-13)
	{
		cout << "E1 is zero!" << endl;
		return -1;
	}
	if(E3/E1 - 0.5*E2*E2/E1/E1 < 0)
	{
		printf("sqrt(negative): %lf - %lf\n", E3/E1 , 0.5*E2*E2/E1/E1);
		return -1;
	}
	// else
		// printf("sqrt(positive): %lf - %lf\n", E3/E1 , 0.5*E2*E2/E1/E1);
	return max(max(
		abs(E2/E1 - sqrt(gamma*(gamma-1)*(E3/E1 - 0.5*E2*E2/E1/E1))),
		abs(E2/E1)),
		abs(E2/E1 + sqrt(gamma*(gamma-1)*(E3/E1 - 0.5*E2*E2/E1/E1)))
		);
}

double F1(double E1, double E2, double E3)
{
	return E2;
}

double F2(double E1, double E2, double E3)
{
	return E2*E2/E1 + (gamma-1)*(E3 - 0.5*E2*E2/E1);
}

double F3(double E1, double E2, double E3)
{
	return 0.5*E2*E2*E2/E1/E1 + gamma * E2/E1*(E3 - 0.5*E2*E2/E1);
}

double R1(double E1, double E2, double E3, double r)
{
	return 0;
}

double R2(double E1, double E2, double E3, double r)
{
	return 2*(gamma-1)*E1/r*(E3/E1 - 0.5*E2*E2/E1/E1);
}

double R3(double E1, double E2, double E3, double r)
{
	return Qn(E1/r/r, (gamma-1)/k_m*(E3/E1 - 0.5*E2*E2/E1/E1))*r*r;
}

double u_(double E1, double E2, double E3, double r)
{
	return E2/E1;
}

double rho_(double E1, double E2, double E3, double r)
{
	return E1/r/r;
}

double T_(double E1, double E2, double E3, double r)
{
	double e = E3/E1 - 0.5*E2*E2/E1/E1;
	return e*(gamma-1)/k_m;
}

int main()
{
	FILE *f = fopen("aa.txt", "w"), *ff = fopen("bb.txt","w"),
		*u0_file = fopen("u0.txt", "w"),
		*u1_file = fopen("u1.txt", "w"),
		*u2_file = fopen("u2.txt", "w"),
		*u3_file = fopen("u3.txt", "w");
	int M = 100; // N - t coordinate; M - r coordinate
	vector<vector<vector<double>>>
		E(3, vector<vector<double>>(1, vector<double>(M, 0))),
		F(3, vector<vector<double>>(1, vector<double>(M, 0))),
		R(3, vector<vector<double>>(1, vector<double>(M, 0)));
	double rho_0 = 1, T_0 = 1, A = 1, r1 = 50, r_eps = 10, SR_max = -1;
	double tau = 1, L = 100, dt, dr = L/double(M), CFL = 0.5;
	double r, e, u, t = 0;
	for(int j = 0; j < M; j++)
	{
		r = j*dr + r_eps;
		e = 1;
		if(j < int(double(M) * r1 / L) || j >= int(double(M) * (r1 + 1) / L))
			u = 0;
		else
			u = A * sin(PI*(r - (r1+r_eps)));
		E[0][0][j] = rho_0*r*r;
		E[1][0][j] = rho_0*u*r*r;
		E[2][0][j] = rho_0*(u*u/2 + e)*r*r;

		F[0][0][j] = rho_0*u*r*r;
		F[1][0][j] = rho_0*(u*u + (gamma-1)*e)*r*r;
		F[2][0][j] = rho_0*u*(u*u/2 + gamma*e)*r*r;

		R[0][0][j] = 0;
		R[1][0][j] = 2*rho_0*(gamma-1)*e*r;
		R[2][0][j] = Qn(rho_0, T_0)*r*r;

		fprintf(u0_file, "%lf %lf\n", j*dr + r_eps, u);
	}

	for(int j = 0; j < M; j++)
		if(SR_max < SR(E[0][0][j], E[1][0][j], E[2][0][j]))
			SR_max = SR(E[0][0][j], E[1][0][j], E[2][0][j]);
	printf("SR_max = %lf\n", SR_max);
	dt = CFL*dr/SR_max;
	printf("dt = %lf\n", dt);

	vector<vector<double>>
		ER(3, vector<double>(M+1, 0)), // -1/2, 1/2, ... , M+1/2
		EL(3, vector<double>(M+1, 0)); // ER[s][k] ~ k - 1/2
	vector<double> fdif(3);

	fprintf(ff, "u[%d,%d] = %lf\n", 0, 0, u_(E[0][0][0],E[1][0][0],E[2][0][0], r_eps));
	int i = 0;
	while(t < tau)
	{
		//predictor
		for(int s = 0; s < 3; s++)
		{
			E[s].push_back(vector<double>(M, 0));
			// E[s][i+1][0] = E[s][i][0] - 0.5*dt/dr*(F[s][i][1] - F[s][i][0]) + 0.5*dt*R[s][i][0];
			// E[s][i+1].push_back(0.0);
			for(int j = 1; j < M - 1; j++)
				E[s][i+1][j] = E[s][i][j] - 0.25*dt/dr*(F[s][i][j+1] - F[s][i][j-1]) + 0.5*dt*R[s][i][j];
			// E[s][i+1].push_back(0.0);
			
			// снос значений на индексы 0 и M-1 (с 1 и M-2 соотв.)
			E[s][i+1][0] = E[s][i+1][1];
			E[s][i+1][M - 1] = E[s][i+1][M - 2];

			for(int j = 1; j < 3; j++)
				if(s==1)
					printf("%d it: %lf = %lf - %lf + %lf\n", j, E[s][i+1][j],E[s][i][j],0.25*dt/dr*(F[s][i][j+1] - F[s][i][j-1]),0.5*dt*R[s][i][j]);
			//streams

			// ER[s][0] = E[s][i+1][0] - 0;
			// EL[s][0] = E[s][i+1][0] + 0;
			// ER[s][1] = E[s][i+1][1] - 0.5*minmod(E[s][i+1][1] - E[s][i+1][0], E[s][i+1][2] - E[s][i+1][1]);
			// EL[s][1] = E[s][i+1][0] + 0;
			for(int j = 2; j < M - 1; j++)
			{
				ER[s][j] = E[s][i+1][j] - 0.5*minmod(E[s][i+1][j] - E[s][i+1][j-1], E[s][i+1][j+1] - E[s][i+1][j]);
				EL[s][j] = E[s][i+1][j-1] + 0.5*minmod(E[s][i+1][j-1] - E[s][i+1][j-2], E[s][i+1][j] - E[s][i+1][j-1]);
				
			}
			// ER[s][M-1] = E[s][i+1][M-1] - 0;
			// EL[s][M-1] = E[s][i+1][M-2] + 0.5*minmod(E[s][i+1][M-2] - E[s][i+1][M-3], E[s][i+1][M-1] - E[s][i+1][M-2]);
			// ER[s][M] = E[s][i+1][M-1] - 0;
			// EL[s][M] = E[s][i+1][M-1] + 0;
		}
		printf("after predictor %lf %lf %lf\n", E[1][i+1][0], E[1][i+1][1], E[1][i+1][2]);

		//corrector
		for(int j = 2; j < M - 2; j++)
		{
			fdif[0] = 0.5 * (F1(ER[0][j+1], ER[1][j+1], ER[2][j+1]) +
				F1(EL[0][j+1], EL[1][j+1], EL[2][j+1]) -
				SR(ER[0][j+1], ER[1][j+1], ER[2][j+1]) * (ER[0][j+1] - EL[0][j+1])) -
				0.5 * (F1(ER[0][j], ER[1][j], ER[2][j]) +
				F1(EL[0][j], EL[1][j], EL[2][j]) -
				SR(ER[0][j], ER[1][j], ER[2][j]) * (ER[0][j] - EL[0][j]));

			fdif[1] = 0.5 * (F2(ER[0][j+1], ER[1][j+1], ER[2][j+1]) +
				F2(EL[0][j+1], EL[1][j+1], EL[2][j+1]) -
				SR(ER[0][j+1], ER[1][j+1], ER[2][j+1]) * (ER[1][j+1] - EL[1][j+1])) -
				0.5 * (F2(ER[0][j], ER[1][j], ER[2][j]) +
				F2(EL[0][j], EL[1][j], EL[2][j]) -
				SR(ER[0][j], ER[1][j], ER[2][j]) * (ER[1][j] - EL[1][j]));

			fdif[2] = 0.5 * (F3(ER[0][j+1], ER[1][j+1], ER[2][j+1]) +
				F3(EL[0][j+1], EL[1][j+1], EL[2][j+1]) -
				SR(ER[0][j+1], ER[1][j+1], ER[2][j+1]) * (ER[2][j+1] - EL[2][j+1])) -
				0.5 * (F3(ER[0][j], ER[1][j], ER[2][j]) +
				F3(EL[0][j], EL[1][j], EL[2][j]) -
				SR(ER[0][j], ER[1][j], ER[2][j]) * (ER[2][j] - EL[2][j]));

			for(int s = 0; s < 3; s++){
				E[s][i+1][j] = E[s][i][j] - dt/dr*fdif[s] + dt*R[s][i][j];
				if(s==1 && j < 3)
					printf("%d it: %lf = %lf - %lf + %lf\n", j, E[s][i+1][j],E[s][i][j],dt/dr*fdif[s],dt*R[s][i][j]);
			}
		}
		// снос значений на 0, 1, M-2, M-1 узлы (с 2 и M-3 соотв.)
		for(int s = 0; s < 3; s++)
		{
			E[s][i+1][0] = E[s][i+1][2];
			E[s][i+1][1] = E[s][i+1][2];
			E[s][i+1][M-2] = E[s][i+1][M-3];
			E[s][i+1][M-1] = E[s][i+1][M-3];
		}
		printf("after corrector %lf %lf %lf\n", E[1][i+1][0], E[1][i+1][1], E[1][i+1][2]);

		//F and R update
		for(int s = 0; s < 3; s++)
		{
			F[s].push_back(vector<double>(M, 0));
			R[s].push_back(vector<double>(M, 0));
		}
		for(int j = 0; j < M; j++)
		{
			F[0][i+1][j] = F1(E[0][i+1][j], E[1][i+1][j], E[2][i+1][j]);
			F[1][i+1][j] = F2(E[0][i+1][j], E[1][i+1][j], E[2][i+1][j]);
			F[2][i+1][j] = F3(E[0][i+1][j], E[1][i+1][j], E[2][i+1][j]);

			R[0][i+1][j] = R1(E[0][i+1][j], E[1][i+1][j], E[2][i+1][j], j*dr + r_eps);
			R[1][i+1][j] = R2(E[0][i+1][j], E[1][i+1][j], E[2][i+1][j], j*dr + r_eps);
			R[2][i+1][j] = R3(E[0][i+1][j], E[1][i+1][j], E[2][i+1][j], j*dr + r_eps);
		}

		
		fprintf(ff, "u[%d,%d] = %lf\n", i+1, 0, u_(E[0][i+1][0],E[1][i+1][0],E[2][i+1][0], 0*dr+r_eps));

		// CFL
		SR_max = -1;
		for(int j = 2; j < M - 1; j++)
			if(SR_max < SR(ER[0][j], ER[1][j], ER[2][j]))
				SR_max = SR(ER[0][j], ER[1][j], ER[2][j]);
		dt = CFL*dr/SR_max;
		t += dt;
		fprintf(f, "%lf %lf %d\n", dt, t, i);
		for(int j = 2; j < M - 1; j++)
		{
			fprintf(f, "SR[%d][%d] = %lf\n", i+1, j, SR(ER[0][j], ER[1][j], ER[2][j]));
			if(SR(ER[0][j], ER[1][j], ER[2][j]) == -1)
			{
				printf("vse ploho\n");
				return -1;
			}
		}
		switch(i){
			case 100:
				for(int j = 0; j < M; j++)
					fprintf(u1_file, "%lf %lf\n", j*dr + r_eps, u_(E[0][i][j], E[1][i][j], E[2][i][j], 0));
				break;
			case 200:
				for(int j = 0; j < M; j++)
					fprintf(u2_file, "%lf %lf\n", j*dr + r_eps, u_(E[0][i][j], E[1][i][j], E[2][i][j], 0));
				break;
			case 300:
				for(int j = 0; j < M; j++)
					fprintf(u3_file, "%lf %lf\n", j*dr + r_eps, u_(E[0][i][j], E[1][i][j], E[2][i][j], 0));
				break;
		}

		i++;
	}
	fclose(f);
	cout << E.size() << " " << E[0].size() << " " << E[0][0].size() << endl;
	return 0;
}