//
//  main.cpp
//  Ising Model
//
//  Created by Íris on 23/06/2018.
//  Copyright © 2018 Íris e João. All rights reserved.
//

#include <iostream>
#include <vector>
#include <sstream>
#include <iostream>
#include <array>
#include <math.h>
#include <iomanip>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <ctime>
#include <fstream>
#include <time.h>
#include <random>

using namespace std;
using std::vector;
using std::cout;
using std::endl;



int main(int argc, const char * argv[]) {
    

    ofstream file1("/Users/Damiao/Desktop/energia1_L18.txt");
    ofstream file2("/Users/Damiao/Desktop/magnetização1_L18.txt");
    ofstream file3("/Users/Damiao/Desktop/calor1_L18.txt");
    ofstream file4("/Users/Damiao/Desktop/qui1_L18.txt");
    
    std::random_device rd;    // Generates a random number
    std::mt19937 gen(rd());
    
    
    int L = 6; // Size of the lattice
    int N = pow (L,2);
    int z = 4; // número de vizinhos
    double r;
    
    
    
    // Inicializa a rede de dimensão N.
    
    int matriz_inicial [N];
    int nn [N][z];
    
    // Helical- bouundary conditions - Identificação dos vizinhos para cada spin i
    
    for (int i=0; i<=(N-1);i++){
        matriz_inicial [i] = 1;
        
        //Vizinhos do spin i
        
        nn[i][0] = (i+1)%(N);
        nn[i][1] = ((i-1)%(N)+(N))%(N);
        nn[i][2] = (i+L)%(N);
        nn[i][3] = ((i-L)%(N) + (N))%(N);
    }
    
    
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // IMPLEMENTAR MÉTODO DE MONTE CARLO METROPOLIS
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    double J = 1;
    
    //Intervalo de temperatura considerado
    
    double tStep=0.05;
    double Tmax=5;
    
    int nmonte = 10000000;
    int nrelax = 10000000;
    
    int r1; //dois números aleatórios
    std::uniform_int_distribution<> dis_int(0,N-1);
    double r2;
    
    int energias [z+1];
    double d_e;
    double A;
    
    
    for(double T=0.1;T<=Tmax;T=T+tStep){
         cout << T << endl;
        // Relaxação
        for (int n=0; n<=nrelax; n++){
           
            d_e = 0;
            r1 = dis_int(gen); // 1 passo gerar um flip aleatório
            d_e = 2*J*matriz_inicial[r1]*(matriz_inicial[nn[r1][0]]+matriz_inicial[nn[r1][1]]+matriz_inicial[nn[r1][2]]+matriz_inicial[nn[r1][3]]);
            if (d_e <= 0){
                matriz_inicial [r1] = (-1)*matriz_inicial[r1];
            }
            else {
                A = exp((-1/T)*(d_e));
                
                r2 = std::generate_canonical<double, 10>(gen);
                
                if (r2<A){
                    matriz_inicial [r1] = (-1)*matriz_inicial[r1];
                }
            }
        }
        
        // Equilibrio - Relaxação - Onde se calculam as variáveis termodinâmicas
        
        double m=0; //Magnetização total
        double h=0; //Energia Total - Hamiltoniano
        
        // Para calcular os valor de energia do sistema inicial assim como de magnetização
        
        for (int i=0; i<=(N-1);i++){
            h = h -((0.5)*matriz_inicial[i]*(matriz_inicial[nn[i][0]]+matriz_inicial[nn[i][1]]+matriz_inicial[nn[i][2]]+matriz_inicial[nn[i][3]]));
            m = m + matriz_inicial[i];
        }
        
    
        double H = 0;
        double M = 0;
        double H2 = 0;
        double M2 = 0;
        
        for (int n=0; n<nmonte; n++){
            d_e=0;
            r1 = dis_int(gen); // 1 passo gerar um flip aleatório
            d_e = 2*J*matriz_inicial[r1]*(matriz_inicial[nn[r1][0]]+matriz_inicial[nn[r1][1]]+matriz_inicial[nn[r1][2]]+matriz_inicial[nn[r1][3]]);
            if (d_e <= 0){
                matriz_inicial [r1] = (-1)*matriz_inicial[r1];
                h = h + d_e;
                 m = m + 2*matriz_inicial[r1];
                
            }
            else {
                    A = exp((-1/T)*(d_e));
                    
                    r2 = std::generate_canonical<double, 10>(gen);
                if (r2<A){
                    
                    matriz_inicial [r1] = (-1)*matriz_inicial[r1];
                    h = h + d_e;
                    m = m + 2*matriz_inicial[r1];
                                    }
            }
            H = (H + h);
            M = (M + fabs(m));
            
            H2 = H2+h*h;
            M2 = M2+m*m;
        }
        
        double E_av = (H/N)/(nmonte); //Energia média por spin
        double M_av = ((M)/N)/(nmonte); // Magnetização média por spin
        
        H2 = H2/(nmonte);
        M2 = M2/(nmonte);
        
        double E1 = (H)/(nmonte);
        double M1 = (M)/(nmonte);
        
        double  c =(H2-(E1*E1))/(T*T*N); // Calor especifico por spin
        double  Xi = (M2-(M1*M1))/(T*N); // Susceptibilidade por spin
        
        
        file1 << T << "  " << E_av << endl;;
        file2 << T << "  " << M_av<< endl;
        file3 << T << "  " << c << endl;
        file4 << T << "  " << Xi << endl;
        
       
        
    }
    return 0;
}

