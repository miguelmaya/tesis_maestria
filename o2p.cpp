/// Para obtener la trayectoria de los pasos, las distancias y ganancias
///  Este es un programa para simular a un agente en busca de 2 cosas (O2P) ///
/// Maneja MILLONES de objetivos usando arreglos globales ///
/// NO usa Clases ni estructuras ///
/// La salida es la dispersion en las distancias y desplazam promediados en un ensamble ///
/// para cada valor de beta. También la distancia (1 y 2) promediado en el ensamble ///
/// Imprime un archivo de salida con las distancias en cada iteración ///
// Me da un archivo de salida de la ganacias para cada valor de beta  //
// En este caso minimizo (l_ij + l_jm)/(k_j + k_m)  //


#include <iostream>  // Para salida a la pantalla
#include <fstream>  // Para la salida de archivos
#include <cstdlib> // para los numeros aleatorios
#include <ctime> // para tomar datos del reloj
#include <cmath>
#include <iomanip> // para dar formato a la salida
#include <algorithm>  //std::fill
#include <vector>  // para asignar la longit de un array con una var
#include <sstream>  // para CONVERTIR NUMBER TO STRING

using namespace std;

// Aquí van las CONSTANTES //
const int num_mon = 1;  
const int num_tree = 10000;     // número de blancos presentes
const int kmax = 1000;          // valor de corte del tamaño de los blancos
const int maxit = 250;          
const int maxdob = 2*maxit;     // número de blancos visitados 2*maxit
const int escala = 10000;
const int ensamb = 50;          // Ensamble sobre el cual promediar
const int ta = 2;               // número de valores de beta a utilizar


float tree_x[num_tree], tree_y[num_tree];
int tree_k[num_tree];

// Aquí van las FUNCIONES a utilizar //

float sqr(float num)
{
	return num*num;	
}

float Alea()
{
	return rand()/(float)RAND_MAX;
}


float Promedio(float  matriz[ ], float limite) // aquí calculo el coeficiente de variacion
{
	 int j; 
	 float prom = 0, promcuad = 0, tot;  

     for ( j = 0; j < limite; j++)
     {
		 prom += matriz[j];
		 promcuad += sqr(matriz[j]);
	 }
	 tot = (promcuad/limite - sqr(prom/limite))/sqr(prom/limite);
	 return tot;
}

float Media (float matr[])
{
	float sum = 0;
	for(int i = 0; i < ensamb; i++)
		sum += matr[i];
	return sum/ensamb;
}

float Suma (float matr[])
{
	float sum = 0;
	for(int i = 0; i < maxit; i++)
		sum += matr[i];
	return sum;
}

float DesvEst (float matriz[ ], float promed)
{
	float tran = 0;
	for (int k = 0; k < ensamb; k++)
		tran += sqr(matriz[k] - promed);
	return sqrt(tran/ensamb);
}

template <typename T>
string NTS ( T Number )
{
    stringstream ss;
    ss << Number;
    return ss.str();
}

int main   ()
{
    //srand48(time(NULL));  //la semilla de la función rand()
    srand48(132);  //la semilla de la función rand()
	float p[kmax];  // puedo teclear kmax desde la linea de comandos
	float c, ran, pe, ratio, ratio_min, b;
    float ratio_1, ratio_2, transit;
	float monkey_x[num_mon], monkey_y[num_mon];
	float beta[ta] ={2.9,3.0};                       /// VALORES DE BETA A UTILIZAR
	float dist1[maxit], dist2[maxit], desp[maxit], dist[maxdob];
	float sigm1[ensamb], sigm2[ensamb], sigm3[ensamb], sigm4[ensamb], distrec1[ensamb], distrec2[ensamb], desprec[ensamb];
	float proms1, proms2, proms3, proms4;
	int cont, limite, it, primpaso, segpaso, fruit1, fruit2;
	float fruta, fruta1, fruta2,frutaobt, l_0, l_1, l_2, l_tot, dp1, dp2, dp3;
	int i, k, j, w, bi;
	
	ofstream f1 ("2promdist.txt");
	ofstream f2 ("2promdesp.txt");
	ofstream f3 ("2distdesprec.dat");
	
	std::string fname1 = "2posi";
	std::string fname2 = "rep";
	std::string ext = ".dat";
	string filename3 = "2ganancia";
	
	for (int val = 0; val < ta; val++)
	{
		b = beta[val];
		bi = ceil(b*100);
		ofstream ganancia;
		ganancia.open((filename3 + NTS(bi) + ext).c_str());
		
		for (int rep = 0; rep < ensamb; rep++)
		{
			ofstream f4;   // imprimir las posiciones en cada iteración
			f4.open((fname1 + NTS(bi) + fname2 + NTS(rep) + ext).c_str());
			monkey_x[0] = 0.5;
			monkey_y[0] = 0.5;
  
			for ( i = 0; i < num_tree; i++)
			{
				tree_x[i] = drand48();
				tree_y[i] = drand48();
			}
			c = 0;
			fruta = 0;
			for ( k = 1; k <= kmax; k++)
				c += pow(k, -b);
         
			pe = 0;
			for ( k = 1; k <= kmax; k++)
			{
				pe += pow(k, -b)/c;
				p[k-1] = pe;
			}
			cont = 0;
			while ( cont < num_tree)
			{
				ran = drand48();
				i = 0;
				while (( ran >= p[i] ) && (i < (kmax-1)))
					i ++;
				tree_k[cont] = i + 1;  // debe empezar en 1 no en 0
				fruta += tree_k[cont];
				cont++;
			}
			fruta1 = 0;
			fruta2 = 0;
			frutaobt = 0;
			l_0 = 1 / sqrt(num_tree);
			
			f4 << setiosflags(ios::fixed);
       
			for (it = 0; it < maxit; it++)
			{
				ratio_min = 1e10;
				transit = 1e10;
				for ( i = 0; i < num_tree; i++)
				{
					if ( tree_k[i] > 0 )
					{
						l_1 = sqrt(sqr(monkey_x[0]-tree_x[i]) + sqr(monkey_y[0]-tree_y[i]));
						ratio_1 = l_1/tree_k[i];
						//ratio_2 = (l_1 + 0.00003)/(tree_k[i] + 1000);
						//if (ratio_2 < transit )  // optimización
						//{
							for ( w = 0; w < num_tree; w++)
							{
					     
								if (tree_k[w] > 0 && w != i)
								{
									l_2 = sqrt(sqr(tree_x[i]-tree_x[w]) + sqr(tree_y[i]-tree_y[w])); 
							
									ratio = (l_1 + l_2)/(tree_k[i] + tree_k[w]);  //minimizo la suma de las distancias
									if (ratio < ratio_min )
									{
										ratio_min = transit = ratio;
										primpaso = i;
										segpaso = w;                                          
									}
								}
							}
						//}
					}     
				}    
           
				l_1 = sqrt(sqr(monkey_x[0]-tree_x[primpaso]) + sqr(monkey_y[0]-tree_y[primpaso]));
				l_2 = sqrt(sqr(tree_x[primpaso]-tree_x[segpaso]) + sqr(tree_y[primpaso]-tree_y[segpaso]));
				l_tot = l_1 + l_2;
                     
				monkey_x[0] = tree_x[segpaso];  // el agente se para en el segundo objetivo
				monkey_y[0] = tree_y[segpaso];  // el agente se para en el segundo objetivo  
    
				fruta1 += tree_k[primpaso];
				fruta2 += tree_k[segpaso];
				frutaobt += tree_k[primpaso];
                f4 << setw(10) << setprecision(6)  << tree_x[primpaso] << '\t';
                f4 << setw(10) << setprecision(6)  << tree_y[primpaso] << '\t';
                f4 << setw(10) << setprecision(6)  << l_1 << '\t';
                f4 << setw(10) << tree_k[primpaso] << '\t';
                f4 << setw(8) << setprecision(0) << frutaobt << endl;
                frutaobt += tree_k[segpaso];
                f4 << setw(10) << setprecision(6)  << tree_x[segpaso] << '\t';
                f4 << setw(10) << setprecision(6)  << tree_y[segpaso] << '\t';
                f4 << setw(10) << setprecision(6)  << l_2 << '\t';
                f4 << setw(10) << tree_k[segpaso] << '\t';      
                f4 << setw(8) << setprecision(0) << frutaobt << endl;          
           				
				
				tree_k[primpaso] = 0;  // no regresa a un lugar ya visitado
				tree_k[segpaso] = 0;  // no regresa a un lugar ya visitado
				
				dist1[it] = l_1;
				dist2[it] = l_2;
				dist[2*it] = l_1;
				dist[2*it+1] = l_2;
				desp[it] = l_tot;

			}
			
			f4.close();
			
			// Ahora calculo la dispersión para un valor de beta y una corrida //
			
			sigm1[rep] = Promedio(dist, maxdob);
			sigm2[rep] = Promedio(dist1, maxit);
			sigm3[rep] = Promedio(dist2, maxit);
			sigm4[rep] = Promedio(desp, maxit);
			
			// Aquí guardo la suma de las distancias y desp recorridas //
			
			distrec1[rep] = Suma(dist1);
			distrec2[rep] = Suma(dist2);
			desprec[rep] = Suma(desp); 		
			
			// Guardo la fruta comida y la distancia recorrida
			ganancia << fruta << '\t' << fruta1 << '\t' << fruta2 << '\t' << distrec1[rep] << '\t' << distrec2[rep] << endl;							
		}
		ganancia.close();
		// Promedio de la Dispersión en las distancias //
		proms1 = Media(sigm1);
		proms2 = Media(sigm2);
		proms3 = Media(sigm3);
		
		f1 << setiosflags(ios::fixed);
		f1 << setw(4) << setprecision(2) << b << "\t";
		f1 << setw(10) << setprecision(6)  << proms1 << "\t";
		f1 << setw(10) << setprecision(6)  << DesvEst(sigm1, proms1) << "\t";
		f1 << setw(10) << setprecision(6)  << proms2 << "\t";
		f1 << setw(10) << setprecision(6)  << DesvEst(sigm2, proms2) << "\t";
		f1 << setw(10) << setprecision(6)  << proms3 << "\t";
		f1 << setw(10) << setprecision(6)  << DesvEst(sigm3, proms3) << "\n";
		// Promedio de la Dispersión en los desplazamientos //
		proms4 = Media(sigm4);
		f2 << setiosflags(ios::fixed);		
		f2 << setw(4) << setprecision(2) << b << "\t";		
		f2 << setw(10) << setprecision(6)  << proms4 << "\t";		
		f2 << setw(10) << setprecision(6)  << DesvEst(sigm4, proms4) << "\n";
		
		// Promedio de las ditacias (1 y 2) y de los desplazamientos recorridos //
		dp1 = Media(distrec1);
		dp2 = Media(distrec2);
		dp3 = Media(desprec);
		f3 << setiosflags(ios::fixed);		
		f3 << setw(4) << setprecision(2) << b << "\t";		
		f3 << setw(10) << setprecision(6)  << dp1 << "\t";		
		f3 << setw(10) << setprecision(6)  << DesvEst(distrec1, dp1) << "\t";		
		f3 << setw(10) << setprecision(6)  << dp2 << "\t";
		f3 << setw(10) << setprecision(6)  << DesvEst(distrec2, dp2) << "\t";
		f3 << setw(10) << setprecision(6)  << dp3 << "\t";
		f3 << setw(10) << setprecision(6)  << DesvEst(desprec, dp3) << "\n";
                              
     }
     f1.close();
     f2.close();
     f3.close();

     //cin.get();
}
 
