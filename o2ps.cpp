/// Programa para obtener la trayectoria del caminante que optimiza a segundo paso (o2p)
///  Este es un programa para simular a un agente en busca de 2 cosas (O2P) ///
/// Maneja MILLONES de objetivos usando arreglos globales ///
/// NO usa Clases ni estructuras ///
/// La salida es la dispersion en las distancias y desplazam promediados en un ensamble ///
/// para cada valor de beta. También la distancia (1 y 2) promediado en el ensamble ///
/// Imprime un archivo de salida con las distancias en cada iteración ///
// Me da un archivo de salida con la "trayectoria.dat que contiene 6 columnas dende:  //
// 1) posición_x, 2)posición_y; 3)distancia_recorrida 4)recurso_captirado 5)recursos acumulados
/// Los renglones van intercalados 1,3,5,7,... información al visitar el primer blanco de cada iteración
/// los 2,4,6,8,... la información al visitar el segundo blanco de cada iteración 
// En este caso minimizo (l_ij + l_jm)/(k_j + k_m)  //

///  Programó: Miguel Angel Maya Contreas, julio de 2014



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
const int num_tree = 1000;     // número de blancos presentes
const int kmax = 1000;          // valor de corte del tamaño de los blancos
const int maxit = 250;          
const int maxdob = 2*maxit;     // número de blancos visitados 2*maxit
const int escala = 10000;

float tree_x[num_tree], tree_y[num_tree];
int tree_k[num_tree];

// Aquí van las FUNCIONES a utilizar //

float sqr(float num)
{
	return num*num;	
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
    srand48(time(NULL));  //la semilla de la función rand()
	float p[kmax];  // puedo teclear kmax desde la linea de comandos
	float c, ran, pe, ratio, ratio_min;
    float ratio_1, ratio_2, transit;
	float monkey_x[num_mon], monkey_y[num_mon];
	float dist1[maxit], dist2[maxit], desp[maxit], dist[maxdob];
	float proms1, proms2, proms3, proms4;
	int cont, limite, it, primpaso, segpaso, fruit1, fruit2;
	float fruta, fruta1, fruta2,frutaobt, l_0, l_1, l_2, l_tot, dp1, dp2, dp3;
	int i, k, j, w, bi;
	float b = 3.0;  /// valor del parámetro beta

	ofstream f4 ("trayectoria.dat");   // imprimir las posiciones en cada iteración

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
		cont++;
	}
	fruta1 = 0;
	fruta2 = 0;
	frutaobt = 0;
	l_0 = 0.5 / sqrt(num_tree);  // distancia caracterítica
	
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
		
	}		
	f4.close();
}
 
