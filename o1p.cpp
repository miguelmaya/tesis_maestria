///  Este es un programa para simular un Agente FORRAJEANDO ///
/// Maneja MILLONES de árboles usando arreglos globales ///
/// NO usa Clases ni estructuras ///
/// optimización a primer paso ///

#include <iostream>
#include <fstream>
#include <cstdlib> // para los numeros aleatorios
#include <time.h> // para tomar datos del reloj
#include <cmath>
#include <iomanip> // para dar formato a la salida
#include <algorithm>  //std::fill
#include <vector>  // para asignar la longit de un array con una var

using namespace std;

const int num_mon = 1;  
const int num_tree = 1000000;    // número de blancos presentes
const int kmax = 1000;           // valor de corte del tamaño de los blancos
const int maxit = 1000;          // número de blancos visitados
float b = 3; 	                 // valor de beta

float tree_x[num_tree], tree_y[num_tree];
int tree_k[num_tree];

float sqr(float num)
{
	return num*num;	
}

int main   ()
{
    srand48(time(NULL));  //la semilla de la función rand()
	float p[kmax];  // puedo teclear kmax desde la linea de comandos
	float c, ran, pe, ratio, ratio_min;

	float monkey_x[num_mon], monkey_y[num_mon];
	int cont, it, k_min;
	float ran2, fruta, l, l_0;
	int i, k, j, optimo;
	ofstream f1 ("power7.txt");
	ofstream f2 ("monos7.txt");
	
    monkey_x[0] = 0.5;
    monkey_y[0] = 0.5;

	for ( i = 0; i < num_tree; i++)
	{
        tree_x[i] = drand48();
		tree_y[i] = drand48();
    }
    cout << " Asignando tamaño de árboles \n";
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
           f1 << setiosflags(ios::fixed);
           f1 << setw(10) << cont+1 << "\t";
           f1 << setw(10) << tree_k[cont] << "\n";
           cont++;
    }
     
    f1.close();
    cout << "...listo !!!";
    cout << " Tamaño de árbol k máximo es: " << kmax << endl;
    cout << " Fruta disponible: " << fruta << endl;
    cout << "Exponente asignado: " << b << endl;
    cout << " Número de iteraciones: " << maxit << endl;
     
    l_0 = 1 / sqrt(num_tree);
    cout << " El valor de l_0 es: " << l_0 << endl;
    cout << " Número de monos: " << num_mon << endl;
    cout << "Número de árboles: " << num_tree << endl;   
    cout << " Comiendo...";
     
    for (it = 0; it < maxit; it++)
    {
         ratio_min = 1e10;
         for ( i = 0; i < num_tree; i++)
         {
             if ( tree_k[i] > 0 )
             {
                  l = sqrt(sqr((monkey_x[0]-tree_x[i])) + sqr((monkey_y[0]-tree_y[i])));
                  ratio = l/(l_0 * tree_k[i]);
                  if (ratio < ratio_min )
                  {
                      ratio_min = ratio;
                      optimo = i;
                  }
             }
         }
         monkey_x[0] = tree_x[optimo];
         monkey_y[0] = tree_y[optimo];

         f2 << setiosflags(ios::fixed);
         f2 << setw(10) << it+1 << "\t";
         f2 << setw(15) << setprecision(7) << monkey_x[0] <<'\t';
         f2 << setw(15) << setprecision(7) << monkey_y[0] << endl;
         
         tree_k[optimo] = 0;
    }
    f2.close();
    cout << " Salida de datos en power7.txt y monos7.txt \n";
     //cin.get();
}
 
