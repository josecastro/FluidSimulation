//-----------------------------------------------------------------------------
// Generador.cpp
//
// Project: Tephra Irazu simulation - (C) CIC-ITCR
// Author : José Castro
// Date   : July 15, 2005
//
// Purpose:
//      C++ simulation generation files
//
//-----------------------------------------------------------------------------

#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <stdlib.h>
#include "Simulacion.h"

using namespace std;

void generateCubeScenario(int scenario, char file[], int sizeCube, double dx, double viscosidad);

void generateCubeScenario1(int sizeCube, ofstream& sw);
void generateCubeScenario2(int sizeCube, ofstream& sw);
void generateCubeScenario3(int sizeCube, ofstream& sw);
void generateCubeScenario4(int sizeCube, ofstream& sw);
void generateCubeScenario5(int sizeCube, ofstream& sw);
void generateCubeScenario6(int sizeCube, ofstream& sw);

int main(int argc, char* argv[])
{
  try
    {
      if (argc == 1)
	{
	  cout << "\nuso : $ " << argv[0]
	       << " <tipo> <tamano> <particulas> <dx> <viscosidad> <archivo>\n\n";
	  cout << "    tipo       : 1 - cubo con agua suspendida en el medio\n";
	  cout << "                 2 - cubo con piramide de agua\n";
	  cout << "                 3 - cubo con cubo de agua suspendida en una esquina superior\n";
	  cout << "                 4 - cubo con agua estable\n";
	  cout << "                 5 - bloque grande de agua suspendida\n";
	  cout << "                 6 - cubo de agua suspendido sobre agua\n";
	  cout << "    tamano     : cantidad de voxels en cada lado de el cubo (total = tamano^3)\n";
	  cout << "    particulas : catidad de particulas por voxel\n";
	  cout << "    dx         : tamano de los voxels\n";
	  cout << "    viscosidad : viscosidad cinematica\n";
	  cout << "    archivo    : prefijo del nombre de archivo en que se guarda escenario\n\n";
	  return 0;
	}

      // Tipo y archivo a generar
      int tipoGen       = atoi(argv[1]);
      int sizeCube      = atoi(argv[2]);
      int particulas    = atoi(argv[3]);
      double dx         = atof(argv[4]);
      double viscosidad = atof(argv[5]);
      char* archivo  = argv[6];

      generateCubeScenario(tipoGen, archivo, sizeCube, dx, viscosidad);

      Simulacion sim(0, archivo);
      sim.fillRandomParticles(particulas);
      sim.saveParticles();
    }
  catch (exception e)
    {
      cerr << "Uso: generador <tama\\o-escenario> <archivo-salida>\n";
    }
  return 0;
}

void generateCubeScenario(int scenario, char file[], int sizeCube, double dx, double viscosidad)
{
  // Archivo de escritura de las pruebas
  char name[30];
  sprintf(name, "%s_0.txt", file);
  ofstream sw(name, ios::out);

  sw << sizeCube << ' ' << sizeCube << ' ' << sizeCube << '\n';
  sw << "1 " << viscosidad << ' ' << dx << " 0.0\n"; // densidad, viscosidad, dx, tiempo1

  for (int i = 0; i < sizeCube; i++)
    for (int j = 0; j < sizeCube; j++)
      {
	sw << i << ' ' << j << " 0 1 0 0 0 0\n";
	sw << i << ' ' << j << ' ' << sizeCube-1 << " 1 0 0 0 0\n";
	sw << "0 " << i << ' ' << j << " 1 0 0 0 0\n";
	sw << sizeCube-1 << ' ' << i << ' ' << j << " 1 0 0 0 0\n";
	sw << i << " 0 " << j << " 1 0 0 0 0\n";
	sw << i << ' ' << sizeCube-1 << ' ' << j << " 1 0 0 0 0\n";
      }

  switch (scenario)
    {
    case 1: generateCubeScenario1(sizeCube, sw); break;
    case 2: generateCubeScenario2(sizeCube, sw); break;
    case 3: generateCubeScenario3(sizeCube, sw); break;
    case 4: generateCubeScenario4(sizeCube, sw); break;
    case 5: generateCubeScenario5(sizeCube, sw); break;
    case 6: generateCubeScenario6(sizeCube, sw); break;
    }

  sw << "-1\n\n\n";
  sw.close();
}

/* Metodo que utiliza el generador NxN */
void generateCubeScenario1(int sizeCube, ofstream& sw)
{
  // Archivo de escritura de las pruebas

  int i0, j0, k0, sizeWater;
  i0 = j0 = k0 = sizeWater = sizeCube / 3;

  /* En este paso, se crean los voxels líquidos con las velocidades
   * iniciales en cero, y las caras que son extremos poseen las
   * velocidades indefinidas
   */

  double velx = 0;
  double vely = 0;
  double velz = 0;

  for (int i = 0; i < sizeWater; i++)
    for (int j = 0; j < sizeWater; j++)
      for (int k = 0; k < sizeWater; k++)
	// Se escriben las velocidades para los voxels. Si son internos
	// poseen velocidades inicializadas en cero. En caso contrario las
	// velocidades están indefinidas
	sw << i+i0 << ' ' << j+j0 << ' ' << k+k0 << " 0 "
	   << velx << ' ' << vely << ' ' << velz << " 0 \n";
}

/* Metodo que utiliza el generador NxN */
void generateCubeScenario2(int sizeCube, ofstream& sw)
{
  int sizeWater = 0;

  for (int j = 1; j < sizeCube-1; ++j)
    for (int i = j; i < sizeCube-1; ++i)
      for (int k = j; k < sizeCube-1; ++k)
	++sizeWater;

  /* En este paso, se crean los voxels líquidos con las velocidades
   * iniciales en cero, y las caras que son extremos poseen las
   * velocidades indefinidas
   */

  double velx = 0;
  double vely = 0;
  double velz = 0;

  for (int j = 1; j < sizeCube-1; ++j)
    for (int i = j; i < sizeCube-1; ++i)
      for (int k = j; k < sizeCube-1; ++k)
	// Se escriben las velocidades para los voxels. Si son internos
	// poseen velocidades inicializadas en cero. En caso contrario las
	// velocidades están indefinidas
	sw << i << ' ' << j << ' ' << k << " 0 "
	   << velx << ' ' << vely << ' ' << velz << " 0 \n";
}

/* Metodo que utiliza el generador NxN */
void generateCubeScenario3(int sizeCube, ofstream& sw)
{

  int sizeWater = (sizeCube-2 - sizeCube/2);

  sizeWater = sizeWater*sizeWater*sizeWater;

  /* En este paso, se crean los voxels líquidos con las velocidades
   * iniciales en cero, y las caras que son extremos poseen las
   * velocidades indefinidas
   */

  double velx    = 0.0;
  double vely    = 0.0;
  double velz    = 0.0;
  double presion = 0.0;

  for (int j = sizeCube/2; j < sizeCube-2; ++j)
    for (int i = sizeCube/2; i < sizeCube-2; ++i)
      for (int k = sizeCube/2; k < sizeCube-2; ++k)
	// Se escriben las velocidades para los voxels. Si son internos
	// poseen velocidades inicializadas en cero. En caso contrario las
	// velocidades están indefinidas
	sw << i << ' ' << j << ' ' << k << " 0 "
	   << velx << ' ' << vely << ' ' << velz << ' ' << presion << '\n';
}

/* Metodo que utiliza el generador NxN */
void generateCubeScenario4(int sizeCube, ofstream& sw)
{
  int i0, j0, k0, sizeWaterX, sizeWaterY, sizeWaterZ;

  i0 = j0 = k0 = 1;
  sizeWaterZ = sizeWaterX = sizeCube - 2;
  sizeWaterY = sizeCube/2;


  /* En este paso, se crean los voxels líquidos con las velocidades
   * iniciales en cero, y las caras que son extremos poseen las
   * velocidades indefinidas
   */

  double velx = 0;
  double vely = 0;
  double velz = 0;

  for (int i = 0; i < sizeWaterX; i++)
    for (int j = 0; j < sizeWaterY; j++)
      for (int k = 0; k < sizeWaterZ; k++)
	// Se escriben las velocidades para los voxels. Si son internos
	// poseen velocidades inicializadas en cero. En caso contrario las
	// velocidades están indefinidas
	sw << i+i0 << ' ' << j+j0 << ' ' << k+k0 << " 0 "
	   << velx << ' ' << vely << ' ' << velz << " 0\n";
}

/* Metodo que utiliza el generador NxN */
void generateCubeScenario5(int sizeCube, ofstream& sw)
{

  int sizeWater = (sizeCube-2 - sizeCube/2);

  sizeWater = sizeWater*sizeWater*sizeWater;

  /* En este paso, se crean los voxels líquidos con las velocidades
   * iniciales en cero, y las caras que son extremos poseen las
   * velocidades indefinidas
   */

  double velx    = 0.0;
  double vely    = 0.0;
  double velz    = 0.0;
  double presion = 0.0;

  for (int j = sizeCube/2; j < sizeCube-2; ++j)
    for (int i = 2; i < sizeCube-2; ++i)
      for (int k = 2; k < sizeCube-2; ++k)
	// Se escriben las velocidades para los voxels. Si son internos
	// poseen velocidades inicializadas en cero. En caso contrario las
	// velocidades están indefinidas
	sw << i << ' ' << j << ' ' << k << " 0 "
	   << velx << ' ' << vely << ' ' << velz << ' ' << presion << '\n';
}

/* Metodo que utiliza el generador NxN */
void generateCubeScenario6(int sizeCube, ofstream& sw)
{

  int sizeWater = (sizeCube-2 - sizeCube/2);

  sizeWater = sizeWater*sizeWater*sizeWater;

  /* En este paso, se crean los voxels líquidos con las velocidades
   * iniciales en cero, y las caras que son extremos poseen las
   * velocidades indefinidas
   */

  double velx    = 0.0;
  double vely    = 0.0;
  double velz    = 0.0;
  double presion = 0.0;

  for (int j = sizeCube/2; j < sizeCube-2; ++j)
    for (int i = sizeCube/2; i < sizeCube-2; ++i)
      for (int k = sizeCube/2; k < sizeCube-2; ++k)
	// Se escriben las velocidades para los voxels. Si son internos
	// poseen velocidades inicializadas en cero. En caso contrario las
	// velocidades están indefinidas
	sw << i << ' ' << j << ' ' << k << " 0 "
	   << velx << ' ' << vely << ' ' << velz << ' ' << presion << '\n';

  int i0, j0, k0, sizeWaterX, sizeWaterY, sizeWaterZ;

  i0 = j0 = k0 = 1;
  sizeWaterZ = sizeWaterX = sizeCube - 2;
  sizeWaterY = sizeCube/6;


  /* En este paso, se crean los voxels líquidos con las velocidades
   * iniciales en cero, y las caras que son extremos poseen las
   * velocidades indefinidas
   */

  for (int i = 0; i < sizeWaterX; i++)
    for (int j = 0; j < sizeWaterY; j++)
      for (int k = 0; k < sizeWaterZ; k++)
	// Se escriben las velocidades para los voxels. Si son internos
	// poseen velocidades inicializadas en cero. En caso contrario las
	// velocidades están indefinidas
	sw << i+i0 << ' ' << j+j0 << ' ' << k+k0 << " 0 "
	   << velx << ' ' << vely << ' ' << velz << " 0\n";

}

