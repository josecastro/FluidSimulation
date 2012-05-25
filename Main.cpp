//-----------------------------------------------------------------------------
// Main.cpp
//
// Project: Tephra Irazu simulation - (C) CIC-ITCR
// Author : Santiago Nu\ez
// Date   : July 15, 2005
//
// Purpose:
//      CPP implementaion of fluid simulation using voxel method
//      Taken from Carlson Lewis CS Georgia Tech Ph.D. dissertation, 2004
//
//-----------------------------------------------------------------------------
#include <assert.h>
#include <time.h>
#include <stdlib.h>
#include "Simulacion.h"

int main(int argc, char* argv[])
{
  if (argc == 1)
    {
      cout << "\nuso: $ " << argv[0] << " <iteraciones> <tiempo> <step> <archivo>\n\n";
      cout << "    iteraciones : cantidad maxima de iteraciones\n";
      cout << "    tiempo      : cantidad maxima de tiempo simulado (en segundos)\n";
      cout << "    step        : etapa inicial de la simulacion\n";
      cout << "    archivo     : prefijo del nombre de archivo de la simulacion\n\n";
      return 0;
    }

  int iteraciones = atoi(argv[1]);
  int tiempo      = atoi(argv[2]);
  int step        = atoi(argv[3]);
  Simulacion simulacion(step, argv[4]);
  simulacion.loadParticles();
  simulacion.write("corrida.dat");
  simulacion.simulate(iteraciones, tiempo);
  return 0;
}
