//-----------------------------------------------------------------------------
// Simulacion.cpp
//
// Project: Tephra Irazu simulation - (C) CIC-ITCR
// Author : José Castro
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

#include "VoxelSpace.cpp"
#include "Particula.cpp"

void Simulacion::read(char* in_file)
{
  cout << "Lectura del archivo.\n";
  // Manejo de archivos

  char nombre[30];

  if (in_file == NULL) in_file = _nombre;
  sprintf(nombre, "%s_%d.txt", in_file, _step);

  fstream archivo(nombre, ios::in);
  // Apertura del archivo de puntos

  try
    {
      /* Leer las dimensiones de la Matriz para I, J, K */
      int I, J, K;
      double dx;

      archivo >> I;
      archivo >> J;
      archivo >> K;

      /* Leer la cantidad de partÌculas, la presiÛn, la densidad
       * y la viscosidad
       */
      archivo >> _densidad;
      archivo >> _viscosidad;
      archivo >> dx;
      archivo >> _time1;
      cout << " Viscosidad = " << _viscosidad << '\n';
      cout << " dx = " << dx << '\n';
      _maxDt =  1.0/30.0;
      cout << " maxDt = " << _maxDt << '\n';

      // Declarar la matriz de tres dimensiones
      _voxels.init(I,J,K,dx);

      // Se inicializa el espacio como uno de aire. Se modifica
      // con los datos del archivo (evita que la lectura sea dependiente
      // del orden). En versiones futuras debe mejorase.
      for(VoxelSpace::Cursor cursor = _voxels.begin(); !cursor.finished(); ++cursor)
	cursor->Init(_presion, _viscosidad);

      /* Se crea la matriz de 3 dimensiones */

      /* Iterativamente desde i = 2 se leen los datos del grid
       * Importante: se asume que los datos est·n ordenados crecientemente
       * de acuerdo a cada dimensiÛn
       */
      int    ipos, jpos, kpos, tipo;
      double xvel, yvel, zvel, presion;

      _nLiquidas = 0;
      while (true)
	{
	  archivo >> ipos; if (ipos == -1) break;
	  archivo >> jpos;
	  archivo >> kpos;
	  archivo >> tipo;
	  archivo >> xvel;
	  archivo >> yvel;
	  archivo >> zvel;
	  archivo >> presion;
	  // Si es una celda de fluido, incrementar la cuenta
	  if (tipo == liquido)
	    _nLiquidas++;
	  // Modifica los valores de la celda

	  VoxelSpace::Cursor cursor = _voxels(ipos, jpos, kpos);
	  cursor->tipo = (Estado)tipo;
	  cursor->u[X] = xvel;
	  cursor->u[Y] = yvel;
	  cursor->u[Z] = zvel;
	  cursor->pressure = presion;
	}
    }
  catch (exception e)
    {
      cerr << "Error de estructura del archivo.\n";
    }
  archivo.close();
}

void Simulacion::write(char* out_file)
{
  cout << "Escritura de archivo.\n";
  // Manejo de archivos
  char nombre[30];

  if (out_file == NULL) out_file = _nombre;
  sprintf(nombre, "%s_%d.txt", out_file, _step);

  fstream archivo(nombre, ios::out);
  // Apertura del archivo de puntos

  try
    {
      /* Leer las dimensiones de la Matriz para I, J, K */
      // actual = Regex.Split(lineas[0], " ");
      archivo << _voxels.I() << ' ' << _voxels.J() << ' ' << _voxels.K() << '\n';

      /* Leer la cantidad de partÌculas, la presiÛn, la densidad
       * y la viscosidad
       */
      archivo << _densidad << ' ' << _viscosidad << ' ' << _voxels.dx() << ' ' << _time1 << '\n';

      /* Iterativamente desde i = 2 se leen los datos del grid
       * Importante: se asume que los datos est·n ordenados crecientemente
       * de acuerdo a cada dimensiÛn
       */

      for(VoxelSpace::Cursor cursor = _voxels.begin(); !cursor.finished(); ++cursor)
	{
	  archivo << cursor.i() << ' ' << cursor.j() << ' ' << cursor.k() << ' ';
	  archivo << (int)cursor->tipo << ' ';
	  archivo << cursor->u[X] << ' ';
	  archivo << cursor->u[Y] << ' ';
	  archivo << cursor->u[Z] << ' ';
	  archivo << cursor->pressure << '\n';
	}

      archivo << "-1\n";
    }
  catch (exception e)
    {
      cerr << "Error de estructura del archivo.\n";
    }
  archivo.close();
}

Simulacion::Simulacion(int step, char in_file_prefix[])
{
  _step      = step;
  strcpy(_nombre, in_file_prefix);
  cout << _nombre << '\n';

  read();
  // Se inicializa el tensor de condiciones de cuerpos
  // Inicialmente se considera como la gravedad en el vector de y.
  _condForma = Vector3D(0, -9.806, 0);
}

/* Esta es la funcion que se encarga de efectuar la simulacion */
void Simulacion::simulate(int iteraciones, double tiempo)
{
  _niteraciones = iteraciones;
  _time         = tiempo;
  ofstream sw("timeslices.txt", ios::out);

  try
    {
      // 1. Se reubican las particulas
      // 2. Se aplican las condiciones de borde antes de iniciar la simulaciÛn
      //    para determinar el dominio de las velocidades
      /* 3. Mientras el tiempo de ejecucion sea mayor que cero, ejecute los siguientes pasos:
       *	a. Calcule el valor de dt
       *	b. Reubique las particulas con los datos actuales de las velocidades.
       *	c. Calcule las condiciones de borde
       *	d. Resuelva el valor de las presiones
       *	e. Guarde la informacion de las particulas actuales
       *	f. Guarde el valor de dt para el paso actual
       *	g. time = time - dt
       */
      while ((_time > 0) && (_niteraciones > 0))
	{
	  _dt    = findDt()*0.9;
          _time0 = _time1;
	  relocateParticles();
	  initNewCells();
	  solidBoundaryConditions();
	  emptyBoundaryConditions();
	  cout << "(best guess).\n";
	  calculateNewVelocitiesStep1();
	  recountLiquidCells();
	  for (int i = 0; i < 3; ++i) solveVelocity(i);
	  solvePressures();
	  ++_step;
	  --_niteraciones;
	  _time  -= _dt;
	  _time1 += _dt;
	  savePovRay();
	  cout << "\n step = " << _step << ", tiempo = " << _time1 << " dt = " << _dt  << "\n\n";
	  sw << "Frame: " << _step << ' ' << _dt << '\n';
 	}
      write(); // guarda estado actual de la simulacion
      saveParticles();
    }
  catch(exception e)
    {
      cerr << "Error: ";
    }
  sw.close();
}


/* Este metodo se encarga de distribuir las particulas dentro de las celdas con fluido. La cantidad
 * de particulas por celda es fija, pero su distribucion espacial es aleatoria */
void Simulacion::fillRandomParticles(int cant)
{
  // Se determina la cantidad de particulas por celda que deben existir
  int intPart = _nparticulas = 0;
  double dx   = _voxels.dx();

  // srand(time(NULL));

  // Se llenan las celdas que son liquidas con las particulas
  for(VoxelSpace::Cursor cursor = _voxels.begin(); !cursor.finished(); ++cursor)
    if (cursor->tipo == liquido)
      for (int p = 0; p < cant; p++)
	{
	  _particulas.push_back(Particula(cursor.i(),cursor.j(),cursor.k(),
					  Vector3D(((1.0*random()/RAND_MAX)-0.5+cursor.i())*dx,
						   ((1.0*random()/RAND_MAX)-0.5+cursor.j())*dx,
						   ((1.0*random()/RAND_MAX)-0.5+cursor.k())*dx)
					  )
				);
	  _nparticulas += cant;
	}
}

/* Este metodo se encarga de distribuir las particulas dentro de las celdas con fluido. La cantidad
 * de particulas por celda es fija, pero su distribucion espacial es aleatoria */
void Simulacion::fillParticles(int cant)
{
  // Se determina la cantidad de particulas por celda que deben existir
  int intPart = _nparticulas = 0;
  double dx   = _voxels.dx();
  double dx1  = dx /cant;
  double dx2  = dx1/2;

  // srand(time(NULL));

  int i,j,k;
  double x,y,z;

  // Se llenan las celdas que son liquidas con las particulas
  for(VoxelSpace::Cursor cursor = _voxels.begin(); !cursor.finished(); ++cursor)
    if (cursor->tipo == liquido)
      for (i = 0, x = dx2; i < cant; ++i, x += dx1)
	for (j = 0, y = dx2; j < cant; ++j, y += dx1)
	  for (k = 0, z = dx2; k < cant; ++k, z += dx1)
	    {
	      _particulas.push_back(Particula(cursor.i(),cursor.j(),cursor.k(),
					  Vector3D((cursor.i()-0.5)*dx+x,
						   (cursor.j()-0.5)*dx+y,
						   (cursor.k()-0.5)*dx+z)
					  )
				);
	      _nparticulas += cant;
	    }
}

/* Metodo que calcula el dt a partir de la seccion 3.2.3.1 */
double Simulacion::findDt()
{
  //Variable que guarda el maximo valor
  double max_uvw = 0;
  double result;

  for(VoxelSpace::Cursor cursor = _voxels.begin(); !cursor.finished(); ++cursor)
    //Verifica que las celdas sean liquidas
    if (cursor->tipo == liquido)
      {
	for (int i = 0; i < 3; ++i)
	  if( max_uvw < abs(cursor->u[i]) ) max_uvw = abs(cursor->u[i]);
      }
    else if (cursor->tipo == vacio)
      {
	if ((cursor.i()>0) && (cursor[IZQUIERDA]->tipo == liquido) && (max_uvw < abs(cursor->u[X])))
	  max_uvw = abs(cursor->u[X]);
	if ((cursor.j()>0) && (cursor[ABAJO]->tipo == liquido) && (max_uvw < abs(cursor->u[Y])))
	  max_uvw = abs(cursor->u[Y]);
	if ((cursor.k()>0) && (cursor[ATRAS]->tipo == liquido) && (max_uvw < abs(cursor->u[Z])))
	  max_uvw = abs(cursor->u[Z]);
      }


  cout << "\n Max(u,v,w) = " << max_uvw << '\n';
  /* Se completa la formula 3.19 con una desigualdad */
  if ((max_uvw == 0) || (_maxDt <  _voxels.dx()/max_uvw))
    result = _maxDt;
  else
    result = _voxels.dx()/max_uvw;

  return result;
}

/* Metodo que calcula la nueva ubicacion de las particulas */
void Simulacion::relocateParticles()
{
  for (vector<Particula>::iterator p = _particulas.begin();
       p != _particulas.end(); ++p)
    p->moveParticle(_voxels, _dt);
}

void Simulacion::initNewCells()
{
  double dx  = _voxels.dx();
  double dx2 = dx/2;

  for(VoxelSpace::Cursor cursor = _voxels.begin(); !cursor.finished(); ++cursor)
    if (cursor->newCelda)
      {
	Vector3D u(cursor.i(), cursor.j(), cursor.k());
	u *= dx;
	for (int i = 0; i < 3; ++i)
	  if (cursor[i]->tipo == vacio)
	    {
	      u[i] -= dx2;
	      cursor->u[i] = getClosestU(u)[i];
	      u[i] += dx2;
	    }
      }
    else if (cursor->tipo == vacio)
      {
	Vector3D u(cursor.i(), cursor.j(), cursor.k());
	u *= dx;
	for (int i = 0; i < 3; ++i)
	  if (cursor[i]->newCelda)
	    {
	      u[i] -= dx2;
	      cursor->u[i] = getClosestU(u)[i];
	      u[i] += dx2;
	    }
      }

  for(VoxelSpace::Cursor cursor = _voxels.begin(); !cursor.finished(); ++cursor)
    {
      if (cursor->newCelda)
	cursor->newCelda = false;
      else if (cursor->tipo==liquido)
	cursor->tipo = vacio;
      cursor->borde = false;
    }

  for (vector<Particula>::iterator p = _particulas.begin();
       p != _particulas.end(); ++p)
    _voxels(p->i, p->j, p->k)->tipo = liquido;
}

Vector3D& Simulacion::getClosestU(Vector3D x)
{
  double   dx2 = _voxels.dx()/2;
  double   maxDistance = dx2*_voxels.I();
  double   minDistance = maxDistance;
  Vector3D xmin;
  Vector3D prom(0,0,0);
  int      count = 0;

  for (vector<Particula>::iterator p = _particulas.begin();
       p != _particulas.end(); ++p)
    {
      double dist = distancia(x, p->pos);
      if (dist < minDistance)
	{
	  xmin = p->u;
	  minDistance = dist;
	}
      if (dist <= dx2)
	{
	  prom += p->u;
	  ++count;
	}
    }
  if (count == 0)
    prom = xmin;
  else
    prom *= (1.0/count);

  return prom+Vector3D(0,0,0);
}

void Simulacion::solidBoundaryConditions()
{
  for(VoxelSpace::Cursor cursor = _voxels.begin(); !cursor.finished(); ++cursor)
    {
      if (cursor->tipo == solido)
	{
	  if ((cursor.i()>0) && (cursor[IZQUIERDA]->tipo == liquido))
	    cursor->u = -cursor[IZQUIERDA]->u;
	  else if ((cursor.j()>0) && (cursor[ABAJO]->tipo == liquido))
	    cursor->u = -cursor[ABAJO]->u;
	  else if ((cursor.k()>0) && (cursor[ATRAS]->tipo == liquido))
	    cursor->u = -cursor[ATRAS]->u;

	  if ((cursor.i()>0) && (cursor[IZQUIERDA]->tipo != solido))
	    cursor->u[X] = 0.0;
	  if ((cursor.j()>0) && (cursor[ABAJO]->tipo != solido))
	    cursor->u[Y] = 0.0;
	  if ((cursor.k()>0) && (cursor[ATRAS]->tipo != solido))
	    cursor->u[Z] = 0.0;
	}
      else if (cursor->tipo == liquido)
	{
	  int i;
	  for (i = 0; i < 3; ++i)
	    {
	      if (cursor[i]->tipo == solido)
		{
		  cursor[i]->u = -cursor->u;
		  break;
		}
	    }

	  for (i = 0; i < 3; ++i)
	    if (cursor[i]->tipo == solido)
	      cursor->u[i] = 0.0;
	}
    }
}

void Simulacion::emptyBoundaryConditions()
{
  for(VoxelSpace::Cursor cursor = _voxels.begin(); !cursor.finished(); ++cursor)
    if (cursor->tipo == liquido)
      emptyBoundaryConditions(
			      cursor->u[X], cursor->u[Y], cursor->u[Z],
			      cursor[DERECHA]->u[X],
			      cursor[ARRIBA]->u[Y],
			      cursor[ADELANTE]->u[Z],
			      cursor[IZQUIERDA]->tipo,
			      cursor[ABAJO]->tipo,
			      cursor[ATRAS]->tipo,
			      cursor[DERECHA]->tipo,
			      cursor[ARRIBA]->tipo,
			      cursor[ADELANTE]->tipo,
			      cursor->borde
			      );
}

void Simulacion::emptyBoundaryConditions(
					 Float izquierda, Float abajo, Float atras,
					 Float derecha, Float arriba, Float adelante,
					 Estado eizquierda, Estado eabajo, Estado eatras,
					 Estado ederecha, Estado earriba, Estado eadelante,
					 bool& borde)
{
  Float  cubo[] = {izquierda, abajo, atras, derecha, arriba, adelante};
  Estado ecubo[] = {eizquierda, eabajo, eatras, ederecha, earriba, eadelante};
  bool   utilizable[6];

  int i;
  double sum = 0.0;
  int vacios = 0;
  int utilizables = 0;
  int liquidos = 0;

  for (i = 0; i < 6; ++i)
    {
      if (ecubo[i] != vacio)
	{
	  utilizable[i] = false;
	  if (ecubo[i] == liquido) ++liquidos;
	}
      else
	{
	  vacios++;
	  if (ecubo[opuesto(i)] == liquido)
	    {
	      utilizable[i] = false;
	      cubo[i] = cubo[opuesto(i)];
	    }
	  else
	    {
	      utilizable[i] = true;
	      ++utilizables;
	    }
	}
    }
  if (vacios == 0) return;
  borde = true;
  if (vacios == 6)
    {
      izquierda = izquierda + _condForma[X]*_dt;
      derecha   = derecha   + _condForma[X]*_dt;
      arriba    = arriba    + _condForma[Y]*_dt;
      abajo     = abajo     + _condForma[Y]*_dt;
      adelante  = adelante  + _condForma[Z]*_dt;
      atras     = atras     + _condForma[Z]*_dt;
      return;
    }
  if (liquidos == 0)
    {
      for (int i = 0; i < 3; ++i)
	{
	  if (ecubo[i] == solido)
	    if ((cubo[i]+_condForma[i]) < 0.0)
	      cubo[i] = 0.0;
	    else
	      cubo[i] = cubo[i]+_condForma[i];
	  int k = opuesto(i);
	  if (ecubo[k] == solido)
	    if ((cubo[k]+_condForma[k]) > 0.0)
	      cubo[k] = 0.0;
	    else
	      cubo[k] = cubo[k]+_condForma[k];
	}
    }

  for (i = 0; i < 3; ++i) sum += cubo[opuesto(i)] - cubo[i];

  switch (vacios)
    {
    case 1:
    case 5:
      for (i = 0; i < 3; ++i)
	if ((ecubo[i] == vacio) && (ecubo[opuesto(i)] != vacio))
	  cubo[i] += sum;
	else if ((ecubo[i] != vacio) && (ecubo[opuesto(i)] == vacio))
	  cubo[opuesto(i)] -= sum;
      break;
    case 2:
      if (utilizables == 2) break;
      sum /= 2;
      for (i = 0; i < 3; ++i)
	{
	  if (ecubo[i] == vacio)
	    cubo[i] += sum;
	  if (ecubo[opuesto(i) == vacio])
	    cubo[opuesto(i)] -= sum;
	}
      break;
    case 3:
      if (utilizables == 0) break;
      for (i = 0; i < 3; ++i)
	if ((ecubo[i] == vacio) && (ecubo[opuesto(i)] != vacio))
	  cubo[i] += sum;
	else if ((ecubo[i] != vacio) && (ecubo[opuesto(i)] == vacio))
	  cubo[opuesto(i)] -= sum;
      break;
    case 4:
      sum /= utilizables;

      for (i = 0; i < 3; ++i)
	{
	  if (utilizable[i])
	    cubo[i] += sum;
	  if (utilizable[opuesto(i)])
	    cubo[opuesto(i)] -= sum;
	}
      break;
    }
}

void Simulacion::calculateNewVelocitiesStep1()
{
  for(VoxelSpace::Cursor c = _voxels.begin(); !c.finished(); ++c)
    if (c->tipo == liquido)
      {
	c->unew = c->u; // por si acaso hay bordes que no son liquidos
	for (int i = 0; i < 3; ++i)
	  if (c[i]->tipo == liquido)
	    c->unew[i] = c.bestGuess(i, _dt, _condForma);
      }
}

// OJO estoy aqui
void Simulacion::recountLiquidCells()
{
  _nLiquidas = _nCoord[0] = _nCoord[1] = _nCoord[2] = 0;

  for(VoxelSpace::Cursor c = _voxels.begin(); !c.finished(); ++c)
    {
      for (int i = 0; i < 4; ++i) c->solverIndex[i] = -1;
      if (c->tipo == liquido)
	{
	  c->u = c->unew;
	  c->solverIndex[CELDA] = _nLiquidas;
	  ++_nLiquidas;
	  for (int i = 0; i < 3; ++i)
	    if (c[i]->tipo == liquido)
	      {
		c->solverIndex[i] = _nCoord[i];
		++_nCoord[i];
	      }
	}
    }
}

void Simulacion::solvePressures()
{
  Matrix A; A = Matrices.newMatrix(_nLiquidas);
  Vector b; b = Vectors.newVector(_nLiquidas);
  Vector x; x = Vectors.newVector(_nLiquidas);
  double c = _voxels.dx()*_voxels.dx()*_densidad/_dt;

  for(VoxelSpace::Cursor cursor = _voxels.begin(); !cursor.finished(); ++cursor)
    if (cursor->tipo == liquido)
      {
	int i = cursor->solverIndex[CELDA];
	double diagValue = 6.0;
	for (int v = 0; v < 6; ++v)
	  switch (cursor[v]->tipo)
	    {
	    case liquido: A[i][cursor[v]->solverIndex[3]] = -1.0; break;
	    case solido : diagValue -= 1.0; break;
	    }
	A[i][i] = diagValue;
	b[i] = -c*cursor.divergenceVelocity();
	x[i] = cursor->pressure;
      }
  ConjugateGradientSolve(A, x, b, 1000, 0.00000001);

  for(VoxelSpace::Cursor cursor = _voxels.begin(); !cursor.finished(); ++cursor)
    if (cursor->tipo == liquido)
      cursor->pressure = x[cursor->solverIndex[CELDA]];

  c = _dt / _densidad;
  for(VoxelSpace::Cursor cursor = _voxels.begin(); !cursor.finished(); ++cursor)
    if (cursor->tipo == liquido)
      cursor->u -= cursor.gradientPressure()*c;
}

void Simulacion::solveVelocity(int dir)
{
   if (_nCoord[dir] == 0) return;
  Matrix  D;    D    = Matrices.newMatrix(_nCoord[dir]);
  Vector  u;    u    = Vectors.newVector(_nCoord[dir]);
  Vector  unew; unew = Vectors.newVector(_nCoord[dir]);
  Vector  M;    M    = Vectors.newVector(_nCoord[dir]);
  double c = (_viscosidad*_dt)/(_voxels.dx()*_voxels.dx());

  for(VoxelSpace::Cursor cursor = _voxels.begin(); !cursor.finished(); ++cursor)
    if ((cursor->tipo == liquido) && (cursor[dir]->tipo == liquido))
      {
	int i   = cursor->solverIndex[dir];
	D[i][i] = 1.0 + 6.0*c;
	M[i]    = 1.0 + 6.0*c;
	unew[i] = cursor->u[dir];
	u[i]    = unew[i];
	for (int v = 0; v < 6; ++v)
	  if (cursor[v]->solverIndex[dir] != -1)
	    D[i][cursor[v]->solverIndex[dir]] = -c;
	  else
	    u[i] += cursor[v].get(dir, unew[i])*c;
      }

  if (unew == 0.0) return; // comparacion de vector con una constante

  ConjugateGradientPreconditioner(D, M, unew, u, 1000, 0.00000001);

  for (VoxelSpace::Cursor cursor = _voxels.begin(); !cursor.finished(); ++cursor)
    if ((cursor->tipo == liquido) && (cursor[dir]->tipo == liquido))
      cursor->u[dir] = unew[cursor->solverIndex[dir]];
}

/* Metodo que guarda una instantanea de la posicion actual de las
 * particulas en un archivo para hacer render en PovRay
 */
void Simulacion::savePovRay()
{
  try
    {
      char fname[30];
      sprintf(fname, "%s_%d.pov",_nombre, _step);
      fstream sw(fname, ios::out);

      sw << "#if (clock >= " << _time0 << ")\n";
      sw << "#if (clock <  " << _time1 << ")\n";

      int count = 0;

      for (vector<Particula>::iterator p = _particulas.begin(); p != _particulas.end(); ++p)
	  	++count;

      int total = 12000;
      int mod = ( count > total ? (count / total) : 1);
      int n = 0;

      if (_step == 8) cout << " total = " << total << ", mod = " << mod << '\n';

      for (vector<Particula>::iterator p = _particulas.begin(); p != _particulas.end(); ++p)
	  if ((n++ % mod) == 0)
	  {
	    sw << "sphere {";
	    sw << " <" << p->pos[X] << ',' << p->pos[Y] << ',' << p->pos[Z] << ">, TAM";
	    sw << " texture { pigment { color Cyan } } }\n";
	  }

      sw << "#end\n";
      sw << "#end\n";
      sw.close();
    }
  //MIND: Definir cual excepcion es la que va aqui y en el metodo sig
  catch (exception e)
    {
      cout << "Error: ";
    }
}

/* Metodo que guarda una instantanea de la posicion actual de las
 * particulas en un archivo
 */
void Simulacion::saveParticles()
{
  try
    {
      char fname[30];
      sprintf(fname, "%s_%d.par",_nombre, _step);
      fstream sw(fname, ios::out);

      sw << _particulas.size() << '\n';
      for (vector<Particula>::iterator p = _particulas.begin(); p != _particulas.end(); ++p)
	sw << p->pos[X] << ' ' << p->pos[Y] << ' ' << p->pos[Z] << ' '
	   << p->i << ' ' << p->j << ' ' << p->k << '\n';

      sw.close();
    }
  //MIND: Definir cual excepcion es la que va aqui y en el metodo sig
  catch (exception e)
    {
      cout << "Error: ";
    }
}

/* Metodo que lee una instantanea de la posicion actual de las
 * particulas en un archivo
 */
void Simulacion::loadParticles()
{
  try
    {
      char fname[30];
      sprintf(fname, "%s_%d.par", _nombre, _step);
      fstream sw(fname, ios::in);

      _particulas.clear();
      sw >> _nparticulas;

      for (int i = 0; i < _nparticulas; ++i)
	{
	  double x, y, z;;
	  int i, j, k;
	  sw >> x;
	  sw >> y;
	  sw >> z;
	  sw >> i;
	  sw >> j;
	  sw >> k;
	  _particulas.push_back(Particula(i,j,k,Vector3D(x,y,z)));
	}
      sw.close();
    }
  //MIND: Definir cual excepcion es la que va aqui y en el metodo sig
  catch (exception e)
    {
      cout << "Error: ";
    }
}

/* Metodo que guarda las velocidades en un instante de la simulacion.
 */

void Simulacion::saveVelocity()
{
  try
    {
      char fname[30];
      sprintf(fname, "velocity_%d.txt",_step);
      fstream sw(fname, ios::out);

      for(VoxelSpace::Cursor c = _voxels.begin(); !c.finished(); ++c)
	if (c->tipo == liquido)

	  sw << "Voxel[" << c.i() << ',' << c.j() << ',' << c.k()
	     << "] -> vx = " << c->u[X] << ", vy = " << c->u[Y]
	     << ", vz = " << c->u[Z] << '\n';

      sw.close();
    }

  catch (exception e)
    {
      cerr << "Error: ";
    }
}
