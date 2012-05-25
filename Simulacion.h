//-----------------------------------------------------------------------------
// Simulacion.h
//
// Project: Tephra Irazu simulation - (C) CIC-ITCR
// Author : José Castro
// Date   : July 15, 2005
//
// Purpose:
//      CPP specification of fluid simulation using voxel method
//      Taken from Carlson Lewis CS Georgia Tech Ph.D. dissertation, 2004
//
//-----------------------------------------------------------------------------

#ifndef SIMULACION
#define SIMULACION 1

#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>

#include "Componentes.h"
#include "Matrix.h"

using namespace std;

const double ERROR = 1.0E-8;

typedef enum { IZQUIERDA=0, ABAJO=1, ATRAS=2, DERECHA=3, ARRIBA=4, ADELANTE=5 } vecinos;

#define opuesto(i) ((i+3)%6)

class SimulationError
{
 public:
  static const int ERROR_SIZE = 64;
  SimulationError( int num, const char* string )
    {
      _errNo = num;
      strncpy( _errString, string, ERROR_SIZE );
    }
  const char* errString(void);
  int errNo(void);
 private:

  int _errNo;
  char _errString[ERROR_SIZE];
};


class VoxelSpace
{
 private:
  Celda* _voxels;
  int    _I, _J, _K;
  int    _numVoxels;
  int    _vecino[6];
  double _dx;

 public:
  class Cursor
    {
    private:
      VoxelSpace* _vSpace;
      Celda*      _celda;
      int         _index;
      int         _i, _j, _k;
    public:
      void init(VoxelSpace& vs, int i, int j, int k) 
	{ _vSpace = &vs; _i = i; _j = j; _k = k; 
	_index = _vSpace->index(i,j,k); _celda = _vSpace->_voxels+_index; }

      void init(VoxelSpace& vs, int index) 
	{ _vSpace = &vs; _index = index; _celda = _vSpace->_voxels+_index; 
	_i = index / (vs._J*vs._K); index %= (vs._J*vs._K); _j = index / vs._K; _k = index % vs._K;}

      Celda*     operator->() { return _celda; }
      operator   Celda*() { return _celda; }
      Cursor&    operator[](int p);
      Cursor&    operator++();
      bool       finished(){ return _index == _vSpace->numVoxels(); }
      Velocidad& u(Vector3D& p);
      Vector3D&  gradientPressure();
      double     divergenceVelocity();
      double bestGuess(int coord, double dt, Vector3D& f);
      double get(int i, double guess);
      int    i() { return _i; }
      int    j() { return _j; }
      int    k() { return _k; }
    };

  inline int index(int i, int j, int k) { return i*(_J*_K)+j*_K+k; }
  inline int numVoxels() { return _numVoxels; }
  inline int I() { return _I; }
  inline int J() { return _J; }
  inline int K() { return _K; }
  inline double dx() { return _dx; }
  VoxelSpace() : _voxels(NULL), _I(0), _J(0), _K(0) {}
  ~VoxelSpace() { if (_voxels != NULL) delete [] _voxels; _voxels = NULL; }
  void init(int I, int J, int K, double dx)
    { 
      _I = I; _J = J; _K = K; 
      _numVoxels = I*J*K;
      _dx = dx;
      _voxels    = new Celda[_numVoxels]; 
      _vecino[ARRIBA]    = index(0,1,0)-index(0,0,0);
      _vecino[DERECHA]   = index(1,0,0)-index(0,0,0);
      _vecino[ADELANTE]  = index(0,0,1)-index(0,0,0);
      _vecino[ABAJO]     = -_vecino[ARRIBA];
      _vecino[ATRAS]     = -_vecino[ADELANTE];
      _vecino[IZQUIERDA] = -_vecino[DERECHA];
    }
  Cursor& operator[](int pos);
  Cursor& operator()(int i, int j, int k);
  Cursor& begin();
};

class Particula
{
 public:
  int i, j, k;
  Vector3D pos;
  Vector3D newPos;
  Vector3D u;
  Particula(int _i, int _j, int _k, Vector3D _pos) :
    i(_i), j(_j), k(_k), pos(_pos) {}
  void relocateParticle(VoxelSpace& space);
  void moveParticle(VoxelSpace& space, double dt);
};

class Float
{
 private:
  double* d;
 public:
  Float(double& x) : d(&x) {}
  operator double&() { return *d; }
  Float& operator=(Float& other) { *d = *other.d; return *this; }
  Float& operator=(double x) { *d = x; return *this; }
};

class Simulacion
{
 private:
  char       _nombre[30];       // nombre de la simulacion
  int        _step;             // etapa actual de la simulacion
  double     _time;             // tiempo que queda por simular
  double     _time0;            // tiempo anterior (segundos)
  double     _time1;            // tiempo proximo  (segundos)
  int        _niteraciones;	// Numero de iteraciones que se van a realizar 
  int        _nparticulas;	// Numero total de particulas                          
  int        _nLiquidas;	// Contador que contabiliza la cantidad de celdas con fluidos
  int	     _nCoord[3];
  double     _dt;		// Dt que se va a utilizar en los calculos
  double     _maxDt;    // Dt maximo permisible
  double     _presion;		// Escalar correspondiente a la presion                
  double     _densidad;		// Densidad valor rho en ecuacion (2)    
  double     _viscosidad;	// Densidad valor rho en ecuacion (2)   
  VoxelSpace _voxels;		// Conjunto de voxels
  vector<Particula> _particulas;
  Vector3D   _condForma;	// Condiciones de forma de 3.21
  double     findDt();
  void       relocateParticles();
  void       initNewCells();
  Vector3D&  getClosestU(Vector3D x);
  void       solidBoundaryConditions();
  void       emptyBoundaryConditions();
  void       emptyBoundaryConditions(
				     Float izquierda, Float abajo, Float atras,
				     Float derecha, Float arriba, Float adelante,
				     Estado eizquierda, Estado eabajo, Estado eatras,
				     Estado ederecha, Estado earriba, Estado eadelante,
				     bool& borde
				     );
  void       calculateNewVelocitiesStep1();
  void       recountLiquidCells();
  void       solvePressures();
  void	     solveVelocity(int dir);
  void       savePovRay();
  void       saveVelocity();

 public:
  Simulacion(int step, char nombre[]);
  double width() { return _voxels.I() * _voxels.dx(); }
  void   read(char* in_file = NULL);
  void   write(char* out_file = NULL);
  void   simulate(int iteraciones, double tiempo);
  void   fillRandomParticles(int cant);
  void   fillParticles(int cant1);
  void   saveParticles();
  void   loadParticles();
};

#endif
