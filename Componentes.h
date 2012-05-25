//-----------------------------------------------------------------------------
// Componentes.h
//
// Project: Tephra Irazu simulation - (C) CIC-ITCR
// Author : José Castro
// Date   : July 15, 2005
//
// Purpose:
//      HPP specification of data for fluid simulation
//
//-----------------------------------------------------------------------------
#ifndef COMPONENTES
#define COMPONENTES 1

#include <vector>
#include <cmath>

using namespace std;

#define VECTOR_BUFFER_SIZE 15
typedef enum {liquido=0, solido=1, vacio=2} Estado;
typedef enum { X = 0, Y = 1, Z = 2, CELDA=3} Coordenada;

class Vector3D
{
 public:
  double coord[3];
  double&   operator [](int i) { return coord[i]; }
  double    operator [](int i) const { return coord[i]; }
  Vector3D( double x=0.0, double y=0.0, double z=0.0 )
	  { coord[X] = x; coord[Y] = y; coord[Z] = z; }
  Vector3D(const Vector3D& otro )
	  { coord[X] = otro[X]; coord[Y] = otro[Y]; coord[Z] = otro[Z]; }

  Vector3D& operator *=(double k) { coord[X] *= k; coord[Y] *= k; coord[Z] *= k; return *this; }
  Vector3D& operator +=(const Vector3D& otro ) { coord[X] += otro[X]; coord[Y] += otro[Y]; coord[Z] += otro[Z]; return *this; }
  Vector3D& operator -=(const Vector3D& otro ) { coord[X] -= otro[X]; coord[Y] -= otro[Y]; coord[Z] -= otro[Z]; return *this; }
  Vector3D& operator * (double k) const;
  Vector3D& operator + (const Vector3D& otro) const;
  Vector3D& operator - (const Vector3D& otro) const;
  Vector3D& operator - () const;
};

double distancia(Vector3D& v1, Vector3D& v2);

typedef Vector3D Velocidad;

class Celda 
{
  /* 
   * El objeto aproximación se requiere para el cálculo del nuevo u, de acuerdo
   * con las ecuaciones 3.28 - 3.31
   */
 public:
  Velocidad u;
  Velocidad unew;
  Estado    tipo;
  bool      newCelda;
  bool      borde;
  double pressure;
  double viscosidad;      // OJO se asume una viscosidad constante (1era Version
  int    solverIndex[4];

  void Init(double pr, double vs);
  Celda() : 
    tipo(vacio), pressure(0.0), viscosidad(0.0),
    u(), newCelda(false) {}
  Celda(double pr, double vs);
  Celda(Velocidad& v, Estado e, double pr, double vs);
};

#endif
