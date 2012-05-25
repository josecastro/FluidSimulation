//-----------------------------------------------------------------------------
// Componentes.cpp
//
// Project: Tephra Irazu simulation - (C) CIC-ITCR
// Author : José Castro
// Date   : July 15, 2005
//
// Purpose:
//      CPP implementation of data for fluid simulation
//
//-----------------------------------------------------------------------------

#include "Componentes.h"

Vector3D vectorBuffer[VECTOR_BUFFER_SIZE];
int      vectorPos = 0;

Vector3D& Vector3D::operator *(double k) const
{
  Vector3D& v = vectorBuffer[vectorPos];
  vectorPos = (vectorPos+1) % VECTOR_BUFFER_SIZE;

  for (int i = 0; i < 3; ++i)
	  v[i] = coord[i]*k;

  return v;
}
 
Vector3D& Vector3D::operator +(const Vector3D& otro) const
{
  Vector3D& v = vectorBuffer[vectorPos];
  vectorPos = (vectorPos+1) % VECTOR_BUFFER_SIZE;

  for (int i = 0; i < 3; ++i)
	  v[i] = coord[i]+otro[i];

  return v;
}

Vector3D& Vector3D::operator -(const Vector3D& otro) const
{
  Vector3D& v = vectorBuffer[vectorPos];
  vectorPos = (vectorPos+1) % VECTOR_BUFFER_SIZE;

  for (int i = 0; i < 3; ++i)
	  v[i] = coord[i]-otro[i];

  return v;
}

double distancia(Vector3D& p1, Vector3D& p2)
{
  double result = 0.0;

	for (int i = 0; i < 3; ++i)
		result = max(result, abs(p1[i]-p2[i]));
  return result;
}

Vector3D& Vector3D::operator - () const
{
  Vector3D& v = vectorBuffer[vectorPos];
  vectorPos = (vectorPos+1) % VECTOR_BUFFER_SIZE;

	  for (int i = 0; i < 3; ++i)
	      v[i] = -coord[i];

  return v;
}
void Celda::Init(double pr, double vs)
{
  tipo = vacio;
  pressure = pr;
  // Constante
  viscosidad = vs;
  solverIndex[0] = solverIndex[1] = solverIndex[2] = solverIndex[3] = 0;
  newCelda = false;
  borde    = false;
}

Celda::Celda(double pr, double vs) { Init(pr, vs); }

Celda::Celda(Velocidad& v, Estado e, double pr, double vs)
{
  u = v;
  tipo = e;
  pressure = pr;
  // Constante
  viscosidad = vs;
  solverIndex[0] = solverIndex[1] = solverIndex[2] = solverIndex[3] = 0;
  newCelda = false;
}
