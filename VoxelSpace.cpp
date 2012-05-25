//-----------------------------------------------------------------------------
// VoxelSpace.cpp
//
// Project: Tephra Irazu simulation - (C) CIC-ITCR
// Author : José Castro
// Date   : July 15, 2005
//
// Purpose:
//      CPP implementaion of voxel space for fluid simulation
//      Taken from Carlson Lewis CS Georgia Tech Ph.D. dissertation, 2004
//
//-----------------------------------------------------------------------------

#define CURSOR_BUFFER_SIZE 35

static  VoxelSpace::Cursor cursors[CURSOR_BUFFER_SIZE];
int     cursor_pos = 0;

extern Vector3D vectorBuffer[];
extern int vectorPos;

//-----------------------------------------------------------------------------
//
// Class Cursor (pertenece a VoxelSpace)
//
//  Utilizada para accesar las Celdas de el espacio
//  el Cursor tiene conocimiento de sus vecinos
//
//  Dado un cursor c se puede accesar
//  c[ARRIBA], c[ABAJO] ...
//  con el operador -> se obtiene un puntero a la Celda
//-----------------------------------------------------------------------------

VoxelSpace::Cursor& VoxelSpace::Cursor::operator[](int p)
{
  int i = cursor_pos;
	
  assert(_index+_vSpace->_vecino[p] >= 0);
  if (_index+_vSpace->_vecino[p] >= _vSpace->numVoxels())
    {
      cout << "i,j,k = " << _i << ' ' << _j << ' ' << _k << '\n';
      cout << " p = " << p << '\n';
      cout << " vecino[p] = " << _vSpace->_vecino[p] << '\n';
    }
  assert(_index+_vSpace->_vecino[p] < _vSpace->numVoxels());
	
  cursor_pos = (cursor_pos + 1) % CURSOR_BUFFER_SIZE;
  cursors[i].init(*_vSpace, _index+_vSpace->_vecino[p]);

  return cursors[i];
}

//---------------------------------------------------------------
// u(p) obtiene la velocidad de la particula p dentro del espacio
//      utilizando las velocidades en las caras del voxel
//---------------------------------------------------------------
Velocidad& VoxelSpace::Cursor::u(Vector3D& p)
{
  Velocidad& v = vectorBuffer[vectorPos];
  VoxelSpace::Cursor& c = *this;
  vectorPos = (vectorPos + 1) % VECTOR_BUFFER_SIZE;
  double dx = _vSpace->dx();
  Vector3D pos((_i-0.5)*dx, (_j-0.5)*dx, (_k-0.5)*dx);

  for (int i = 0; i < 3; ++i)
    v[i] = (pos[i]+dx-p[i])*(_celda->u[i])+(p[i]-pos[i])*(c[opuesto(i)]->u[i]);

  v *= (1/dx);

  return v;
}


//-------------------------------------------------------
// Diferencial (derivada) de la presion
//   NOTA: Se utilizan condiciones de borde de
//         Dirichlet y de Newmann Pearson
//
//       1.si la celda contigua es solida la derivada
//         de la presion en esa frontera es 0
//       2.si la celda contigua es vacia el diferencial
//         es igual a la presion (presion del vacio es 0)
//         
//-------------------------------------------------------
Vector3D& VoxelSpace::Cursor::gradientPressure()
{
  Vector3D& dp = vectorBuffer[vectorPos];
  VoxelSpace::Cursor& c = *this;

  vectorPos = (vectorPos+1) % VECTOR_BUFFER_SIZE;
  double p = (*this)->pressure;

  for (int i = 0; i < 3; ++i)
    switch (c[i]->tipo)
      {
      case liquido:
	dp[i] = p - c[i]->pressure; break;
      case solido:
	dp[i] = 0; break;
      case vacio:
	dp[i] = p;
      }

  dp *= (1.0/_vSpace->dx());

  return dp;
}

double VoxelSpace::Cursor::get(int dir, double guess)
{
  if (_celda->tipo != vacio)
    return _celda->u[dir];
  else if ((dir == X) && (_i > 0) && ((*this)[dir]->tipo != vacio))
    return _celda->u[dir];
  else if ((dir == Y) && (_j > 0) && ((*this)[dir]->tipo != vacio))
    return _celda->u[dir];
  else if ((dir == Z) && (_j > 0) && ((*this)[dir]->tipo != vacio))
    return _celda->u[dir];
  else
    return guess; // 0.0;
}

double VoxelSpace::Cursor::bestGuess(int coord, double dt, Vector3D& f)
{
  int coord2 = (coord+1) % 3;
  int coord3 = (coord+2) % 3;
  VoxelSpace::Cursor& c = *this;
  double dx     = _vSpace->dx();
  Vector3D u    = _celda->u;
  Vector3D uant = c[coord]->u;
  Vector3D usig;

  usig[X] = c[DERECHA ]->u[X];
  usig[Y] = c[ARRIBA  ]->u[Y];
  usig[Z] = c[ADELANTE]->u[Z];

  double u90  = c[coord2         ].get(coord,u[coord]);
  double u180 = c[coord3         ].get(coord,u[coord]);
  double u270 = c[opuesto(coord2)].get(coord,u[coord]);
  double u360 = c[opuesto(coord3)].get(coord,u[coord]);
	
  double vdes = c[opuesto(coord2)][coord]->u[coord2];
  double wdes = c[opuesto(coord3)][coord]->u[coord3];
	
  double result = _celda->u[coord] + dt/(4.0*dx)*
    (   (uant[coord]+u[coord])*(uant[coord]+u[coord]) 
	- (u[coord]+usig[coord])*(u[coord]+usig[coord])
	+ (u90 +u[coord])*(uant[coord2]+u[coord2]) 
	- (u[coord]+u270)*(usig[coord2]+vdes)
	+ (u180+u[coord])*(uant[coord3]+u[coord3]) 
	- (u[coord]+u360)*(usig[coord3]+wdes) )
		
    + dt*f[coord];

  return result;
}

double VoxelSpace::Cursor::divergenceVelocity()
{
  VoxelSpace::Cursor& c = *this;
	
  return 
    (c[DERECHA ]->u[X] - c->u[X] + 
     c[ARRIBA  ]->u[Y] - c->u[Y] + 
     c[ADELANTE]->u[Z] - c->u[Z]
     ) /_vSpace->dx();
}

VoxelSpace::Cursor& VoxelSpace::Cursor::operator++()
{
  ++_index;
  ++_k;
  if (_k == _vSpace->_K)
    {
      _k = 0;
      ++_j;
      if (_j == _vSpace->_J)
	{
	  _j = 0;
	  ++_i;
	}
    }
  ++_celda;
  return *this;
}

VoxelSpace::Cursor& VoxelSpace::begin()
{
  Cursor& cursor = cursors[cursor_pos];
  cursor_pos = (cursor_pos + 1) % CURSOR_BUFFER_SIZE;
  cursor.init(*this, 0, 0, 0);

  return cursor; 
}

VoxelSpace::Cursor& VoxelSpace::operator()(int i, int j, int k)
{ 
  Cursor& cursor = cursors[cursor_pos];
  cursor_pos = (cursor_pos + 1) % CURSOR_BUFFER_SIZE;
  cursor.init(*this, i, j, k);

  return cursor; 
}

VoxelSpace::Cursor& VoxelSpace::operator[](int pos)
{
  Cursor& c = cursors[cursor_pos];
  cursor_pos = (cursor_pos + 1) % CURSOR_BUFFER_SIZE;
  c.init(*this, pos);

  return c;
}
