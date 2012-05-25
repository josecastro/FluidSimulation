//-----------------------------------------------------------------------------
// Particula.cpp
//
// Project: Tephra Irazu simulation - (C) CIC-ITCR
// Author : Santiago Nu\ez
// Date   : July 15, 2005
//
// Purpose:
//      CPP implementaion of particle movement for fluid simulation using MAC
//      Taken from Carlson Lewis CS Georgia Tech Ph.D. dissertation, 2004
//
//-----------------------------------------------------------------------------

void Particula::moveParticle(VoxelSpace& space, double dt)
{
  int i0 = i;
  int j0 = j;
  int k0 = k;
  VoxelSpace::Cursor celda, newCelda;
	
  do
    {
      i = i0; j = j0; k = k0;
      celda = space(i,j,k);

      u = celda.u(pos);
      newPos = pos + u*dt;
      relocateParticle(space);
	
      // i,j,k pueden haber cambiado
      newCelda = space(i,j,k);
	
      if (newCelda->tipo == liquido)
	{
	  u = (u + newCelda.u(newPos))*0.5;
	  newPos = pos + u*dt;
	  relocateParticle(space);
	}
      newCelda = space(i,j,k);
      switch (newCelda->tipo)
	{
    	case solido  : dt /= 2.0; break;
    	case vacio   : newCelda->newCelda = true; break;
	case liquido : u = newCelda.u(newPos); break;
	}
    }
  while (newCelda->tipo == solido);

  assert(abs(newCelda.i()-celda.i()) <= 1);
  assert(abs(newCelda.j()-celda.j()) <= 1);
  assert(abs(newCelda.k()-celda.k()) <= 1);
  assert(newCelda.i() > 0);
  assert(newCelda.j() > 0);
  assert(newCelda.k() > 0);
  assert(newCelda.i() < space.I()-1);
  assert(newCelda.j() < space.J()-1);
  assert(newCelda.k() < space.K()-1);

  pos = newPos;
}

void Particula::relocateParticle(VoxelSpace& space)
{
  double dx   = space.dx();
  double xpos = (i-0.5)*dx;
  double ypos = (j-0.5)*dx;
  double zpos = (k-0.5)*dx;

  while (newPos[X] > xpos+dx) { ++i; xpos += dx; }
  while (newPos[X] < xpos)    { --i; xpos -= dx; }

  while (newPos[Y] > ypos+dx) { ++j; ypos += dx; }
  while (newPos[Y] < ypos)    { --j; ypos -= dx; }

  while (newPos[Z] > zpos+dx) { ++k; zpos += dx; }
  while (newPos[Z] < zpos)    { --k; zpos -= dx; }
}
