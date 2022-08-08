/*

   Compute pair correlation function for arbitrary box shapes under 
   periodic boundary conditions.
                        July 2022, Slovenia

*/

/*
  Author:           Anze Hubman
                    National Institute of Chemistry/Theory department

  Disclaimer:       This code can be used and modified freely.
                    I am not taking responsibility for any mistakes that might occur.
                    For details contact me e.g. on LinkedIn.

  Code description:

                    A C++ code to compute g(r) for an arbitrary simulation box shape.
                    Primarily intended, but not restricted to  LAMMPS. Compile with g++ -O3.
  
                    For LAMMPS:

                    - to create 'shape.dat':
                      1. add thermo_style  custom atoms step lx ly lz xy xz yz
                      2. from log.script_name extract atoms step lx ly lz xy xz yz to a
                         continuous file using e.g. grep

                    - to create 'shape_i.dat' and 'shape_j.dat'
                      1. add dump        trj all custom n_samp trajectory_file element id x y z
                                         dump_modify trj sort id element name_of_i name_of_j name_of_k ...
                      2. grep 'name_of_i' trajectory > coor_i.dat   (repeat for all desired types)

                    General remark: to compute g(r) between i and i use coor_i.dat for i and j case.
*/


#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
using namespace std;

// global variables
const int n_i      = 96;                   // number of atoms of type i
const int n_j      = 96;                   // number of atoms of type j
const int n_step   = 500000;               // number of MD steps
const int n_samp   = 1000;                 // sampling frequency
const int n_frames = (n_step/n_samp)+1;    // LAMMPS writes frame at t=0
string coor_i      = "coor_i.dat";         // trajectory of type i
string coor_j      = "coor_j.dat";         // trajectroy of type j
string shape       = "shape.dat";          // h matrix elements for each t
const int n_bid    = 100;                  // discretization of r
const float r_cut  = 12.0;                 // cutoff radius
const float pi     = 4.0*atan(1.0);

float xi[n_frames][n_i];
float yi[n_frames][n_i];
float zi[n_frames][n_i];
float xj[n_frames][n_j];
float yj[n_frames][n_j];
float zj[n_frames][n_j];
float s[n_frames][6];
float g[n_frames][n_bid];

// main loop
int main()
{
  string dummy1, fname;
  int    dummy2, dummy3, row, ii, col, bid;
  float  h[3][3], hinv[3][3], deth, invdeth, V0, V1;
  float  xk, yk, zk, lxk, lyk, lzk, xyk, xzk, yzk, delta, volume;
  float  dsx, dsy, dsz, rx, ry, rz, r, ri, g_smooth[n_bid];

  // read and rearrange trajectory
  ii  = 0;
  col = 0;
  fstream infile_i(coor_i);
  while (infile_i >> dummy1 >> dummy2 >> xk >> yk >> zk)
    {
      row = ii/n_i;
      xi[row][col] = xk;
      yi[row][col] = yk;
      zi[row][col] = zk;
      ii  += 1;
      col += 1;
      if (col == n_i) col = 0;
    }
  cout << "Coordinates of i read." << endl;
  infile_i.close();

  ii  = 0;
  col = 0;
  fstream infile_j(coor_j);
  while (infile_j >> dummy1 >> dummy2 >> xk >> yk >> zk)
    {
      row = ii/n_j;
      xj[row][col] = xk;
      yj[row][col] = yk;
      zj[row][col] = zk;
      ii  += 1;
      col += 1;
      if (col == n_j) col = 0;
    }
  cout << "Coordinates of j read." << endl;

  // read cell shape
  row = 0;
  fstream infile_s(shape);
  while (infile_s >> dummy2 >> dummy3 >> lxk >> lyk >> lzk >> xyk >> xzk >> yzk)
    {
      s[row][0] = lxk;
      s[row][1] = lyk;
      s[row][2] = lzk;
      s[row][3] = xyk;
      s[row][4] = xzk;
      s[row][5] = yzk;
      row += 1;
    }
  cout << "Cell shape read." << endl;
  infile_s.close();

  // initialise g(r) for each frame
  for (int i = 0; i < n_frames; i++) {
    for (int j = 0; j < n_bid; j++) {
      g[i][j] = 0.0;
    }
  }

  // loop over frames
  delta = r_cut/n_bid;
  cout << "Computing g(r)..." << endl;

  for (int k = 0; k < n_frames; k++)
    {

      if ((k % 500) == 0) cout << "Current frame: " << k << endl;

      // generate h matrix at k-th step
      for (int i = 0; i < 3; i++) {
	for (int j = 0; j < 3; j++) {
	  h[i][j] = 0.0;
	}
      }

      h[0][0] = s[k][0];
      h[0][1] = s[k][3];
      h[0][2] = s[k][4];
      h[1][1] = s[k][1];
      h[1][2] = s[k][5];
      h[2][2] = s[k][2];

      // determinant of h
      deth = h[0][0] * (h[1][1] * h[2][2] - h[1][2] * h[2][1])
	- h[0][1] * (h[1][0] * h[2][2] - h[1][2] * h[2][0])
	+ h[0][2] * (h[1][0] * h[2][1] - h[1][1] * h[2][0]);

      // inverse of determinant (V = deth)
      invdeth = 1.0/deth;
      volume  = deth;

      // inverse of h
      hinv[0][0] = (h[1][1] * h[2][2] - h[2][1] * h[1][2]) * invdeth;
      hinv[0][1] = (h[0][2] * h[2][1] - h[0][1] * h[2][2]) * invdeth;
      hinv[0][2] = (h[0][1] * h[1][2] - h[0][2] * h[1][1]) * invdeth;
      hinv[1][0] = (h[1][2] * h[2][0] - h[1][0] * h[2][2]) * invdeth;
      hinv[1][1] = (h[0][0] * h[2][2] - h[0][2] * h[2][0]) * invdeth;
      hinv[1][2] = (h[1][0] * h[0][2] - h[0][0] * h[1][2]) * invdeth;
      hinv[2][0] = (h[1][0] * h[2][1] - h[2][0] * h[1][1]) * invdeth;
      hinv[2][1] = (h[2][0] * h[0][1] - h[0][0] * h[2][1]) * invdeth;
      hinv[2][2] = (h[0][0] * h[1][1] - h[1][0] * h[0][1]) * invdeth;

      // compute distances between atoms (i,j) at k-th step
      for (int i = 0; i < n_i; i++)  {
	for (int j = 0; j < n_j; j++) {

	  dsx = hinv[0][0]*(xi[k][i] - xj[k][j]) + hinv[0][1]*(yi[k][i] - yj[k][j]) + hinv[0][2]*(zi[k][i] - zj[k][j]);
	  dsy = hinv[1][0]*(xi[k][i] - xj[k][j]) + hinv[1][1]*(yi[k][i] - yj[k][j]) + hinv[1][2]*(zi[k][i] - zj[k][j]);
	  dsz = hinv[2][0]*(xi[k][i] - xj[k][j]) + hinv[2][1]*(yi[k][i] - yj[k][j]) + hinv[2][2]*(zi[k][i] - zj[k][j]);

	  dsx = dsx - nearbyint(dsx);
	  dsy = dsy - nearbyint(dsy);
	  dsz = dsz - nearbyint(dsz);

	  rx = h[0][0]*dsx + h[0][1]*dsy + h[0][2]*dsz;
	  ry = h[1][0]*dsx + h[1][1]*dsy + h[1][2]*dsz;
	  rz = h[2][0]*dsx + h[2][1]*dsy + h[2][2]*dsz;

	  r  = sqrt(rx*rx + ry*ry + rz*rz);   // this is valid for orthorombic only, to be changed !!!

	  //	  cout << i << " " << j << " " << r << endl; 

	  if ((r >= 0.001) && (r < r_cut))
	    {
	      bid = int(r/delta);
	      g[k][bid] += 1.0;
	    }
	}
      }

      // normalise k-th g(r)
      for (int i = 0; i < n_bid; i++)
	{
	  ri = (i+0.5)*delta;
	  V0 = (4.0/3.0)*pi*(i*delta)*(i*delta)*(i*delta);
	  V1 = (4.0/3.0)*pi*((i+1)*(i+1)*(i+1)*delta*delta*delta);
	  g[k][i] = g[k][i]*volume/(n_i*n_j*(V1-V0));
	}

      // output g(r)
      if ((k % 500) == 0) {
	ofstream outfile("gr_ij_"+to_string(k));
	for (int j = 0; j < n_bid; j++)
	  {
	    outfile << (j+0.5)*delta << " " << g[k][j] << "\n";
	  }
	outfile.close();
      }

    }

  return 0;
}
