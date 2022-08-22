/*

   Compute velocity autocorrelation function from an MD simulation
                    Anze Hubman, Aug. 2022
*/


/*

  For compilation use g++ -O3 -mcmodel=medium

*/


#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
using namespace std;

// global
const int n_atom   = 576;                   // number of atoms
const int n_step   = 500000;                // number of MD steps
const int n_samp   = 10;                    // sampling frequency
const int n_frames = (n_step/n_samp)+1;     // number of frames, LAMMPS also stores frame at t=0
string trj_file    = "velocity.dat";        // trajectory file (format specification to be added)

float vx[n_frames][n_atom];
float vy[n_frames][n_atom];
float vz[n_frames][n_atom];

// main
int main()
{
  int ii, col, dm2, row;
  string dm1;
  float dmx, dmy, dmz, vx_k, vy_k, vz_k, average_origin, average_atoms;
  float vacf[n_frames];
  

  // read trajectory and rearrange velocities
  ii  = 0;
  col = 0;
  fstream infile(trj_file);
  while (infile >> dm1 >> dm2 >> dmx >> dmy >> dmz >> vx_k >> vy_k >> vz_k)
    {
      row = ii/n_atom;
      vx[row][col] = vx_k;
      vy[row][col] = vy_k;
      vz[row][col] = vz_k;
      ii += 1;
      col += 1;
      if (col == n_atom) col = 0;
    }
  cout << "Trajectory read." << endl;
  infile.close();

  // compute VACF; averaging over all atoms and over specified time origins
  for (int tau = 0; tau < n_frames; tau++)
    {

      if ((tau % 100) == 0) cout << "Current tau: " << n_samp*tau << endl;
      
      average_origin = 0.0;
      for (int i = 0; i < n_frames-tau; i++)
	{
	  average_atoms = 0.0;
	  for (int j = 0; j < n_atom; j++)
	    {
	      average_atoms += vx[i][j]*vx[i+tau][j] + vy[i][j]*vy[i+tau][j] + vz[i][j]*vz[i+tau][j];
	    }

	  average_atoms = average_atoms/n_atom;

	  average_origin += average_atoms;
	}
      average_origin = average_origin/(n_frames-tau);
      vacf[tau] = average_origin;
    }

  // write VACF to file
  ofstream outfile("vacf.dat");

  for (int i = 0; i < n_frames; i++)
  {
    outfile << i*n_samp << " " << vacf[i]/vacf[0] << "\n";   // normalised
  }
  outfile.close();

  return 0;
}

