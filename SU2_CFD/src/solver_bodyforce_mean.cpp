/*!
 * \file solver_bodyforce_mean.cpp
 * \brief Main subroutines for solving bodyforce turbomachinery problems.
 * \author E.C.Bunschoten
 * \version 6.2.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */
#include "../include/solver_structure.hpp"
#include <vector>

CBodyForceModelSolver::CBodyForceModelSolver(void) : CSolver() {

}
CBodyForceModelSolver::CBodyForceModelSolver(CGeometry *geometry, CConfig *config, CSolver *fluidsolver, unsigned short iMesh) : CSolver() {
	/*
	BodyForceModel constructor
	All BFM settings defined in the configuration file are stored in the right locations
	*/
	unsigned short iZone = config->GetiZone();
	unsigned short BFM_Formulation = config->GetBFM_Formulation();
	SetBFM_Formulation(BFM_Formulation);	// Storing BFM formulation

	cout << "Initializing body-force model" << endl;
	// Storing rotation vector components
	Rotation_vector[0] = config->GetRotation_Rate_X(iZone);
	Rotation_vector[1] = config->GetRotation_Rate_Y(iZone);
	Rotation_vector[2] = config->GetRotation_Rate_Z(iZone);
	Rotation_rate = config->GetBody_Force_Rotation()* 2 * M_PI/60.0; 	// Storing rotation rate
	nPoint_geom = geometry->GetnPoint();	// Storing geometry cell count

	
	
}

void CBodyForceModelSolver::LoadRestart(CGeometry **geometry, CSolver ***fluidsolver, CConfig *config, int val_iter, bool val_update_geo){
	/*
	This function currently doesn't do anything. Upon restarting, the blade geometry is interpolated onto the mesh nodes and 
	the body-force sources are computed based on the restart solution. 
	*/
	cout <<"Loading restart from body-force model" << endl;
}

void CBodyForceModelSolver::PreprocessBFMParams(CGeometry *geometry, CConfig *config, CSolver *fluidSolver){
	/*
	This function performs all the preprocessing functionalities required for the BFM to run. This includes 
	reading the BFM geometry input file, interpolating the blade geometry to the mesh and computing the 
	cylindrical coordinates and projections of the mesh cells.
	*/
	cout << "Preprocessing function called" << endl;
	// Reading BFM input file
    ReadInputFile(config, geometry);
	SU2_MPI::Barrier(MPI_COMM_WORLD);
	
	// Computing cylindrical projections and coordinates of mesh cells
	ComputeCylProjections(geometry, config, fluidSolver);

	// Interpolating blade geometry to mesh
    InterpolateGeometry(geometry, config, fluidSolver);
	SU2_MPI::Barrier(MPI_COMM_WORLD);

	// In case of Thollets BFM, the metal blockage gradient field is computed
	if(kind_bfm==1){
    	ComputeBlockageGradient(geometry, config, fluidSolver);
		SU2_MPI::Barrier(MPI_COMM_WORLD);
	}

	
}
void CBodyForceModelSolver::ComputeRelVelocity(CGeometry *geometry, CConfig *config, CSolver *fluidsolver){
	/*
	Body-force models require the cylindrical, relative velocity. This function computes the relative velocity in the cylindrical
	coordinate system. This function therefore has to be called each iteration prior to the BFM itself.
	*/
	su2double *Coord, *U_i,*Geometric_Parameters;
	su2double W_ax, W_r, W_th, rotFac;
	for(unsigned long iPoint=0; iPoint<geometry->GetnPoint(); iPoint++){

		// Obtaining node coordinates
		Coord = geometry->node[iPoint]->GetCoord();

		// Obtaining solution flow variables
		U_i = fluidsolver->node[iPoint]->GetSolution();

		// Computing relative velocity components in axial, radial and tangential direction
		W_ax = 0;	// Relative, axial velocity
		W_r = 0;	// Relative, radial velocity
		W_th = 0;	// Relative, tangential velocity

		// Obtaining BFM geometric parameters from node
		Geometric_Parameters = fluidsolver->node[iPoint]->GetBodyForceParameters();
		rotFac = Geometric_Parameters[6];	// Rotation factor, distinguish between rotor and stator

		// Adding to the relative velocity components through dot product of absolute velocity and respective unit vector
		for(int iDim = 0; iDim < nDim; iDim++){
			W_ax += proj_vector_axial[iPoint][iDim] * U_i[iDim+1]/U_i[0];
			W_th += proj_vector_tangential[iPoint][iDim] * U_i[iDim+1]/U_i[0];
			W_r += proj_vector_radial[iPoint][iDim] * U_i[iDim+1]/U_i[0];
		}
		// Induced velocity due to rotation is subtracted from tangential velocity
		W_th -= rotFac * Rotation_rate * cyl_coordinates[iPoint][1];

		// Storing relative velocity vector in class data structure
		relative_vel[iPoint][0] = W_ax;
		relative_vel[iPoint][1] = W_th;
		relative_vel[iPoint][2] = W_r;
	}
}

void CBodyForceModelSolver::ComputeCylProjections(CGeometry *geometry, CConfig *config, CSolver *fluidsolver){
	/*
	This function computes the cylindrical coordinates and projections of the mesh cells based on the rotation 
	vector components defined in the configuration file.
	*/
	cout << "Computing cylindrical unit vectors and coordinates" << endl;
	su2double *Coord;
	su2double rot_dot_x, rot_dot_rot;
	su2double ax, radius;
	su2double axial_vector[nDim], radial_vector[nDim];

	for(unsigned long iPoint=0; iPoint<nPoint_geom; iPoint++){
		// Getting node Cartesian coordinates
		Coord = geometry->node[iPoint]->GetCoord();

		rot_dot_x = 0;	// dot product between axial unit vector and coordinate vector
		rot_dot_rot = 0;	// dot product of axial unit vector
		radius = 0;	// Radial coordinate
		ax = 0;		// Axial coordinate

		// Computing axial dot products
		for(int iDim = 0; iDim < nDim; iDim ++){
			rot_dot_x += Rotation_vector[iDim]*Coord[iDim];
			rot_dot_rot += Rotation_vector[iDim]* Rotation_vector[iDim];
		}
		//
		for(int iDim = 0; iDim < nDim; iDim ++){
			axial_vector[iDim] = (rot_dot_x/rot_dot_rot)* Rotation_vector[iDim]; // Parallel projection of coordinate vector onto rotation axis
			radial_vector[iDim] = Coord[iDim] - axial_vector[iDim];	// Orthogonal projection of coordinate vector onto rotation axis
			ax += axial_vector[iDim]*axial_vector[iDim];	
			radius += radial_vector[iDim]*radial_vector[iDim];
			proj_vector_axial[iPoint][iDim] = Rotation_vector[iDim];	// Appending axial unit vector component
		}

		// Computing absolute values for cylindrical coordinates
		radius = sqrt(radius);	// Computing radial coordinate
		ax = sqrt(ax);	// Computing axial coordinate
		if(rot_dot_x < 0.0){ax = -ax;}	// In case the coordinate parallel projection is negative, the axial coordinate is flipped

		// Storing cylindrical coordinates in class property
		cyl_coordinates[iPoint][0] = ax;
		cyl_coordinates[iPoint][1] = radius;

		// Computation of tangential and radial projection vectors
		int i_left, i_right;
		for(int iDim=0; iDim<nDim; iDim++){

			// Computing tangential unit vector components through cross product of radial and axial vectors
			i_left = (iDim + 1) % nDim;
			i_right = (iDim + 2) % nDim;
			proj_vector_tangential[iPoint][iDim] = (proj_vector_axial[iPoint][i_left]*radial_vector[i_right] 
			- proj_vector_axial[iPoint][i_right]*radial_vector[i_left])/radius;

			//Radial projection unit vector is normalized radial coordinate vector
			proj_vector_radial[iPoint][iDim] = radial_vector[iDim]/radius;
		}
		
			
	}
}
void CBodyForceModelSolver::BFM_Hall(CConfig *config, CGeometry *geometry, CSolver *fluidsolver){
	/*
	This function describes Hall's BFM formulation. It computes the normal body-force magnitude, transforms it to
	a 3D, Cartesian force vector and stores it in the node as a body-force source vector.
	*/
	unsigned long iPoint;
	su2double *U_i, *V_i, *Coord_i, *Geometric_Parameters;
	su2double pi = M_PI, pitch;
	su2double bfFac, b, Nx, Nt, Nr, x_le, rotFac, BF_blades, radius;
	su2double *W;
	su2double WdotN, W_nx, W_nr, W_nth, W_px, W_pth, W_pr, W_p, W_mag;
	su2double delta, F_n;
	su2double F_ax, F_th, F_r, e_source;
	su2double BF_source[nDim+2], F[nDim];


	for(iPoint=0; iPoint<geometry->GetnPoint(); iPoint++){
		// Getting solution variables from the flow solver
		U_i = fluidsolver->node[iPoint]->GetSolution();
		// Getting body-force parameters stored in the node
		Geometric_Parameters = fluidsolver->node[iPoint]->GetBodyForceParameters();
		bfFac = Geometric_Parameters[0]; 	// Body-force factor. Multiplies all body-forces.
		b = Geometric_Parameters[1];			// Metal blockage factor
		Nx = Geometric_Parameters[2]; 		// Camber normal component in axial direction
		Nt = Geometric_Parameters[3]; 			// Camber normal component in tangential direction
		Nr = Geometric_Parameters[4];			// Camber normal component in radial direction
		x_le = Geometric_Parameters[5]; 		// Axial distance from leading edge
		rotFac = Geometric_Parameters[6];	// Rotation factor. Multiplies the body-force rotation value.
		BF_blades = Geometric_Parameters[7];				// Blade row blade count

		// Getting cylindrical, relative velocity vector
		W = relative_vel[iPoint];

		radius = cyl_coordinates[iPoint][1];
		pitch = 2 * pi * radius / BF_blades;	// Computing pitch

		WdotN = W[0] * Nx + W[2] * Nr + W[1] * Nt;		// Dot product of relative velocity and camber normal vector
		W_nx = WdotN * Nx, W_nr = WdotN * Nr, W_nth = WdotN * Nt;		// Relative velocity components normal to the blade
		W_px = W[0] - W_nx, W_pr = W[2] - W_nr, W_pth = W[1] - W_nth;  // Relative velocity components parallel to the blade 
		
		W_p = sqrt(W_px * W_px + W_pr * W_pr + W_pth * W_pth);	// Parallel relative velocity magnitude
		// Relative velocity magnitude
		W_mag = 0;
		for(unsigned short iDim=0; iDim<nDim; iDim++){
			W_mag += W[iDim]*W[iDim];
		}
		W_mag = sqrt(W_mag);
		// Calculating the deviation angle
		delta = asin(WdotN / W_mag);

		// Computing normal body force magnitude
		F_n = -bfFac * pi * delta * (1 / pitch) * (1 / abs(Nt)) * W_mag * W_mag;

		// Transforming the normal and force component to cyllindrical coordinates
			
		F_ax = F_n * (cos(delta) * Nx - sin(delta) * (W_px / W_p));		// Axial body-force component
		F_r = F_n * (cos(delta) * Nr - sin(delta) * (W_pr / W_p));		// Radial body-force component
		F_th = F_n * (cos(delta) * Nt - sin(delta) * (W_pth / W_p));	// Tangential body-force component
		e_source = rotFac * Rotation_rate * radius * F_th;				// Energy source term
		
		// Appending Cartesial body-forces to body-force vector
		for(int iDim=0; iDim<nDim; iDim++){
			F[iDim] = U_i[0] * (F_ax*proj_vector_axial[iPoint][iDim] + F_th*proj_vector_tangential[iPoint][iDim] + F_r*proj_vector_radial[iPoint][iDim]);
		}

		// Appending source term values to the body-force source term vector
		BF_source[0] = 0.0;	// empty density component

		// Appending momentum source terms
		for(int iDim = 0; iDim < nDim; iDim ++) {
			BF_source[iDim + 1] = F[iDim];
		}
		// Appending energy source term
		BF_source[nDim+1] = U_i[0] * e_source;

		// Storing body-force source term at node
		fluidsolver->node[iPoint]->SetBody_Force_Source(BF_source);
	}
}
void CBodyForceModelSolver::BFM_Thollet(CConfig *config, CGeometry *geometry, CSolver *fluidsolver){
	/*
	This function computes the body-force source term according to Thollet's BFM formulation
	*/
    /*--- Load all relevant values from config file ---*/
    unsigned long iPoint;
    su2double *U_i, *V_i, *Coord_i, *Geometric_Parameters;
    su2double gamma = config->GetGamma(), R_gas = config->GetGas_Constant(), BF_radius = config->GetBody_Force_Radius();

    /*--- Initialize common variables ---*/
    su2double pi = M_PI, pitch;
	su2double radius, ax;
	su2double W_mag;
	su2double *W;
	su2double bfFac, b, Nx, Nt, Nr, x_le, rotFac, BF_blades;
	su2double WdotN, W_nx, W_nth, W_nr, W_px, W_pth, W_pr, W_p;
	su2double delta, V_sound, M_rel;
	su2double F_n_inc, F_n, F_p;
	su2double F_ax, F_r, F_th, F_x, F_y, F_z, e_source, F[nDim];
	su2double C_f, Re_x, mu=1.716E-5;
	su2double Kprime, K=1;
	su2double BF_source[nDim+2];

    for ( iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++ ) {
		
		// Extracting flow variables and primitive variables.
        U_i = fluidsolver->node[iPoint]->GetSolution();
        V_i = fluidsolver->node[iPoint]->GetPrimitive();
		
		// Extracting node coordinates 
        Coord_i = geometry->node[iPoint]->GetCoord();
			
		// Getting interpolated geometric parameters from node
		Geometric_Parameters = fluidsolver->node[iPoint]->GetBodyForceParameters();

		bfFac = Geometric_Parameters[0]; 	// Body-force factor. Multiplies all body-forces.
		b = Geometric_Parameters[1];			// Metal blockage factor
		Nx = Geometric_Parameters[2]; 		// Camber normal component in axial direction
		Nt = Geometric_Parameters[3]; 			// Camber normal component in tangential direction
		Nr = Geometric_Parameters[4];			// Camber normal component in radial direction
		x_le = Geometric_Parameters[5]; 		// Axial distance from leading edge
		rotFac = Geometric_Parameters[6];	// Rotation factor. Multiplies the body-force rotation value.
		BF_blades = Geometric_Parameters[7];				// Blade row blade count
		
		// Getting relative velocity vector
		W = relative_vel[iPoint];

		// Getting cylindrical coordinates
		ax = cyl_coordinates[iPoint][0];
		radius = cyl_coordinates[iPoint][1];

		WdotN = W[0] * Nx + W[2] * Nr + W[1] * Nt;		// Dot product of relative velocity and camber normal vector
		W_nx = WdotN * Nx, W_nr = WdotN * Nr, W_nth = WdotN * Nt;		// Relative velocity components normal to the blade
		W_px = W[0] - W_nx, W_pr = W[2] - W_nr, W_pth = W[1] - W_nth;  // Relative velocity components parallel to the blade 
		
		W_p = sqrt(W_px * W_px + W_pr * W_pr + W_pth * W_pth);	// Parallel relative velocity magnitude
		W_mag = sqrt(W[0] * W[0] + W[1] * W[1] + W[2] * W[2]);					// Relative velocity magnitude
		
		// Calculating the deviation angle
		delta = asin(WdotN / W_mag);
		// Calculating te local, relative Mach number
		V_sound = sqrt(gamma * R_gas * V_i[0]);
		M_rel = W_mag / V_sound;
		
		// Blade pitch 
		pitch = 2 * pi * radius / BF_blades;
	
		
		// Incompressible normal force magnitude 
		F_n_inc = pi * delta * (1 / pitch) * (1 / abs(Nt)) * (1 / b) * W_mag * W_mag;
		
		
		// Axial Reynolds number, based on leading edge distance.
		Re_x = (abs(ax - x_le) * W_mag * U_i[0]) / mu;
		if(Re_x == 0.0){
			Re_x = (0.001 * W_mag * U_i[0]) / mu;
		}
		// Parallel force friction factor.
		C_f = 0.0592 * pow(Re_x, -0.2) ;

		// Calculating the normal force compressibility factor
		if (M_rel < 1) {
			Kprime = 1 / (sqrt(1 - (M_rel * M_rel)));
			if (Kprime <= 3) {
				K = Kprime;
			}
			if (Kprime > 3) {
				K = 3;
			}
		}
		if (M_rel > 1) {
			Kprime = 2 / (pi * sqrt((M_rel * M_rel) - 1));
			if (Kprime <= 3) {
				K = Kprime;
			}
			if (Kprime > 3) {
				K = 3;
			}
		}
	
		
		// Calculating the final values for the normal and parallel body-force
		F_n = -bfFac * K * F_n_inc;
		F_p = -bfFac * C_f * W_mag * W_mag * (1 / pitch) * (1 / abs(Nt)) * (1 / b);
		
		
		// Transforming the normal and parallel force components to cyllindrical coordinates
		
		F_ax = F_n * (cos(delta) * Nx - sin(delta) * (W_px / W_p)) + F_p * W[0] / W_mag;
		F_r = F_n * (cos(delta) * Nr - sin(delta) * (W_pr / W_p)) + F_p * W[2] / W_mag;
		F_th = F_n * (cos(delta) * Nt - sin(delta) * (W_pth / W_p)) + F_p * W[1] / W_mag;
		e_source = rotFac * Rotation_rate * radius * F_th;
		
		// Appending Cartesial body-forces to body-force vector
		for(int iDim=0; iDim<nDim; iDim++){
			F[iDim] = U_i[0] * (F_ax*proj_vector_axial[iPoint][iDim] + F_th*proj_vector_tangential[iPoint][iDim] + F_r*proj_vector_radial[iPoint][iDim]);
		}

		// Appending source term values to the body-force source term vector
		BF_source[0] = 0.0;	// empty density component

		// Appending momentum source terms
		for(int iDim = 0; iDim < nDim; iDim ++) {
			BF_source[iDim + 1] = F[iDim];
		}

		// Appending energy source term
		BF_source[nDim+1] = U_i[0] * e_source;

		// Storing BFM source term in node
		fluidsolver->node[iPoint]->SetBody_Force_Source(BF_source);
	};
	
};
void CBodyForceModelSolver::AllocateMemory(CGeometry *geometry){
	/*
	This function allocates memory for all the relevant BFM properties. Mainly for the cylindrical projection vectors and interpolation data
	*/
	unsigned long n_points = geometry->GetnPoint();
	unsigned short nDim = geometry->GetnDim();
	
	// Allocating memory for cylindrical projection vectors and coordinates
	proj_vector_axial  = new su2double*[n_points];
	proj_vector_tangential  = new su2double*[n_points];
	proj_vector_radial  = new su2double*[n_points];
	relative_vel = new su2double*[n_points];
	cyl_coordinates = new su2double*[n_points];
	for(unsigned long iPoint=0; iPoint<n_points; iPoint++){
		proj_vector_axial[iPoint] = new su2double[nDim];
		proj_vector_tangential[iPoint] = new su2double[nDim];
		proj_vector_radial[iPoint] = new su2double[nDim];
		relative_vel[iPoint] = new su2double[nDim];
		cyl_coordinates[iPoint] = new su2double[2];
	}

	// Allocating memory for the BFM interpolation data
	row_blade_count = new su2double[n_rows];
	row_rotation = new su2double[n_rows];

	ax_coordinate = new su2double**[n_rows];
	rad_coordinate = new su2double**[n_rows];
	camber_normal_data = new su2double***[n_rows];
	blockage_data = new su2double**[n_rows];
	X_le_data = new su2double**[n_rows];
	chord_data = new su2double**[n_rows];

	for(unsigned long i_p=0; i_p<n_rows; i_p++){
		ax_coordinate[i_p] = new su2double*[n_sec];
		rad_coordinate[i_p] = new su2double*[n_sec];
		camber_normal_data[i_p] = new su2double**[n_sec];
		blockage_data[i_p] = new su2double*[n_sec];
		X_le_data[i_p] = new su2double*[n_sec];
		chord_data[i_p] = new su2double*[n_sec];
		for(unsigned long i_sec=0; i_sec<n_sec; i_sec++){
			ax_coordinate[i_p][i_sec] = new su2double[n_axial];
			rad_coordinate[i_p][i_sec] = new su2double[n_axial];
			camber_normal_data[i_p][i_sec] = new su2double*[n_axial];
			blockage_data[i_p][i_sec] = new su2double[n_axial];
			X_le_data[i_p][i_sec] = new su2double[n_axial];
			chord_data[i_p][i_sec] = new su2double[n_axial];
			for(unsigned long i_ax=0; i_ax<n_axial; i_ax++){
				camber_normal_data[i_p][i_sec][i_ax] = new su2double[3];
			}
		}
	}
}
void CBodyForceModelSolver::ReadInputFile(CConfig *config, CGeometry *geometry){
	/*
	This function reads the BFM input file and stores the interpolation data in the respective data structures accordingly. 
	*/
    cout << "Reading Body Force Model geometry input file" <<endl;

	// Problem dimension count
    nDim = geometry->GetnDim();
    
	string line;
	// Opening body-force model input file
	std::ifstream inputFile((config->GetBFM_inputName()).c_str());
	
	// Reading number of blade rows, number of radial blade sections and axial points from first line
	inputFile >> n_rows >> n_sec >> n_axial;

	// Allocating memory based on interpolation data size
	AllocateMemory(geometry);

	su2double axial_coordinate, radial_coordinate, N_axial, N_tangential, N_radial,
	blockage_factor, axial_c_LE, axial_chord, rotation_factor, blade_count;

	// Looping over blade rows, radial sections and axial points to store interpolation data
    for (int i_row=0; i_row < n_rows; i_row++){
		for (int i_sec=0; i_sec < n_sec; i_sec++){
			// Skipping radial section header line
			getline(inputFile, line);
			for (int i_ax=0; i_ax<n_axial; i_ax++){
				getline(inputFile, line);
				// Storing values read on the line
				inputFile >> ax_coordinate[i_row][i_sec][i_ax] 
				>> rad_coordinate[i_row][i_sec][i_ax] 
				>> camber_normal_data[i_row][i_sec][i_ax][0]
				>> camber_normal_data[i_row][i_sec][i_ax][1]
				>> camber_normal_data[i_row][i_sec][i_ax][2]
				>> blockage_data[i_row][i_sec][i_ax]
				>> X_le_data[i_row][i_sec][i_ax]
				>> chord_data[i_row][i_sec][i_ax]
				>> rotation_factor 
				>> blade_count;
				
				if(axial_coordinate <= x_min){
					x_min = axial_coordinate;
				}
					
				
				
			}
		}
		row_rotation[i_row] = rotation_factor;	// Blade row rotation factor
		row_blade_count[i_row] = blade_count;	// Blade row blade count

	}
	inputFile.close();
}
void CBodyForceModelSolver::InterpolateGeometry(CGeometry *geometry, CConfig *config, CSolver *fluidSolver){
	/*
	This function interpolates the BFM interpolation data onto the mesh nodes and stores the data locally. 
	*/
	cout << "Interpolating geometric parameters to mesh" << endl;
  	unsigned long iPoint;
  	unsigned short iMesh, iDim;
  	unsigned short iZone = config->GetiZone();
	su2double BF_radius = config->GetBody_Force_Radius();	
	bool RC_interp = true;
	bool S_projection = false;

	// Subtracting 1 meter from minimum axial coordinate to ensure it's located outside of the domain
	x_min = x_min - 1.0;
	
	// Initializing values for the BFM geometric parameters
	su2double b=1.0, Nx = 0.0, Nt = 1.0, Nr = 0.0, bfFac = 0.0, x_le = 0.0, rotFac = 0.0, bladeCount = 1, chord = 1.0;
	// BFM parameter array
	su2double BodyForceParams[9] = {bfFac, b, Nx, Nt, Nr, x_le, rotFac, bladeCount, chord};
	
	su2double x1, x2, r1, r2;
	int nInt=0;
	bool inside = false;
	su2double dist = 0, deNom=0;
	su2double eNum_b = 0, eNum_Nx = 0, eNum_Nt = 0, eNum_Nr = 0, eNum_x_le = 0, eNum_chord = 0;
	su2double cell_axial_coordinate, cell_radial_coordinate;
	su2double rot_dot_x, rot_dot_rot;
	su2double axial[nDim], radial[nDim];
	//su2double x_side [4][2] {{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}};
	vector <vector<su2double>> x_side {{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}};
	vector<su2double> blockage_side, Nx_side, Nt_side, Nr_side, x_le_side, chord_side;
	vector<su2double> edge {0, 0}, v_x {0, 0}, v_opp {0, 0}, proj_x {0, 0}, proj_opp {0, 0};
	su2double orientation;
	// Looping over points in the domain to interpolate body-force parameters onto the respective cells
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++){

		// Getting cylindrical coordinates for current cell
		cell_axial_coordinate = cyl_coordinates[iPoint][0];
		cell_radial_coordinate = cyl_coordinates[iPoint][1];

		// Setting default BodyForceParams vector with farfield values
		b = 1.0;
		Nx = 0.0;
		Nt = 1.0;
		Nr = 0.0;
		bfFac = 0.0;
		x_le = 0.0;
		rotFac = 0.0;
		bladeCount = 1;
		chord = 1.0;
		BodyForceParams[0] = bfFac;
		BodyForceParams[1] = b;
		BodyForceParams[2] = Nx;
		BodyForceParams[3] = Nt;
		BodyForceParams[4] = Nr;
		BodyForceParams[5] = x_le;
		BodyForceParams[6] = rotFac;
		BodyForceParams[7] = bladeCount;
		BodyForceParams[8] = chord;
		
		// Looping over blade rows
		for (int q=0; q < n_rows; q++){
			// Looping over radial sections
			for(int j=0; j < n_sec-1; j++){
				// Looping over axial points
				for(int n = 0; n < n_axial - 1; n++){
					// Appending interpolation values to sides loop
					x_side[0][0] = ax_coordinate[q][j][n];
					x_side[0][1] = rad_coordinate[q][j][n];
					x_side[1][0] = ax_coordinate[q][j+1][n];
					x_side[1][1] = rad_coordinate[q][j+1][n];
					x_side[2][0] = ax_coordinate[q][j+1][n+1];
					x_side[2][1] = rad_coordinate[q][j+1][n+1];
					x_side[3][0] = ax_coordinate[q][j][n+1];
					x_side[3][1] = rad_coordinate[q][j][n+1];
		
					// Interpolation data along interpolation cell
					blockage_side = {blockage_data[q][j][n], blockage_data[q][j+1][n], blockage_data[q][j+1][n+1], blockage_data[q][j][n+1]};
					Nx_side = {camber_normal_data[q][j][n][0], camber_normal_data[q][j+1][n][0], camber_normal_data[q][j+1][n+1][0], camber_normal_data[q][j][n+1][0]};
					Nt_side = {camber_normal_data[q][j][n][1], camber_normal_data[q][j+1][n][1], camber_normal_data[q][j+1][n+1][1], camber_normal_data[q][j][n+1][1]};
					Nr_side = {camber_normal_data[q][j][n][2], camber_normal_data[q][j+1][n][2], camber_normal_data[q][j+1][n+1][2], camber_normal_data[q][j][n+1][2]};
					x_le_side = {X_le_data[q][j][n], X_le_data[q][j+1][n], X_le_data[q][j+1][n+1], X_le_data[q][j][n+1]};
					chord_side = {chord_data[q][j][n], chord_data[q][j+1][n], chord_data[q][j+1][n+1], chord_data[q][j][n+1]};

					// Checking whether point is inside of the current input node.
					inside = RC_inclusion(cell_axial_coordinate, cell_radial_coordinate, x_side);

					// If the point is inside, the body-force parameters are interpolated through bilinear interpolation from the corner points.
					if (inside){
						
						bfFac = 1.0;
						// b = DW_average(x_side, blockage_side, cell_axial_coordinate, cell_radial_coordinate);
						// Nx = DW_average(x_side, Nx_side, cell_axial_coordinate, cell_radial_coordinate);
						// Nt = DW_average(x_side, Nt_side, cell_axial_coordinate, cell_radial_coordinate);
						// Nr = DW_average(x_side, Nr_side, cell_axial_coordinate, cell_radial_coordinate);
						// x_le = DW_average(x_side, x_le_side, cell_axial_coordinate, cell_radial_coordinate);
						// chord = DW_average(x_side, chord_side, cell_axial_coordinate, cell_radial_coordinate);
						// rotFac = row_rotation[q];
						// bladeCount = row_blade_count[q];

						b = Bilinear_Interp(x_side, blockage_side, cell_axial_coordinate, cell_radial_coordinate);
						Nx = Bilinear_Interp(x_side, Nx_side, cell_axial_coordinate, cell_radial_coordinate);
						Nt = Bilinear_Interp(x_side, Nt_side, cell_axial_coordinate, cell_radial_coordinate);
						Nr = Bilinear_Interp(x_side, Nr_side, cell_axial_coordinate, cell_radial_coordinate);
						x_le = Bilinear_Interp(x_side, x_le_side, cell_axial_coordinate, cell_radial_coordinate);
						chord = Bilinear_Interp(x_side, chord_side, cell_axial_coordinate, cell_radial_coordinate);
						rotFac = row_rotation[q];
						bladeCount = row_blade_count[q];

						// Setting loop variables to their maximum to exit the loop
						q = n_rows;
						j = n_sec;
						n = n_axial;
					}
					//Storing interpolated geometric parameters into vector
					BodyForceParams[0] = bfFac;
					BodyForceParams[1] = b;
					BodyForceParams[2] = Nx;
					BodyForceParams[3] = Nt;
					BodyForceParams[4] = Nr;
					BodyForceParams[5] = x_le;
					BodyForceParams[6] = rotFac;
					BodyForceParams[7] = bladeCount;
					BodyForceParams[8] = chord;
				}
			}
		}
		
		// Storing body-force parameter vector onto node
		fluidSolver->node[iPoint]->SetBodyForceParameters(BodyForceParams);
	}
}

void CBodyForceModelSolver::ComputeBlockageGradient(CGeometry *geometry, CConfig *config, CSolver *fluidSolver){
	/*
	This function computes the spatial metal blockage gradient field, which is used by Thollet's BFM to account for obstruction due to blade thickness
	The gradient method used is similar to the Green Gauss method used for computing the gradient of the primary flow variables
	*/
	cout << "Computing metal blockage factor gradient field" << endl;
	unsigned short iDim, jDim, iNeigh;
	unsigned long iPoint, jPoint;
	su2double Blockage_i, Blockage_j, *BFMVec, *Coord_i, *Coord_j, r11, r12, r13, r22, r23, r23_a,
	r23_b, r33, weight, product, z11, z12, z13, z22, z23, z33, detR2, cvector[nDim], smatrix[nDim][nDim];
	bool singular;
	
	/*--- Loop over points of the grid ---*/
	
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
		
		/*--- Set the value of the singular ---*/
		singular = false;
		
		/*--- Get coordinates ---*/
		
		Coord_i = geometry->node[iPoint]->GetCoord();
		/*--- Get primitives from CVariable ---*/
		BFMVec = fluidSolver->node[iPoint]->GetBodyForceParameters();
		Blockage_i = BFMVec[1];
		/*--- Inizialization of variables ---*/
		
		
		for (iDim = 0; iDim < nDim; iDim++){
			cvector[iDim] = 0.0;
		}
		r11 = 0.0; r12 = 0.0;   r13 = 0.0;    r22 = 0.0;
		r23 = 0.0; r23_a = 0.0; r23_b = 0.0;  r33 = 0.0;
		
		/*
		AD::StartPreacc();
		AD::SetPreaccIn(PrimVar_i, nPrimVarGrad);
		AD::SetPreaccIn(Coord_i, nDim);
		*/
		for (iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnPoint(); iNeigh++) {
		jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
		Coord_j = geometry->node[jPoint]->GetCoord();
		
		BFMVec = fluidSolver->node[jPoint]->GetBodyForceParameters();
		Blockage_j = BFMVec[1];
		
		/*
		AD::SetPreaccIn(Coord_j, nDim);
		AD::SetPreaccIn(PrimVar_j, nPrimVarGrad);
		*/
		
		weight = 0.0;
		for (iDim = 0; iDim < nDim; iDim++)
			weight += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);
		
		/*--- Sumations for entries of upper triangular matrix R ---*/
		if (weight != 0.0) {
			
			r11 += (Coord_j[0]-Coord_i[0])*(Coord_j[0]-Coord_i[0])/weight;
			r12 += (Coord_j[0]-Coord_i[0])*(Coord_j[1]-Coord_i[1])/weight;
			r22 += (Coord_j[1]-Coord_i[1])*(Coord_j[1]-Coord_i[1])/weight;
			
			if (nDim == 3) {
			r13 += (Coord_j[0]-Coord_i[0])*(Coord_j[2]-Coord_i[2])/weight;
			r23_a += (Coord_j[1]-Coord_i[1])*(Coord_j[2]-Coord_i[2])/weight;
			r23_b += (Coord_j[0]-Coord_i[0])*(Coord_j[2]-Coord_i[2])/weight;
			r33 += (Coord_j[2]-Coord_i[2])*(Coord_j[2]-Coord_i[2])/weight;
			}
			
			/*--- Entries of c:= transpose(A)*b ---*/
			
			
			for (iDim = 0; iDim < nDim; iDim++){
				cvector[iDim] += (Coord_j[iDim]-Coord_i[iDim])*(Blockage_j-Blockage_i)/weight;
			}
			
		}
		
		}
		//cout << " Cvector[0]: " << cvector[0] << " Cvector[1]: " << cvector[1] << endl;
		/*--- Entries of upper triangular matrix R ---*/
		
		if (r11 >= 0.0) r11 = sqrt(r11); else r11 = 0.0;
		if (r11 != 0.0) r12 = r12/r11; else r12 = 0.0;
		if (r22-r12*r12 >= 0.0) r22 = sqrt(r22-r12*r12); else r22 = 0.0;
		if (nDim == 3) {
		if (r11 != 0.0) r13 = r13/r11; else r13 = 0.0;
		if ((r22 != 0.0) && (r11*r22 != 0.0)) r23 = r23_a/r22 - r23_b*r12/(r11*r22); else r23 = 0.0;
		if (r33-r23*r23-r13*r13 >= 0.0) r33 = sqrt(r33-r23*r23-r13*r13); else r33 = 0.0;
		}
		
		/*--- Compute determinant ---*/
		
		if (nDim == 2) detR2 = (r11*r22)*(r11*r22);
		else detR2 = (r11*r22*r33)*(r11*r22*r33);
		
		/*--- Detect singular matrices ---*/
		
		if (abs(detR2) <= EPS) { detR2 = 1.0; singular = true; }
		//cout << "Determinant: " << detR2 << endl;
		/*--- S matrix := inv(R)*traspose(inv(R)) ---*/
		
		if (singular) {
		for (iDim = 0; iDim < nDim; iDim++)
			for (jDim = 0; jDim < nDim; jDim++)
			smatrix[iDim][jDim] = 0.0;
		}
		else {
		if (nDim == 2) {
			smatrix[0][0] = (r12*r12+r22*r22)/detR2;
			smatrix[0][1] = -r11*r12/detR2;
			smatrix[1][0] = smatrix[0][1];
			smatrix[1][1] = r11*r11/detR2;
		}
		else {
			z11 = r22*r33; z12 = -r12*r33; z13 = r12*r23-r13*r22;
			z22 = r11*r33; z23 = -r11*r23; z33 = r11*r22;
			smatrix[0][0] = (z11*z11+z12*z12+z13*z13)/detR2;
			smatrix[0][1] = (z12*z22+z13*z23)/detR2;
			smatrix[0][2] = (z13*z33)/detR2;
			smatrix[1][0] = smatrix[0][1];
			smatrix[1][1] = (z22*z22+z23*z23)/detR2;
			smatrix[1][2] = (z23*z33)/detR2;
			smatrix[2][0] = smatrix[0][2];
			smatrix[2][1] = smatrix[1][2];
			smatrix[2][2] = (z33*z33)/detR2;
		}
		}
		/*--- Computation of the gradient: S*c ---*/
		
		for (iDim = 0; iDim < nDim; iDim++) {
			product = 0.0;
			for (jDim = 0; jDim < nDim; jDim++) {
			product += smatrix[iDim][jDim]*cvector[jDim];
			}
			fluidSolver->node[iPoint]->SetGradient_Blockage(iDim, product);
		}
		
		/*
		AD::SetPreaccOut(node[iPoint]->GetGradient_Primitive(), nPrimVarGrad, nDim);
		AD::EndPreacc();
		*/
	}
	

}

void CBodyForceModelSolver::ComputeBFMSources(CConfig *config, CGeometry *geometry, CSolver *fluidsolver){
	ComputeBodyForce_Source(config, geometry, fluidsolver);
	SU2_MPI::Barrier(MPI_COMM_WORLD);
	if(kind_bfm==1){
		ComputeBlockage_Source(config, geometry, fluidsolver);
		SU2_MPI::Barrier(MPI_COMM_WORLD);
	}
}
void CBodyForceModelSolver::ComputeBodyForce_Source(CConfig *config, CGeometry *geometry, CSolver *fluidsolver){
	bool thollet=true;
	ComputeRelVelocity(geometry, config, fluidsolver);
	if(kind_bfm==1){
		BFM_Thollet(config, geometry, fluidsolver);
	}
	if(kind_bfm==0){
		BFM_Hall(config, geometry, fluidsolver);
	}
	
}

void CBodyForceModelSolver::ComputeBlockage_Source(CConfig *config, CGeometry *geometry, CSolver *fluidsolver){
	/*
	This function computes the local residual vector due to metal blockage for the body-force Euler solver.
	It uses the blockage factor and its derivatives in x- and r-direction provided by the interpolation function which is 
	performed during the solver initialization.
	*/
	
	// Obtaining the dimension count, flow and primitive variables and the interpolated body-force parameters.
	unsigned long iPoint;
	unsigned long nDim = geometry->GetnDim();
	su2double *U_i, *V_i, enthalpy, *Coord_i, Density, *Geometric_Parameters;
	int iDim;
	
	//  Defining the blockage residual vector.
	su2double Blockage_Vector[nDim + 2];
    //cout << "Blockage function is being called" <<endl;
	// The function loops over all points in the zone, calculating and storing the residual blockage vector for each
	// respective node.
	su2double b = 1.0;
	su2double bffac = 0.0;
	su2double BGradient[nDim], Blockage_Div = 0.0;
	su2double MomGradient[nDim];
	su2double rho;
	su2double Velocity[nDim];
	for ( iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++ ) {
		// Getting the solution and primitive variables of the current node.
		U_i = fluidsolver->node[iPoint]->GetSolution();
        V_i = fluidsolver->node[iPoint]->GetPrimitive();
		enthalpy = fluidsolver->node[iPoint]->GetEnthalpy();
		// Getting the node coordinates.
		Coord_i = geometry->node[iPoint]-> GetCoord();
		
		// The blockage factor is extracted at the node, resulting from the interpolation during initialization.
		Geometric_Parameters = fluidsolver->node[iPoint]->GetBodyForceParameters();
		bffac = Geometric_Parameters[0];
		b = Geometric_Parameters[1];
		
		
		// Obtaining blockage derivative from gradient field, calculated during initialization.
		for(iDim=0; iDim < nDim; iDim++){
			BGradient[iDim] = fluidsolver->node[iPoint]->GetGradient_Blockage(iDim);
		}
		
		// Fluid density
		rho = U_i[0];
		
		// Initializing blockage divergence with zero.
		Blockage_Div = 0.0;
		// Calculating blockage divergence term.
		for(iDim=0; iDim < nDim; iDim ++){
			Velocity[iDim] = U_i[iDim + 1]/rho;
			Blockage_Div += rho * Velocity[iDim] * BGradient[iDim];
		}
		
		// Inserting values into the blockage residual vector.
		if(nDim == 2){
			Blockage_Vector[0] = -(1 / b) * Blockage_Div;
			Blockage_Vector[1] = -(1 / b) * Velocity[0] * Blockage_Div;
			Blockage_Vector[2] = -(1 / b) * Velocity[1] * Blockage_Div;
			Blockage_Vector[3] = -(1 / b) * enthalpy * Blockage_Div;
		}else{

			Blockage_Vector[0] = -bffac * (1 / b) * Blockage_Div;
			Blockage_Vector[1] = -bffac * (1 / b) * Velocity[0] * Blockage_Div;
			Blockage_Vector[2] = -bffac * (1 / b) * Velocity[1] * Blockage_Div;
			Blockage_Vector[3] = -bffac * (1 / b) * Velocity[2] * Blockage_Div;
			Blockage_Vector[4] = -bffac * (1 / b) * enthalpy * Blockage_Div;
		}
		
		// Storing the blockage residual vector at the current node.
		fluidsolver->node[iPoint]->SetBlockage_Source(Blockage_Vector);
	};
}

bool CBodyForceModelSolver::RC_inclusion(su2double axial, su2double radial, vector<vector<su2double>>x_side){

	unsigned short n_int {0};
	bool inside = false;
	unsigned short n_nodes = x_side.size();
	unsigned short i_next;
	su2double ax_min {0};
	su2double determinant, S, R;
	for(unsigned short i=0; i<n_nodes; i++){
		if(x_side[i][0] < ax_min) ax_min = x_side[1][0];
	}
	ax_min -= 1;
	for(unsigned short i=0; i<n_nodes; i++){
		i_next = (i+1) % n_nodes;
		determinant = (axial - ax_min)*(x_side[i_next][1] - x_side[i][1]);
		if(determinant != 0){
			S = (1 / determinant) * (axial - ax_min)*(radial - x_side[i][1]);
			R = (1 / determinant) * ((x_side[i][1] - x_side[i_next][1])*(ax_min - x_side[i][0]) + 
			(x_side[i_next][0] - x_side[i][0])*(radial - x_side[i][1]));
			if((S >= 0 && S <= 1) && (R >= 0 && R <= 1)){
				n_int ++;
			}
		}
	}
	if(n_int % 2 != 0){
		inside = true;
	}
	return inside;
}
su2double CBodyForceModelSolver::Bilinear_Interp(vector<vector<su2double>> x_side, vector<su2double> y_side, su2double axial, su2double radial){
	// This function allows for bilinear interpolation over a given set of corner points stored in x_side, with values stored in y_side. 
	// This is done through transforming the corner point coordinates to a unit square.

	su2double m, l, interp_solution;

	// Axial, unit distance
	m = (axial - x_side.at(0).at(0))/(x_side.at(3).at(0) - x_side.at(0).at(0));

	// Lower and upper radii
	su2double r_low = x_side.at(0).at(1) + (x_side.at(3).at(1) - x_side.at(0).at(1))*(axial - x_side.at(0).at(0));
	su2double r_high = x_side.at(1).at(1) + (x_side.at(2).at(1) - x_side.at(1).at(1))*(axial - x_side.at(0).at(0));

	// Radial, unit distance
	l = (radial - r_low)/(r_high - r_low);

	// Interpolating solution based on unit coordinates
	interp_solution = (1 - m)*(1 - l)*y_side.at(0) +(1 - m)*l*y_side.at(1) + m*(1 - l)*y_side.at(3) + m*l*y_side.at(2);
	return interp_solution;
}
su2double CBodyForceModelSolver::DW_average(vector<vector<su2double>> x_side, vector<su2double> y_side, su2double axial, su2double radial){
	su2double eNum{0};
	su2double deNom{0};
	su2double dist{0};
	su2double y_interp;
	bool intersect = false;
	unsigned short n_nodes = y_side.size();
	for(unsigned short i=0; i<n_nodes; i++){
		dist = sqrt(pow(axial - x_side[i][0], 2) + pow(radial - x_side[i][1], 2));
		if(dist == 0){
			y_interp = y_side[i];
			i = n_nodes;
			intersect = true;
		}else{
			eNum += y_side[i]/dist;
			deNom += 1 / dist;
		}
	}
	if(!intersect){
		y_interp = eNum / deNom;
	}
	return y_interp;
}
CBodyForceModelSolver::CBodyForceModelSolver::~CBodyForceModelSolver(){
	for(unsigned long iPoint=0; iPoint<nPoint_geom; iPoint++){
		delete [] proj_vector_axial[iPoint];
		delete [] proj_vector_tangential[iPoint];
		delete [] proj_vector_radial[iPoint];
		delete [] relative_vel[iPoint];
		delete [] cyl_coordinates[iPoint];
	}
	delete [] cyl_coordinates;
	delete [] proj_vector_axial;
	delete [] proj_vector_tangential;
	delete [] proj_vector_radial;
	delete [] relative_vel;

	for(size_t i_row=0; i_row<n_rows; i_row++){
		for(size_t i_sec=0; i_sec<n_sec; i_sec++){
			for(size_t i_ax=0; i_ax<n_axial; i_ax++){
				delete [] camber_normal_data[i_row][i_sec][i_ax];
			}
			delete [] ax_coordinate[i_row][i_sec];
			delete [] rad_coordinate[i_row][i_sec];
			delete [] camber_normal_data[i_row][i_sec];
			delete [] blockage_data[i_row][i_sec];
			delete [] X_le_data[i_row][i_sec];
		}
		delete [] ax_coordinate[i_row];
		delete [] rad_coordinate[i_row];
		delete [] camber_normal_data[i_row];
		delete [] blockage_data[i_row];
		delete [] X_le_data[i_row];
	}
	delete [] ax_coordinate;
	delete [] rad_coordinate;
	delete [] camber_normal_data;
	delete [] blockage_data;
	delete [] X_le_data;
	delete [] row_rotation;
	delete [] row_blade_count;
	cout << "Deleted CBodyForceModel container." << endl;
}

su2double CBodyForceModelSolver::vector_dot_product(vector<su2double> v_1, vector<su2double> v_2){

	int n = v_1.size();
	su2double dot_product {0};
	for(int i=0; i<n; i++){
		dot_product += v_1.at(i) * v_2.at(i);
	}
	return dot_product;
}
