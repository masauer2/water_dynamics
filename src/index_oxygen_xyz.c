// I have a trajectory of H2O molecules with an H3O+ ion. I want to track the extra proton.

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>

#define NUM_OF_FRAMES 90677 // Total frames are there in our simulation
#define NUM_OF_ATOMS 385 // Total number of atoms in each of the frames
#define BOX_SIZE 15.8122 // Simulation box size


struct BookKeeping{
	int o[3];
};
/*
Stores the ids of the hydronium ion from the id file
*/

struct hydronium_idx{
	int i[NUM_OF_FRAMES]; //index [0-90676]
	int idx[NUM_OF_FRAMES]; //hydronium oxygen id [0-385]
};

/*
Stores an instance of H3O molecule
*/
struct H3O{
	double o1[3]; //Stores coordinates (x,y,z) of hydronium oxygen 1
	double h1[3]; //Stores coordinates (x,y,z) of hydronium hydrogen 1 
	double h2[3]; //Stores coordinates (x,y,z) of hydronium hydrogen 2
	double h3[3]; //Stores coordinates (x,y,z) of hydronium hydrogen 3 

	// Stores index of the atoms
	int index_o1;
	int index_h1; 
	int index_h2;
	int index_h3;

	// We need to find trajectories where the index of the hydronium ion is constant
	// Store the index according to hydronium_idx[frame_number].idx[n]
	int ion_index;

	// What are the three closest water molecule ids?
	int closest_waters[3];
};

/*
Stores an instance of H2O molecule
*/
struct H2O{
	double o1[3]; //Stores coordinates (x,y,z) of water oxygen 1
	double h1[3]; //Stores coordinates (x,y,z) of water hydrogen 1
	double h2[3]; //Stores coordinates (x,y,z) of water hydrogen 2

	// Stores index of the atoms
	int index_o1;
	int index_h1;
	int index_h2;
};

/*
Stores the coordinates of the simulation timesteps.
Here, the size of each array is the number of atoms per simulation frame.
*/
struct Frame{
	char atom_name[NUM_OF_ATOMS]; //Stores all of the atom names for a frame (either O or H)
	double x[NUM_OF_ATOMS]; //Stores x coordinates for a frame
	double y[NUM_OF_ATOMS]; //Stores y coordinates for a frame
	double z[NUM_OF_ATOMS]; //Stores z coordinates for a frame
};

struct Stack{
	struct H3O ion;
	struct H2O water[127];
};

struct oxygens{
	double coords[3];
};

/*
Function to read the id file and store the hydronium ion index at each timestep.
Note that regardless of how many frames we want to read in, we always read in all of the indexes.

NOTE: THE PROVIDED INDEX FILE HAS THE OXYGEN NUMBER, WE CONVERT FROM OXYGEN NUMBER TO ATOM NUMBER.
*/
void read_ids(struct hydronium_idx *arr){
	FILE *fp2; //Pointer to file
    int i = 0; 
	int FILELENGTH = NUM_OF_FRAMES;
	fp2 = fopen("../id_hydronium.dat", "r");

	// Read in all of the indeces and store in struct
	while(i < FILELENGTH){
		fscanf(fp2, "%d %d", &arr->i[i], &arr->idx[i]);
		if(arr->idx[i] != 0){
			arr->idx[i] = arr->idx[i]*3+1;
		}
		i++;
	}
	fclose(fp2);
}

void apply_pbc_coordinates(double *vec){

	if(vec[0] < -BOX_SIZE * 0.5 ){
		vec[0] = vec[0] + BOX_SIZE;
	}
	else if(vec[0] >= BOX_SIZE * 0.5 ){
	        vec[0] = vec[0] - BOX_SIZE;
	}

	// Check coordinate in y-direction
	if( vec[1] < -BOX_SIZE * 0.5 ){
	        vec[1] = vec[1] + BOX_SIZE;
	}
	else if(vec[1] >= BOX_SIZE * 0.5 ){
	        vec[1] = vec[1] - BOX_SIZE;
	}
	// Check coordinate in z-direction
	if( vec[2] < -BOX_SIZE * 0.5 ){
	        vec[2] = vec[2] + BOX_SIZE;
	}
	else if(vec[2] >= BOX_SIZE * 0.5 ){
	        vec[2] = vec[2] - BOX_SIZE;
	}

}

void apply_pbc_distances(double *vec1, double *vec2, double *store_here){
	double x_rsize = 1.0 / BOX_SIZE;
	double distance_vector[3];
	// Compute the distance vector in the x-direction
	distance_vector[0] = vec1[0] - vec2[0];
	distance_vector[1] = vec1[1] - vec2[1];
	distance_vector[2] = vec1[2] - vec2[2];
	distance_vector[0] -= BOX_SIZE * nearbyint(distance_vector[0] * x_rsize);
	distance_vector[1] -= BOX_SIZE * nearbyint(distance_vector[1] * x_rsize);
	distance_vector[2] -= BOX_SIZE * nearbyint(distance_vector[2] * x_rsize);
	store_here[0] = distance_vector[0];
	store_here[1] = distance_vector[1];
	store_here[2] = distance_vector[2];
}
/*
Using the minimum image convention and pbc conditions, compute the distance between two atoms.
We will use this function to compute the distance (bond length) between a given oxygen and all possible hydrogens.
*/
double compute_distance(double vec1[], double vec2[]){
	double distance_vector[3];
	double distance_magnitude;
	double coordinates[3] = {vec1[0], vec1[1], vec1[2]};
	double coordinates_2[3] = {vec2[0], vec2[1], vec2[2]};
	apply_pbc_distances(coordinates, coordinates_2, distance_vector);
	
	
	distance_magnitude = sqrt(pow(distance_vector[0],2) + pow(distance_vector[1],2) + pow(distance_vector[2],2));	
	return distance_magnitude;
}

double compute_hbond_angle(double h1[], double o1[], double o2[]){
	double distance_vector[3];
	double distance_vector_2[3];
	double distance_magnitude;
	
	apply_pbc_distances(h1, o1, distance_vector);
	apply_pbc_distances(o2, o1, distance_vector_2);
	distance_magnitude = sqrt(pow(distance_vector[0],2) + pow(distance_vector[1],2) + pow(distance_vector[2],2));	
	double angle = (180/3.14)*acos((distance_vector[0]*distance_vector_2[0]+distance_vector[1]*distance_vector_2[1]+distance_vector[2]*distance_vector_2[2])/(sqrt(distance_vector[0]*distance_vector[0] + distance_vector[1]*distance_vector[1] + distance_vector[2]*distance_vector[2])*sqrt(distance_vector_2[0]*distance_vector_2[0] + distance_vector_2[1]*distance_vector_2[1] + distance_vector_2[2]*distance_vector_2[2])));
	return angle;
}

/*
Assign the water molecule coordinates to the respective struct.

water_id is the unique water id number as stored in the struct water.
However, the water id as stored in the original coordinate file is actually oxygen_id,
since oxygen_id stores the raw index of the water oxygen from the coordinate file.
We match the indexing of the water molecules to the indexing found in id_hydronium.dat.
*/
void assign_waters(int water_id, int oxygen_id, int* hydrogen_index, struct Frame *analysis_frame, struct H2O *waters){

				//Assign the oxygen coordinates
				waters[water_id].o1[0] = analysis_frame->x[oxygen_id];
				waters[water_id].o1[1] = analysis_frame->y[oxygen_id];
				waters[water_id].o1[2] = analysis_frame->z[oxygen_id];

				//Assign the hydrogen coordinates
				waters[water_id].h1[0] = analysis_frame->x[hydrogen_index[0]];
				waters[water_id].h1[1] = analysis_frame->y[hydrogen_index[0]];
				waters[water_id].h1[2] = analysis_frame->z[hydrogen_index[0]];

				waters[water_id].h2[0] = analysis_frame->x[hydrogen_index[1]];
				waters[water_id].h2[1] = analysis_frame->y[hydrogen_index[1]];
				waters[water_id].h2[2] = analysis_frame->z[hydrogen_index[1]];
				//Assign the oxygen index
				waters[water_id].index_o1 = oxygen_id;

				//Assign the hydrogen indeces
				waters[water_id].index_h1 = hydrogen_index[0];
				waters[water_id].index_h2 = hydrogen_index[1];

}


/*
Assign the hydronium ion coordinates to the respective struct.
Print out some info about the hydronium ion.
*/
void assign_ion(int oxygen_id, int* hydrogen_index, struct Frame *analysis_frame, struct H3O *hydronium){

	//Assign the oxygen coordinates
	hydronium->o1[0] = analysis_frame->x[oxygen_id];
    hydronium->o1[1] = analysis_frame->y[oxygen_id];
    hydronium->o1[2] = analysis_frame->z[oxygen_id];

	//Assign the hydrogen coordinates
	hydronium->h1[0] = analysis_frame->x[hydrogen_index[0]];
	hydronium->h1[1] = analysis_frame->y[hydrogen_index[0]];
	hydronium->h1[2] = analysis_frame->z[hydrogen_index[0]];

	hydronium->h2[0] = analysis_frame->x[hydrogen_index[1]];
	hydronium->h2[1] = analysis_frame->y[hydrogen_index[1]];
	hydronium->h2[2] = analysis_frame->z[hydrogen_index[1]];

	hydronium->h3[0] = analysis_frame->x[hydrogen_index[2]];
	hydronium->h3[1] = analysis_frame->y[hydrogen_index[2]];
	hydronium->h3[2] = analysis_frame->z[hydrogen_index[2]];

	//Assign the oxygen index
	hydronium->index_o1 = oxygen_id;

	//Assign the hydrogen indeces
    hydronium->index_h1 = hydrogen_index[0];
	hydronium->index_h2 = hydrogen_index[1];
	hydronium->index_h3 = hydrogen_index[2];

	hydronium->ion_index= oxygen_id;
	//Print out some basic info about the hydronium ion

}

void assign_closest_waters(int water_idx[], struct BookKeeping *closestOs){
	closestOs -> o[0] = water_idx[0];
	closestOs -> o[1] = water_idx[1];
	closestOs -> o[2] = water_idx[2];

}
/*
Use the distance formula in order to find the three closest hydrogens to the oxygen molecule based off of
the provided id.

int oxygen_id = the hydronium oxygen
struct Frame *analysis_frame = the frame we are analyzing
struct H3O *hydronium = the hydronium coordinates

Start with a list of distances that are larger than the box length. If we compute a distance that is smaller than
one of these elements, replace that list element. We will repeat this process until we are left with the three
shortest bond lengths. These three hydrogens will define the hydronium ion hydrogens.
*/
void identify_ion(int oxygen_id, struct Frame *analysis_frame, struct H3O *hydronium){

	double distance; // Current distance that is being tested
	double smallest_distances[] = {16,17,18}; // List containing shortest hydrogens to oxygen
    int hydrogen_idx[] = {-1,-1,-1}; // Index list corresponding to smallest_distances[]
	double ion_vector[3];
	double hydrogen_vector[3];
	/*
	Loop over all of the atoms in a given frame and look for hydrogen atoms. Once we find a hydrogen atom,
	we will test to see if the hydrogen belongs to the hydronium oxygen.
	*/
	for(int i = 1; i < NUM_OF_ATOMS+1; i++){
		if(analysis_frame->atom_name[i] == 'H'){
			ion_vector[0] = analysis_frame->x[oxygen_id];
			ion_vector[1] = analysis_frame->y[oxygen_id];
			ion_vector[2] = analysis_frame->z[oxygen_id];
			hydrogen_vector[0] = analysis_frame->x[i];
			hydrogen_vector[1] = analysis_frame->y[i];
			hydrogen_vector[2] = analysis_frame->z[i];

			distance = compute_distance(ion_vector, hydrogen_vector); //Compute the distance

			if(distance < smallest_distances[2] && distance > smallest_distances[1]){
				smallest_distances[2] = distance; 
				hydrogen_idx[2] = i;
			}

			else if(distance < smallest_distances[1] && distance > smallest_distances[0]){
				smallest_distances[2] = smallest_distances[1]; 
				smallest_distances[1] = distance; 
				hydrogen_idx[2] = hydrogen_idx[1];
				hydrogen_idx[1] = i;
			}

			else if(distance < smallest_distances[0]){
				smallest_distances[2] = smallest_distances[1]; 
				smallest_distances[1] = smallest_distances[0]; 
				smallest_distances[0] = distance;
			    hydrogen_idx[2] = hydrogen_idx[1];
				hydrogen_idx[1] = hydrogen_idx[0];	
				hydrogen_idx[0] = i;
			}
		}
	}
	
	//Assign the coordinates of the hydrogen to the hydronium oxygen and print out information about the ion
	assign_ion(oxygen_id, hydrogen_idx, analysis_frame, hydronium);
}

/*
Use the distance formula in order to find the two closest hydrogens to each of the oxygen molecules that
are not the hydronium ion oxygen.

The usage of this function is similat to identify_ion(). However, we are not given the indexes of all oxygens.
Therefore, we loop to find all the oxygens, and then re-loop to find its corresponding hydrogens. 
Then, we build our water molecule.
*/
void identify_waters(int oxygen_id, struct Frame *analysis_frame, struct H2O *waters){

	//Test distances
	double distance;
	double distances[] = {16,17};
	int indexe[] = {-1,-1};
	double ion_vector[3];
	double hydrogen_vector[3];

	int k = -1;
	for(int i = 0; i < NUM_OF_ATOMS; i++){
		if(analysis_frame->atom_name[i] == 'O' && i != oxygen_id){
			k++;
			double distances[] = {16,17}; //dup?
			int indexe[] = {-1,-1};
			for(int j = 0; j < NUM_OF_ATOMS; j++){
				if(analysis_frame->atom_name[j] == 'H'){
					ion_vector[0] = analysis_frame->x[i];
					ion_vector[1] = analysis_frame->y[i];
					ion_vector[2] = analysis_frame->z[i];
					hydrogen_vector[0] = analysis_frame->x[j];
					hydrogen_vector[1] = analysis_frame->y[j];
					hydrogen_vector[2] = analysis_frame->z[j];
					distance = compute_distance(ion_vector, hydrogen_vector);
					if(distance < distances[1] && distance > distances[0]){ 
						distances[1] = distance; 
						indexe[1] = j;
					}
					else if(distance < distances[0]){
						distances[1] = distances[0]; 
						distances[0] = distance;
						indexe[1] = indexe[0];	
						indexe[0] = j;
					}
				}
			}

			// Build water molecules
			assign_waters(k, i, indexe, analysis_frame, waters);
			}
		}
}


/*
USAGE: ./index_waters begin_frame end_frame
*/

void identify_closest_oxygens(struct Stack *stack, struct BookKeeping *closestOs){
	// Test distances
	// We seek the three closest oxygens to our hydronium ion
	// Is this in eigen/zundel form? long lived pair = 1-2ps
	double distance = 18.0;
	double distanceh1;
	int index_close;
	double vec1[3];
	double vec2[3];
	double angle;
	double smallest_distances[] = {16,17,18};
	int water_idx[] = {-1,-1,-1};
	double angles[] = {-1,-1,-1};
	int k = -1;


	//Iterate through all the possible water molecules. Remember that we have 128 oxygens in our system, but one is the hydronium ion
		for(int i = 0; i < 127; i++){
			distanceh1 = compute_distance(stack->ion.h1,stack->water[i].o1);
			angle = compute_hbond_angle(stack->ion.h1, stack->ion.o1, stack->water[i].o1);
			if(distanceh1 < distance && angle < 45){
				distance = distanceh1;
				index_close = i;
				angles[0] = angle;
			}
		}

		closestOs -> o[0] = index_close;
		distance = 16.0;

		for(int i = 0; i < 127; i++){
			distanceh1 = compute_distance(stack->ion.h2,stack->water[i].o1);
			angle = compute_hbond_angle(stack->ion.h2, stack->ion.o1, stack->water[i].o1);
			if(distanceh1 < distance && angle < 45){
				distance = distanceh1;
				index_close = i;
				angles[1] = angle;
			}
		}
		closestOs -> o[1] = index_close;
		distance = 16.0;

		for(int i = 0; i < 127; i++){
			distanceh1 = compute_distance(stack->ion.h3,stack->water[i].o1);
			angle = compute_hbond_angle(stack->ion.h3, stack->ion.o1, stack->water[i].o1);
			if(distanceh1 < distance && angle < 45){
				distance = distanceh1;
				index_close = i;
				angles[2] = angle;
			}
		}
		closestOs -> o[2] = index_close;
	}
	



void frame_to_stack(struct Stack *stack, struct H3O *hydronium, struct H2O *waters){
	stack->ion = *hydronium;
	for(int i = 0; i < 127; i++){
		stack->water[i] = waters[i];
	}
}

int test_lists(struct BookKeeping *keep, struct BookKeeping *old_keep){
	int flag = 0;
	for(int i = 0; i< 3; i++){
			if(old_keep->o[i] == keep->o[0] || old_keep->o[i] == keep->o[1] || old_keep->o[i] == keep->o[2]){
				flag = 0;
			}
			else{
				flag = 1;
				return flag;
			}
		}
		return flag;
	}

void identify_index(struct Stack *last_stack, struct Stack *current_stack, struct H3O *hydronium, struct H2O *waters){
	struct Stack *copy_stack = current_stack;
	current_stack->ion = *hydronium;
	double checker;
	double checker2;
	double checker3;
	double smallest_distance = BOX_SIZE;
	double smallest_distance_2 = BOX_SIZE;
	double smallest_distance_3 = BOX_SIZE;
	int smallest_distance_idx;
	for(int i = 0; i < 127; i++){
		smallest_distance = BOX_SIZE;
		smallest_distance_2 = BOX_SIZE;
		smallest_distance_3 = BOX_SIZE;
		smallest_distance_idx = 0;
		for(int j = 0; j < 127; j++){
			checker = compute_distance(current_stack->water[i].o1, last_stack->water[j].o1);
			checker2 = compute_distance(current_stack->water[i].h1, last_stack->water[j].h1);
			checker3 = compute_distance(current_stack->water[i].h2, last_stack->water[j].h2);
			if (checker < smallest_distance && checker2 < smallest_distance_2 && checker3 < smallest_distance_3){
				smallest_distance = checker;
				smallest_distance_idx = j;
			}
		}
	//current_stack->water[i] = copy_stack->water[smallest_distance_idx];
	}
}

int identify_traj(struct Stack *current_stack, double os[][3][4], double traj[3][4]){
	double checker;
	double smallest_distance = BOX_SIZE;
	double o1[] = {os[0][0][0], os[0][1][0],os[0][2][0]};
	double o2[] = {os[0][0][1],os[0][1][1],os[0][2][1]};
	double o3[] = {os[0][0][2],os[0][1][2],os[0][2][2]};
	double o4[] = {os[0][0][3],os[0][1][3],os[0][2][3]};

	int o1_index;
	int o2_index;
	int o3_index;
	int o4_index;
	int id = 0;

	//O1
	for(int i = 0; i < 128; i++){
			if(i == 0){
				checker = compute_distance(current_stack->ion.o1, o1);
				if (checker < smallest_distance){
					smallest_distance = checker;
					o1_index = i;
					traj[0][0] = current_stack->ion.o1[0];
					traj[1][0] = current_stack->ion.o1[1];
					traj[2][0] = current_stack->ion.o1[2];
				}	
			}

			else if(i != 0){
				checker = compute_distance(current_stack->water[i-1].o1, o1);
				if (checker < smallest_distance){
					smallest_distance = checker;
					o1_index = i;
					traj[0][0] = current_stack->water[i-1].o1[0];
					traj[1][0] = current_stack->water[i-1].o1[1];
					traj[2][0] = current_stack->water[i-1].o1[2];
				}
			}
		
	}

	if(o1_index == 0){
		id = 1;
	}

	smallest_distance = BOX_SIZE;
	//O2
	for(int i = 0; i < 128; i++){
		if(i == 0){
			checker = compute_distance(current_stack->ion.o1, o2);
			if (checker < smallest_distance){
				smallest_distance = checker;
				o2_index = i;
				traj[0][1] = current_stack->ion.o1[0];
				traj[1][1] = current_stack->ion.o1[1];
				traj[2][1] = current_stack->ion.o1[2];
			}
		}

		else if(i != 0){
			checker = compute_distance(current_stack->water[i-1].o1, o2);
			if (checker < smallest_distance){
				smallest_distance = checker;
				o2_index = i;
				traj[0][1] = current_stack->water[i-1].o1[0];
				traj[1][1] = current_stack->water[i-1].o1[1];
				traj[2][1] = current_stack->water[i-1].o1[2];
			}
		}
		
	}
	if(o2_index == 0){
			id = 2;
		}
	

	smallest_distance = BOX_SIZE;
	for(int i = 0; i < 128; i++){
		if(i == 0 ){
			checker = compute_distance(current_stack->ion.o1, o3);
			if (checker < smallest_distance){
				smallest_distance = checker;
				o3_index = i;
				traj[0][2] = current_stack->ion.o1[0];
				traj[1][2] = current_stack->ion.o1[1];
				traj[2][2] = current_stack->ion.o1[2];
			}
		}

		else if(i != 0 ){
			checker = compute_distance(current_stack->water[i-1].o1, o3);
			if (checker < smallest_distance){
				smallest_distance = checker;
				o3_index = i;
				traj[0][2] = current_stack->water[i-1].o1[0];
				traj[1][2] = current_stack->water[i-1].o1[1];
				traj[2][2] = current_stack->water[i-1].o1[2];
			}
		}
	}
	if(o3_index == 0){
		id = 3;
	}
	

	smallest_distance = BOX_SIZE;
	for(int i = 0; i < 128; i++){
		if(i == 0  ){
				checker = compute_distance(current_stack->ion.o1, o4);
				if (checker < smallest_distance){
					smallest_distance = checker;
					o4_index = i;
					traj[0][3] = current_stack->ion.o1[0];
					traj[1][3] = current_stack->ion.o1[1];
					traj[2][3] = current_stack->ion.o1[2];
				}
		}

		else if(i != 0){
			checker = compute_distance(current_stack->water[i-1].o1, o4);
			if (checker < smallest_distance){
				smallest_distance = checker;
				o4_index = i;
				traj[0][3] = current_stack->water[i-1].o1[0];
				traj[1][3] = current_stack->water[i-1].o1[1];
				traj[2][3] = current_stack->water[i-1].o1[2];
			}
		}
	}
	if(o4_index == 0){
		id = 4;
	}
	
	printf("ID: %d\n", id);
	return id - 1;
}

void check_errors(int begin_frame, int end_frame){
	int error_state = 1; // Error state
	if(end_frame < begin_frame){
		error_state = error_state * 3;
	}

	if(begin_frame < 0 || begin_frame > 90676){
		error_state = error_state * 5;
	}
	
	if(end_frame < 0 || end_frame > 90676 || end_frame < begin_frame){
		error_state = error_state * 7;
	}

	if(error_state % 3 == 0){
		printf("\n-----------------------------------------\n\nEnd frame is before than starting frame.\n\n-----------------------------------------\n");
	}
	if(error_state % 5 == 0){
		printf("\n------------------------------------------------------------------\n\nStarting frame is out of range. \nThe beginning frame needs to take on a value between 0 and 90676.\n\n------------------------------------------------------------------\n");
	}
	if(error_state % 7 == 0){
		printf("\n------------------------------------------------------------------\n\nEnd frame is out of range. \nThe end frame needs to take on a value between 0 and 90676.\n\n------------------------------------------------------------------\n");
	}
	if(error_state != 1){
		printf("\n------------------------------------------------\n\nThe correct usage of this program is as follows.\n\n ./index_waters begin_frame end_frame.\n\n------------------------------------------------\n");
		exit(0);
	}		
}

void calculate_OH3_distance(struct Frame *analysis_frame, double oxyCoords[], double hCoords[][3][3]){
	
	double distance; // Current distance that is being tested
	double smallest_distances[] = {16,17,18}; // List containing shortest hydrogens to oxygen
	double ion_vector[3];
	double hydrogen_vector[3];
	int hydrogen_idx[3];
	for(int i = 1; i < NUM_OF_ATOMS+1; i++){
		if(analysis_frame->atom_name[i] == 'H'){
			ion_vector[0] = oxyCoords[0];
			ion_vector[1] = oxyCoords[1];
			ion_vector[2] = oxyCoords[2];
			hydrogen_vector[0] = analysis_frame->x[i];
			hydrogen_vector[1] = analysis_frame->y[i];
			hydrogen_vector[2] = analysis_frame->z[i];

			distance = compute_distance(ion_vector, hydrogen_vector); //Compute the distance

			if(distance < smallest_distances[2] && distance > smallest_distances[1]){
				smallest_distances[2] = distance; 
				hydrogen_idx[2] = i;
			}

			else if(distance < smallest_distances[1] && distance > smallest_distances[0]){
				smallest_distances[2] = smallest_distances[1]; 
				smallest_distances[1] = distance; 
				hydrogen_idx[2] = hydrogen_idx[1];
				hydrogen_idx[1] = i;
			}

			else if(distance < smallest_distances[0]){
				smallest_distances[2] = smallest_distances[1]; 
				smallest_distances[1] = smallest_distances[0]; 
				smallest_distances[0] = distance;
			    hydrogen_idx[2] = hydrogen_idx[1];
				hydrogen_idx[1] = hydrogen_idx[0];	
				hydrogen_idx[0] = i;
			}
		}
	}
	hCoords[1][0][0] = analysis_frame->x[hydrogen_idx[0]];
	hCoords[1][1][0] =analysis_frame->y[hydrogen_idx[0]];
	hCoords[1][2][0] = analysis_frame->z[hydrogen_idx[0]];
	hCoords[1][0][1] = analysis_frame->x[hydrogen_idx[1]];
	hCoords[1][1][1] =analysis_frame->y[hydrogen_idx[1]];
	hCoords[1][2][1] = analysis_frame->z[hydrogen_idx[1]];
	hCoords[1][0][2] = analysis_frame->x[hydrogen_idx[2]];
	hCoords[1][1][2] =analysis_frame->y[hydrogen_idx[2]];
	hCoords[1][2][2] = analysis_frame->z[hydrogen_idx[2]];
}

void calculate_OH2_distance(struct Frame *analysis_frame, double oxyCoords[], double hCoords[][3][2]){
	
	double distance; // Current distance that is being tested
	double smallest_distances[] = {16,17}; // List containing shortest hydrogens to oxygen
	double ion_vector[3];
	double hydrogen_vector[3];
	int hydrogen_idx[2];
	for(int i = 1; i < NUM_OF_ATOMS+1; i++){
		if(analysis_frame->atom_name[i] == 'H'){

			ion_vector[0] = oxyCoords[0];
			ion_vector[1] = oxyCoords[1];
			ion_vector[2] = oxyCoords[2];
			
			hydrogen_vector[0] = analysis_frame->x[i];
			hydrogen_vector[1] = analysis_frame->y[i];
			hydrogen_vector[2] = analysis_frame->z[i];

			distance = compute_distance(ion_vector, hydrogen_vector); //Compute the distance

			if(distance < smallest_distances[1] && distance > smallest_distances[0]){
				smallest_distances[1] = distance; 
				hydrogen_idx[1] = i;
			}
			else if(distance < smallest_distances[0]){
				smallest_distances[1] = smallest_distances[0]; 
				smallest_distances[0] = distance;
				hydrogen_idx[1] = hydrogen_idx[0];	
				hydrogen_idx[0] = i;
			}
		}
	}
	hCoords[1][0][0] = analysis_frame->x[hydrogen_idx[0]];
	hCoords[1][1][0] =analysis_frame->y[hydrogen_idx[0]];
	hCoords[1][2][0] = analysis_frame->z[hydrogen_idx[0]];
	hCoords[1][0][1] = analysis_frame->x[hydrogen_idx[1]];
	hCoords[1][1][1] =analysis_frame->y[hydrogen_idx[1]];
	hCoords[1][2][1] = analysis_frame->z[hydrogen_idx[1]];
}

int main(int argc, char* argv[]){

	/* Checking for errors */
	int reference_stack[3];
	if(argc != 3){
		printf("\n------------------------------------------------\n\nWrong number of inputs. \nThe correct usage of this program is as follows.\n\n ./index_waters begin_frame end_frame.\n\n------------------------------------------------\n");
		exit(0);
	}

	int begin_frame;
	int end_frame;
	int reference_index;
	
	/* Setting the frames to be reading */

	sscanf (argv[1], "%d",&begin_frame);
	sscanf (argv[2], "%d",&end_frame);

	double oxy[2][3][4];
	

	check_errors(begin_frame, end_frame);
	// Frames to read in with this code
	int frame_number = end_frame-begin_frame;

	clock_t begin; // Variable to time the entire code segment
	clock_t one_iteration; // Variable to see how long it takes to run one frame
	struct hydronium_idx ids; // Struct that will store the hydronium oxygen ids in an array

	read_ids(&ids); // Read in the hydronium oxygen ids
	
	reference_index = ids.idx[begin_frame];

	FILE *fp;
	FILE *filewriter;
	FILE *framewriter;
	FILE *traj;
	FILE *ratiofile;

	char *token;
	char buffer[200];
	
	fp = fopen("../water_scan-pos-1.xyz", "r");
	filewriter = fopen("../output.txt", "w");
	framewriter = fopen("../frame_change.txt", "w");
	traj = fopen("../traj_seg_fixed.xyz","w");
	ratiofile = fopen("../ratio.txt","w");
	begin= clock(); // Start timing the code

	/*
	When we are reading frames not starting from the beginning, fgets() the lines that
	we do not need
	*/
	
	if(begin_frame != 0){
		for(int garb = 0; garb < 387*begin_frame; garb++){
			fgets(buffer, 200, fp);
		}
	}
	
	/* Iterate over all of the frames that we want to read and store all of the
	atom coordinates in a Frame struct 
	
	To prevent stack overflow, we read in one frame at a time. This program will need
	to be implemented with heap storage to read in large amounts of data at once.
	*/

	int new = 1;
	int old = 0;

	int i; //Counter for current atom number in the frame

	struct Frame *current_frame =  malloc(sizeof(struct Frame));  // Frame to store coordinates	
	struct Stack *current_stack = malloc(sizeof(struct Stack)); 
	struct Stack *last_stack = malloc(sizeof(struct Stack));

	int j = 0;

	for(int k = begin_frame; k < end_frame+1; k++){
		if(k != begin_frame){
			last_stack = current_stack;
			current_stack = malloc(sizeof(struct Stack));
		}

		i = 0; //Start at the first atom
		fgets(buffer, 100, fp); // First line -> Number of atoms
		fgets(buffer, 100, fp); // Second line -> info about timestep

		while((i < NUM_OF_ATOMS)){

			// Read in one line from the file and tokenize
			fgets(buffer, 100, fp);
            token = strtok(buffer, " \t");    
			
			// First column -> atom name 
			current_frame->atom_name[i] = token[0];
			token = strtok(NULL, " \t");
			
			// Second column -> x coordinate
			current_frame->x[i] = atof(token);
			token = strtok(NULL, " \t");
			
			// Third column -> y coordinate
			current_frame->y[i] = atof(token);
			token = strtok(NULL, " \t");
			
			// Fourth column -> z coordinate
			current_frame->z[i] = atof(token);

			i++; // Read in the next atom
		}

		printf("======================================================\n");
		printf("Finished reading frame %d\n", k);

		one_iteration = clock(); // Begin timing one iteration of calculations


		int current_H3O = ids.idx[k]; //The current H3O oxygen index from the index array
 		struct H3O *hydronium = malloc (sizeof(struct H3O)); // Struct to store the atom coordinates of the hydronium ion
		struct H2O waters[127]; // Struct to store the 127 water molecule atom coordinates
		
		
		identify_ion(current_H3O, current_frame, hydronium); // Identify the hydronium ion
		identify_waters(current_H3O, current_frame, waters); // Identify all of the water molecules
		frame_to_stack(current_stack, hydronium, waters);

		
		if(k != begin_frame){
			identify_index(last_stack, current_stack, hydronium, waters);
		}


		struct BookKeeping *keep = malloc(sizeof(struct BookKeeping));
		struct BookKeeping *old_keep = malloc(sizeof(struct BookKeeping));

		identify_closest_oxygens(current_stack, keep); // Identify all of the water molecules
		identify_closest_oxygens(last_stack, old_keep); // Identify all of the water molecules

		one_iteration = clock() - one_iteration; // Stop timing the iteration

		int tracker;
		current_stack->ion.ion_index = (current_stack->ion.ion_index-1)/3;

		printf("j = %d\n", j);
		double h2c[2][3][2];
		double h3c[2][3][3];

		if(j == 0){
			oxy[old][0][0] = current_stack->ion.o1[0];
			oxy[old][1][0] = current_stack->ion.o1[1];
			oxy[old][2][0] = current_stack->ion.o1[2];
			oxy[old][0][1] = current_stack->water[keep->o[0]].o1[0];
			oxy[old][1][1] = current_stack->water[keep->o[0]].o1[1];
			oxy[old][2][1] = current_stack->water[keep->o[0]].o1[2];
			oxy[old][0][2] = current_stack->water[keep->o[1]].o1[0];
			oxy[old][1][2] = current_stack->water[keep->o[1]].o1[1];
			oxy[old][2][2] = current_stack->water[keep->o[1]].o1[2];
			oxy[old][0][3] = current_stack->water[keep->o[2]].o1[0];
			oxy[old][1][3] = current_stack->water[keep->o[2]].o1[1];
			oxy[old][2][3] = current_stack->water[keep->o[2]].o1[2];

			fprintf(traj, "%d\nt = %d\n", 13, j);

			for (int index = 0; index < 4; index++){
				if(index == 0){
					fprintf(traj, "O\t%f\t%f\t%f\n", oxy[old][0][index], oxy[old][1][index], oxy[old][2][index]);
					
					double oxygenCoordinates[] = {oxy[old][0][index], oxy[old][1][index], oxy[old][2][index]};
					calculate_OH3_distance(current_frame, oxygenCoordinates, h3c);

					for(int i = 0; i <3;i++){
						for(int j = 0; j < 3;j++){
							h3c[old][i][j] = h3c[new][i][j];
						}
					}
						
					fprintf(traj, "H\t%f\t%f\t%f\n", h3c[0][0][0], h3c[0][1][0], h3c[0][2][0]);
					fprintf(traj, "H\t%f\t%f\t%f\n", h3c[0][0][1], h3c[0][1][1], h3c[0][2][1]);
					fprintf(traj, "H\t%f\t%f\t%f\n", h3c[0][0][2], h3c[0][1][2], h3c[0][2][2]);
					
				}
				
				if(index != 0){
					fprintf(traj, "O\t%f\t%f\t%f\n", oxy[old][0][index], oxy[old][1][index], oxy[old][2][index]);
					double oxygenCoordinates[] = {oxy[old][0][index], oxy[old][1][index], oxy[old][2][index]};
					calculate_OH2_distance(current_frame, oxygenCoordinates, h2c);
					for(int i = 0; i <3;i++){
						for(int j = 0; j < 2;j++){
							h2c[old][i][j] = h2c[new][i][j];
						}
					}
					fprintf(traj, "H\t%f\t%f\t%f\n", h2c[0][0][0], h2c[0][1][0], h2c[0][2][0]);
					fprintf(traj, "H\t%f\t%f\t%f\n", h2c[0][0][1], h2c[0][1][1], h2c[0][2][1]);
				}
				
			}
			printf("O1: [%f %f %f]\n", oxy[old][0][0], oxy[old][1][0], oxy[old][2][0]);
			printf("O2: [%f %f %f]\n", oxy[old][0][1], oxy[old][1][1], oxy[old][2][1]);
			printf("O3: [%f %f %f]\n", oxy[old][0][2], oxy[old][1][2], oxy[old][2][2]);
			printf("O4: [%f %f %f]\n", oxy[old][0][3], oxy[old][1][3], oxy[old][2][3]);
			}
			



		if(j != 0){
			double trajs[3][4];
			tracker = identify_traj(current_stack, oxy, trajs);
			printf("%d\n", tracker);
			fprintf(traj, "%d\nt = %d\n", 13, j);
			for (int index = 0; index < 4; index++){
				if(index == tracker){

					oxy[new][0][index] = trajs[0][index];
					oxy[new][1][index] = trajs[1][index];
					oxy[new][2][index] = trajs[2][index];

					fprintf(traj, "O\t%f\t%f\t%f\n", oxy[new][0][index], oxy[new][1][index], oxy[new][2][index]);
			
					double oxygenCoordinates[] = {oxy[new][0][index], oxy[new][1][index], oxy[new][2][index]};
					if(index == tracker){
						calculate_OH3_distance(current_frame, oxygenCoordinates, h3c);
						fprintf(traj, "H\t%f\t%f\t%f\n", h3c[new][0][0], h3c[new][1][0], h3c[new][2][0]);
						fprintf(traj, "H\t%f\t%f\t%f\n", h3c[new][0][1], h3c[new][1][1], h3c[new][2][1]);
						fprintf(traj, "H\t%f\t%f\t%f\n", h3c[new][0][2], h3c[new][1][2], h3c[new][2][2]);
					}

					oxy[old][0][index] = oxy[new][0][index];
					oxy[old][1][index] = oxy[new][1][index];
					oxy[old][2][index] = oxy[new][2][index];
					h3c[old][0][index] = h3c[new][0][index];
					h3c[old][1][index] = h3c[new][1][index];
					h3c[old][2][index] = h3c[new][2][index];
				}
			}

			for (int index = 0; index < 4; index++){
				if(index != tracker){

					oxy[new][0][index] = trajs[0][index];
					oxy[new][1][index] = trajs[1][index];
					oxy[new][2][index] = trajs[2][index];

					fprintf(traj, "O\t%f\t%f\t%f\n", oxy[new][0][index], oxy[new][1][index], oxy[new][2][index], oxy[new][0][index]-oxy[old][0][index], oxy[new][1][index]-oxy[old][1][index], oxy[new][2][index]-oxy[old][2][index]);
				
					double oxygenCoordinates[] = {oxy[new][0][index], oxy[new][1][index], oxy[new][2][index]};
					calculate_OH2_distance(current_frame, oxygenCoordinates, h2c);
					fprintf(traj, "H\t%f\t%f\t%f\n", h2c[new][0][0], h2c[new][1][0], h2c[new][2][0]);
					fprintf(traj, "H\t%f\t%f\t%f\n", h2c[new][0][1], h2c[new][1][1], h2c[new][2][1]);

					oxy[old][0][index] = oxy[new][0][index];
					oxy[old][1][index] = oxy[new][1][index];
					oxy[old][2][index] = oxy[new][2][index];

					h2c[old][0][index] = h2c[new][0][index];
					h2c[old][1][index] = h2c[new][1][index];
					h2c[old][2][index] = h2c[new][2][index];
					}
				}
			
			printf("O1: [%f %f %f]\n", oxy[old][0][0], oxy[old][1][0], oxy[old][2][0]);
			printf("O2: [%f %f %f]\n", oxy[old][0][1], oxy[old][1][1], oxy[old][2][1]);
			printf("O3: [%f %f %f]\n", oxy[old][0][2], oxy[old][1][2], oxy[old][2][2]);
			printf("O4: [%f %f %f]\n", oxy[old][0][3], oxy[old][1][3], oxy[old][2][3]);
		}
		j++;
	}
	begin = clock() - begin; //Stop timing the code
	
	printf("Processing frames %d-%d took %f seconds.\n", begin_frame, end_frame, ((double)begin)/CLOCKS_PER_SEC);
	fclose(fp); //Close file
	fclose(filewriter);
	fclose(framewriter);
	fclose(traj);
	fclose(ratiofile);
}
