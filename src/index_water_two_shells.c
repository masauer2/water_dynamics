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
	float o1[3]; //Stores coordinates (x,y,z) of hydronium oxygen 1
	float h1[3]; //Stores coordinates (x,y,z) of hydronium hydrogen 1 
	float h2[3]; //Stores coordinates (x,y,z) of hydronium hydrogen 2
	float h3[3]; //Stores coordinates (x,y,z) of hydronium hydrogen 3 

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
	float o1[3]; //Stores coordinates (x,y,z) of water oxygen 1
	float h1[3]; //Stores coordinates (x,y,z) of water hydrogen 1
	float h2[3]; //Stores coordinates (x,y,z) of water hydrogen 2

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
	float x[NUM_OF_ATOMS]; //Stores x coordinates for a frame
	float y[NUM_OF_ATOMS]; //Stores y coordinates for a frame
	float z[NUM_OF_ATOMS]; //Stores z coordinates for a frame
};

struct Stack{
	struct H3O ion;
	struct H2O water[127];
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

/*
Using the minimum image convention and pbc conditions, compute the distance between two atoms.
We will use this function to compute the distance (bond length) between a given oxygen and all possible hydrogens.
*/
float compute_distance(int oxygen_id, int i, struct Frame *analysis_frame){
	float distance_vector[3];
	float distance_magnitude;

	// Check coordinate in x-direction
	if( analysis_frame->x[i] < -BOX_SIZE * 0.5 ){
		analysis_frame -> x[i] = analysis_frame->x[i] + BOX_SIZE;
	}
	else if(analysis_frame->x[i] >= BOX_SIZE * 0.5 ){
	        analysis_frame->x[i] = analysis_frame->x[i] - BOX_SIZE;
	}

	// Check coordinate in y-direction
	if( analysis_frame->y[i] < -BOX_SIZE * 0.5 ){
	        analysis_frame->y[i] = analysis_frame->y[i] + BOX_SIZE;
	}
	else if(analysis_frame->y[i] >= BOX_SIZE * 0.5 ){
	        analysis_frame->y[i] = analysis_frame->y[i] - BOX_SIZE;
	}
	// Check coordinate in z-direction
	if( analysis_frame->z[i] < -BOX_SIZE * 0.5 ){
	        analysis_frame->z[i] = analysis_frame->z[i] + BOX_SIZE;
	}
	else if(analysis_frame->z[i] >= BOX_SIZE * 0.5 ){
	        analysis_frame->z[i] = analysis_frame->z[i] - BOX_SIZE;
	}

	// Compute the distance vector in the x-direction
	distance_vector[0] = analysis_frame->x[i] - analysis_frame->x[oxygen_id];

	// Distance vector must also obey this convention
	// Check distance vector in x-direction
	if(distance_vector[0] > BOX_SIZE * 0.5){
		distance_vector[0] = distance_vector[0] - BOX_SIZE;
	}
	else if(distance_vector[0] <= -BOX_SIZE * 0.5){
		distance_vector[0] = distance_vector[0] + BOX_SIZE;
	}

	// Compute the distance vector in the y-direction
	distance_vector[1] = analysis_frame->y[i] - analysis_frame->y[oxygen_id];

	// Check distance vector in y-direction
	if(distance_vector[1] > BOX_SIZE * 0.5){
	        distance_vector[1] = distance_vector[1] - BOX_SIZE;
	}
	else if(distance_vector[1] <= -BOX_SIZE * 0.5){
	        distance_vector[1] = distance_vector[1] + BOX_SIZE;
	}
	
	// Compute the distance vector in the z-direction
	distance_vector[2] = analysis_frame->z[i] - analysis_frame->z[oxygen_id];	

	// Check distance vector in z-direction
	if(distance_vector[2] > BOX_SIZE * 0.5){
	        distance_vector[2] = distance_vector[2] - BOX_SIZE;
	}
	else if(distance_vector[2] <= -BOX_SIZE * 0.5){
	        distance_vector[2] = distance_vector[2] + BOX_SIZE;
	}

	// Compute the distance between the two atoms and return this value
	distance_magnitude = sqrt(pow(distance_vector[0],2) + pow(distance_vector[1],2) + pow(distance_vector[2],2));	
	return distance_magnitude;
}

float compute_oxygen_distances(float o1[], float o2[]){
	float distance_vector[3];
	float distance_magnitude;

	// Check coordinate in x-direction
	/*
	if( o1[0] < -BOX_SIZE * 0.5 ){
		 o1[0] =  o1[0] + BOX_SIZE;
	}
	else if( o1[0] >= BOX_SIZE * 0.5 ){
	         o1[0] =  o1[0] - BOX_SIZE;
	}

	// Check coordinate in y-direction
	if( o1[1] < -BOX_SIZE * 0.5 ){
	        o1[1] = o1[1] + BOX_SIZE;
	}
	else if(o1[1] >= BOX_SIZE * 0.5 ){
	        o1[1] = o1[1] - BOX_SIZE;
	}
	// Check coordinate in z-direction
	if(o1[2]< -BOX_SIZE * 0.5 ){
	        o1[2] = o1[2] + BOX_SIZE;
	}
	else if(o1[2] >= BOX_SIZE * 0.5 ){
	        o1[2] = o1[2]- BOX_SIZE;
	}

	if( o2[0] < -BOX_SIZE * 0.5 ){
		 o2[0] =  o2[0] + BOX_SIZE;
	}
	else if( o2[0] >= BOX_SIZE * 0.5 ){
	         o2[0] =  o2[0] - BOX_SIZE;
	}

	// Check coordinate in y-direction
	if( o2[1] < -BOX_SIZE * 0.5 ){
	        o2[1] = o2[1] + BOX_SIZE;
	}
	else if(o2[1] >= BOX_SIZE * 0.5 ){
	        o2[1] = o2[1] - BOX_SIZE;
	}
	// Check coordinate in z-direction
	if(o2[2]< -BOX_SIZE * 0.5 ){
	        o2[2] = o2[2] + BOX_SIZE;
	}
	else if(o2[2] >= BOX_SIZE * 0.5 ){
	        o2[2] = o2[2]- BOX_SIZE;
	}
*/
	// Compute the distance vector in the x-direction
	distance_vector[0] = o1[0] - o1[0];

	// Distance vector must also obey this convention
	// Check distance vector in x-direction
	if(distance_vector[0] > BOX_SIZE * 0.5){
		distance_vector[0] = distance_vector[0] - BOX_SIZE;
	}
	else if(distance_vector[0] <= -BOX_SIZE * 0.5){
		distance_vector[0] = distance_vector[0] + BOX_SIZE;
	}

	// Compute the distance vector in the y-direction
	distance_vector[1] = o1[1] - o2[1];

	// Check distance vector in y-direction
	if(distance_vector[1] > BOX_SIZE * 0.5){
	        distance_vector[1] = distance_vector[1] - BOX_SIZE;
	}
	else if(distance_vector[1] <= -BOX_SIZE * 0.5){
	        distance_vector[1] = distance_vector[1] + BOX_SIZE;
	}
	
	// Compute the distance vector in the z-direction
	distance_vector[2] = o1[2] - o2[2];

	// Check distance vector in z-direction
	if(distance_vector[2] > BOX_SIZE * 0.5){
	        distance_vector[2] = distance_vector[2] - BOX_SIZE;
	}
	else if(distance_vector[2] <= -BOX_SIZE * 0.5){
	        distance_vector[2] = distance_vector[2] + BOX_SIZE;
	}

	// Compute the distance between the two atoms and return this value
	distance_magnitude = sqrt(pow(distance_vector[0],2) + pow(distance_vector[1],2) + pow(distance_vector[2],2));	
	return distance_magnitude;
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
	printf("Oxygen ID IN FRAME: %d\n", oxygen_id);
	printf("HYDRONIUM HYDROGEN IDS: %d %d %d\n", hydrogen_index[0], hydrogen_index[1], hydrogen_index[2]);
}

void assign_closest_waters(int idx, int water_idx[], struct BookKeeping *closestOs){
	closestOs -> o[0] = water_idx[0];
	closestOs -> o[1] = water_idx[1];
	closestOs -> o[2] = water_idx[2];
	printf("Three Closest Oxygens to ID %d: [%d %d %d].\n", idx, water_idx[0], water_idx[1],water_idx[2]);

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

	float distance; // Current distance that is being tested
	float smallest_distances[] = {16,17,18}; // List containing shortest hydrogens to oxygen
    int hydrogen_idx[] = {-1,-1,-1}; // Index list corresponding to smallest_distances[]

	/*
	Loop over all of the atoms in a given frame and look for hydrogen atoms. Once we find a hydrogen atom,
	we will test to see if the hydrogen belongs to the hydronium oxygen.
	*/
	for(int i = 1; i < NUM_OF_ATOMS+1; i++){
		if(analysis_frame->atom_name[i] == 'H'){
			distance = compute_distance(oxygen_id, i, analysis_frame); //Compute the distance

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
	float distance;
	float distances[] = {16,17};
	int indexe[] = {-1,-1};
	int k = -1;
	for(int i = 0; i < NUM_OF_ATOMS; i++){
		if(analysis_frame->atom_name[i] == 'O' && i != oxygen_id){
			k++;
			float distances[] = {16,17}; //dup?
			int indexe[] = {-1,-1};
			for(int j = 0; j < NUM_OF_ATOMS; j++){
				if(analysis_frame->atom_name[j] == 'H'){
					distance = compute_distance(i, j, analysis_frame);
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

void identify_first_shell(struct Stack *stack, int idx_to_search, struct BookKeeping *closestOs){
	// Test distances
	// We seek the three closest oxygens to our hydronium ion
	// Is this in eigen/zundel form? long lived pair = 1-2ps
	float distance;
	float distanceh1;
	float distanceh2;
	float distanceh3;
	float smallest_distances[] = {16,17,18};
	int water_idx[] = {-1,-1,-1};
	int k = -1;

	//Iterate through all the possible water molecules. Remember that we have 128 oxygens in our system, but one is the hydronium ion
	for(int i = 0; i < 129; i++){
			k++;
			if(k != idx_to_search && idx_to_search != 0){
				distanceh1 = compute_oxygen_distances(stack->water[idx_to_search].h1, stack->water[i].o1);
				distanceh2 = compute_oxygen_distances(stack->water[idx_to_search].h2, stack->water[i].o1);
				if(distanceh1 < distanceh2){
					distance = distanceh1;
				}
				else if(distanceh2 < distanceh1){
					distance = distanceh2;
				}
				
				
				if(distance < smallest_distances[2] && distance > smallest_distances[1]){
					smallest_distances[2] = distance; 
					water_idx[2] = i;
				}

				else if(distance < smallest_distances[1] && distance > smallest_distances[0]){
					smallest_distances[2] = smallest_distances[1]; 
					smallest_distances[1] = distance; 
					water_idx[2] = water_idx[1];
					water_idx[1] = i;
				}

				else if(distance < smallest_distances[0]){
					smallest_distances[2] = smallest_distances[1]; 
					smallest_distances[1] = smallest_distances[0]; 
					smallest_distances[0] = distance;
					water_idx[2] = water_idx[1];
					water_idx[1] = water_idx[0];	
					water_idx[0] = i;
				}
			}
			else if(k != idx_to_search && idx_to_search == 0){
				distanceh1 = compute_oxygen_distances(stack->ion.h1, stack->water[i].o1);
				distanceh2 = compute_oxygen_distances(stack->ion.h2, stack->water[i].o1);
				distanceh3 = compute_oxygen_distances(stack->ion.h3, stack->water[i].o1);
				if(distanceh1 < distanceh2 && distanceh1 < distanceh3){
					distance = distanceh1;
				}
				else if(distanceh2 < distanceh1 && distanceh2 < distanceh3){
					distance = distanceh2;
				}
				else if(distanceh3 < distanceh1 && distanceh3 < distanceh2){
					distance = distanceh3;
				}
				if(distance < smallest_distances[2] && distance > smallest_distances[1]){
					smallest_distances[2] = distance; 
					water_idx[2] = i;
				}

				else if(distance < smallest_distances[1] && distance > smallest_distances[0]){
					smallest_distances[2] = smallest_distances[1]; 
					smallest_distances[1] = distance; 
					water_idx[2] = water_idx[1];
					water_idx[1] = i;
				}

				else if(distance < smallest_distances[0]){
					smallest_distances[2] = smallest_distances[1]; 
					smallest_distances[1] = smallest_distances[0]; 
					smallest_distances[0] = distance;
					water_idx[2] = water_idx[1];
					water_idx[1] = water_idx[0];	
					water_idx[0] = i;
				}
			}
	}
	assign_closest_waters(idx_to_search, water_idx, closestOs);
	// Build water molecules
	//assign_waters(k, i, water_idx, analysis_frame, waters);

}

void frame_to_stack(struct Stack *stack, struct H3O *hydronium, struct H2O *waters){
	stack->ion = *hydronium;
	for(int i = 0; i < 128; i++){
		stack->water[i] = waters[i];
	}
}


void identify_index(struct Stack *last_stack, struct Stack *current_stack, struct H3O *hydronium, struct H2O *waters){
	current_stack->ion = *hydronium;
	float checker;
	float smallest_distance = BOX_SIZE;
	int smallest_distance_idx;
	for(int i = 0; i < 128; i++){
		smallest_distance = BOX_SIZE;
		smallest_distance_idx = 0;
		for(int j = 0; j < 128; j++){
			checker = compute_oxygen_distances(current_stack->water[i].o1, last_stack->water[j].o1);
			if (checker < smallest_distance){
				smallest_distance = checker;
				smallest_distance_idx = j;
			}
		}
	current_stack->water[i] = waters[smallest_distance_idx];
	}
	printf("OXYGEN ID IN STACK: %d\n", current_stack->ion.index_o1);
}

int main(int argc, char* argv[]){

	/* Checking for errors */

	if(argc != 3){
		printf("\n------------------------------------------------\n\nWrong number of inputs. \nThe correct usage of this program is as follows.\n\n ./index_waters begin_frame end_frame.\n\n------------------------------------------------\n");
		exit(0);
	}

	int begin_frame;
	int end_frame;

	/* Setting the frames to be reading */

	sscanf (argv[1],"%d",&begin_frame);
	sscanf (argv[2],"%d",&end_frame);
	

	/* Checking for more errors */
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
	
	// Frames to read in with this code
	int frame_number = end_frame-begin_frame;

	clock_t begin; // Variable to time the entire code segment
	clock_t one_iteration; // Variable to see how long it takes to run one frame
	struct hydronium_idx ids; // Struct that will store the hydronium oxygen ids in an array
	read_ids(&ids); // Read in the hydronium oxygen ids
	
	FILE *fp;
	FILE *filewriter;
	FILE *framewriter;
	char *token;
	char buffer[200];
	
	fp = fopen("../water_scan-pos-1.xyz", "r");
	filewriter = fopen("../output.txt", "w");
	framewriter = fopen("../frame_change.txt", "w");
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

	int i; //Counter for current atom number in the frame
	struct Frame *current_frame =  malloc(sizeof(struct Frame));  // Frame to store coordinates	

	struct Stack *current_stack = malloc(sizeof(struct Stack)); 
	struct Stack *last_stack = malloc(sizeof(struct Stack));

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
		current_stack->ion.ion_index = (current_stack->ion.ion_index-1)/3;

		struct BookKeeping *keep = malloc(sizeof(struct BookKeeping));
		struct BookKeeping *keep2 = malloc(sizeof(struct BookKeeping));
		struct BookKeeping *keep3 = malloc(sizeof(struct BookKeeping));
		struct BookKeeping *keep4 = malloc(sizeof(struct BookKeeping));
		struct BookKeeping *old_keep = malloc(sizeof(struct BookKeeping));
		identify_first_shell(current_stack, current_stack->ion.ion_index, keep); // Identify all of the water molecules
		identify_first_shell(current_stack,keep->o[0], keep2);
		identify_first_shell(current_stack,keep->o[1], keep3);
		identify_first_shell(current_stack,keep->o[2], keep4);
		identify_first_shell(last_stack, last_stack->ion.ion_index, old_keep); // Identify all of the water molecules

		one_iteration = clock() - one_iteration; // Stop timing the iteration
		
		
		if(last_stack->ion.ion_index != current_stack->ion.ion_index || old_keep->o[0]  != keep->o[0] || old_keep->o[1]  != keep->o[1] || old_keep->o[2]  != keep->o[2]){
		fprintf(framewriter, "%d\n", k);
		}
		fprintf(filewriter, "Frame %d\nH3O ID: %d\n", k, current_stack->ion.ion_index);
		fprintf(filewriter, "Closest Water Molecule IDS: [%d %d %d]\n\n", keep->o[0], keep->o[1], keep->o[2]);
		printf("Finished processing frame %d.\nProcessing took %f seconds.\n", k, ((double)one_iteration)/CLOCKS_PER_SEC);
		printf("======================================================\n\n");
	
		
	}
	
	begin = clock() - begin; //Stop timing the code
	
	printf("Processing frames %d-%d took %f seconds.\n", begin_frame, end_frame, ((double)begin)/CLOCKS_PER_SEC);
	fclose(fp); //Close file
	fclose(filewriter);
	fclose(framewriter);
}
