#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>

#define NUM_OF_FRAMES 90677
#define NUM_OF_ATOMS 385
#define BOX_SIZE 15.8122

/*
Stores the ids of the hydronium ion from the id file
*/
struct hydronium_idx{
	int i[NUM_OF_FRAMES];
	int idx[NUM_OF_FRAMES];
};

/*
Stores an instance of H3O molecule
*/
struct H3O{
	float o1[3];
	float h1[3];
	float h2[3];
	float h3[3];
	int index_h1;
	int index_h2;
	int index_h3;
};

/*
Stores an instance of H2O molecule
*/
struct H2O{
	float o1[3];
	float h1[3];
	float h2[3];
	int index_h1;
	int index_h2;
};

/*
Stores the coordinates of the simulation timesteps
*/
struct Frame{
	char atom_name[NUM_OF_ATOMS];
	float x[NUM_OF_ATOMS];
	float y[NUM_OF_ATOMS];
	float z[NUM_OF_ATOMS];
};


void print_array(int* array_to_print, int SIZE){
	for(int i = 0; i < SIZE; i++){
		printf("%d -> %d \n", i, array_to_print[i]);
	}
}

void print_frame(struct Frame frame_to_print, int SIZE){
	for(int i = 0; i < SIZE; i++){
		printf("%d -> %c: %f %f %f\n", i, frame_to_print.atom_name[i], frame_to_print.x[i], frame_to_print.y[i], frame_to_print.z[i]);
	}
}


/*
Function to read the id file and store the hydronium ion index at each timestep
*/
void read_ids(struct hydronium_idx *arr){
	FILE *fp2;
    int i = 0;
	int FILELENGTH = NUM_OF_FRAMES;
	fp2 = fopen("id_hydronium.dat", "r");
	while(i < FILELENGTH){
		fscanf(fp2, "%d %d", &arr->i[i], &arr->idx[i]);
		i++;
	      }
	fclose(fp2);
}

/*
Using the minimum image convention and pbc conditions, compute the distance between two atoms.
Generally, should be used to compute distances between an oxygen atom and hydrogen atoms.
*/
float compute_distance(int oxygen_id, int i, struct Frame *analysis_frame){
	float distance_vector[3];
	float distance_magnitude;
	if( analysis_frame->x[i] < -BOX_SIZE * 0.5 ){
		analysis_frame -> x[i] = analysis_frame->x[i] + BOX_SIZE;
	}
	else if(analysis_frame->x[i] >= BOX_SIZE * 0.5 ){
	        analysis_frame->x[i] = analysis_frame->x[i] - BOX_SIZE;
	}

	if( analysis_frame->y[i] < -BOX_SIZE * 0.5 ){
	        analysis_frame->y[i] = analysis_frame->y[i] + BOX_SIZE;
	}
	else if(analysis_frame->y[i] >= BOX_SIZE * 0.5 ){
	        analysis_frame->y[i] = analysis_frame->y[i] - BOX_SIZE;
	}

	if( analysis_frame->z[i] < -BOX_SIZE * 0.5 ){
	        analysis_frame->z[i] = analysis_frame->z[i] + BOX_SIZE;
	}
	else if(analysis_frame->z[i] >= BOX_SIZE * 0.5 ){
	        analysis_frame->z[i] = analysis_frame->z[i] - BOX_SIZE;
	}
	distance_vector[0] = analysis_frame->x[i] - analysis_frame->x[oxygen_id];
	if(distance_vector[0] > BOX_SIZE * 0.5){
		distance_vector[0] = distance_vector[0] - BOX_SIZE;
	}
	else if(distance_vector[0] <= -BOX_SIZE * 0.5){
		distance_vector[0] = distance_vector[0] + BOX_SIZE;
	}
	
	distance_vector[1] = analysis_frame->y[i] - analysis_frame->y[oxygen_id];
	if(distance_vector[1] > BOX_SIZE * 0.5){
	        distance_vector[1] = distance_vector[1] - BOX_SIZE;
	}
	else if(distance_vector[1] <= -BOX_SIZE * 0.5){
	        distance_vector[1] = distance_vector[1] + BOX_SIZE;
	}
	
	distance_vector[2] = analysis_frame->z[i] - analysis_frame->z[oxygen_id];	
	if(distance_vector[2] > BOX_SIZE * 0.5){
	        distance_vector[2] = distance_vector[2] - BOX_SIZE;
	}
	else if(distance_vector[2] <= -BOX_SIZE * 0.5){
	        distance_vector[2] = distance_vector[2] + BOX_SIZE;
	}
	distance_magnitude = sqrt(pow(distance_vector[0],2) + pow(distance_vector[1],2) + pow(distance_vector[2],2));	
	return distance_magnitude;
}

/*
Assign the hydronium ion coordinates to the respective struct.
*/
void assign_ion(int oxygen_id, int* hydrogen_index, struct Frame *analysis_frame, struct H3O *hydronium){
	hydronium->o1[0] = analysis_frame->x[oxygen_id];
    hydronium->o1[1] = analysis_frame->y[oxygen_id];
    hydronium->o1[2] = analysis_frame->z[oxygen_id];
	hydronium->h1[0] = analysis_frame->x[hydrogen_index[0]];
	hydronium->h1[1] = analysis_frame->y[hydrogen_index[0]];
	hydronium->h1[2] = analysis_frame->z[hydrogen_index[0]];
	hydronium->h2[0] = analysis_frame->x[hydrogen_index[1]];
	hydronium->h2[1] = analysis_frame->y[hydrogen_index[1]];
	hydronium->h2[2] = analysis_frame->z[hydrogen_index[1]];
	hydronium->h3[0] = analysis_frame->x[hydrogen_index[2]];
	hydronium->h3[1] = analysis_frame->y[hydrogen_index[2]];
	hydronium->h3[2] = analysis_frame->z[hydrogen_index[2]];
    hydronium->index_h1 = hydrogen_index[0];
	hydronium->index_h2 = hydrogen_index[1];
	hydronium->index_h3 = hydrogen_index[2];
	printf("Oxygen ID: %d\n", oxygen_id);
	printf("%d %d %d\n", hydrogen_index[0], hydrogen_index[1], hydrogen_index[2]);
}

/*
Use the distance formula in order to find the three closest hydrogens to the oxygen molecule based off of
the provided id.
*/
void identify_ion(int oxygen_id, struct Frame *analysis_frame, struct H3O *hydronium){
	float distance;
	float distances[] = {16,17,18};
    int indexe[] = {-1,-1,-1};
	for(int i = 0; i < NUM_OF_ATOMS; i++){
		if(analysis_frame->atom_name[i] == 'H'){
			distance = compute_distance(oxygen_id, i, analysis_frame);
			if(distance < distances[2] && distance > distances[1]){
				distances[2] = distance; 
				indexe[2] = i;
			}
			else if(distance < distances[1] && distance > distances[0]){
				distances[2] = distances[1]; 
				distances[1] = distance; 
				indexe[2] = indexe[1];
				indexe[1] = i;
			}
			else if(distance < distances[0]){
				distances[2] = distances[1]; 
				distances[1] = distances[0]; 
				distances[0] = distance;
			    indexe[2] = indexe[1];
				indexe[1] = indexe[0];	
				indexe[0] = i;
			}
		}
	}
	
	assign_ion(oxygen_id, indexe, analysis_frame, hydronium);
	
}

/*
Use the distance formula in order to find the two closest hydrogens to each of the oxygen molecules that
are not the hydronium ion oxygen.
*/
void identify_waters(int oxygen_id, struct Frame *analysis_frame, struct H2O *waters){
		float distance;
		float distances[] = {16,17};
		int indexe[] = {-1,-1};
		int k = -1;
		for(int i = 0; i < NUM_OF_ATOMS; i++){
			if(analysis_frame->atom_name[i] == 'O' && i != oxygen_id){
				k++;
				float distances[] = {16,17};
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
				waters[k].o1[0] = analysis_frame->x[i];
				waters[k].o1[1] = analysis_frame->y[i];
				waters[k].o1[2] = analysis_frame->z[i];
				waters[k].h1[0] = analysis_frame->x[indexe[0]];
				waters[k].h1[1] = analysis_frame->y[indexe[0]];
				waters[k].h1[2] = analysis_frame->z[indexe[0]];
				waters[k].h2[0] = analysis_frame->x[indexe[1]];
				waters[k].h2[1] = analysis_frame->y[indexe[1]];
				waters[k].h2[2] = analysis_frame->z[indexe[1]];
				waters[k].index_h1 = indexe[0];
				waters[k].index_h2 = indexe[1];
			}
		
		//assign_ion(oxygen_id, indexe, analysis_frame, hydronium);
		}
		/*
		for(int i = 0; i < 126; i++){
			printf("Oxygen ID: %d\n", i);
			printf("%d %d\n", waters[i].index_h1 , waters[i].index_h2);
			printf("======================================================\n\n");		
		}*/
}



int main(int argc, char* argv[]){

	clock_t begin;
	clock_t one_iteration;
	struct hydronium_idx ids;
	read_ids(&ids);
	
	FILE *fp;
	int i = 0;
	char buffer[200];
	char *token;
	int begin_frame;
	int end_frame;
	sscanf (argv[1],"%d",&begin_frame);
	sscanf (argv[2],"%d",&end_frame);


	int frame_number = end_frame-begin_frame;
	fp = fopen("water_scan-pos-1.xyz", "r");
	begin= clock();
	if(begin_frame != 0){
		for(int garb = 0; garb < 387*begin_frame; garb++){
			fgets(buffer, 200, fp);
		}
	}
	
	for(int k = begin_frame; k < end_frame + 1; k++){
		struct Frame frame;
		i = 0;
		fgets(buffer, 100, fp);
		fgets(buffer, 100, fp);
		while((i < NUM_OF_ATOMS)){
			fgets(buffer, 100, fp);
            token = strtok(buffer, " \t");    
			
			frame.atom_name[i] = token[0];
			token = strtok(NULL, " \t");
			
			frame.x[i] = atof(token);
			token = strtok(NULL, " \t");
			
			frame.y[i] = atof(token);
			token = strtok(NULL, " \t");
			
			frame.z[i] = atof(token);
			i++;
		}
		printf("======================================================\n");
		printf("Finished reading frame %d\n", k);
		one_iteration = clock();
		int current_H3O = ids.idx[k];
 		struct H3O hydronium;
		struct H2O waters[127];
		identify_ion(current_H3O, &frame, &hydronium);
		identify_waters(current_H3O, &frame, waters);	
		one_iteration = clock() - one_iteration;
		printf("Finished processing frame %d.\nProcessing took %f seconds.\n", k, ((double)one_iteration)/CLOCKS_PER_SEC);
		printf("======================================================\n\n");
	}
	begin = clock() - begin;
	printf("Processing frames %d-%d took %f seconds.\n", begin_frame, end_frame, ((double)begin)/CLOCKS_PER_SEC);
	fclose(fp);
}
