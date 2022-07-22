#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

#define NUM_OF_FRAMES 50000
#define NUM_OF_ATOMS 385
#define BOX_SIZE 15.8122

struct hydronium_idx{
	int i[NUM_OF_FRAMES];
	int idx[NUM_OF_FRAMES];
};

struct H3O{
	float o1[3];
	float h1[3];
	float h2[3];
	float h3[3];
};

struct H2O{
	float o1[3];
	float h1[3];
	float h2[3];
};

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

	printf("Oxygen ID: %d\n", oxygen_id);
	printf("%d %d %d\n", hydrogen_index[0], hydrogen_index[1], hydrogen_index[2]);
	printf("======================================================\n\n");
}

void identify_ion(int oxygen_id, struct Frame *analysis_frame, struct H3O *hydronium){
	float distance;
	float distances[] = {16,17,18};
    int indexe[] = {-1,-1,-1};
	for(int i = 0; i < NUM_OF_ATOMS; i++){
		if(analysis_frame->atom_name[i] == 'H'){
			distance = compute_distance(oxygen_id, i, analysis_frame);
			if(distance < distances[2] && distance > distances[1]){
				distances[2] = distance; 
				printf("%d\n", i);
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


int main(){
	struct hydronium_idx ids;
	read_ids(&ids);
	
	FILE *fp;
	int i = 0;
	char buffer[200];
	char *token;
	int frame_number = NUM_OF_FRAMES;
	fp = fopen("water_scan-pos-1.xyz", "r");
	
	for(int k = 0; k < frame_number; k++){
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
		int current_H3O = ids.idx[k];
 		struct H3O hydronium;
		identify_ion(current_H3O, &frame, &hydronium);	
	}
	fclose(fp);
}
