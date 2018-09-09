#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#define grid_size 512
#define total_entries (grid_size * grid_size)

#define top 0
#define right 1
#define bottom 2
#define left 3

#define rho 0.5
#define eta 2E-4
#define G 0.75

double  startTime, finishTime;

struct element_s {
    float u;
    float u1;
    float u2;
};

// Create a grid
int createMatrix(struct element_s *matrix, int size, int hit) {
    //struct element_s matrix[total_entries] = {0};
    for (int i = 0; i < size; i++) {
        matrix[i].u = 0;
        matrix[i].u1 = 0;
        matrix[i].u2 = 0;
    }
    if (hit == 1) {
        int hit_idx = grid_size * (grid_size / 2) + (grid_size / 2);
        matrix[hit_idx].u = 1;
        matrix[hit_idx].u1 = 1;
    }
    return 0;
}

int get_array_idx(int i, int j) {
	return (grid_size * i + j);
}

int idx_to_ij(int idx, int position) {
    int i = idx / grid_size;
    int j = idx % grid_size;

    switch(position){
        case 0:
            return (grid_size * (i - 1) + j);

        case 1:
            return (grid_size * i + (j + 1));

        case 2:
            return (grid_size * (i + 1) + j);

        case 3:
            return (grid_size * i + (j - 1));

        default:
            return 0;
    }
}

int store_prev(struct element_s *entry, float new_value) {
    entry->u2 = entry->u1;
    entry->u1 = entry->u;
    entry->u = new_value;

    return 0;
}

int printMatrix(struct element_s *matrix, int u_idx) {
    for (int i = 0; i < grid_size; i++) {
        for (int j = 0; j < grid_size; j++) {
        	if (u_idx == 0) {
        		printf("(%d, %d) = %.6f\t", i, j, matrix[grid_size * i + j].u);
        	}
            else if (u_idx == 1) {
            	printf("(%d, %d) = %.6f\t", i, j, matrix[grid_size * i + j].u1);
            }
            else if (u_idx == 2) {
            	printf("(%d, %d) = %.6f\t", i, j, matrix[grid_size * i + j].u2);
            }
        }
        printf("\n");
    }
    printf("\n");
    return 0;
}

int printOutput(struct element_s *matrix) {
    int output_idx = grid_size * (grid_size / 2) + (grid_size / 2);
    //printf("(%d, %d) %.6f\n", grid_size / 2, grid_size / 2, matrix[output_idx].u);
    printf("%.6f\n", matrix[output_idx].u);

    return 0;
}

int printLocal(struct element_s element, int rank) {
    printf("ID [%d] has data:\n", rank);
    printf("u = %.6f\t", element.u);
    printf("u1 = %.6f\t", element.u1);
    printf("u2 = %.6f\n", element.u2);

    return 0;
}

// Create a MPI datatype of an element.
int mpi_element_init(MPI_Datatype *mpi_element) {
    // instance of structure
    //struct element_s element;
    // temporary loop indexer
    //int i = 0;
    // number of blocks in the struct
    int count = 3;
    // set up 3 blocks
    int blocks[3] = {1, 1, 1};
    // element internal types
    MPI_Datatype types[3] = {
        MPI_FLOAT,
        MPI_FLOAT,
        MPI_FLOAT
    };
    // internal displacements
    MPI_Aint dis[4] = {
        offsetof(struct element_s, u),
        offsetof(struct element_s, u1),
        offsetof(struct element_s, u2)
    };

    MPI_Type_create_struct(count, blocks, dis, types, mpi_element);
    MPI_Type_commit(mpi_element);

    return(EXIT_SUCCESS);
}

int main(int argc, char **argv) {

	if(argc < 1)
    {
        printf("Missing argument for iteration\n");
        return 0;
    }
    int T = atoi(argv[1]);

    MPI_Datatype mpi_node;

    // Initialize the MPI environment
    MPI_Init(&argc, &argv);
    startTime = MPI_Wtime();

    // Get the number of processes
    int num_proc;
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);

    // Get the id of the process
    int id;
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    // Create an MPI element datatype
    mpi_element_init(&mpi_node);

    int chunk_size = total_entries / num_proc;
    int first_idx = id * chunk_size;
    int end_idx = (id == num_proc - 1) ? total_entries : first_idx + chunk_size;

    //create a local array for each processor
    int local_length = end_idx - first_idx;

    struct element_s sample[total_entries] = {0};
    struct element_s *local_data; 
    createMatrix(sample, total_entries, 1);
    local_data = malloc(local_length * sizeof(struct element_s));

    while (T --> 0) {
    	MPI_Bcast(&chunk_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
        //Split the work to each processor
        MPI_Scatter(sample, local_length, mpi_node, local_data, local_length, mpi_node, 0, MPI_COMM_WORLD);

        for (int local_idx = 0; local_idx < local_length; local_idx++) {

        	int entry_idx = id * chunk_size + local_idx;
        	if (entry_idx % grid_size == 0 || entry_idx % grid_size == (grid_size - 1) ||
                        entry_idx < get_array_idx(0, grid_size - 1) || entry_idx > get_array_idx(grid_size - 2, grid_size - 1)) {
        	    //ignore elements located on the edges
        	}
        	//interior
        	else {
        	    float temp;
        	    temp = (rho*(sample[idx_to_ij(entry_idx, top)].u1 + sample[idx_to_ij(entry_idx, bottom)].u1 + sample[idx_to_ij(entry_idx, left)].u1 + sample[idx_to_ij(entry_idx, right)].u1 - 4 * local_data[local_idx].u1)
        	                + 2 * local_data[local_idx].u1 - (1 - eta) * local_data[local_idx].u2) / (1 + eta);
        	    //printLocal(local_data[local_idx], id);
                store_prev(&local_data[local_idx], temp);
        	}
        }

        MPI_Barrier(MPI_COMM_WORLD);
        //Collect the results from each processor
        MPI_Gather(local_data, local_length, mpi_node, sample, local_length, mpi_node, 0, MPI_COMM_WORLD);
        MPI_Bcast(sample, total_entries, mpi_node, 0, MPI_COMM_WORLD);
        //Split the work to each processor
        MPI_Scatter(sample, local_length, mpi_node, local_data, local_length, mpi_node, 0, MPI_COMM_WORLD);

        for (int local_idx = 0; local_idx < local_length; local_idx++) {
            int entry_idx = id * chunk_size + local_idx;
        	if (entry_idx == get_array_idx(0, 0) || entry_idx == get_array_idx(0, grid_size - 1) ||
                        entry_idx == get_array_idx(grid_size - 1, 0) || entry_idx == get_array_idx(grid_size - 1, grid_size - 1)) {
        	    //ignore elements located at the corner
        	}
        	//top edge
        	else if (entry_idx < get_array_idx(0, grid_size)) {
        	    float temp;
        	    temp = G * sample[entry_idx + grid_size].u;
        	    store_prev(&local_data[local_idx], temp);
        	}
        	//bottom edge
        	else if (entry_idx > get_array_idx(grid_size - 2, grid_size - 1)) {
        	    float temp;
        	    temp = G * sample[entry_idx - grid_size].u;
        	    store_prev(&local_data[local_idx], temp);
        	}
        	//left edge
        	else if (entry_idx % grid_size == 0) {
        	    float temp;
        	    temp = G * sample[entry_idx + 1].u;
        	    store_prev(&local_data[local_idx], temp);
        	}
        	//right edge
        	else if (entry_idx % grid_size == (grid_size - 1)) {
        	    float temp;
        	    temp = G * sample[entry_idx - 1].u;
        	    store_prev(&local_data[local_idx], temp);
        	}
        }
        MPI_Barrier(MPI_COMM_WORLD);
        //Collect the results from each processor
        MPI_Gather(local_data, local_length, mpi_node, sample, local_length, mpi_node, 0, MPI_COMM_WORLD);
        MPI_Bcast(sample, total_entries, mpi_node, 0, MPI_COMM_WORLD);
        //Split the work to each processor
        MPI_Scatter(sample, local_length, mpi_node, local_data, local_length, mpi_node, 0, MPI_COMM_WORLD);

        for (int local_idx = 0; local_idx < local_length; local_idx++) {
            int entry_idx = id * chunk_size + local_idx;
            //top left corner
        	if (entry_idx == get_array_idx(0, 0)) {
        	    float temp;
        	    temp = G * sample[entry_idx + grid_size].u;
        	    store_prev(&local_data[local_idx], temp);
        	}
        	//top right and bottom right
        	else if (entry_idx == get_array_idx(0, grid_size - 1) || entry_idx == get_array_idx(grid_size - 1, grid_size - 1)) {
        	    float temp;
        	    temp = G * sample[entry_idx - 1].u;
        	    store_prev(&local_data[local_idx], temp);
        	}
        	//bottom left
        	else if (entry_idx == get_array_idx(grid_size - 1, 0)) {
        	    float temp;
        	    temp = G * sample[entry_idx - grid_size].u;
        	    store_prev(&local_data[local_idx], temp);
        	}
        }
        MPI_Barrier(MPI_COMM_WORLD);

        //Collect the results from each processor
        MPI_Gather(local_data, local_length, mpi_node, sample, local_length, mpi_node, 0, MPI_COMM_WORLD);

        for (int i = 0; i < total_entries; i++){
            store_prev(&sample[i], sample[i].u);
        }

        if (id == 0) {
            //printMatrix(sample, 0);
            printOutput(sample);
        }

    }

    finishTime = MPI_Wtime();
    // Finalize the MPI environment.
    MPI_Finalize();

    printf("Parallel Elapsed time: %f seconds\n", (finishTime - startTime));
}