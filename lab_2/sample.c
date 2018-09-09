#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define grid_size 4
#define rho 0.5
#define eta 2E-4
#define G 0.75

float** initialize(int size, int hit) {
	float **grid = (float **)malloc(grid_size * sizeof(float*));
	int i;
	for (i = 0; i < grid_size; i++) {
		grid[i] = (float *)malloc(grid_size * sizeof(float));
	}
	if (hit == 1) {
		grid[grid_size / 2][grid_size / 2] = 1;
	}

	return grid;
}

int store_prev(int size, float **previous, float **current) {
	int i, j;
	for (i = 0; i < grid_size; i++) {
		for (j = 0; j < grid_size; j++) {
			previous[i][j] = current[i][j];
		}
	}
	return 0;
}

int print_grid(int size, float **grid) {
	int i, j;
	for (i = 0; i < size; i++) {
		for (j = 0; j < size; j++) {
			printf("(%d, %d): %.6f  ", i, j, grid[i][j]);
		}
		printf("\n");
	}
	return 0;
}

int print_output(float **grid) {
	printf("(%d, %d) %.6f\n", grid_size / 2, grid_size / 2, grid[grid_size / 2][grid_size / 2]);
}

int main(int argc, char *argv[]) {

	if(argc < 1)
    {
        printf("Missing argument for iteration\n");
        return 0;
    }
    int T = atoi(argv[1]);

	float **grid;
	//printf("Create An Empty Grid\n");
	grid = initialize(grid_size, 1);
	//print_grid(grid_size, grid);

	float **u = initialize(grid_size, 0);
	float **u1 = initialize(grid_size, 0);
	float **u2 = initialize(grid_size, 0);

	store_prev(grid_size, u1, grid);

	while (T --> 0) {
		/*
		printf("u:\n");
		print_grid(grid_size, u);
		printf("u1:\n");
		print_grid(grid_size, u1);
		printf("u2:\n");
		print_grid(grid_size, u2);
		*/
		//printf("Interior\n");
		for (int i = 1; i < grid_size - 1; i++) {
			for (int j = 1; j <grid_size - 1; j++) {
				u[i][j] = (rho*(u1[i-1][j] + u1[i+1][j] + u1[i][j-1] + u1[i][j+1] - 4 * u1[i][j])
							+ 2 * u1[i][j] - (1 - eta) * u2[i][j]) / (1 + eta);
				/*printf("Index: %d\n", (grid_size*i + j));
				printf("u1[i-1][j]=%.6f\n", u1[i-1][j]);
				printf("u1[i+1][j]=%.6f\n", u1[i+1][j]);
				printf("u1[i][j-1]=%.6f\n", u1[i][j-1]);
				printf("u1[i][j+1]=%.6f\n", u1[i][j+1]);
				printf("u1[i][j]=%.6f\n", u1[i][j]);
				printf("u2[i][j]=%.6f\n", u2[i][j]);
				*/
			}
		}
		//printf("Done Interior\n");
		//print_grid(grid_size, u);

		//printf("Side\n");
		for (int i = 1; i < grid_size; i++) {
			u[0][i] = G * u [1][i];
			u[grid_size - 1][i] = G * u[grid_size - 2][i];
			u[i][0] = G * u[i][1];
			u[i][grid_size - 1] = G * u[i][grid_size - 2];
		}
		//printf("Done Side\n");
		//print_grid(grid_size, u);

		//printf("Corner\n");
		u[0][0] = G * u[1][0];
		u[grid_size - 1][0] = G * u[grid_size - 2][0];
		u[0][grid_size - 1] = G * u[0][grid_size - 2];
		u[grid_size - 1][grid_size - 1] = G * u[grid_size - 1][grid_size - 2];
		//printf("Done Corner\n");

		store_prev(grid_size, u2, u1);
		store_prev(grid_size, u1, u);

		printf("u:\n");
		print_grid(grid_size, u);
		/*
		printf("u1:\n");
		print_grid(grid_size, u1);
		printf("u2:\n");
		print_grid(grid_size, u2);
		*/
		print_output(u);
	}

}
