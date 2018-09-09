#include "lodepng.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <time.h>

#define bytes_per_pixel 4
#define max(x, y)  ((x) > (y) ? (x) : (y))

/*
	Pooling Performance
*/
void pool_process(char* input_filename, char* output_filename, int num_of_threads)
{
	//declare variables to load images
	unsigned error;
	unsigned char *image, *new_image;
	unsigned width, height;
	unsigned new_width, new_height;
	error = lodepng_decode32_file(&image, &width, &height, input_filename);
	if(error) printf("error %u: %s\n", error, lodepng_error_text(error));
	//both width and height become half of their original value
	new_width = (unsigned)((int)width / 2);
 	new_height = (unsigned)((int)height / 2);
	new_image = malloc(new_width * new_height * bytes_per_pixel * sizeof(unsigned char));

	//start timer
 	double start = omp_get_wtime( );
	//pool image
	int i = 0, j = 0;
 	#pragma omp parallel for num_threads(num_of_threads) private(j) collapse(2)
 	for (i = 0; i < height; i+=2)
 	{
 		for (j = 0; j < width; j+=2)
 		{
 			float max = 0.0;
 			//get pixel pointer of input image
 			/*
 			indices of a 2x2 squares
 			---------------------------------
 			|	1 (i,j)		|	2 (i,j+1)	|
 			---------------------------------
 			|	3 (i+1,j+1)	|	4 (i+1,j+1)	|
 			---------------------------------
 			*/
 			int idx_1 = bytes_per_pixel * width * i + bytes_per_pixel * j;
 			int idx_2 = bytes_per_pixel * width * i + bytes_per_pixel * (j + 1);
 			int idx_3 = bytes_per_pixel * width * (i + 1) + bytes_per_pixel * j;
 			int idx_4 = bytes_per_pixel * width * (i + 1) + bytes_per_pixel * (j + 1);
 			//get pixel pointer of output image
 			int pointer = bytes_per_pixel * new_width * (i / 2) + bytes_per_pixel * (j / 2);
 			//find the maximum value by comparing each byte of each pixel within the square
 			for (int rgba = 0; rgba < bytes_per_pixel; rgba++)
 			{
 				max = max(max(image[idx_1 + rgba], image[idx_2 + rgba]),max(image[idx_3 + rgba], image[idx_4 + rgba]));
 				new_image[pointer + rgba] = (unsigned char)max;
 			}
 		}
 	}
 	// end the timer
	double end = omp_get_wtime();
	// stores the difference in diff
	double diff = end - start;
	printf("Time: %f seconds\n", diff);

 	//load the new image produced
	lodepng_encode32_file(output_filename, new_image, new_width, new_height);

	free(image);
	free(new_image);
}

int main(int argc, char *argv[]) {

	if (argc < 4)
	{
		printf("Error: missing arguments\nRun the command as following:\n./pool <name of input png> <name of output png> <# threads>\n");
		return 0;
	}
	//retrieve each arguments
	char* input_filename = argv[1];
	char* output_filename = argv[2];
	int num_of_threads = atoi(argv[3]);

	pool_process(input_filename, output_filename, num_of_threads);

    // done
    return 0;
}
