#include "lodepng.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define bytes_per_pixel 4

/*
	Rectification Performance
*/
void rectify_process(char* input_filename, char* output_filename, int num_of_threads)
{
	//declare variables to load images
	unsigned error;
	unsigned char *image, *new_image;
	unsigned width, height;
	error = lodepng_decode32_file(&image, &width, &height, input_filename);
	if(error) printf("error %u: %s\n", error, lodepng_error_text(error));
	new_image = malloc(width * height * bytes_per_pixel * sizeof(unsigned char));

	//start timer
 	double start = omp_get_wtime( );
	//rectify image
	int i = 0, j = 0;
	#pragma omp parallel for num_threads(num_of_threads) collapse(2) private(j)
	for (i = 0; i < height; i++)
	{
		for (j = 0; j < width; j++)
		{
			//get pixel pointer of images
			int idx = bytes_per_pixel * width * i + bytes_per_pixel * j;
			//compare the pixel values to 127, and set equal to 127 if lower
			signed int value_R = (int)image[idx] -127;
			signed int value_G = (int)image[idx + 1] - 127;
			signed int value_B = (int)image[idx + 2] - 127;
			value_R = (value_R >= 0) ? value_R : 0;
			value_G = (value_G >= 0) ? value_G : 0;
			value_B = (value_B >= 0) ? value_B : 0;
			new_image[idx] = (unsigned char)(value_R + 127);
			new_image[idx + 1] = (unsigned char)(value_G + 127);
			new_image[idx + 2] = (unsigned char)(value_B + 127);
			new_image[idx + 3] = image[idx + 3];
		}
	}
	// end the timer
	double end = omp_get_wtime();
	// stores the difference in diff
	double diff = end - start;
	printf("Time: %f seconds\n", diff);

	//load the new image produced
	lodepng_encode32_file(output_filename, new_image, width, height);

	free(image);
	free(new_image);
}

int main(int argc, char *argv[]) {

	if (argc < 4)
	{
		printf("Error: missing arguments\nRun the command as following:\n./rectify <name of input png> <name of output png> <# threads>\n");
		return 0;
	}
	//retrieve each arguments
	char* input_filename = argv[1];
	char* output_filename = argv[2];
	int num_of_threads = atoi(argv[3]);

	rectify_process(input_filename, output_filename, num_of_threads);

    // done
    return 0;
}
