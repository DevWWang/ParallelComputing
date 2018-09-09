#include "lodepng.h"
#include "wm.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define bytes_per_pixel 4
#define matrices_size 3

/*
	Convolution Performance
*/
void convolve_process(char* input_filename, char* output_filename, int num_of_threads)
{
	//declare variables to load images
 	unsigned error;
 	unsigned char *image, *new_image;
 	unsigned width, height;
 	unsigned new_width, new_height;
 	error = lodepng_decode32_file(&image, &width, &height, input_filename);
 	if(error) printf("error %u: %s\n", error, lodepng_error_text(error));
 	//both width and height of new images are 2 pixels less than the dimension of the original image
 	new_width = (unsigned)((int)width - 2);
 	new_height = (unsigned)((int)height - 2);
 	new_image = malloc(new_width * new_height * bytes_per_pixel * sizeof(unsigned char));

 	//start timer
 	double start = omp_get_wtime( );
 	// convolve image
 	int i = 0, j = 0;
 	#pragma omp parallel for num_threads(num_of_threads) private(j) collapse(2)
 	for (i = 1; i < (height - 1); i++)
 	{
 		for (j = 1; j < (width - 1); j++)
 		{
 			float sum_R =0.0, sum_G = 0.0, sum_B = 0.0;
 			//compute the corresponding pixel in the output image using a weighted sum of the input pixel and its neighbors
 			for (int ii = 0; ii < matrices_size; ii++)
 			{
 				for (int jj = 0; jj < matrices_size; jj++)
 				{
 					//get pixel pointer of input image
 					int idx = bytes_per_pixel * (i + ii - 1) * width + bytes_per_pixel * (j + jj - 1);
 					sum_R += image[idx] * w[ii][jj];
 					sum_G += image[idx + 1] * w[ii][jj];
 					sum_B += image[idx + 2] * w[ii][jj];
 				}
 			}
 			//clamp for R
 			sum_R = (sum_R < 0) ? 0 : sum_R;
 			sum_R = (sum_R > 255) ? 255 : sum_R;
 			//clamp for G
 			sum_G = (sum_G < 0) ? 0 : sum_G;
 			sum_G = (sum_G > 255) ? 255 : sum_G;
 			//clamp for B
 			sum_B = (sum_B < 0) ? 0 : sum_B;
 			sum_B = (sum_B > 255) ? 255 : sum_B;
 			//get pixel pointer of output image
 			int pointer = bytes_per_pixel * new_width * (i - 1) + bytes_per_pixel * (j - 1);
 			//restore the pixel for new images
 			new_image[pointer] = (unsigned char)sum_R;
			new_image[pointer + 1] = (unsigned char)sum_G;
			new_image[pointer + 2] = (unsigned char)sum_B;
			new_image[pointer + 3] = 255;
			//clear the variables
			sum_R = 0;
			sum_G = 0;
			sum_B = 0;
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
		printf("Error: missing arguments\nRun the command as following:\n./convolve <name of input png> <name of output png> <# threads>\n");
		return 0;
	}
	//retrieve each arguments
	char* input_filename = argv[1];
	char* output_filename = argv[2];
	int num_of_threads = atoi(argv[3]);

	convolve_process(input_filename, output_filename, num_of_threads);

    // done
    return 0;
}
