#include <math.h>
#include <inttypes.h>
#include <algorithm>
#include <unordered_set>



#define IB_Z 0
#define IB_Y 1
#define IB_X 2



void CppDownsampleMapping(const char *prefix, int64_t *segmentation, float input_resolution[3], int64_t output_resolution[3], int64_t input_grid_size[3])
{
    // get the number of entries
    int64_t input_nentries = input_grid_size[IB_Z] * input_grid_size[IB_Y] * input_grid_size[IB_X];

    // get downsample ratios
    float zdown = ((float) output_resolution[IB_Z]) / input_resolution[IB_Z];
    float ydown = ((float) output_resolution[IB_Y]) / input_resolution[IB_Y];
    float xdown = ((float) output_resolution[IB_X]) / input_resolution[IB_X];

    // get the output resolution size
    int64_t output_grid_size[3];
    output_grid_size[IB_Z] = (int64_t) ceil(input_grid_size[IB_Z] / zdown);
    output_grid_size[IB_Y] = (int64_t) ceil(input_grid_size[IB_Y] / ydown);
    output_grid_size[IB_X] = (int64_t) ceil(input_grid_size[IB_X] / xdown);
    int64_t output_sheet_size = output_grid_size[IB_Y] * output_grid_size[IB_X];
    int64_t output_row_size = output_grid_size[IB_X];

    int64_t max_segment = 0;
    for (int64_t iv = 0; iv < input_nentries; ++iv)
        if (segmentation[iv] > max_segment) max_segment = segmentation[iv];
    max_segment++;

    // create a set for each segment of downsampled locations
    std::unordered_set<int64_t> *downsample_sets = new std::unordered_set<int64_t>[max_segment];
    for (int64_t iv = 0; iv < max_segment; ++iv)
        downsample_sets[iv] = std::unordered_set<int64_t>();

    int64_t index = 0;
    for (int64_t iz = 0; iz < input_grid_size[IB_Z]; ++iz) {
        for (int64_t iy = 0; iy < input_grid_size[IB_Y]; ++iy) {
            for (int64_t ix = 0; ix < input_grid_size[IB_X]; ++ix, ++index) {
                int64_t segment = segmentation[index];
                if (!segment) continue;

                int64_t iw = (int64_t) (iz / zdown);
                int64_t iv = (int64_t) (iy / ydown);
                int64_t iu = (int64_t) (ix / xdown);

                int64_t downsample_index = iw * output_sheet_size + iv * output_row_size + iu;
                downsample_sets[segment].insert(downsample_index);
            }
        }
    }

    // write the downsampling information
    char downsample_filename[4096];
    sprintf(downsample_filename, "skeletons/%s/downsample-%03ldx%03ldx%03ld.bytes", prefix, output_resolution[IB_X], output_resolution[IB_Y], output_resolution[IB_Z]);

    // open the output file
    FILE *dfp = fopen(downsample_filename, "wb");
    if (!dfp) { fprintf(stderr, "Failed to write to %s\n", downsample_filename); exit(-1); }

    // write the upsampling information
    char upsample_filename[4096];
    sprintf(upsample_filename, "skeletons/%s/upsample-%03ldx%03ldx%03ld.bytes", prefix, output_resolution[IB_X], output_resolution[IB_Y], output_resolution[IB_Z]);

    // open the output file
    FILE *ufp = fopen(upsample_filename, "wb");
    if (!ufp) { fprintf(stderr, "Failed to write to %s\n", upsample_filename); exit(-1); }

    // write the number of segments
    fwrite(&output_grid_size[IB_Z], sizeof(int64_t), 1, dfp);
    fwrite(&output_grid_size[IB_Y], sizeof(int64_t), 1, dfp);
    fwrite(&output_grid_size[IB_X], sizeof(int64_t), 1, dfp);
    fwrite(&max_segment, sizeof(int64_t), 1, dfp);

    // write the output file size of the upsample version
    fwrite(&(input_grid_size[IB_Z]), sizeof(int64_t), 1, ufp);
    fwrite(&(input_grid_size[IB_Y]), sizeof(int64_t), 1, ufp);
    fwrite(&(input_grid_size[IB_X]), sizeof(int64_t), 1, ufp);
    fwrite(&max_segment, sizeof(int64_t), 1, ufp);

    // output values for downsampling
    for (int64_t label = 0; label < max_segment; ++label) {
        // write the size for this set
        int64_t nelements = downsample_sets[label].size();
        fwrite(&nelements, sizeof(int64_t), 1, dfp);
        fwrite(&nelements, sizeof(int64_t), 1, ufp);
        for (std::unordered_set<int64_t>::iterator it = downsample_sets[label].begin(); it != downsample_sets[label].end(); ++it) {
            int64_t element = *it;
            fwrite(&element, sizeof(int64_t), 1, dfp);

            int64_t iz = element / (output_grid_size[IB_Y] * output_grid_size[IB_X]);
            int64_t iy = (element - iz * output_grid_size[IB_Y] * output_grid_size[IB_X]) / output_grid_size[IB_X];
            int64_t ix = element % output_grid_size[IB_X];

            int64_t zmin = (int64_t) (zdown * iz);
            int64_t ymin = (int64_t) (ydown * iy);
            int64_t xmin = (int64_t) (xdown * ix);

            int64_t zmax = std::min((int64_t) ceil(zdown * (iz + 1) + 1), input_grid_size[IB_Z]);
            int64_t ymax = std::min((int64_t) ceil(ydown * (iy + 1) + 1), input_grid_size[IB_Y]);
            int64_t xmax = std::min((int64_t) ceil(xdown * (ix + 1) + 1), input_grid_size[IB_X]);

            double closest_to_center = input_grid_size[IB_Z] * input_grid_size[IB_Y] * input_grid_size[IB_X];
            int64_t upsample_index = -1;

            int64_t zcenter = (zmax + zmin) / 2;
            int64_t ycenter = (ymax + ymin) / 2;
            int64_t xcenter = (xmax + xmin) / 2;

            for (int64_t iw = zmin; iw < zmax; ++iw) {
                for (int64_t iv = ymin; iv < ymax; ++iv) {
                    for (int64_t iu = xmin; iu < xmax; ++iu) {
                        int64_t linear_index = iw * input_grid_size[IB_Y] * input_grid_size[IB_X] + iv * input_grid_size[IB_X] + iu;


                        // find the closest point to the center
                        if (segmentation[linear_index] != label) continue;

                        double distance = abs(iw - zcenter) + abs(iv - ycenter) + abs(iu - xcenter);
                        if (distance < closest_to_center) {
                            closest_to_center = distance;
                            upsample_index = linear_index;
                        }
                    }
                }
            }

            fwrite(&upsample_index, sizeof(int64_t), 1, ufp);
        }
    }

    // close the file
    fclose(dfp);
    fclose(ufp);

    // free memory
    delete[] downsample_sets;
}
