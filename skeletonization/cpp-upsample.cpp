/* c++ file to upsample the skeletons to full resolution */

#include <math.h>
#include <unordered_set>
#include <map>
#include <set>
#include <vector>
#include "cpp-generate_skeletons.h"



// global variables for upsampling operation

static std::map<int64_t, int64_t> *down_to_up;
static int64_t *segmentation;
static unsigned char *skeleton;
static std::set<std::pair<int64_t, int64_t> > connected_joints;


// convenient variables for moving between high and low resolutions

static float zdown;
static float ydown;
static float xdown;

static int64_t up_grid_size[3];
static int64_t up_nentries;
static int64_t up_sheet_size;
static int64_t up_row_size;

static int64_t down_grid_size[3];
static int64_t down_nentries;
static int64_t down_sheet_size;
static int64_t down_row_size;



// conver the index to indices
static void IndexToIndices(int64_t iv, int64_t &ix, int64_t &iy, int64_t &iz)
{
    iz = iv / down_sheet_size;
    iy = (iv - iz * down_sheet_size) / down_row_size;
    ix = iv % down_row_size;
}




static int MapDown2Up(const char *prefix, int64_t skeleton_resolution[3])
{
    // get the downsample filename
    char downsample_filename[4096];
    sprintf(downsample_filename, "skeletons/%s/downsample-%03ldx%03ldx%03ld.bytes", prefix, skeleton_resolution[IB_X], skeleton_resolution[IB_Y], skeleton_resolution[IB_Z]);

    FILE *dfp = fopen(downsample_filename, "rb");
    if (!dfp) { fprintf(stderr, "Failed to read %s\n", downsample_filename); return 0; }

    // get the upsample filename
    char upsample_filename[4096];
    sprintf(upsample_filename, "skeletons/%s/upsample-%03ldx%03ldx%03ld.bytes", prefix, skeleton_resolution[IB_X], skeleton_resolution[IB_Y], skeleton_resolution[IB_Z]);

    FILE *ufp = fopen(upsample_filename, "rb");
    if (!ufp) { fprintf(stderr, "Failed to read %s\n", upsample_filename); return 0; }

    // read downsample header
    int64_t down_max_segment;
    if (fread(&(down_grid_size[IB_Z]), sizeof(int64_t), 1, dfp) != 1) { fprintf(stderr, "Failed to read %s\n", downsample_filename); return 0; }
    if (fread(&(down_grid_size[IB_Y]), sizeof(int64_t), 1, dfp) != 1) { fprintf(stderr, "Failed to read %s\n", downsample_filename); return 0; }
    if (fread(&(down_grid_size[IB_X]), sizeof(int64_t), 1, dfp) != 1) { fprintf(stderr, "Failed to read %s\n", downsample_filename); return 0; }
    if (fread(&down_max_segment, sizeof(int64_t), 1, dfp) != 1) { fprintf(stderr, "Failed to read %s\n", downsample_filename); return 0; }

    // read upsample header
    int64_t up_max_segment;
    if (fread(&(up_grid_size[IB_Z]), sizeof(int64_t), 1, ufp) != 1) { fprintf(stderr, "Failed to read %s\n", upsample_filename); return 0; }
    if (fread(&(up_grid_size[IB_Y]), sizeof(int64_t), 1, ufp) != 1) { fprintf(stderr, "Failed to read %s\n", upsample_filename); return 0; }
    if (fread(&(up_grid_size[IB_X]), sizeof(int64_t), 1, ufp) != 1) { fprintf(stderr, "Failed to read %s\n", upsample_filename); return 0; }
    if (fread(&up_max_segment, sizeof(int64_t), 1, ufp) != 1) { fprintf(stderr, "Failed to read %s\n", upsample_filename); return 0; }

    down_to_up = new std::map<int64_t, int64_t>[up_max_segment];
    for (int64_t label = 0; label < up_max_segment; ++label) {
        down_to_up[label] = std::map<int64_t, int64_t>();

        int64_t down_nelements, up_nelements;
        if (fread(&down_nelements, sizeof(int64_t), 1, dfp) != 1) { fprintf(stderr, "Failed to read %s\n", downsample_filename); return 0; }
        if (fread(&up_nelements, sizeof(int64_t), 1, ufp) != 1) { fprintf(stderr, "Failed to read %s\n", upsample_filename); return 0; }

        int64_t *down_elements = new int64_t[down_nelements];
        int64_t *up_elements = new int64_t[up_nelements];
        if (fread(down_elements, sizeof(int64_t), down_nelements, dfp) != (uint64_t)down_nelements) { fprintf(stderr, "Failed to read %s\n", downsample_filename); return 0; }
        if (fread(up_elements, sizeof(int64_t), up_nelements, ufp) != (uint64_t)up_nelements) { fprintf(stderr, "Failed to read %s\n", upsample_filename); return 0; }

        for (int64_t ie = 0; ie < down_nelements; ++ie)
            down_to_up[label][down_elements[ie]] = up_elements[ie];
    }

    fclose(dfp);
    fclose(ufp);

    return 1;
}



static void FindEndpointVector(int64_t index, double &vx, double &vy, double &vz)
{
    std::vector<int64_t> path_from_endpoint = std::vector<int64_t>();
    path_from_endpoint.push_back(index);

    while (path_from_endpoint.size() < 4) {
        short nneighbors = 0;
        int64_t only_neighbor = -1;

        int64_t ix, iy, iz;
        IndexToIndices(index, ix, iy, iz);

        for (int64_t iw = iz - 1; iw <= iz + 1; ++iw) {
            if (iw < 0 || iw >= down_grid_size[IB_Z]) continue;
            for (int64_t iv = iy - 1; iv <= iy + 1; ++iv) {
                if (iv < 0 || iv >= down_grid_size[IB_Y]) continue;
                for (int64_t iu = ix - 1; iu <= ix + 1; ++iu) {
                    if (iu < 0 || iu >= down_grid_size[IB_X]) continue;

                    int64_t neighbor_index = iw * down_grid_size[IB_Y] * down_grid_size[IB_X] + iv * down_grid_size[IB_X] + iu;
                    if (!skeleton[neighbor_index]) continue;

                    // skip if the neighbor is this index (i.e., it is not a neighbor)
                    if (neighbor_index == index) continue;

                    nneighbors += 1;
                    only_neighbor = neighbor_index;
                }
            }
        }

        // if there were no neighbors break since there are no more endpoints
        if (!nneighbors) break;
        // if there are two neighbors break since there is a split
        else if (nneighbors > 1) break;
        else {
            // mask out this skeleton point so next iteration works
            skeleton[index] = 0;
            // reset the index to the neighbors value
            index = only_neighbor;
            // add this neighbor to the path
            path_from_endpoint.push_back(index);
        }
    }

    // reset the skeletons
    for (uint64_t iv = 0; iv < path_from_endpoint.size(); ++iv) {
        skeleton[path_from_endpoint[iv]] = 1;
    }

    // find the vector
    if (path_from_endpoint.size() == 1) {
        vx = 0.0;
        vy = 0.0;
        vz = 0.0;

        // cannot normalize
        return;
    }
    else {
        int64_t ix, iy, iz, ii, ij, ik;
        IndexToIndices(path_from_endpoint[0], ix, iy, iz);
        IndexToIndices(path_from_endpoint[path_from_endpoint.size() - 1], ii, ij, ik);

        vx = ix - ii;
        vy = iy - ij;
        vz = iz - ik;
    }

    // we do not change the coordinate system for the vectors from anisotropic to isotropic
    // this allows us to easily compute edges because we can multiply by the resolutions to convert
    // anisotropic coordinates to isotropic ones

    double normalization = sqrt(vx * vx + vy * vy + vz * vz);
    vx = vx / normalization;
    vy = vy / normalization;
    vz = vz / normalization;
}



void CppFindEndpointVectors(const char *prefix, int64_t skeleton_resolution[3], float output_resolution[3])
{
    // get the mapping from downsampled locations to upsampled ones
    if (!MapDown2Up(prefix, skeleton_resolution)) return;

    // get downsample ratios
    zdown = ((float) skeleton_resolution[IB_Z]) / output_resolution[IB_Z];
    ydown = ((float) skeleton_resolution[IB_Y]) / output_resolution[IB_Y];
    xdown = ((float) skeleton_resolution[IB_X]) / output_resolution[IB_X];

    // set global variables
    up_nentries = up_grid_size[IB_Z] * up_grid_size[IB_Y] * up_grid_size[IB_X];
    up_sheet_size = up_grid_size[IB_Y] * up_grid_size[IB_X];
    up_row_size = up_grid_size[IB_X];

    down_nentries = down_grid_size[IB_Z] * down_grid_size[IB_Y] * down_grid_size[IB_X];
    down_sheet_size = down_grid_size[IB_Y] * down_grid_size[IB_X];
    down_row_size = down_grid_size[IB_X];

    // I/O filenames
    char input_filename[4096];
    sprintf(input_filename, "skeletons/%s/thinning-%03ldx%03ldx%03ld-downsample-skeleton.pts", prefix, skeleton_resolution[IB_X], skeleton_resolution[IB_Y], skeleton_resolution[IB_Z]);

    char output_filename[4096];
    sprintf(output_filename, "skeletons/%s/thinning-%03ldx%03ldx%03ld-endpoint-vectors.vec", prefix, skeleton_resolution[IB_X], skeleton_resolution[IB_Y], skeleton_resolution[IB_Z]);

    // open files for read/write
    FILE *rfp = fopen(input_filename, "rb");
    if (!rfp) { fprintf(stderr, "Failed to read %s\n", input_filename); return; }

    FILE *wfp = fopen(output_filename, "wb");
    if (!wfp) { fprintf(stderr, "Failed to write %s\n", output_filename); return; }

    // read header
    int64_t max_label;
    int64_t input_grid_size[3];
    if (fread(&(input_grid_size[IB_Z]), sizeof(int64_t), 1, rfp) != 1) { fprintf(stderr, "Failed to read %s\n", input_filename); return; }
    if (fread(&(input_grid_size[IB_Y]), sizeof(int64_t), 1, rfp) != 1) { fprintf(stderr, "Failed to read %s\n", input_filename); return; }
    if (fread(&(input_grid_size[IB_X]), sizeof(int64_t), 1, rfp) != 1) { fprintf(stderr, "Failed to read %s\n", input_filename); return; }
    if (fread(&max_label, sizeof(int64_t), 1, rfp) != 1) { fprintf(stderr, "Failed to read %s\n", input_filename); return; }

    // write the header
    if (fwrite(&(up_grid_size[IB_Z]), sizeof(int64_t), 1, wfp) != 1) { fprintf(stderr, "Failed to write %s\n", output_filename); return; }
    if (fwrite(&(up_grid_size[IB_Y]), sizeof(int64_t), 1, wfp) != 1) { fprintf(stderr, "Failed to write %s\n", output_filename); return; }
    if (fwrite(&(up_grid_size[IB_X]), sizeof(int64_t), 1, wfp) != 1) { fprintf(stderr, "Failed to write %s\n", output_filename); return; }
    if (fwrite(&max_label, sizeof(int64_t), 1, wfp) != 1) { fprintf(stderr, "Failed to write %s\n", output_filename); return; }

    for (int64_t label = 0; label < max_label; ++label) {
        int64_t nelements;
        if (fread(&nelements, sizeof(int64_t), 1, rfp) != 1) { fprintf(stderr, "Failed to read %s\n", input_filename); return; }

        skeleton = new unsigned char[down_nentries];
        for (int64_t iv = 0; iv < down_nentries; ++iv) skeleton[iv] = 0;

        // find all of the downsampled elements
        int64_t *down_elements = new int64_t[nelements];
        if (fread(down_elements, sizeof(int64_t), nelements, rfp) != (uint64_t)nelements) { fprintf(stderr, "Failed to read %s\n", input_filename); return; }

        int64_t nendpoints = 0;
        for (int64_t ie = 0; ie < nelements; ++ie) {
            if (down_elements[ie] < 0) {
                skeleton[-1 * down_elements[ie]] = 1;
                nendpoints++;
            }
            else skeleton[down_elements[ie]] = 1;
        }
        if (fwrite(&nendpoints, sizeof(int64_t), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); return; }

        // go through all down elements to find endpoints
        for (int64_t ie = 0; ie < nelements; ++ie) {
            if (down_elements[ie] >= 0) continue;

            double vx, vy, vz;
            FindEndpointVector(-1 * down_elements[ie], vx, vy, vz);

            // get the corresponding up element for this endpoint
            int64_t up_element = down_to_up[label][-1 * down_elements[ie]];

            // save the up element with the vector
            if (fwrite(&up_element, sizeof(int64_t), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); return; }
            if (fwrite(&vz, sizeof(double), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); return; }
            if (fwrite(&vy, sizeof(double), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); return; }
            if (fwrite(&vx, sizeof(double), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); return; }
        }


        delete[] skeleton;
    }

    // close the file
    fclose(rfp);
    fclose(wfp);

    delete[] down_to_up;
}



// operation that takes skeletons and
void CppApplyUpsampleOperation(const char *prefix, int64_t *input_segmentation, int64_t skeleton_resolution[3], float output_resolution[3])
{
    // get the mapping from downsampled locations to upsampled ones
    if (!MapDown2Up(prefix, skeleton_resolution)) return;

    // get a list of labels for each downsampled index
    segmentation = input_segmentation;

    // get downsample ratios
    zdown = ((float) skeleton_resolution[IB_Z]) / output_resolution[IB_Z];
    ydown = ((float) skeleton_resolution[IB_Y]) / output_resolution[IB_Y];
    xdown = ((float) skeleton_resolution[IB_X]) / output_resolution[IB_X];

    // set global variables
    up_nentries = up_grid_size[IB_Z] * up_grid_size[IB_Y] * up_grid_size[IB_X];
    up_sheet_size = up_grid_size[IB_Y] * up_grid_size[IB_X];
    up_row_size = up_grid_size[IB_X];

    down_nentries = down_grid_size[IB_Z] * down_grid_size[IB_Y] * down_grid_size[IB_X];
    down_sheet_size = down_grid_size[IB_Y] * down_grid_size[IB_X];
    down_row_size = down_grid_size[IB_X];

    // I/O filenames
    char input_filename[4096];
    sprintf(input_filename, "skeletons/%s/thinning-%03ldx%03ldx%03ld-downsample-skeleton.pts", prefix, skeleton_resolution[IB_X], skeleton_resolution[IB_Y], skeleton_resolution[IB_Z]);


    char output_filename[4096];
    sprintf(output_filename, "skeletons/%s/thinning-%03ldx%03ldx%03ld-upsample-skeleton.pts", prefix, skeleton_resolution[IB_X], skeleton_resolution[IB_Y], skeleton_resolution[IB_Z]);

    // open files for read/write
    FILE *rfp = fopen(input_filename, "rb");
    if (!rfp) { fprintf(stderr, "Failed to read %s\n", input_filename); return; }

    FILE *wfp = fopen(output_filename, "wb");
    if (!wfp) { fprintf(stderr, "Failed to write %s\n", output_filename); return; }

    // read header
    int64_t max_label;
    int64_t input_grid_size[3];
    if (fread(&(input_grid_size[IB_Z]), sizeof(int64_t), 1, rfp) != 1) { fprintf(stderr, "Failed to read %s\n", input_filename); return; }
    if (fread(&(input_grid_size[IB_Y]), sizeof(int64_t), 1, rfp) != 1) { fprintf(stderr, "Failed to read %s\n", input_filename); return; }
    if (fread(&(input_grid_size[IB_X]), sizeof(int64_t), 1, rfp) != 1) { fprintf(stderr, "Failed to read %s\n", input_filename); return; }
    if (fread(&max_label, sizeof(int64_t), 1, rfp) != 1) { fprintf(stderr, "Failed to read %s\n", input_filename); return; }

    // write the header
    if (fwrite(&(up_grid_size[IB_Z]), sizeof(int64_t), 1, wfp) != 1) { fprintf(stderr, "Failed to write %s\n", output_filename); return; }
    if (fwrite(&(up_grid_size[IB_Y]), sizeof(int64_t), 1, wfp) != 1) { fprintf(stderr, "Failed to write %s\n", output_filename); return; }
    if (fwrite(&(up_grid_size[IB_X]), sizeof(int64_t), 1, wfp) != 1) { fprintf(stderr, "Failed to write %s\n", output_filename); return; }
    if (fwrite(&max_label, sizeof(int64_t), 1, wfp) != 1) { fprintf(stderr, "Failed to write %s\n", output_filename); return; }

    // go through all skeletons
    for (int64_t label = 0; label < max_label; ++label) {
        int64_t nelements;
        if (fread(&nelements, sizeof(int64_t), 1, rfp) != 1) { fprintf(stderr, "Failed to read %s\n", input_filename); return; }
        if (fwrite(&nelements, sizeof(int64_t), 1, wfp) != 1) { fprintf(stderr, "Failed to write %s\n", output_filename); return; }

        // just run naive method where endpoints in downsampled are transfered
        int64_t *down_elements = new int64_t[nelements];
        if (fread(down_elements, sizeof(int64_t), nelements, rfp) != (uint64_t)nelements) { fprintf(stderr, "Failed to read %s\n", input_filename); return; }

        int64_t *up_elements = new int64_t[nelements];
        for (int64_t ie = 0; ie < nelements; ++ie) {
            int64_t down_index = down_elements[ie];

            if (down_index < 0) {
                down_index = -1 * down_index;
                up_elements[ie] = -1 * down_to_up[label][down_index];
            }
            else {
                up_elements[ie] = down_to_up[label][down_index];
            }
        }

        if (fwrite(up_elements, sizeof(int64_t), nelements, wfp) != (uint64_t)nelements) { fprintf(stderr, "Failed to write %s\n", output_filename); return; }

        // free memory
        delete[] down_elements;
        delete[] up_elements;
    }

    // free memory
    delete[] down_to_up;

    // close the files
    fclose(rfp);
    fclose(wfp);
}
