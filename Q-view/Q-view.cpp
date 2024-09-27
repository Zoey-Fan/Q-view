/*
 * date: 2024-04-26
 * author: Mengyu Ma@National University of Defense Technology
 *         Yifan Zhang@National University of Defense Technology
 * e-mail: mamengyu10@nudt.edu.cn
 *         zhangyifan.000@nudt.edu.cn
 * description:  3d viewshed query.
 * run:  ./Q-VIEW --index ../mergeindex/indexes/test_area.3idx
 * input:radius 250; cordinate for sensor 114.1670 22.2806; height for Sensor 100
 */
#include <iostream>
#include <fstream>
#include <sstream>
#include <omp.h>
#include <sys/time.h>
#include <png.h>
#include <ogr_spatialref.h>
#include "query3dRtree.h"
#include "kernel.h"
#include <proj.h>

using namespace std;

typedef tuple<double, double, double, double, double, double, double, double, double> ValueType;
typedef RTree<ValueType, double, 3> RTree3d;

struct Rect3d
{
    Rect3d() {}

    Rect3d(double minX, double minY, double minZ, double maxX, double maxY, double maxZ)
    {
        min[0] = minX;
        min[1] = minY;
        min[2] = minZ;
        max[0] = maxX;
        max[1] = maxY;
        max[2] = maxZ;
    }
    double min[3];
    double max[3];
};
// various operations
class Vec3d
{
public:
    double x, y, z;
    // Calculate the dot product of this vector with another vector b.
    double dot(const Vec3d &b)
    {
        return Vec3d::x * b.x + Vec3d::y * b.y + Vec3d::z * b.z;
    }
    // Calculate the cross product of this vector with another vector b.
    Vec3d cross(const Vec3d &b)
    {
        return Vec3d(
            Vec3d::y * b.z - Vec3d::z * b.y,
            Vec3d::z * b.x - Vec3d::x * b.z,
            Vec3d::x * b.y - Vec3d::y * b.x);
    }
    // Normalize this vector to make its length equal to 1.
    Vec3d normalize()
    {
        const double s = 1.0f / sqrtf(Vec3d::x * Vec3d::x + Vec3d::y * Vec3d::y + Vec3d::z * Vec3d::z);
        return Vec3d(Vec3d::x * s, Vec3d::y * s, Vec3d::z * s);
    }
    Vec3d operator+(const Vec3d &b)
    {
        return Vec3d(
            Vec3d::x + b.x,
            Vec3d::y + b.y,
            Vec3d::z + b.z);
    }
    Vec3d operator+=(const Vec3d &b)
    {
        *this = Vec3d::operator+(b);
        return *this;
    }

    Vec3d operator-(const Vec3d &b)
    {
        return Vec3d(
            Vec3d::x - b.x,
            Vec3d::y - b.y,
            Vec3d::z - b.z);
    }
    Vec3d operator-=(const Vec3d &b)
    {
        *this = Vec3d::operator-(b);
        return *this;
    }
    Vec3d operator*(const Vec3d &b)
    {
        return Vec3d(
            Vec3d::x * b.x,
            Vec3d::y * b.y,
            Vec3d::z * b.z);
    }
    Vec3d operator*=(const Vec3d &b)
    {
        *this = Vec3d::operator*(b);
        return *this;
    }
    Vec3d operator*(double b)
    {
        return Vec3d(
            Vec3d::x * b,
            Vec3d::y * b,
            Vec3d::z * b);
    }
    Vec3d operator*=(double b)
    {
        *this = Vec3d::operator*(b);
        return *this;
    }
    Vec3d operator/(const Vec3d &b)
    {
        return Vec3d(
            Vec3d::x / b.x,
            Vec3d::y / b.y,
            Vec3d::z / b.z);
    }
    Vec3d operator/=(const Vec3d &b)
    {
        *this = Vec3d::operator/(b);
        return *this;
    }
    Vec3d operator/(double b)
    {
        return Vec3d(
            Vec3d::x * b,
            Vec3d::y * b,
            Vec3d::z * b);
    }
    Vec3d operator/=(double b)
    {
        *this = Vec3d::operator/(b);
        return *this;
    }
    Vec3d(double x, double y, double z)
    {
        Vec3d::x = x;
        Vec3d::y = y;
        Vec3d::z = z;
    }
    Vec3d(double x)
    {
        Vec3d::x = x;
        Vec3d::y = x;
        Vec3d::z = x;
    }
    Vec3d()
    {
        //
    }
    ~Vec3d()
    {
        //
    }
};

#define EPSILON 0.000001f
// Check if a line segment intersects with a triangle.
bool lineSegIntersectTri(Vec3d line[2], Vec3d tri[3], Vec3d *point)
{
    // Calculate the two edges of the triangle and the direction vector of the line segment.
    Vec3d e0 = tri[1] - tri[0];
    Vec3d e1 = tri[2] - tri[0];
    Vec3d dir = line[1] - line[0];
    Vec3d dir_norm = dir.normalize();
    Vec3d h = dir_norm.cross(e1);
    const double a = e0.dot(h);

    if (a > -EPSILON && a < EPSILON)
    {
        return false;
    }

    Vec3d s = line[0] - tri[0];
    const double f = 1.0f / a;
    const double u = f * s.dot(h);

    if (u < 0.0f || u > 1.0f)
    {
        return false;
    }

    Vec3d q = s.cross(e0);
    const double v = f * dir_norm.dot(q);

    if (v < 0.0f || u + v > 1.0f)
    {
        return false;
    }

    const double t = f * e1.dot(q);

    if (t > EPSILON && t < sqrtf(dir.dot(dir))) // segment intersection
    {
        if (point)
        {
            *point = line[0] + dir_norm * t;
        }

        return true;
    }

    return false;
}
// Find the maximum and minimum values in an array.
void find(double *arr, int nSize, double &aMax, double &aMin)
{
    aMax = *arr;
    aMin = *arr;
    for (int i = 1; i < nSize; i++)
    {
        if (*(arr + i) > aMax)
            aMax = *(arr + i);
        if (*(arr + i) < aMin)
            aMin = *(arr + i);
    }
}

#define MAXH 100000
#define MINH -100000

const double PI = 3.14159265358979323846;
// Degrees to radians function.
double degToRad(double degree){
    return degree * PI / 180.0;
}

// This function is a callback used to determine the height at a given (x, y) location.
double HeightCallback(ValueType id, const double x, const double y){
    // Create a triangle from the given id, which contains vertex positions.
    Vec3d tri[3] =
        {
            {get<0>(id), get<1>(id), get<2>(id)},
            {get<3>(id), get<4>(id), get<5>(id)},
            {get<6>(id), get<7>(id), get<8>(id)},
        };
    // Define a line segment from (x, y, MAXH) to (x, y, MINH).
    Vec3d line[2] =
        {
            {x, y, MAXH},
            {x, y, MINH},
        };
    Vec3d *point = new Vec3d();
    // Check if the line segment intersects with the triangle.
    if (lineSegIntersectTri(line, tri, point))
    {
        return point->z;
    }
    // If there is no intersection, return a value indicating to continue searching.
    return MINH - 1;
}

// This function is a callback used to check for intersection between a line segment and a triangle.
bool IntersectCallback(ValueType id, const double a[3], const double b[3])
{
    Vec3d tri[3] =
        {
            {get<0>(id), get<1>(id), get<2>(id)},
            {get<3>(id), get<4>(id), get<5>(id)},
            {get<6>(id), get<7>(id), get<8>(id)},
        };
    Vec3d line[2] =
        {
            {a[0], a[1], a[2]},
            {b[0], b[1], b[2]},
        };
    // If an intersection is found, return false to indicate the intersection.
    if (lineSegIntersectTri(line, tri, NULL))
    {
        return false;
    }
    // If no intersection is found, return true to indicate continue searching.
    return true;
}

void Usage()
{
    printf("Usage:   [--index:     input file path  ]\n");
}

// Write viewshed results into csv.
void writeViewshedToCsv(vs_t **viewshed, int ncols, int nrows, char *fileName)
{
    std::ofstream myfile;
    myfile.open(fileName);
    int count_0 = 0;
    int count_1 = 0;
    for (int i = 0; i < nrows; i++)
    {
        for (int j = 0; j < ncols; j++)
        {
            // int index = i * ncols + j;
            if (viewshed[i][j] == 1)
            {
                myfile << "1,";
                count_1 += 1;
            }
            else
            {
                myfile << "0,";
                count_0 += 1;
            }
        }
        myfile << "\n";
    }
    cout << "the invisible area is " << count_0 << " square meter" << endl;
    cout << "the visibile area is " << count_1 << " square meter" << endl;
    myfile.close();
}

int main(int nArgc, char **papszArgv)
{
    // Load the 3D Rtree index
    const char *index = NULL;
    for (int iArg = 1; iArg < nArgc; iArg++)
    {
        if (strcmp(papszArgv[iArg], "--index")==0)
        {
            index = papszArgv[iArg + 1];
        }
    }
    if (index == NULL)
    {
        cout << "[ERROR] Wrong input!" << endl;
        Usage();
        return 0;
    }

    struct timeval t1, t2, t3, t4;
    gettimeofday(&t1, NULL);
    RTree3d tree;
    tree.Load(index);
    gettimeofday(&t2, NULL);
    cout << "load" << (t2.tv_sec - t1.tv_sec) + (double)(t2.tv_usec - t1.tv_usec) / 1000000.0 << endl;

    vector<ValueType> itemList;
    // Define the size of image.
    int tile_size = 500;
    int min_x = 835000;
    int min_y = 815500;

    // Input the raidus
    int radius = 0;
    cout << "\nEnter Sensor Radius(in Pixels): ";
    cin >> radius;
    // enter the sensor cordinate
    double sensor_x = 0; // longtitude
    double sensor_y = 0; // latitude
    cout << "Enter cordinate for sensor: ";
    cin >> sensor_x >> sensor_y;

    PJ_CONTEXT *context = proj_context_create();
    PJ *proj_from = proj_create_crs_to_crs(context, "EPSG:4326", "EPSG:2326", NULL);
    PJ_COORD coord = proj_coord(sensor_y, sensor_x, 0, 0);
    PJ_COORD result = proj_trans(proj_from, PJ_FWD, coord);
    // delete
    proj_destroy(proj_from);
    proj_context_destroy(context);

    sensor_x = result.xy.y-min_x;
    sensor_y = result.xy.x-min_y;
    //std::cout<<"sensor_x: "<<sensor_x<<" and sensor_y: "<<sensor_y<<std::endl;

    // spatial index leftdown cordinate.
    double indexRange_x = -1975;
    double indexRange_y = 43;

    // Calculate area of interest.
    int minX = max(indexRange_x, indexRange_x + sensor_x - radius);
    int minY = max(indexRange_y, indexRange_y + sensor_y - radius);
    int maxX = min(indexRange_x + tile_size - 1, indexRange_x + sensor_x + radius);
    int maxY = min(indexRange_y + tile_size - 1, indexRange_y + sensor_y + radius);

    //cout << "minX: " << minX << " minY: " << minY << " maxX: " << maxX << " maxY: " << maxY << endl;
    int width = ((maxX - minX) + 1); //*2 high resolution
    int height = ((maxY - minY) + 1);
    //cout << "width*height= " << width << " * " << height << endl;
    // Calculate the actural height.
    int terrainHeight = tree.Getheight3d(indexRange_x + sensor_x, indexRange_y + sensor_y, HeightCallback) + 0.1;
    cout << "the actural height of oberversion is " << terrainHeight << endl;
    cout << "Enter Height for Sensor: ";
    int sensor_h = 0;
    cin >> sensor_h;
    sensor_h = terrainHeight + sensor_h;

    // Writting into png
    png_bytep *rowPointers = (png_bytep *)malloc(height * sizeof(png_bytep));
    for (int i = 0; i < height; i++)
        rowPointers[i] = (png_bytep)malloc(width * 4);
    png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    png_infop info_ptr = png_create_info_struct(png_ptr);
    char fullPath[256];
    sprintf(fullPath, "viewshed_%d_%d_%d_%d.png", static_cast<int>(sensor_x), static_cast<int>(sensor_y), sensor_h, radius);
    FILE *temp_png = fopen(fullPath, "wb");
    png_init_io(png_ptr, temp_png);
    png_set_IHDR(png_ptr, info_ptr, width, height, 8, PNG_COLOR_TYPE_RGB_ALPHA, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
    png_write_info(png_ptr, info_ptr);

    // Create viewshed array
    vs_t **viewshed = new vs_t *[height];
    for (int i = 0; i < height; i++)
        viewshed[i] = new vs_t[width];

    gettimeofday(&t3, NULL);
    // viewshed query
    #pragma omp parallel for num_threads(8) schedule(dynamic)
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            int ii = height - 1 - i;
            double z = tree.Getheight3d(j + minX, i + minY, HeightCallback) + 0.1;
            Rect3d search_rect(j + minX, i + minY, z, minX + radius, minY + radius, sensor_h); // 250*250 500*500

            if (z < 0)
                z = 0;
            if (tree.Intersect3d(search_rect.min, search_rect.max, IntersectCallback)) // true
            {
                // invisible
                rowPointers[ii][4 * j] = 255;   // red
                rowPointers[ii][4 * j + 1] = 0; // green
                rowPointers[ii][4 * j + 2] = 0; // blue
                rowPointers[ii][4 * j + 3] = 255;
                viewshed[i][j] = 0;
            }
            else
            {
                // visible
                rowPointers[ii][4 * j] = 0;       // red
                rowPointers[ii][4 * j + 1] = 255; // green
                rowPointers[ii][4 * j + 2] = 0;   // blue
                rowPointers[ii][4 * j + 3] = 255;
                viewshed[i][j] = 1;
            }
            if ((minX + radius == j) && (minY + radius == i))
            {
                rowPointers[ii][4 * j] = 0;     // red
                rowPointers[ii][4 * j + 1] = 0; // green
                rowPointers[ii][4 * j + 2] = 0; // blue
                rowPointers[ii][4 * j + 3] = 255;
            }
        }
    }
    gettimeofday(&t4, NULL);
    printf("%lfs\n", (t4.tv_sec - t3.tv_sec) + (double)(t4.tv_usec - t3.tv_usec) / 1000000.0);

    // writing into csv
    char buffer[256];
    sprintf(buffer, "viewshed_%d_%d_%d_%d.csv", static_cast<int>(sensor_x), static_cast<int>(sensor_y), sensor_h, radius);
    writeViewshedToCsv(viewshed, width, height, buffer);
    // writing into png
    png_write_image(png_ptr, rowPointers);
    png_write_end(png_ptr, NULL);
    fclose(temp_png);
    // delete
    for (int i = 0; i < height; i++)
    {
        delete[] viewshed[i];
    }
    delete[] viewshed;

    return 0;
}
