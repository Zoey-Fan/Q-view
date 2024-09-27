/*
 * date: 2024-04-26
 * author: Mengyu Ma@National University of Defense Technology
 *         Yifan Zhang@National University of Defense Technology
 * e-mail: mamengyu10@nudt.edu.cn
 *         zhangyifan.000@nudt.edu.cn
 * description:  3d curve query.
 * run: ./Q-CURVE --index ../mergeindex/indexes/test_area.3idx --positionStart 114.166 22.2825 1 --positionEnd 114.166 22.280 1 --velocity 50 --angle 50
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <omp.h>
#include <sys/time.h>
#include <math.h>
#include <iomanip>
#include "query3dRtree.h"
#include "kernel.h"
#include <cstdlib>
#include <proj.h>
using namespace std;

typedef tuple<double, double, double, double, double, double, double, double, double> ValueType;
typedef RTree<ValueType, double, 3> RTree3d;

int TILE_SIZE = 500;
constexpr double DEG_TO_RAD_LOCAL = 3.1415926535897932 / 180.0;
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

class Vec3d
{
public:
    double x, y, z;

    double dot(const Vec3d &b)
    {
        return Vec3d::x * b.x + Vec3d::y * b.y + Vec3d::z * b.z;
    }

    Vec3d cross(const Vec3d &b)
    {
        return Vec3d(
            Vec3d::y * b.z - Vec3d::z * b.y,
            Vec3d::z * b.x - Vec3d::x * b.z,
            Vec3d::x * b.y - Vec3d::y * b.x);
    }

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

struct Point3D
{
    double start_x, start_y, start_z;
    double end_x, end_y, end_z;
    double positionXAtTHalf, positionYAtTHalf;
    int radius;
};

struct Velocity3D
{
    double vz, vxy, angle, v;
};

struct Object
{
    Point3D position;
    Velocity3D velocity;
};

#define EPSILON 0.000001f
// Check if a line segment intersects with a triangle.
bool lineSegIntersectTri(Vec3d line[2], Vec3d tri[3], Vec3d *point)
{
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

// This function is a callback used to determine the height at a given (x, y) location.
double HeightCallback(ValueType id, const double x, const double y)
{
    Vec3d tri[3] =
        {
            {get<0>(id), get<1>(id), get<2>(id)},
            {get<3>(id), get<4>(id), get<5>(id)},
            {get<6>(id), get<7>(id), get<8>(id)},
        };
    Vec3d line[2] =
        {
            {x, y, MAXH},
            {x, y, MINH},
        };
    Vec3d *point = new Vec3d();
    if (lineSegIntersectTri(line, tri, point))
    {
        return point->z;
    }
    return MINH - 1; // keep going
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
    if (lineSegIntersectTri(line, tri, NULL))
    {
        //~ printf("tri1:%lf %lf %lf \n", get<0>(id) , get<1>(id) , get<2>(id) );
        //~ printf("tri2:%lf %lf %lf \n",get<3>(id) , get<4>(id) , get<5>(id));
        //~ printf("tri3:%lf %lf %lf \n", get<6>(id) , get<7>(id) , get<8>(id) );
        //~ printf("a:%lf %lf %lf \n",a[0],a[1],a[2]);
        //~ printf("b:%lf %lf %lf \n",b[0],b[1],b[2]);
        return false;
    }
    return true; // keep going
}

void Usage()
{
    printf("Usage:           [--index:     path for 3d model index            ]\n"
           "                 [--positionStart:  start position                ]\n"
           "                 [--positionEnd:    end position                  ]\n"
           "                 [--velocity:       velocity                      ]\n"
           "                 [--angle:          vertical angle                ]\n"
          );
}

// Write viewshed results into csv.
void writeViewshedToCsv(vs_t **viewshed, int ncols, int nrows, char *fileName)
{
    std::ofstream myfile_viewshed;
    myfile_viewshed.open(fileName);

    for (int i = 0; i < nrows; i++)
    {
        for (int j = 0; j < ncols; j++)
        {
            // int index = i * ncols + j;
            if (viewshed[i][j] == 1)
            {
                myfile_viewshed << "1,";
            }
            else
            {
                myfile_viewshed << "0,";
            }
        }
        myfile_viewshed << "\n";
    }
    myfile_viewshed.close();
}

void solveQuadraticEquation(double x1, double y1, double x2, double y2, double x3, double y3, double *a, double *b, double *c){

    double denom = (x1 - x2) * (x1 - x3) * (x2 - x3);

    *a = (y1 / ((x1 - x2) * (x1 - x3))) + (y2 / ((x2 - x1) * (x2 - x3))) + (y3 / ((x3 - x1) * (x3 - x2)));

    *b = (-y1 * (x2 + x3) / ((x1 - x2) * (x1 - x3))) - (y2 * (x1 + x3) / ((x2 - x1) * (x2 - x3))) - (y3 * (x1 + x2) / ((x3 - x1) * (x3 - x2)));

    *c = (y1 * x2 * x3 / ((x1 - x2) * (x1 - x3))) + (y2 * x1 * x3 / ((x2 - x1) * (x2 - x3))) + (y3 * x1 * x2 / ((x3 - x1) * (x3 - x2)));
}

// Function to calculate the trajectory of an object and write it to a CSV file.
std::vector<std::vector<double>> calculateTrajectory(Object object, RTree3d tree)
{
    char filePath[256];
    sprintf(filePath, "trajectory_%d_%d_%d_allPoints.csv",
            static_cast<int>(object.position.start_x),
            static_cast<int>(object.position.start_y),
            static_cast<int>(object.position.start_z));

    // Calculate the direction vector.
    double directionX = object.position.end_x - object.position.start_x;
    double directionY = object.position.end_y - object.position.start_y;
    double distance = sqrt(directionX * directionX + directionY * directionY);
    // Calculate the time of flight.
    double timeOfFlight = distance / object.velocity.vxy;
    // Calculate the position at time t/2.
    object.position.positionYAtTHalf = (-9.81 * 0.5) * (0.5 * timeOfFlight) * (0.5 * timeOfFlight) + object.velocity.vz * (0.5 * timeOfFlight);
    object.position.positionXAtTHalf = object.velocity.vxy * (timeOfFlight * 0.5);

    double a, b, c;
    // Solve the quadratic equation to find the trajectory.
    solveQuadraticEquation(0.0, 0.0, distance, object.position.end_z - object.position.start_z, object.position.positionXAtTHalf, object.position.positionYAtTHalf, &a, &b, &c);
    // wheather the curve is suitable
    if (object.position.positionYAtTHalf < 0.5 * (object.position.end_z - object.position.start_z))
    {
        cout << "Cannot generate the correct projectile trajectory! Exiting program.";
        exit(EXIT_FAILURE);
    }

    // Discrete the trajectory
    int stepSize = 10;                     // step of distance is 10 meters
    int numPoints = static_cast<int>(ceil(distance / stepSize)) + 1; // we got n points

    // 3d_points
    vector<double> trajectoryZ(numPoints);
    vector<double> trajectoryX(numPoints);
    vector<double> trajectoryY(numPoints);

    double* discretizedPoints = new double[numPoints];
    for (int i = 0; i < numPoints; i++)
    {
        discretizedPoints[i] = (i * distance) / (numPoints - 1); 
        // cout << discretizedPoints[i] << ", ";
    }
    // cout << endl;
    double direction = atan2(object.position.end_y - object.position.start_y, object.position.end_x - object.position.start_x);
   
    for (int i = 0; i < numPoints; i++) {
        trajectoryZ[i] = a * discretizedPoints[i] * discretizedPoints[i] + b * discretizedPoints[i] + c + object.position.start_z;
        trajectoryX[i] = discretizedPoints[i] * cos(direction) + object.position.start_x;
        trajectoryY[i] = discretizedPoints[i] * sin(direction) + object.position.start_y;
        // std::cout << trajectoryX[i] << ", " << trajectoryY[i] << ", " << trajectoryZ[i] << std::endl;
    }

    
    std::ofstream outFile;
    outFile.open(filePath);
    std::vector<std::vector<double>> mergedCoordinates;

    // Ensure all vectors have the same size
    if (trajectoryX.size() == trajectoryY.size() && trajectoryY.size() == trajectoryZ.size()) {
        for (size_t i = 0; i < trajectoryX.size(); ++i) {
            std::vector<double> xyz = {trajectoryX[i], trajectoryY[i], trajectoryZ[i]};  // Create a new row with xyz coordinates
            mergedCoordinates.push_back(xyz);  // Add the new row to the merged container
        }
    }

     // Write the trajectory to the CSV file
    for (const auto& xyz : mergedCoordinates) {
        for (size_t i = 0; i < xyz.size(); ++i) {
            outFile << xyz[i];  // Write the xyz coordinate value
            if (i < xyz.size() - 1) {
                outFile << ",";  // Add a comma separator
            }
        }
        outFile << "\n";  // New line for the next point
    }

    return mergedCoordinates;
}

int main(int nArgc, char **papszArgv)
{
    // Load index data
    const char *index = NULL;
    Object ob;
    double degree = 0.0;
    ob.position.radius = 100;
    for (int iArg = 1; iArg < nArgc; iArg++)
    {
        if (strcmp(papszArgv[iArg], "--index")==0)
        {
            index = papszArgv[iArg + 1];
        }
        else if (strcmp(papszArgv[iArg], "--positionStart")==0)
        {
            ob.position.start_x = atof(papszArgv[iArg + 1]); // longitude
            ob.position.start_y = atof(papszArgv[iArg + 2]); // latitude
            ob.position.start_z = atof(papszArgv[iArg + 3]); // height
        }
        else if (strcmp(papszArgv[iArg], "--positionEnd")==0)
        {
            ob.position.end_x = atof(papszArgv[iArg + 1]);
            ob.position.end_y = atof(papszArgv[iArg + 2]);
            ob.position.end_z = atof(papszArgv[iArg + 3]);
        }
        else if (strcmp(papszArgv[iArg], "--velocity")==0)
        {
            ob.velocity.v = atof(papszArgv[iArg + 1]);
        }

        else if (strcmp(papszArgv[iArg], "--angle")==0)
        {
            degree = atof(papszArgv[iArg + 1]);
        }
    }
    if (index == NULL)
    {
        cout << "[ERROR] Wrong input!" << endl;
        Usage();
        return 0;
    }

    // Convert the angle value to radians.
    ob.velocity.angle = degree * M_PI / 180.0;
    ob.velocity.vxy = ob.velocity.v * cos(ob.velocity.angle);
    ob.velocity.vz = ob.velocity.v * sin(ob.velocity.angle);

    PJ_CONTEXT *context = proj_context_create();
    PJ *proj_from = proj_create_crs_to_crs(context, "EPSG:4326", "EPSG:2326", NULL);
    PJ_COORD start_coord = proj_coord(ob.position.start_y, ob.position.start_x, 0, 0);
    PJ_COORD end_coord = proj_coord(ob.position.end_y, ob.position.end_x, 0, 0);
    PJ_COORD start_result = proj_trans(proj_from, PJ_FWD, start_coord);
    PJ_COORD end_result = proj_trans(proj_from, PJ_FWD, end_coord);

    // delete
    proj_destroy(proj_from);
    proj_context_destroy(context);
    // The minimum value of the index range.
    double indexRange_x = -1975;
    double indexRange_y = 43;
    //Translate the coordinates to the Hong Kong coordinate system.
    ob.position.start_x = start_result.xy.y - 835000; 
    ob.position.start_y =start_result.xy.x - 815500; 
    ob.position.end_x = end_result.xy.y - 835000;   
    ob.position.end_y =end_result.xy.x - 815500;  
    // std::cout << "start_Y: " << ob.position.start_y << " and start_X: " << ob.position.start_x << std::endl;
    // std::cout << "end_Y: " << ob.position.end_y << " and end_X: " << ob.position.end_x << std::endl;

    struct timeval t1, t2, t3, t4;
    gettimeofday(&t1, NULL);
    RTree3d tree;
    tree.Load(index);
    gettimeofday(&t2, NULL);
    std::cout << "The time taken to load index: " << (t2.tv_sec - t1.tv_sec) + (double)(t2.tv_usec - t1.tv_usec) / 1000000.0 << "s"<<std::endl;

    vector<ValueType> itemList;
    //The terrain height of the start point and end point.
    double startTerrainHeight = tree.Getheight3d(indexRange_x + ob.position.start_x, indexRange_y + ob.position.start_y, HeightCallback) + 0.1;
    double endTerrainHeight = tree.Getheight3d(indexRange_x + ob.position.end_x, indexRange_y + ob.position.end_y, HeightCallback) + 0.1;
    // std::cout << "the obstacle height of the input start corrdinate is " << startTerrainHeight << endl;
    // std::cout << "the obstacle height of the input end corrdinate is " << endTerrainHeight << endl;
    
    //The actural height of the start point and end point.
    ob.position.start_z = startTerrainHeight + ob.position.start_z;
    ob.position.end_z = endTerrainHeight + ob.position.end_z;

    // Return curve trajectory points.
    std::vector<std::vector<double>> trajectory = calculateTrajectory(ob, tree);
    //cout << "trajectory: " << trajectory.size() << endl;

    std::ofstream myFile_visible, myFile_invisible;
    char visiblePointBuffer[256];
    char invisiblePointBuffer[256];
    sprintf(visiblePointBuffer, "curve_%d_%d_%d_visiblePoints.csv", static_cast<int>(ob.position.start_x), static_cast<int>(ob.position.start_y), static_cast<int>(ob.position.start_z));
    sprintf(invisiblePointBuffer, "curve_%d_%d_%d_invisiblePoints.csv", static_cast<int>(ob.position.start_x), static_cast<int>(ob.position.start_y), static_cast<int>(ob.position.start_z));
    myFile_visible.open(visiblePointBuffer);
    myFile_invisible.open(invisiblePointBuffer);
    //the number of visible points.
    int numCount = 1;
    gettimeofday(&t3, NULL);
    for (size_t i = 1; i < trajectory.size(); ++i)
    {
        const std::vector<double> &currentPoint = trajectory[i];
        const std::vector<double> &previousPoint = trajectory[i - 1];

        Rect3d search_rect(indexRange_x + previousPoint[0], indexRange_y + previousPoint[1], previousPoint[2], indexRange_x + currentPoint[0], indexRange_y + currentPoint[1], currentPoint[2]); // 遍历两个相邻的点
        if (tree.Intersect3d(search_rect.min, search_rect.max, IntersectCallback))
        {
            numCount += 1; // invisible, until the point of intersection that obstructs the view is found, return the number of points.
            break;
        }
        numCount = i + 1; //fully visible, obtain the number of points for the desired viewshed.
    }

    for (int m = 0; m < numCount; m++) // 遍历计算每个点的可视域
    {
        // Calculate interest area
        int minX = max(indexRange_x, indexRange_x + static_cast<int>(trajectory[m][0]) - ob.position.radius);
        int minY = max(indexRange_y, indexRange_y + static_cast<int>(trajectory[m][1]) - ob.position.radius);
        int maxX = min(indexRange_x + TILE_SIZE - 1, indexRange_x + static_cast<int>(trajectory[m][0]) + ob.position.radius);
        int maxY = min(indexRange_y + TILE_SIZE - 1, indexRange_y + static_cast<int>(trajectory[m][1]) + ob.position.radius);
        int width = (maxX - minX) + 1;
        int height = (maxY - minY) + 1;
        // Create viewshed array
        vs_t **viewshed = new vs_t *[height];
        for (int i = 0; i < height; i++)
            viewshed[i] = new vs_t[width];
        // start to analysis
        #pragma omp parallel for num_threads(16) schedule(dynamic)
        for (int i = 0; i < height; i++) 
        {
            for (int j = 0; j < width; j++) 
            {
                double z = tree.Getheight3d(j + minX, i + minY, HeightCallback) + 0.1;
                Rect3d search_rect(j + minX, i + minY, z, minX + ob.position.radius, minY + ob.position.radius, trajectory[m][2]);
                if (z < 0)
                    z = 0;
                if (tree.Intersect3d(search_rect.min, search_rect.max, IntersectCallback)) 
                {
                    viewshed[i][j] = 0;// invisible
                }
                else
                {
                    viewshed[i][j] = 1;// visible
                }
            }
        }

        // Writing viewshed analysis results into csv
        char viewshedBuffer[256];
        sprintf(viewshedBuffer, "./csv-curve/viewshed_%d_%d_%d.csv", static_cast<int>(trajectory[m][0]), static_cast<int>(trajectory[m][1]), static_cast<int>(trajectory[m][2]));
        writeViewshedToCsv(viewshed, width, height, viewshedBuffer);
        for (int i = 0; i < height; i++)
        {
            delete[] viewshed[i];
        }
        delete[] viewshed;
    }
    PJ_CONTEXT *context1 = proj_context_create();
    PJ * transformation =  proj_create_crs_to_crs(context1, "EPSG:2326", "EPSG:4326", NULL);// Transform the coordinates from the spatial reference system defined by EPSG:2326 to the one defined by EPSG:4326.
    gettimeofday(&t4, NULL);
    for (int i = 0; i < trajectory.size(); ++i)
    {
        trajectory[i][1] += 815500;
        trajectory[i][0] += 835000;
        PJ_COORD coordTrans=proj_coord(trajectory[i][1],trajectory[i][0],0,0);
        PJ_COORD coordTrans_result = proj_trans(transformation, PJ_FWD, coordTrans);
        //cout<<"trajectory_y: "<<coordTrans_result.xy.x<<" trajectory_x: "<<coordTrans_result.xy.y<<endl;
        trajectory[i][1]=coordTrans_result.xy.x;
        trajectory[i][0]=coordTrans_result.xy.y;
        if (i < numCount)
        {
            for (size_t j = 0; j < trajectory[i].size(); ++j)
            {
                myFile_visible << std::fixed << std::setprecision(13) << trajectory[i][j]; // 写入xyz坐标值
                if (j < trajectory[i].size() - 1)
                {
                    myFile_visible << ","; // 添加逗号分隔符
                }
            }
            myFile_visible << "\n"; // 换行以表示下一行
        }
        if (i >= numCount - 1)
        {
            for (size_t j = 0; j < trajectory[i].size(); ++j)
            {
                myFile_invisible << std::fixed << std::setprecision(13) << trajectory[i][j]; // 写入xyz坐标值
                if (j < trajectory[i].size() - 1)
                {
                    myFile_invisible << ","; // 添加逗号分隔符
                }
            }
            myFile_invisible << "\n"; // 换行以表示下一行
        }
    }
    // 清理资源
    proj_destroy(transformation);
    proj_context_destroy(context1);
    printf("The time to calculate: %lfs\n", (t4.tv_sec - t3.tv_sec) + (double)(t4.tv_usec - t3.tv_usec) / 1000000.0);
    return 0;
}
