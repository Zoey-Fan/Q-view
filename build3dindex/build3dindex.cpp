/*
 * date: 2024-04-17
 * author: Mengyu Ma@National University of Defense Technology
 *         Yifan Zhang@National University of Defense Technology
 * e-mail: mamengyu10@nudt.edu.cn
 *         zhangyifan.000@nudt.edu.cn
 * description:  create 3D indexes.
 * run: run.py
 * result: tile_23_26.3idx (Total faces:1774685  50.218376s, minx -1974.635391 miny 43.427783 maxx -1724.635391 maxy 293.427783)
 *         tile_23_27.3idx (Total faces:2001545  56.762683s, minx -1974.635391 miny 293.427783 maxx -1724.635391 maxy 543.427783)
 *         tile_24_26.3idx (Total faces:1843398  51.952792s, minx -1724.635391 miny 43.427783 maxx -1474.635391 maxy 293.427783)
 *         tile_24_27.3idx (Total faces:1157073  31.960064s, minx -1724.635391 miny 293.427783 maxx -1474.635391 maxy 543.42778)
 */

#include <string>
#include <iostream>
#include <sys/time.h>
#include <math.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/boost/graph/graph_traits_Delaunay_triangulation_2.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include "query3dRtree.h"

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Projection_traits = CGAL::Projection_traits_xy_3<Kernel>;
using Point_2 = Kernel::Point_2;
using Point_3 = Kernel::Point_3;
using Segment_3 = Kernel::Segment_3;
// Triangulated Irregular Network
using TIN = CGAL::Delaunay_triangulation_2<Projection_traits>;
using Mesh = CGAL::Surface_mesh<Point_3>;

using namespace std;

typedef tuple<double, double, double, double, double, double, double, double, double> ValueType;
typedef RTree<ValueType, double, 3> RTree3d;

// A 3D rectangular bounding box defined by minimum and maximum coordinates along the x, y, and z axes.
struct Rect3d
{
    Rect3d() {}
    // Constructor to initialize the bounding box with specific minimum and maximum coordinates.
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
// Find the maximum and minimum values.
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

void Usage()
{
    printf("Usage:           [--threshold:  Threshold]\n"
           "                 [--input:      input file path                  ]\n"
           "                 [--output:     index file path                  ]\n");
}

int main(int nArgc, char **papszArgv)
{
    int format = 0;
    const char *input = NULL;
    const char *output = NULL;
    for (int iArg = 1; iArg < nArgc; iArg++)
    {
        if (strcmp(papszArgv[iArg], "--input") == 0)
        {
            input = papszArgv[iArg + 1];
        }
        else if (strcmp(papszArgv[iArg], "--output") == 0)
        {
            output = papszArgv[iArg + 1];
        }
    }
    if (input == NULL || output == NULL)
    {
        std::cout << "[ERROR] Wrong input!" << endl;
        Usage();
        return 0;
    }

    ifstream fs(input);
    string sline, s0;
    int vIndex = 0;
    vector<vector<double>> posList;
    RTree3d tree;
    int faceCount = 0;
    double minx, miny, maxx, maxy;
    minx = 1000000;
    miny = 1000000;
    maxx = -1000000;
    maxy = -1000000;
    struct timeval t1, t2;
    gettimeofday(&t1, NULL);
    // Read an OBJ file.
    while (getline(fs, sline))
    {
        if (sline[0] == 'v')
        {
            if (sline[1] == 'n')
            {
                continue;
            }
            else if (sline[1] == 't')
            {
                continue;
            }
            else
            {
                istringstream ins(sline);
                vector<double> posIter(3);
                ins >> s0 >> posIter[0] >> posIter[1] >> posIter[2];
                posIter[0] = posIter[0];
                posIter[1] = posIter[1];
                posIter[2] = posIter[2];
                // Find bounding of these verties.
                if (minx > posIter[0])
                    minx = posIter[0];
                if (maxx < posIter[0])
                    maxx = posIter[0];
                if (miny > posIter[1])
                    miny = posIter[1];
                if (maxy < posIter[1])
                    maxy = posIter[1];
                posList.push_back(posIter);
                vIndex++;
            }
        }
        // Build a bounding box for each triangular face and add it to the index tree.
        if (sline[0] == 'f')
        {
            istringstream inFace(sline);
            inFace >> s0;
            int vIndex;
            vector<int> faceIter(3);

            for (int i = 0; i < 3; i++)
            {
                s0.clear();
                vIndex = 0;
                inFace >> s0;
                for (int j = 0; s0[j] != '/'; j++)
                    vIndex = vIndex * 10 + (s0[j] - 48);
                faceIter[i] = vIndex;
            }
            for (int i = 0; i < 3; i++)
                faceIter[i]--;
            double Px[3], Py[3], Pz[3];
            double Pminx, Pminy, Pminz, Pmaxx, Pmaxy, Pmaxz;
            for (int i = 0; i < 3; i++)
            {
                Px[i] = posList[faceIter[i]][0];
                Py[i] = posList[faceIter[i]][1];
                Pz[i] = posList[faceIter[i]][2];
            }
            find((double *)Px, 3, Pmaxx, Pminx);
            find((double *)Py, 3, Pmaxy, Pminy);
            find((double *)Pz, 3, Pmaxz, Pminz);
            Rect3d Prect(Pminx, Pminy, Pminz, Pmaxx, Pmaxy, Pmaxz);
            tree.Insert(Prect.min, Prect.max, make_tuple(Px[0], Py[0], Pz[0], Px[1], Py[1], Pz[1], Px[2], Py[2], Pz[2]));
            faceCount++;
        }
    }

    printf("Output Index: %s\n", output);
    printf("Range :minx %lf miny %lf maxx %lf maxy %lf\n", minx, miny, maxx, maxy);

    std::cout << "Total faces:" << faceCount << "\n";
    gettimeofday(&t2, NULL);
    printf("%lfs \n", (t2.tv_sec - t1.tv_sec) + (double)(t2.tv_usec - t1.tv_usec) / 1000000.0);
    tree.Save(output);

    return 0;
}
