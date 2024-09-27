/*
 * date: 2024-04-26
 * author: Mengyu Ma@National University of Defense Technology
 *         Yifan Zhang@National University of Defense Technology
 * e-mail: mamengyu10@nudt.edu.cn
 *         zhangyifan.000@nudt.edu.cn
 * description:  merge 3d indexes.
 * run:   ./Mergeindex --inputdir ../build3dindex/test_indexes --output ./indexes/test_area.3idx
 * results: output: /indexes/test_area.3idx, Total faces: 6776701, minx -1974.635391 miny 43.427783 maxx -1474.635391 maxy 543.42778
 */

#include <string>
#include <iostream>
#include <dirent.h>
#include <sys/time.h>
#include <math.h>
#include <sys/stat.h>
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
/// construct a 3D bounding box
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

// This function browses through the specified directory and collects the paths of all non-directory files.
// It takes a pointer to a directory path (dir) and a reference to a vector of strings (vFiles) to store the file paths.
void BrowseDir(char *dir, vector<string> &vFiles)
{
    // Pointer to the directory structure and entry structure.
    DIR *pDir = NULL;
    struct dirent *ent;
    // Buffer to store the path of the child (file or directory).
    char childpath[512];
    // Open the directory specified by 'dir'.
    pDir = opendir((char *)dir);
    if (pDir == NULL)
    {
        printf("Open %s failed\n", dir);
        return;
    }
    // Clear the childpath buffer.
    memset(childpath, 0, sizeof(childpath));
    while ((ent = readdir(pDir)) != NULL)
    {
        struct stat fileStat;

        sprintf(childpath, "%s/%s", dir, ent->d_name);
        int status = stat(childpath, &fileStat);
        // If the entry is not a directory, add its path to the vector of files.
        if (!S_ISDIR(fileStat.st_mode))
            vFiles.push_back(childpath);
    }

    closedir(pDir); // Close the directory
}

void Usage()
{
    printf("Usage:           [--inputdir:      input file dir            ]\n"
           "                 [--output:        index file path           ]\n");
}

int main(int nArgc, char **papszArgv)
{
    // Pointers to store the input directory and output file path.
    char *inputdir = NULL;
    char *output = NULL;
    for (int iArg = 1; iArg < nArgc; iArg++)
    {
        if (strcmp(papszArgv[iArg], "--inputdir") == 0)
        {
            inputdir = papszArgv[iArg + 1];
        }
        else if (strcmp(papszArgv[iArg], "--output") == 0)
        {
            output = papszArgv[iArg + 1];
        }
    }
    if (inputdir == NULL || output == NULL)
    {
        cout << "[ERROR] Wrong input!" << endl;
        Usage();
        return 0;
    }

    RTree3d treeMerge;
    int fcount = 0;
    // Initialize variables to store the minimum and maximum bounds.
    double minx, miny, maxx, maxy;
    minx = 100000000;
    miny = 100000000;
    maxx = -100000000;
    maxy = -100000000;
    // Vector to store the list of files in the input directory.
    vector<string> fileList;
    BrowseDir(inputdir, fileList);
    for (int i = 0; i < fileList.size(); i++)
    {
        RTree3d rtree;
        rtree.Load((fileList[i]).c_str());
        printf("%s\n", (fileList[i]).c_str());

        RTree3d::Iterator it;
        for (rtree.GetFirst(it);
             !rtree.IsNull(it);
             rtree.GetNext(it))
        {
            ValueType value = rtree.GetAt(it);

            double boundsMin[3] = {0, 0, 0};
            double boundsMax[3] = {0, 0, 0};
            it.GetBounds(boundsMin, boundsMax); // Get the bounds of the current entry.
            if (minx > boundsMin[0])
                minx = boundsMin[0];
            if (maxx < boundsMax[0])
                maxx = boundsMax[0];
            if (miny > boundsMin[1])
                miny = boundsMin[1];
            if (maxy < boundsMax[1])
                maxy = boundsMax[1];
            treeMerge.Insert(boundsMin, boundsMax, value); // Insert the current entry into the treeMerge RTree3d instance.
            fcount++;
        }
    }

    printf("Output Index: %s\n", output);
    treeMerge.Save(output);
    printf("Range :minx %lf miny %lf maxx %lf maxy %lf\n", minx, miny, maxx, maxy);
    printf("Total faces:%d\n", fcount);

    return 0;
}
