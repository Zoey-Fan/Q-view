/*
 * date: 2024-04-26
 * author: Mengyu Ma@National University of Defense Technology
 *         Yifan Zhang@National University of Defense Technology
 * e-mail: mamengyu10@nudt.edu.cn
 *         zhangyifan.000@nudt.edu.cn
 * description:  3d index query function.
 */
#ifndef QUERY3DRTREE_H
#define QUERY3DRTREE_H

#include <functional>
#include "RTree.h"

//#define RTREE_TEMPLATE template<typename DATATYPE, typename ELEMTYPE, int NUMDIMS, typename ELEMTYPEREAL, int TMAXNODES, int TMINNODES>
#define RTREE_TEMPLATE template<class DATATYPE, class ELEMTYPE, int NUMDIMS, class ELEMTYPEREAL, int TMAXNODES, int TMINNODES>
#define RTREE_QUAL RTree<DATATYPE, ELEMTYPE, NUMDIMS, ELEMTYPEREAL, TMAXNODES, TMINNODES>
#define MAXH 1000000000
#define MINH -1000000000

// This function checks for intersection between a 3D segment and elements within a R-tree.
// It takes the start and end points of the segment as arrays of three elements (ELEMTYPE),
// and a callback function that is called when an intersection is detected.
RTREE_TEMPLATE
bool RTREE_QUAL::Intersect3d(const ELEMTYPE s_start[3], const ELEMTYPE s_end[3], 
    std::function<bool (const DATATYPE&, const ELEMTYPE [3], const ELEMTYPE [3])> callback) const
{
    Seg3d seg;
    for (int axis = 0; axis < NUMDIMS; ++axis) {
        seg.s_start[axis] = s_start[axis];
        seg.s_end[axis] = s_end[axis];
    }
    bool IntersectFlag = false;// Initialize a flag to track whether an intersection has been found.

    // Recursively check for intersections starting from the root of the 3D R-tree.
    Intersect3d(m_root, &seg, IntersectFlag, callback);
    
    return IntersectFlag;
}

// This function retrieves the height at a specific 3D location using a 3D R-tree structure.
// It takes the x and y coordinates of the location and a callback function that processes the data.
// The function returns the height value at the given location.
RTREE_TEMPLATE
ELEMTYPE RTREE_QUAL::Getheight3d(const ELEMTYPE x, const ELEMTYPE y, 
    std::function<ELEMTYPE (const DATATYPE&, const ELEMTYPE, const ELEMTYPE)> callback) const
{
    double Height = MINH - 1;
     // Recursively search for the height at the given x, y coordinates starting from the root of the 3D R-tree.
    Getheight3d(m_root, x, y, Height, callback);
    return Height;
}

// This function checks for overlap between two 2D line segments.
// It takes the coordinates of the endpoints of the first segment (x1, y1, x2, y2)
// and the coordinates of the endpoints of the second segment (tx1, ty1, tx2, ty2).
// The function returns true if the segments overlap, otherwise false.
RTREE_TEMPLATE
inline bool RTREE_QUAL::OverlapSeg2d(ELEMTYPE x1, ELEMTYPE y1, ELEMTYPE x2, ELEMTYPE y2, 
    ELEMTYPE tx1, ELEMTYPE ty1, ELEMTYPE tx2, ELEMTYPE ty2) const
{
    if ((x1 > x2 ? x1 : x2) < (tx1 < tx2 ? tx1 : tx2) ||
        (y1 > y2 ? y1 : y2) < (ty1 < ty2 ? ty1 : ty2) ||
        (tx1 > tx2 ? tx1 : tx2) < (x1 < x2 ? x1 : x2) ||
        (ty1 > ty2 ? ty1 : ty2) < (y1 < y2 ? y1 : y2))
        return false;
    // Check for overlap using the cross product method.
    if ((((x1 - tx1) * (ty2 - ty1) - (y1 - ty1) * (tx2 - tx1)) *
        ((x2 - tx1) * (ty2 - ty1) - (y2 - ty1) * (tx2 - tx1))) > 0 ||
        (((tx1 - x1) * (y2 - y1) - (ty1 - y1) * (x2 - x1)) *
        ((tx2 - x1) * (y2 - y1) - (ty2 - y1) * (x2 - x1))) > 0)
        return false;

     // If none of the above conditions are met, the segments overlap.
    return true;
}

// This function checks for an overlap between a 2D line segment and a rectangle.
// It takes the coordinates of the endpoints of the line segment (x1, y1, x2, y2)
// and the minimum and maximum x and y coordinates of the rectangle (minx, miny, maxx, maxy).
// The function returns true if there is an overlap, otherwise false.

RTREE_TEMPLATE
inline bool RTREE_QUAL::OverlapSegRect2d(ELEMTYPE x1, ELEMTYPE y1, ELEMTYPE x2, ELEMTYPE y2, 
    ELEMTYPE minx, ELEMTYPE miny, ELEMTYPE maxx, ELEMTYPE maxy) const
{
    // Check if any endpoint of the segment is inside the rectangle.
    if ((x1 >= minx && x1 <= maxx && y1 >= miny && y1 <= maxy) ||
        (x2 >= minx && x2 <= maxx && y2 >= miny && y2 <= maxy))
        return true;

    // Check for overlap between the segment and the rectangle using the OverlapSeg2d function.
    // This checks for proper overlap, not just touching at the corners.
    return OverlapSeg2d(x1, y1, x2, y2, minx, miny, maxx, maxy) || 
           OverlapSeg2d(x1, y1, x2, y2, minx, maxy, maxx, miny);
}

// This function checks for overlap between a 3D segment and a 3D rectangle.
// It takes pointers to a 3D segment (a_segA) and a 3D rectangle (a_rectB).
// The function returns true if there is an overlap, otherwise false.
RTREE_TEMPLATE
bool RTREE_QUAL::Overlap3d(Seg3d* a_segA, Rect* a_rectB) const
{
    ASSERT(a_segA && a_rectB);// Ensure that the pointers to the segment and rectangle are valid.
    
    // Check for overlap between the segment's XY plane projection and the rectangle's XY plane projection.
    if (!(OverlapSegRect2d(a_segA->s_start[0], a_segA->s_start[1], 
        a_segA->s_end[0], a_segA->s_end[1], 
        a_rectB->m_min[0], a_rectB->m_min[1], 
        a_rectB->m_max[0], a_rectB->m_max[1])))
        return false;
    // Check for overlap between the segment's XZ plane projection and the rectangle's XZ plane projection.
    if (!(OverlapSegRect2d(a_segA->s_start[0], a_segA->s_start[2], 
        a_segA->s_end[0], a_segA->s_end[2], 
        a_rectB->m_min[0], a_rectB->m_min[2], 
        a_rectB->m_max[0], a_rectB->m_max[2])))
        return false;
    // Check for overlap between the segment's YZ plane projection and the rectangle's YZ plane projection.
    if (!(OverlapSegRect2d(a_segA->s_start[1], a_segA->s_start[2], 
        a_segA->s_end[1], a_segA->s_end[2], 
        a_rectB->m_min[1], a_rectB->m_min[2], 
        a_rectB->m_max[1], a_rectB->m_max[2])))
        return false;
    // If all projections overlap, then the 3D segment and rectangle overlap.
    return true;
}

// This function checks for intersection between a 3D segment and the nodes of an R-tree.
RTREE_TEMPLATE
bool RTREE_QUAL::Intersect3d(Node* a_node, Seg3d* a_seg, bool& IntersectFlag, 
    std::function<bool (const DATATYPE&, const ELEMTYPE [NUMDIMS], const ELEMTYPE [NUMDIMS])> callback) const
{
    // Ensure that the node and segment pointers are valid.
    ASSERT(a_node);
    ASSERT(a_node->m_level >= 0);
    ASSERT(a_seg);
    // If the current node is an internal node (not a leaf), recursively check its children.
    if (a_node->IsInternalNode()) {
        // Check for overlap between the segment and the bounding rectangle of the current branch.
        for (int index = 0; index < a_node->m_count; ++index) {
            if (Overlap3d(a_seg, &a_node->m_branch[index].m_rect)) {
                if (!Intersect3d(a_node->m_branch[index].m_child, a_seg, IntersectFlag, callback)) {
                    return false;// If an intersection is found, stop the search.
                }
            }
        }
    } else {
        // If the current node is a leaf node, check for intersection with the data it contains.
        for (int index = 0; index < a_node->m_count; ++index) {
            // Check for overlap between the segment and the bounding rectangle of the current branch.
            if (Overlap3d(a_seg, &a_node->m_branch[index].m_rect)) {
                // Retrieve the data associated with the current branch.
                DATATYPE& id = a_node->m_branch[index].m_data;
                // Call the callback function with the data, start, and end points of the segment.
                if (callback && !callback(id, a_seg->s_start, a_seg->s_end)) {
                    IntersectFlag = true;// Set the intersection flag if the callback returns false (indicating an intersection).
                    return false;// Stop the search.
                }
            }
        }
    }
    // If no intersections were found, continue the search.
    return true;
}

// This function retrieves the maximum height at a specific 2D location (x, y) within an 3D R-tree.
RTREE_TEMPLATE
bool RTREE_QUAL::Getheight3d(Node *a_node, ELEMTYPE x, ELEMTYPE y, ELEMTYPE &Height,
                             std::function<ELEMTYPE(const DATATYPE &, const ELEMTYPE, const ELEMTYPE)> callback) const
{
    ASSERT(a_node);
    ASSERT(a_node->m_level >= 0);
    bool found = false;
    if (a_node->IsInternalNode())
    {
        // If the current node is an internal node (not a leaf), recursively check its children.
        for (int index = 0; index < a_node->m_count; ++index)
        {
            // Check if the location (x, y) is within the bounding rectangle of the current branch.
            if (!(x < (&a_node->m_branch[index].m_rect)->m_min[0]) &&
                !(y < (&a_node->m_branch[index].m_rect)->m_min[1]) &&
                !(x > (&a_node->m_branch[index].m_rect)->m_max[0]) &&
                !(y > (&a_node->m_branch[index].m_rect)->m_max[1]))
            {
                if (!Getheight3d(a_node->m_branch[index].m_child, x, y, Height, callback))
                {
                    found = true;
                }
            }
        }
    }
    else
    {
        // If the current node is a leaf node, check for the maximum height at the location (x, y).
        for (int index = 0; index < a_node->m_count; ++index)
        {
            //~ printf("adadfadsfasdffafd\n");
            //~ printf("%lf %lf a_rectB %lf  %lf  %lf  %lf  \n",a_seg->s_start[0],a_seg->s_start[1],(&a_node->m_branch[index].m_rect)->m_min[0],(&a_node->m_branch[index].m_rect)->m_min[1],(&a_node->m_branch[index].m_rect)->m_max[0],(&a_node->m_branch[index].m_rect)->m_max[1]);

            if (!(x < (&a_node->m_branch[index].m_rect)->m_min[0]) &&
                !(y < (&a_node->m_branch[index].m_rect)->m_min[1]) &&
                !(x > (&a_node->m_branch[index].m_rect)->m_max[0]) &&
                !(y > (&a_node->m_branch[index].m_rect)->m_max[1]))
            {
                // Retrieve the data associated with the current branch.
                DATATYPE &id = a_node->m_branch[index].m_data;
                if (callback)
                {
                    ELEMTYPE h = callback(id, x, y); // h is the actural height of xy
                    // If the height returned by the callback is greater than the current maximum, update the maximum height.
                    if (h > Height)
                    {
                        Height = h;
                        found = true;// Stop the search
                    }
                }
            }
        }
    }
    return found;
}
#endif // QUERY3DRTREE_H
