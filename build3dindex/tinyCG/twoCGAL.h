#ifndef TWOCGAL_H
#define TWOCGAL_H

#include "tinyCG_global.h"
#include "Vec2.hpp"

//��ά���㼸���㷨

namespace tinyCG
{

class TINYCG_EXPORT twoCGAL
{
public:
    twoCGAL();

    //��֪�߶���ĳ�������ľ��룬��õ������
    static void CalPointFromLineWithDistance_D(const Vec2d & O, const Vec2d & E, double d, Vec2d& P);
    static void CalPointFromLineWithDistance_F(const Vec2f & O, const Vec2f & E, float d, Vec2f& P);


    //�������Ƿ��ཻ
    static bool RectIntersect(double minx1, double miny1,double maxx1, double maxy1,
        double minx2, double miny2, double maxx2, double maxy2);
};

}

#endif // TWOCGAL_H
