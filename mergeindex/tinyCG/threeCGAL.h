#ifndef THREECGAL_H
#define THREECGAL_H

#include "tinyCG_global.h"
#include "Vec3.hpp"

//��ά���㼸���㷨

namespace tinyCG
{

class TINYCG_EXPORT threeCGAL
{
public:
    threeCGAL();

    //�������ֱ�ߵľ���
    static double CalDistancePointAndLine(Vec3d &point, Vec3d &lineBegin, Vec3d &lineEnd);

    //�����������ķ�����
    static void CalNormal(const Vec3d& v1, const Vec3d& v2, const Vec3d& v3, Vec3d &vn);
};

}

#endif // THREECGAL_H
