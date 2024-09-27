#include "threeCGAL.h"

namespace tinyCG
{

threeCGAL::threeCGAL()
{

}

double threeCGAL::CalDistancePointAndLine(Vec3d &point, Vec3d &lineBegin, Vec3d &lineEnd)
{
    //ֱ�߷�������
    Vec3d n = lineEnd -lineBegin;

    //ֱ����ĳһ����������������
    Vec3d m = point - lineBegin;

    return (n ^ m).length() / n.length();
}

//�����������ķ�����
void threeCGAL::CalNormal(const Vec3d& v1, const Vec3d& v2, const Vec3d& v3, Vec3d &vn)
{
    //v1(n1,n2,n3);
    //ƽ�淽��: na * (x �C n1) + nb * (y �C n2) + nc * (z �C n3) = 0 ;
    double na = (v2.y()-v1.y())*(v3.z()-v1.z())-(v2.z()-v1.z())*(v3.y()-v1.y());
    double nb = (v2.z()-v1.z())*(v3.x()-v1.x())-(v2.x()-v1.x())*(v3.z()-v1.z());
    double nc = (v2.x()-v1.x())*(v3.y()-v1.y())-(v2.y()-v1.y())*(v3.x()-v1.x());

    //ƽ�淨����
    vn.set(na,nb,nc);
}

}
