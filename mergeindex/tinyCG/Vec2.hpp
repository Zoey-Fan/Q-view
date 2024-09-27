#ifndef VEC2_HPP
#define VEC2_HPP

#include <iostream>

//��άʸ��

namespace tinyCG
{

template <class T>
class Vec2
{
public:
    //ʸ������
    typedef T value_type;
    value_type _v[2];

    //���캯��
    Vec2() {_v[0]=T(0); _v[1]=T(0);}
    Vec2(value_type x,value_type y) { _v[0]=x; _v[1]=y; }

    //��ֵ
    inline void set( value_type x, value_type y ) { _v[0]=x; _v[1]=y; }

    //���ط����жϴ�С
    inline bool operator == (const Vec2& v) const { return _v[0]==v._v[0] && _v[1]==v._v[1]; }
    inline bool operator != (const Vec2& v) const { return _v[0]!=v._v[0] || _v[1]!=v._v[1]; }
    inline bool operator <  (const Vec2& v) const
    {
        if (_v[0]<v._v[0]) return true;
        else if (_v[0]>v._v[0]) return false;
        else return (_v[1]<v._v[1]);
    }

    //������
    inline value_type& operator [] (int i) { return _v[i]; }
    inline value_type operator [] (int i) const { return _v[i]; }

    //�ж�ֵ�Ƿ�����
    inline bool valid() const { return !isNaN(); }
    inline bool isNaN() const { return _isnan(_v[0]) || _isnan(_v[1]); }

    //���
    inline value_type operator * (const Vec2& rhs) const
    {
        return _v[0]*rhs._v[0]+_v[1]*rhs._v[1];
    }

    //��������
    inline const Vec2 operator * (value_type rhs) const
    {
        return Vec2(_v[0]*rhs, _v[1]*rhs);
    }
    //��������
    inline Vec2& operator *= (value_type rhs)
    {
        _v[0]*=rhs;
        _v[1]*=rhs;
        return *this;
    }

    //����һ������
    inline const Vec2 operator / (value_type rhs) const
    {
        return Vec2(_v[0]/rhs, _v[1]/rhs);
    }
    //����һ������
    inline Vec2& operator /= (value_type rhs)
    {
        _v[0]/=rhs;
        _v[1]/=rhs;
        return *this;
    }

    //ʸ�����
    inline const Vec2 operator + (const Vec2& rhs) const
    {
        return Vec2(_v[0]+rhs._v[0], _v[1]+rhs._v[1]);
    }
    //ʸ�����
    inline Vec2& operator += (const Vec2& rhs)
    {
        _v[0] += rhs._v[0];
        _v[1] += rhs._v[1];
        return *this;
    }

    //ʸ�����
    inline const Vec2 operator - (const Vec2& rhs) const
    {
        return Vec2(_v[0]-rhs._v[0], _v[1]-rhs._v[1]);
    }
    //ʸ�����
    inline Vec2& operator -= (const Vec2& rhs)
    {
        _v[0]-=rhs._v[0];
        _v[1]-=rhs._v[1];
        return *this;
    }

    //�������
    inline const Vec2 componentMultiply(const Vec2& rhs) const
    {
        return Vec2(_v[0]*rhs[0], _v[1]*rhs[1]);
    }

    //�������
    inline const Vec2 componentDivide( const Vec2& rhs) const
    {
        return Vec2(_v[0]/rhs[0], _v[1]/rhs[1]);
    }

    //ȡ��
    inline const Vec2 operator - () const
    {
        return Vec2 (-_v[0], -_v[1]);
    }

    //ʸ���ĳ���
    inline value_type length() const
    {
        return sqrt( _v[0]*_v[0] + _v[1]*_v[1] );
    }

    //ʸ�����ȵ�ƽ��
    inline value_type length2( void ) const
    {
        return _v[0]*_v[0] + _v[1]*_v[1];
    }

    //��һ��������������������ǰ���ȡ�
    inline value_type normalize()
    {
        value_type norm = length();
        if (norm>0.0)
        {
            value_type inv = 1.0/norm;
            _v[0] *= inv;
            _v[1] *= inv;
        }
        return( norm );
    }

    //getter and setter
    inline value_type x() const { return _v[0]; }
    inline value_type y() const { return _v[1]; }
    inline value_type u() const { return _v[0]; }
    inline value_type v() const { return _v[1]; }

	inline value_type& x() { return _v[0]; }
	inline value_type& y() { return _v[1]; }
	inline value_type& u() { return _v[0]; }
	inline value_type& v() { return _v[1]; }
	
    //���
    friend std::ostream & operator<<(std::ostream & out, Vec2 & A)
    {
        out <<"("<< A.x() <<","<< A.y() <<")";
        return out;
    }
};


typedef Vec2<float> Vec2f;
typedef Vec2<double> Vec2d;
typedef Vec2<int> Vec2i;
typedef Vec2<size_t> Vec2ul;

}

#endif // VEC2_HPP
