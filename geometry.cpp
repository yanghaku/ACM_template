#include<cstdio>
#include<cmath>
#include<algorithm>
using namespace std;

/** 1. 浮点比较等相关函数
 */
const int maxn=1e5+10;
const double eps=1e-6;
const double INF=1e20;
const double pi=acos(-1.0);
// 1.1 浮点判号, -1为负数, 1为正数, 0表示x==0
inline int dcmp(double x){
	if(fabs(x)<eps)return 0;
	return x<0? -1: 1;
}
// 1.2 平方函数, 返回x的平方
inline double sqr(double x){ return x*x; }
// 1.3 角度与弧度相互转化
inline double degree_rad(double ang){return ang/180*pi; }
inline double rad_degree(double rad){return rad/pi*180; }

/** 2.计算几何点类与相关函数
 */
/** 2.1 点类(也是起点为原点的向量)
 *  其中除了运算符重载还有 旋转,求模,求单位向量等函数
 */
struct Point{
	double x,y;
	Point(double _x=0,double _y=0):x(_x),y(_y){}
	void input() {scanf("%lf%ld",&x,&y); }
	void output(){printf("%.2f %.2f\n",x,y); }
	bool operator==(const Point& b)const{
		return (dcmp(x-b.x)==0 && dcmp(y-b.y)==0);
	}
	bool operator<(const Point& b)const{
		return (dcmp(x-b.x)==0? dcmp(y-b.y)<0 : x<b.x);
	}
	Point operator+(const Point& b)const{
		return Point(x+b.x,y+b.y);
	}
	Point operator-(const Point& b)const{
		return Point(x-b.x,y-b.y);
	}
	Point operator*(double a){
		return Point(x*a,y*a);
	}
	Point operator/(double a){
		return Point(x/a,y/a);
	}
	double len2(){ // 返回长度的平方(与原点距离)
		return sqr(x)+sqr(y);
	}
	double len(){ // 返回与原点的距离
		return sqrt(len2());
	}
	Point change_len(double r){ // 转化为长度为r的向量
		double l=len();
		if(dcmp(l)==0)return *this; //零向量返回自身
		r /= l;
		return Point(x*r,y*r);
	}
	Point rorate_left(){ // 顺时针旋转90度
		return Point(-y,x);
	}
	Point rorate_right(){ //逆时针旋转90度
		return Point(y,-x);
	}
	Point rorate(const Point& p,double ang){ //绕点P逆时针旋转ang
		Point v=(*this)-p;
		double c=cos(ang),s=sin(ang);
		return Point(p.x+v.x*c-v.y*s,p.y+v.x*s+v.y*c);
	}
	Point unit(){ //单位法向量
		double l=len();
		return Point(-y/l,x/l);
	}
};

// 2.2 叉积
inline double det(const Point& a,const Point& b){
	return a.x*b.y-a.y*b.x;
}
// 2.3 点积
inline double dot(const Point& a,const Point& b){
	return a.x*b.x+a.y*b.y;
}
// 2.4 两点距离
inline double dis(const Point& a,const Point& b){
	return (a-b).len();
}
// 2.5 两个向量的夹角 oa 与 ob
inline double rad(Point a,Point b){
	return fabs( atan2(fabs(det(a,b)),dot(a,b)) );
}
// 2.6 判断向量平行
bool parallel(Point a,Point b){
	double p=rad(a,b);
	return dcmp(p)==0 || dcmp(p-pi)==0;
}
// 2.7 计算 pa与pb的夹角 
inline double rad(Point p,Point a,Point b){
	return fabs( atan2(fabs(det(a-p,b-p)), dot(a-p,b-p)) );
}


/** 3.直线,线段,向量
 */
struct Line{
	Point s,e;//向量起点, 终点
	double k;
};




/*******____ 三维几何____*********/
/** 
 *  1.三维点(原点为起点的向量)
 */
struct Point3{
	double x,y,z;
	Point3(double _x=0,double _y=0,double _z=0):x(_x),y(_y),z(_z){}
	void input(){
		scanf("%lf%lf%lf",&x,&y,&z);
	}
	void output(){
		printf("%.2f %.2f %.2f\n",x,y,z);
	}
	bool operator==(const Point3& b)const{
		return dcmp(x-b.x)==0&&dcmp(y-b.y)==0&&dcmp(z-b.z)==0;
	}
	bool operator<(const Point3& b)const{
		return dcmp(x-b.x)==0?(dcmp(y-b.y)==0?dcmp(z-b.z)<0:y<b.y):x<b.x;
	}
	Point3 operator-(const Point3& b)const{
		return Point3(x-b.x,y-b.y,z-b.z);
	}
	Point3 operator+(const Point3& b)const{
		return Point3(x+b.x,y+b.y,z+b.z);
	}
	Point3 operator*(const double& k)const{
		return Point3(x*k,y*k,z*k);
	}
	Point3 operator/(const double& k)const{
		return Point3(x/k,y/k);
	}
	double len(){return sqrt(x*x+y*y+z*z); }
	double len2(){return x*x+y*y+z*z; }
	// 两点距离
	double distance(const Point3& b)const{
		return sqrt(sqr(x-b.x)+sqr(y-b.y)+sqr(z-b.z));
	}
	// 点乘
	double operator * (const Point3& b)const{
		return x*b.x+y*b.y+z*b.z;
	}
	// 叉乘
	Point3 operator ^ (const Point3& b)const{
		return Point3(y*b.z-z*b.y,z*b.x-x*b.z,x*b.y-y*b.x);
	}
	// pa 与 pb 的夹角
	double rad(Point3 a,Point3 b){
		Point3 p=(*this);
		return acos( ((a-p)*(b-p)) / (a.distance(p)*b.distance(p)) );
	}
	// 将长度变为 r
	Point3 trunc(double r){
		double l=len();
		if(dcmp(l)==0)return *this;
		r /= l;
		return Point3(x*r,y*r,z*r);
	}
	// 返回单位向量
	Point3 unit(){  return *this/len(); }
	// 混合积
	double mix(const Point3& a,const Point3& b){
		return (*this)*(a^b);
	}
};

/** 2. 三维直线(线段,向量)
 */
struct Line3{
	Point3 s,e; //start -> end
	Line3(){}
	Line3(Point3 _s,Point3 _e):s(_s),e(_e){}
	bool operator==(const Line3& v)const{
		return (s==v.s)&&(e==v.e);
	}
	void input(){s.input();e.input();}
	void output(){s.output();e.output();}
	double length(){return s.distance(e);}
	// 点到直线的距离
	double dis_point_to_line(Point3 p){
		return ((e-s)^(p-s)).len()/s.distance(e);
	}
	// 点到线段的距离
	double dis_point_to_seg(Point3 p){
		if(dcmp((p-s)*(e-s))<0 || dcmp((p-e)*(s-e))<0 )
			return min(p.distance(s),e.distance(p));
		return dis_point_to_line(p);
	}
	// 返回点p到直线上的投影
	Point3 line_prog(Point3 p){
		return s+( ((e-s)*((e-s)*(p-s)))/((e-s).len2()) );
	}
	// p绕此向量旋转 ang
	Point3 rorate(Point3 p,double ang){
		if(dcmp(((s-p)^(e-p)).len()) ==0)return p;
		Point3 f1= (e-s)^(p-s);
		Point3 f2= (e-s)^(f1);
		double len=((s-p)^(e-p)).len()/s.distance(e);
		f1 = f1.trunc(len); f2=f2.trunc(len);
		Point3 h=p+f2;
		Point3 pp=h+f1;
		return h+((p-h)*cos(ang))+((pp-h)*sin(ang));
	}
	// 点在直线上
	bool point_on_seg(const Point3& p){
		return dcmp(((s-p)^(e-p)).len())==0 && dcmp((s-p)*(e-p))==0;
	}

};


/**  3.平面
 */
struct Plane{
	Point3 a,b,c,o; //平面上三个点, 和法向量
	Plane(){}
	Plane(Point3 _a,Point3 _b,Point3 _c):a(_a),b(_b),c(_c){
		o=pvec();
	}
	//  返回法向量
	Point3 pvec(){ return (b-a)^(c-a); }
	// ax+by+cz+d=0;
	Plane(double _a,double _b,double _c,double _d){
		o=Point3(_a,_b,_c);
		if(dcmp(_a) != 0)
			a=Point3((-_d-_c-_b)/_a,1,1);
		else if(dcmp(_b) !=0)
			a=Point3(1,(-_d-_c-_a)/_b,1);
		else if(dcmp(_c) !=0)
			a=Point3(1,1,(_d-_a-_b)/_c);
	}
	// 点在平面上
	bool point_on_plane(Point3 p){
		return dcmp((p-a)*o)==0;
	}
	// 两平面夹角
	double plane_angle(Plane f){
		return acos(o*f.o/(o.len()*f.o.len()));
	}
	// 平面和直线的交点, 返回值是交点个数
	int crossline(Line3 u,Point3& p){
		double x=o*(u.e-a);
		double y=o*(u.s-a);
		double d=x-y;
		if(dcmp(d)==0)return 0;
		p=((u.s*x)-(u.e*y))/d;
		return 1;
	}
	// 点到平面最近点(投影)
	Point3 point_to_plane(Point3 p){
		Line3 u=Line3(p,p+o);
		crossline(u,p);
		return p;
	}
	// 平面和平面的交线
	int crossplane(Plane f,Line3& u){
		Point3 oo=o^f.o;
		Point3 v=o^oo;
		double d=fabs(f.o*v);
		if(dcmp(d)==0)return 0;
		Point3 q=a+(v*(f.o*(f.a-a))/d);
		u=Line3(q,q+oo);
		return 1;
	}
	// 三角形的面积(需要a,b,c都被初始化)
	double area(){
		double x1=(a-b).len(),x2=(a-c).len(),x3=(b-c).len();
		double p=(x1+x2+x3)/2;
		return sqrt(p*(p-x1)*(p-x2)*(p-x3));
	}
	// 平面沿p方向前进x距离 (需要a,b,c都初始化)
	Plane go(Point3 p,double x){
		o=o.trunc(x);
		double ang=a.rad(a+o,p);
		if(dcmp(ang-pi/2) >=0)o=o*(-1);
		return Plane(a+o,b+o,c+o);
	}
};