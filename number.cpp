/**
 *	1. 素数相关
 */

/** 1.1埃氏筛 O(nlogn)
 *	notPrime[] ture表示非素数
 */
const int maxn=1e5+10;
bool notPrime[maxn]; 
void getPrime(){
	//memset(notPrime,false,sizeof(notPrime));
	notPrime[0]=notPrime[1]=true;
	for(int i=2;i<maxn;++i){
		if(!notPrime[i]){
			if(i>maxn/i)continue;	//防止溢出
			for(int j=i*i;j<maxn;j+=i)notPrime[j]=true;
		}
	}
}

/** 1.2线性筛 O(n)
 *	得到升序的素数数组, 个数保存在 tot 中
 *  并且可以得到 notPrime[]布尔数组, true表示非素数
 */
const int maxn=1e5+10;
bool notPrime[maxn];
int tot, prime[maxn];
void getPrime(){
	//memset(notPrime,0,sizeof(prime));
	//memset(prime,0,sizeof(prime));
	tot=0;
	for(int i=2;i<maxn;++i){
		if(!notPrime[i])prime[tot++] = i;
		for(int j=0;j<tot && i*prime[j]<maxn ;++j){
			notPrime[ i*prime[j] ]=true;
			if(i%prime[j]==0)break;
		}
	}
}

/** 1.3区间筛 (poj2689)
 *	需要1.2的线性筛预处理出来prime[]数组 和 tot
 *	复杂度 O(R-L)
 *	notPrime2[i] 保存区间上的 L+i 是否为素数
 *  prime2[i] 保存第i个素数, cnt为总个数
 */
bool notPrime2[1000010];
int cnt, prime2[1000010];
void getPrime2(int L,int R){
	memset(notPrime2,0,sizeof(notPrime2));
	if(L<2)L=2;
	for(int i=0; i<tot && (long long)prime[i]*prime[i]<=R;++i){
		int s=L/prime[i]+(L%prime[i]>0);
		if(1==s)s=2;
		for(int j=s;(long long)j*prime[i]<=R;++j)
			if((long long)j*prime[i]>=L)
				notPrime2[j*prime[i]-L]=true;
	}
	cnt=0;
	for(int i=0;i<=R-L;++i)
		if(!notPrime2[i])
			prime2[cnt++]=i+L;
}

/** 1.4.1 miller_rabin素性检验 (<1e8时)
	只需2,3,5和一些特例就能判断1e8以内的素数
	960946321是1e9时的特例
*/
bool miller_rabin(int n){
    if(2==n || 3==n || 5==n || 7==n)return true;
    if(1==n|| 0==(n&1) ||0==n%3||0==n%5||4097==n
        ||1048577==n||16777217==n||25326001==n)return false;
    long long x,d;
    long long s[]={2,3,5};  
    for(int k=0;k<3;++k){
        int i=ceil( log(n-1.0)/log(2.0) )-1;
        d=1;
        for(;i>=0;--i){
            x=d;
            d=(d*d)%n;
            if(d==1 && x!=1 && x!=n-1) return false;
            if(((n-1) & (1<<i)) >0) d=(d*s[k])%n;
        }
        if(d!=1)return false;
    }
    return true;
}

/** 1.4.2 miller_rabin素性检验 复杂度O(logn) (n<2^63)
	是素数返回True, (有可能是伪素数)
	不是素数就返回false (一定不是素数)
	对于任意奇数n>2和正整数s. 出错的概率<=2^(-S)
	(需要调用快速幂q_pow和快速乘q_mul)
*/
bool miller_rabin(long long n,int S=50){
	if(n==2)return true;
	if(n<2 || !(n&1))return false;//偶数或者小于2
	long long a,x,y,u=n-1;
	int t=0;
	while( (u&1)==0 ){ u>>=1; ++t; }
	srand(time(NULL)); // 一定不能忘记置随机数种子
	for(int i=0;i<S;++i){
		a= rand()%(n-1)+1;
		x=q_pow(a,u,n);
		for(int j=0;j<t;++j){
			y=q_mul(x,x,n);
			if(y==1 && x!=1 && x!=n-1)return false;
			x=y;
		}
		if(x!=1)return false;
	}
	return true;
}
/** 1.4.3 大数素数测试
 *  最好用java的BigInteger, bool isProbablePrime()
 *  当java没法过的时候再用这个
 */
const int MAXL=4;
#define M10 1000000000
#define Z10 9
const int zero[MAXL-1]={0};
struct bnum{
	int data[MAXL];
	void read(){
		memset(data,0,sizeof(data));
		char buf[32];
		scanf("%s",buf);
		int len=(int)strlen(buf);
		int i=0,k;
		while(len >= Z10){
			for(k=len-Z10;k<len;++k)data[i]=data[i]*10+buf[k]-'0';
			++i;
			len -= Z10;
		}
		if(len>0){
			for(k=0;k<len;++k)data[i]=data[i]*10+buf[k]-'0';
		}
	}
	bool operator==(const bnum& x){ return memcmp(data,x.data,sizeof(data))==0; }
	bnum& operator=(const int x){
		memset(data,0,sizeof(data));
		data[0]=x;
		return *this;
	}
	bnum operator+(const bnum& x){
		int i,carry=0;
		bnum ans;
		for(int i=0;i<MAXL;++i){
			ans.data[i]=data[i]+x.data[i]+carry;
			carry = ans.data[i]/M10;
			ans.data[i] %= M10;
		}
		return ans;
	}
	bnum operator-(const bnum& x){
		int i,carry=0;
		bnum ans;
		for(i=0;i<MAXL;++i){
			ans.data[i]=data[i]-x.data[i]-carry;
			if(ans.data[i]<0){
				ans.data[i] += M10;
				carry=1;
			}else carry=0;
		}
		return ans;
	}
	bnum operator%(const bnum& x){
		for(int i=MAXL-1;i>-1;--i){
			if(data[i]<x.data[i])return *this;
			else if(data[i]>x.data[i])break;
		}
		return ((*this) - x);
	}
	bnum& div2(){
		int carry=0,tmp;
		for(int i=MAXL-1;i>-1;--i){
			tmp=data[i]&1;
			data[i]=(data[i]+carry)>>1;
			carry=tmp*M10;
		}
		return *this;
	}
	bool is_odd(){
		return (data[0]&1)==1;
	}
	bool is_zero(){
		for(int i=0;i<MAXL;++i){
			if(data[i])return false;
		}
		return true;
	}
};
void mul_mod(bnum& a0,bnum& b0,bnum& p,bnum& ans){
	bnum tmp=a0,b=b0;
	ans=0;
	while(!b.is_zero()){
		if(b.is_odd())ans=(ans+tmp)%p;
		tmp=(tmp+tmp)%p;
		b.div2();
	}
}
void pow_mod(bnum& a0,bnum& b0,bnum& p,bnum& ans){
	bnum tmp=a0,b=b0;
	ans=1;
	while(!b.is_zero()){
		if(b.is_odd())mul_mod(ans,tmp,p,ans);
		mul_mod(tmp,tmp,p,tmp);
		b.div2();
	}
}
bool MillerRabinTest(bnum& p,int iter){
	int i,small=0,j,d=0;
	for(i=1;i<MAXL;++i){
		if(p.data[i])break;
	}
	if(i==MAXL){
		if(p.data[0]<2)return false;
		if(p.data[0]==2)return true;
		small=1;
	}
	if(!p.is_odd())return false;
	bnum a,s,m,one,pd1;
	one=1;
	s=pd1=p-one;
	while(!s.is_odd()){
		s.div2(); ++d;
	}
	for(int i=0;i<iter;++i){
		a=rand();
		if(small)a.data[0]=a.data[0]%(p.data[0]-1)+1;
		else{
			a.data[1]=a.data[0]/M10;
			a.data[0]%=M10;
		}
		if(a==one)continue;
		pow_mod(a,s,p,m);
		for(j=0;j<d && !(m==one) && !(m==pd1);++j)mul_mod(m,m,p,m);
		if(!(m==pd1) && j>0)return false;
	}
	return true;
}

/**
 *  2.欧拉函数
 */

/** 2.1 求欧拉函数单个值 O(n^0.5)
 */
long long phi(long long n){
	long long ans=n;
	int m=sqrt(n+0.5);
	for(int i=2;i<=m;++i){
		if(n%i==0){
			ans-=ans/i;
			while(n%i==0)n /= i;
		}
	}
	if(n>1)ans -= ans/n;
	return ans;
}

/** 2.2 筛法求欧拉函数 O(nlogn)
 *  筛除maxn以内的所有欧拉函数值 
 */
const int maxn=1e5+10;
int phi[maxn];
void getPhi(){
	memset(phi,0,sizeof(phi));
	phi[1]=1;
	for(int i=2;i<maxn;++i)if(!phi[i])
		for(int j=i;j<maxn;j += i){
			if(!phi[j])phi[j]=j;
			phi[j]=phi[j]/i*(i-1);
		}
}

/** 2.3 线性筛 O(n)
 *  同时得到欧拉函数phi[]和素数数组
 *  notPrime[] 记录了每一个值的素性
 *  tot 是prime[] 数组的个数
 */
const int maxn=1e5+10;
bool notPrime[maxn];
int phi[maxn];
int tot, prime[maxn];
void phi_and_prime(){
	memset(notPrime,0,sizeof(notPrime));
	phi[1]=1;
	tot=0;
	for(int i=2;i<maxn;++i){
		if(!notPrime[i]){
			prime[ tot++ ]=i;
			phi[i]=i-1;
		}
		for(int j=0;j<tot;++j){
			if(i*prime[j]>=maxn)break;
			notPrime[i*prime[j]]=true;
			if(i%prime[j]==0){
				phi[i*prime[j]]=phi[i]*prime[j];
			}
			else{
				phi[i*prime[j]]=phi[i]*(prime[j]-1);
			}
		}
	}
}


/**
 *  3.合数分解
 */
/** 3.1 朴素分解算法 O(n^0.5)
 *  factor第一个保存的素因子,第二个保存的次数
 *  可以用prime试除优化
 */
long long factor[100][2];
int fatCnt;
int getFactor(long long x){
	fatCnt=0;
	for(int i=2;i<=x/i;++i){
		factor[fatCnt][1]=0;
		if(x%i==0){
			factor[fatCnt][0]=i;
			while(x%i==0){
				++factor[fatCnt][1];
				x /= i;
			}
			++fatCnt;
		}
	}
	if(x != 1){
		factor[fatCnt][0]=x;
		factor[fatCnt++][1]=1;
	}
	return fatCnt;
}

/** 3.2 pollard_rho算法进行质因素分解
 *  结果保存在fac[]里 (返回的是无序的!)
 *  需要q_mul,miller_rabin;
 *  调用之前findfac()之前要将tot变为0;
 */
long long fac[200];
int tot;
// 找出一个因子
long long pollard_rho(long long x,long long c){
	long long i=1,k=2;
	long long x0=rand()%(x-1)+1;
	long long y=x0;
	while(1){
		++i;
		x0=(q_mul(x0,x0,x)+c)%x;
		long long d=abs(__gcd(y-x0,x));
		if(d!=1 && d!=x)return d;
		if(y==x0)return x;
		if(i==k){ y=x0; k += k; }
	}
}
//对n进行素因子分解, 存入fac
void findfac(long long n){
	if(n==1)return;
	if(miller_rabin(n)){
		fac[tot++]=n;
		return;
	}
	long long p=n;
	while(p>=n)p=pollard_rho(p,rand()%(n-1)+1);
	findfac(p);	findfac(n/p);
}


/**
 *	4.gcd, 扩展gcd, 快速乘, 快速幂
 */

// 4.1 辗转相除法(可以用algorithm里的__gcd() )
int gcd(int a,int b){
	return b==0? a : gcd(b,a%b);
}
int lcm(int a,int b){
	return a/gcd(a,b)*b;
}

/** 4.2 扩展gcd,  x*a + y*b = d = gcd(a,b)
 *  数据很有可能会溢出, 最好用long long
 */
void gcd(int a,int b,int& d,int& x,int& y){
	if(b){ gcd(b,a%b,d,y,x); y-=x*(a/b); }
	else { d=a; x=1; y=0; }
}

// 4.3 快速gcd (在位数比较大时的更相减损法)
int qgcd(int a,int b){
	if(a==0)return b;
	if(b==0)return a;
	if(!(a&1) && !(b&1))
		return qgcd(a>>1,b>>1)<<1;
	else if(!(b&1))
		return qgcd(a,b>>1);
	else if(!(b&1))
		return qgcd(a>>1,b);
	else
		return qgcd(abs(a-b),min(a,b));
}

// 4.4 快速幂 O(logb)
// !!!!! 当p为long long的时候, 要把里面的乘法改为快速乘
long long q_pow(long long a,long long b,long long p){
	long long ans=1;
	while(b){
		if(b&1)ans=ans*a%p;
		a=a*a%p;
		b>>=1;
	}
	return ans;
}

// 4.5 log快速乘取模, O(logb)
long long q_mul(long long a,long long b,long long p){
	a%=p;
	b%=p;
	ll ans=0;
	while(b){
		if(b&1){
			ans+=a;
			if(ans>=p)ans-=p;
		}
		a<<=1;
		if(a>=p)a-=p;
		b>>=1;
	}
	return ans;
}

/** 4.6 O(1) 快速乘
 *  (常数很大, 而且特别大数会溢出)
 */
inline long long mul_O1(long long x,long long y,long long p){
	long long tmp=(x*y-(long long)((long double)x/p*y+1.0e-8)*p);
	return tmp<0? tmp+p : tmp;
}


/** 4.7 O1 快速乘 (常数极小)
	只适用于 X86_64 位的评测机
	需要关闭编译器优化
	并且只能运算 (unsigned long long)
*/
#pragma GCC optimize("O0")
typedef unsigned long long ull;
inline ull q_mul(ull a,ull b,ull p){
    ull ans;
    __asm__ __volatile__( 
        "movq %1,%%rax\n\t"
        "mulq %2\n\t"
        "pushq %%rax\n\t"
        "movq %%rdx,%%rax\n\t"
        "xorq %%rdx,%%rdx\n\t"
        "divq %3\n\t"
        "popq %%rax\n\t"
        "divq %3\n\t"
        "movq %%rdx,%0\n\t"
        :"=r"(ans)
        :"r"(a),"r"(b),"r"(p)
    );
    return ans;
}

/**
 *  5.逆元
 */

// 5.1 扩展欧几里德
long long inv(long long a,long long mod){
	long long x,y,d;
	gcd(a,mod,d,x,y);
	return d==1? (x%n+n)%n : -1;
}

// 5.2 欧拉定理
long long inv(long long a,long long p){// p需要是素数
	return q_pow(a,p-2,p);
}

// 5.3 递归
long long inv(long long a,long long m){
	if(a==1)return 1;
	return inv(m%a,m)*(m-m/a)%m;
}

// 5.4 线性递推
void inverse(){
	memset(inv,0,sizeof(inv));
	inv[1]=1;
	for(int i=2;i<=maxn;++i){
		inv[i]=inv[mod%i]*(mod-mod/i)%mod;
	}
}

//5.5 阶乘的逆元
const int maxn=1e5+10;
long long fac[maxn],inv[maxn];
void init(long long p){
	fac[0]=1;
	for(int i=1;i<maxn;++i)fac[i]=fac[i-1]*i%p;
	inv[maxn-1]=q_pow(fac[maxn-1],p-2,p);
	for(int i=maxn-2;i>=0;--i)inv[i]=inv[i+1]*(i+1)%p;
}


/**
 *  6.模线性方程(组)
 */

/** 6.1 求模线性方程的特解 O(logn)
 * a*x = b (mod m)的一个特解
 * 如果没有解会返回m
 */
long long solve(long long a,long long b,long long m){
	long long x,y,d;
	gcd(a,m,d,x,y);  //扩展gcd求逆元和最大公因数
	if(b%d==0){
		x=(x%m+m)%m;
		return x*(b/d)%(m/d);
	}
	else return m;
}

/** 6.2 求模线性方程 [0,m) 的区间所有解 O(logn)
 *  结果存在vector里, 
 */
vector<long long>ans;
void solve(long long a,long long b,long long m){
	ans.clear();
	long long x,y,d;
	gcd(a,m,d,x,y);
	if(b%d==0){
		x=(x%m+m)%m;
		y=m/d;
		ans.push_back(x*(b/d)%y);
		for(int i=1;i<d;++i)
			ans.push_back((ans[i-1]+y)%m);
	}
}

/** 6.3 中国剩余定理(互素时) O(nlogM)
 *  求解 x= a[i] (mod m[i]) 的方程组
 *  要求m[i]两两互素 , 返回一个特解
 */
long long China(long long a[],long long m[],int n){
	long long M=1,ans=0;
	long long x,y,d,tm;
	for(int i=0;i<n;++i)M*=m[i];
	for(int i=0;i<n;++i){
		tm=M/m[i];
		gcd(tm,m[i],d,x,y);
		ans = (ans + tm*x*a[i])%M;
	}
	return (ans+M)%M;
}

/** 6.4 中国剩余定理(不互素时) O(nlogM)
 *  求解 x = a[i] (nod m[i]), m[i]可以不互素
 *  如果无解就返回-1 ;
 */
long long China_ex(long long a[],long long m[],int n){
	if(1==n&&0==a[0])return m[0];
	long long ans=a[0],_lcm=m[0];
	bool flag=true;
	long long x,y,_gcd,tmp;
	for(int i=1;i<n;++i){
		gcd(_lcm,m[i],_gcd,x,y);
		if((a[i]-ans)%_gcd){ flag=false; break; }
		tmp=_lcm*( ((a[i]-ans)/_gcd*x)%(m[i]/_gcd) );
		_lcm=_lcm/_gcd*m[i];
		ans=(ans+tmp)%_lcm;
		if(ans<0)ans+=_lcm;
	}
	return flag? ans : -1;
}

/** 7.数论其他知识
 */

/** 7.1 原根
 *  需要预先筛除prime[]数组,
 *  getFactor()求出质因子
 *  接口: primitive_root(P) 返回P的最小原根
 */
int primitive_root(int P){
	if(2==P)return 1;
	getFactor(P-1);
	for(int g=2;g<P;++g){
		bool flag=true;
		for(int i=0;i<fatCnt;++i){
			if(q_pow(g,(P-1)/factor[i][0],P)==1){
				flag=false;
				break;
			}
		}
		if(flag)return g;
	}
	return -1;
}

/** 7.2 勒让德符号 (p需要是奇素数) O(logn)
 *  当为1 : d是模 p的平方剩余
 *  当为-1: d不是模p的平常剩余
 *  当为0: d整除p
 */
int Legendre(long long d,long long p){
	int coef= (d>0)? 1: (((p-1)%4==0)? 1:-1);
	d=(d>0)? d: -d;
	d%=p;
	if(q_pow(d,(p-1)>>1,p)==1)return coef;
	else return -coef;
}

/** 7.3 平方剩余  O((logn)^2)
 *  求 x^2=a(mod m)的最小整数解
 *  无解返回 -1
 */
long long modsqr(long long a,long long m){
	long long b,k,i,x;
	a%=m;
	if(2==m)return a%m;
	if(q_pow(a,(m-1)>>1,m)==1){
		if(m%4 != 3){
			for(b=1;q_pow(b,(m-1)>>1,m)==1;++b);
			i=(m-1)>>1;k=0;
			while(true){
				i /= 2; k /= 2;
				long long h1=q_pow(a,i,m);
				long long h2=q_pow(b,k,m);
				if((h1*h2+1)%m==0)k+=(m-1)>>1;
				if(i&1)break;
			}
			long long t1=q_pow(a,(i+1)>>1,m);
			long long t2=q_pow(b,k>>1,m);
			x=q_mul(t1,t2,m);
		}
		else x=q_pow(a,(m+1)>>2,m);
		if(x*2>m)x=m-x;
		return x;
	}
	return -1;
}

/** 7.4离散对数 O(n^0.5logn)
 *  求a^x=b(mod n)的最小整数解x, 无解返回-1
 *  (n为素数)
 */
long long log_mod(long long a,long long b,long long n){
	long long e=1;
	int m=(int)sqrt(n+0.5);
	long long v=inv(q_pow(a,m,n),n);
	map<long long,int>x;
	x[1]=0;
	for(int i=1;i<m;++i){
		e=q_mul(e,a,n);
		if(!x.count(e))x[e]=i;
	}
	for(int i=0;i<m;++i){
		if(x.count(b))return i*m+x[b];
		b=q_mul(b,v,n);
	}
	return -1;
}

/** 7.5 N次剩余
 *  求解 x^N=a(mod P)的所有整数解x, 其中P为素数
 *  找到原根g后, 所有P-1都可以与 g^i 建立一一对应关系
 *  然后转换成 g^(N*y) = g^t (mod (p-1))
 *  就变为 y*N=t(mod (p-1))   (t可以用离散对数求解)
 *
 *  需要写好的函数: 4.2扩展gcd, 4.4q_pow, 4.5q_mul
 *        3.1质因子分解, 7.1原根, 7.4离散对数
 */
vector<long long>residue(long long N,long long a,long long p){
	long long g=primitive_root(p);
	long long m=log_mod(g,a,p);
	vector<long long>ret;
	if(0==a){ ret.push_back(0); return ret; }
	if(-1==m)return ret;
	long long A=N,B=p-1,C=m,x,y;
	long long d;
	gcd(A,B,d,x,y);
	if(C%d != 0)return ret;
	x=x*(C/d)%B;
	long long delta=B/d;
	for(int i=0;i<d;++i){
		x=((x+delta)%B+B)%B;
		ret.push_back((q_pow(g,x,p)));
	}
	sort(ret.begin(),ret.end());
	ret.erase(unique(ret.begin(),ret.end()), ret.end());
	return ret;
}

/** 
 * 8. mobius 莫比乌斯函数
 */
// 8.1 单独求值 O(n^0.5)
int mu(int n){
	int cnt,k=0;
	for(int i=2;i*i<=n;++i){
		if(n%i)continue;
		cnt=0;
		++k;
		while(!(n%i)){
			n/=i;
			++cnt;
		}
		if(cnt>1)return 0;
	}
	if(n != 1)++k;
	return (k&1)? -1:1;
}

// 8.2 莫比乌斯函数递推 O(nlogn)
const int maxn=1e5+10;
int mu[maxn];
void mobius(){
	int delta, target;
	for(int i=1;i<maxn;++i){
		target= (i==1)? 1:0;
		delta=target - mu[i];
		mu[i]=delta;
		for(int j=i<<1;j<maxn;j += i)
			mu[j] += delta;
	}
}

// 8.3 线性筛(与素数筛差不多)
const int maxn=1e5+10;
bool notPrime[maxn];
int cnt, prime[9600], mu[maxn];
void mobius()
{
    memset(vis,0,sizeof(vis));
    mu[1]=1;
    cnt=0;
    for(int i=2;i<maxn;++i)
    {
        if(!vis[i])
        {
            prime[cnt++]=i;
            mu[i]=-1;
        }
        for(int j=0;j<cnt&&i*prime[j]<maxn;++j)
        {
            vis[i*prime[j]]=1;
            if(i%prime[j]==0)
            {
                mu[i*prime[j]]=0;
                break;
            }
            else mu[i*prime[j]]= -mu[i];
        }
    }
}


/**
 *  9.组合数取模
 */

//  9.1. 当0<=n,m<=1000时,靠杨辉三角递推
const int maxn=1001;
long long C[maxn][maxn];
void init(){
	C[0][0]=1;
    for(int i=1;i<maxn;++i)
        for(int j=0;j<maxn;++j){
            C[i][j]=C[i-1][j]+C[i-1][j-1];
            if(C[i][j]>=P)C[i][j] -= P;
        }
}

// 9.2. 当0<=n,m<2e6时 用阶乘的逆元求解
//		阶乘的逆元可以预处理出来, 也可以每次查询的时候再计算    
const int maxn=1e5+10;
long long fac[maxn],inv[maxn];
void init(){
    fac[0]=fac[1]=1;
    for(int i=2;i<maxn;++i)fac[i]=fac[i-1]*i%P;
}
long long C(int n,int k){
    long long x,y,d;
    long long tmp=(fac[k]*fac[n-k])%P;
    gcd(tmp,P,d,x,y);
    y=fac[n]*(x%P)%P;
    if(y<0)y += P;
    return y;
}

// 9.3 当1<=m,n<=1e18, 2<=P<=1e6时(P为素数) Lucas定理
const long long P=100003;
long long fac[P];
long long q_pow(long long a){
    long long ans=1;
    long long b=P-2;
    while(b){
        if(b&1)ans=ans*a%P;
        a=a*a%P;
        b>>=1;
    }
    return ans;
}
void init(){
    fac[0]=fac[1]=1;
    for(int i=2;i<P;++i)
        fac[i]=fac[i-1]*i%P;
}
long long C(long long n,long long m){
    if(m>n)return 0;
    return q_pow(fac[m]*fac[n-m]%P)*fac[n]%P;
}
long long Lucas(long long n,long long m){
    if(m==0)return 1;
    return (C(n%P,m%P)*Lucas(n/P,m/P))%P;
}


// 9.3. 当1<=m,n<=1e18, 2<=P<=1e6(P为合数时) 拓展Lucas定理加中国剩余定理
//      需要gcd,q_pow,  调用ll ans=solve(n,m,P);
typedef long long ll;
ll inv(ll a,ll mod){
	if(!a)return 0;
	ll x,d,y;
	gcd(a,mod,d,x,y);
	x%=mod;
	if(x<=0)x+=mod;
	return x;
}
ll Mul(ll n,ll pi,ll pk){
	if(!n)return 1LL;
	ll ans=1LL;
	if(n/pk){
		for(ll i=2;i<=pk;++i)
			if(i%pi)ans=ans*i%pk;
		ans=q_pow(ans,n/pk,pk);
	}
	for(ll i=2;i<=n%pk;++i)
		if(i%pi)ans=ans*i%pk;
	return ans*Mul(n/pi,pi,pk)%pk;
}
ll C(ll n,ll m,ll Mod,ll pi,ll pk){
	if(m>n)return 0;
	ll a=Mul(n,pi,pk),b=Mul(m,pi,pk),c=Mul(n-m,pi,pk);
	ll k=0,ans;
	for(ll i=n;i;i /= pi)k += i/pi;
	for(ll i=m;i;i /= pi)k -= i/pi;
	for(ll i=n-m;i;i /= pi)k -= i/pi;
	ans=a*inv(b,pk)%pk*inv(c,pk)%pk*q_pow(pi,k,pk)%pk;
	return ans*(Mod/pk)%Mod*inv(Mod/pk,pk)%Mod;
}
ll solve(long long n,long long m,long long P){
	long long ans=0;
	for(ll x=P,i=2;i<=P;++i)
		if(x%i==0){
			ll pk=1;
			while(x%i==0)pk *= i,x /= i;
			ans=(ans+C(n,m,P,i,pk))%P;
		}
	return ans;
}

/** 10. 矩阵类
 */

// 10.1. 矩阵类, 使用前务必调用clear()清零
const int maxN=1010;
const int maxM=1010;
struct Matrix{
	int n,m;
	int a[maxN][maxM];
	inline void clear(){ n=m=0; memset(a,0,sizeof(a)); }
	Matrix operator+(const Matrix& b)const{
		Matrix tmp;tmp.n=n;tmp.m=m;
		for(int i=0;i<n;++i)
			for(int j=0;j<m;++j)
				tmp.a[i][j]=a[i][j]+b.a[i][j];
		return tmp;
	}
	Matrix operator-(const Matrix& b)const{
		Matrix tmp;tmp.n=n;tmp.m=m;
		for(int i=0;i<n;++i)
			for(int j=0;j<m;++j)
				tmp.a[i][j]=a[i][j]-b.a[i][j];
		return tmp;
	}
	Matrix operator*(const Matrix& b)const{
		Matrix tmp;tmp.clear();
		tmp.n=n;tmp.m=b.m;
		for(int i=0;i<n;++i)
			for(int j=0;j<b.m;++j)
				for(int k=0;k<m;++k)
					tmp.a[i][j] += a[i][k]*b.a[k][j];
		return tmp;
	}
};

// 10.2矩阵的逆 O(n^3)
// 输入 A原矩阵, C逆矩阵 n矩阵的阶数
// 当矩阵不是满秩的, 返回0
// 当算整数取模的时候,只需要改动前两个内联函数和变成乘逆元即可
inline vector<double>operator*(const vector<double>&a,double b){
	int n=a.size();
	vector<double> res(n,0);
	for(int i=0;i<n;++i)res[i]=a[i]*b;
	return res;
}
inline vector<double>operator-(const vector<double>&a,const vector<double>& b){
	int n=a.size();
	vector<double> res(n,0);
	for(int i=0;i<n;++i)res[i]=a[i]-b[i];
	return res;
}
const double eps=1e-8;
inline int inverse(vector<double>A[],vector<double>C[],int n){
	for(int i=0;i<n;++i){
		C[i]=vector<double>(n,0);
		C[i][i]=1;
	}
	for(int i=0;i<n;++i){
		for(int j=i;j<n;++j)if(fabs(A[j][i]>eps)){
			swap(A[i],A[j]);
			swap(C[i],C[j]);
			break;
		}
		if(fabs(A[i][i])<eps)return 0;//矩阵不是满秩的
		C[i]=C[i]*(1.0/A[i][i]);
		A[i]=A[i]*(1.0/A[i][i]);
		for(int j=0;j<n;++j)
			if(j!=i&&fabs(A[j][i]>eps)){
				C[j]=C[j]-C[i]*A[j][i];
				A[j]=A[j]-A[i]*A[j][i];
			}
	}
	return 1;
}

/** 10.3矩阵快速幂加速递推
 *  考虑到递推式 f[x]=a[n-1]*f[x-1]+a[n-2]*f[x-2]+....+a[0]*f[x-n]
 *  可以变为:
 *    [ 0  ,  1 ,  0 ,...... 0 ]   [f[x-n]  ]
 *    [ 0  ,  0 ,  1 ,...... 0 ]   [f[x-n+1]]
 *  A=[ ...................... ] B=[........]
 *    [ 0  ,  0 ,  0 ,...... 1 ]   [ f[x-2] ]
 *    [a[0],a[1],a[2],...a[n-1]]   [ f[x-1] ]
 *  构造出A,B矩阵用快速幂计算即可
 *  (需要矩阵类)
 */
// a[],b[]为递推的系数, n为矩阵大小, t为递推次数 O(n^3logt)
int solve(int a[],int b[],int n,int t){
	Matrix M,F,E;
	M.clear(); M.n=M.m=n;
	E.clear(); E.n=E.m=n;
	F.clear(); F.n=n;F.m=1;
	for(int i=0;i<n-1;++i)M.a[i][i+1]=1;
	for(int i=0;i<n;++i){
		M.a[n-1][i]=a[i];
		F.a[i][0]=b[i];
		E.a[i][1]=1;
	}
	if(t<n)return F.a[t][0];
	for(t -= n-1;t;t>>=1){
		if(t&1)E=M*E;
		M=M*M;
	}
	F=E*F;
	return F.a[n-1][0];
}

/** 
 *  11. 高斯消元
 */

/** 11.1 高斯消元Gauss O(n^3)
 *  a[][maxn]为方程组对应的矩阵, n为未知数上的个数
 *  l,ans存储解, l[]表示是否为自由元 True表示不是自由元
 *  返回解空间的维数
 */
const int maxn=105;
const double eps=1e-8;
int Gauss(double a[][maxn],bool l[],double ans[],const int& n){
	int res=0,r=0;
	for(int i=0;i<n;++i)l[i]=false;
	for(int i=0;i<n;++i){
		for(int j=r;j<n;++j)if(fabs(a[j][i])>eps){
			for(int k=i;k<=n;++k)swap(a[j][k],a[r][k]);
			break;
		}
		if(fabs(a[r][i])<eps){ ++res; continue; }
		for(int j=0;j<n;++j)
			if(j!=r && fabs(a[j][i])>eps){
				double tmp=a[j][i]/a[r][i];
				for(int k=i;k<=n;++k)a[j][k] -= tmp*a[r][k];
			}
		l[i]=true;
		++r;
	}
	for(int i=0;i<n;++i)if(l[i])
		for(int j=0;j<n;++j)if(fabs(a[j][i])>eps)
			ans[i]=a[j][n]/a[j][i];
	return res;
}
/** 11.2 列主元高斯消元
 *  求解 a[][]*x[]=b[]
 *  返回是否有唯一解,若有解就存在b[]中
 */
const int maxn=100;
int Gauss(int n,double a[][maxn],double b[]){
	int i,j,row=0;
	double maxP,tmp;
	for(int k=0;k<n;++k){ //枚举每一列
		for(maxP=0,i=k;i<n;++i){ //对于第k列,往下找最大的一行
			if(fabs(a[i][k]) > fabs(maxP))maxP=a[row=i][k];
		}
		if(fabs(maxP)<eps)return 0;//如果下面全是0,那么就没有唯一解
		if(row != k){
			for(j=k;j<n;++j)swap(a[k][j],a[row][j]);
			swap(b[k],b[row]);
		}
		for(j=k+1;j<n;++j){ //枚举k行后的每一行
			a[k][j] /= maxP;
			for(i=k+1;i<n;++i)a[i][j] -= a[i][k]*a[k][j];
		}
		b[k]/=maxP;
		for(i=n-1;i>-1;--i)
			for(j=i+1;j<n;++j)b[i] -= a[i][j]*b[j];
	}
	return 1;
}



/** 12.分数类
 *  通过分子,分母进行构造
 *  重载了 +,-,*,/,<,==
 */
struct Fraction{
	long long num,den;
	Fraction(long long n=0,d=0){
		if(d<0){ n=-n; d=-d;}
		assert(d!=0);
		long long g=__gcd(abs(n),d);
		num=n/g;den=d/g;
	}
	Fraction operator+(const Fraction& b)const{
		return Fraction(num*b.den+den*b.num,den*b.den);
	}
	Fraction operator-(const Fraction& b)const{
		return Fraction(num*b.den-den*b.num,den*b.den);
	}
	Fraction operator*(const Fraction& b)const{
		return Fraction(num*b.num,den*b.den);
	}
	Fraction operator/(const Fraction& b)const{
		return Fraction(num*b.den,den*b.num);
	}
	bool operator<(const Fraction& b)const{
		return num*b.den < den*b.num;
	}
	bool operator==(const Fraction& b)const{
		return num*b.den == den*b.num;
	}
};


/** 13. 二分、三分计算
 */
// 13.1 二分计算, test()为自己的判断函数
const double eps=1e-15;
double Bsearch(double l,double r){
	while(r-l > eps){
		double mid=(l+r)/2;
		if(test(mid))l=mid;
		else r=mid;
	}
	return (l+r)/2;
}

/** 13.2 三分
 */
// 13.2.1 三等分法  求F(x)最大值, 返回x
double Tsearch(double l,double r){
	double midl=0,midr=0;
	while(r-l > eps){
		midl=(2*l + r)/3;
		midr=(2*r + l)/3;
		if(F(midl)<F(midr))l=midl;
		else r=midr;
	}
	return midl;
}

// 13.2.2 midmid法
double Tsearch(double l,double r){
	double mid=0,midmid=0;
	while(r-l > eps){
		mid=(r+l)/2;
		midmid=(mid+r)/2;
		if(F(mid)>F(midmid))r=midmid;
		else l=mid;
	}
	return mid;
}

// 13.2.3 优选法  (求最大值)
const double eps=1e-9;
const double cef=(sqrt(5.0)-1.0)*0.5;

double Tsearch_s(double l,double r){
	double midl=r-(r-l)*cef;
	double midr=l+(r-l)*cef;
	while(r-l > eps){
		if(F(midl)<F(midr)){
			l=midl; midl=midr;
			midr=l+ (r-l)*cef;
		}
		else{
			r=midr; midr=midl;
			midl=r- (r-l)*cef;
		}
	}
	return midr;
}

/** 14. 求积分
 */

/** 14.1 Simpson 自适应辛普森方法
 *  eps为精度要求, F()为要积分的函数
 *  直接调用 asr(积分下限,积分上限) 即可 
 */
const double eps=1e-10;
inline double simpon(double l,double r){
	double mid=(l+r)/2;
	return (F(l)+4*F(mid)+F(r))*(r-l)/6;
}
double asr(double l,double r,double A){
	double mid=(l+r)/2;
	double L=simpon(l,mid),R=simpon(mid,r);
	if(fabs(L+R-A)<=15*eps)return L+R+(L+R-A)/15.0;
	else return asr(l,mid,L)+asr(mid,r,R);
}
double asr(double l,double r){
	return asr(l,r,simpon(l,r));
}

/** 14.2 Romberg 龙贝格方法 O((logK)^2)
 */
const double eps=1e-8;
double Romberg(double L,double R){
	vector<double> t;
	double h=R-L,last,curr;
	int k=1,i=1;
	t.push_back(h* (F(L)+F(R))/2 );
	while(true){
		last=t.back();
		curr=0;
		double x=L+h/2;
		for(int j=0;j<k;++j,x += h)curr += F(x);
		curr=(t[0]+h*curr)/2;
		double k1=4.0/3.0, k2=1.0/3.0;
		for(int j=0;j<i;++j){
			double tmp=k1*curr-k2*t[j];
			t[j]=curr;curr=tmp;
			k2 /= 4*k1-k2;
			k1=k2+1;
		}
		t.push_back(curr);
		k *= 2;
		h /= 2;
		++i;
		if(fabs(last-curr)<eps)break; 
	}
	return t.back();
}

/** 15. FFT 快速傅里叶变换
 *  例:hdu1402高精乘
 */
const double PI=acos(-1.0);
struct Complex{
	double x,y; //x+yi
	Complex(double _x=0.0,double _y=0.0):x(_x),y(_y){}
	Complex operator-(const Complex& b)const{
		return Complex(x-b.x,y-b.y);
	}
	Complex operator+(const Complex& b)const{
		return Complex(x+b.x,y+b.y);
	}
	Complex operator*(const Complex& b)const{
		return Complex(x*b.x-y*b.y,x*b.y+y*b.x);
	}
};
/* 进行FFT与IFFT前的反转变换
 * 位置i与i二进制反转后的位置互换
 * len 必须为2的幂
 */
void change(Complex y[],int len){
	int i,j,k;
	for(i=1,j=len/2;i<len-1;++i){
		if(i<j)swap(y[i],y[j]);
		k=len/2;
		while(j>=k){
			j -= k;
			k /= 2;
		}
		if(j<k)j += k;
	}
}
/* len必须为2的幂
 * on==1 :DFT, on==-1: IDFT
 */
void fft(Complex y[],int len,int on){
	change(y,len);
	for(int h=2;h<=len;h<<=1){
		Complex wn(cos(-on*2*PI/h),sin(-on*2*PI/h));
		for(int j=0;j<len;j+=h){
			Complex w(1,0);
			for(int k=j;k<j+h/2;++k){
				Complex u=y[k];
				Complex t=w*y[k+h/2];
				y[k]=u+t;
				y[k+h/2]=u-t;
				w=w*wn;
			}
		}
	}
	if(on == -1)
		for(int i=0;i<len;++i)
			y[i].x /= len;
}
const int maxn=200010;
Complex x1[maxn],x2[maxn];
char str1[maxn/2],str2[maxn/2];
int sum[maxn];
int main(){
	while(scanf("%s%s",str1,str2)!=EOF){
		int len1=strlen(str1);
		int len2=strlen(str2);
		int len=1;
		while(len<len1*2 || len<len2*2)len<<=1;
		for(int i=0;i<len1;++i)
			x1[i]=Complex(str1[len1-i-1]-'0',0);
		for(int i=len1;i<len;++i)
			x1[i]=Complex(0,0);
		for(int i=0;i<len2;++i)
			x2[i]=Complex(str2[len2-i-1]-'0',0);
		for(int i=len2;i<len;++i)
			x2[i]=Complex(0,0);
		fft(x1,len,1);
		fft(x2,len,1);
		for(int i=0;i<len;++i)
			x1[i]=x1[i]*x2[i];
		fft(x1,len,-1);
		for(int i=0;i<len;++i)
			sum[i]=(int)(x1[i].x+0.5);
		for(int i=0;i<len;++i){
			sum[i+1] += sum[i]/10;
			sum[i] %= 10;
		}
		len=len1+len2-1;
		while(sum[len]<=0&&len>0)--len;
		for(int i=len;i>=0;--i)
			printf("%c",sum[i]+'0');
		printf("\n");
	}
	return 0;
}

/** 16 NTT 快速数论变换, 



/** 
 *   18. 高阶方程求根 O(N^3*logK)
 *  对于方程 a[n]*x^n+a[n-1]*x^n-1+....+a[1]*x+a[0]=0,求出方程的所有实数解
 *  因为利用了递归和vector,时间效率不高,
 */
const double eps=1e-12;
const double inf=1e+12;
inline int sign(double x){ return (x< -eps)? -1: x>eps; } 
inline double get(const vector<double>& coef,double x){
	double e=1,s=0;
	for(int i=0;i<coef.size();++i,e*=x)s+=coef[i]*e;
	return s;
}
double find(const vector<double>& coef,int n,double lo,double hi){
	double sign_lo,sign_hi;
	if((sign_lo=sign(get(coef,lo))) == 0)return lo;
	if((sign_hi=sign(get(coef,hi))) == 0)return hi;
	if( sign_lo*sign_hi>0 )return inf;
	for(int step=0;step<100&&hi-lo>eps;++step){
		double m=(lo+hi)*0.5;
		int sign_mid=sign(get(coef,m));
		if(sign_mid==0)return m;
		if(sign_lo*sign_mid<0)hi=m;
		else lo=m;
	}
	return (lo+hi)*0.5;
}
vector<double> solve(vector<double>coef,int n){
	vector<double>ret;
	if(n==1){
		if(sign(coef[1]))ret.push_back(-coef[0]/coef[1]);
		return ret;
	}
	vector<double> dcoef(n); //求导
	for(int i=0;i<n;++i)dcoef[i]=coef[i+1]*(i+1);
	vector<double> droot=solve(dcoef,n-1);
	droot.insert(droot.begin(),-inf);
	droot.push_back(inf);
	for(int i=0;i+1<droot.size();++i){
		double tmp=find(coef,n,droot[i],droot[i+1]);
		if(tmp<inf)ret.push_back(tmp);
	}
	return ret;
}

/** 19. 多项式求根(牛顿法)
 *  c[]多项式系数,n为多项式度数, 求在[a,b]的根
 *  输出根 (要保证[a,b]有根)
 */
double fabs(double x){ return (x<0)? -x:x; }
double F(int m,double c[],double x){
	double p=c[m];
	for(int i=m;i>0;--i)p=p*x+c[i-1];
	return p;
}
int newton(double x0,double *r,double c[],double cp[],int n,double a,double b,double eps){
	int MAX_ITERATION=1000;
	int i=1;
	double x1,x2,fp,eps2=eps/10.0;
	x1=x0;
	while(i<MAX_ITERATION){
		x2=F(n,c,x1);
		fp=F(n-1,cp,x1);
		if(fabs(fp)<1e-9 && fabs(x2)>1.0 )return 0;
		x2=x1-x2/fp;
		if(fabs(x1-x2)<eps2){
			if(x2<a || x2>b)return 0;
			*r=x2;
			return 1;
		}
		x1=x2;
		++i;
	}
	return 0;
}
double polynomial_root(double c[],int n,double a,double b,double eps){
	double *cp;
	int i;
	double root;
	cp=(double*)calloc(n,sizeof(double));
	for(i=n-1;i>-1;--i)cp[i]=(i+1)*c[i+1];
	if(a>b){
		root=a;
		a=b;
		b=root;
	}
	if((!newton(a,&root,c,cp,n,a,b,eps)) && (!newton(b,&root,c,cp,n,a,b,eps)))
		newton((a+b)*0.5,&root,c,cp,n,a,b,eps);
	free(cp);
	if(fabs(root)<eps)return fabs(root);
	else return root;
}







/**
 *  19.其他
 */
// 19.1 进制转换
// 将x进制的串s转换为y进制的串
string transform(int x,int y,string s){
	int len=s.size(),sum=0;
	string res="";
	for(int i=0;i<len;++i){
		if(s[i]=='-')continue;
		if(s[i]>='0'&&s[i]<='9')
			sum=sum*x+s[i]-'0';
		else if(s[i]>='A'&&s[i]<='Z')
			sum=sum*x+s[i]-'A'+10;
		else sum=sum*x+s[i]-'a'+10+26;
	}
	while(sum){
		char tmp=sum%y;
		sum /= y;
		if(tmp<=9)tmp += '0';
		else if(tmp<=36)tmp +='A'-10;
		else tmp += 'A'-10-26;
		res=tmp+res;
	}
	if(res.size()==0)res="0";
	if(s[0]=='-')res='-'+res;
	return res;
}

// 19.2 格雷码 O(2^n)
// 给一个n, 求一个0~2^n-1的排列, 使得相邻两项(包括首尾)的二进制只有一位不同
vector<int> initGray(int n){
	vector<int>res;res.resize(1<<n);
	for(int i=0;i<(1<<n);++i)
		res[i]=(i^(i>>1));
	return res;
}

