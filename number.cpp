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
	long long tmp=x;
	for(int i=2;i<=tmp/i;++i){
		factor[fatCnt][1]=0;
		if(tmp%i==0){
			factor[fatCnt][0]=i;
			while(tmp%i==0){
				++factor[fatCnt][1];
				tmp /= i;
			}
			++fatCnt;
		}
	}
	if(tmp != 1){
		factor[fatCnt][0]=tmp;
		factor[fatCnt][1]=1;
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


