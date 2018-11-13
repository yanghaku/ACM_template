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
			if(y==1 && x!=1 && x!=n-1) return false;
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
 *  需要gcd(),q_mul,miller_rabin;
 *  调用之前findfac()之前要将tot变为0;
 */
long long fac[400];
int tot;
// 找出一个因子
long long pollard_rho(long long x,long long c){
	long long i=1,k=2;
	srand(time(NULL));
	long long x0=rand()%(x-1)+1;
	long long y=x0;
	while(1){
		++i;
		x0=(q_mul(x0,x0,x)+c)%x;
		long long d=gcd(y-x0,x);
		if(d!=1 && d!=x)return d;
		if(y==x0)return x;
		if(i==k){ y=x0; k += k; }
	}
}
//对n进行素因子分解, 存入fac
void findfac(long long n,int k=107){
	if(n==1)return;
	if(miller_rabin(n)){
		fac[tot++]=n;
		return;
	}
	long long p=n;
	int c=k;
	while(p>=n)p=pollard_rho(p,c--);
	findfac(p,k);
	findfac(n/p,k);
}








// 辗转相除法(可以用algorithm里的__gcd() )
int gcd(int a,int b){
	return b==0? a : gcd(b,a%b);
}
int lcm(int a,int b){
	return a/gcd(a,b)*b;
}

/** 扩展gcd,  x*a + y*b = d = gcd(a,b)
 *  数据很有可能会溢出, 最好用long long
 */
void gcd(int a,int b,int& d,int& x,int& y){
	if(b){ gcd(b,a%b,d,y,x); y-=x*(a/b); }
	else { d=a; x=1; y=0; }
}


typedef long long ll;
ll q_pow(ll a,ll b,ll p){
	ll ans=1;
	while(b){
		if(b&1)ans=ans*a%p;
		a=a*a%p;
		b>>=1;
	}
	return ans;
}


ll q_mul(ll a,ll b,ll p){
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

/**
	O(1) 快速乘
	(常数很大, 而且特别大数会溢出)
*/
inline ll mul_O1(ll x,ll y){
	long long tmp=(x*y-(long long)((long double)x/p*y+1.0e-8)*p);
	return tmp<0? tmp+p : tmp;
}


/**
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


/** 逆元四种求法
  1. 扩展欧几里德
  2. 欧拉定理(或费马小定理)
  3. 线性递推打表
  4. 递归求法
 */

long long inv(long long a,long long mod){
	long long x,y,d;
	gcd(a,mod,d,x,y);
	return d==1? (x%n+n)%n : -1;
}
long long inv(long long a,long long p){// p需要是素数
	return q_pow(a,p-2,p);
}
long long inv(long long a,long long m){
	if(a==1)return 1;
	return inv(m%a,m)*(m-m/a)%m;
}
void inverse(){
	memset(inv,0,sizeof(inv));
	inv[1]=1;
	for(int i=2;i<=maxn;++i){
		inv[i]=inv[mod%i]*(mod-mod/i)%mod;
	}
}

