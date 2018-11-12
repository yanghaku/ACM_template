// 辗转相除法(可以用algorithm里的__gcd() )
int gcd(int a,int b){
	return b==0? a : gcd(b,a%b);
}
int lcm(int a,int b){
	return a/gcd(a,b)*b;
}

// 扩展gcd,  x*a + y*b = d = gcd(a,b)
void gcd(int a,int b,int& d,int& x,int& y){
	if(b){ gcd(b,a%b,d,y,x); y-=x*(a/b); }
	else { d=a; x=1; y=0; }
}


const int maxn=1e5+10;
bool vis[maxn];
void init(){
	//memset(vis,0,sizeof(vis));
	int m=sqrt(maxn+0.5);
	for(int i=2;i<=m;++i)if(!vis[i])
		for(int j=i*i;j<maxn;j+=i)vis[j]=1;
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
	ll ans=0;
	while(b){
		if(b&1)ans=(ans+a)%p;
		a=(a<<1)%p; //a=(a+a)%p;
		b>>=1;
	}
	return ans;
}


