#include<cstdio>
#include<cmath>
#include<algorithm>
using namespace std;

typedef long long ll;
const ll P = 1000000007LL;

int main() {
	double ans = (double)P*P;
	ll y;
	for(ll x=sqrt(P); x<P; x++) {
		double na = max((double)x*x, (double)P/x/2*P/x/2*(P-(ll)x*x%P));
		if(na < ans) {
			ans = na;
			y = x;
		}
		if((ll)x*x > ans) break;
	}
	printf("%lld\n%lld\n", y, P-(ll)y*y%P);
	return 0;
}
