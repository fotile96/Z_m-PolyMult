/* Author : Jiatai Huang
 * Email  : fotile13@gmail.com
 * 
 * */
#include<cstdio>
#include<cassert>
#include<complex>
#include<cstring>
#include<algorithm>
#include<cstdlib>
#include<ctime>
#include<cmath>
#include<map>

#define M_PIl          3.141592653589793238462643383279502884L

#define ASSERT_ON

/* #define USE_LONG_DOUBLE */

using namespace std;

#ifdef USE_LONG_DOUBLE
typedef long double llf;
#else
typedef double llf;
#endif
typedef complex<llf> comp;

typedef long long ll;
const ll P = 2147483647, M = 364889, __I2 = P - M * M % P, P2 = P / 2;
#ifdef USE_LONG_DOUBLE
const llf __I = sqrtl(__I2);
#else
const llf __I = sqrt(__I2);
#endif

void butterfly(comp *a, int n) {
	for(int i=0, j=0; i<n; i++) {
		if(j>i) swap(a[i], a[j]);
		int k=n;
		while(j & (k>>=1)) j &= ~k;
		j |= k;
	}
}

const int MAXN = 8401000;

int FFT(comp *a, int n, int f) { 
	int t=1;

	// for polynomial multiplication
	while(t < n+n) t <<= 1;
	if(f == 1) memset(a+n, 0, sizeof(comp) * (t-n));
	n = t;

	butterfly(a, n);

	static comp* roots[32];
	static comp cPool[2*MAXN], *cp=cPool;

	for(int hi=1,h=2; h<=n; h=1<<(++hi)) {
		if(!roots[hi]) {
			// For n>=2^20, using the powers of the h-th primitive unity root to obtain all h-th unity roots does not give enough precision.
			roots[hi] = cp;
			for(int k=0; k<h; k++)
#ifdef USE_LONG_DOUBLE
				cp[k] = comp(cosl(-2.*M_PIl*k/h), sinl(-2.*M_PIl*k/h));
#else
				cp[k] = comp(cos(-2.*M_PIl*k/h), sin(-2.*M_PIl*k/h));
#endif
			cp += h;
		}
		comp *r = roots[hi];
		if(f==1) for(int j=0; j<n; j+=h) {
			for(int k=j; k<j+h/2; k++) {
				comp u = a[k];
				comp t = r[k-j] * a[k+h/2];
				a[k] = u + t;
				a[k+h/2] = u - t;
			}
		}
		else for(int j=0; j<n; j+=h) {
			for(int k=j; k<j+h/2; k++) {
				comp u = a[k];
				comp t = conj(r[k-j]) * a[k+h/2];
				a[k] = u + t;
				a[k+h/2] = u - t;
			}
		}
	}
	if(f == -1)
		for(int i=0; i<n; i++)
			a[i] /= n;
	return n;
}

void DAC(comp *a, const ll *d, int n) {
	for(int i=0; i<n; i++) {
#ifdef ASSERT_ON
		assert(d[i]>=0 && d[i] < P);
#endif
		ll whiten = d[i] >= P - P2 ? d[i] - P : d[i];
		ll t = whiten / M, r = whiten - M * t;
/* #define STRICT_SAFETY */
#ifdef STRICT_SAFETY
		// make t and r independent (assuming d[] random) to bound the result of DFT, slightly slower
		if(rand()&1)
			if(r < 0) {
				r += M;
				t--;
			}
			else {
				r -= M;
				t++;
			}
#endif
		a[i] = comp(r, __I * t);
	}
}

void ADC(ll *d, const comp *a, int n) {
	for(int i=0; i<n; i++) {
#ifdef USE_LONG_DOUBLE
		ll x = llroundl(real(a[i]));
		ll y = llroundl(imag(a[i]) / __I) % P;
#else
		ll x = llround(real(a[i]));
		ll y = llround(imag(a[i]) / __I) % P;
#endif
		d[i] = (y * M + x) % P;
		if(d[i] < 0) d[i] += P;
	}
}


void Mult(ll *c, ll *a, ll *b, int n) {
	static ll noise[MAXN];
	static int lastn, m;
	static comp noise_hat[MAXN];

	if(lastn != n) {
		// may use a dictionary to hold random vectors and their FFT results.
		for(int i=0; i<n; i++) noise[i] = rand() % P;
		DAC(noise_hat, noise, n);
		m = FFT(noise_hat, n, 1);
		lastn = n;
	}

	static ll a_noise[MAXN];
	static comp tmp[2][MAXN];
	for(int i=0; i<n; i++) a_noise[i] = (a[i] + noise[i]) % P;
	DAC(tmp[0], a_noise, n);
	DAC(tmp[1], b, n);
	FFT(tmp[0], n, 1);
	FFT(tmp[1], n, 1);
	for(int i=0; i<m; i++) tmp[0][i] *= tmp[1][i];
	FFT(tmp[0], n, -1);
	ADC(c, tmp[0], n*2);
	for(int i=0; i<m; i++) tmp[0][i] = tmp[1][i] * noise_hat[i];
	FFT(tmp[0], n, -1);
	ADC(a_noise, tmp[0], n*2-1);
	for(int i=0; i<n*2-1; i++) {
		c[i] -= a_noise[i];
		if(c[i] < 0) c[i] += P;
	}
}


int n = 4000000;
ll a[MAXN], b[MAXN], c[MAXN];

ll EvalPoly(const ll *a, ll x, int n) {
	ll ret = 0;
	for(int i=n; i>=0; i--) ret = (ret * x + a[i]) % P;
	return ret;
}

int main() {
	srand(time(0));
	int lastclk = 0;

	for(int ti=1; ; ti++) {

	for(int i=0; i<n; i++) a[i] = rand() % P;
	/* for(int i=0; i<n; i++) b[i] = rand() % P; */
	for(int i=0; i<n; i++) b[i] = rand()&1 ? P2 : rand() % P;
	Mult(c, a, b, n);
	printf("#%d: %d\n", ti, clock() - lastclk);
	/* for(int i=0; i<n; i++) { */
	/* 	ll t = 0; */
	/* 	for(int k=0; k<=i; k++) t = (t + a[k] * b[i-k]) % P; */
	/* 	assert(t == c[i]); */
	/* } */
	for(int k=1; k<=10; k++) {
		ll x = rand() % P;
		ll fa = EvalPoly(a, x, n-1), fb = EvalPoly(b, x, n-1), fc = EvalPoly(c, x, 2*n-2);
		assert(fa * fb % P == fc);
	}
	lastclk = clock();
	}
	return 0;
}
