#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <map>
#include <cmath>
#include <queue>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

using namespace std;

const int hmin = 2;
const int hmax = 6 + 6;

const int numDice = 100;
const int hminCoin = 0;
const int hmaxCoin = numDice;
const int numCoins = numDice;
const int Nsamples = 1000000;

int coin() {
  int s = (rand() % 100) < 50;
  return s;
}

int dice6() {
  int s = rand() % 6 + 1;
  return s;
}

void sampleDice() {
  map<char, int> histogram;

  for (int i = hmin; i <= hmax; i++) {
    histogram[i] = 0;
  }

  for (int i = 0; i < Nsamples; i++) {
    const int d1 = dice6();
    const int d2 = dice6();
    histogram[d1 + d2]++;
    //cout << dice6() << " ";
  }

  for (int i = hmin; i <= hmax; i++) {
    cout << i << ": " << histogram[i] * 100 / Nsamples << "%" << endl;
  }

  cout << endl;
}

double choose(int n, int k) {
    if (k == 0) {
      return 1;
    }
    return (n * choose(n - 1, k - 1)) / k;
}

void sampleCoins() {
  map<char, int> histogram;
  map<char, int> importanceHistogram;
  int* samples = new int[Nsamples];
  double* W = new double[Nsamples];
  double* cumulativeW = new double[Nsamples];
  double proposal = 1.0 / Nsamples;

  for (int i = hminCoin; i <= hmaxCoin; i++) {
    histogram[i] = 0;
    importanceHistogram[i] = 0;
  }

  double sum = 0.0;
  for (int i = 0; i < Nsamples; i++) {
    int k = 0;
    for (int j = 0; j < numCoins; j++) {
      k += coin();
    }
    samples[i] = k;
    W[i] = (double) k / proposal;
    histogram[k]++;
    sum += W[i];
  }

  //cout << "sum: " << sum << endl;
  for (int i = 0; i < Nsamples; i++) {
    W[i] /= sum;
    if (i > 0) {
      cumulativeW[i] = cumulativeW[i - 1] + W[i];
      if (i % 1000 == 1000 - 1) {
        //cout << "CW" << i << ": " << cumulativeW[i] << endl;
      }
    } else {
      cumulativeW[0] = W[0];
    }
  }

  for (int i = 0; i < Nsamples; i++) {
    //int r = rand() % Nsamples;
    double r = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
    int pos = upper_bound(cumulativeW, cumulativeW + Nsamples, r) - cumulativeW + 1;
    if (!pos) {
      cerr << "not found: " << r << endl;
      exit(0);
    }
    int sampleK = samples[pos];
    importanceHistogram[sampleK]++;
  }

  double td1 = 0.0;
  double td2 = 0.0;
  for (int i = hminCoin; i <= hmaxCoin; i++) {
    double nCk = choose(numCoins, i);
    double pct = nCk * pow(0.5, numCoins);
    double N = pct * Nsamples;
    cout << i << ": " << N << " ";
    cout << i << ": " << histogram[i] * 100 / Nsamples << "%: " << abs(N - histogram[i]) << " ";
    cout << i << ": " << importanceHistogram[i] * 100 / Nsamples << "%: " << abs(N - importanceHistogram[i]) << endl;
    td1 += abs(N - histogram[i]);
    td2 += abs(N - importanceHistogram[i]);
  }
  cout << "TD1: " << td1 << endl;
  cout << "TD2: " << td2 << endl;

  cout << endl;
}

double f(double x) {
  return sqrt(1 - x * x);
}

double uniform(double x) {
  return x;
}

double randX(double x) {
  double r = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
  return r;
}

double p(double x) {
  return sqrt(1 - x * x);
}

void integral() {
  double s = 0.0;
  double s2 = 0.0;
  double s3 = 0.0;
  for (int i = 0; i < Nsamples; i++) {
    double x = (double) i / Nsamples;
    double x2 = uniform(x);
    double x3 = randX(x);
    double x4 = p(x);
    s += f(x2);
    s2 += f(x3);
    s3 += f(x4);
  }
  s = 4.0 * s / Nsamples;
  s2 = 4.0 * s2 / Nsamples;
  s3 = 4.0 * s3 / Nsamples;

  double d = abs(s - M_PI);
  double d2 = abs(s2 - M_PI);
  double d3 = abs(s3 - M_PI);

  cout << "sum1: " << s << " " << setprecision(9) << d << endl;
  cout << "sum2: " << s2 << " " << setprecision(9) << d2 << endl;
  cout << "sum3: " << s3 << " " << setprecision(9) << d3 << endl;
}

double riemann(double a, double b) {
  return (f(a) + f(b)) * 0.5 * (b - a);
}

void bisectionIntegral() {
  double s = 0.0;

  double* starts = new double[Nsamples];
  double* ends = new double[Nsamples];

  auto cmp = [&](int left, int right) {
    //return (f(ends[left]) - f(starts[left])) < (f(ends[right]) - f(starts[right]));
    //double r1 = (f(ends[left]) - f(starts[left])) * (ends[left] - starts[left]);
    //double r2 = (f(ends[right]) - f(starts[right])) * (ends[right] - starts[right]);
    double r1 = riemann(ends[left], starts[left]);
    double r2 = riemann(ends[right], starts[right]);
    //return ends[left] - starts[left] < ends[right] - starts[right];
    return f(ends[left]) - f(starts[left]) > f(ends[right]) - f(starts[right]);
    //return r1 > r2;
  };
  std::priority_queue<int, std::vector<int>, decltype(cmp)> q(cmp);

  starts[0] = 0.0;
  ends[0] = 1.0;
  q.push(0);
  s += riemann(starts[0], ends[0]);

  for (int i = 0; i < Nsamples - 1; i++) {
    int index = q.top();
    q.pop();

    double start = starts[index];
    double end = ends[index];
    double mid = (start + end) * 0.5;

    if (i % 1000 == 0) {
      // cout << "m: " << mid << endl;
    }

    s -= riemann(start, end);

    starts[index] = start;
    ends[index] = mid;
    q.push(index);
    s += riemann(start, mid);

    starts[i + 1] = mid;
    ends[i + 1] = end;
    q.push(i + 1);
    s += riemann(mid, end);

    //s += riemann((double) i / Nsamples, (double) (i + 1) / Nsamples);
  }

  s = 4.0 * s; // (Nsamples - 1);

  double d = abs(s - M_PI);

  cout << "sum4: " << s << " " << setprecision(9) << d << endl;
}

int main() {
  srand(1337);

  integral();
  bisectionIntegral();

  // sampleCoins();
}
