#include <stdlib.h>
#include <iostream>

using namespace std;

int dice6() {
  int s = rand() % 6;
  return s;
}

void sample() {
  for (int i = 0; i < 10; i++) {
    cout << dice6() << endl;
  }
}

int main() {
  srand(1);

  sample();
}
