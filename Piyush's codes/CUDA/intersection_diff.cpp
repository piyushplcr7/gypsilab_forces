#include <iostream>
#include <vector>

using namespace std;
void test(vector<int> vec1, vector<int> vec2) {
  int N = vec1.size();
  vector<bool> vec1tracker(N, false);
  vector<bool> vec2tracker(N, false);
  vector<int> intersection{}, diff1{}, diff2{};
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      if (vec1[i] == vec2[j]) {
        vec1tracker[i] = true;
        vec2tracker[j] = true;
      }
    }
  }
  // constructing intersection and diffs
  for (int i = 0; i < N; ++i) {
    if (vec1tracker[i])
      intersection.push_back(vec1[i]);
    else
      diff1.push_back(vec1[i]);

    if (!vec2tracker[i])
      diff2.push_back(vec2[i]);
  }

  cout << "Intersection = " << endl;
  for (auto el : intersection)
    cout << el << ",";
  cout << endl;

  cout << "diff1 = " << endl;
  for (auto el : diff1)
    cout << el << ",";
  cout << endl;

  cout << "diff2 = " << endl;
  for (auto el : diff2)
    cout << el << ",";
  cout << endl;
}

int main() {
  vector<int> vec1 = {1, 2, 5};
  vector<int> vec2 = {2, 3, 5};
  test(vec1, vec2);
  return 0;
}