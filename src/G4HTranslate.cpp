#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector G4HTranslate(std::string sequence) {
  int n = sequence.size();
  NumericVector score(n);
  int i = 0;

  while (i < n) {
    char base = sequence[i];
    if (base == 'G' || base == 'C') {
      int run_length = 0;
      while (i + run_length < n && sequence[i + run_length] == base) {
        run_length++;
      }
      if (base == 'G') {
        for (int j = i; j < i + run_length; ++j) {
          score[j] = std::min(run_length, 4);
        }
      } else {
        for (int j = i; j < i + run_length; ++j) {
          score[j] = -std::min(run_length, 4);
        }
      }
      i += run_length;
    } else {
      i++;
    }
  }

  return score;
}
