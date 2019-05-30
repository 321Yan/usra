#include <Rcpp.h>
#include <iostream>
#include <vector>

using namespace std;
using namespace Rcpp;





// [[Rcpp::export]]
List genos_to_sparse_prep (CharacterMatrix genos) {
  
  int nrow = genos.nrow();
  int ncol = genos.ncol();
  
  std::vector<int> ind_i;
  std::vector<int> ind_j;
  
  for (int j = 0; j < ncol; j ++) {
    for (int i = 0; i < nrow; i++) {
      if (genos(i, j)[0] == '1') {
        ind_i.push_back(j*2+1);
        ind_j.push_back(i+1);
      }
      
      if (genos(i, j)[2] == '1') {
        ind_i.push_back(j*2+2);
        ind_j.push_back(i+1);
      }
      
    }
  }
  
  IntegerVector r_i = wrap(ind_i);
  IntegerVector r_j = wrap(ind_j);
  IntegerVector x(ind_i.size(), 1);
  
  return List::create(r_i, r_j, x);
  
  
}
