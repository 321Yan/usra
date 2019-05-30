#include <Rcpp.h>
#include <iostream>
#include <vector>

using namespace std;
using namespace Rcpp;



// this function does not check for missing values

// [[Rcpp::export]]
S4 genos_to_sparse (CharacterMatrix genos) {
  
  int nrow = genos.nrow();
  int ncol = genos.ncol();
  
  std::vector<int> ind_i;
  std::vector<int> int_p(1, 0);
  
  for (int i = 0; i < nrow; i++) {
    int sum = 0;
    for (int j = 0; j < ncol; j++) {
      if (genos(i, j)[0] == '1') {
        ind_i.push_back(j*2);
		sum++;
      }
      
      if (genos(i, j)[2] == '1') {
        ind_i.push_back(j*2+1);
		sum++;
      }
    }
    int_p.push_back(int_p.back() + sum);
    
  }
  
  NumericVector x(ind_i.size(), 1);
  
  S4 smat("dgCMatrix");
  smat.slot("i") = ind_i;
  smat.slot("p") = int_p;
  smat.slot("x") = x;
  smat.slot("Dim") = IntegerVector::create(ncol*2, nrow);
  
  
  // return List::create(ind_i, int_p, x);
  return smat;
  
}