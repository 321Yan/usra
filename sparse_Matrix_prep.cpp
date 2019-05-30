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


// problematic
// List sparse_Matrix(NumericVector i, NumericVector j, NumericVector x) {
//   
//   if (Rf_isNull(i) || Rf_isNull(j) || Rf_isNull(x)) stop("All of i, j, x must be provided");
//   Rcout << "0" << endl;
//   
//   int n = i.size();
//   IntegerVector dims((int)(max(i)), (int)(max(j)));
//   NumericVector new_i(n);
//   NumericVector new_j(n);
//   
//   Rcout << "1" << endl;
//   
//   for (int k = 0; k < n; k++) {
//     new_i[k] = i[k] - 1;
//     new_j[k] = j[k] - 1;
//     Rcout << "new_i is "<< new_i[k] << endl;
//     Rcout << "new_j is "<< new_j[k] << endl;
//     Rcout << "k is " << k << endl;
//     
//     
//     
//   }
//   
//   Rcout << "2" << endl;
//   
//   S4 r("dgTMatrix");
//   r.slot("Dim") = dims;
//   r.slot("x") = x;
//   r.slot("i") = new_i;
//   r.slot("j") = new_j;
//   
//   Rcout << "3" << endl;
//   
//   return r; 
// }
//   
  
