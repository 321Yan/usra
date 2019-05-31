// [[Rcpp::depends(RcppEigen)]]

#include <RcppEigen.h>
#include <Rcpp.h>


using namespace Rcpp;


// [[Rcpp::export]]
S4 genos_to_sparse3 (CharacterMatrix genos) {
  
  typedef Eigen::SparseMatrix<double> SpMat; 
  typedef Eigen::Triplet<double> T;
  
  int nrow = genos.nrow();
  int ncol = genos.ncol();
  std::vector<T> tripletList;
  tripletList.reserve(0.05*ncol*nrow);
  
  for (int j = 0; j < ncol; j ++) {
    for (int i = 0; i < nrow; i++) {
      if (genos(i, j)[0] == '1') {
        tripletList.push_back(T(j*2,i,1));
      }
      
      if (genos(i, j)[2] == '1') {
        tripletList.push_back(T(j*2+1,i,1));
      }
      
    }
  }
  
  SpMat mat(2*ncol,nrow);
  mat.setFromTriplets(tripletList.begin(), tripletList.end());
  
  S4 dgcmat(wrap(mat));
  
  
  return dgcmat;
  
  
}
