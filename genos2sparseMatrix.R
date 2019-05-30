get_SMindex_by_person <- function(person_index, person_genos){
  #j_1 stores the positions where the
  #individual with person_index carries an SNV
  #on their 1st haplotype
  #
  #These will become the non-zero column entries
  #in our sparseMatrix (i.e. argument j to sparseMatrix)
  j_1 <- which(person_genos %in% c("1|0", "1|1"))

  #j_2 stores the positions where the
  #individual with person_index carries an SNV
  #on their 2nd haplotype
  #
  #These will also become non-zero column entries
  #in our sparseMatrix (i.e. argument j to sparseMatrix)
  j_2 <- which(person_genos %in% c("0|1", "1|1"))


  # create the vector of row positions for person with person_index,
  # which will be supplied to sparseMatrix
  i_pos <- c(rep((2*person_index - 1), length(j_1)),
             rep((2*person_index), length(j_2)))

  return(list(i_pos = i_pos,
              j_pos = c(j_1, j_2)))

}


genos2sparseMatrix <- function(genos){
  # Get the row and column location of each mutation
  # person-by-person (i.e. row-by-row) from the genotypes
  # matrix returned by read.vcfR
  #
  # NOTE: index_by_person, below, is a list of lists
  # the first item in the first list are the row postions of SNVs for the first person
  # and the second item in the first list are the column postions of SNVs for the first person
  index_by_person <- lapply(1:ncol(genos), function(x){
    get_SMindex_by_person(x, genos[, x])
  })

  #input the row and column data into the sparseMatrix
  SM_format <- sparseMatrix(i = unlist(lapply(index_by_person, `[[`, 1)),
                            j = unlist(lapply(index_by_person, `[[`, 2)),
                            x = rep(1, length(unlist(lapply(index_by_person, `[[`, 1)))))

  return(SM_format)
}


#------------------------#
#      Examples          #
#------------------------#
library(vcfR)
library(Matrix)
setwd("D:/sfu/usra")
vcf <- read.vcfR("chr22exons.vcf.gz")
# vcf <- read.vcfR("D:/sfu/usra/LDheatmap/vignettes/snp_in_vcf.vcf")
# vcf <- read.vcfR(system.file("vcf/CEU.exon.2010_09.genotypes.vcf.gz", package="GGtools"))

#NOTE: removing the first column since it does not hold genotype data
my_genos <- vcf@gt[, -1]

#genotypes for the first 30 mutations and first 5 people
my_genos[1:30, 1:5]

#run convert to sparseMatrix
haps <- genos2sparseMatrix(my_genos)
myhaps <- GT_to_numeric(my_genos)
#view the first 10 mutations for the first 5 people
#REMEMBER in sparseMatrix format people are rows and mutations are columns
haps[1:10, 1:10]


library(rbenchmark)
benchmark(
  genos2sparseMatrix(my_genos),
  GT_to_numeric(my_genos),
  replications = 2
)

# with dimnames
#                           test replications elapsed relative user.self sys.self user.child sys.child
# 1 genos2sparseMatrix(my_genos)           20    47.8    1.000     46.34     1.39         NA        NA
# 2      GT_to_numeric(my_genos)           20   157.2    3.289    150.60     6.22         NA        NA

# without dimnames
# same
# space used also same

mat11 = matrix(rep("1|1",35146*2548),nrow=35146)

mat11s = mat11[1:6000,1:400]

mat00 = matrix(rep("0|0",35146*2548),nrow=35146)



#                        test replications elapsed relative user.self sys.self user.child sys.child
# 1 genos2sparseMatrix(mat11)            2   82.95    6.578     34.94    47.86         NA        NA
# 2      GT_to_numeric(mat11)            2   12.61    1.000     12.23     0.30         NA        NA



benchmark(
  lapply(1:ncol(my_genos), function(x){
    get_SMindex_by_person(x, my_genos[, x])
  }),
  sparseMatrix(i = unlist(lapply(index_by_person, `[[`, 1)),
               j = unlist(lapply(index_by_person, `[[`, 2)),
               x = rep(1, length(unlist(lapply(index_by_person, `[[`, 1))))),
  replications = 20
)

#  replications elapsed relative user.self sys.self user.child sys.child
# 1           20   38.86   11.103     38.50     0.20         NA        NA
# 2           20    3.50    1.000      2.42     1.08         NA        NA




index11 = lapply(1:ncol(mat11), function(x){
  get_SMindex_by_person(x, mat11[, x])
})


benchmark(
  lapply(1:ncol(mat11), function(x){
    get_SMindex_by_person(x, mat11[, x])
  }),
  sparseMatrix(i = unlist(lapply(index11, `[[`, 1)),
               j = unlist(lapply(index11, `[[`, 2)),
               x = rep(1, length(unlist(lapply(index11, `[[`, 1))))),
  replications = 2
)

# mat11s replication 20, elapsed, relative
# 2.24    1.000
# 6.08    2.714

# mat11 replication 2, elapsed, relative
# 33.81    1.000
# 43.33    1.282




## Now in ./Tsparse.R
## setAs("dgTMatrix", "dgCMatrix",
##       function(from) .Call(Tsparse_to_Csparse, from, FALSE)
##       )

library(Rcpp)
sourceCpp("D:/sfu/usra/sparse_Matrix_prep.cpp")
sourceCpp("D:/sfu/usra/sparse_Matrix.cpp")

index_by_person <- lapply(1:ncol(my_genos), function(x){
  get_SMindex_by_person(x, my_genos[, x])
})

ori_mat <- sparseMatrix(i = unlist(lapply(index_by_person, `[[`, 1)),
                          j = unlist(lapply(index_by_person, `[[`, 2)),
                          x = rep(1, length(unlist(lapply(index_by_person, `[[`, 1)))))

mylist = genos_to_sparse_prep(my_genos)


myfunction = function(my_genos){
  List = genos_to_sparse_prep(my_genos)
  sparseMatrix(i = List[[1]], j = List[[2]], x = List[[3]])
}

benchmark(
  genos2sparseMatrix(my_genos),
  myfunction(my_genos),
  replications = 20
)

#                           test replications elapsed relative user.self sys.self user.child sys.child
# 1 genos2sparseMatrix(my_genos)           20   68.92    4.537     67.29     1.50         NA        NA
# 2         myfunction(my_genos)           20   15.19    1.000     14.19     0.99         NA        NA



################## more optimazatioin #######################
# extrect essentials from sparseMatrix

benchmark(
  GT_to_numeric(mat11),
  myfunction(mat11),
  genos2sparseMatrix(mat11),
  replications = 2
)


#                        test replications elapsed relative user.self sys.self user.child sys.child
# 3 genos2sparseMatrix(mat11)            2   62.95    5.512     27.82    33.94         NA        NA
# 1      GT_to_numeric(mat11)            2   11.42    1.000     11.17     0.18         NA        NA
# 2         myfunction(mat11)            2   38.72    3.391     19.03    19.32         NA        NA


myfunction_missing = function(my_genos){
  List = GT_to_sparse_prep(my_genos)
  sparseMatrix(i = List[[1]], j = List[[2]], x = List[[3]])
}


benchmark(
  myfunction(my_genos),
  myfunction_missing(my_genos),
  replications = 20
  
)

smat1 = myfunction(my_genos)
smat2 = myfunction_missing(my_genos)
#                           test replications elapsed relative user.self sys.self user.child sys.child
# 1         myfunction(my_genos)           20   15.58    1.000     13.96     1.63         NA        NA
# 2 myfunction_missing(my_genos)           20  132.65    8.514    131.47     0.84         NA        NA


unique(as.vector(my_genos))




# ---------------------------------------------
directCpp = genos_to_sparse(my_genos)
smat1 = myfunction(my_genos)
identical(directCpp, smat1)

benchmark(
  genos_to_sparse(my_genos),
  genos_to_sparse_prep(my_genos),
  sparseMatrix(i = List[[1]], j = List[[2]], x = List[[3]]),
  replications = 20
  
)
# s4 here, List is produced by genos_to_sparse_prep(my_genos)
#                                                        test replications elapsed relative
# 1                                 genos_to_sparse(my_genos)           20   26.84   10.284
# 2                            genos_to_sparse_prep(my_genos)           20   13.99    5.360
# 3 sparseMatrix(i = List[[1]], j = List[[2]], x = List[[3]])           20    2.61    1.000


benchmark(
  genos_to_sparse(mat11),
  genos_to_sparse_prep(mat11),
  replications = 20
  
)
# !!! here genos_to_sparse returns a list of elements that will be filled in to the slots
#                          test replications elapsed relative user.self sys.self user.child
# 1      genos_to_sparse(mat11)           20   63.08    1.000     41.58    21.48         NA
# 2 genos_to_sparse_prep(mat11)           20   69.24    1.098     37.39    31.72         NA
#   sys.child
# 1        NA
# 2        NA



benchmark(
  genos_to_sparse(mat11),
  myfunction(mat11),
  replications = 5
  
)
# !!! here genos_to_sparse returns a dgcMatrix
#                     test replications elapsed relative user.self sys.self user.child sys.child
# 1 genos_to_sparse(mat11)            5   14.84    1.000     10.34      4.5         NA        NA
# 2      myfunction(mat11)            5   87.05    5.866     56.32     30.0         NA        NA



