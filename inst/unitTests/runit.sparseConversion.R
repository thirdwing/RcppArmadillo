#!/usr/bin/r -t
#
# Copyright (C) 2014  Dirk Eddelbuettel
#
# This file is part of RcppArmadillo. It is based on the documentation of package Matrix, 
# slam, SparseM, spam and SciPy, which are respectively created by developers: of the packages:
# Douglas Bates, Martin Maechler; Kurt Hornik, David Meyer, Christian Buchta; 
# Roger Koenker, Pin Ng; Reinhard Furrer, Florian Gerber, Daniel Gerber, 
# Kaspar Moesinger, Youcef Saad, Esmond G. Ng, Barry W. Peyton, Joseph W.H. Liu,
# Alan D. George; the developers of SciPy. It is also modified by Binxiang Ni
#
# RcppArmadillo is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# RcppArmadillo is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with RcppArmadillo.  If not, see <http://www.gnu.org/licenses/>.

.runThisTest <- suppressMessages(require(Matrix))

if (.runThisTest) {
  
    .setUp <- RcppArmadillo:::unit_test_setup("sparse.cpp") 
  
    ## setting up an example matrix -- using the fact that the as<sp_mat>
    ## converter prefers sparse matrix objects create by the Matrix package
    suppressMessages(require(Matrix))

    # Matrix
    # https://cran.r-project.org/web/packages/Matrix/Matrix.pdf
    
    ## p10 (dgCMatrix)
    set.seed(7)
    m <- matrix(0, 5, 5)
    m[sample(length(m), size = 14)] <- rep(1:9, length=14)
    mm <- as(m, "CsparseMatrix")
    checkEquals(mm, asSpMat(mm), msg="dgC2dgC")
    
    ## p35 (ddiMatrix)
    d2 <- Diagonal(x = c(10, 1))
    mtxt <- c("10 0",
              "0 1")
    M <- as.matrix(read.table(textConnection(mtxt)))
    dimnames(M) <- NULL
    dgc <- as(M, "dgCMatrix")
    checkEquals(dgc, asSpMat(d2), msg="ddi2dgC")
    
    cd2 <- chol(d2)
    dgc <- as(chol(dgc), "dgCMatrix")
    checkEquals(dgc, asSpMat(cd2), msg="ddi2dgC")
    
    ## p36 (dgCMatrix)
    m <- Matrix(c(0,0,2:0), 3,5)
    checkEquals(m, asSpMat(m), msg="dgC2dgC")
    
    ## p40 (dgTMatrix)
    m <- Matrix(0+1:28, nrow = 4)
    m[-3,c(2,4:5,7)] <- m[ 3, 1:4] <- m[1:3, 6] <- 0
    mT <- as(m, "dgTMatrix")
    dgc <- as(m, "dgCMatrix")
    checkEquals(dgc, asSpMat(mT), msg="dgT2dgC")
    
    ## p40 (dgTMatrix)
    T2 <- new("dgTMatrix",
              i = as.integer(c(1,1,0,3,3)),
              j = as.integer(c(2,2,4,0,0)), x=10*1:5, Dim=4:5)
    dgc <- as(T2, "dgCMatrix")
    checkEquals(dgc, asSpMat(T2), msg="dgT2dgC")
    
    # p42 (ddiMatrix)
    mtxt <- c("1 0 0",
              "0 1 0",
              "0 0 1")
    M <- as.matrix(read.table(textConnection(mtxt)))
    dimnames(M) <- NULL
    dgc <- as(M, "dgCMatrix")
    checkEquals(dgc, asSpMat(Diagonal(3)), msg="ddi2dgC")
    
    mtxt <- c("1000 0   0",
              "0    100 0",
              "0    0   10")
    M <- as.matrix(read.table(textConnection(mtxt)))
    dimnames(M) <- NULL
    dgc <- as(M, "dgCMatrix")
    checkEquals(dgc, asSpMat(Diagonal(x = 10^(3:1))), msg="ddi2dgC")
    
    ## p42 (dgTMatrix)
    M1 <- Matrix(0+0:5, 2,3)
    M <- kronecker(Diagonal(3), M1)
    dgc <- as(M, "dgCMatrix")
    checkEquals(dgc, asSpMat(M), msg="dgT2dgC")
    
    ## p42 (ldiMatrix) (To be continued)
    ## Diagonal(x = (1:4) >= 2)
    
    ## p42 (dsCMatrix)
    S <- crossprod(Matrix(rbinom(60, size=1, prob=0.1), 10,6))
    dgc <- as(S, "dgCMatrix")
    checkEquals(dgc, asSpMat(S), msg="dsC2dgC")
    
    SI <- S + 10*.symDiagonal(6)
    dgc <- as(SI, "dgCMatrix")
    checkEquals(dgc, asSpMat(SI), msg="dgT2dgC")
    
    ## p42 (dtCMatrix)
    I4 <- .sparseDiagonal(4, shape="t")
    dgc <- as(I4, "dgCMatrix")
    checkEquals(dgc, asSpMat(I4), msg="dtC2dgC")
    
    ## p43 (ddiMatrix)
    I5 <- Diagonal(5)
    mtxt <- c("1 0 0 0 0",
              "0 1 0 0 0",
              "0 0 1 0 0",
              "0 0 0 1 0",
              "0 0 0 0 1")
    M <- as.matrix(read.table(textConnection(mtxt)))
    dimnames(M) <- NULL
    dgc <- as(M, "dgCMatrix")
    checkEquals(dgc, asSpMat(I5), msg="ddi2dgC")
    
    D5 <- Diagonal(x = 10*(1:5))
    mtxt <- c("10 0 0 0 0",
              "0 20 0 0 0",
              "0 0 30 0 0",
              "0 0 0 40 0",
              "0 0 0 0 50")
    M <- as.matrix(read.table(textConnection(mtxt)))
    dimnames(M) <- NULL
    dgc <- as(M, "dgCMatrix")
    checkEquals(dgc, asSpMat(D5), msg="ddi2dgC")
    
    # mtxt <- c("0.1 0    0          0     0",
    #           "0   0.05 0          0     0",
    #           "0   0    0.03333333 0     0",
    #           "0   0    0          0.025 0",
    #           "0   0    0          0     0.02")
    # M <- as.matrix(read.table(textConnection(mtxt)))
    # dimnames(M) <- NULL
    # dgc <- as(M, "dgCMatrix")
    # checkEquals(dgc, asSpMat(solve(D5)), msg="ddi2dgC")
    
    ## p49 (dgTMatrix)
    m <- spMatrix(10,20, i= 1:8, j=2:9, x = c(0:2,3:-1))
    dgc <- as(m, "dgCMatrix")
    checkEquals(dgc, asSpMat(m), msg="dgT2dgC")
    checkEquals(drop0(m), asSpMat(drop0(m)), msg="dgT2dgC")
    
    ## p49 (dtCMatrix)
    t5 <- new("dtCMatrix", Dim = c(5L, 5L), uplo = "L",
              x = c(10, 1, 3, 10, 1, 10, 1, 10, 10),
              i = c(0L,2L,4L, 1L, 3L,2L,4L, 3L, 4L),
              p = c(0L, 3L, 5L, 7:9))
    dgc <- as(t5, "dgCMatrix")
    checkEquals(dgc, asSpMat(t5), msg="dtC2dgC")
    
    ## p50 (dsCMatrix)
    mm <- Matrix(toeplitz(c(10, 0, 1, 0, 3)), sparse = TRUE)
    dgc <- as(mm, "dgCMatrix")
    checkEquals(dgc, asSpMat(mm), msg="dsC2dgC")
    
    ## p51 (dgTMatrix)
    mT <- as(mm, "dgTMatrix")
    checkEquals(dgc, asSpMat(mT), msg="dgT2dgC")
    
    ## p51 (dsTMatrix)
    symM <- as(mT, "symmetricMatrix")
    checkEquals(dgc, asSpMat(symM), msg="dsT2dgC")
    
    ## p51 (dsCMatrix)
    symC <- as(symM, "CsparseMatrix")
    checkEquals(dgc, asSpMat(symC), msg="dsC2dgC")
    
    sC <- Matrix(mT, sparse=TRUE, forceCheck=TRUE)
    checkEquals(dgc, asSpMat(sC), msg="dsC2dgC")
    
    ## p51 (dsTMatrix)
    sym2 <- as(symC, "TsparseMatrix")
    checkEquals(dgc, asSpMat(sym2), msg="dsT2dgC")
    
    ## p53 (dsRMatrix)
    m2 <- new("dsRMatrix", Dim = c(2L,2L),
              x = c(3,1), j = c(1L,1L), p = 0:2)
    mtxt <- c("0 3",
              "3 1")
    M <- as.matrix(read.table(textConnection(mtxt)))
    dimnames(M) <- NULL
    dgc <- as(M, "dgCMatrix")
    checkEquals(dgc, asSpMat(m2), msg="dsR2dgC")
    
    ds2 <- forceSymmetric(diag(2))
    dR <- as(ds2, "RsparseMatrix")
    mtxt <- c("1 0",
              "0 1")
    M <- as.matrix(read.table(textConnection(mtxt)))
    dimnames(M) <- NULL
    dgc <- as(M, "dgCMatrix")
    checkEquals(dgc, asSpMat(dR), msg="dsR2dgC")
    
    ## p56 (dtTMatrix)
    t1 <- new("dtTMatrix", x= c(3,7), i= 0:1, j=3:2, Dim= as.integer(c(4,4)))
    mtxt <- c("0 0 0 3",
             "0 0 7 0",
             "0 0 0 0",
             "0 0 0 0")
    M <- as.matrix(read.table(textConnection(mtxt)))
    dimnames(M) <- NULL
    dgc <- as(M, "dgCMatrix")
    checkEquals(dgc, asSpMat(t1), msg="dtT2dgC")
    
    tu <- t1 ; tu@diag <- "U"
    mtxt <- c("1 0 0 3",
              "0 1 7 0",
              "0 0 1 0",
              "0 0 0 1")
    M <- as.matrix(read.table(textConnection(mtxt)))
    dimnames(M) <- NULL
    dgc <- as(M, "dgCMatrix")
    checkEquals(dgc, asSpMat(tu), msg="dtT2dgC")
    
    ## p56 (dtCMatrix)
    cu <- as(tu, "dtCMatrix")
    dgc <- as(cu, "dgCMatrix")
    checkEquals(dgc, asSpMat(cu), msg="dtC2dgC")
    
    ## p56 (dtTMatrix)
    t1[1,2:3] <- -1:-2
    mtxt <- c("0 -1 -2 3",
              "0  0  7 0",
              "0  0  0 0",
              "0  0  0 0")
    M <- as.matrix(read.table(textConnection(mtxt)))
    dimnames(M) <- NULL
    dgc <- as(M, "dgCMatrix")
    checkEquals(dgc, asSpMat(t1), msg="dtT2dgC")
    
    diag(t1) <- 10*c(1:2,3:2)
    mtxt <- c("10 -1 -2 3",
              "0  20  7 0",
              "0  0  30 0",
              "0  0  0 20")
    M <- as.matrix(read.table(textConnection(mtxt)))
    dimnames(M) <- NULL
    dgc <- as(M, "dgCMatrix")
    checkEquals(dgc, asSpMat(t1), msg="dtT2dgC")
    
    ## p56 (dtCMatrix)
    U5 <- new("dtCMatrix", i= c(1L, 0:3), p=c(0L,0L,0:2, 5L), Dim = c(5L, 5L),
              x = rep(1, 5), diag = "U")
    dgc <- as(U5, "dgCMatrix")
    checkEquals(dgc, asSpMat(U5), msg="dtC2dgC")
    
    ## p59 (dtRMatrix)
    m2 <- new("dtRMatrix", Dim = c(2L,2L),
               x = c(5, 1:2), p = c(0L,2:3), j= c(0:1,1L))
    mtxt <- c("5 1",
              "0 2")
    M <- as.matrix(read.table(textConnection(mtxt)))
    dimnames(M) <- NULL
    dgc <- as(M, "dgCMatrix")
    checkEquals(dgc, asSpMat(m2), msg="dtR2dgC")
    
    m3 <- as(Diagonal(2), "RsparseMatrix")
    mtxt <- c("1 0",
              "0 1")
    M <- as.matrix(read.table(textConnection(mtxt)))
    dimnames(M) <- NULL
    dgc <- as(M, "dgCMatrix")
    checkEquals(dgc, asSpMat(m3), msg="dtR2dgC")
    
    ## p74 (pMatrix)
    p1 <- as(c(2,3,1), "pMatrix")
    mtxt <- c("0 1 0",
              "0 0 1",
              "1 0 0")
    M <- as.matrix(read.table(textConnection(mtxt)))
    dimnames(M) <- NULL
    dgc <- as(M, "dgCMatrix")
    checkEquals(dgc, asSpMat(p1), msg="p2dgC")
    
    ## p74 (indMatrix) (To be continued)
    sm1 <- as(rep(c(2,3,1), e=3), "indMatrix")
    set.seed(27)
    s10 <- as(sample(10, 30, replace=TRUE),"indMatrix")
    set.seed(27)
    IM1 <- as(sample(1:20, 100, replace=TRUE), "indMatrix")
    IM2 <- as(sample(1:18, 100, replace=TRUE), "indMatrix")
    
    ## p74 (dgCMatrix)
    c12 <- crossprod(IM1,IM2)
    checkEquals(c12, asSpMat(c12), msg="dgC2dgC")
    
    ## p74 (indMatrix) (To be continued)
    # as(2:4, "indMatrix")
    # as(list(2:4, 5), "indMatrix")
    
    ### p74 (ngTMatrix) (To be continued)
    # ngt <- as(sm1, "ngTMatrix")
    # mtxt <- c("0 1 0",
    #           "0 1 0",
    #           "0 1 0",
    #           "0 0 1",
    #           "0 0 1",
    #           "0 0 1",
    #           "1 0 0",
    #           "1 0 0",
    #           "1 0 0")
    # M <- as.matrix(read.table(textConnection(mtxt)))
    # dimnames(M) <- NULL
    # dgc <- as(M, "dgCMatrix")
    # checkEquals(dgc, asSpMat(ngt), msg="ngT2dgC")
    # 
    # ngt <- s10[1:7, 1:4] 
    # mtxt <- c("0 0 0 0",
    #           "1 0 0 0",
    #           "0 0 0 0",
    #           "0 0 0 1",
    #           "0 0 1 0",
    #           "0 0 0 0",
    #           "1 0 0 0")
    # M <- as.matrix(read.table(textConnection(mtxt)))
    # dimnames(M) <- NULL
    # dgc <- as(M, "dgCMatrix")
    # checkEquals(dgc, asSpMat(ngt), msg="ngT2dgC")
    
    ## p74 (indMatrix) (To be continued)
    # ind <- s10[1:4, ]
    # dgc <- as(as(ind, "matrix"), "dgCMatrix")
    # checkEquals(dgc, asSpMat(ind), msg="ind2dgC")
    # 
    # I1 <- as(c(5:1,6:4,7:3), "indMatrix")
    # dgc <- as(as(I1, "matrix"), "dgCMatrix")
    # checkEquals(dgc, asSpMat(I1), msg="ind2dgC")
    
    ## p74 (pMatrix)
    I2 <- as(7:1, "pMatrix")
    dgc <- as(as(I2, "matrix"), "dgCMatrix")
    checkEquals(dgc, asSpMat(I2), msg="p2dgC")
    
    ## p77 (dgTMatrix)
    A <- spMatrix(10,20, i = c(1,3:8),
                  j = c(2,9,6:10),
                  x = 7 * (1:7))
    dgc <- as(A, "dgCMatrix")
    checkEquals(dgc, asSpMat(A), msg="dgT2dgC")
    
    ## p85 (ldiMatrix) (To be continued)
    lM <- Diagonal(x = c(TRUE,FALSE,FALSE))
    
    ## p85 (ddiMatrix)
    # ddi <- crossprod(lM) 
    # mtxt <- c("1 0 0",
    #           "0 0 0",
    #           "0 0 0")
    # M <- as.matrix(read.table(textConnection(mtxt)))
    # dimnames(M) <- NULL
    # dgc <- as(M, "dgCMatrix")
    # checkEquals(dgc, asSpMat(ddi), msg="ddi2dgC")
    
    # ## p85 (ntTMatrix) (To be continued)
    # nM <- as(lM, "nMatrix")
    # checkEquals(dgc, asSpMat(nM), msg="ntT2dgC")
    
    # ## p85 (nsTMatrix) (To be continued)
    # nsc <- crossprod(nM)
    # checkEquals(dgc, asSpMat(nsc), msg="nsC2dgC")
    
    ## p87 (dgCMatrix)
    m <- Matrix(c(0,0,2:0), 3,5, dimnames=list(LETTERS[1:3],NULL))
    dimnames(m) <- NULL
    checkEquals(m, asSpMat(m), msg="dgC2dgC")
    
    # ## p87 (lgCMatrix) (To be continued)
    # lm <- (m > 1)
    
    ## p111 (ngCMatrix) (To be continued)
    # m <- Matrix(c(0,0,2:0), 3,5, dimnames=list(LETTERS[1:3],NULL))
    # dimnames(m) <- NULL
    # nm <- as(m, "nsparseMatrix")
    
    # ## p111 (lgCMatrix) (To be continued)
    # nnm <- !nm     # no longer sparse
    # nnm <- as(nnm, "sparseMatrix")
    
    ## p116 (pMatrix)
    pm1 <- as(as.integer(c(2,3,1)), "pMatrix")
    dgc <- as(as(pm1, "matrix"), "dgCMatrix")
    checkEquals(dgc, asSpMat(pm1), msg="p2dgC")
    
    pm2 <- t(pm1)
    dgc <- as(as(pm2, "matrix"), "dgCMatrix")
    checkEquals(dgc, asSpMat(pm2), msg="p2dgC")
    
    ## p116 (pMatrix)
    set.seed(11)
    ## random permutation matrix :
    p10 <- as(sample(10),"pMatrix")
    dgc <- as(as(p10, "matrix"), "dgCMatrix")
    checkEquals(dgc, asSpMat(p10), msg="p2dgC")
    
    # ## p116 (ngTMatrix) (To be continued)
    # as(pm1, "ngTMatrix")
    # p10[1:7, 1:4]
    
    ## p116 (indMatrix)
    # ind <- p10[1:3, ]
    # dgc <- as(as(ind, "matrix"), "dgCMatrix")
    # checkEquals(dgc, asSpMat(ind), msg="ind2dgC")
    
    ## p118 (dgCMatrix)
    f1 <- gl(5, 3, labels = LETTERS[1:5])
    X <- as(f1, "sparseMatrix")
    dimnames(X) <- NULL
    checkEquals(X, asSpMat(X), msg="dgC2dgC")
    
    ## p129 (dsCMatrix)
    S9 <- rsparsematrix(9, 9, nnz = 10, symmetric=TRUE)
    dgc <- as(S9, "dgCMatrix")
    checkEquals(dgc, asSpMat(S9), msg="dsC2dgC")
    
    # ## p129 (ngCMatrix) (To be continued)
    # (n7 <- rsparsematrix(5, 12, nnz = 10, rand.x = NULL))
    
    ## p129 (dgTMatrix)
    set.seed(129)
    T2 <- rsparsematrix(40, 12, nnz = 99, giveCsparse=FALSE)
    dgc <- as(T2, "dgCMatrix")
    checkEquals(dgc, asSpMat(T2), msg="dgT2dgC")
    
    ## p142 (dgCMatrix)
    i <- c(1,3:8); j <- c(2,9,6:10); x <- 7 * (1:7)
    A <- sparseMatrix(i, j, x = x)
    checkEquals(A, asSpMat(A), msg="dgC2dgC")
    
    ## p142 (dsCMatrix)
    sA <- sparseMatrix(i, j, x = x, symmetric = TRUE)
    dgc <- as(sA, "dgCMatrix")
    checkEquals(dgc, asSpMat(sA), msg="dsC2dgC")
    
    ## p142 (dtCMatrix)
    tA <- sparseMatrix(i, j, x = x, triangular= TRUE)
    dgc <- as(tA, "dgCMatrix")
    checkEquals(dgc, asSpMat(tA), msg="dtC2dgC")
    
    ## p152 (dgTMatrix)
    A <- spMatrix(10,20, i = c(1,3:8),
                  j = c(2,9,6:10),
                  x = 7 * (1:7))
    dgc <- as(A, "dgCMatrix")
    checkEquals(dgc, asSpMat(A), msg="dgT2dgC")
    
    ## p152 (lgTMatrix) (To be continued)
    # L <- spMatrix(9, 30, i = rep(1:9, 3), 1:27,
    #               (1:27) %% 4 != 1)
    
    # slam
    # https://cran.r-project.org/web/packages/slam/slam.pdf
    
    ## p4 (dgCMatrix)
    x <- matrix(c(1, 0, 0, 2, 1, 0), nrow = 3)
    SM <- Matrix(x, sparse = TRUE)
    checkEquals(SM, asSpMat(SM), msg="dgC2dgC")
    
    ## p9 (dgCMatrix)
    x <- matrix(c(1, 0, 0, 2, 1, NA), nrow = 2)
    SM <- Matrix(x, sparse = TRUE)
    checkEquals(SM, asSpMat(SM), msg="dgC2dgC")
    
    ## p12 (dgCMatrix)
    x <- matrix(c(1, 0, 0, 2), nrow = 2)
    SM <- Matrix(x, sparse = TRUE)
    dgc <- as(SM, "dgCMatrix")
    checkEquals(dgc, asSpMat(SM), msg="dgC2dgC")
    
    # SparseM
    # https://cran.r-project.org/web/packages/SparseM/SparseM.pdf
    
    ## p14
    
    ## p21 (dgTMatrix)
    set.seed(21)
    a <- rnorm(20*5)
    A <- matrix(a,20,5)
    A[row(A)>col(A)+4|row(A)<col(A)+3] <- 0
    TA <- as(A, "dgTMatrix")
    CA <- as(A, "dgCMatrix")
    checkEquals(CA, asSpMat(TA), msg="dgT2dgC")
    
    set.seed(22)
    b <- rnorm(20*5)
    B <- matrix(b,20,5)
    B[row(A)>col(A)+2|row(A)<col(A)+2] <- 0
    TB <- as(B, "dgTMatrix")
    CB <- as(B, "dgCMatrix")
    checkEquals(CB, asSpMat(TB), msg="dgT2dgC")
    
    # (dgCMatrix)
    cp <- crossprod(CA, CB)
    checkEquals(cp, asSpMat(cp), msg="dgC2dgC")
    
    ## p23 (dgRMatrix)
    n1 <- 10
    p <- 5
    a <- rnorm(n1*p)
    a[abs(a)<0.5] <- 0
    A <- matrix(a,n1,p)
    RA <- as(A, "dgRMatrix")
    CA <- as(A, "dgCMatrix")
    checkEquals(CA, asSpMat(RA), msg="dgR2dgC")
    
    ## p25 (dgRMatrix)
    n1 <- 10
    n2 <- 10
    p <- 6
    set.seed(16)
    a <- rnorm(n1*p)
    a[abs(a) < 0.5] <- 0
    A <- matrix(a,n1,p)
    RA <- as(A, "dgRMatrix")
    CA <- as(A, "dgCMatrix")
    checkEquals(CA, asSpMat(RA), msg="dgR2dgC")
    
    set.seed(25)
    b <- rnorm(n2*p)
    b[abs(b)<1.0] <- 0
    B <- matrix(b,n2,p)
    RB <- as(B, "dgRMatrix")
    CB <- as(B, "dgCMatrix")
    checkEquals(CB, asSpMat(RB), msg="dgR2dgC")
    
    ## (dgTMatrix)
    dgt <- RA %x% matrix(1:4,2,2)
    dgc <- as(dgt, "dgCMatrix")
    checkEquals(dgc, asSpMat(dgt), msg="dgT2dgC")
    
    # spam
    # https://cran.r-project.org/web/packages/spam/spam.pdf
    
    ## p4 (dgTMatrix)
    x <- matrix(c(1, 0, 0, 2, 1, 0), nrow = 3)
    dgt <- as(x, "dgTMatrix")
    dgc <- as(x, "dgCMatrix")
    checkEquals(dgc, asSpMat(dgt), msg="dgT2dgC")
    
    ## p12 (dgTMatrix)
    x <- matrix(c(1, 0, 0, 2), nrow = 2)
    dgt <- as(x, "dgTMatrix")
    dgc <- as(x, "dgCMatrix")
    checkEquals(dgc, asSpMat(dgt), msg="dgT2dgC")
    
    ## p14 (dgTMatrix)
    x <- matrix(c(1, 0, 0, 2, 1, NA), nrow = 3)
    dgt <- as(x, "dgTMatrix")
    dgc <- as(x, "dgCMatrix")
    checkEquals(dgc, asSpMat(dgt), msg="dgT2dgC")
    
    
    # SciPy
    # https://docs.scipy.org/doc/scipy/reference/sparse.html
    
    ## (dgCMatrix)
    mtxt <- c("1 2 0",
              "0 0 3",
              "4 0 5")
    M <- as.matrix(read.table(textConnection(mtxt)))
    dimnames(M) <- NULL
    A <- Matrix(M, sparse = TRUE)
    checkEquals(A, asSpMat(A), msg="dgC2dgC")
    
    V <- c(1, 0, -1)
    checkEquals((A * V), asSpMat(A * V), msg="dgC2dgC")
    
    ## (dgTMatrix)
    dgt <- new("dgTMatrix", x = c(4, 5, 7, 9), i = as.integer(c(0, 3, 1, 0)), j = as.integer(c(0, 3, 1, 2)), Dim= as.integer(c(4,4)))
    dgc <- as(dgt, "dgCMatrix")
    checkEquals(dgc, asSpMat(dgt), msg="dgT2dgC")
    
    ## (dgRMatrix)
    dgt <- new("dgTMatrix", x = rep(1, 7), i = as.integer(c(0, 0, 1, 3, 1, 0, 0)), j = as.integer(c(0, 2, 1, 3, 1, 0, 0)), Dim= as.integer(c(4, 4)))
    dgr <- as(as(dgt, "matrix"), "dgRMatrix")
    dgc <- as(dgt, "dgCMatrix")
    checkEquals(dgc, asSpMat(dgt), msg="dgT2dgC")
    checkEquals(dgc, asSpMat(dgr), msg="dgR2dgC")
}



