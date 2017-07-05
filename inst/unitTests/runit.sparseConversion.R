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
    test.as.dgC2dgC <- function() {
        set.seed(7)
        m <- matrix(0, 5, 5)
        m[sample(length(m), size = 14)] <- rep(1:9, length=14)
        mm <- as(m, "CsparseMatrix")
        checkEquals(mm, asSpMat(mm), msg="dgC2dgC")
    }
    
    ## p35 (ddiMatrix)
    test.as.ddi2dgC <- function() {
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
    }
   
    ## p36 (dgCMatrix)
    test.as.dgC2dgC <- function() {
        m <- Matrix(c(0,0,2:0), 3,5)
        checkEquals(m, asSpMat(m), msg="dgC2dgC")
    }
    
    ## p40 (dgTMatrix)
    test.as.dgT2dgC <- function() {
        m <- Matrix(0+1:28, nrow = 4)
        m[-3,c(2,4:5,7)] <- m[ 3, 1:4] <- m[1:3, 6] <- 0
        mT <- as(m, "dgTMatrix")
        dgc <- as(m, "dgCMatrix")
        checkEquals(dgc, asSpMat(mT), msg="dgT2dgC")
    }
    
    ## p40 (dgTMatrix)
    test.as.dgT2dgC <- function() {
        T2 <- new("dgTMatrix",
                  i = as.integer(c(1,1,0,3,3)),
                  j = as.integer(c(2,2,4,0,0)), x=10*1:5, Dim=4:5)
        dgc <- as(T2, "dgCMatrix")
        checkEquals(dgc, asSpMat(T2), msg="dgT2dgC")
    }
    
    # p42 (ddiMatrix)
    test.as.ddi2dgC <- function() {
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
    }
    
    ## p42 (dgTMatrix)
    test.as.dgT2dgC <- function() {
        M1 <- Matrix(0+0:5, 2,3)
        M <- kronecker(Diagonal(3), M1)
        dgc <- as(M, "dgCMatrix")
        checkEquals(dgc, asSpMat(M), msg="dgT2dgC")
    }
    
}



