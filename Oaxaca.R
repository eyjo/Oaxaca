oaxaca <- function(m1, m2, FE = FALSE) UseMethod("oaxaca")

oaxaca.default <- function(m1, m2, FE = FALSE) {
    # Accuire the full model matrices
    X1 <- model.matrix(m1)
    X2 <- model.matrix(m2)
    # Accuire the parameters
    B1 <- as.matrix(coef(m1))
    B2 <- as.matrix(coef(m2))
    B1[is.na(B1)] <- 0
    B2[is.na(B2)] <- 0
    # Check if we are using FE model
    if (FE)
    {
        B1[1,1] <- m1$intercept
        B2[1,1] <- m2$intercept
    }
    # Variance-covariance matrice of the parameters
    VB1 <- vcov(m1)
    VB2 <- vcov(m2)
    # Check for model mismatch
    allvars <- union(rownames(B1), rownames(B2))
    if (rownames(B1) != allvars || rownames(B2) != allvars) {
        warning("Model mismatch!")
        X1 <- expandMatrix(X1, allvars)
        X2 <- expandMatrix(X2, allvars)
        B1 <- expandMatrix(B1, allvars, TRUE)
        B2 <- expandMatrix(B2, allvars, TRUE)
        VB1 <- expandMatrix(VB1, allvars)
        VB2 <- expandMatrix(VB2, allvars)
        VB1 <- expandMatrix(VB1, allvars, TRUE)
        VB2 <- expandMatrix(VB2, allvars, TRUE)
    }
    k <- NROW(B1)
    # Create the variance-covariance matrix of mean explanatory variables
    VX1 <- cov(X1) / nrow(X1)
    VX2 <- cov(X2) / nrow(X2)
    # Create model weights:
    W <- list("W = 1 (Oaxaca, 1973)" = diag(1, k),
            "W = 0 (Blinder, 1973)" = diag(0, k),
            "W = 0.5 (Reimers 1983)" = diag(0.5, k),
            "W = Omega (Neumark 1988)" = solve(crossprod(X1) + crossprod(X2)) %*% crossprod(X1))
    # The estimate relies on mean values
    X1 <- as.matrix(colMeans(X1, na.rm = TRUE))
    X2 <- as.matrix(colMeans(X2, na.rm = TRUE))
    # Mean difference
    ans <- list(R = mean(model.response(m1$model)) - mean(model.response(m2$model)),
        VR = crossprod(X1, VB1) %*% X1 + crossprod(B1, VX1) %*% B1 + sum(diag(VX1 %*% VB1)) +
            crossprod(X2, VB2) %*% X2 + crossprod(B2, VX2) %*% B2 + sum(diag(VX2 %*% VB2)))
    # The Blinder-Oaxaca variance decomposition
    Q <- sapply(W, function(w) crossprod(X1 - X2, w %*% B1 + (diag(1, k) - w) %*% B2))
    U <- sapply(W, function(w) (crossprod(X1, diag(1, k) - w) + crossprod(X2, w)) %*% (B1 - B2))
    VQ <- sapply(W, function(w)
        sum(diag((VX1 + VX2) %*% (w %*% tcrossprod(VB1, w) + (diag(1, k) - w) %*% tcrossprod(VB2, diag(1, k) - w)))) +
        crossprod(X1 - X2, w %*% tcrossprod(VB1, w) + (diag(1, k) - w) %*% tcrossprod(VB2, diag(1, k) - w)) %*% (X1 - X2) +
        crossprod(w %*% B1 + (diag(1, k) - w) %*% B2, VX1 + VX2) %*% (w %*% B1 + (diag(1, k) - w) %*% B2))
    VU <- sapply(W, function(w)
        sum(diag((crossprod(diag(1, k) - w, VX1) %*% (diag(1, k) - w) + crossprod(w, VX2) %*% w) %*% (VB1 + VB2))) +
        crossprod(crossprod(diag(1, k) - w, X1) + crossprod(w, X2), VB1 + VB2) %*% (crossprod(diag(1, k) - w, X1) + crossprod(w, X2)) +
        crossprod(B1 - B2, crossprod(diag(1, k) - w, VX1) %*% (diag(1, k) - w) + crossprod(w, VX2) %*% w) %*% (B1 - B2))
    # Prep resulting list
    ans$Q <- Q
    ans$U <- U
    ans$VQ <- VQ
    ans$VU <- VU
    ans$W <- W
    ans$call <- match.call()
    class(ans) <- "oaxaca"
    ans
}

print.oaxaca <- function(x) {
    se <- sqrt(x$VQ)
    zval <- x$Q / se
    decomp <- cbind(Explained = x$Q, StdErr = se, "z-value" = zval, "Pr(>|z|)" = 2*pnorm(-abs(zval)))
    cat("\nBlinder-Oaxaca decomposition\n\nCall:\n")
    print(x$call)
    decomp <- matrix(nrow = 1, ncol = 4)
    decomp[1, c(1, 2)] <- c(x$R, sqrt(x$VR))
    decomp[1, 3] <- c(decomp[, 1] / decomp[, 2])
    decomp[1, 4] <- c(2 * pnorm(-abs(decomp[, 3])))
    colnames(decomp) <- c("Difference", "StdErr", "z-value", "Pr(>|z|)")
    rownames(decomp) <- "Mean"
    cat("\n")
    print(zapsmall(decomp))
    cat("\nLinear decomposition:\n")
    for (i in 1:4) {
        cat(paste("\nWeight: ", names(x$W[i]), "\n"))
        decomp <- cbind(Difference = c(x$Q[i], x$U[i]), StdErr = c(sqrt(x$VQ[i]), sqrt(x$VU[i])))
        decomp <- cbind(decomp, "z-value" = decomp[, 1] / decomp[, 2])
        decomp <- cbind(decomp, "Pr(>|z|)" = 2 * pnorm(-abs(decomp[, 3])))
        rownames(decomp) <- c("Explained", "Unexplained")
        print(zapsmall(decomp))
    }
    cat("\n")
}

expandMatrix <- function(x, index, transpose = FALSE) {
    if (transpose)
        x <- t(x)
    tmp <- matrix(0, nrow = nrow(x), ncol = length(index))
    colnames(tmp) <- index
    rownames(tmp) <- rownames(x)
    tmp[, colnames(x)] <- x
    if (transpose)
        tmp <- t(tmp)
    tmp
}

oaxacaOmega <- function(m1, m2, insert = TRUE, FE = FALSE) UseMethod("oaxacaOmega")

oaxacaOmega.default <- function(m1, m2, insert = TRUE, FE = FALSE) {
    X1 <- model.matrix(m1)
    X2 <- model.matrix(m2)
    # Accuire the parameters
    B1 <- as.matrix(coef(m1))
    B2 <- as.matrix(coef(m2))
    B1[is.na(B1)] <- 0
    B2[is.na(B2)] <- 0
    # Check if we are using FE model
    if (FE)
    {
        B1[1,1] <- m1$intercept
        B2[1,1] <- m2$intercept
    }
    # Check for model mismatch
    # if ((NROW(B1) != NROW(B2)) || (rownames(B1) != rownames(B2))) {
    if (insert)
    {
        allvars <- union(rownames(B1), rownames(B2))
        # warning("Model mismatch!")
        X1 <- expandMatrix(X1, allvars)
        X2 <- expandMatrix(X2, allvars)
        B1 <- expandMatrix(B1, allvars, TRUE)
        B2 <- expandMatrix(B2, allvars, TRUE)
    # }
    } else {
        allvars <- intersect(rownames(B1), rownames(B2))
        X1 <- X1[,allvars]
        X2 <- X2[,allvars]
        B1 <- B1[allvars,]
        B2 <- B2[allvars,]
    }
    # The Neumark Omega
    Omega <- solve(crossprod(X1) + crossprod(X2)) %*% crossprod(X1)
    b <- Omega %*% B1 + (diag(1, NROW(B1)) - Omega) %*% B2
    # The estimate relies on mean values
    X1 <- as.matrix(colMeans(X1, na.rm = TRUE))
    X2 <- as.matrix(colMeans(X2, na.rm = TRUE))
    # The Neumark decomposition:
    dw <- list(Explained = crossprod(X1 - X2, b),
        Unexplained = crossprod(X1, B1 - b) + crossprod(X2, b - B2))
    class(dw) <- "oaxacaOmega"
    dw
}

print.oaxacaOmega <- function(x) {
    cat("\nOaxaca linear decomposition (Neumark, 1988):\n")
    cat(paste("Explained = ", x$Explained, ",\tUnexplained = ", x$Unexplained, "\n\n", sep = ""))
}


