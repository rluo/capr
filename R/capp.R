capReg <- function(Y, X, nD = 1, method = c("CAP", "CAP-C"), CAP.OC = FALSE, max.itr = 1000, tol = 1e-4, trace = FALSE, score.return = TRUE, gamma0.mat = NULL, ninitial = NULL) {
    n <- length(Y)
    q <- ncol(X)
    p <- ncol(Y[[1]])

    if (is.null(colnames(X))) {
        colnames(X) <- paste0("X", 1:q)
    }

    if (method[1] == "CAP") {
        re1 <- MatReg_QC_opt(Y, X, method = method, max.itr = max.itr, tol = tol, trace = FALSE, score.return = score.return, gamma0.mat = gamma0.mat, ninitial = ninitial)

        Phi.est <- matrix(re1$gamma, ncol = 1)
        beta.est <- matrix(re1$beta, ncol = 1)

        if (score.return) {
            score <- matrix(re1$score, ncol = 1)
        }

        if (nD > 1) {
            Phi.est <- matrix(Phi.est, ncol = 1)
            for (j in 2:nD)
            {
                re.tmp <- NULL
                try(re.tmp <- MatReg_QC_opt2(Y, X, Phi0 = Phi.est, method = method, CAP.OC = CAP.OC, max.itr = max.itr, tol = tol, trace = FALSE, score.return = score.return, gamma0.mat = gamma0.mat, ninitial = ninitial))

                if (is.null(re.tmp) == FALSE) {
                    Phi.est <- cbind(Phi.est, re.tmp$gamma)
                    beta.est <- cbind(beta.est, re.tmp$beta)

                    if (score.return) {
                        score <- cbind(score, re.tmp$score)
                    }
                } else {
                    break
                }
            }
        }

        colnames(Phi.est) <- colnames(beta.est) <- paste0("D", 1:ncol(Phi.est))
        rownames(Phi.est) <- paste0("V", 1:p)
        rownames(beta.est) <- colnames(X)

        if (ncol(Phi.est) > 1) {
            DfD <- diag.level(Y, Phi.est)
        } else {
            DfD <- 1
        }

        if (score.return) {
            colnames(score) <- paste0("D", 1:ncol(Phi.est))
            re <- list(gamma = Phi.est, beta = beta.est, orthogonality = t(Phi.est) %*% Phi.est, DfD = DfD, score = score)
        } else {
            re <- list(gamma = Phi.est, beta = beta.est, orthogonality = t(Phi.est) %*% Phi.est, DfD = DfD)
        }

        return(re)
    } else {
        # ================================================================
        # find common eigenvectors and subject-specific eigenvalues
        Ymat <- NULL
        Group <- NULL
        Tvec <- rep(NA, n)
        for (i in 1:n)
        {
            Tvec[i] <- nrow(Y[[i]])

            Ymat <- rbind(Ymat, Y[[i]])

            Group <- c(Group, rep(i, Tvec[i]))
        }
        re.FCPCA <- FCPCA(Data = Ymat, Group = Group)
        # eigenvalues
        lambda <- re.FCPCA$lambda
        # common eigenvectors
        phi <- re.FCPCA$loadings.common
        for (j in 1:p)
        {
            if (phi[1, j] < 0) {
                phi[, j] <- -phi[, j]
            }
        }
        # ================================================================

        if (method[1] == "CAP-C") {
            optmat <- matrix(NA, p, p)
            colnames(optmat) <- paste0("Dim", 1:p)
            rownames(optmat) <- paste0("BetaDim", 1:p)
            beta.CPC <- matrix(NA, q, p)
            colnames(beta.CPC) <- paste0("Dim", 1:p)
            rownames(beta.CPC) <- colnames(X)
            for (j in 1:p)
            {
                beta.tmp <- MatReg_QC_beta(Y, X, gamma = phi[, j])$beta
                optmat[j, ] <- (apply(lambda, 2, function(x) {
                    return(sum(x * Tvec * exp(-X %*% beta.tmp)))
                }) / apply(lambda, 2, sum)) * n / 2 + sum(Tvec * X %*% beta.tmp) / 2

                beta.CPC[, j] <- beta.tmp
            }

            min.idx <- apply(optmat, 1, which.min)
            sidx <- which(apply(cbind(min.idx, 1:p), 1, function(x) {
                x[1] == x[2]
            }) == TRUE)
            if (length(sidx) > 0) {
                svar <- rep(NA, length(sidx))
                for (j in 1:length(sidx))
                {
                    svar[j] <- sum(exp(X %*% beta.CPC[, sidx[j]]))
                }
                order.idx <- sidx[sort(svar, decreasing = TRUE, index.return = TRUE)$ix]
                if (nD <= length(order.idx)) {
                    Phi.est <- phi[, order.idx[1:nD]]
                    beta.est <- beta.CPC[, order.idx[1:nD]]

                    PC.idx <- order.idx[1:nD]

                    if (score.return) {
                        score <- lambda[, order.idx[1:nD]]
                    }
                } else {
                    Phi.est <- phi[, order.idx]
                    beta.est <- beta.CPC[, order.idx]

                    PC.idx <- order.idx

                    if (score.return) {
                        score <- lambda[, order.idx]
                    }
                }

                if (score.return) {
                    re <- list(gamma = Phi.est, beta = beta.est, orthogonality = t(Phi.est) %*% Phi.est, PC.idx = PC.idx, aPC.idx = order.idx, minmax = TRUE, score = score)
                } else {
                    re <- list(gamma = Phi.est, beta = beta.est, orthogonality = t(Phi.est) %*% Phi.est, PC.idx = PC.idx, aPC.idx = order.idx, minmax = TRUE)
                }
            } else {
                re1 <- MatReg_QC_opt(Y, X, method = method, max.itr = max.itr, tol = tol, trace = FALSE, score.return = score.return, gamma0.mat = gamma0.mat, ninitial = ninitial)

                Phi.est <- matrix(re1$gamma, ncol = 1)
                beta.est <- matrix(re1$beta, ncol = 1)

                if (score.return) {
                    score <- matrix(re1$score, ncol = 1)
                }

                if (nD > 1) {
                    Phi.est <- matrix(Phi.est, ncol = 1)
                    for (j in 2:nD)
                    {
                        re.tmp <- MatReg_QC_opt2(Y, X, Phi0 = Phi.est, method = method, CAP.OC = CAP.OC, max.itr = max.itr, tol = tol, trace = FALSE, score.return = score.return, gamma0.mat = gamma0.mat, ninitial = ninitial)

                        Phi.est <- cbind(Phi.est, re.tmp$gamma)
                        beta.est <- cbind(beta.est, re.tmp$beta)

                        if (score.return) {
                            score <- cbind(score, re.tmp$score)
                        }
                    }
                }

                PC.idx <- rep(NA, ncol(Phi.est))
                for (j in 1:ncol(Phi.est))
                {
                    PC.idx[j] <- which.min(apply(phi, 2, function(x) {
                        return(sqrt(sum((x - Phi.est[, j])^2)))
                    }))
                }
                colnames(Phi.est) <- colnames(beta.est) <- paste0("Dim", PC.idx)
                rownames(Phi.est) <- paste0("V", 1:p)
                rownames(beta.est) <- colnames(X)

                if (score.return) {
                    colnames(score) <- paste0("Dim", PC.idx)
                }

                if (score.return) {
                    re <- list(gamma = Phi.est, beta = beta.est, orthogonality = t(Phi.est) %*% Phi.est, PC.idx = PC.idx, minmax = FALSE, score = score)
                } else {
                    re <- list(gamma = Phi.est, beta = beta.est, orthogonality = t(Phi.est) %*% Phi.est, PC.idx = PC.idx, minmax = FALSE)
                }
            }
        }

        return(re)
    }
}
