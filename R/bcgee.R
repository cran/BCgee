derivative_tD_n_W_inv_n <- function(r1, mu_n, eta_n, X_n, tD_n, W_n_inv, A_n_inv, A_n_inv_d, Working_corr, linkd1, linkd2, varfd1, cluster_condition, cluster_size)
{
	if(cluster_condition==1)
		Working_corr_inv <- solve(Working_corr)
	else
		Working_corr_inv <- solve(Working_corr[1:cluster_size, 1:cluster_size])

	tmp <- (A_n_inv_d*varfd1(mu_n)*linkd1(eta_n)*X_n[,r1])%*%Working_corr_inv%*%A_n_inv
	##r1: component of beta wrt derivative has been taken
	t(linkd2(eta_n)*X_n[,r1]*X_n)%*%W_n_inv - .5*tD_n%*%( tmp+t(tmp) )
} 

derivative_D_n <- function(r1, eta_n, X_n, linkd2)
{
	##r1: component of beta wrt derivative has been taken
	linkd2(eta_n)*X_n[,r1]*X_n
} 

k_r1r2_r3_fun <- function(p, mu_n, eta_n, X_n, tD_n, Cov_n, tD_W_inv_n, W_n_inv, A_n_inv, A_n_inv_d, Working_corr, linkd1, linkd2, varfd1, cluster_condition, cluster_size)
{

	k_r1r2_r3_vals <- lapply(1:p, function(PARS) Reduce("+", lapply(1:length(eta_n), 
																														function(x) derivative_tD_n_W_inv_n(
																														r1=PARS, 			 	
																														mu_n=mu_n[[x]], 
																														eta_n=eta_n[[x]], 
																														X_n=X_n[[x]], 
																														tD_n=tD_n[[x]], 
																														W_n_inv=W_n_inv[[x]], 
																														A_n_inv=A_n_inv[[x]], 
																														A_n_inv_d=A_n_inv_d[[x]], 
																														Working_corr=Working_corr,
																														linkd1=linkd1, linkd2=linkd2, varfd1=varfd1,
																														cluster_condition=cluster_condition,
																														cluster_size=cluster_size[x])%*%Cov_n[[x]]%*%t(tD_W_inv_n[[x]]))))

	k_r1r2_r3_vals <- lapply(1:p, function(z) Reduce("rbind",lapply(k_r1r2_r3_vals,"[",z,)))
	k_r1r2_r3_vals
}



k_r1r2r3_fun <- function(p, mu_n, eta_n, X_n, D_n, tD_n, Cov_n, tD_W_inv_n, W_n_inv, A_n_inv, A_n_inv_d, Working_corr, linkd1, linkd2, varfd1, cluster_condition, cluster_size=NULL)
{
	k_r1r2r3_vals1 <- lapply(1:p, function(PARS) Reduce("+", lapply(1:length(eta_n), 
																														function(x) -derivative_tD_n_W_inv_n(
																														r1=PARS, 			 	
																														mu_n=mu_n[[x]], 
																														eta_n=eta_n[[x]], 
																														X_n=X_n[[x]], 
																														tD_n=tD_n[[x]], 
																														W_n_inv=W_n_inv[[x]], 
																														A_n_inv=A_n_inv[[x]], 
																														A_n_inv_d=A_n_inv_d[[x]], 
																														Working_corr=Working_corr,
																														linkd1=linkd1, linkd2=linkd2, varfd1=varfd1,
																														cluster_condition=cluster_condition,
																														cluster_size=cluster_size[x])%*%D_n[[x]])))

	k_r1r2r3_vals2 <- lapply(1:p, function(PARS) Reduce("+", lapply(1:length(eta_n), 
																														function(x) -tD_W_inv_n[[x]]%*%derivative_D_n(
																														r1=PARS, 			 	
																														eta_n=eta_n[[x]], 
																														X_n=X_n[[x]], linkd2=linkd2))))

	k_r1r2r3_vals1 <- lapply(1:p, function(z){ tmp <- Reduce("rbind",lapply(k_r1r2r3_vals1,"[",z,)); tmp+t(tmp)})
	k_r1r2r3_vals2 <- lapply(1:p, function(z) Reduce("rbind",lapply(k_r1r2r3_vals2,"[",z,)))
	k_r1r2r3_vals <- mapply("+", k_r1r2r3_vals1, k_r1r2r3_vals2, SIMPLIFY = FALSE)

	k_r1r2r3_vals
}


##Link functions and their derivatives
##Supported functions are the same as in package gee (to date 4/14/17)
link_funcs <- function(link_name)
{
	##link function
	link <- switch(link_name,
									Identity = function(x) x,
									Logarithm = function(x) exp(x),
									Logit = function(x) exp(x)/(1+exp(x)),
									Reciprocal = function(x) 1/x,
									Probit = function(x) pnorm(x),
									Cloglog = function(x) 1-exp(-exp(x))
								)

	##first derivative of the link function
	linkd1 <- switch(link_name,
									Identity = function(x) rep(1, length(x)),
									Logarithm = function(x) exp(x),
									Logit = function(x) exp(x)/(1+exp(x))^2,
									Reciprocal = function(x) -1/x^2,
									Probit = function(x) dnorm(x),
									Cloglog = function(x) exp(x-exp(x))
									)

	##second derivative of the link function
	linkd2 <- switch(link_name,
									Identity = function(x) rep(0, length(x)),
									Logarithm = function(x) exp(x),
									Logit = function(x) -exp(x)*(exp(x)-1)/(1+exp(x))^3,
									Reciprocal = function(x) 2/x^3,
									Probit = function(x) -dnorm(x)*x,
									Cloglog = function(x) exp(x-exp(x))*(1-exp(x))
									)

	list(link=link, linkd1=linkd1, linkd2=linkd2)
}


##Variance functions and their derivatives
##Supported functions are the same as in package gee (to date 4/14/17)
var_funcs <- function(var_name)
{
	##variance function
	varf <- switch(var_name,
									Gaussian = function(x) rep(1, length(x)),
									Poisson = function(x) x,
									Binomial = function(x) x*(1-x),
									Gamma = function(x) x^2
								)

	##first derivative of the variance function
	varfd1 <- switch(var_name,
									Gaussian = function(x) rep(0, length(x)),
									Poisson = function(x) rep(1, length(x)),
									Binomial = function(x) 1-2*x,
									Gamma = function(x) 2*x
									)

	list(varf=varf, varfd1=varfd1)
}




BCgee <- function(fit)
{
	##check for class of fit
	if(!any(class(fit)%in%"gee"))
		stop("fit must be retrieved from function gee in package gee.")

		links <- link_funcs(fit$model$link)
		linkf <- links$link
		linkd1 <- links$linkd1
		linkd2 <- links$linkd2
		vars <- var_funcs(fit$model$varfun)
		varf <- vars$varf
		varfd1 <- vars$varfd1

    mnames <- c("", "formula", "data", "offset", "subset", "na.action", "id")
    cnames <- names(fit$call)
    cnames <- cnames[match(mnames, cnames, 0)]
    mcall <- fit$call[cnames]
		mcall[[1]] <- as.name("model.frame")
		model.f <- eval(mcall, parent.frame())
		Terms <- attributes(model.f)$terms
		contrasts <- model.extract(model.f, contrasts)
		X <- model.matrix(Terms, model.f, contrasts)
		offset <- model.offset(model.f)
			if(is.null(offset))
				offset <- 0

		y <- as.matrix(model.extract(model.f, "response"))
		weights <- rep(1, length(y))
		flag_weights <- 0
		if(ncol(y)==2)
		{
			flag_weights <- 1
			weights <- as.vector(y%*%c(1, 1))
			y <- y[,1]
		}

		##linear predictor
		eta_n <- tapply(fit$linear.predictors + offset, fit$id, function(x) x )
		##number of clusters
		n_cluster <- length(eta_n)
		##fitted values
		mu_n <- lapply(eta_n, "linkf" )
		##
		weights_n <- by(weights, fit$id, function(y) y )

		##design matrix
		X_n <- by(X, fit$id, function(y) as.matrix(y) )
		##D matrix
		D_n <- lapply(1:n_cluster, function(x) linkd1(eta_n[[x]])*X_n[[x]] )
		tD_n <- lapply(D_n, "t")
		#
		if(!flag_weights)
		{
			A_n_inv <- lapply(mu_n, function(y) (fit$scale^(-1/2))*diag(varf(y)^(-1/2)))
			A_n_inv_d <- lapply(mu_n, function(y) (fit$scale^(-1/2))*diag(varf(y)^(-3/2)))
		}else{		
			lst_mu_wei <- mapply(function(x,y) list(mu=x,w=y), mu_n, weights_n, SIMPLIFY=FALSE)
			A_n_inv <- lapply(lst_mu_wei, function(Z) (fit$scale^(-1/2))*diag(Z$w^(1/2)*varf(Z$mu)^(-1/2)))
			A_n_inv_d <- lapply(lst_mu_wei, function(Z) (fit$scale^(-1/2))*diag(Z$w^(1/2)*varf(Z$mu)^(-3/2)))
		}

		Working_corr <- fit$working.correlation
		##account for non-constant cluster size
		cluster_size <- table(fit$id)
		cluster_condition <- length(unique(cluster_size))
		if(cluster_condition==1)
		{
			Working_corr_inv <- solve(Working_corr)
			W_n_inv <- lapply(1:n_cluster, function(x) A_n_inv[[x]]%*%Working_corr_inv%*%A_n_inv[[x]])
		}else{
			W_n_inv <- lapply(1:n_cluster, function(x) A_n_inv[[x]]%*%solve(Working_corr[1:cluster_size[x],1:cluster_size[x]])%*%A_n_inv[[x]])
		}

		tD_W_inv_n <- lapply(1:n_cluster, function(x) tD_n[[x]]%*%W_n_inv[[x]])
		residuals_n <- (fit$y-fit$fitted.values*weights)/weights
		residuals_n <- tapply(residuals_n, fit$id, function(x) x )
		Cov_n <- lapply(1:n_cluster, function(x) residuals_n[[x]]%*%t(residuals_n[[x]]))

		###The sandwich (asymptotic variance of estimator of regression parameters)
		Sandwich <- fit$robust.variance
		#The bread in the sandwich
		B <- fit$naive.variance

		p <- length(coefficients(fit))
		k_r1r2_r3_vals <- k_r1r2_r3_fun(p, mu_n, eta_n, X_n, tD_n, Cov_n, tD_W_inv_n, W_n_inv, A_n_inv, A_n_inv_d, 
																		Working_corr, linkd1, linkd2, varfd1, cluster_condition, cluster_size)
		k_r1r2r3_vals <- k_r1r2r3_fun(p, mu_n, eta_n, X_n, D_n, tD_n, Cov_n, tD_W_inv_n, W_n_inv, A_n_inv, A_n_inv_d, 
																		Working_corr, linkd1, linkd2, varfd1, cluster_condition, cluster_size)

		#first piece of the bias
		b1 <- B%*%unlist( lapply(k_r1r2_r3_vals, function(x) sum(x*B)) )
		#second piece of the bias
		b2 <- .5*B%*%unlist( lapply(k_r1r2r3_vals, function(x) sum(x*Sandwich)) )

		result <- fit
		##bias estimates
		result$bias <- as.numeric(b1+b2)
		attr(result$bias,"names") <- attributes(coefficients(fit))$names
		##bias corrected estimates
		result$coefficients <- result$coefficients - result$bias
		##updated linear predictors
		result$linear.predictors <- as.numeric(X%*%result$coefficients)# + offset)
		##updated fitted
		result$fitted.values <- as.numeric(linkf(result$linear.predictors))
		##updated residuals
		result$residuals <- as.numeric(result$y-result$fitted.values*weights)/weights

		##updated naive variance
		eta_n_upd <- tapply(result$linear.predictors, result$id, function(x) x )
		D_n_upd <- lapply(1:n_cluster, function(x) linkd1(eta_n_upd[[x]])*X_n[[x]] )
		tD_n_upd <- lapply(D_n_upd, "t")

		if(!flag_weights)
		{
			A_n_inv_upd <- by(result$fitted.values, fit$id, function(y) (fit$scale^(-.5))*diag(varf(y)^(-.5)))		
		}else{		
			mu_n_upd <- by(result$fitted.values, fit$id, function(y) y )
			lst_mu_wei_upd <- mapply(function(x,y) list(mu=x,w=y), mu_n_upd, weights_n, SIMPLIFY=FALSE)
			A_n_inv_upd <- lapply(lst_mu_wei_upd, function(Z) (fit$scale^(-1/2))*diag(Z$w^(1/2)*varf(Z$mu)^(-1/2)))
		}

		if(cluster_condition==1)
		{
			W_n_inv_upd <- lapply(1:n_cluster, function(x) A_n_inv_upd[[x]]%*%Working_corr_inv%*%A_n_inv_upd[[x]])
		}else{
			W_n_inv_upd <- lapply(1:n_cluster, function(x) A_n_inv_upd[[x]]%*%solve(Working_corr[1:cluster_size[x],1:cluster_size[x]])%*%A_n_inv_upd[[x]])
		}

		tD_W_inv_n_upd <- lapply(1:n_cluster, function(x) tD_n_upd[[x]]%*%W_n_inv_upd[[x]])
		result$naive.variance <- solve(Reduce("+", lapply(1:n_cluster, function(x) tD_W_inv_n_upd[[x]]%*%D_n_upd[[x]])))
		##updated robust variance
		residuals_n_upd <- tapply(result$residuals, result$id, function(x) x )
		Cov_n_upd <- lapply(1:n_cluster, function(x) residuals_n_upd[[x]]%*%t(residuals_n_upd[[x]]))
		V_upd <- Reduce("+", lapply(1:n_cluster, function(x) tD_W_inv_n_upd[[x]]%*%Cov_n_upd[[x]]%*%t(tD_W_inv_n_upd[[x]])))
		result$robust.variance <- result$naive.variance%*%V_upd%*%result$naive.variance

    attr(result, "class") <- c("BCgee")
		result
}


print.BCgee <- function(x, digits = NULL, quote = FALSE, prefix = "", ...)
{
    if(is.null(digits)) digits <- options()$digits else options(digits =
                                                                digits)
#    cat("\n BCgee: BIAS CORRECTED ESTIMATES OF REGRESSION PARAMATERS OF ")
#    cat("\n        GENERALIZED LINEAR MODELS FOR DEPENDENT DATA ")
    cat("\n BCgee: Bias correction based on estimates retrieved from gee package")
    cat("\n       ", x$version, "\n")
    cat("\nModel:\n")
    cat(" Link:                     ", x$model$link, "\n")
    cat(" Variance to Mean Relation:", x$model$varfun, "\n")
    if(!is.null(x$model$M))
        cat(" Correlation Structure:    ", x$model$corstr, ", M =", x$
            model$M, "\n")
    else cat(" Correlation Structure:    ", x$model$corstr, "\n")
    cat("\nCall:\n")
    dput(x$call)                        #       cat("\nTerms:\n")
###        ys <- matrix(rep(as.matrix(x$id, ncol = 1), 5), ncol = 5)
    ys <- matrix(rep(matrix(x$id, ncol = 1), 5), ncol = 5)
    ys[, 2] <- x$y
    ys[, 3] <- x$linear.predictors
    ys[, 4] <- x$fitted.values
    ys[, 5] <- x$residuals
    dimnames(ys) <- list(1:length(x$y), c("ID", "Y", "LP", "fitted",
                                          "Residual")) #       cat("\nFitted Values:\n")
    cat("\nNumber of observations : ", x$nobs, "\n")
    cat("\nMaximum cluster size   : ", x$max.id, "\n")
    nas <- x$nas
    if(any(nas))
        cat("\n\nCoefficients: (", sum(nas),
            " not defined because of singularities)\n", sep = "")
    else cat("\n\nBias-corrected coefficients:\n")
    print(x$coefficients, digits = digits)
    cat("\nEstimated Scale Parameter: ", format(round(x$scale, digits)))
    cat("\nNumber of Iterations: ", x$iterations)
    nc <- min(4, nrow(x$working.correlation))
    cat("\n\nWorking Correlation\n")
    print(x$working.correlation[1:nc, 1:nc], digits = digits)
    cat("\n\nReturned Error Value:\n")
    print(x$error)
    invisible(x)
}


print.summary.BCgee <- function(x, digits = NULL, quote = FALSE, prefix = "", ... )
{
    if(is.null(digits))
        digits <- options()$digits
    else options(digits = digits)
#    cat("\n",x$title)
#    cat("\n",x$version,"\n")
    cat("\n BCgee: Bias correction based on estimates retrieved from gee package")
    cat("\n       ", x$version, "\n")
    cat("\nModel:\n")
    cat(" Link:                     ",x$model$link,"\n")
    cat(" Variance to Mean Relation:",x$model$varfun,"\n")
    if(!is.null(x$model$M))
        cat(" Correlation Structure:    ",x$model$corstr,", M =",x$model$M,"\n")
    else 	cat(" Correlation Structure:    ",x$model$corstr,"\n")
    cat("\nCall:\n")
    dput(x$call)
    cat("\nSummary of Residuals:\n")
    print(x$residual.summary, digits = digits)
    nas <- x$nas
###	if(any(nas))
    if(!is.null(nas) && any(nas))
        cat("\n\nCoefficients: (", sum(nas),
            " not defined because of singularities)\n", sep = "")
    else cat("\n\nCoefficients:\n")
    print(x$coefficients, digits = digits)
    cat("\nEstimated Scale Parameter: ", format(round(x$scale,digits)))
    cat("\nNumber of Iterations: ", x$iterations)
    cat("\n\nWorking Correlation\n")
    print(x$working.correlation,digits=digits)
    if(!is.null(x$naive.correlation)) {
        cat("\n\nNaive Correlation of Estimates:\n")
        print(x$naive.correlation)
    }
    if(!is.null(x$robust.correlation)) {
        cat("\n\nRobust Correlation of Estimates:\n")
        print(x$robust.correlation)
    }
    invisible(x)
}

summary.BCgee <- function(object, correlation = TRUE, ...)
{
    coef <- object$coefficients
    resid <- object$residuals
    p <- object$rank
    if(is.null(p))
        p <- sum(!is.na(coef))
    if(!p) {
        warning("This model has zero rank --- no summary is provided")
        return(object)
    }
    nas <- is.na(coef)
    cnames <- names(coef[!nas])
    coef <- matrix(rep(coef[!nas], 5), ncol = 5)
    dimnames(coef) <- list(cnames, c("Estimate",
                                     "Naive S.E.",  "Naive z",
                                     "Robust S.E.", "Robust z"))
    rse <- sqrt(diag(object$robust.variance))
    nse <- sqrt(diag(object$naive.variance))
    coef[,2] <- nse
    coef[,3] <- coef[,1]/coef[,2]
    coef[,4] <- rse
    coef[,5] <- coef[,1]/coef[,4]
    summary <- list()
    summary$call <- object$call
    summary$version <- object$version
    summary$nobs <- object$nobs
    summary$residual.summary <- quantile(as.vector(object$residuals))
    names(summary$residual.summary) <- c("Min", "1Q", "Median", "3Q", "Max")
    summary$model<- object$model
    summary$title <- object$title
    summary$coefficients <- coef
    summary$working.correlation <- object$working.correlation
    summary$scale <- object$scale
    summary$error <- paste("Error code was", object$error)
    summary$working.correlation <- object$working.correlation
    summary$iterations <- object$iterations
    if ( correlation ) {
        ##	rob.var <- object$robust.variance
        ##	nai.var <- object$naive.variance
        ##	summary$robust.correlation <- rob.var /
        ##	outer(sqrt(diag(rob.var)),sqrt(diag(rob.var)))
        ##	dimnames(summary$robust.correlation) <- list(object$xnames,object$xnames)
        ##	summary$naive.correlation <- nai.var /
        ##	outer(sqrt(diag(nai.var)),sqrt(diag(nai.var)))
        ##	dimnames(summary$naive.correlation) <- list(object$xnames,object$xnames)
    }
    attr(summary,"class") <- "summary.BCgee"
    summary
}

