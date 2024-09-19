#' Multi-Level Optimal Bayes Function (MLOB)
#'
#' This function estimates between-group parameter beta_b of multilevel latent variable model with control variables
#' using a multivariate regularized Bayesian approach.
#' It is primarily designed for cases where data can be divided into multiple groups
#' and ensures that the data is balanced across groups.
#'
#' @param formula A formula specifying the model (e.g., \code{Y ~ X + C...}) where Y is the dependent variable, X is the target variable (always), and C... represents control variables.
#' @param data A data frame containing the variables in the formula.
#' @param group A positive numeric value representing the number of groups (J) in the data. The number of rows in \code{data} must be divisible by \code{group}.
#' @param balancing Logical. If \code{TRUE}, the data will be balanced across groups. Defaults to \code{FALSE}.
#' @param alpha A numeric value representing the confidence level used to calculate confidence intervals for the estimators. Defaults to \code{0.05}.
#' @param jackknife Logical. If \code{TRUE}, the jackknife resampling method for calculating of the standard error of between-group parameter beta_b will be applied. Defaults to \code{FALSE}.
#' @param ... Additional arguments passed to the function.
#'
#' @details
#' This function checks whether the data is balanced (i.e., whether the same number of individuals are present
#' in each group). If the data is unbalanced and `balancing = FALSE`, the function will stop 
#' and throw an error message. 
#'
#' MultiLevelOptimalBayes (MLOB) is designed for small-sample two-level latent variable models, commonly used in
#' psychology, education, and other disciplines that deal with hierarchical data structures.
#'
#' @return A list containing the results of the Bayesian estimation including model formula, response and predictor variables, and other relevant details.
#'
#' @examples
#' # Example usage with the iris dataset
#' result <- mlob(Sepal.Length ~ Sepal.Width + Petal.Length + Petal.Width, data = iris, group = 3)
#' 
#' # View summary statistics (similar to summary of a linear model)
#' summary(result)
#' ## 
#' ## Call:
#' ## mlob(Sepal.Length ~ Sepal.Width + Petal.Length + Petal.Width, data = data, group = 3) 
#' ## 
#' #' ## Summary of Coefficients:
#' ##                     Estimate Std. Error Lower CI (95%) Upper CI (95%)    Z value  Pr(>|z|) Significance
#' ## beta_b              0.5018467  0.8777217     -1.2184563       2.222150  0.5717605 0.5674842
#' ## gamma_Petal.Length  0.7117494  0.6214602     -0.5062901       1.929789  1.1452857 0.2520908
#' ## gamma_Petal.Width  -0.5576712  1.6291191     -3.7506860       2.635344 -0.3423145 0.7321142
#' ## ---
#' ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#'
#' @export
mlob <- function(formula, data, group, balancing = FALSE, alpha = 0.05, jackknife = FALSE, ...) {
  # Ensure data is a data frame
  if (!is.data.frame(data)) {
    stop("The 'data' argument must be a data frame.")
  }
  
  # Check that 'group' is a number representing the number of groups (J)
  if (!is.numeric(group) || length(group) != 1 || group <= 0) {
    stop("'group' must be a positive numeric value representing the number of groups (J).")
  }
  
  # Ensure the data length is divisible by J (i.e., balanced data)
  n_rows <- nrow(data)
  if (balancing == FALSE){
    if (n_rows %% group != 0) {
      stop(paste("Data is not balanced: the number of rows in the data is not divisible by", group))
    }
  }
  
  
  # Parse the formula (Y ~ X + C...)
  all_vars <- all.vars(formula)
  response_var <- all_vars[1]  # The response variable (Y)
  predictor_var <- all_vars[2]  # The first predictor variable (X)
  control_vars <- all_vars[-c(1, 2)]  # Control variables (C...)
  
  # Extract relevant data columns
  y <- data[[response_var]]
  x <- data[[predictor_var]]
  C <- data[control_vars]
  C <- as.matrix(C)
  
  # Calculate group size (n) and total observations (kn)
  n <- n_rows / group  # Size of each group
  kn <- n * group  # Total observations
  
  # Prepare data_CV list
  data_CV <- list(
    y = y,
    x = x,
    C = C,
    k = group,
    n = n,
    kn = kn,
    kc = length(control_vars)# Number of control variables
  )
  
  ML <- estimate_ML_CV(data_CV) # run ML estimator and get a preliminary estimation of b_b for estimate_Bay_CV
  data_CV$b_b = ML$beta_b_ML_CV # dummy real value of beta_b
  
  # Call estimate_Bay_CV function with data_CV
  Bay <- estimate_Bay_CV(data_CV)
  
  # recalculate SE of Bayesian estimator with jackknife if a
  if (jackknife == TRUE){
    Bay_jackknife <- estimate_Bay_CV_SE_jackknife_individual(data_CV)
    Bay$SE_beta_Bay <-Bay_jackknife$SE_beta_Bay_ML_jackknife_individual
  }
  
  # Generate the result output
  
  # Number of control variables (kc)
  kc <- data_CV$kc
  
  # Create the list of gamma values dynamically
  Coefficients <- data.frame(
    beta_b = Bay$beta_b_Bay,
    t(sapply(1:kc, function(i) Bay$gamma[i]))   # Dynamic number of gamma columns
  )
  
  colnames(Coefficients) <- c("beta_b", paste0("gamma_", control_vars))  # Adjust gamma column names
  
  Standard_Error <- data.frame(
    beta_b = Bay$SE_beta_Bay,
    t(sapply(1:kc, function(i) Bay$SE_gamma[i]))
  )
  
  colnames(Standard_Error) <- c("beta_b", paste0("gamma_", control_vars))
  
  Confidence_Interval <- data.frame(
    Lower = c(Bay$beta_b_Bay - qnorm(1-alpha/2) * Bay$SE_beta_Bay, Bay$gamma - qnorm(1-alpha/2) * Bay$SE_gamma),
    Upper = c(Bay$beta_b_Bay + qnorm(1-alpha/2) * Bay$SE_beta_Bay, Bay$gamma + qnorm(1-alpha/2) * Bay$SE_gamma)
  )
  
  rownames(Confidence_Interval) <- c("beta_b", paste0("gamma_", control_vars))
  
  Confidence_level <- paste0((1 - alpha) * 100, "%")
  
  Z_value <- data.frame(
    beta_b = Bay$beta_b_Bay / Bay$SE_beta_Bay,
    t(sapply(1:kc, function(i) Bay$gamma[i] / Bay$SE_gamma[i]))
  )
  
  colnames(Z_value) <- c("beta_b", paste0("gamma_", control_vars))
  
  p_value <- data.frame(
    beta_b = 2 * (1 - pnorm(abs(Bay$beta_b_Bay / Bay$SE_beta_Bay))),
    t(sapply(1:kc, function(i) 2 * (1 - pnorm(abs(Bay$gamma[i] / Bay$SE_gamma[i])))))
  )
  
  colnames(p_value) <- c("beta_b", paste0("gamma_", control_vars))
  
  
  # Create the dynamic call_info string
  call_info <- paste0("mlob(", deparse(formula), ", data = data, group = ", data_CV$k)
  
  # Conditionally add `balancing` if it's TRUE
  if (!missing(balancing) && balancing) {
    call_info <- paste0(call_info, ", balancing = ", balancing)
  }
  
  # Conditionally add `alpha` if it's not the default value
  if (!missing(alpha) && alpha != 0.05) {
    call_info <- paste0(call_info, ", alpha = ", alpha)
  }
  
  # Conditionally add `jackknife` if it's not TRUE
  if (!missing(jackknife) && !jackknife) {
    call_info <- paste0(call_info, ", jackknife = ", jackknife)
  }
  
  call_info <- paste0(call_info, ")")
  
  
  result <- list(
    Coefficients = Coefficients,
    Standard_Error = Standard_Error,
    Confidence_Interval = Confidence_Interval,
    Confidence_level = Confidence_level,
    Z_value = Z_value,
    p_value = p_value,
    call_info = call_info
    
  )
  
  class(result) <- "mlob_result"  # Assign custom class
  
  return(result)
}


#' @exportS3Method
print.mlob_result <- function(x, ...) {
  # Custom print method for mlob_result
  # Extract call information
  cat("Call:\n", x$call_info, "\n\n")
  
  # Print each section without row names
  cat("Coefficients\n")
  print(x$Coefficients, row.names = FALSE)
  
  cat("\nStandard_Error\n")
  print(x$Standard_Error, row.names = FALSE)
  
  cat("\nConfidence_Interval (", x$Confidence_level, ")\n", sep = "")
  print(x$Confidence_Interval)
  
  cat("\nZ value\n")
  print(x$Z_value, row.names = FALSE)
  
  cat("\np value\n")
  print(x$p_value, row.names = FALSE)
}


#' @exportS3Method
summary.mlob_result <- function(object, ...) {
  # Custom print method for summary_mlob_result
  # Function to assign stars based on p-values
  signif_stars <- function(pval) {
    if (pval < 0.001) {
      return("***")
    } else if (pval < 0.01) {
      return("**")
    } else if (pval < 0.05) {
      return("*")
    } else if (pval < 0.1) {
      return(".")
    } else {
      return(" ")
    }
  }
  
  # Prepare the summary table with significance stars
  summary_table <- data.frame(
    Estimate = as.numeric(object$Coefficients),
    `Std. Error` = as.numeric(object$Standard_Error),
    `Lower CI` = as.numeric(object$Confidence_Interval$Lower),
    `Upper CI` = as.numeric(object$Confidence_Interval$Upper),
    `Z value` = as.numeric(object$Z_value),
    `Pr(>|z|)` = as.numeric(object$p_value),
    Significance = sapply(object$p_value, signif_stars)
  )
  
  # Extract control variable names from the Coefficients (removing the 'beta_b' column)
  control_vars <- colnames(object$Coefficients)[-1]  # Get control variable names (gamma_Age, gamma_Education, etc.)
  
  # Update the row names for the summary table using the actual control variable names
  rownames(summary_table) <- c("beta_b", control_vars)
  
  # Redact column names for consistency
  colnames(summary_table) <- c("Estimate", "Std. Error", paste0("Lower CI (", object$Confidence_level, ")"),
                               paste0("Upper CI (", object$Confidence_level, ")"), "Z value", "Pr(>|z|)", "Significance")
  
  
  # Print the summary header
  cat("Call:\n", object$call_info, "\n\n")
  cat("Summary of Coefficients:\n")
  
  # Print the summary table without row names for a cleaner display
  print(summary_table, row.names = TRUE)
  
  # Add significance codes
  cat("\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  
  invisible(summary_table)  # Return the table invisibly so it doesn't clutter the output
}


estimate_ML_CV <- function(data) {
  
  library(pracma)
  
  ##
  
  tog <- 0 ## Matlab load toggle and change signature data - used for debug
  
  if (tog>0) {
    
    library(R.matlab)
    data <- readMat("data_CV_matlab.mat") # change file name
    #print(data)
    
    data1 <- list(
      k=(data$data[[1]])[1],
      n=(data$data[[2]])[1],
      kc=(data$data[[3]])[1],
      ICC_x=(data$data[[4]])[1],
      ICC_y=(data$data[[5]])[1],
      ICC_C =(data$data[[6]])[1],
      b0=(data$data[[7]])[1],
      b_w=(data$data[[8]])[1],
      b_b=(data$data[[9]])[1],
      gamma=(data$data[[10]])[1],
      kn=(data$data[[11]])[1],
      m_x=(data$data[[12]])[1],
      var_x1=(data$data[[13]])[1],
      var_x2=(data$data[[14]])[1],
      var_e1=(data$data[[15]])[1],
      var_e2=(data$data[[16]])[1],
      
      cov_mat =(data$data[[17]]),
      cov_mat_b =(data$data[[18]]),
      
      x2=(data$data[[19]]),
      x=(data$data[[20]]),
      e2=(data$data[[21]]),
      m_C =(data$data[[22]])[1],
      var_C1 = (data$data[[23]])[1],
      var_C2 =(data$data[[24]])[1],
      
      C2=(data$data[[25]]),
      C=(data$data[[26]]),
      y=(data$data[[27]]),
      x2D=(data$data[[28]])
    )
    
    
    data <-data1
    #print(data)
    
  }
  
  
  data$j <- data$k
  # Reshape the data.x into an n x j matrix
  x <- matrix(data$x, nrow = data$n, ncol = data$k)
  
  av_x <- mean(x)
  av_x_j <- colMeans(x)
  
  SSA <- data$n * sum(av_x_j^2) - data$n * data$k * (av_x^2)
  SSD <- sum(data$x * data$x) - data$n * sum(av_x_j * av_x_j)
  
  
  MSA <- SSA / (data$j - 1)
  MSD <- SSD / ((data$n - 1) * data$j)
  
  tau_x2 <- (MSA - MSD) / data$n
  sigma_x2 <- MSD
  
  # Reshape the data.y into an n x j matrix
  y <- matrix(data$y, nrow = data$n, ncol = data$j)
  
  av_y <- mean(y)
  av_y_j <- colMeans(y)
  
  SPA <- data$n * sum(av_y_j * av_x_j) - data$n * data$j * (av_x * av_y)
  SPD <- sum(data$x *data$y) - data$n * sum(av_y_j * av_x_j)
  
  MPA <- SPA / (data$j - 1)
  MPD <- SPD / ((data$n - 1) * data$j)
  
  tau_yx <- (MPA - MPD) / data$n
  sigma_yx <- MPD
  
  # ML estimation of beta_b
  beta_b_ML <- tau_yx / tau_x2
  
  # Multiple control variables part
  beta_b_ML_tilde <- NULL
  data$C = matrix(data$C)
  if (!is.null(data$kc) && data$kc != 0) {
    
    phi <- matrix(0, nrow = 3, ncol = data$kc)
    data$C = matrix(data$C, nrow = data$k*data$n, ncol = data$kc)
    for (i in 1:data$kc) {
      
      if (data$kc>1){
        # Reshape data$C[,i] into a matrix with n rows and J columns
        C_i_C <- matrix(data$C[, i], nrow = data$n, ncol = data$k)
      }else{
        C_i_C <- matrix(data$C, nrow = data$n, ncol = data$k)
      }
      
      av_C <- mean(C_i_C)
      av_C_j <- colMeans(C_i_C)
      
      SPA_C <- data$n * sum(av_C_j * av_x_j) - data$n * data$j* (av_x * av_C)
      SPD_C <- sum(x * data$C[, i]) - data$n * sum(av_C_j * av_x_j)
      
      MPA_C <- SPA_C / (data$j - 1)
      MPD_C <- SPD_C / ((data$n - 1) * data$j)
      
      tau_C_x <- (MPA_C - MPD_C) / data$n
      sigma_C_x <- MPD_C
      
      phi[2, i] <- tau_C_x / tau_x2
      phi[3, i] <- sigma_C_x / sigma_x2
    }
    
    C <- data$C
    x_b <- rep(av_x_j, each = data$n)
    #x_w <- as.vector(t(x)) - x_b 
    x_w <- as.vector(data$x) - x_b
    C_tilde <- C - cbind(x_b, x_w) %*% (phi[2:3, ])
    
    C_tilde_mat <- cbind(1, C_tilde)
    
    gamma_all <- solve(t(C_tilde_mat) %*% C_tilde_mat) %*% t(C_tilde_mat) %*% data$y
    gamma <- gamma_all[-1]
    
    y_tilde <- data$y - C %*% gamma
    y_tilde <- matrix(y_tilde, nrow = data$n, ncol = data$j)
    
    av_y_tilde <- mean(y_tilde)
    av_y_tilde_j <- colMeans(y_tilde)
    
    SPA_tilde <- data$n * sum(av_y_tilde_j * av_x_j) - data$n * data$j * (av_x * av_y_tilde)
    #SPD_tilde <- sum(x * as.vector(t(y_tilde))) - data$n * sum(av_y_tilde_j * av_x_j)
    SPD_tilde <- sum(as.vector(data$x) * as.vector(y_tilde)) - data$n * sum(av_y_tilde_j * av_x_j)
    
    MPA_tilde <- SPA_tilde / (data$j - 1)
    MPD_tilde <- SPD_tilde / ((data$n - 1) * data$j)
    
    tau_y_tilde_x <- (MPA_tilde - MPD_tilde) / data$n
    sigma_y_tilde_x <- MPD_tilde
    
    # ML estimation of beta_b for the model with control variables
    beta_b_ML_tilde <- tau_y_tilde_x / tau_x2
  }
  
  return(list(beta_b_ML_CV = beta_b_ML_tilde, beta_b_ML_no_CV = beta_b_ML, gamma = gamma, tau_x2 = tau_x2, sigma_x2 = sigma_x2, tau_yx = tau_yx, sigma_yx = sigma_yx))
}


estimate_Bay_CV <- function(data) {
  
  library(pracma)
  
  tog <- 0 ## Matlab load toggle and change signature data - used for debug
  
  if (tog>0) {
    
    library(R.matlab)
    data <- readMat("data_CV_matlab.mat") #change file name
    #print(data)
    
    data1 <- list(
      k=(data$data[[1]])[1],
      n=(data$data[[2]])[1],
      kc=(data$data[[3]])[1],
      ICC_x=(data$data[[4]])[1],
      ICC_y=(data$data[[5]])[1],
      ICC_C =(data$data[[6]])[1],
      b0=(data$data[[7]])[1],
      b_w=(data$data[[8]])[1],
      b_b=(data$data[[9]])[1],
      gamma=(data$data[[10]])[1],
      kn=(data$data[[11]])[1],
      m_x=(data$data[[12]])[1],
      var_x1=(data$data[[13]])[1],
      var_x2=(data$data[[14]])[1],
      var_e1=(data$data[[15]])[1],
      var_e2=(data$data[[16]])[1],
      
      cov_mat =(data$data[[17]]),
      cov_mat_b =(data$data[[18]]),
      
      x2=(data$data[[19]]),
      x=(data$data[[20]]),
      e2=(data$data[[21]]),
      m_C =(data$data[[22]])[1],
      var_C1 = (data$data[[23]])[1],
      var_C2 =(data$data[[24]])[1],
      
      C2=(data$data[[25]]),
      C=(data$data[[26]]),
      y=(data$data[[27]]),
      x2D=(data$data[[28]])
    )
    
    
    data <-data1
    #print(data)
    
  }
  
  
  J <- data$k
  n <- data$n
  
  a <- rep(0, n * J + J + 1)
  
  for (i in 1:(n * J)) {
    a[i] <- -1 / (n * (n - 1) * J)
  }
  
  for (i in (n * J + 1):(n * J + J)) {
    a[i] <- (n * J - 1) / ((n - 1) * (J - 1) * J)
  }
  
  a[n * J + J + 1] <- -J / (J - 1)
  
  A <- diag(a)
  
  # B1 = [A, zeros(n*J+J+1); zeros(n*J+J+1,2*(n*J+J+1))];
  
  # B2 = [zeros(n*J+J+1), A/2; A/2, zeros(n*J+J+1)];
  
  # calculate covariance matrices:
  
  #cov(x):
  
  x <- matrix(data$x, n, J)
  
  av_x <- mean(x)
  av_x_j <- colMeans(x)
  
  SSA <- n * sum(av_x_j^2) - n * J * (av_x^2)
  SSD <- sum(data$x^2) - n * sum(av_x_j^2)
  
  MSA <- SSA / (J - 1)
  MSD <- SSD / ((n - 1) * J)
  
  tau_x2 <- (MSA - MSD) / n
  sigma_x2 <- MSD
  
  #Multiple control variables part
  # checks the existence of control variables and their number is not 0
  
  #Find 3 phi's for every of kc control variables from equation 
  # C = phi1 + phi2*X_b + phi3*X_w + eps. Our aim is phi2 and phi3  
  
  if (!is.null(data$kc) && data$kc != 0) { 
    phi <- matrix(0, 3, data$kc)
    
    C <- list()
    
    for (i in 1:data$kc) {
      C[[i]] <- list()
      
      # 2 cases for number of control variables
      if (data$kc>1){
        # Reshape data$C[,i] into a matrix with n rows and J columns
        C[[i]]$C <- matrix(data$C[, i], nrow = n, ncol = J)
      }else{
        C[[i]]$C <- matrix(data$C, nrow = n, ncol = J)
      }
      
      # Compute the average of the entire C matrix
      C[[i]]$av_C <- mean(C[[i]]$C)
      
      # Compute the average of each column in the C matrix
      C[[i]]$av_C_j <- colMeans(C[[i]]$C)
      
      # Compute SPA_C and SPD_C
      C[[i]]$SPA_C <- n * sum(C[[i]]$av_C_j * av_x_j) - n * J * (av_x * C[[i]]$av_C)
      # 2 cases for number of control variables
      if (data$kc>1){
        C[[i]]$SPD_C <- sum(data$x * data$C[, i]) - n * sum(C[[i]]$av_C_j * av_x_j)
      }else{
        C[[i]]$SPD_C <- sum(data$x * data$C) - n * sum(C[[i]]$av_C_j * av_x_j) 
      }
      # Compute MPA_C and MPD_C
      C[[i]]$MPA_C <- C[[i]]$SPA_C / (J - 1)
      C[[i]]$MPD_C <- C[[i]]$SPD_C / ((n - 1) * J)
      
      # Compute tau_C_x and sigma_C_x
      C[[i]]$tau_C_x <- (C[[i]]$MPA_C - C[[i]]$MPD_C) / n
      C[[i]]$sigma_C_x <- C[[i]]$MPD_C
      
      # Initialize phi as a vector of zeros
      C[[i]]$phi <- numeric(3)
      
      # Compute phi values
      C[[i]]$phi[2] <- C[[i]]$tau_C_x / tau_x2
      C[[i]]$phi[3] <- C[[i]]$sigma_C_x / sigma_x2
      
      # Store the phi vector in the phi matrix
      phi[, i] <- C[[i]]$phi
    }
    
    C <- data$C
    x_b <- matrix(rep(av_x_j, each = n), nrow = n * J, ncol = 1)
    #x_b <- matrix(rep(av_x_j, each = n), nrow = n, byrow = TRUE)
    x_w <- as.vector(matrix(data$x, nrow = n * J, ncol = 1)) - x_b
    C_tilde <- C - cbind(x_b, x_w) %*% phi[2:nrow(phi), ]
    C_tilde_mat <- cbind(rep(1, data$kn), C_tilde)
    
    #Estimate gamma from equation y = gamma1 + gamma*C_tilde + eps. Our
    #interest is vestor of parameters gamma without intersection gamma1
    gamma_all <- solve(t(C_tilde_mat) %*% C_tilde_mat) %*% t(C_tilde_mat) %*% data$y
    
    SE_gamma <- numeric(data$kc)
    
    for (i in 1:data$kc) {
      var_res_gamma <- sum((data$y - C_tilde_mat[, i + 1] * gamma_all[i + 1])^2) / (data$kn - 2)
      SE_gamma[i] <- sqrt(var_res_gamma / ((data$kn - 1) * var(C_tilde_mat[, i + 1])))
    }
    
    gamma <- gamma_all[-1]
    # Find y_tilde - the difference between real and estimated 
    # y (without intersection)
    if (data$kc>1){
      y_tilde <- data$y - C %*% gamma
    }else{
      y_tilde <- data$y - C * gamma
    }
    data$y <- y_tilde
  }
  #cov(x,y)
  y <- matrix(data$y, n, J)
  
  av_y <- mean(y)
  av_y_j <- colMeans(y)
  
  SPA <- n * sum(av_y_j * av_x_j) - n * J * (av_x * av_y)
  SPD <- sum(data$x * data$y) - n * sum(av_y_j * av_x_j)
  
  MPA <- SPA / (J - 1)
  MPD <- SPD / ((n - 1) * J)
  
  tau_yx <- (MPA - MPD) / n
  sigma_yx <- MPD
  
  #cov(y):
  SOA <- n * sum(av_y_j^2) - n * J * (av_y^2)
  SOD <- sum(data$y^2) - n * sum(av_y_j^2)
  
  MOA <- SOA / (J - 1)
  MOD <- SOD / ((n - 1) * J)
  
  tau_y2 <- (MOA - MOD) / n
  sigma_y2 <- MOD
  
  #create Sigma_x matrix of covariances of x_ij, mean(x_j), mean(x)
  
  # create Sigma blockvise with 9 blocks
  
  # block 1: cov(x_ij,x_ij): 
  a1 <- diag(sigma_y2 * rep(1, n))
  #a1 <- diag(sigma_x2, n)
  
  a2 <- matrix(tau_x2, n, n)
  
  a <- a1 + a2
  block1 <- kronecker(diag(J), a)
  # Block 2: cov(y_ij, mean(y_j))
  
  block2 <- kronecker(diag(J), (tau_x2 + sigma_x2 / n) * matrix(1, n, 1))
  #block2 <- kronecker(diag(J), (tau_y2 + sigma_y2 / n) * rep(1, n))
  
  block3 <- matrix((tau_x2 / J + sigma_x2 / (n * J)), n * J, 1)
  #block3 <- (tau_y2 / J + sigma_y2 / (n * J)) * rep(1, n * J)
  block4 <- t(block2)
  
  block5 <- (tau_x2 + sigma_x2 / n) * diag(J)
  
  block6 <- (tau_x2 / J + sigma_x2 / (n * J)) * matrix(1, J, 1)
  #block6 <- (tau_y2 / J + sigma_y2 / (n * J)) * rep(1, J
  block7 <- t(block3)
  
  block8 <- t(block6)
  
  block9 <- tau_x2 / J + sigma_x2 / (n * J)
  
  Sigma_x <- rbind(cbind(block1, block2, block3), cbind(block4, block5, block6), cbind(block7, block8, block9))
  
  a1 <- diag(sigma_yx, n)
  a2 <- matrix(tau_yx, n, n)
  
  a <- a1 + a2
  block1 <- kronecker(diag(J), a)
  
  block2 <- kronecker(diag(J), (tau_yx + sigma_yx / n) * matrix(1, n, 1))
  
  block3 <- matrix((tau_yx / J + sigma_yx / (n * J)), n * J, 1)
  
  block4 <- t(block2)
  
  block5 <- (tau_yx + sigma_yx / n) * diag(J)
  
  block6 <- (tau_yx / J + sigma_yx / (n * J)) * matrix(1, J, 1)
  
  block7 <- t(block3)
  
  block8 <- t(block6)
  
  block9 <- tau_yx / J + sigma_yx / (n * J)
  
  Sigma_yx <- rbind(cbind(block1, block2, block3), cbind(block4, block5, block6), cbind(block7, block8, block9))
  
  a1 <- diag(sigma_y2, n)
  a2 <- matrix(tau_y2, n, n)
  
  a <- a1 + a2
  block1 <- kronecker(diag(J), a)
  
  block2 <- kronecker(diag(J), (tau_y2 + sigma_y2 / n) * matrix(1, n, 1))
  
  block3 <- matrix((tau_y2 / J + sigma_y2 / (n * J)), n * J, 1)
  
  block4 <- t(block2)
  
  block5 <- (tau_y2 + sigma_y2 / n) * diag(J)
  
  block6 <- (tau_y2 / J + sigma_y2 / (n * J)) * matrix(1, J, 1)
  
  block7 <- t(block3)
  
  block8 <- t(block6)
  
  block9 <- tau_y2 / J + sigma_y2 / (n * J)
  
  Sigma_y <- rbind(cbind(block1, block2, block3), cbind(block4, block5, block6), cbind(block7, block8, block9))
  
  #spectral decomposition of Sigma_x
  
  #compute eigenvalues matrix D_x of Sigma_x
  
  # Initialize D_x
  D_x <- rep(0, n * J + J + 1)
  
  # Set values for specific indices
  D_x[(J + 2):(n * J + 1)] <- sigma_x2
  D_x[(n * J + 2):(n * J + J)] <- (n + 1) * (tau_x2 + sigma_x2 / n)
  D_x[n * J + J + 1] <- (n * J + J + 1) * (tau_x2 + sigma_x2 / n) / J
  
  # Create diagonal matrix D_x
  D_x <- diag(D_x)
  
  # Initialize V_yx matrix
  V_yx <- matrix(0, n * J + J + 1, n * J + J + 1)
  
  # lambda_i = 0 ((J+1) pieces)
  for (i in 1:J) {
    V_yx[((n * (i - 1) + 1):(n * i)), i] <- -1 / sqrt(n * (n + 1))
    V_yx[(n * J + i), i] <- sqrt(n / (n + 1))
  }
  
  V_yx[1:(n * J + J), J + 1] <- -1 / sqrt((n * J + J) * (n * J + J + 1))
  V_yx[(n * J + J + 1), J + 1] <- sqrt((n * J + J) / (n * J + J + 1))
  
  # lambda_i = sigma_yx ((n-1)J pieces)
  AU <- matrix(0, n, n - 1)
  for (i in 1:(n - 1)) {
    AU[1:i, i] <- -1 / sqrt(i * (i + 1))
    AU[(i + 1), i] <- i / sqrt(i * (i + 1))
  }
  
  for (j in 1:J) {
    V_yx[(n * (j - 1) + 1):(n * j), (J + 2 + (n - 1) * (j - 1)):(J + 1 + (n - 1) * j)] <- AU
  }
  
  # lambda_i = (n + 1)(tau_yx + sigma_yx / n) ((J - 1) pieces)
  for (i in 1:(J - 1)) {
    V_yx[1:(i * n), (n * J + 1 + i)] <- -1 / sqrt(i * (i + 1) * (n + 1))
    V_yx[(n * J + 1):(n * J + i), (n * J + 1 + i)] <- -1 / sqrt(i * (i + 1) * (n + 1))
    V_yx[(i * n + 1):((i + 1) * n), (n * J + 1 + i)] <- i / sqrt(i * (i + 1) * (n + 1))
    V_yx[(n * J + i + 1), (n * J + 1 + i)] <- i / sqrt(i * (i + 1) * (n + 1))
  }
  
  # lambda_nJ + J + 1 = (n * J + J + 1)(tau_yx + sigma_yx / n) / J (1 piece)
  V_yx[1:(n * J + J + 1), (n * J + J + 1)] <- 1 / sqrt(n * J + J + 1)
  
  # Spectral decomposition of Sigma_y
  
  # Initialize D_y
  D_y <- rep(0, n * J + J + 1)
  
  # Set values for specific indices
  D_y[(J + 2):(n * J + 1)] <- sigma_y2
  D_y[(n * J + 2):(n * J + J)] <- (n + 1) * (tau_y2 + sigma_y2 / n)
  D_y[n * J + J + 1] <- (n * J + J + 1) * (tau_y2 + sigma_y2 / n) / J
  
  # Create diagonal matrix D_y
  D_y <- diag(D_y)
  
  # Initialize V_y matrix
  V_y <- matrix(0, n * J + J + 1, n * J + J + 1)
  
  # lambda_i = 0 ((J + 1) pieces)
  for (i in 1:J) {
    V_y[((n * (i - 1) + 1):(n * i)), i] <- -1 / sqrt(n * (n + 1))
    V_y[(n * J + i), i] <- sqrt(n / (n + 1))
  }
  
  V_y[1:(n * J + J), J + 1] <- -1 / sqrt((n * J + J) * (n * J + J + 1))
  V_y[(n * J + J + 1), J + 1] <- sqrt((n * J + J) / (n * J + J + 1))
  
  # lambda_i = sigma_yx ((n - 1)J pieces)
  for (i in 1:(n - 1)) {
    AU[1:i, i] <- -1 / sqrt(i * (i + 1))
    AU[(i + 1), i] <- i / sqrt(i * (i + 1))
  }
  
  for (j in 1:J) {
    V_y[(n * (j - 1) + 1):(n * j), (J + 2 + (n - 1) * (j - 1)):(J + 1 + (n - 1) * j)] <- AU
  }
  V_x <- matrix(0, nrow = n*sum(J)+sum(J)+1, ncol = n*sum(J)+sum(J)+1)
  
  # lambda_i = (n + 1)(tau_yx + sigma_yx / n) ((J - 1) pieces)
  for (i in 1:(J - 1)) {
    V_x[1:(i * n), (n * J + 1 + i)] <- -1 / sqrt(i * (i + 1) * (n + 1))
    V_x[(n * J + 1):(n * J + i), (n * J + 1 + i)] <- -1 / sqrt(i * (i + 1) * (n + 1))
    V_x[(i * n + 1):((i + 1) * n), (n * J + 1 + i)] <- i / sqrt(i * (i + 1) * (n + 1))
    V_x[(n * J + i + 1), (n * J + 1 + i)] <- i / sqrt(i * (i + 1) * (n + 1))
  }
  
  # lambda_nJ + J + 1 = (n * J + J + 1)(tau_yx + sigma_yx / n) / J (1 piece)
  V_x[1:(n * J + J + 1), (n * J + J + 1)] <- 1 / sqrt(n * J + J + 1)
  # check the correctness of spectral decomposition:
  
  # diff_x = sum(sum((Sigma_x - V_x*D_x*V_x').^2));
  
  # disp(['difference between covariance matrix for x and its decomposition is ', num2str(diff_x)]);
  
  
  # Spectral decomposition of Sigma_yx
  
  # compute eigenvalues matrix D_yx of Sigma_yx
  D_yx <- matrix(0, n * sum(J) + sum(J) + 1, 1)
  
  for (i in ((J+2 ): ((n*J)+1))){
    D_yx[i,] <- sigma_yx
  }
  
  for (i in (((n*J)+2) : ((n*J)+J))){
    D_yx[i,] <- (n+1)*(tau_yx + sigma_yx/n)
  }
  
  D_yx[(n*J)+J+1,] <- ((n*J)+J+1)*(tau_yx + sigma_yx/n)/J
  
  D_yx <- diag(D_yx[,1])
  
  # compute eigenvectors of Sigma_yx
  
  # create matrix of eigenvectors V_yx for covariance matrix Sigma_yx
  
  V_yx <- matrix(0, nrow = n*sum(J)+sum(J)+1, ncol = n*sum(J)+sum(J)+1)
  
  
  # lambda_i = 0 ((J+1) pieces)
  
  for (i in 1:(length(seq_len(J)))){
    V_yx[((n*(i-1))+1):(n*i),i] <- -1/sqrt(n*(n+1))
    V_yx[(n*J)+i,i] <- sqrt(n/(n+1))
  }
  
  V_yx[1:((n*J)+J),J+1] <- -1/sqrt((n*J+J)*(n*J+J+1))
  V_yx[(n*J)+J+1,J+1] <- sqrt((n*J+J)/(n*J+J+1))
  
  # lambda_i = sigma_yx ((n-1)J pieces)
  # auxiliary matrix AU
  
  AU <- matrix(0, nrow = n, ncol = n - 1)
  
  for (i in 1:(n-1)){
    AU[1:i,i] <- -1/sqrt(i*(i+1))
    AU[i+1,i] <- i/sqrt(i*(i+1))
  }
  
  for (j in 1:length(seq_len(J))){
    V_y[(n*(j-1)+1):(n*j),(J+2+((n-1)*(j-1))):(J+1+((n-1)*j))] <- AU
  }
  # lambda_i = (n+1)(tau_yx + sigma_yx/n) ((J-1) pieces)
  
  for (i in 1:(J-1)){
    V_y[1:(i*n),(n*J)+1+i] <- -1/sqrt(i*(i+1)*(n+1))
    V_y[((n*J)+1):((n*J)+i),(n*J)+1+i] <- -1/sqrt(i*(i+1)*(n+1))
    V_y[(i*n+1):((i+1)*n),(n*J)+1+i] <- i/sqrt(i*(i+1)*(n+1))
    V_y[(n*J)+i+1,(n*J)+1+i] <- i/sqrt(i*(i+1)*(n+1))
  }
  
  # lambda_nJ+J+1 = (nJ+J+1)(tau_yx + sigma_yx/n)/J (1 piece)
  
  V_yx[1:((n*J)+J+1),(n*J)+J+1] <- 1/sqrt(n*J+J+1)
  
  # check the correctness of spectral decomposition:
  
  # diff_yx = sum(sum((Sigma_yx - V_yx*D_yx*V_yx').^2));
  
  # disp(['difference between covariance matrix for x and y and its decomposition is ', num2str(diff_yx)]);
  
  
  # Spectral decomposition of Sigma_y
  
  # compute eigenvalues matrix D_y of Sigma_y
  
  D_y <- matrix(0, nrow = n*sum(J)+sum(J)+1, ncol = 1)
  
  for (i in ((J+2) : ((n*J)+1))){
    D_y[i,] <- sigma_y2
  }
  
  for (i in ((n*J+2) : ((n*J)+J))){
    D_y[i,] <- (n+1)*(tau_y2 + sigma_y2/n)
  }
  
  D_y[(n*J)+J+1,] <- ((n*J)+J+1)*(tau_y2 + sigma_y2/n)/J
  
  D_y <- diag(D_y[,1])
  
  # compute eigenvectors of Sigma_y
  
  # create matrix of eigenvectors V_y for covariance matrix Sigma_y
  
  V_y <- matrix(0, nrow = n*sum(J)+sum(J)+1, ncol = n*sum(J)+sum(J)+1)
  
  # lambda_i = 0 ((J+1) pieces)
  
  for (i in 1:(length(seq_len(J)))){
    V_y[((n*(i-1))+1):(n*i),i] <- -1/sqrt(n*(n+1))
    V_y[(n*J)+i,i] <- sqrt(n/(n+1))
  }
  
  V_y[1:((n*J)+J),J+1] <- -1/sqrt((n*J+J)*(n*J+J+1))
  V_y[(n*J)+J+1,J+1] <- sqrt((n*J+J)/(n*J+J+1))
  
  # lambda_i = sigma_yx ((n-1)J pieces)
  
  # auxiliary matrix AU
  
  AU <- matrix(0, nrow = n, ncol = n - 1)
  for (i in 1:(n-1)){
    AU[1:i,i] <- -1/sqrt(i*(i+1))
    AU[i+1,i] <- i/sqrt(i*(i+1))
  }
  
  for (j in 1:J){
    V_y[(n*(j-1)+1):(n*j),(J+2+((n-1)*(j-1))):(J+1+((n-1)*j))] <- AU
  }
  
  # lambda_i = (n+1)(tau_yx + sigma_yx/n) ((J-1) pieces)
  
  for (i in 1:(J-1)){
    V_y[1:(i*n),(n*J)+1+i] <- -1/sqrt(i*(i+1)*(n+1))
    V_y[((n*J)+1):((n*J)+i),(n*J)+1+i] <- -1/sqrt(i*(i+1)*(n+1))
    V_y[(i*n+1):((i+1)*n),(n*J)+1+i] <- i/sqrt(i*(i+1)*(n+1))
    V_y[(n*J)+i+1,(n*J)+1+i] <- i/sqrt(i*(i+1)*(n+1))
  }
  
  # lambda_nJ+J+1 = (nJ+J+1)(tau_yx + sigma_yx/n)/J (1 piece)
  
  V_y[1:((n*J)+J+1),(n*J)+J+1] <- 1/sqrt(n*J+J+1)
  
  # check the correctness of spectral decomposition:
  
  # diff_y = sum(sum((Sigma_y - V_y*D_y*V_y').^2));
  
  # disp(['difference between covariance matrix for y and its decomposition is ', num2str(diff_y)]);
  
  
  # calculate the eigenvalues of transformed sum Z_x'*A*Z_x =
  # W_x'*(Sigma_x^0.5)*A*(Sigma_x^0.5)*W_x = {use spectral decomposition and
  # substitute sqrt(D) with S=sqrt(D), because D is diagonal} =
  # = W_x*V_x*S_x*V_x'*A*V_x*S_x*V_x'*W_x = W_x*C_x*W_x
  # C_x = V_x*S_x*V_x'*A*V_x*S_x*V_x' with internal part
  # L1 = S_x*V_x'*A*V_x*S_x  - diagonal -> gives eigenvalues of C_x, which we use
  # to calculate MSE:
  
  S_x <- sqrt(D_x)
  t_v_x <-t(V_x)
  L1 <- S_x %*% t_v_x %*% A %*% V_x %*% S_x 
  
  L1 <- diag(L1)
  
  L1 <- Re(L1)# could be complex insignificant tails, because of computer precision
  
  for (i in 1:(n*J+J+1)){ # delete computer precision tails
    if (abs(L1[i])<0.0000001){
      L1[i] <- 0
    }
  }
  
  # calculate the eigenvalues of transformed sum Z_yx'*A*Z_yx =
  # W_x'*(Sigma_x^0.5)*A*(Sigma_x^0.5)*W_y = {use spectral decomposition and
  # substitute sqrt(D) with S=sqrt(D), because D is diagonal} =
  # {W_x, W_y ~N(0,I)} = W_x'*V_x*S_x*V_x'*A*V_x*S_y*V_y'*W_y =
  # {H_x = V_x'*W_x ~ N(0,I), H_y = V_y'*W_y ~ N(0,I)} =
  # = H_x'*S_x*V_x'*A*V_x*S_y*H_y = H_x'*Q*H_y = {H = [H_x, H_y]} = H'*Q1*H1
  # = H' * EH^(-0.5) * EH^(0.5) * Q1 * EH^(0.5) * EH^(-0.5) *H =
  # = {EH^(-0.5) * H = H1 ~ N(0,I)} = H1' *EH^(0.5) * Q1 * EH^(0.5) * H1 =
  # {use spectral decomposition and substitute: EH = V_H*D_H*V_H' ->
  # EH^(0.5) = V_H*S_H*V_H' with with S_H=sqrt(D_H), because D_H is diagonal}
  # H1' * V_H * S_H * V_H' * Q1 * V_H * S_H * V_H' * H1 =
  # {H1_tilde = V_H' * H1 ~ N(0,I)} =
  # = H1_tilde' * S_H * V_H' * Q1 * V_H * S_H * H1_tilde, with internal part
  # L2 = S_H * V_H' * Q1 * V_H * S_H  - diagonal -> gives eigenvalues of L2,
  # which we use to calculate MSE:
  
  S_y <- sqrt(D_y)
  
  Q <- S_x %*% t_v_x %*% A %*% V_y %*% S_y ## wrong
  
  #
  upper_left <- matrix(0, nrow(Q), ncol(Q))
  lower_right <- matrix(0, nrow(Q), ncol(Q))
  
  upper_right <- Q / 2
  lower_left <- Q / 2
  
  Q1 <- rbind(cbind(upper_left, upper_right), cbind(lower_left, lower_right))
  
  #
  
  
  #
  upper_left <- diag(nrow(D_yx))
  lower_right <- diag(nrow(D_yx))
  
  pinv_S_x_S_y_D_yx <- pinv(S_x) %*% pinv(S_y) %*% D_yx
  
  
  upper_right <- pinv_S_x_S_y_D_yx
  lower_left <- pinv_S_x_S_y_D_yx
  
  EH <- rbind(cbind(upper_left, upper_right), cbind(lower_left, lower_right))
  #
  
  
  D1 <- EH[(1:(n*J+J+1)), ((n*J+J+2):(2*(n*J+J+1)))]
  D_H <- matrix(0, nrow = nrow(EH), ncol = ncol(EH))
  for (i in 1:(n*J+J+1)){
    D_H[2*i-1,2*i-1] <- 1 + D1[i,i]
    D_H[2*i,2*i] <- 1 - D1[i,i]
  }
  S_H <- ifelse(D_H>0,sqrt(D_H),0)
  rm(D1)
  
  ##-- 
  
  
  V_H <- matrix(0, nrow = nrow(D_H), ncol = ncol(D_H))
  
  for (i in 1:(n*J+J+1)){
    V_H[i,(2*i)-1] <- 1/sqrt(2)
    V_H[i, 2*i] <- 1/sqrt(2)
    V_H[(n*J)+J+1 + i, (2*i)-1] <- 1/sqrt(2)
    V_H[(n*J)+J+1 + i, 2*i] <- -1/sqrt(2)
  }
  
  L2 <- S_H %*% t(V_H) %*% Q1 %*% V_H %*% S_H
  
  L2 <- diag(L2)
  L2 <- Re(L2)# could be complex insignificant tails, because of computer precision
  
  for (i in 1:(2*(n*J+J+1))){# delete computer precision tails
    if (abs(L2[i])<0.0000001){
      L2[i] <- 0
    }
  }
  
  
  
  
  # Use L1 and L2 and ML estimator of beta_b as coefficients for computing optimal MSE
  
  beta_b_ML <- tau_yx / tau_x2
  
  K_sum1 <- sum(L1)^2/(2*(t(L1)%*%L1))# unused
  
  T_sum1 <- (t(L1)%*%L1)/sum(L1)# unused
  
  K_sum2 <- sum(L2)^2/(2*(t(L2)%*%L2))
  
  T_sum2 <- (t(L2)%*%L2)/sum(L2)
  
  beta_b <- data$b2# real value of between parameter ##check if b_b is same as b2
  
  
  
  # MSE of ML estimated parameter beta_b_ML
  
  MSE_add <- (((K_sum2 * T_sum2)/ (T_sum1*(K_sum1-1))) - beta_b)^2 +
    ((K_sum2 * T_sum2^2 * (K_sum1 + K_sum2 - 1)) / (T_sum1^2 * (K_sum1-1)^2 * (K_sum1-2)))
  
  # MSE_add = (K_sum1^2 * T_sum1/ (K_sum2 * T_sum2*(K_sum1-1)) - beta_b)^2 +...
  #    K_sum1^4 * T_sum1^2 * (K_sum1+K_sum2-1)/ (K_sum2^3 * T_sum2^2 * (K_sum1-1)^2 * (K_sum1-2)); # original
  
  
  # optimize MSE with beta_b_ML
  
  
  Tau02 <- kronecker(    matrix(1,  nrow = length(seq(0, 1, by = 0.01)), ncol = 1)   , t( matrix(seq(0.05, 10, by = 0.05), nrow = 1) )  )# restrict search interval of tau02 to 10
  
  W <-kronecker(t(matrix(seq(0, 1, by = 0.01), nrow=1)), matrix(1,nrow = length(seq(0.5, by = 0.05)), ncol=1))
  ## recheck???  
  
  MSE_ML <- matrix(0, nrow = length(W), ncol = 1)
  
  
  
  # use grid search to find w and tau02 that give smallest MSE
  for (i in 1:length(W)){
    MSE_ML[i] <- ((K_sum2 * T_sum2^2 * (K_sum2+1) * ((1-W[i])*Tau02[i] + W[i]*sum(L1)))/
                    ((((1-W[i])*Tau02[i] + W[i]*sum(L1))^2 - 2*W[i]^2*(t(L1)%*%L1)) *
                       (((1-W[i])*Tau02[i] + W[i]*sum(L1))^2 - 4*W[i]^2*(t(L1)%*%L1))) - 
                    (2*beta_b_ML *K_sum2 * T_sum2 * ((1-W[i])*Tau02[i] + W[i]*sum(L1)))/
                    ((((1-W[i])*Tau02[i] + W[i]*sum(L1))^2) - (2*W[i]^2*(t(L1)%*%L1)))+ beta_b_ML)
    
  }
  
  # optimize MSE with beta_b
  
  MSE <- matrix(0, nrow = length(W), ncol = 1)
  
  
  # use grid search to find w and tau02 that give smallest MSE
  for (i in 1:length(W)){
    MSE[i] <- ((K_sum2 * T_sum2^2 * (K_sum2+1) * ((1-W[i])*Tau02[i] + W[i]*sum(L1)))/
                 ((((1-W[i])*Tau02[i] + W[i]*sum(L1))^2 - 2*W[i]^2*(t(L1)%*%L1)) *
                    (((1-W[i])*Tau02[i] + W[i]*sum(L1))^2 - 4*W[i]^2*(t(L1)%*%L1))) - 
                 (2*data$b_b *K_sum2 * T_sum2 * ((1-W[i])*Tau02[i] + W[i]*sum(L1)))/
                 ((((1-W[i])*Tau02[i] + W[i]*sum(L1))^2) - (2*W[i]^2*(t(L1)%*%L1)))+data$b_b)
    
    
    
  }
  
  
  # use MSE to find beta_b_Bay with smallest MSE:
  
  #
  ind_MSEbetabBay <- which.min(MSE)
  MSE_beta_b_Bay <- MSE[ind_MSEbetabBay]
  #
  
  
  w_opt <- W[ind_MSEbetabBay]
  tau02_opt <- Tau02[ind_MSEbetabBay]
  
  beta_b_Bay <- tau_yx/((1-w_opt)*tau02_opt + w_opt*tau_x2)
  
  
  # use MSE_ML to find beta_b_Bay_ML with smallest MSE:
  
  #
  ind_MSEbetabBayML <- which.min(MSE_ML)
  MSE_beta_b_Bay_ML <- MSE_ML[ind_MSEbetabBayML]
  #
  
  w_opt_ML <- W[ind_MSEbetabBayML]
  tau02_opt_ML <- Tau02[ind_MSEbetabBayML]
  
  beta_b_Bay_ML <- tau_yx/((1-w_opt_ML)*tau02_opt_ML + w_opt_ML*tau_x2)
  
  
  
  # Compute standard errors of beta_b_Bay, beta_b_Bay_ML, beta_b_ML using their distributions
  
  if (w_opt != 0) {
    K_B_Bay <- (w_opt * T_sum1 * K_sum1 + (1 - w_opt) * tau02_opt)^2 / (w_opt^2 * T_sum1^2 * K_sum1)
    
    T_B_Bay <- (w_opt^2 * T_sum1^2 * K_sum1) / (w_opt * T_sum1 * K_sum1 + (1 - w_opt) * tau02_opt)
  } else {
    K_B_Bay <- 4  # to avoid negative values in SE_beta_Bay
    
    T_B_Bay <- 0.25
  }
  
  if (w_opt_ML != 0) {
    K_B_Bay_ML <- (w_opt_ML * T_sum1 * K_sum1 + (1 - w_opt_ML) * tau02_opt_ML)^2 / (w_opt_ML^2 * T_sum1^2 * K_sum1)
    
    T_B_Bay_ML <- (w_opt_ML^2 * T_sum1^2 * K_sum1) / (w_opt_ML * T_sum1 * K_sum1 + (1 - w_opt_ML) * tau02_opt_ML)
  } else {
    K_B_Bay_ML <- 4  # to avoid negative values in SE_beta_Bay_ML
    
    T_B_Bay_ML <- 0.25
  }
  
  
  
  SE_beta_ML <- abs((T_sum2 / (T_sum1 * (K_sum1 - 1))) * sqrt(abs((K_sum2 * (K_sum1 + K_sum2 - 1)) / (K_sum1 - 2))))
  
  SE_beta_Bay <- abs((T_sum2 / (T_B_Bay * (K_B_Bay - 1))) * sqrt(abs((K_sum2 * (K_B_Bay + K_sum2 - 1)) / (K_B_Bay - 2))))
  
  SE_beta_Bay_ML <- abs((T_sum2 / (T_B_Bay_ML * (K_B_Bay_ML - 1))) * sqrt(abs((K_sum2 * (K_B_Bay_ML + K_sum2 - 1)) / (K_B_Bay_ML - 2))))
  
  
  #Bay <- 0.1
  Bay <- list()
  Bay$A <- A
  #Bay$B2 = B2;
  Bay$Sigma_x <- Sigma_x
  Bay$Sigma_yx <- Sigma_yx
  Bay$Sigma_y <- Sigma_y
  Bay$tau_x2 <- tau_x2
  Bay$sigma_x2 <- sigma_x2
  Bay$tau_yx <- tau_yx
  Bay$sigma_yx <- sigma_yx
  Bay$tau_y2 <- tau_y2
  Bay$sigma_y2 <- sigma_y2
  Bay$D_x <- D_x
  Bay$V_x <- V_x
  Bay$D_yx <- D_yx
  Bay$V_yx <- V_yx
  Bay$D_y <- D_y
  Bay$V_y <- V_y
  Bay$L1 <- L1
  Bay$L2 <- L2
  Bay$MSE <- MSE
  Bay$MSE_ML <- MSE_ML
  Bay$W <- W
  Bay$Tau02 <- Tau02
  Bay$beta_b_ML <- beta_b_ML
  Bay$beta_b_Bay <- beta_b_Bay
  Bay$beta_b_Bay_ML <- beta_b_Bay_ML
  Bay$gamma <- gamma
  Bay$MSE_beta_b_Bay <- MSE_beta_b_Bay
  Bay$MSE_beta_b_Bay_ML <- MSE_beta_b_Bay_ML
  Bay$MSE_add <- MSE_add
  Bay$SE_beta_ML = SE_beta_ML
  Bay$SE_beta_Bay = SE_beta_Bay
  Bay$SE_beta_Bay_ML = SE_beta_Bay_ML
  Bay$SE_gamma = SE_gamma
  
  Bay_list <- list( 'A'=A,  'Sigma_x'=Sigma_x, 'Sigma_y'=Sigma_y, 'Sigma_yx'=Sigma_yx, 'tau_x2'=tau_x2, 
                    'sigma_x2'=sigma_x2, 'tau_y2'=tau_y2, 'sigma_y2'=sigma_y2,  'tau_yx'=tau_yx, 
                    'sigma_yx'=sigma_yx,'D_x'=D_x,  'D_y'=D_y, 'D_yx'=D_yx,'V_x'=V_x, 'V_y'=V_y, 
                    'V_yx'=V_yx,'L1'=L1, 'L2'=L2,  'MSE'=MSE,  'MSE_ML'=MSE_ML,'W'=W,'Tau02'=Tau02, 
                    'beta_b_ML'=beta_b_ML, 'beta_b_Bay'=beta_b_Bay,'beta_b_Bay_ML'=beta_b_Bay_ML, 
                    'gamma' = gamma, 'MSE_beta_b_Bay'=MSE_beta_b_Bay,'MSE_beta_b_Bay_ML'=MSE_beta_b_Bay_ML, 
                    'MSE_add'= MSE_add, 'SE_beta_ML' = SE_beta_ML, 'SE_beta_Bay' = SE_beta_Bay,
                    'SE_beta_Bay_ML' = SE_beta_Bay_ML, 'SE_gamma' = SE_gamma)
  
  #print(Bay_list)
}


estimate_Bay_CV_SE_jackknife_individual <- function(data) {
  library(pracma)
  
  # Original estimator calculation (to get initial results)
  original_result <- estimate_Bay_CV(data)
  
  # Initialize parameters
  n <- data$n
  J <- data$k
  
  # Initialize vectors to store jackknife estimates for individual deletion
  jackknife_beta_b_Bay_individual <- numeric(n)
  jackknife_beta_b_Bay_ML_individual <- numeric(n)
  
  r = n # set up the number of replications for jackknife, r from 1 to n^J
  # Jackknife by deleting one individual from each group simultaneously
  for (i in 1:r) {
    # Generate J random indices (one for each group) from 1 to n
    random_indices <- sample(1:n, J, replace = TRUE)
    
    # Initialize a list to store indices to keep for each group
    indices_to_keep <- vector("list", J)
    
    for (j in 1:J) {
      indices_to_keep[[j]] <- (n * (j - 1) + 1):(n * j)
      indices_to_keep[[j]] <- indices_to_keep[[j]][-random_indices[j]]
    }
    
    # Create new indices after deleting one individual from each group
    indices <- unlist(indices_to_keep)
    
    # Create a new dataset excluding the selected individuals
    data_jackknife <- data
    data_jackknife$x <- data$x[indices]
    data_jackknife$y <- data$y[indices]
    if (data$kc>1){
      data_jackknife$C <- data$C[indices, , drop = FALSE]
    } else {
      data_jackknife$C <- data$C[indices, drop = FALSE]
    }
    data_jackknife$n <- n - 1
    data_jackknife$kn <- (n-1)*J
    
    # Recalculate estimators without the selected individuals
    result_jackknife <- estimate_Bay_CV(data_jackknife)
    
    # Store jackknife estimates
    jackknife_beta_b_Bay_individual[i] <- result_jackknife$beta_b_Bay
    jackknife_beta_b_Bay_ML_individual[i] <- result_jackknife$beta_b_Bay_ML
  }
  
  # Calculate jackknife means for individual deletion
  jackknife_mean_beta_b_Bay_individual <- mean(jackknife_beta_b_Bay_individual)
  jackknife_mean_beta_b_Bay_ML_individual <- mean(jackknife_beta_b_Bay_ML_individual)
  
  # Calculate jackknife standard errors for individual deletion
  SE_beta_Bay_jackknife_individual <- sqrt((n - 1) / n * sum((jackknife_beta_b_Bay_individual - jackknife_mean_beta_b_Bay_individual)^2))
  SE_beta_Bay_ML_jackknife_individual <- sqrt((n - 1) / n * sum((jackknife_beta_b_Bay_ML_individual - jackknife_mean_beta_b_Bay_ML_individual)^2))
  
  # Output updated results with jackknife standard errors
  list(
    beta_b_Bay = original_result$beta_b_Bay,
    beta_b_Bay_ML = original_result$beta_b_Bay_ML,
    SE_beta_Bay_jackknife_individual = SE_beta_Bay_jackknife_individual,
    SE_beta_Bay_ML_jackknife_individual = SE_beta_Bay_ML_jackknife_individual
    #original_result = original_result
  )
}

# Example usage:
# You would use your actual `estimate_Bay_CV` function here, which returns a complete `Bay` list.
# result <- mlob(y ~ X + Educ1 + Educ2, data = data, group = 15, alpha = 0.01, jackknife = FALSE)
# print(result)
# summary(result)
