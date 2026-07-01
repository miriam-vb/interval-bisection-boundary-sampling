# ----------------------------------------------------------------------------
# Interval Bisection Sampling of Bias Thresholds (2D)
# 
# This function uses interval bisection to approximate decision-invariant
# bias adjustment thresholds for network meta-analysis, with options to use 
# either preset or user-defined decision functions. Thresholds represent 
# the amount of simultaneous adjustment needed in up to two individual data 
# points before the treatment decision changes.
#
# @param data  Data frame containing the NMA data that is passed as an argument
#    to decision_function
# @param decision_function  Function accepting NMA data and bias adjustment
#    used to implement the decision rule at each step of the boundary finding 
#    method
# @param ind1  Numerical vector indicating the indices of the sequential list
#    of data points for which the first generic bias adjustment should be 
#    applied
# @param ind2  Numerical vector indicating the indices of the sequential list
#    of data points for which the second generic bias adjustment should be 
#    applied
# @param admin  Administrative cutoff value for bias adjustment beyond which 
#    decision invariance will not be assessed
# @param tol  Tolerance for the absolute difference between converging boundary
#    estimates
# @param rad_jump  Angle (in radians) by which the angle of bias assessment in 
#    the polar framework will increase for each sequential iteration 
# @param dist_tol  Euclidean distance tolerance between any two sequential points
# @param plot  Boolean determining whether the function call should also output
#    a plot of the invariant region
# @param preset  Numeric value determining whether a specific preset 
#    decision_function should be implemented rather than a user-supplied function
# @param parallel  Boolean determining whether to parallelize the threshold 
#    convergence method using all available cores (as opposed to sequential 
#    evaluation)
#
# @return  List containing thresh.df, a data frame of thresholds and new 
#    recommended treatments with columns \code{Bias_Ind_1}, \code{Bias_Ind_2}, 
#    and \code{New_Rec}, and args, a list of the arguments defined in the 
#    original function call
# ----------------------------------------------------------------------------

bias_thresh_2D <- function(data, decision_function, ind1, ind2, admin = 5, 
                          tol = 10**(-3), rad_jump = pi/90, dist_tol = 0.5, 
                          plot = TRUE, preset = 1, parallel = FALSE){
  
  # set decision_function to frequentist threshold analysis using the 
  # projection matrix with max efficacy as default
  if (preset == 1) {
    n <- length(data$vec)
    decision_function <- function(data, bias = rep(0,n)) {
      if ("C" %in% names(data)) {
        # allow for arm-level bias assessment using the matrix mapping arms to
        # contrasts
        H <- solve(t(data$X)%*%(data$W)%*%(data$X))%*%t(data$X)%*%(data$W)%*%
          (data$C)
      } else {
        H <- solve(t(data$X)%*%(data$W)%*%(data$X))%*%t(data$X)%*%(data$W)
      }
      
      trtm_estimates <- H%*%(data$vec+bias)
      if (all(trtm_estimates < 0)) {
        best <- c(0)
      } else {
        best <- c(which.max(trtm_estimates))
      }
      return(best)
    }
  }
  
  # output warning if an invalid preset was selected and no decision function
  # was supplied
  if (is.null(decision_function)) {
    stop("Invalid preset selected with no alternative decision function provided")
  }
  
  best <- decision_function(data)
  
  
  # function for evaluating Euclidean distance between vectors
  eucDist <- function(vec1, vec2) {
    return(sqrt(sum((vec1 - vec2)**2)))
  }
  
  library(DescTools)
  
  # implement random 2D bias adjustment at each selected data point and use
  # IVT and interval bisection method to compute thresholds
  library(doFuture)
  n_cores <- availableCores()
  
  # define function for implementing threshold convergence
  thresh_conv <- function(core, n_cores) {
    thresh.df <- data.frame(matrix(ncol=4,nrow=0))
    theta <- 2*pi*(core-1)/n_cores
    grain <- 0
    while (theta <= 2*pi*(core/n_cores)) {
      # ensure initial bias is greater than tol
      r <- runif(1,min = tol, max = 1)
      r0 <- 0
      r1 <- r
      r2 <- admin
      best0 <- best
      
      bvec <- rep(0,n)
      for (i in ind1) {
        bvec[i] <- PolToCart(r,theta)$x
      }
      for (i in ind2) {
        bvec[i] <- PolToCart(r,theta)$y
      }
      best1 <- decision_function(data, bias = bvec)
      
      for (i in ind1) {
        bvec[i] <- PolToCart(admin,theta)$x
      }
      for (i in ind2) {
        bvec[i] <- PolToCart(admin,theta)$y
      }
      best2 <- decision_function(data, bias = bvec)
      
      if (setequal(best2, best)) {
        # record administrative threshold if recommendation doesn't shift
        x <- PolToCart(admin,theta)$x
        y <- PolToCart(admin,theta)$y
        if (nrow(thresh.df) != 0) {
          if (eucDist(c(x,y),c(as.numeric(thresh.df[nrow(thresh.df),1]),
                               as.numeric(thresh.df[nrow(thresh.df),2]))) > dist_tol) {
            # if sequential points are not within dist_tol of each other by 
            # Euclidean distance, decrease rad_jump
            grain <- grain + 1
            theta <- theta - rad_jump/(2*grain)
            next
          }
        }
        trt <- "Admin"
      } else {
        # iterate until biases are within tolerance
        while(abs(r2 - r0) > tol) {
          # select interval [r0,r1] or [r1,r2], then obtain midpoint of 
          # chosen interval and update variables
          if (!setequal(best0,best1)) {
            min <- r0
            max <- r1
            mid <- min + (r1 - r0)/2
            r0 <- min
            r1 <- mid
            r2 <- max
          } else if (!setequal(best1,best2)) {
            min <- r1
            max <- r2
            mid <- min + (r2 - r1)/2
            r0 <- min
            r1 <- mid
            r2 <- max
          }
          # update treatment recommendations for each point of bias
          for (i in ind1) {
            bvec[i] <- PolToCart(r0,theta)$x
          }
          for (i in ind2) {
            bvec[i] <- PolToCart(r0,theta)$y
          }
          best0 <- decision_function(data, bias = bvec)
          
          for (i in ind1) {
            bvec[i] <- PolToCart(r1,theta)$x
          }
          for (i in ind2) {
            bvec[i] <- PolToCart(r1,theta)$y
          }
          best1 <- decision_function(data, bias = bvec)
          
          for (i in ind1) {
            bvec[i] <- PolToCart(r2,theta)$x
          }
          for (i in ind2) {
            bvec[i] <- PolToCart(r2,theta)$y
          }
          best2 <- decision_function(data, bias = bvec)
        }
        x <- PolToCart(r0,theta)$x
        y <- PolToCart(r0,theta)$y
        trt <- paste0(best2, collapse = " ")
        if (nrow(thresh.df) != 0) {
          if (eucDist(c(x,y),c(as.numeric(thresh.df[nrow(thresh.df),1]),
                               as.numeric(thresh.df[nrow(thresh.df),2]))) > dist_tol) {
            # if sequential points are not within dist_tol of each other by 
            # Euclidean distance, decrease rad_jump
            grain <- grain + 1
            theta <- theta - rad_jump/(2*grain)
            next
          }
        }
      }
      # if the last point has met dist_tol with the endpoint, cease iteration
      if (theta == 2*pi*(core/n_cores)) {
        break
      }
      
      # store bias threshold and admin indicator/new superior treatment
      # report the point just inside the invariant region
      row <- c(x,y,trt,theta)
      thresh.df <- rbind(thresh.df, row)
      
      # update theta
      theta <- theta + rad_jump
      grain <- 0
      
      # cap evaluation of dist_tol at end point
      if (theta >= 2*pi*(core/n_cores)) {
        theta <- 2*pi*(core/n_cores)
      }
    }
    colnames(thresh.df) <- c("Bias_Ind_1", "Bias_Ind_2", "New_Rec","Theta")
    return(thresh.df) 
  }
  
  # allow for parallelized option
  if (parallel) {
    plan(multisession)
    thresh <- foreach(core = 1:n_cores, .options.future = 
                        list(seed = TRUE)) %dofuture% {
                          thresh_conv(core, n_cores)
                        }
    # reform the data.frame using futures
    thresh.df <- Reduce(function(x, y) merge(x, y, all=TRUE), thresh)
    rownames(thresh.df) <- NULL
    
  } else {
    # evaluate once sequentially
    thresh.df <- thresh_conv(1, 1)  
  }
  
  # sort and convert bias columns to numeric data type
  thresh.df[, c(1,2,4)] <- apply(thresh.df[, c(1,2,4)], 2, 
                                 function(x) as.numeric(as.character(x)))
  thresh.df <- thresh.df[order(thresh.df$Theta),]
  thresh.df <- thresh.df[,1:3]
  
  
  # return bias values and index of treatment that became superior at each switch 
  # (or that it was admin cutoff, indicating a potential invariant region)
  thresh_obj <- list(thresh.df = thresh.df, args = list(data = data, 
                   decision_function = decision_function, ind1 = ind1, 
                   ind2 = ind2, admin = admin, tol = tol, rad_jump = rad_jump, 
                   dist_tol = dist_tol, plot = plot, preset = preset, 
                   parallel = parallel))
  
  if (plot == TRUE) {
    # print out a graph of the boundary points of the region
    print_thresh_2D(thresh_obj)
  }
  
  return(thresh_obj)
}


# ----------------------------------------------------------------------------
# Visualization of Bias Thresholds (2D)
# 
# This function allows for repeated plotting of the decision-invariant bias 
# adjustment thresholds and invariant region for two-dimensional threshold 
# analysis, and is called by \code{bias-thresh-2D} automatically.
#
# @param thresh_obj  List object obtained from \code{bias-thresh-2D} function 
#    call containing the estimated decision-invariant bias
#    adjustment thresholds and a list of the arguments defined in the 
#    original function call 
# ----------------------------------------------------------------------------

print_thresh_2D <- function(thresh_obj){
  # load packages
  library(ggplot2)
  library(ggforce)
  
  # import arguments from bias_thresh_2D call
  thresh.df <- thresh_obj$thresh.df
  ind1 <- thresh_obj$args$ind1
  ind2 <- thresh_obj$args$ind2
  admin <- thresh_obj$args$admin
  
  # define labels for sets or individual indices of bias adjustment
  biasX <- paste0("Bias (", ifelse(length(ind1) > 1, "Indices ", "Index "), 
                  paste0(ind1, collapse = ", "), ")")
  biasY <- paste0("Bias (", ifelse(length(ind2) > 1, "Indices ", "Index "), 
                  paste0(ind2, collapse = ", "), ")")
  
  plt <- ggplot(data = thresh.df, aes(x = Bias_Ind_1, y = Bias_Ind_2)) +
    geom_hline(yintercept = 0, linewidth = 0.2, color = "grey") +
    geom_vline(xintercept = 0, linewidth = 0.2, color = "grey") +
    
    # fill the approximate invariant region
    geom_polygon(alpha = 0.5, fill = "grey") +
    
    # colour the threshold in accordance with the recommendation shift
    geom_point(aes(colour = as.factor(New_Rec))) +
    
    # visualize the boundary of the administrative cutoff
    geom_circle(data = data.frame(null = c(0)), aes(x0=0,y0=0,r=admin), 
                inherit.aes=FALSE, linetype=2) +
    
    # label the plot in accordance with the settings of the boundary finding
    # function call
    labs(x = biasX, y=biasY, color = "New Recommendation") +
    coord_fixed() +
    theme_classic()
  
  # print ggplot object
  print(plt)
}

