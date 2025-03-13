#' Maximum Pseudolikelihood estimation for a Markov Random Field model for
#' ordinal variables with the possibility to supply the adjacency matrix.
#'
#' The function \code{mple} estimates the parameters for the ordinal MRF
#' by optimizing the joint pseudolikelihood with Newton-Raphson.
#'
#' @param x An \code{no_persons} by \code{no_nodes} matrix containing the
#' categories coded as non-negative integers (i.e., coded
#' \code{0, 1, ..., no_categories}) for \code{no_persons} independent
#' observations on \code{no_nodes} variables in the network or graph.
#'
#' @param no_categories The maximum category.
#'
# @param precision A number between zero and one. The prior precision that is
# desired for edge selection. Equal to one minus the desired type-1 error.
# Defaults to \code{.975}.
#'
#' @param convergence_criterion The convergence criterion for the
#' pseudoposterior values in the EM algorithm. Defaults to
#' \code{sqrt(.Machine$double.eps)}.
#'
#' @param maximum_iterations The maximum number of EM iterations used. Defaults
#' to \code{1e3}. A warning is issued if procedure has not converged in
#' \code{maximum_iterations} iterations.
#'
#' @param thresholds A \code{no_nodes} by \code{no_categories} matrix
#' \code{thresholds}. Used as starting values in the Newton-Raphson procedure.
#' Optional.
#'
#' @param interactions A \code{no_nodes} by \code{no_nodes} matrices
#' \code{interactions}. Used as starting values in the Newton-Raphson procedure.
#' Optional.
#'
#' @param G A \code{no_nodes} by \code{no_nodes} adjacency matrix.
#'
#' @return A list containing the \code{no_nodes} by \code{no_nodes} matrices
#' \code{interactions} and the \code{no_nodes} by
#' \code{no_categories} matrix \code{thresholds}. The matrix \code{interactions}
#' is a numeric matrix which contains the maximum pseudolikelihood estimates of
#' the pairwise association parameters. The matrix \code{thresholds} contains
#' the maximum pseudolikelihood estimates of the category thresholds parameters.
#' @export

mple_G = function(x,
                  convergence_criterion = sqrt(.Machine$double.eps),
                  maximum_iterations = 1e3,
                  thresholds,
                  interactions,
                  G) {
  # Check data input ------------------------------------------------------------
  if(!inherits(x, what = "matrix") && !inherits(x, what = "data.frame"))
    stop("The input x needs to be a matrix or dataframe.")
  if(inherits(x, what = "data.frame"))
    x = data.matrix(x)
  if(ncol(x) < 2)
    stop("The matrix x should have more than one variable (columns).")
  if(nrow(x) < 2)
    stop("The matrix x should have more than one observation (rows).")

  # Format the data input -------------------------------------------------------
  data = reformat_data(x = x)
  x = data$x
  no_categories = data$no_categories
  no_nodes = ncol(x)
  no_interactions = no_nodes * (no_nodes - 1) / 2
  no_thresholds = sum(no_categories)
  no_parameters = no_thresholds + no_interactions

  # Check NR input --------------------------------------------------------------
  if(convergence_criterion <= 0)
    stop("Parameter ``convergence_criterion'' needs to be positive.")
  if(maximum_iterations <= 0 ||
     abs(maximum_iterations - round(maximum_iterations)) >
     sqrt(.Machine$double.eps))
    stop("Parameter ``maximum_iterations'' needs to be a positive integer.")

  # Starting values -----------------------------------------------------------
  if(!hasArg("thresholds")) {
    thresholds = matrix(0,
                        nrow = no_nodes,
                        ncol = max(no_categories))
  }
  if(!hasArg("interactions")) {
    interactions = matrix(0,
                          nrow = no_nodes,
                          ncol = no_nodes)
  }

  # UPDATE: Create the mask vector v based on G  -------------------------------
  v = c(rep(1, no_thresholds), G[lower.tri(G)])

  # Newton-Raphson ------------------------------------------------------------
  log_pl = log_pseudolikelihood(interactions * G,
                                thresholds,
                                observations = x,
                                no_categories)

  hessian = matrix(data = NA,
                   nrow = no_parameters,
                   ncol = no_parameters)
  gradient = matrix(data = NA,
                    nrow = 1,
                    ncol = no_parameters)

  for(iteration in 1:maximum_iterations) {
    old_log_pl = log_pl

    # Compute gradient vector (first order derivatives) ------------------------
    gradient[1:no_thresholds] =
      gradient_thresholds_pseudolikelihood(interactions = G * interactions,
                                           thresholds = thresholds,
                                           observations = x,
                                           no_categories)
    gradient[-c(1:no_thresholds)] =
      gradient_interactions_pseudolikelihood(interactions = G * interactions,
                                             thresholds = thresholds,
                                             observations = x,
                                             no_categories)

    # UPDATE: Apply the adjacency mask to the gradient -------------------------
    gradient = gradient * matrix(v, nrow = 1)

    # Compute Hessian matrix (second order partial derivatives) ---------------
    hessian[1:no_thresholds, 1:no_thresholds] =
      hessian_thresholds_pseudolikelihood(interactions = G * interactions,
                                          thresholds = thresholds,
                                          observations = x,
                                          no_categories)

    hessian[-(1:no_thresholds), -(1:no_thresholds)] =
      hessian_interactions_pseudolikelihood(interactions = G * interactions,
                                            thresholds = thresholds,
                                            observations = x,
                                            no_categories)

    hessian[-(1:no_thresholds), 1:no_thresholds] =
      hessian_crossparameters(interactions = G * interactions,
                              thresholds = thresholds,
                              observations = x,
                              no_categories)

    hessian[1:no_thresholds, -(1:no_thresholds)] =
      t(hessian[-(1:no_thresholds), 1:no_thresholds])

    # Apply the adjacency mask to the Hessian
    hessian = hessian * (matrix(v, ncol = 1) %*% matrix(v, nrow = 1))

    # Update parameter values (Newton-Raphson step) ---------------------------
    Delta = gradient %*% solve(hessian + diag(1e-6, nrow(hessian))) # issue arises here!!!!!
    if(any(is.nan(Delta)) || any(is.infinite(Delta)))
      stop("log_pseudolikelihood optimization failed. Please check the data. If the data checks out, please try different starting values.")

    cntr = 0
    for(node in 1:no_nodes) {
      for(category in 1:no_categories[node]) {
        cntr = cntr + 1
        thresholds[node, category] = thresholds[node, category] - Delta[cntr]
      }
    }
    for(node in 1:(no_nodes - 1)) {
      for(node_2 in (node + 1):no_nodes) {
        cntr = cntr + 1
        if(G[node, node_2] == 1) {  # Ensure updates only for allowed interactions
          interactions[node, node_2] = interactions[node, node_2] - Delta[cntr]
          interactions[node_2, node] = interactions[node, node_2]
        }
      }
    }

    # Check for convergence ---------------------------------------------------
    log_pl = log_pseudolikelihood(interactions * G,
                                  thresholds,
                                  observations = x,
                                  no_categories)

    if(abs(log_pl - old_log_pl) < convergence_criterion)
      break
  }

  if(abs(log_pl - old_log_pl) >= convergence_criterion &&
     iteration == maximum_iterations)
    warning(paste("The optimization procedure did not converge in",
                  maximum_iterations, "iterations.",
                  sep = " "),
            call. = FALSE)

  # Preparing the output --------------------------------------------------------
  if(is.null(colnames(x))){
    data_columnnames = paste0("node ", 1:no_nodes)
    colnames(interactions) = data_columnnames
    rownames(interactions) = data_columnnames
    rownames(thresholds) = data_columnnames
  } else {
    data_columnnames <- colnames(x)
    colnames(interactions) = data_columnnames
    rownames(interactions) = data_columnnames
    rownames(thresholds) = data_columnnames
  }
  colnames(thresholds) = paste0("category ", 1:max(no_categories))

  return(list(interactions = interactions, thresholds = thresholds))
}
