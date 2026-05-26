# Example code to run the DT trans-prob function. 
# Following contains bunch of functions. Example is on bottom of code-file. Hope it works as desired :) 


### --- Functions
# 1) Function to call within code. Note good name, but used (initially) to plot the trans probs. Hence 'get_plot_df':
get_plot_df = function(ndf, mdl, transitions, ci_hazard = FALSE, ci_trans = FALSE, nsim = 100, alpha = 0.05, obs_start_day = 4){
    get_named_transition_matrix = function(transitions){
        # Function used to create identity matrix on which transitions probabilities are based on (see diag_transition_matrix in add_trans_probs)
        transitions_from = sub("->(.*)", "", transitions) |> unique()
        transitions_to = sub("(.*)->", "", transitions) |> unique()
        observed_states = c(transitions_from, transitions_to) |> unique() |> sort()
        transition_matrix = diag(1, nrow = length(observed_states))
        colnames(transition_matrix) = rownames(transition_matrix) = observed_states

        return(transition_matrix)
    }

  # Define identity matrix for trans_probs (see diag_transition_matrix in add_trans_probs)
  diag_transition_matrix = get_named_transition_matrix(transitions)

  # Use provided model to predict hazard
  ndf$hazard = predict(mdl, type="response", newdata = ndf)
  ndf = ndf |> dplyr::mutate(from = sub("->(.*)", "", transition), to = sub("(.*)->", "", transition)) 

  # Group ndf. So far, I hard coded this. Great Code, much wow. 
  ndf = ndf |> 
    dplyr::group_by(transition) |>
    dplyr::arrange(transition, tend) 

  # Add trans probs and re-assign trans prob to df
  ndf = ndf |>
    add_trans_prob(diag_transition_matrix)

  if (ci_hazard){
    ndf = ndf |> 
        # NOTE: I cleaned the code and originally, add_trans_prob did not remove the grouping. Now it does :D Sorry for that, but towards the end of this file, I pasted all the original code - which is ugly
        dplyr::group_by(transition) |>
        dplyr::arrange(transition, tend) |> 
        add_hazard_ci(object = mdl, nsim = nsim, alpha = alpha)
  }
  if (ci_trans) {
    ndf = ndf |> 
        dplyr::group_by(transition) |>
        dplyr::arrange(transition, tend) |> 
        add_trans_ci(diag_transition_matrix = diag_transition_matrix, object = mdl, nsim = nsim, alpha = alpha)
  }

  return(ndf |> dplyr::ungroup())
}

add_trans_prob = function(data, diag_transition_matrix){
    # Add trans_prob NA here to ensure I misstakes if they happen --> trans_prob = NA --> something went wrong
    data = data |> dplyr::ungroup(transition) |> dplyr::mutate(trans_prob = NA)

    # List of hazards matrices at each time for each grouping factor
    # i.e. for every grouping factor, we obtain t_max number of hazard matrices, e.g.:
        # [[1]][[56]]
        #           1          2         3          4            5
        # 1 0.9748542 0.00000000 0.0000000 0.01962785 5.517984e-03
        # 2 0.1089366 0.89106325 0.0000000 0.00000000 1.686735e-07
        # 3 0.0000000 0.02118114 0.9657027 0.00000000 1.311617e-02
        # 4 0.0000000 0.00000000 0.0000000 1.00000000 0.000000e+00
        # 5 0.0000000 0.00000000 0.0000000 0.00000000 1.000000e+00

        # [[1]][[57]]
        #           1          2         3         4            5
        # 1 0.9752644 0.00000000 0.0000000 0.0193041 5.431470e-03
        # 2 0.1073383 0.89266157 0.0000000 0.0000000 1.354789e-07
        # 3 0.0000000 0.02064386 0.9662805 0.0000000 1.307566e-02
        # 4 0.0000000 0.00000000 0.0000000 1.0000000 0.000000e+00
        # 5 0.0000000 0.00000000 0.0000000 0.0000000 1.000000e+00
    # The indices of the list are: [[grouping_factor_index]][[time_point]]
    hazard_matrix_list =  data |> 
        dplyr::group_map(.f = ~ get_transition_matrices(data = .x, diag_transition_matrix = diag_transition_matrix, tend_min = min(data$tend), tend_max = max(data$tend)))

    # For every grouping_factor, build the matrix product up to time point t
    for (group_idx in unique(dplyr::group_indices(data))) {
        # Get current list of hazard matrices:
        hazard_matrices = hazard_matrix_list[[group_idx]]
        # Reduce is like comprod, but for matrices. Think its pretty fast, too: 
        trans_probs = Reduce("%*%", hazard_matrices, accumulate = TRUE)
        # Replace NAs in trans_prob column where possible (i.e. for current grouping)
        data[data |> dplyr::group_indices() == group_idx, ] = data[data |> dplyr::group_indices() == group_idx, ] |> 
        dplyr::rowwise() |> 
        dplyr::mutate(
            trans_prob = trans_probs[[tend]][from, to]
        ) 
    }

    return(data)
}

get_transition_matrices = function(data, diag_transition_matrix, tend_min, tend_max, grp = NA){
    # Called within add_trans_prob. Given a start and an end time point, return all hazard matrices for these time points:
    lapply(
        tend_min:tend_max, 
        function(t) get_transition_matrix(t = t, data = data, diag_transition_matrix = diag_transition_matrix, grp = grp)
    )
}

get_transition_matrix = function(t, data, diag_transition_matrix, grp){
    # Get the hazard matrix for a current time point 
    
    # Grab currently relevant hazards:
    hazards_df = data |> 
        dplyr::filter(tend == t) |> 
        dplyr::select(from, to, hazard) |> 
        unique()

    # Create matrix of 0s to plug the hazards into. "Inherit" from transition matrix
    # This is basically 1 - 1, but in desired dimensions because diag_transition_matrix has correct dimensions. 
    # I think it wouldve been a hussle to get the right dimensions in this function-context, but I do not fully remember why I did this instead of a new matrix. 
    hazard_matrix = diag_transition_matrix - diag_transition_matrix 

    # Shit code to replace hazards in matrix: 
    for (row_idx in 1:nrow(hazards_df)){
        from = hazards_df[row_idx, "from"] |> unlist(use.names = FALSE)
        to = hazards_df[row_idx, "to"] |> unlist(use.names = FALSE)
        hazard = hazards_df[row_idx, "hazard"] |> unlist(use.names = FALSE)

        hazard_matrix[from, to] = hazard
    }

    # 1 - hazard_rowsums 
    hazards_t = diag(1 - rowSums(hazard_matrix)) + hazard_matrix
    
    return(hazards_t)
}

add_hazard_ci = function(newdata, object, nsim, alpha, obs_start_day = 4){
  # Basically copy pasta from what u've done
  X     <- mgcv::predict.gam(object, newdata = newdata, type = "lpmatrix")
  coefs <- coef(object)
  V     <- object$Vp

  # define groups: 1. all grouping variables -> cumu hazards, 2. all but transition -> trans_prob
  groups_array <- dplyr::group_indices(newdata)
  groups_trans <- newdata |> dplyr::ungroup(transition) |> dplyr::group_indices()

  # Use object to identify used family which contains the link / respose function
  # NOTE: If you use the object here, any link function is possible, not just poisson
  fam = object$family
  sim_coef_mat <- mvtnorm::rmvnorm(nsim, mean = coefs, sigma = V)
  sim_fit_mat <- apply(sim_coef_mat, 1, function(z) fam$linkinv(X %*% z))

  # create list with replicated newdata
  nlst <- as.list(replicate(nsim, newdata[,c("tend", "transition", "from")], simplify=F))

  # add cumu-hazard in each element and calculate trans_prob with perturbed hazards
  nlst <- lapply(1:nsim, function(i) {
    nlst[[i]] <- cbind(nlst[[i]], hazard = sim_fit_mat[, i]) # add hazard

    nlst[[i]] = nlst[[i]] |> dplyr::mutate(from = sub("->(.*)", "", transition), to = sub("(.*)->", "", transition)) 

    nlst[[i]]
  })

  sim_haz <- do.call(cbind, lapply(nlst, function(df) df$hazard))
  newdata$hazard_lower <- apply(sim_haz, 1, quantile, probs = alpha / 2)
  newdata$hazard_upper <- apply(sim_haz, 1, quantile, probs = 1 - alpha / 2)

  newdata

}


### Add CIs for transition probabilitiy
# Based on pammtools: https://github.com/adibender/pammtools/blob/master/R/add-functions.R
add_trans_ci <- function(newdata, object, diag_transition_matrix, nsim, alpha, obs_start_day = 4){
  newdata = newdata |> 
    dplyr::mutate(from = sub("->(.*)", "", transition), to = sub("(.*)->", "", transition)) 

  X     <- mgcv::predict.gam(object, newdata = newdata, type = "lpmatrix")
  coefs <- coef(object)
  V     <- object$Vp

  # define groups: 1. all grouping variables -> cumu hazards, 2. all but transition -> trans_prob
  groups_array <- dplyr::group_indices(newdata)
  groups_trans <- newdata |> dplyr::ungroup(transition) |> dplyr::group_indices()

  # Use object to identify used family which contains the link / respose function
  fam = object$family
  sim_coef_mat <- mvtnorm::rmvnorm(nsim, mean = coefs, sigma = V)
  sim_fit_mat <- apply(sim_coef_mat, 1, function(z) fam$linkinv(X %*% z))

  # create list with replicated newdata
  nlst <- as.list(replicate(nsim, newdata[,c("tend", "transition", "completed_mv_duration", "from")], simplify=F))

  # add cumu-hazard in each element and calculate trans_prob with perturbed hazards
  nlst <- lapply(1:nsim, function(i) {
    nlst[[i]] <- cbind(nlst[[i]], hazard = sim_fit_mat[, i]) # add hazard

    nlst[[i]] = nlst[[i]] |> dplyr::mutate(from = sub("->(.*)", "", transition), to = sub("(.*)->", "", transition)) 

    # Get transition probabilities based on new hazards
    nlst[[i]] <- split(nlst[[i]], groups_trans) |> purrr::map_dfr(add_trans_prob, diag_transition_matrix = diag_transition_matrix)

    nlst[[i]]
  })

  sim_trans_probs <- do.call(cbind, lapply(nlst, function(df) df$trans_prob))
  newdata$trans_lower <- apply(sim_trans_probs, 1, quantile, probs = alpha / 2)
  newdata$trans_upper <- apply(sim_trans_probs, 1, quantile, probs = 1 - alpha / 2)

  newdata
}


### --- Code Example
# (Use the .Rds-object added to the mail. It contains ped data frame for first 100 IDs in original data.)
# Read data
ped = readRDS("test_ped.Rds")
# Fit model - should take ~1sec
mdl <- mgcv::bam(
    ped_status ~ transition + s(tend, by = transition),
    data = ped,
    family   = binomial(),          
    discrete = TRUE,
    method   = "fREML"
)
# Create new df
ndf <- pammtools::make_newdata(
  ped,
  transition = unique(transition),
  tend = unique(tend),
) 
# Get transition probabilities
plot_df = get_plot_df(
  ndf = ndf, 
  mdl = mdl, 
  transitions = unique(ped$transition),
  ci_hazard = TRUE, 
  ci_trans = TRUE,
  nsim = 5, 
  alpha = 0.05
)

plot_df |> 
    ggplot(aes(x = tend, y = hazard)) + 
    geom_ribbon(aes(ymax = hazard_upper, ymin = hazard_lower), xmin = 0, xmax = Inf, alpha = 0.1, fill = "blue") +
    geom_line() + 
    facet_wrap(~transition)

plot_df |> 
    ggplot(aes(x = tend, y = trans_prob)) + 
    geom_ribbon(aes(ymax = trans_upper, ymin = trans_lower), xmin = 0, xmax = Inf, alpha = 0.1, fill = "blue") +
    geom_line() + 
    facet_wrap(~transition)






# # Here you go fill all the original code - have fun...
# enforce_deterministic_mvyes_hazards = function(ndf, obs_start_day, mvyes_state_idx = "3"){
#     grab_relevant_hazards = function(data){
#         data |> dplyr::filter(from != mvyes_state_idx, completed_mv_duration <= tend) |> dplyr::pull(hazard) 
#     }
#     # If completed_mv_duration is used as covariate in the model, the inital patient jounrey ins acutally deterministic
#     # i.e. we enforce everyone with a given mv_duration to stay exactly that duration in mv:yes
#     # Thus: 
#         # - set hazards to 0 while patient has to remain in mv:yes
#         # - set hazards to 1 while patient left mv:yes 
#         # - use conditional cause specific hazards [conditioned on "any cause"] when tend == completed_mv_duration 
#             # --> i.e. divide cause specific hazards by all cause hazard at t

#     data_grouping = dplyr::group_vars(ndf)
#     ndf = ndf |> dplyr::ungroup()

#     init_ndf = ndf
#     # Save the currently calculated hazards to plot them as sanity display
#     ndf$ghost_hazard = ndf$hazard
#     # Hazard to 0 where patient must still remain in mv:yes bc completed_mv_duration is larger than time spent in state
#     # ndf[ndf$completed_mv_duration > ndf$tend + (obs_start_day - 1) & ndf$from == mvyes_state_idx, "hazard"] = 0
#     ndf[ndf$completed_mv_duration > ndf$tend + (obs_start_day - 1), "hazard"] = 0
#     # Hazard to 1 as survival in state must be set to 0
#     ndf[(ndf$completed_mv_duration < ndf$tend + (obs_start_day - 1) & ndf$from == mvyes_state_idx), "hazard"] = 0
#     # # Use conditional cause specific hazard to identify the transitions
#     ndf[
#         (ndf$completed_mv_duration == ndf$tend + (obs_start_day - 1) & ndf$from == mvyes_state_idx), 
#         "hazard"
#     ] = ndf |> 
#         dplyr::filter(completed_mv_duration == tend + (obs_start_day - 1), from == mvyes_state_idx) |> 
#         dplyr::select(transition, from, completed_mv_duration, hazard) |> 
#         dplyr::mutate(.by = c(from, completed_mv_duration), all_cause_hazard = sum(hazard)) |> 
#         dplyr::mutate(cond_hazard = hazard / all_cause_hazard) |> 
#         dplyr::pull(cond_hazard)

#     # Check that I only adjusted the hazards leaving mv:yes
#     other_hazards_init_value = init_ndf |> grab_relevant_hazards()
#     other_hazards_post_value = ndf |> grab_relevant_hazards()
#     other_hazards_unchanged = other_hazards_init_value == other_hazards_post_value

#     if (!all(other_hazards_unchanged)){
#         browser()
#         # Affected hazards
#         before = init_ndf |> dplyr::filter(from != mvyes_state_idx)
#         after = ndf |> dplyr::filter(from != mvyes_state_idx) 
#         before[!other_hazards_unchanged, ] |> dplyr::select(transition, tend, completed_mv_duration)
#     }
#     ndf = ndf |> dplyr::mutate(ghost_hazard = dplyr::if_else(from == mvyes_state_idx, hazard, ghost_hazard))
#     return(ndf |> dplyr::group_by(dplyr::across(dplyr::all_of(data_grouping))))
# }

# # Add hazard, cumu hazard and trans probs
# get_plot_df = function(ndf, mdl, transitions, ci_hazard = FALSE, ci_trans = FALSE, nsim = 100, alpha = 0.05, obs_start_day = 4){
#   # Allow a Flag to not adjust the hazard, i.e. an unconditional view on completed_mv_duration
#   # This is done using the flag -1 in completed_mv_duration: 
#   ignore_completed_mv_duration = all(unique(ndf$completed_mv_duration) == -1)
#   # Get Identity matrix with column / row names corresponding to all observed transitions
#   diag_transition_matrix = get_named_transition_matrix(transitions)

#   ndf$hazard = predict(mdl, type="response", newdata = ndf)
#   ndf = ndf |> dplyr::mutate(from = sub("->(.*)", "", transition), to = sub("(.*)->", "", transition)) 
#   if (!ignore_completed_mv_duration) {
#     ndf = enforce_deterministic_mvyes_hazards(ndf, obs_start_day = obs_start_day)
#   }

#   ndf = ndf |> 
#     dplyr::group_by(
#         age, apache_score, bmi, year, primary_diagnosis, admission_type, sex,  
#         n_days_oral_nutrition_d1_4, n_days_low_protein_intake_d1_4, n_days_low_calorie_intake_d1_4, n_days_propofol_d1_4, 
#         completed_mv_duration, transition
#     ) |>
#     dplyr::arrange(
#         age, apache_score, bmi, year, year, primary_diagnosis, admission_type, sex,  
#         n_days_oral_nutrition_d1_4, n_days_low_protein_intake_d1_4, n_days_low_calorie_intake_d1_4, n_days_propofol_d1_4, 
#         completed_mv_duration, transition, tend
#     ) 
#   ndf = ndf |>
#     # Dont need it
#     # dplyr::mutate(
#     #   cumu_hazard = cumsum(hazard)
#     # ) |> 
#     add_trans_prob(diag_transition_matrix)

#   if (ci_hazard){
#     ndf = ndf |> 
#         add_hazard_ci(object = mdl, nsim = nsim, alpha = alpha)
#   }
#   if (ci_trans) {
#     ndf = ndf |> 
#         add_trans_ci(diag_transition_matrix = diag_transition_matrix, object = mdl, nsim = nsim, alpha = alpha)
#   }

#   return(ndf |> dplyr::ungroup())
# }

# add_hazard_ci = function(newdata, object, nsim, alpha, obs_start_day = 4){
#   ignore_completed_mv_duration = all(unique(newdata$completed_mv_duration) == -1)

#   X     <- mgcv::predict.gam(object, newdata = newdata, type = "lpmatrix")
#   coefs <- coef(object)
#   V     <- object$Vp

#   # define groups: 1. all grouping variables -> cumu hazards, 2. all but transition -> trans_prob
#   groups_array <- dplyr::group_indices(newdata)
#   groups_trans <- newdata |> dplyr::ungroup(transition) |> dplyr::group_indices()

#   # Use object to identify used family which contains the link / respose function
#   fam = object$family
#   sim_coef_mat <- mvtnorm::rmvnorm(nsim, mean = coefs, sigma = V)
#   sim_fit_mat <- apply(sim_coef_mat, 1, function(z) fam$linkinv(X %*% z))

#   # create list with replicated newdata
#   nlst <- as.list(replicate(nsim, newdata[,c("tend", "transition", "completed_mv_duration", "from")], simplify=F))

#   # add cumu-hazard in each element and calculate trans_prob with perturbed hazards
#   nlst <- lapply(1:nsim, function(i) {
#     nlst[[i]] <- cbind(nlst[[i]], hazard = sim_fit_mat[, i]) # add hazard
#     # Set the hazards where completed_mv_duration > tend to 0! 
#         # bzw. we start observations LATER than day 1 --> on time 0, the mv duration may be alrger than 0

#     if (!ignore_completed_mv_duration) {
#         nlst[[i]] = enforce_deterministic_mvyes_hazards(nlst[[i]], obs_start_day = obs_start_day)
#     }
#     # nlst[[i]][nlst[[i]]$completed_mv_duration > nlst[[i]]$tend + (obs_start_day - 1), "hazard"] = 0
#     nlst[[i]] = nlst[[i]] |> dplyr::mutate(from = sub("->(.*)", "", transition), to = sub("(.*)->", "", transition)) 

#     nlst[[i]]
#   })

#   sim_haz <- do.call(cbind, lapply(nlst, function(df) df$hazard))
#   newdata$hazard_lower <- apply(sim_haz, 1, quantile, probs = alpha / 2)
#   newdata$hazard_upper <- apply(sim_haz, 1, quantile, probs = 1 - alpha / 2)

#   newdata

# }


#   add_trans_prob = function(data, diag_transition_matrix){
#     data = data |> dplyr::ungroup(transition) |> dplyr::mutate(trans_prob = NA)
#     # List of hazards matrices at each time for each grouping factor
#     hazard_matrix_list =  data |> 
#       dplyr::group_map(.f = ~ get_transition_matrices(data = .x, diag_transition_matrix = diag_transition_matrix, tend_min = min(data$tend), tend_max = max(data$tend)))

#     for (group_idx in unique(dplyr::group_indices(data))) {
#       hazard_matrices = hazard_matrix_list[[group_idx]]
#       trans_probs = Reduce("%*%", hazard_matrices, accumulate = TRUE)
#       data[data |> dplyr::group_indices() == group_idx, ] = data[data |> dplyr::group_indices() == group_idx, ] |> 
#         dplyr::rowwise() |> 
#         dplyr::mutate(
#           trans_prob = trans_probs[[tend]][from, to]
#         ) 
#     }

#     return(data)
#   }

#   get_transition_matrices = function(data, diag_transition_matrix, tend_min, tend_max, grp = NA){
#     lapply(
#       tend_min:tend_max, 
#       function(t) get_transition_matrix(t = t, data = data, diag_transition_matrix = diag_transition_matrix, grp = grp)
#     )
#   }

#   get_transition_matrix = function(t, data, diag_transition_matrix, grp){
#     hazards_df = data |> 
#       dplyr::filter(tend == t) |> 
#       dplyr::select(from, to, hazard) |> 
#       unique()
    
#     # Create matrix of 0s to plug the hazards into. "Inherit" from transition matrix
#     hazard_matrix = diag_transition_matrix - diag_transition_matrix 

#     for (row_idx in 1:nrow(hazards_df)){
#       from = hazards_df[row_idx, "from"] |> unlist(use.names = FALSE)
#       to = hazards_df[row_idx, "to"] |> unlist(use.names = FALSE)
#       hazard = hazards_df[row_idx, "hazard"] |> unlist(use.names = FALSE)

#       hazard_matrix[from, to] = hazard
#     }

#     hazards_t = diag(1 - rowSums(hazard_matrix)) + hazard_matrix
#     return(hazards_t)
#   }

#   get_named_transition_matrix = function(transitions){
#     transitions_from = sub("->(.*)", "", transitions) |> unique()
#     transitions_to = sub("(.*)->", "", transitions) |> unique()
#     observed_states = c(transitions_from, transitions_to) |> unique() |> sort()
#     transition_matrix = diag(1, nrow = length(observed_states))
#     colnames(transition_matrix) = rownames(transition_matrix) = observed_states

#     return(transition_matrix)
#   }



# ### Add CIs for transition probabilitiy
# # Based on pammtools: https://github.com/adibender/pammtools/blob/master/R/add-functions.R
# add_trans_ci <- function(newdata, object, diag_transition_matrix, nsim, alpha, obs_start_day = 4){
#   ignore_completed_mv_duration = all(unique(newdata$completed_mv_duration) == -1)

#   newdata = newdata |> 
#     dplyr::mutate(from = sub("->(.*)", "", transition), to = sub("(.*)->", "", transition)) 
#     # dplyr::group_by(year, primary_diagnosis, admission_type, sex,  n_days_oral_nutrition_d1_4, n_days_low_protein_intake_d1_4, n_days_low_calorie_intake_d1_4, n_days_propofol_d1_4, completed_mv_duration, transition) |>
#     # dplyr::arrange(year, primary_diagnosis, admission_type, sex,  n_days_oral_nutrition_d1_4, n_days_low_protein_intake_d1_4, n_days_low_calorie_intake_d1_4, n_days_propofol_d1_4, completed_mv_duration, transition, tend) 

#   X     <- mgcv::predict.gam(object, newdata = newdata, type = "lpmatrix")
#   coefs <- coef(object)
#   V     <- object$Vp

#   # define groups: 1. all grouping variables -> cumu hazards, 2. all but transition -> trans_prob
#   groups_array <- dplyr::group_indices(newdata)
#   groups_trans <- newdata |> dplyr::ungroup(transition) |> dplyr::group_indices()

#   # Use object to identify used family which contains the link / respose function
#   fam = object$family
#   sim_coef_mat <- mvtnorm::rmvnorm(nsim, mean = coefs, sigma = V)
#   sim_fit_mat <- apply(sim_coef_mat, 1, function(z) fam$linkinv(X %*% z))

#   # create list with replicated newdata
#   nlst <- as.list(replicate(nsim, newdata[,c("tend", "transition", "completed_mv_duration", "from")], simplify=F))

#   # add cumu-hazard in each element and calculate trans_prob with perturbed hazards
#   nlst <- lapply(1:nsim, function(i) {
#     nlst[[i]] <- cbind(nlst[[i]], hazard = sim_fit_mat[, i]) # add hazard
#     # Set the hazards where completed_mv_duration > tend to 0! 
#         # bzw. we start observations LATER than day 1 --> on time 0, the mv duration may be alrger than 0

#     if (!ignore_completed_mv_duration) {
#         nlst[[i]] = enforce_deterministic_mvyes_hazards(nlst[[i]], obs_start_day = obs_start_day)
#     }
#     # nlst[[i]][nlst[[i]]$completed_mv_duration > nlst[[i]]$tend + (obs_start_day - 1), "hazard"] = 0
#     nlst[[i]] = nlst[[i]] |> dplyr::mutate(from = sub("->(.*)", "", transition), to = sub("(.*)->", "", transition)) 

#     # Get transition probabilities based on hazards
#     nlst[[i]] <- split(nlst[[i]], groups_trans) |> purrr::map_dfr(add_trans_prob, diag_transition_matrix = diag_transition_matrix)

#     nlst[[i]]
#   })

#   sim_trans_probs <- do.call(cbind, lapply(nlst, function(df) df$trans_prob))
#   newdata$trans_lower <- apply(sim_trans_probs, 1, quantile, probs = alpha / 2)
#   newdata$trans_upper <- apply(sim_trans_probs, 1, quantile, probs = 1 - alpha / 2)

#   newdata
# }


# get_hazard_ratios = function(ndf, ref_model, cmp_model, ci = FALSE, nsim = 100, alpha = 0.05, pal_bias_days = 2, obs_start_day = 4){
#   get_hazard_ci = function(newdata, ref_model, cmp_model, nsim, alpha){
#       get_sim_fit_mat = function(newdata, object, nsim){
#         X     <- mgcv::predict.gam(object, newdata = newdata, type = "lpmatrix")
#         coefs <- coef(object)
#         V     <- object$Vp

#         # Use object to identify used family which contains the link / respose function
#         fam = object$family
#         sim_coef_mat <- mvtnorm::rmvnorm(nsim, mean = coefs, sigma = V)
#         sim_fit_mat <- apply(sim_coef_mat, 1, function(z) fam$linkinv(X %*% z))

#         return(sim_fit_mat)
#       }
#       ref_sim_fit_mat = get_sim_fit_mat(newdata, ref_model, nsim)
#       cmp_sim_fit_mat = get_sim_fit_mat(newdata, cmp_model, nsim)

#       # create list with replicated newdata
#       nlst <- as.list(replicate(nsim, newdata[,c("tend", "transition", "completed_mv_duration", "from", "to")], simplify=F))

#       # add cumu-hazard in each element and calculate trans_prob with perturbed hazards
#       nlst <- lapply(1:nsim, function(i) {
#         ref_sim_ndf = cbind(nlst[[i]], hazard = ref_sim_fit_mat[, i]) 
#         cmp_sim_ndf = cbind(nlst[[i]], hazard = cmp_sim_fit_mat[, i]) 

#         ref_sim_ndf  = enforce_deterministic_mvyes_hazards(ref_sim_ndf, obs_start_day = obs_start_day)
#         cmp_sim_ndf  = enforce_deterministic_mvyes_hazards(cmp_sim_ndf, obs_start_day = obs_start_day)
#         # For the palliative bias model, additionally enforce hazard of 0 the first two days someone is discharged from mv:yes
#         cmp_sim_ndf[cmp_sim_ndf$from == "2" & cmp_sim_ndf$to == "5" & 2 + cmp_sim_ndf$completed_mv_duration >= cmp_sim_ndf$tend + (obs_start_day - 1), "hazard"] = 0
#         cmp_sim_ndf[cmp_sim_ndf$from == "1" & cmp_sim_ndf$to == "5" & 1 + cmp_sim_ndf$completed_mv_duration >= cmp_sim_ndf$tend + (obs_start_day - 1), "hazard"] = 0

#         nlst[[i]] <- cbind(nlst[[i]], hazard_ratio = cmp_sim_ndf$hazard / ref_sim_ndf$hazard) # add hazard
#         nlst[[i]]
#       })

#       sim_haz <- do.call(cbind, lapply(nlst, function(df) df$hazard_ratio))
#       # NAs may occur in ratio if both hazards were 0
#       newdata$hr_lower <- apply(sim_haz, 1, quantile, probs = alpha / 2, na.rm = T)
#       newdata$hr_upper <- apply(sim_haz, 1, quantile, probs = 1 - alpha / 2, na.rm = T)

#       newdata
#       }
#   ndf = ndf |> dplyr::mutate(from = sub("->(.*)", "", transition), to = sub("(.*)->", "", transition)) 
#   ref_ndf = cmp_ndf = ndf
#   ref_ndf$hazard = predict(ref_model, newdata = ndf, type = "response") 
#   cmp_ndf$hazard = predict(cmp_model, newdata = ndf, type = "response") 

#   ref_ndf = enforce_deterministic_mvyes_hazards(ref_ndf, obs_start_day = obs_start_day)
#   cmp_ndf = enforce_deterministic_mvyes_hazards(cmp_ndf, obs_start_day = obs_start_day)
#   # For the palliative bias model, additionally enforce hazard of 0 the first two days someone is discharged from mv:yes
#   cmp_ndf[cmp_ndf$from == "2" & cmp_ndf$to == "5" & 2 + cmp_ndf$completed_mv_duration >= cmp_ndf$tend + (obs_start_day - 1), "hazard"] = 0
#   cmp_ndf[cmp_ndf$from == "1" & cmp_ndf$to == "5" & 1 + cmp_ndf$completed_mv_duration >= cmp_ndf$tend + (obs_start_day - 1), "hazard"] = 0

#   ndf = ndf |> 
#     dplyr::left_join(ref_ndf |> dplyr::select(transition, tend, hazard), by = c("transition", "tend")) |> 
#     dplyr::left_join(cmp_ndf |> dplyr::select(transition, tend, hazard), by = c("transition", "tend"), suffix = c("", "_cmp")) |> 
#     dplyr::mutate(hazard_ratio = hazard_cmp / hazard)

#   if (ci){
#     ndf = get_hazard_ci(ndf, ref_model, cmp_model, nsim, alpha)
#   }
#   return(ndf)
# }