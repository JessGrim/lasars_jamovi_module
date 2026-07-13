lasarsmodelClass <- if (requireNamespace('jmvcore', quietly=TRUE)) R6::R6Class(
  "lasarsmodelClass",
  inherit = lasarsmodelBase,
  private = list(
    
    .run = function() {
      
      #LOADING CHECKS/ERRORS####
      if (check_required_inputs(self$options))
        return()
      if (!check_response_options(self))
        return()
      if (!check_direction_requirements(self,
                                       direction_required = self$options$estDirectPref == TRUE,
                                       dir_var_name = 'direction',
                                       rev_label_name = 'rev_label'))
        return()
      if (!check_latent_states(self))
        return()
      
      #THE LL FUNCTION####
      ll_func <- function(x, data, subscale, response, direction = NULL, rev_label, resp_opts, sample = FALSE) {
        
        #Fix sd
        sd = 1
        
        ###DEFINE RESP PARS
        centrePref <- x["centrePref"]
        
        #get the oddsPref estimate (if it's being estimated)
        if("oddsPref" %in% names(x)){
          oddsPref <- x["oddsPref"]
        } else{
          oddsPref <- 0
        }
        
        #Define directionPref based on whether directionPref is in the pars vector
        if("directionPref" %in% names(x)){
          directionPref <- x["directionPref"]
        } else{
          directionPref <- 0
        }
        
        
        cutpoints <- make_thresholds(resp_opts, centrePref, oddsPref, directionPref)
        
        
        #Mus
        scales <- unique(data[[subscale]])
        
        mu <- c()
        for (i in 1:length(scales)){
          mu[i] <- x[paste0("mu.", scales[i])]
        }
        
        names(mu)<- scales
        
        
        if(sample){
          #replace response column with NAs
          data[[response]] <- rep(NA, nrow(data))
          samples <- rep(NA, nrow(data))
          
          for(i in scales) {
            for (j in unique(data[[direction]])) {
              if (j != rev_label){ 
                #use_trials is a vector which indicates which trials to use this time around
                use_trials <- data[[subscale]] == i & data[[direction]] == j
                #samples gets a random sample from a norm dist for all trials which satisfy this loop's conditions
                samples[use_trials] <- rnorm(sum(use_trials), mu[i], sd)
              } else {
                use_trials <- data[[subscale]] == i & data[[direction]] == j
                samples[use_trials] <- rnorm(sum(use_trials), -mu[i], sd)
              }
            }
          }
          
          #the cut function identifies which bin the sample falls into (eg. between cut.1 and cut.2)
          #as.numeric turns these bins into numbers (eg. 2) which correspond with response options
          samples <- as.numeric(cut(samples, c(-Inf, cutpoints, Inf)))
          data$response <- samples
          return(data)
          
        } else{
          #GET LL
          like <- rep(NA, length(data[[response]]))
          
          if(self$options$estDirectPref == TRUE){
            for(i in scales) {
              for (j in unique(data[[direction]])){
                if(j != rev_label){
                  use_trials <- data[[subscale]] == i & data[[direction]] == j
                  p <- diff(c(0, pnorm(c(cutpoints, Inf), mu[i], sd)))
                  like[use_trials] <- p[data[[response]][use_trials]]
                } else if (j == rev_label) {
                  use_trials <- data[[subscale]] == i & data[[direction]] == j
                  p <- diff(c(0, pnorm(c(cutpoints, Inf), -mu[i], sd)))
                  like[use_trials] <- p[data[[response]][use_trials]]
                }
              }
            }
          } else{
            for(i in scales){
              use_trials <- data[[subscale]] == i
              p <- diff(c(0, pnorm(c(cutpoints, Inf), mu[i], sd)))
              like[use_trials] <- p[data[[response]][use_trials]]
            }
          }
          
          out <- sum(log(pmax(like, 1e-10)))
          return(-out)
        }
      }
      
      #CAPTURE INPUT VARIABLES####
      data <- self$data
      response <- self$options$response
      subject <- self$options$subject
      subscale <- self$options$subscale
      direction <- self$options$direction
      estDirectPref <- self$options$estDirectPref
      rev_label <- self$options$rev_label
      resp_opts <- self$options$resp_opts
      subj_est_table <- self$options$subj_est_table
      
      #DEFINE STARTPOINTS (dynamically)####
      #mu pars
      mu_starts <- c(rep(0, length(unique(data[[subscale]]))))
      names(mu_starts) <- c(paste0('mu.', unique(data[[subscale]])))
      
      #Response pars
      if(estDirectPref == TRUE){
        resp_starts <- c("centrePref" = 0,
                         "oddsPref" = 0,
                         "directionPref" = 0)
      } else{
        resp_starts <- c("centrePref" = 0,
                         "oddsPref" = 0)
      }
      
      start_vals <- c(mu_starts, resp_starts)
      
      #RUN OPTIM####
      #Check if there is anything in the results table yet
      if(is.null(self$results$resultsTable$state)){
        
        #Set upper and lower bounds
        upper_bound <- 20
        lower_bound <- -20
        
        #remove NA values before processing
        data <- data[!is.na(data[[response]]), ]
        
        #Run optim over each subject's data and add to a list
        results_list <- lapply(unique(data[[subject]]), function(subj) {
          
          this_subj_data <- data[data[[subject]] == subj, ]
          
          # Early skip for empty / bad data
          if (nrow(this_subj_data) == 0) {
            return(list(
              par = rep(NA_real_, length(start_vals)),
              convergence = 999,
              error = "Empty subject data"
            ))
          }
          
          tryCatch({
            fit <- optim(
              par = start_vals,
              fn = ll_func,
              data = this_subj_data,
              subscale = subscale,
              response = response,
              direction = direction,
              rev_label = rev_label,
              resp_opts = resp_opts,
              method = "L-BFGS-B",
              lower = rep(lower_bound, length(start_vals)),
              upper = rep(upper_bound, length(start_vals))
            )
            
            # normal return
            fit
            
          }, error = function(e) {
            
            # return a fake "failed" result
            list(
              par = setNames(rep(NA_real_, length(start_vals)), names(start_vals)),
              convergence = 999,
              error = e$message
            )
          })
        })
        
        results_df <- process_results(results_list, data, subject, mu_indices = NULL,
                                      upper_bound, lower_bound)
        self$results$resultsTable$setState(results_df)
        
      } else{
        results_df <- self$results$resultsTable$state
      }
      
      #SUBJ-LEVEL ESTS####
      if(subj_est_table){
        fill_subject_results_table(self$results$resultsTable, results_df)
      }
      
      
      ##DESCRIPTIVES####
      #Remove non-converging subjects
      nc_ps <- results_df$subject[results_df$convergence != 'Success']
      results_df <- results_df[!results_df$subject %in% nc_ps, ]
      
      if(length(nc_ps) != 0){
        nonconverg_notice <- jmvcore::Notice$new(options=self$options, name='dynamicNotice')
        nonconverg_notice$setContent(paste0('Estimates from ', length(nc_ps), ' subject(s) removed from Descriptive and Correlation analyses due to convergence concerns.'))
        self$results$insert(1, nonconverg_notice)
      }
      
      #Latent State Descripts
      mu_rows <- colnames(results_df)[grepl("mu", colnames(results_df))]
      fill_descriptive_table(
        self$results$mu_desc,
        results_df,
        mu_rows,
        note_key = "mu_note",
        note_text = "Group-level estimates calculated from subject-level fits"
      )
       
      #Response Style Descripts
      rs_rows <- setdiff(
        colnames(results_df),
        c("subject", "convergence", grep("mu", colnames(results_df), value = TRUE))
      )
      fill_descriptive_table(
        self$results$rs_desc,
        results_df,
        rs_rows,
        note_key = "rs_note",
        note_text = "Group-level estimates calculated from subject-level fits"
      )
      
      #GG_PAIRS PLOT####
      #Make plot data
      plotData <- results_df[, c("subject", grep("mu.", colnames(results_df), value = TRUE))]
      
      #Do some jamovi stuff
      image <- self$results$pairs_plot
      image$setState(plotData)
      
      if(self$options$ls_plot){
        self$results$pairs_plot$setVisible(visible=TRUE)
      }
      
      #SUBJ-LEVEL JAMOVI FILE####
      if (self$options$get_subj_ests) {
        #Retrieve the option object
        option <- self$options$option('get_subj_ests')
        
        #Check for compatibility
        if (is.null(option$perform))
          return()
        
        #Call $perform() with a callback
        option$perform(function(action) {
          list(
            data = results_df,
            title = 'Subject-level Results'
          )
        })
      }
      
    }, #<- end of run brackets!
    
    #PLOTTING FUNCTION
    .pairsplot=function(image, ...) {
      plotData <- image$state
      if (is.null(plotData))
        return()
      build_pairs_plot(plotData)
    })
)
