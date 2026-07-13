origmodelClass <- if (requireNamespace('jmvcore', quietly=TRUE)) R6::R6Class(
  "origmodelClass",
  inherit = origmodelBase,
  private = list(
    .run = function() {
      
      #LOADING CHECKS/ERRORS####
      if (check_required_inputs(self$options))
        return()
      if (!check_response_options(self))
        return()
      if (!check_direction_requirements(self, direction_required = TRUE,
                                       dir_var_name = 'direction',
                                       rev_label_name = 'dirRev'))
        return()
      if (!check_latent_states(self))
        return()
      
      #THE LL FUNCTION####
      ll_func <- function(x, data, subscale, response, direction = NULL, rev_label = NULL, resp_opts) {
        
        no_of_cuts <- resp_opts - 1
        fixed_cut <- ifelse(no_of_cuts%% 2 == 0, (ceiling(no_of_cuts/2) + 1), ceiling(no_of_cuts/2))
        
        
        #Define sd
        sd = exp(x["sd"])
        
        
        ###DEFINE CUTPOINTS
        #assign cut.1
        cut.1 <- x["cut.1"]
        cutpoints <- cut.1
        
        #define cuts up to fixed cut
        if(fixed_cut != 2){
          #loop through until you reach the fixed cut
          for (point in 2:(fixed_cut-1)) {
            #and define cut points
            assign(paste0("cut.", point), get(paste0("cut.", (point-1))) + exp(x[paste0("cut.", point)]))
            cutpoints <- c(cutpoints, get(paste0("cut.", point)))
            names(cutpoints)[point] <- paste0("cut.", point)
          }
        }
        
        #define fixed cut
        assign(paste0("cut.", fixed_cut), get(paste0("cut.", (fixed_cut -1))) + 1)
        cutpoints <- c(cutpoints, get(paste0("cut.", fixed_cut)))
        names(cutpoints)[fixed_cut] <- paste0("cut.", fixed_cut)
        
        
        #loop through remaining cuts and define them
        for (point in (fixed_cut+1):no_of_cuts) {
          assign(paste0("cut.", point), get(paste0("cut.", (point-1))) + exp(x[paste0("cut.", point)]))
          cutpoints <- c(cutpoints, get(paste0("cut.", point)))
          names(cutpoints)[point] <- paste0("cut.", point)
        }
        
        #Mus
        scales <- unique(data[[subscale]])
        
        mu <- c()
        for (i in 1:length(scales)){
          mu[i] <- x[paste0("mu.", scales[i])]
        }
        
        names(mu)<- scales
        
        
        #GET LL
        like <- rep(NA, length(data[[response]]))
        
        if(!is.null(direction)){
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
      
      
      #CAPTURE INPUT VARIABLES####
      data <- self$data
      response <- self$options$response
      subject <- self$options$subject
      subscale <- self$options$subscale
      direction <- self$options$direction
      rev_label <- self$options$dirRev
      resp_opts <- self$options$resp_opts
      subj_est_table <- self$options$subj_est_table
      
      #DEFINE STARTPOINTS (dynamically)####
      #mu pars
      mu_starts <- c(rep(0, length(unique(data[[subscale]]))))
      names(mu_starts) <- c(paste0('mu.', unique(data[[subscale]])))
      
      #Response pars
      no_of_cuts <- resp_opts - 1
      cut_pars <- paste0("cut.", 1:no_of_cuts)
      
      #remove the middle one
      fixed_cut <- ifelse(no_of_cuts%% 2 == 0, (ceiling(no_of_cuts/2) + 1), ceiling(no_of_cuts/2))
      cut_pars <- cut_pars[-fixed_cut]
      
      resp_starts <- rep(0, length(cut_pars))
      names(resp_starts) <- cut_pars
      resp_starts <- c(resp_starts, 'sd' = 1)
      
      start_vals <- c(mu_starts, resp_starts)
      
      #RUN OPTIM####
      #Check if anything is already in the state
      if(is.null(self$results$resultsTable$state)){
        
        #Set upper and lower bounds
        upper_bound <- 20
        lower_bound <- -20
         
        #remove NA values before processing
        data <- data[!is.na(data[[response]]), ]
         
        #Run optim over each subject's data and add to a list
        results_list <- lapply(unique(data[[subject]]), function(subj) {
           
          this_subj_data <- data[data[[subject]] == subj, ]
           
          # Optional: early skip for empty / bad data
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
         
        results_df <- process_results(results_list, data, subject,
                                      mu_indices = grep("mu.", names(results_list[[1]]$par)),
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
    }
  )
)

