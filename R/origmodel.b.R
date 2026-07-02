OrigmodelClass <- if (requireNamespace('jmvcore', quietly=TRUE)) R6::R6Class(
  "OrigmodelClass",
  inherit = OrigmodelBase,
  private = list(
    .run = function() {
      
      #LOADING CHECKS/ERRORS####
      #Don't run the analysis until something is in the boxes
      if (is.null(self$options$subject))
        return()
      if (is.null(self$options$response))
        return()
      if (is.null(self$options$subscale))
        return()
      
      #wait until there is something in the direction box
      if (is.null(self$options$direction)){
        
        dirMissing_notice <- jmvcore::Notice$new(options=self$options, name='dirMissing')
        dirMissing_notice$setContent("This analysis requires a directional variable.
                                       \nIf none exists, use the Latent States and Response Style Analysis")
        self$results$insert(1, dirMissing_notice)
        
        return()
      }
      
      #Ask for a way to identify reverse coded items
      if(!is.null(self$options$direction) &&
         (is.null(self$options$dirRev) ||
          is.na(self$options$dirRev) ||
          self$options$dirRev == "")
      )
      {
        negNeeded_notice <- jmvcore::Notice$new(options=self$options, name='negNeeded')
        negNeeded_notice$setContent(paste0("Please enter the reverse-coding identifying label.
                                            \nThis is the string or character from the '", self$options$direction, "' column 
                                            \nused to indicate that an item is reverse-coded."))
        self$results$insert(1, negNeeded_notice)
        
        return()
      }
      
      #Check for more than 1 latent state
      #wait until there is something in the direction box
      scale_levels <- unique(self$data[[self$options$subscale]])
      if (!is.null(self$options$subscale) && length(scale_levels) <= 1){
        
        tooFewLS_notice <- jmvcore::Notice$new(options=self$options, name='tooFewLS')
        tooFewLS_notice$setContent(
          paste0(
            "This analysis requires more than 1 latent state.\n",
            "The '", self$options$subscale, "' column only contains levels: ",
            paste(scale_levels, collapse = ", ")
          )
        )
        
        self$results$insert(1, tooFewLS_notice)
        
        return()
      }
      
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
        
        #Get subj IDs
        subjects <- unique(data[[subject]])
        
        #Bind all the results together
        results_df <- do.call(rbind, lapply(seq_along(results_list), function(i) {
          c(
            round(results_list[[i]]$par[grep("mu.", names(results_list[[i]]$par))], 3),
            convergence = results_list[[i]]$convergence
          )
        }))
        
        #MANAGE RESULTS####
        results_df <- as.data.frame(results_df)
        
        #Add subject IDs
        rownames(results_df) <- unique(data[[subject]])
        
        #Record convergence code
        results_df$convergence <- ifelse(results_df$convergence == 0, "Success",
                                         ifelse(results_df$convergence == 1, "Error - Limit Reached", paste0("Error - optim code: ", results_df$convergence)))
        
        #Override convergence code if bound reached
        check_cols <- setdiff(names(results_df), "convergence")
        
        has_upper <- apply(results_df[check_cols], 1, function(x) any(x == upper_bound))
        has_lower <- apply(results_df[check_cols], 1, function(x) any(x == lower_bound))
        
        results_df$convergence <- ifelse(
          has_upper & has_lower, "Both bounds reached",
          ifelse(
            has_upper, "Upper Bound Reached",
            ifelse(
              has_lower, "Lower Bound Reached",
              results_df$convergence
            )
          )
        )
        
        results_df$subject <- rownames(results_df)
        results_df <- results_df[, c("subject", setdiff(names(results_df), "subject"))]
        
        #Set the state to save object
        self$results$resultsTable$setState(results_df)
        
      } else{
        results_df <- self$results$resultsTable$state
      }
      
      #SUBJ-LEVEL ESTS####
      if(subj_est_table){
        #FILL SUBJ-RESULTS TABLE####
        for(add_col in colnames(results_df)){
          self$results$resultsTable$addColumn(name = add_col, format = "zto,3")
        }
        for(add_row in 1:(nrow(results_df)-1)){
          self$results$resultsTable$addRow(add_row)
        }
        
        #Fill the table
        for (i in seq_len(nrow(results_df))) {
          self$results$resultsTable$setRow(
            rowNo = i,
            values = as.list(results_df[i, , drop = FALSE])
          )
        }
        
        self$results$resultsTable$setVisible(visible=TRUE)
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
      mu_desc_table <- self$results$mu_desc
      
      for(i in 1:(length(mu_rows) - 1)){
        mu_desc_table$setRow(
          rowNo = i, values = list(
            par = mu_rows[i],
            num = length(results_df[[mu_rows[i]]]),
            "mean" = mean(results_df[[mu_rows[i]]]),
            med = median(results_df[[mu_rows[i]]]),
            "sd" = sd(results_df[[mu_rows[i]]]),
            se = sd(results_df[[mu_rows[i]]], na.rm = TRUE) / sqrt(sum(!is.na(results_df[[mu_rows[i]]])))
          )
        )
        mu_desc_table$addRow(i)
      }
      final_mu_row <- length(mu_rows)
      mu_desc_table$setRow(
        rowNo = final_mu_row, values = list(
          par = mu_rows[final_mu_row],
          num = length(results_df[[mu_rows[final_mu_row]]]),
          mean = mean(results_df[[mu_rows[final_mu_row]]]),
          med = median(results_df[[mu_rows[final_mu_row]]]),
          sd = sd(results_df[[mu_rows[final_mu_row]]]),
          se = sd(results_df[[mu_rows[final_mu_row]]], na.rm = TRUE) / sqrt(sum(!is.na(results_df[[mu_rows[final_mu_row]]])))
        )
      )
      mu_desc_table$setNote(key = "mu_note", note = "Group-level estimates calculated from subject-level fits")
      
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
      
      if(is.null(plotData)){
        return()
      }
      
      calculate_font_size <- function(x) {
        if (x == 0) {
          return(3)  # Font size for 0
        } else {
          return(3 + 17 * (abs(x) - 0) / (1 - 0))  # Font size scaled from 3 to 20
        }
      }
      
      just_scale_mus <- plotData[, 2:ncol(plotData)]
      
      #Call to ggplot etc.
      p_model <- GGally::ggpairs(plotData, 
                                 columns = 2:ncol(plotData),
                                 upper = 'blank', #set the upper quadrant to be blank
                                 #diag = 'blank', #set diagonal to be blank
                                 lower = list(continuous = GGally::wrap("points", size = 0.6, alpha = 1)), #make scatterplots in lower quadrant
                                 axisLabels = "none", #remove axis labels/ticks
                                 columnLabels = NULL #remove the column and row (subscale) labels
      )
      
      #Print scale names in diagonal cells
      for(i in 1:ncol(just_scale_mus)){
        this_scale <- colnames(just_scale_mus)[i]
        p_model[i,i] <- GGally::ggally_text(this_scale,
                                            color = "black",
                                            size = 3.5)
      }
      
      #print cor coeffs
      for (i in 1:(ncol(just_scale_mus) - 1)){
        for(j in (i + 1):length(just_scale_mus)){
          this_cor <- round(cor(just_scale_mus[i], just_scale_mus[j], method = 'spearman'), 2)
          font_size <- as.numeric(calculate_font_size(this_cor))
          this_cor_char <- as.character(this_cor)
          #if the first character in the string is '-'
          if(substr(this_cor_char, 1, 1) == '-'){
            #start from the 3rd char in the string (the full-stop)
            #and paste '-' in front of it 
            if(this_cor_char != "-1" ){
              this_cor_char <- paste0("-", substring(this_cor_char, first = 3))
            }
          } else{
            #start from the 2nd char in the string (essentially, remove the first)
            if(this_cor_char != "1" && this_cor_char != "0"){
              this_cor_char <- substring(this_cor_char, first = 2)
            }
          }
          p_model[i,j] <- GGally::ggally_text(this_cor_char,
                                              color = "black",
                                              size = font_size)
        }
      }
      
      #Remove the grid lines and background colours
      p_model <- p_model + ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                                          panel.grid.major = ggplot2::element_blank(),
                                          panel.background = ggplot2::element_rect(fill = "white"),
                                          panel.border = ggplot2::element_rect(colour = "black", fill = NA))
      
      return(p_model)
    }
  )
)

