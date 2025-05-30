
#' @title cheesepop: Disaggregating population counts for a single level of demographics
#' (e.g., age groups only or sex group only) - no geospatial covariates.
#' Please use 'slices' if you want covariates used.
#'
#' @description This function disaggregates population estimates by a single demographic (age or sex or religion, etc)
#'
#' @param df A data frame object containing sample data (often partially observed) on age or sex groups population data
#' as well as the estimated overall total counts per administrative unit.
#'
#' @param output_dir This is the directory with the name of the output folder where the
#' disaggregated population proportions and population totals are
#' automatically saved.
#'
#' @return Data frame objects of the output files including the disaggregated population proportions and population totals
#' along with the corresponding measures of uncertainties (lower and upper bounds of 95-percent credible intervals) for each demographic characteristic.
#' In addition, a file containing the model performance/model fit evaluation metrics is also produced.
#'
#'@examples
#'data(toydata)
#' result <- scissors(df = toydata, output_dir = tempdir())
#' @export
#' @importFrom dplyr "%>%"
#'
#'
scissors <-function(df,output_dir)# disaggregates by age only - no covariates
{

  # Check if the output directory exists, if not, create it
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message(paste("Directory", output_dir, "created successfully."))
  } else {
    message(paste("Directory", output_dir, "already exists."))
  }

  # df = dat
  # Partially observed age data
  age_df <- df %>% dplyr::select(starts_with("age_"))
  pred_dt <- prop_dt <- age_df
  pred_dtL <- prop_dtL <- age_df
  pred_dtU <- prop_dtU <- age_df


  age_classes <- names(age_df)
  # Add the total population to the age data
  age_df$total <- df$total
  age_df$total[age_df$total==0] = NA

  age_df$ID <- 1:nrow(age_df) # add IID

  for(i in 1:length(age_classes))
  {

    #  i=2
    # 1) Disaggregate by age - estimate missing age group proportion for each admin unit
    prior.prec <- list(prec = list(prior = "pc.prec",
                                   param = c(1, 0.01))) # using PC prior
    print(paste(paste0("(",i,")"),paste0(age_classes[i], " model is running")))

    age_df[,colnames(age_df)[i]] <- round(age_df[,i])

    form_age <- as.formula(paste0(colnames(age_df)[i], " ~ 1"))
    mod_age  <- inla(form_age,
                     data = age_df,
                     family = "binomial", Ntrials = total,
                     control.predictor = list(compute = TRUE),
                     control.compute = list(dic = TRUE, cpo = TRUE)
    )
    prop_dt[,i] = round(plogis(mod_age$summary.linear.predictor$mean),4)
    prop_dtL[,i] = round(plogis(mod_age$summary.linear.predictor$'0.025quant'),4)
    prop_dtU[,i] = round(plogis(mod_age$summary.linear.predictor$'0.975quant'),4)

    pred_dt[,i] = round(prop_dt[,i]*age_df$total)
    pred_dtL[,i] = round(prop_dtL[,i]*age_df$total)
    pred_dtU[,i] = round(prop_dtU[,i]*age_df$total)
  }

  ## Rename Columns
  # age combined
  # count
  age_classes_pop = paste0("pp_",age_classes)
  colnames(pred_dt) <- age_classes_pop # predicted population counts
  age_classes_popL = paste0(age_classes_pop,"L")
  colnames(pred_dtL) <- age_classes_popL # lower bound of predicted population counts
  age_classes_popU = paste0(age_classes_pop,"U")
  colnames(pred_dtU) <- age_classes_popU # upper bound of predicted population counts

  # proportion
  age_classes_prop = paste0("prp_",age_classes)
  colnames(prop_dt) <- age_classes_prop # predicted population counts
  age_classes_propL = paste0(age_classes_prop,"L")
  colnames(prop_dtL) <- age_classes_propL # lower bound of predicted population counts
  age_classes_propU = paste0(age_classes_prop,"U")
  colnames(prop_dtU) <- age_classes_propU # upper bound of predicted population counts


  ## predictive ability checksall_prop <- f.prop_dt + m.prop_dt # should be a matrix of 1s only!
  all_pop <- apply(pred_dt, 1, sum)# NO NA

  png(paste0(output_dir,"/model_validation_scatter_plot.png"))
  plot(all_pop, df$total,
       xlab="Observed population", ylab = "Predicted population",
       main= "Scatter plot of \n observed versus predicted")
  abline(0,1, col=2, lwd=2) # should be a straight perfect fit
  dev.off()


  residual = all_pop-df$total
  print(mets <- t(c(MAE = mean(abs(residual), na.rm=T),#MAE
                    MAPE = (1/length(df$total))*sum(abs((df$total-all_pop)/df$total))*100,#MAPE
                    RMSE = sqrt(mean(residual^2, na.rm=T)),
                    corr = cor(df$total[!is.na(df$total)],all_pop[!is.na(df$total)]))))# should be with at least 95% correlation

  write.csv(t(mets), paste0(output_dir,"/fit_metrics.csv"),row.names = F)


  # join all data
  full_dat <- cbind(df,
                    pred_dt, pred_dtL,pred_dtU) # everything

  # save the datasets
  write.csv(full_dat, paste0(output_dir,"/full_disaggregated_data.csv"),row.names = F)
  write.csv(pred_dt, paste0(output_dir,"/age_disaggregated_data.csv"),row.names = F)
  write.csv(prop_dt, paste0(output_dir,"/age_proportions.csv"),row.names = F)

  # return output as a list
  return(out <- list(age_pop = pred_dt,
                     age_popL = pred_dtL,
                     age_popU = pred_dtU,
                     age_prop = prop_dt,
                     full_data = full_dat))

}

#githublink <- "https://raw.github.com/wpgp/CheeseCake/main/example_data.Rdata"
#load(url(githublink))
#dat <- read.csv("C:/Users/ccn1r22/OneDrive - University of Southampton/Documents/packages/example_data.csv")
#output_dir <- "C:/Users/ccn1r22/OneDrive - University of Southampton/Documents/packages/output5"
#output_dir <- tempdir()# please specify your output folder here
#system.time(age_dis <- scissors(example_data,output_dir))
