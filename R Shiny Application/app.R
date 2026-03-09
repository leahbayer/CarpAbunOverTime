# Shiny Web Application for SAW/DAW 
# Lets user choose gear type, abundance level, effect size, and time span
# Fits NB intercept + offset model to data
# Runs simulation across sample sizes given a specified effect size
# Runs repeated fits to test if trend (Year) is significant
# Computes power as proportion of significant results
# Returns an interactive plot of power vs sample size for chosen scenario
# Returns table with minimum sample size for each target power
# Returns table with summary of key results
# Returns raw power-by-sample size results downloadable as csv, excel, or pdf
# Includes help tab with detailed information on how to use and interpret app

library(shiny) # web app framework
library(tidyverse) # data manipulation
library(MASS) # nb modeling
library(tools) 
library(DT) # interactive data tables
library(shinyWidgets) # widgets
library(plotly) # interactive ggplot
library(shinyBS) # help buttons


# ------------------------ Helper functions  -------------------------------

# Simulate negative binomial counts over multiple years
simulate_counts <- function(n, params, eff, theta, offset, time = 1:3) {
  # n - sample size per year
  # params - log-scale intercepts
  # eff - effect size => slope coefficient on log scale divided by time
  # theta - dispersion
  # offset - log(sampling time) (5 min for DT; 15 min for EF)
  # time - years spent sampling (1:3, 1:5, c(1, 3, 5))
  
  counts <- matrix(NA, nrow = n, ncol = length(time))
  
  for (yr in seq_along(time)) {
    exp_count <- exp(params + offset + eff * time[yr])
    
    counts[, yr] <- rnbinom(n, mu = exp_count, size = theta)
  }
  
  sim_data <- data.frame(y = as.vector(counts),
                         Year = rep(time, each = n),
                         offset = offset)
  return(sim_data)
}

# Run repeated simulations for a given sample size & estimate power
simulated_power <- function(sample_size, n_sim = 100, params, eff, theta,
                            alpha = 0.05, offset, time){
  p_val <- numeric(n_sim)
  
  for(i in seq_len(n_sim)){
    count <- simulate_counts(n = sample_size,
                             params = params,
                             eff = eff,
                             theta = theta,
                             offset = offset,
                             time = time)
    
    # testing for warnings or errors from the fitted model
    test <- try(assertCondition(glm.nb(y ~ Year + offset(offset),
                                       data = count)), silent = T)
    
    # simulating data until no warnings
    while(any(class(test[[1]]) %in%
              c('simpleWarning', 'simpleError'))){
      count <- simulate_counts(n = sample_size,
                               params = params,
                               eff = eff,
                               theta = theta,
                               offset = offset,
                               time = time)
      
      # testing for warnings or errors from the fitted model
      test <- try(assertCondition(glm.nb(y ~ Year + offset(offset),
                                         data = count)), silent = T)
    }
    
    fit <- glm.nb(y ~ Year + offset(offset), data = count)
    
    p_val[i] <- coef(summary(fit))["Year", "Pr(>|z|)"]
  }
  
  return(mean(p_val < alpha, na.rm = TRUE))
}


# ------------------------ User interface  -------------------------------

ui <- fluidPage(
  titlePanel("Carp Power Analysis Tool (Sample Size Estimation)"),
  sidebarLayout(
    sidebarPanel(
      selectInput("gear", "Gear type:",
                  choices = c("Electrified Dozer Trawl (DT)" = "DT",
                              "Boat Electrofishing (EF)" = "EF")),
      
      selectInput("abundance", "Abundance level:",
                  choices = c("High" = "High", 
                              "Moderate" = "Moderate", 
                              "Rare" = "Rare", 
                              "Custom (enter below)" = "Custom"), 
                  selected = "High"),
      # had to add diff way to show help buttons since Shiny no longer supports
      conditionalPanel(
        condition = "input.abundance == 'Custom'",
        
        numericInput(
          "custom_abundance",
          "Custom abundance (carp per minute):",
          value = 1,
          min = 0,
          step = 1
        ),
        
        numericInput(
          inputId = "custom_time_unit",
          label = tags$span(
            "Sampling time per sample (minutes): ",
            tags$button(
              "?",
              id = "custom_time_help",
              type = "button",
              class = "btn btn-info btn-xs",
              `data-bs-toggle` = "tooltip",
              title = paste(
                "This is the amount of time (in minutes) spent on a single sampling event.",
                "It is used to convert the custom abundance rate (carp per minute)",
                "into an expected mean count per sample.",
                "For example, 1 carp/min over 5 minutes implies an expected mean of 5 carp per sample."
              )
            )
          ),
          value = 1,
          min = 1,
          step = 1
        ),
        
        numericInput(
          inputId = "custom_sd",
          label = tags$span(
            "Custom standard deviation (count scale): ",
            tags$button(
              "?",
              id = "custom_sd_help",
              type = "button",
              class = "btn btn-info btn-xs",
              `data-bs-toggle` = "tooltip",
              title = paste(
                "This is the standard deviation of carp counts per sample (count scale).",
                "It controls how variable catches are around the mean.",
                "Variance must exceed the mean for the negative binomial model.",
                "Larger values represent more heterogeneous or patchy catches.",
                "If unsure, a standard deviation approximately equal to the expected mean is a reasonable starting assumption for field-based fish sampling data."
              )
            )
          ),
          value = 1,
          min = 1,
          step = 1
        ),
        uiOutput("custom_diagnostic")
      ),
      
      selectInput("trend_choice", "Choose trend (reduction in abundance):",
                  choices = c("15% Reduction" = "0.15",
                              "25% Reduction" = "0.25",
                              "50% Reduction" = "0.50",
                              "Custom (enter percent below)" = "custom")),
      # had to add diff way to show help buttons since Shiny no longer supports
      conditionalPanel(
        condition = "input.trend_choice == 'custom'",
        numericInput(
          inputId = "custom_red",
          label = tags$span(
            "Custom percent reduction (0-0.99): ",
            tags$button(
              "?", 
              id = "custom_red_help",
              type = "button",
              class = "btn btn-info btn-xs",
              `data-bs-toggle` = "tooltip",
              title = "Enter the percent reduction per year as a decimal (e.g., 0.25 = 25%)."
            )
          ),
          value = 0.25, min = 0, max = 0.99, step = 0.01
        )
      ),
      
      radioButtons("years_sampled", "Sampling schedule:",
                   choices = c("Annually for 3 years" = "1:3",
                               "Annually for 5 years" = "1:5",
                               "Year 1, 3, and 5" = "c(1,3,5)")),
      
      numericInput(inputId = "alpha",
                   label = tags$span(
                     "Type I error (alpha): ",
                     tags$button(
                       "?",
                       id = "alpha_help",
                       type = "button",
                       class = "btn btn-info btn-xs",
                       `data-bs-toggle` = "tooltip",
                       title = paste(
                         "Alpha is the probability of falsely detecting a trend when none exists.",
                         "A common choice is 0.05 (5%).",
                         "Lower alpha values make the test more conservative and reduce power;",
                         "higher values increase power but raise false-positive risk."
                       )
                     )
                   ),
                   value = 0.05,
                   min = 0.001,
                   max = 0.2,
                   step = 0.005
      ),
      
      sliderInput("sample_range", "Sample size range per year:",
                  min = 10, max = 500, value = c(25, 200), step = 5),
      
      sliderInput(
        inputId = "sample_by",
        label = tags$span("Run simulations for every ____ number of samples: ",
          bsButton("sample_help", label = "?", style = "info", 
                   size = "extra-small")),
        min = 5, max = 50, value = 25, step = 5
      ),
      bsPopover(
        id = "sample_help",
        title = "Note:",
        content = "Step size controls which sample sizes to simulate. Smaller step sizes improve precision but increase runtime.",
        placement = "right",
        trigger = "hover"
      ),
      
      sliderInput(
        inputId = "n_sim",
        label = tags$span("Number of simulations per sample size: ",
          bsButton("sample_help2", label = "?", style = "info", 
                   size = "extra-small")),
        min = 200, max = 1000, value = 500,step = 100
      ),
      bsPopover(
        id = "sample_help2",
        title = "Note:",
        content = "More simulations improve accuracy but increase runtime.",
        placement = "right",
        trigger = "hover"
      ),
      
      actionButton("run_btn", "Run Simulation", class="btn-primary")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Help / About", br(),
                 h3("How to Use the Carp Power Analysis Tool"),
                 p("This tool estimates how many samples are needed to detect a 
                 decline in Asian carp abundance using power analysis with 
                 simulated data. Initial estimates of relative abundance and 
                 its associated uncertainty are obtained from negative binomial 
                 (NB) models fitted to existing data for each gear type and 
                 abundance level. Using these empirical estimates of mean and 
                 variance, the app simulates repeated sampling across years to 
                 estimate statistical power under different sample size and 
                   relative abundance scenarios."),
                 
                 h4("⚙️ How the Simulation Works"),
                 tags$ul(
                   tags$li("A negative binomial model is fitted to observed 
                           data for the selected gear and abundance level. 
                           Output from this fitted model is used to obtain 
                           estimates of initial relative abundance and its 
                           associated variance that are used in subsequent 
                           simulations."),
                   tags$li("An effect size (annual percent reduction in 
                           abundance) and sampling schedule is specified by 
                           the user."),
                   tags$li("Data are simulated based on the selected sampling 
                           schedule and sample size range."),
                   tags$li("For each simulated dataset, a NB regression tests 
                           whether the null hypothesis of no trend in relative 
                           abundance in rejected. This test assumes a type-I 
                           error probability of 0.05 (i.e., \alpha = 0.05)."),
                   tags$li("Power is estimated as the proportion of simulations 
                           where the null hypothesis is correctly rejected.")
                 ),
                 
                 h4("🗂️️ Tabs Explained"),
                 tags$ul(
                   tags$li(tags$b("Power Curve:"), 
                           "Shows the relationship between sample size and 
                           statistical power. Dashed lines mark 60%, 70%, and 
                           80% power thresholds, and red vertical lines 
                           indicate the minimum sample sizes to reach them."),
                   tags$li(tags$b("Data Table:"), 
                           "Lists the smallest sample size needed to achieve 
                           60%, 70%, and 80% power."),
                   tags$li(tags$b("Data Summary:"), 
                           "Displays key model parameters, simulation inputs, 
                           and a table of minimum sample sizes. You can also 
                           download the full raw power-by-sample-size 
                           results."),
                   tags$li(tags$b("Model Fit:"), 
                           "Shows the fitted negative binomial model summary 
                           used for simulations.")
                 ),
                 
                 h4("📝 Input Options Explained"),
                 tags$ul(
                   tags$li(tags$b("Gear type:"), 
                           "Selects between Dozer Trawl (DT, 5-minute sampling 
                           runs) and Electrofishing (EF, 15-minute sampling 
                           runs)."),
                   tags$li(tags$b("Abundance level:"), 
                           "Determines the initial relative abundance (High, 
                           Moderate, Rare, or Custom). Simulated trends in 
                           relative abundance are calculated relative to this 
                           abundance level."),
                   tags$li(tags$b("Custom abundance level:"), 
                           "Allows the user to define a custom relative 
                           abundance instead of using preset data categories. 
                           This is useful if expected catch rates differ from 
                           the historical data. The value entered represents 
                           the expected number of carp per minute of sampling. 
                           The sampling time (in minutes) is then used to 
                           convert this rate into an expected mean count per 
                           sample. Internally, the app calculates the espected 
                           mean count per sample as expected mean = (carp per 
                           minute) * (minutes per sample). This expected mean 
                           is then converted to a log-scale intercept and used 
                           directly in the simulations. Unrealistic values 
                           (very low or very high expected counts) may result 
                           in unstable power estimates."),
                   tags$li(tags$b("Custom standard deviation:"), 
                           "Specifies the expected variability in carp counts 
                           per sample on the count scale. This value is used to 
                           calculate the dispersion parameter (theta) for the 
                           negative binomial distribution. The variance must 
                           exceed the mean for the model to be valid. Larger 
                           standard deviations correspond to more variable 
                           catches."),
                   tags$li(tags$b("Custom sampling time:"), 
                           "Specifies the duration (in minutes) of a single 
                           sampling event when a custom abundance is used. This 
                           value is combined with the custom abundance rate 
                           (carp per minute) to calculate the expected mean 
                           number of carp per sample used in the simulations. 
                           For example, an abundance of 0.5 carp per minute 
                           with a 10-minute sampling time corresponds to an 
                           expected mean of 5 carp per sample."),
                   tags$li(tags$b("Reduction in abundance:"), 
                           "Specifies the expected yearly decline in mean 
                           catch (e.g., a 25% reduction = 0.75× realtive 
                           abundance the previous year)."),
                   tags$li(tags$b("Sampling schedule:"), 
                           "Defines how many years are sampled (e.g., every 
                           year for 3 consecutive years)."),
                   tags$li(tags$b("Type I error (alpha):"),
                           "Defines the probability of falsely detecting a 
                           trend when no real change exists. A common choice is 
                           alpha = 0.05, meaning a 5% chance of a false 
                           positive. Lower alpha values make the test more 
                           conservative and generally reduce power."),
                   tags$li(tags$b("Sample size range per year:"), 
                           "Sets the range of total samples per year to 
                           evaluate. E.g., by default, the app simulates data 
                           from 25 - 200 sampling events annually."),
                   tags$li(tags$b("Run simulations for every __ number of 
                                  samples:"), 
                           "Determines the increment of sample sizes tested 
                           (smaller steps = more detailed but slower 
                           simulations). E.g., by default, the app will 
                           simulate sample sizes of 25, 50, 75, ..., 200."),
                   tags$li(tags$b("Number of simulations per sample size:"), 
                           "Each sample size is simulated this many times to 
                           estimate power (more simulations = smoother results, 
                           longer runtime).")
                 ),
                 
                 h4("💡 Notes"),
                 tags$ul(
                   tags$li("Power is the probability that the analysis will 
                           detect a real decline in abundance."),
                   tags$li("An 80% power level is typically considered 
                           acceptable for study design."),
                   tags$li("If power never exceeds 0.8 within the sample range, 
                   consider increasing sample size, 
                    number of years, or effect size (larger declines are easier 
                           to detect)."),
                   tags$li("Simulation results may vary slightly between runs 
                           because of random variation in simulated data."),
                   tags$li("When using a custom abundance level, note that 
                   unrealistic values (too high or too low) can lead to 
                           unstable model results or misleading power 
                           estimates.")
                 ),
                 
                 h4("📁 Output Files"),
                 p("Use the download buttons in the Data 
                 Summary tab to export raw results in .csv, .xlsx, or .pdf, 
                 including the power value 
                   for each sample size tested.")
        ),
        
        tabPanel("Power Curve", br(),
                 h3("Power Curve"),
                 helpText("This plot shows the relationship between sample 
                 size and statistical power for the chosen scenario. Each point 
                 represents the estimated probability of detecting a decline 
                 in carp abundance under the selected trend, gear type, and 
                          abundance level."),
                 helpText("Horizontal dashed lines mark 60%, 70%, and 80% power 
                 thresholds. Red vertical lines indicate the smallest sample 
                 size needed to reach each threshold."),
                 plotlyOutput("power_plot", height = "600px", width = "1200px")),
        tabPanel("Data Table", br(),
                 h3("Minimum Sample Sizes"),
                 helpText("This table lists the smallest sample size (per year) 
                 needed to achieve 80% statistical power. If 
                 power never exceeds one of these thresholds within the chosen 
                          sample range, it will display ‘NA’."),
                 helpText("This table can be used to guide planning for the 
                 number of samples required in your monitoring design."),
                 DTOutput("powerTable")),
        tabPanel("Simulation Summary", br(),
                 h3("Simulation Summary"),
                 helpText("This tab summarizes the key parameters and results 
                 from your power analysis simulation. It provides both the 
                 model parameters used and the smallest sample sizes needed to 
                          achieve specified power levels."),
                 
                 helpText("Below is an explanation of each metric shown in the 
                          summary table:"),
                 
                 tags$ul(
                   tags$li(tags$b("Intercept (log):"), 
                           " The log expected carp count under 
                           baseline conditions. This represents the log 
                           abundance before any decline is applied and is 
                           estimated from data."),
                   tags$li(tags$b("Theta (dispersion):"), 
                           " The dispersion parameter from the negative 
                           binomial model, controlling how variable the counts 
                           are around the mean. Smaller theta values indicate 
                           higher dispersion."),
                   tags$li(tags$b("Effect size (log scale):"), 
                           " The log-scale slope corresponding to the 
                           expected annual reduction in abundance (e.g., a 25% 
                           reduction corresponds to log(0.75) = -0.288)."),
                   tags$li(tags$b("Reduction in abundance:"), 
                           " The user-defined or preset annual decline in mean 
                           abundance, expressed as a percentage (e.g., 0.25 = 
                           25% reduction per year)."),
                   tags$li(tags$b("Minimum sample size for 60% power, minimum 
                                  sample size for 70% power, minimum sample 
                                  size for 80% power:"), 
                           " The smallest sample size per year required to 
                           reach 60%, 70%, and 80% statistical power, 
                           respectively.")
                 ),
                 
                 helpText("Use these results to determine the sample size 
                 required to reliably detect changes in carp abundance over 
                          time under the specified conditions."),
                 
                 helpText("The full power-by-sample-size results can be 
                 downloaded below for further exploration or record-keeping."),
                 DTOutput("summary_table")
        ),
        tabPanel("Model Fit", br(),
                 h3("Model Fit"),
                 helpText("This tab displays the fitted negative binomial (NB) 
                 model used to estimate baseline abundance and variance for the 
                 selected gear type and abundance level."),
                 helpText("The ‘Model Summary’ output shows the estimated 
                 intercept (average log-count) and dispersion (theta) 
                 parameter. These values are used in the simulation stage."),
                 verbatimTextOutput("model_summary")
      )
    )
  )
),

# JavaScript for tooltip toggle help buttons to create hovering functionality
tags$script(HTML("
  var tooltipTriggerList = [].slice.call(document.querySelectorAll('[data-bs-toggle=\"tooltip\"]'))
  var tooltipList = tooltipTriggerList.map(function (tooltipTriggerEl) {
    return new bootstrap.Tooltip(tooltipTriggerEl)
  })
")),

# style for regular help icons
tags$style(HTML("
  .help-icon {
    display: inline-block;
    padding: 3px 6px;
    background: #0b3d91;
    color: white;
    border-radius: 50%;
    font-weight: bold;
    cursor: pointer;
  }
"))
)



# --------------------------- Server ----------------------------------

server <- function(input, output, session) {
  
  # Load data fits
  saved_fits <- readRDS("saved_fits.rds")
  
  # Load both datasets
  source('DTdata_Clean_ShinyFriendly.R')
  source('EFdata_Clean_ShinyFriendly.R')
  
  DTdata <- UMRdata_clean
  EFdata <- ILRdataEF_clean
  
  # Separate DT analyses by abundance level
  DTdata_high <- DTdata %>% filter(Abundance == 'High')
  DTdata_mod <- DTdata %>% filter(Abundance == 'Moderate')
  DTdata_rare <- DTdata %>% filter(Abundance == 'Rare')
  
  # Update abundance choices based on gear
  observeEvent(input$gear, {
    if (input$gear == "EF") {
      updateSelectInput(session, "abundance",
                        choices = c("Moderate"),
                        selected = "Moderate")
    } else {
      updateSelectInput(session, "abundance",
                        choices = c("High", "Moderate", "Rare"),
                        selected = "High")
    }
  }, ignoreInit = TRUE)
  
  # Returns correct dataset based on inputs
  data_for_model <- reactive({
    if (input$gear == "DT") {
      if (input$abundance == "Custom") {
        return(DTdata) 
      }
      
      switch(input$abundance,
             "High"   = DTdata_high,
             "Moderate" = DTdata_mod,
             "Rare"   = DTdata_rare)
    } else {
      EFdata 
    }
  })
  
  # Fit NB model
  fitted_model <- reactive({
    
    # load from saved_fits
    if (!is.null(saved_fits)) {
      
      if (input$gear == "DT") {
        
        if (input$abundance == "Custom") {
          
          # change to include custom abundance on count scale
          lambda <- input$custom_abundance * input$custom_time_unit
          var_count <- input$custom_sd ^ 2
          
          validate(
            need(
              var_count > lambda,
              "Variance must be greater than the mean for the negative binomial model. 
     Increase the SD or reduce expected abundance."
            )
          )
          
          theta <- lambda ^ 2 / (var_count - lambda)
          custom_coef <- log(lambda)
          
          return(list(model = NULL, 
                      coef = custom_coef, 
                      theta = theta, 
                      custom = TRUE))
          } 
        
        else {
          fit <- saved_fits$DT[[input$abundance]]
          return(list(model = fit, 
                      coef = as.numeric(coef(fit)[1]), 
                      theta = as.numeric(fit$theta), 
                      custom = FALSE))
          }
        } 
      
      else {
        fit <- saved_fits$EF$Moderate
        return(list(model = fit, 
                    coef = as.numeric(coef(fit)[1]), 
                    theta = as.numeric(fit$theta), 
                    custom = FALSE))
        }
      } 
    
    else {
      
      # fallback to original fitting
      df <- data_for_model()
      offset_vec <- rep(ifelse(input$gear == "DT", log(5), log(15)), nrow(df))
      fit <- try(glm.nb(CountLarge ~ 1 + offset(offset_vec), data = df), 
                 silent = TRUE)
      
      if (input$abundance == "Custom") {
        custom_coef <- log(input$custom_abundance)
        fitted_theta <- as.numeric(fit$theta)
        return(list(model = fit, 
                    coef = custom_coef, 
                    theta = fitted_theta, 
                    custom = TRUE))
      }
      
      list(model = fit, 
           coef = as.numeric(coef(fit)[1]), 
           theta = as.numeric(fit$theta), 
           custom = FALSE)
    }
  })
  
  # Display fitted model summary in Model Fit tab
  output$model_summary <- renderPrint({
    fit <- fitted_model()
    req(fit)
    
    cat("Intercept (log-scale):", round(fit$coef, 3), "\n")
    cat("Theta (dispersion):", round(fit$theta, 3), "\n\n")
    
    if (fit$custom) {
      cat("Custom intercept applied.\n\n")
    }
    
    if (!is.null(fit$model)) {
      print(summary(fit$model))
    } else {
      cat("No fitted model: custom abundance values were used.\n")
    }
  })
  
  # Run simulations
  sim_results <- eventReactive(input$run_btn, {
    
    # extract theta & int from fit
    fit <- fitted_model()
    intercept <- fit$coef 
    theta <- fit$theta
    
    # compute effect size
    eff <- if (input$trend_choice == "custom") {
      red <- as.numeric(input$custom_red)
      validate(need(!is.na(red) && red >= 0 && red < 1,
                    "Custom reduction must be a number >= 0 and < 1 
                    (e.g. 0.25 = 25%)."))
      log(1 - red)
    } else {
      red_preset <- as.numeric(input$trend_choice)
      log(1 - red_preset)
    }
    
    # choose sample time based on input
    time <- switch(input$years_sampled,
                   "1:3" = 1:3,
                   "1:5" = 1:5,
                   "c(1, 3, 5)" = c(1, 3, 5))
    
    # choose offset based on input
    offset_sim <- ifelse(input$gear == "DT", log(5), log(15))
    
    # loop over chosen sample range
    samp_seq <- seq(input$sample_range[1],
                    input$sample_range[2],
                    by = input$sample_by)
    
    powers <- numeric(length(samp_seq))
    
    # display a progress bar as sims run
    withProgress(message = "Running simulations...", {
      for (i in seq_along(samp_seq)) {
        incProgress(i / length(samp_seq),
                    detail = paste("Sample size =", samp_seq[i]))
        
        powers[i] <- simulated_power(
          sample_size = samp_seq[i],
          n_sim = input$n_sim,
          params = intercept,
          eff = eff,
          theta = theta,
          offset = offset_sim,
          time = time,
          alpha = input$alpha
        )
      }
    })
    
    df_out <- data.frame(SampleSize = samp_seq, Power = powers)
    
    thresholds <- c(0.6, 0.7, 0.8)
    mins <- sapply(thresholds, function(th) {
      valid_samples <- df_out$SampleSize[df_out$Power >= th]
      if(length(valid_samples) == 0) return(NA)  # Not achieved
      min(valid_samples)
    })
    names(mins) <- paste0("Power_", thresholds)
    
    list(df = df_out, 
         mins = mins,
         intercept = intercept, 
         theta = theta,
         eff = eff, 
         offset = offset_sim)
  })#, ignoreNULL = FALSE) # this command allows sim to run when opened
  
  # Display summary table in Summary tab
  output$summary_table <- renderDT({
    res <- sim_results()
    req(res)
    
    is_custom <- fitted_model()$custom
    
    mu_sample <- if (is_custom) {
      exp(res$intercept)
    } else {
      exp(res$intercept + res$offset)
    }
    
    var_sample <- mu_sample + (mu_sample ^ 2 / res$theta)
    sd_sample  <- sqrt(var_sample)
    
    mins <- res$mins
    red_text <- if (input$trend_choice == "custom") {
      paste0(round(input$custom_red * 100, 1), "% (custom)")
    } else {
      paste0(round(as.numeric(input$trend_choice) * 100, 1), "%")
    }
    
    power_labels <- c(
      "Minimum sample size for 60% power",
      "Minimum sample size for 70% power",
      "Minimum sample size for 80% power"
    )[seq_along(mins)] 
    
    datatable(
      tibble(
        Metric = c("Intercept (log)",
                   "Expected mean per sample (count scale)",
                   "SD per sample (count scale)",
                   "Theta (dispersion)",
                   "Effect size (log scale)",
                   "Reduction in abundance",
                   power_labels),
        Value = c(round(res$intercept, 3),
                  round(mu_sample, 3),
                  round(sd_sample, 3),
                  round(res$theta, 3),
                  round(res$eff, 3),
                  red_text,
                  ifelse(is.na(mins["Power_0.6"]), "Not achieved", 
                         as.character(mins["Power_0.6"])),
                  ifelse(is.na(mins["Power_0.7"]), "Not achieved", 
                         as.character(mins["Power_0.7"])),
                  ifelse(is.na(mins["Power_0.8"]), "Not achieved", 
                         as.character(mins["Power_0.8"]))
        )
      ),
      rownames = FALSE,
      options = list(dom = 't')
    )
  })
  
  # Display power by sample size plot in Power Curve tab
  output$power_plot <- renderPlotly({
    res <- sim_results()
    req(res)
    df <- res$df
    
    mins <- res$mins[!is.na(res$mins)]
    
    if (length(mins) == 0) {
      label_df <- data.frame(
        x = numeric(0),
        y = numeric(0),
        label = character(0)
      )
    } else {
      label_df <- data.frame(
        x = mins + 15,
        y = rep(0.4, length(mins)),
        label = paste("Minimum sample \nsize =", mins)
      )
    }
    
    trend_text <- if (input$trend_choice == "custom") {
      paste0("-", round(input$custom_red * 100, 1), "% (custom)")
    } else {
      paste0("-", round(as.numeric(input$trend_choice) * 100, 1), "%")
    }
    
    abund_text <- if (input$abundance == "Custom") {
      paste0(
        "Custom (",
        input$custom_abundance, " / min × ",
        input$custom_time_unit, " min)"
      )
    } else input$abundance
    
    title_full <- paste0("Power vs Sample Size (Gear: ", input$gear,
                         " | Abundance: ", abund_text,
                         " | Trend: ", trend_text, ")")
    
    p <- ggplot(df, aes(x = SampleSize, y = Power)) +
      geom_line(color = "#1f77b4", linewidth = 0.8) +
      geom_point(color = "#0b3d91", size = 2.5) + 
      geom_hline(yintercept = c(0.6, 0.7, 0.8), linetype = "dashed", 
                 color = "gray60", linewidth = 0.5) +
      geom_vline(data = data.frame(mins = mins), aes(xintercept = mins),
        color = "red", linetype = "dotted", linewidth = 0.8) +
      geom_text(data = label_df, aes(x = x, y = y, label = label),
                color = "red", vjust = -0.5, hjust = 0, size = 3) + 
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
      labs(
        title = title_full,
        x = "Sample Size per Year",
        y = "Estimated Power"
      ) +
      theme_classic() + 
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5, 
                                  color = "#2c3e50", size = 18),
        axis.title = element_text(face = "bold", size = 14),
        axis.text = element_text(color = "black", size = 12),
        panel.grid = element_blank() 
      )
    ggplotly(p) %>%
      layout(
        title = list(
          text = paste0("Power vs Sample Size \n(Gear: ", input$gear,
                        " | Abundance: ", input$abundance, ")")
        ),
        xaxis = list(title = "Sample Size per Year"),
        yaxis = list(title = "Estimated Power")
      )
  })

# Display data table in Data Summary tab
output$powerTable <- renderDT({
  res <- sim_results()
  req(res)
  
  mins <- res$mins
  df <- res$df
  
  # If no sample size reaches power — show warning message
  if (all(is.na(mins))) {
    DT::datatable(
      data.frame(Message = "⚠️ Target power not achieved with selected settings."),
      options = list(dom = 't'),
      rownames = FALSE
    )
  } 
  
  else {
    min_80 <- if (any(df$Power >= 0.8)) min(df$SampleSize[df$Power >= 0.8]) 
    else NA
    
    results_table <- data.frame(
      Gear = input$gear,
      Trend = ifelse(input$trend_choice == "custom",
                     paste0(round(input$custom_red * 100), "% Custom Reduction"),
                     paste0(as.numeric(input$trend_choice) * 100, "% Reduction")),
      Abundance = input$abundance,
      Minimum_Samples = min_80,
      Achieves_80 = ifelse(min_80 <= max(df$SampleSize), 
                              "<span style='color:green;font-weight:bold;'>Yes</span>", 
                              "No")
    )
    
    colnames(results_table)[colnames(results_table) == "Minimum_Samples"] <- "Minimum Number of Samples per Site"
    colnames(results_table)[colnames(results_table) == "Achieves_80"] <- "Achieves 80% Power?"
    
    DT::datatable(
      results_table,
      escape = FALSE,
      rownames = FALSE,
      options = list(dom = 't')
    )
  }
})

output$custom_diagnostic <- renderUI({
  mu  <- input$custom_abundance * input$custom_time_unit
  sd  <- input$custom_sd
  var <- sd ^ 2
  
  msgs <- c()
  
  if (mu < 0.01)
    msgs <- c(msgs, "⚠️ Expected abundance is very low; model power estimates may be unstable.")
  if (mu > 50)
    msgs <- c(msgs, "⚠️ Expected abundance is very high; simulated counts may exceed realistic biological levels.")
  if (var <= mu)
    msgs <- c(msgs, "⚠️ Variance must exceed mean for the negative binomial model.")
  
  if (length(msgs) == 0) return(NULL)
  
  tagList(
    strong("Diagnostics:"),
    lapply(msgs, function(m) p(style = "color:#b30000;", m))
  )
})

}

shinyApp(ui, server)