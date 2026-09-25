library(shiny)
library(deSolve)
library(ggplot2)
library(tidyverse)
library(bslib)
library(plotly)

# Explicit, deterministic load order (does not rely on Shiny's implicit
# alphabetical R/-directory auto-sourcing, though that would also work here
# since Extended_Model_Population.r happens to sort before the other two).
# popgen_consts.R/popgen_port.R/estimated_parameters.csv must live at the
# app root, not in R/: Extended_Model_Population.r's own internal source()/
# read.csv() calls use paths relative to the working directory, not to its
# own file location, and Shiny keeps the working directory at the app root.
source("popgen_consts.R")
source("popgen_port.R")
source("R/Extended_Model_Population.r")
source("R/Forward_Dosimetry_Stochastic.r")
source("R/Reverse_Dosimetry_Stochastic.r")

# ── UI ────────────────────────────────────────────────────────────────────────

# Shared population-generation inputs (forward & reverse tabs each get their
# own copy with a distinct id prefix, mirroring the original app's rev_*
# convention) -- see Extended_Model_Population.r's .generate_population()
# for what each one controls.
population_inputs <- function(prefix) {
  tagList(
    h4("Virtual Population"),
    numericInput(paste0(prefix, "_n"), "Population size (n)", value = 20, min = 1, max = 10000),
    helpText("Larger populations take proportionally longer to simulate."),
    selectInput(paste0(prefix, "_dataset"), "Reference dataset",
                choices = c("ICRP", "P3M")),
    selectInput(paste0(prefix, "_population_type"), "Population variability",
                choices = c("Realistic", "HighVariation")),
    sliderInput(paste0(prefix, "_age_range"), "Age range (years)",
                min = 18, max = 90, value = c(30, 60)),
    sliderInput(paste0(prefix, "_bmi_range"), "BMI range (kg/m2)",
                min = 18.5, max = 40, value = c(20, 28)),
    sliderInput(paste0(prefix, "_height_range"), "Height range (cm)",
                min = 140, max = 210, value = c(155, 195)),
    sliderInput(paste0(prefix, "_prob_male"), "Probability male", min = 0, max = 1,
                value = 0.5, step = 0.05),
    h5("Ethnicity probabilities"),
    helpText("Must sum to ≤ 1; the remainder is treated as \"Other\"."),
    fluidRow(
      column(4, numericInput(paste0(prefix, "_eth_white"), "White", value = 1, min = 0, max = 1, step = 0.05)),
      column(4, numericInput(paste0(prefix, "_eth_black"), "Black", value = 0, min = 0, max = 1, step = 0.05)),
      column(4, numericInput(paste0(prefix, "_eth_nbh"), "Non-Black Hispanic", value = 0, min = 0, max = 1, step = 0.05))
    ),
    numericInput(paste0(prefix, "_seed"), "Random seed", value = 123, min = 1, step = 1),
    helpText("Same seed creates the same virtual population each run; change it for a fresh draw.")
  )
}

ui <- navbarPage("PFAS PBK Model — Population",
theme = bs_theme(bootswatch = "minty"),

  tabPanel("Forward Dosimetry",

    sidebarLayout(
      sidebarPanel(
        h4("Simulation Parameters"),
        numericInput("duration", "Duration", value = 24),
        selectInput("time_scale", "Time scale of simulation",
                    choices = c("minutes", "hours", "days", "weeks", "months", "years")),
        numericInput("time_step", "Time step", value = 0.1, min = 0.01, step = 0.01),
        selectInput("chemical", "PFAS",
                    choices = c("PFHpA", "PFOA", "PFNA", "PFDA", "PFBS",
                                "PFHxS", "PFOS", "DONA", "HFPO_DA", "PFBA", "PFHxA")),
        selectInput("exp_type", "Exposure type",
                    choices = c("continuous", "pharmacokinetics")),
        helpText("Continuous: steady daily per-kg-BW intake rate. ",
                 "Pharmacokinetics: one-off per-kg-BW dose(s), absorbed from the GI tract."),
        numericInput("n_ingestion", "Number of ingestion steps", value = 1, min = 1),
        helpText("Each step defines a daily intake rate (ng/kg BW/day) starting at a given time; ",
                 "for \"pharmacokinetics\", each step is instead a one-off dose given at that time."),
        uiOutput("ingestion_inputs"),
        hr(),
        population_inputs("fwd"),
        br(),
        actionButton("run", "Run Simulation", class = "btn-primary")
      ),
      mainPanel(
        downloadButton("dl_fwd", "Download CSV (all individuals)"),
        br(), br(),
        h5("Population summary"),
        tableOutput("fwd_pop_summary"),
        br(),
        plotlyOutput("fwd_Cplasma_plot", height = "400px"),
        plotlyOutput("fwd_allcomp_plot", height = "500px"),
        plotlyOutput("fwd_elimination_plot", height = "350px"),
        plotlyOutput("fwd_ingestion_plot", height = "300px")
      )
    )
  ),

  tabPanel("Reverse Dosimetry",
    sidebarLayout(
      sidebarPanel(
        h4("Simulation Parameters"),
        numericInput("rev_duration", "Duration", value = 5),
        selectInput("rev_time_scale", "Time scale of simulation",
                    choices = c("minutes", "hours", "days", "weeks", "months", "years"),
                    selected = "years"),
        selectInput("rev_chemical", "PFAS",
                    choices = c("PFHpA", "PFOA", "PFNA", "PFDA", "PFBS",
                                "PFHxS", "PFOS", "DONA", "HFPO_DA", "PFBA", "PFHxA")),
        selectInput("rev_exp_type", "Exposure type",
                    choices = c("continuous", "pharmacokinetics")),
        helpText("Continuous: search for the steady per-kg-BW intake rate that reproduces the POD. ",
                 "Pharmacokinetics: search for a one-off per-kg-BW dose instead."),
        numericInput("rev_ingestion_time", "Dose/rate start time", value = 0, min = 0),
        hr(),
        h4("Point of Departure"),
        numericInput("rev_POD", "POD (μg/L)", value = 1, min = 0),
        helpText("Target tissue concentration to back-calculate exposure from, per individual."),
        numericInput("rev_POD_sd", "POD standard deviation (μg/L)", value = 0, min = 0),
        helpText("Measurement/estimate uncertainty in the POD. Leave at 0 for an exact POD ",
                 "(same target for every individual); above 0, each individual instead searches ",
                 "against its own random POD draw (lognormal, mean = POD above, this SD), so the ",
                 "exposure distribution reflects POD uncertainty as well as physiology."),
        selectInput("rev_compartment", "Tissue",
                    choices = c("Serum"   = "serum",
                                "Liver"   = "liver",
                                "Adipose" = "adipose",
                                "Brain"   = "brain",
                                "Gonads"  = "gonads",
                                "Gut"     = "gut",
                                "Heart"   = "heart",
                                "Lung"    = "lung",
                                "Muscle"  = "muscle",
                                "Skin"    = "skin",
                                "Kidney"  = "kidney")),
        hr(),
        population_inputs("rev"),
        helpText(strong("Runtime note:"), " each individual runs a full numerical search ",
                 "(up to ~1000 ODE solves), so this is much slower per-individual than Forward ",
                 "Dosimetry. Start with n ≤ 15 here."),
        br(),
        actionButton("rev_run", "Run Reverse Dosimetry", class = "btn-primary")
      ),
      mainPanel(
        downloadButton("dl_rev", "Download CSV (per-individual estimates)"),
        br(), br(),
        h5("Population summary"),
        tableOutput("rev_pop_summary"),
        br(),
        h5("Exposure estimate distribution"),
        tableOutput("rev_results_table"),
        plotlyOutput("rev_exposure_density", height = "400px"),
        br(),
        h5("Per-individual estimates"),
        tableOutput("rev_individual_table")
      )
    )
  ),

  tabPanel("Help",
    fluidPage(
      fluidRow(
        column(8, offset = 2,

          h3("How to Use This App"),
          p("This app runs a population (virtual-individual) version of the physiologically
            based kinetic (PBK) model for PFAS compounds: instead of one fixed body weight,
            each run generates ", strong("n"), " virtual individuals with popgen-derived
            organ mass/flow, sex-specific hematocrit, and age-corrected renal clearance
            (see the project's ", code("Extended_Model_Population.r"), "), and solves the
            model once per individual. Forward Dosimetry simulates the resulting
            ", em("distribution"), " of tissue concentrations from a known exposure;
            Reverse Dosimetry estimates the ", em("distribution"), " of exposures that
            would produce a target tissue concentration, one exposure estimate per
            individual."),

          hr(),
          h3("Virtual Population"),
          tags$dl(
            tags$dt("Population size (n)"),
            tags$dd("Number of virtual individuals to simulate."),
            tags$dt("Reference dataset"),
            tags$dd("\"ICRP\" or \"P3M\" -- which reference anthropometric/physiological
                     dataset popgen draws individuals from."),
            tags$dt("Population variability"),
            tags$dd("\"Realistic\" draws organ mass/flow from literature reference
                     distributions; \"HighVariation\" instead samples more broadly across
                     the given ranges, for stress-testing rather than a realistic cohort."),
            tags$dt("Age / BMI / Height range"),
            tags$dd("Bounds each individual's age, BMI, and height are drawn within
                     (height range only weakly constrains \"Realistic\" draws -- see the
                     project notes on popgen's own sampling behaviour)."),
            tags$dt("Probability male"),
            tags$dd("Chance any given individual is male rather than female; sex affects
                     hematocrit, kidney mass, and several organ volumes."),
            tags$dt("Ethnicity probabilities"),
            tags$dd("Probability of White / Black / Non-Black Hispanic; must sum to
                     ≤ 1, with the remainder treated as \"Other\"."),
            tags$dt("Random seed"),
            tags$dd("Fixes the population draw for reproducibility. Change it to sample a
                     different virtual population under the same settings.")
          ),

          hr(),
          h3("Forward Dosimetry"),

          h4("Simulation Parameters"),
          tags$dl(
            tags$dt("Duration & Time scale"),
            tags$dd("Total length of the simulation and its unit (e.g. 24 hours,
                     30 days, 1 year). The time step controls output resolution."),
            tags$dt("PFAS compound"),
            tags$dd("Selects the chemical-specific parameters (e.g. protein
                     binding, renal clearance, half-life) used in the model.")
          ),

          h4("Exposure Type"),
          tags$dl(
            tags$dt("Continuous"),
            tags$dd("The daily per-kg-BW ingestion rate is applied as a constant, ongoing
                     input -- each individual's own body weight converts it to an
                     absolute intake. Appropriate for chronic dietary exposure scenarios."),
            tags$dt("Pharmacokinetics"),
            tags$dd("Each ingestion step is instead a one-off per-kg-BW dose absorbed from
                     the GI tract. Use for single- or repeated-dose pharmacokinetic studies.")
          ),

          h4("Multiple Exposure Steps"),
          p("Define multiple ingestion steps to model changing exposure over time: e.g. a
            non-zero rate starting at time 0, followed by a second step with rate 0
            starting at the switchover time."),

          h4("Outputs"),
          tags$ul(
            tags$li(tags$b("Population summary:"), " how many individuals were generated
                    (and discarded by popgen's own rejection rules)."),
            tags$li(tags$b("Plasma concentration:"), " median trajectory with a
                    2.5th-97.5th percentile band across the population (µg/L)."),
            tags$li(tags$b("All compartments:"), " median concentration per tissue over
                    time, on a log scale."),
            tags$li(tags$b("Cumulative elimination:"), " median + percentile band for
                    total PFAS excreted via urine and faeces (µg)."),
            tags$li(tags$b("Ingestion rate:"), " median + percentile band for the
                    (BW-scaled, hence individual-specific) absolute ingestion applied to
                    the model (µg/time).")
          ),
          p("The downloaded CSV contains every individual's full time-course, tagged by
            an \"id\" column."),

          hr(),
          h3("Reverse Dosimetry"),

          h4("Simulation Parameters"),
          tags$dl(
            tags$dt("Duration & Time scale"),
            tags$dd("The simulation is run until the end of this period; the estimated
                     exposure is the per-kg-BW rate or dose that produces the target
                     tissue concentration at the end of the duration, for each
                     individual."),
            tags$dt("PFAS compound"),
            tags$dd("Selects the chemical-specific parameters used in the model."),
            tags$dt("Exposure type"),
            tags$dd("\"Continuous\" searches for a steady per-kg-BW intake rate;
                     \"pharmacokinetics\" searches for a one-off per-kg-BW dose instead.
                     Both are BW-scaled per individual.")
          ),

          h4("Point of Departure"),
          tags$dl(
            tags$dt("POD (μg/L)"),
            tags$dd("The target internal concentration (point of departure) in
                     the selected tissue."),
            tags$dt("POD standard deviation (μg/L)"),
            tags$dd("Measurement/estimate uncertainty in the POD (e.g. the SD across replicate
                     study measurements). At 0 (default), every individual is fit to the exact
                     same POD, so only physiology drives the exposure distribution. Above 0,
                     each individual instead searches against its own random POD draw
                     (lognormal, mean = POD, this SD), so the resulting exposure distribution
                     reflects both physiology AND POD uncertainty."),
            tags$dt("Tissue"),
            tags$dd("The compartment in which the POD is defined.")
          ),

          h4("Outputs"),
          tags$ul(
            tags$li(tags$b("Population summary:"), " how many individuals were generated."),
            tags$li(tags$b("Exposure estimate distribution:"), " mean/SD/percentiles of
                    the per-kg-BW exposure needed to reach the POD, across the
                    population, plus a density plot (with a rug showing each
                    individual's own estimate)."),
            tags$li(tags$b("Per-individual estimates:"), " each individual's own POD
                    target (identical for everyone unless POD SD > 0), estimated
                    exposure, and the optimizer's relative error (how exactly it hit
                    the target).")
          ),
          p("Unlike Forward Dosimetry, this tab does not plot concentration-time curves:
            the underlying search only records each individual's final estimated exposure
            and how well it matched the POD, not the full time-course at that exposure.
            The downloaded CSV contains one row per individual."),
          p("Results can be downloaded as a CSV file using the Download button.")
        )
      )
    )
  )
)

# ── Shared plot theme ─────────────────────────────────────────────────────────
theme_pfas <- function(base_size = 13) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title       = element_text(face = "bold", size = base_size + 1,
                                      margin = margin(b = 10)),
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey92"),
      axis.title       = element_text(colour = "grey30"),
      axis.text        = element_text(colour = "grey40"),
      legend.position  = "right",
      legend.key.size  = unit(0.5, "cm"),
      plot.margin      = margin(12, 12, 12, 12)
    )
}

# Percentile summary (2.5/50/97.5) of one column, grouped by time -- the
# population analogue of a single trajectory.
percentile_summary <- function(df, value_col) {
  df %>%
    group_by(time) %>%
    summarise(
      p2.5  = quantile(.data[[value_col]], 0.025, na.rm = TRUE),
      p50   = quantile(.data[[value_col]], 0.5,   na.rm = TRUE),
      p97.5 = quantile(.data[[value_col]], 0.975, na.rm = TRUE),
      .groups = "drop"
    )
}

plot_percentile_band <- function(summary_df, y_lab, title, colour = "#2E9E7A") {
  ggplot(summary_df, aes(x = time)) +
    geom_ribbon(aes(ymin = p2.5, ymax = p97.5), fill = colour, alpha = 0.15) +
    geom_line(aes(y = p50), colour = colour, linewidth = 1.1) +
    labs(y = y_lab, title = title) +
    theme_pfas()
}

# Read a set of "<prefix>_<name>" inputs into the population-generation
# fields forward_dosimetry_stochastic()/reverse_dosimetry_stochastic() need.
read_population_inputs <- function(input, prefix) {
  list(
    n               = input[[paste0(prefix, "_n")]],
    dataset         = input[[paste0(prefix, "_dataset")]],
    population_type = input[[paste0(prefix, "_population_type")]],
    age_range       = input[[paste0(prefix, "_age_range")]],
    bmi_range       = input[[paste0(prefix, "_bmi_range")]],
    height_range    = input[[paste0(prefix, "_height_range")]],
    prob_male       = input[[paste0(prefix, "_prob_male")]],
    ethnicity_probs = c(input[[paste0(prefix, "_eth_white")]],
                        input[[paste0(prefix, "_eth_black")]],
                        input[[paste0(prefix, "_eth_nbh")]]),
    seed            = input[[paste0(prefix, "_seed")]]
  )
}

# ── Server ────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  output$ingestion_inputs <- renderUI({
    n <- input$n_ingestion
    dose_label <- if (identical(input$exp_type, "pharmacokinetics")) "Dose" else "Ingestion"
    lapply(1:n, function(i) {
      fluidRow(
        column(6, numericInput(paste0("ingestion_", i),
                               paste0(dose_label, " ", i, " (ng/kg BW)"), value = 10)),
        column(6, numericInput(paste0("ingestion_time_", i),
                               paste("Start time", i), value = i - 1))
      )
    })
  })

  fwd_data <- eventReactive(input$run, {
    n <- input$n_ingestion
    ingestion_vec      <- sapply(1:n, function(i) input[[paste0("ingestion_", i)]])
    ingestion_time_vec <- sapply(1:n, function(i) input[[paste0("ingestion_time_", i)]])

    pop_in <- read_population_inputs(input, "fwd")

    do.call(forward_dosimetry_stochastic, c(pop_in, list(
      duration       = input$duration,
      time_step      = input$time_step,
      chemical       = input$chemical,
      ingestion      = ingestion_vec,
      ingestion_time = ingestion_time_vec,
      exp_type       = input$exp_type,
      time_scale     = input$time_scale
    )))
  })

  output$fwd_pop_summary <- renderTable({
    pop <- fwd_data()$population
    data.frame(
      "Individuals simulated" = nrow(pop),
      "Mean age"    = round(mean(pop$age), 1),
      "% male"      = round(100 * mean(pop$sex == "Male"), 1),
      "Mean BW (kg)" = round(mean(pop$body_mass_kg), 1),
      check.names = FALSE
    )
  })

  output$fwd_Cplasma_plot <- renderPlotly({
    df <- fwd_data()$results
    p <- plot_percentile_band(
      percentile_summary(df, "CA"),
      y_lab = "Plasma concentration (µg/L)",
      title = paste0(input$chemical, " — Plasma concentration vs time (median ± 95% band)")
    ) + labs(x = paste0("Time (", input$time_scale, ")"))
    ggplotly(p) %>% layout(hovermode = "x unified")
  })

  output$fwd_allcomp_plot <- renderPlotly({
    df <- fwd_data()$results
    compartments <- c("CR", "CAdi", "CBra", "CGon", "CHea", "CLun",
                      "CMus", "CSki", "CSpl", "Cpan", "CKb", "Cfil",
                      "CL", "CGI", "CA")
    df_long <- df %>%
      select(time, all_of(compartments)) %>%
      pivot_longer(-time, names_to = "compartment", values_to = "concentration") %>%
      group_by(time, compartment) %>%
      summarise(median_conc = median(abs(concentration), na.rm = TRUE), .groups = "drop")

    p <- ggplot(df_long, aes(x = time, y = median_conc, colour = compartment)) +
      geom_line(linewidth = 0.8, alpha = 0.85) +
      scale_y_log10(labels = scales::label_scientific()) +
      scale_colour_viridis_d(option = "turbo") +
      labs(
        x      = paste0("Time (", input$time_scale, ")"),
        y      = "Median concentration (µg/L, log scale)",
        colour = NULL,
        title  = "All compartment concentrations (population median)"
      ) +
      theme_pfas()
    ggplotly(p)
  })

  output$fwd_elimination_plot <- renderPlotly({
    df <- fwd_data()$results
    aurine_summary <- percentile_summary(df, "Aurine") %>% mutate(route = "Aurine")
    afeces_summary <- percentile_summary(df, "Afeces") %>% mutate(route = "Afeces")
    df_long <- bind_rows(aurine_summary, afeces_summary)

    p <- ggplot(df_long, aes(x = time, y = p50, colour = route, fill = route)) +
      geom_ribbon(aes(ymin = p2.5, ymax = p97.5), alpha = 0.12, colour = NA) +
      geom_line(linewidth = 1) +
      scale_colour_manual(values = c(Aurine = "#3A86FF", Afeces = "#FF6B6B")) +
      scale_fill_manual(values  = c(Aurine = "#3A86FF", Afeces = "#FF6B6B")) +
      labs(
        x      = paste0("Time (", input$time_scale, ")"),
        y      = "Cumulative amount (µg)",
        colour = NULL, fill = NULL,
        title  = "Cumulative elimination (median ± 95% band)"
      ) +
      theme_pfas()
    ggplotly(p) %>% layout(hovermode = "x unified")
  })

  output$fwd_ingestion_plot <- renderPlotly({
    df <- fwd_data()$results
    p <- plot_percentile_band(
      percentile_summary(df, "ingestion"),
      y_lab = "Ingestion rate (µg/time, absolute, individual-specific)",
      title = "Ingestion rate over time (median ± 95% band)",
      colour = "#E07B39"
    ) + labs(x = paste0("Time (", input$time_scale, ")"))
    ggplotly(p) %>% layout(hovermode = "x unified")
  })

  output$dl_fwd <- downloadHandler(
    filename = function() paste0("forward_population_", Sys.Date(), ".csv"),
    content  = function(file) write.csv(fwd_data()$results, file, row.names = FALSE)
  )

  rev_data <- eventReactive(input$rev_run, {
    pop_in <- read_population_inputs(input, "rev")

    do.call(reverse_dosimetry_stochastic, c(pop_in, list(
      chemical       = input$rev_chemical,
      duration       = input$rev_duration,
      time_scale     = input$rev_time_scale,
      POD            = input$rev_POD,
      POD_sd         = input$rev_POD_sd,
      compartment    = input$rev_compartment,
      exp_type       = input$rev_exp_type,
      ingestion_time = input$rev_ingestion_time
    )))
  })

  output$rev_pop_summary <- renderTable({
    pop <- rev_data()$population
    data.frame(
      "Individuals simulated" = nrow(pop),
      "Mean age"    = round(mean(pop$age), 1),
      "% male"      = round(100 * mean(pop$sex == "Male"), 1),
      "Mean BW (kg)" = round(mean(pop$body_mass_kg), 1),
      check.names = FALSE
    )
  })

  output$rev_results_table <- renderTable({
    est <- rev_data()$exposure_estimates
    exposure_label <- if (identical(input$rev_exp_type, "pharmacokinetics"))
      "Estimated dose (ng/kg BW)" else "Estimated intake rate (ng/kg BW/time)"
    df <- data.frame(
      Chemical = input$rev_chemical,
      Tissue   = input$rev_compartment,
      POD      = input$rev_POD,
      "POD SD" = input$rev_POD_sd,
      Duration = paste(input$rev_duration, input$rev_time_scale),
      Mean     = round(mean(est$exposure), 4),
      SD       = round(sd(est$exposure), 4),
      P2.5     = round(unname(quantile(est$exposure, 0.025)), 4),
      Median   = round(median(est$exposure), 4),
      P97.5    = round(unname(quantile(est$exposure, 0.975)), 4),
      check.names = FALSE
    )
    names(df)[names(df) == "Mean"] <- paste0(exposure_label, " — mean")
    df
  })

  output$rev_exposure_density <- renderPlotly({
    est <- rev_data()$exposure_estimates
    # Label uses input$rev_POD/rev_POD_sd (the target actually asked for), not
    # est$POD[1] -- with POD_sd > 0 that column holds each individual's own
    # random draw, not the shared target.
    pod_label <- if (input$rev_POD_sd > 0)
      sprintf("%.3g ± %.3g (mean ± sd)", input$rev_POD, input$rev_POD_sd)
      else as.character(input$rev_POD)
    p <- ggplot(est, aes(x = exposure)) +
      geom_density(fill = "#2a78d6", colour = "#2a78d6", alpha = 0.45, linewidth = 1) +
      geom_rug(alpha = 0.6, colour = "#2a78d6") +
      geom_vline(aes(xintercept = median(exposure)), color = "#0b0b0b",
                 linetype = "dashed", linewidth = 0.8) +
      labs(
        x = paste0("Estimated exposure reproducing POD (", est$compartment[1],
                    " = ", pod_label, " at t = ", est$time[1], " ", input$rev_time_scale, ")"),
        y = "Density",
        title = paste0("Stochastic reverse dosimetry — ", input$rev_chemical,
                        " (n = ", nrow(rev_data()$population), ")")
      ) +
      theme_pfas()
    ggplotly(p)
  })

  output$rev_individual_table <- renderTable({
    est <- rev_data()$exposure_estimates
    pop <- rev_data()$population %>% select(id, age, sex, body_mass_kg)
    est %>%
      left_join(pop, by = "id") %>%
      transmute(
        id, age = round(age, 1), sex, "BW (kg)" = round(body_mass_kg, 1),
        "POD target" = round(POD, 4),
        exposure = round(exposure, 4), rel_error = signif(rel_error, 4)
      )
  })

  output$dl_rev <- downloadHandler(
    filename = function() paste0("reverse_population_", Sys.Date(), ".csv"),
    content  = function(file) write.csv(rev_data()$exposure_estimates, file, row.names = FALSE)
  )
}

shinyApp(ui = ui, server = server)
