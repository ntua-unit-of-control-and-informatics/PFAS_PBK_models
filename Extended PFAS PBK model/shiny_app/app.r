library(shiny)
library(deSolve)
library(ggplot2)
library(tidyverse)
library(bslib)

source("R/Extended_Model.r")

# ── UI ────────────────────────────────────────────────────────────────────────
ui <- navbarPage("PFAS PBK Model",
theme = bs_theme(bootswatch = "minty"),

  tabPanel("Forward Dosimetry",

    sidebarLayout(
      sidebarPanel(
        h4("Simulation Parameters"),
        numericInput("BW", "Bodyweight (kg)", value = 70, min = 1),
        numericInput("duration", "Duration", value = 24),
        selectInput("time_scale", "Time scale of simulation",
                    choices = c("minutes", "hours", "days", "weeks", "months", "years")),
        numericInput("time_step", "Time step", value = 0.1, min = 0.01, step = 0.01),
        selectInput("chemical", "PFAS",
                    choices = c("PFHpA", "PFOA", "PFNA", "PFDA", "PFBS",
                                "PFHxS", "PFOS", "DONA", "HFPO_DA", "PFBA", "PFHxA")),
        selectInput("admin_type", "Administration type",
                    choices = c("Intravenous" = "iv", "Oral" = "oral")),
        conditionalPanel(
          condition = "input.admin_type == 'iv'",
          numericInput("n_admin", "Number of IV doses", value = 1, min = 0),
          helpText("Each dose is defined by an amount (μg) and a time point."),
          uiOutput("admin_inputs")
        ),
        conditionalPanel(
          condition = "input.admin_type == 'oral'",
          numericInput("n_ingestion", "Number of ingestion steps", value = 1, min = 0),
          helpText("Each step defines a daily intake rate starting at a given time."),
          uiOutput("ingestion_inputs"),
          selectInput("exp_type", "Exposure type",
                      choices = c("continuous", "pharmacokinetics")),
          helpText("Continuous: steady daily intake. Pharmacokinetics: single-dose absorption.")
        ),
        br(),
        actionButton("run", "Run Simulation", class = "btn-primary")
      ),
      mainPanel(
        downloadButton("dl_fwd", "Download CSV"),
        br(), br(),
        plotOutput("fwd_Cplasma_plot", height = "400px"),
        plotOutput("fwd_allcomp_plot", height = "500px"),
        plotOutput("fwd_elimination_plot", height = "350px"),
        conditionalPanel(
          condition = "input.admin_type == 'oral'",
          plotOutput("fwd_ingestion_plot", height = "300px")
        )
      )
    )
  ),

  tabPanel("Reverse Dosimetry",
    sidebarLayout(
      sidebarPanel(
        h4("Simulation Parameters"),
        numericInput("rev_BW", "Bodyweight (kg)", value = 70, min = 1),
        numericInput("rev_duration", "Duration", value = 24),
        selectInput("rev_time_scale", "Time scale of simulation",
                    choices = c("minutes", "hours", "days", "weeks", "months", "years"),
                    selected = "years"),
        numericInput("rev_time_step", "Time step", value = 0.1, min = 0.01, step = 0.1),
        selectInput("rev_chemical", "PFAS",
                    choices = c("PFHpA", "PFOA", "PFNA", "PFDA", "PFBS",
                                "PFHxS", "PFOS", "DONA", "HFPO_DA", "PFBA", "PFHxA")),
        selectInput("rev_admin_type", "Administration type",
                    choices = c("Oral" = "oral", "Intravenous" = "iv")),
        conditionalPanel(
          condition = "input.rev_admin_type == 'iv'",
          numericInput("rev_admin_time", "Dose time", value = 0, min = 0),
          helpText("Time point at which the IV dose is administered.")
        ),
        hr(),
        h4("Point of Departure"),
        numericInput("rev_POD", "POD (μg/L)", value = 1, min = 0),
        helpText("Target tissue concentration to back-calculate exposure from."),
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
        br(),
        actionButton("rev_run", "Run Reverse Dosimetry", class = "btn-primary")
      ),
      mainPanel(
        downloadButton("dl_rev", "Download CSV"),
        br(), br(),
        tableOutput("rev_results_table"),
        br(),
        plotOutput("rev_Cplasma_plot", height = "400px"),
        plotOutput("rev_allcomp_plot", height = "500px"),
        plotOutput("rev_elimination_plot", height = "350px"),
        conditionalPanel(
          condition = "input.rev_compartment != 'serum'",
          plotOutput("rev_tissue_plot", height = "350px")
        )
      )
    )
  ),

  tabPanel("Help",
    fluidPage(
      fluidRow(
        column(8, offset = 2,

          h3("How to Use This App"),
          p("This app runs a physiologically based kinetic (PBK) model for PFAS
            compounds. Two modes are available: Forward Dosimetry simulates tissue
            concentrations from a known exposure; Reverse Dosimetry estimates the
            exposure that would produce a target tissue concentration."),

          hr(),
          h3("Forward Dosimetry"),

          h4("Simulation Parameters"),
          tags$dl(
            tags$dt("Bodyweight (kg)"),
            tags$dd("Subject bodyweight used to scale physiological parameters
                     and to convert ingestion rates from ng/kg/day to absolute
                     amounts."),
            tags$dt("Duration & Time scale"),
            tags$dd("Total length of the simulation and its unit (e.g. 24 hours,
                     30 days, 1 year). The time step controls output resolution —
                     smaller values give smoother curves but take longer."),
            tags$dt("PFAS compound"),
            tags$dd("Selects the chemical-specific parameters (e.g. protein
                     binding, renal clearance, half-life) used in the model.")
          ),

          h4("Administration Type"),
          tags$dl(
            tags$dt("Intravenous (IV)"),
            tags$dd("One or more bolus doses injected directly into the bloodstream.
                     Each dose is defined by an amount (μg) and the time at which
                     it is given (in the selected time unit). Use this route to
                     simulate pharmacokinetic studies or to derive bioavailability
                     by comparing with oral results."),
            tags$dt("Oral"),
            tags$dd("Dietary or drinking-water ingestion expressed as ng/kg body
                     weight per day. Multiple ingestion steps allow the intake rate
                     to change over time — each step defines a new daily rate that
                     starts at the specified time and continues until the next step
                     (or the end of the simulation).")
          ),

          h4("Exposure Types (Oral only)"),
          tags$dl(
            tags$dt("Continuous"),
            tags$dd("The daily ingestion rate is applied as a constant, ongoing
                     input. Appropriate for chronic dietary exposure scenarios."),
            tags$dt("Pharmacokinetics"),
            tags$dd("The ingestion amount is administered as a discrete bolus
                     absorbed from the GI tract. Use for single- or repeated-dose
                     pharmacokinetic studies.")
          ),

          h4("Multiple Exposure Steps"),
          p("Both IV and oral routes support multiple steps to model changing
            exposure over time. For example, define two oral steps to simulate
            a high-exposure period followed by a clean period: the first with a
            non-zero intake starting at time 0, and the second with intake 0
            starting at the switchover time."),

          h4("Outputs"),
          tags$ul(
            tags$li(tags$b("Plasma concentration:"), " blood concentration over
                    time (μg/L)."),
            tags$li(tags$b("All compartments:"), " concentrations in all modelled
                    tissues on a log scale."),
            tags$li(tags$b("Cumulative elimination:"), " total PFAS excreted via
                    urine and faeces (μg)."),
            tags$li(tags$b("Ingestion rate:"), " (oral only) the ingestion input
                    function applied to the model.")
          ),
          p("Results can be downloaded as a CSV file using the Download button."),

          hr(),
          h3("Reverse Dosimetry"),

          h4("Simulation Parameters"),
          tags$dl(
            tags$dt("Bodyweight (kg)"),
            tags$dd("Subject bodyweight, used to scale physiological parameters
                     and convert the estimated exposure back to ng/kg/day."),
            tags$dt("Duration & Time scale"),
            tags$dd("The simulation is run until the end of this period. The
                     estimated exposure is the constant intake that produces the
                     target tissue concentration at the end of the duration.
                     The time step controls the resolution of the output plots."),
            tags$dt("PFAS compound"),
            tags$dd("Selects the chemical-specific parameters used in the model.")
          ),

          h4("Administration Type"),
          tags$dl(
            tags$dt("Oral"),
            tags$dd("The optimizer searches for the constant daily oral intake
                     (ng/kg/day) that results in the target tissue concentration
                     at the end of the simulation duration."),
            tags$dt("Intravenous (IV)"),
            tags$dd("The optimizer searches for the single IV bolus dose (μg)
                     administered at the specified dose time that results in the
                     target tissue concentration at the end of the simulation.")
          ),

          h4("Point of Departure"),
          tags$dl(
            tags$dt("POD (μg/L)"),
            tags$dd("The target internal concentration (point of departure) in
                     the selected tissue. This is the value the model is
                     back-calculated from."),
            tags$dt("Tissue"),
            tags$dd("The compartment in which the POD is defined. Supported
                     tissues: serum, liver, adipose, brain, gonads, gut, heart,
                     lung, muscle, skin, kidney.")
          ),

          h4("Outputs"),
          tags$ul(
            tags$li(tags$b("Results table:"), " estimated exposure (ng/kg/day
                    for oral, μg for IV) and relative error between the achieved
                    and target concentration."),
            tags$li(tags$b("Plasma concentration:"), " blood concentration over
                    time under the estimated exposure."),
            tags$li(tags$b("All compartments:"), " concentrations in all modelled
                    tissues on a log scale."),
            tags$li(tags$b("Cumulative elimination:"), " total PFAS excreted via
                    urine and faeces (μg)."),
            tags$li(tags$b("Target tissue:"), " (when tissue is not serum)
                    dedicated concentration-time plot for the selected tissue.")
          ),
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

# ── Server ────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  output$admin_inputs <- renderUI({
    n <- input$n_admin
    lapply(1:n, function(i) {
      fluidRow(
        column(6, numericInput(paste0("admin_dose_", i),
                               paste("Dose", i, "(ug)"), value = 0)),
        column(6, numericInput(paste0("admin_time_", i),
                               paste("Dose time", i), value = i - 1))
      )
    })
  })

  output$ingestion_inputs <- renderUI({
    n <- input$n_ingestion
    lapply(1:n, function(i) {
      fluidRow(
        column(6, numericInput(paste0("ingestion_", i),
                               paste("Ingestion", i, "(ng/kg/day)"), value = 10)),
        column(6, numericInput(paste0("ingestion_time_", i),
                               paste("Start time", i), value = i - 1))
      )
    })
  })

  fwd_data <- eventReactive(input$run, {
    n <- input$n_ingestion
    ingestion_vec      <- sapply(1:n, function(i) input[[paste0("ingestion_", i)]])
    ingestion_time_vec <- sapply(1:n, function(i) input[[paste0("ingestion_time_", i)]])

    n_admin        <- input$n_admin
    admin_dose_vec <- sapply(1:n_admin, function(i) input[[paste0("admin_dose_", i)]])
    admin_time_vec <- sapply(1:n_admin, function(i) input[[paste0("admin_time_", i)]])

    forward_dosimetry(
      BW             = input$BW,
      duration       = input$duration,
      time_step      = input$time_step,
      chemical       = input$chemical,
      ingestion      = ingestion_vec,
      ingestion_time = ingestion_time_vec,
      admin_dose     = admin_dose_vec,
      admin_time     = admin_time_vec,
      admin_type     = input$admin_type,
      exp_type       = input$exp_type,
      time_scale     = input$time_scale
    )
  })

  output$fwd_Cplasma_plot <- renderPlot({
    df <- fwd_data()
    ggplot(df, aes(x = time, y = CA)) +
      geom_line(colour = "#2E9E7A", linewidth = 1.1) +
      geom_area(fill = "#2E9E7A", alpha = 0.08) +
      labs(
        x     = paste0("Time (", input$time_scale, ")"),
        y     = "Plasma concentration (μg/L)",
        title = paste0(input$chemical, " — Plasma concentration vs time")
      ) +
      theme_pfas()
  })

  output$fwd_allcomp_plot <- renderPlot({
    df <- fwd_data()
    compartments <- c("CR", "CAdi", "CBra", "CGon", "CHea", "CLun",
                      "CMus", "CSki", "CSpl", "Cpan", "CKb", "Cfil",
                      "CL", "CGI", "CA")
    df_long <- df %>%
      select(time, all_of(compartments)) %>%
      pivot_longer(-time, names_to = "compartment", values_to = "concentration")

    ggplot(df_long, aes(x = time, y = abs(concentration), colour = compartment)) +
      geom_line(linewidth = 0.8, alpha = 0.85) +
      scale_y_log10(labels = scales::label_scientific()) +
      scale_colour_viridis_d(option = "turbo") +
      labs(
        x      = paste0("Time (", input$time_scale, ")"),
        y      = "Concentration (μg/L, log scale)",
        colour = NULL,
        title  = "All compartment concentrations"
      ) +
      theme_pfas()
  })

  output$fwd_elimination_plot <- renderPlot({
    df <- fwd_data()
    df_long <- df %>%
      select(time, Aurine, Afeces) %>%
      pivot_longer(-time, names_to = "route", values_to = "amount")

    ggplot(df_long, aes(x = time, y = amount, colour = route, fill = route)) +
      geom_line(linewidth = 1) +
      geom_area(alpha = 0.08, position = "identity") +
      scale_colour_manual(values = c(Aurine = "#3A86FF", Afeces = "#FF6B6B")) +
      scale_fill_manual(values  = c(Aurine = "#3A86FF", Afeces = "#FF6B6B")) +
      labs(
        x      = paste0("Time (", input$time_scale, ")"),
        y      = "Cumulative amount (μg)",
        colour = NULL, fill = NULL,
        title  = "Cumulative elimination"
      ) +
      theme_pfas()
  })

  output$fwd_ingestion_plot <- renderPlot({
    req(input$admin_type == "oral")
    df <- fwd_data()
    ggplot(df, aes(x = time, y = ingestion)) +
      geom_step(colour = "#E07B39", linewidth = 1) +
      geom_area(fill = "#E07B39", alpha = 0.08) +
      labs(
        x     = paste0("Time (", input$time_scale, ")"),
        y     = "Ingestion rate (μg/time)",
        title = "Ingestion rate over time"
      ) +
      theme_pfas()
  })

  output$dl_fwd <- downloadHandler(
    filename = function() paste0("forward_", Sys.Date(), ".csv"),
    content  = function(file) write.csv(fwd_data(), file, row.names = FALSE)
  )

  rev_results <- eventReactive(input$rev_run, {
    admin_time <- if (input$rev_admin_type == "iv") input$rev_admin_time else 0

    res <- .reverse_dosimetry(
      chemical    = input$rev_chemical,
      BW          = input$rev_BW,
      duration    = input$rev_duration,
      time_scale  = input$rev_time_scale,
      POD         = input$rev_POD,
      compartment = input$rev_compartment,
      admin_type  = input$rev_admin_type,
      admin_time  = admin_time
    )

    exposure_time <- sort(unique(c(
      seq(0, input$rev_duration, by = input$rev_time_step),
      input$rev_duration
    )))

    if (input$rev_admin_type == "oral") {
      ingestion_input <- switch(input$rev_time_scale,
        "minutes" = (1/24/60) * res$exposure * input$rev_BW / 1000,
        "hours"   = (1/24)    * res$exposure * input$rev_BW / 1000,
        "days"    =             res$exposure * input$rev_BW / 1000,
        "weeks"   = 7         * res$exposure * input$rev_BW / 1000,
        "months"  = 30        * res$exposure * input$rev_BW / 1000,
        "years"   = 365       * res$exposure * input$rev_BW / 1000
      )
      user_input <- list(
        BW             = input$rev_BW,
        exposure_time  = exposure_time,
        chemical       = input$rev_chemical,
        ingestion      = ingestion_input,
        ingestion_time = 0,
        admin_dose     = 0,
        admin_time     = 0,
        admin_type     = "oral",
        exp_type       = "continuous",
        time_scale     = input$rev_time_scale
      )
    } else {
      user_input <- list(
        BW             = input$rev_BW,
        exposure_time  = exposure_time,
        chemical       = input$rev_chemical,
        ingestion      = 0,
        ingestion_time = 0,
        admin_dose     = res$exposure,
        admin_time     = admin_time,
        admin_type     = "iv",
        exp_type       = "continuous",
        time_scale     = input$rev_time_scale
      )
    }

    list(result = res, solution = .run_extended(user_input))
  })

  output$rev_results_table <- renderTable({
    res <- rev_results()$result
    exposure_label <- if (input$rev_admin_type == "oral") "Estimated intake (ng/kg/day)"
                      else "Estimated IV dose (μg)"
    df <- data.frame(
      Chemical = input$rev_chemical,
      Tissue   = input$rev_compartment,
      POD      = res$POD,
      Duration = paste(input$rev_duration, input$rev_time_scale),
      round(res$exposure, 4),
      `Relative error` = round(res$rel_error, 6),
      check.names = FALSE
    )
    colnames(df)[5] <- exposure_label
    df
  })

  output$rev_Cplasma_plot <- renderPlot({
    df <- rev_results()$solution
    ggplot(df, aes(x = time, y = CA)) +
      geom_line(colour = "#2E9E7A", linewidth = 1.1) +
      geom_area(fill = "#2E9E7A", alpha = 0.08) +
      labs(x     = paste0("Time (", input$rev_time_scale, ")"),
           y     = "Plasma concentration (μg/L)",
           title = paste0(input$rev_chemical, " — Plasma concentration vs time")) +
      theme_pfas()
  })

  output$rev_allcomp_plot <- renderPlot({
    df <- rev_results()$solution
    compartments <- c("CR", "CAdi", "CBra", "CGon", "CHea", "CLun",
                      "CMus", "CSki", "CSpl", "Cpan", "CKb", "Cfil",
                      "CL", "CGI", "CA")
    df_long <- df %>%
      select(time, all_of(compartments)) %>%
      pivot_longer(-time, names_to = "compartment", values_to = "concentration")
    ggplot(df_long, aes(x = time, y = abs(concentration), colour = compartment)) +
      geom_line(linewidth = 0.8, alpha = 0.85) +
      scale_y_log10(labels = scales::label_scientific()) +
      scale_colour_viridis_d(option = "turbo") +
      labs(x      = paste0("Time (", input$rev_time_scale, ")"),
           y      = "Concentration (μg/L, log scale)",
           colour = NULL,
           title  = "All compartment concentrations") +
      theme_pfas()
  })

  output$rev_elimination_plot <- renderPlot({
    df <- rev_results()$solution %>%
      select(time, Aurine, Afeces) %>%
      pivot_longer(-time, names_to = "route", values_to = "amount")
    ggplot(df, aes(x = time, y = amount, colour = route, fill = route)) +
      geom_line(linewidth = 1) +
      geom_area(alpha = 0.08, position = "identity") +
      scale_colour_manual(values = c(Aurine = "#3A86FF", Afeces = "#FF6B6B")) +
      scale_fill_manual(values  = c(Aurine = "#3A86FF", Afeces = "#FF6B6B")) +
      labs(x      = paste0("Time (", input$rev_time_scale, ")"),
           y      = "Cumulative amount (μg)",
           colour = NULL, fill = NULL,
           title  = "Cumulative elimination") +
      theme_pfas()
  })

  output$rev_tissue_plot <- renderPlot({
    req(input$rev_compartment != "serum")
    col_map <- c(liver = "CL", adipose = "CAdi", brain = "CBra", gonads = "CGon",
                 gut = "CGI", heart = "CHea", lung = "CLun", muscle = "CMus",
                 skin = "CSki", kidney = "CKb")
    col  <- col_map[input$rev_compartment]
    df   <- rev_results()$solution
    ggplot(df, aes(x = time, y = .data[[col]])) +
      geom_line(colour = "#9B5DE5", linewidth = 1.1) +
      geom_area(fill = "#9B5DE5", alpha = 0.08) +
      labs(x     = paste0("Time (", input$rev_time_scale, ")"),
           y     = "Concentration (μg/L)",
           title = paste0(input$rev_compartment, " concentration vs time")) +
      theme_pfas()
  })

  output$dl_rev <- downloadHandler(
    filename = function() paste0("reverse_", Sys.Date(), ".csv"),
    content  = function(file) {
      res <- rev_results()$result
      write.csv(data.frame(
        chemical       = input$rev_chemical,
        tissue         = input$rev_compartment,
        POD            = res$POD,
        duration       = input$rev_duration,
        time_scale     = input$rev_time_scale,
        admin_type     = input$rev_admin_type,
        exposure       = res$exposure,
        relative_error = res$rel_error
      ), file, row.names = FALSE)
    }
  )
}

shinyApp(ui = ui, server = server)

