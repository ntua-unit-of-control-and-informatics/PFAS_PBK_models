library(shiny)
library(deSolve)
library(ggplot2)
library(tidyverse)
library(bslib)
library(nloptr)

source("R/Extended_Model.r")
source("R/forward_dosimetry.r")
source("R/reverse_dosimetry.r")

# ── UI ────────────────────────────────────────────────────────────────────────
ui <- navbarPage("PFAS PBK Model",
theme = bs_theme(bootswatch = "minty"),

  tabPanel("Forward Dosimetry",

    sidebarLayout(
      sidebarPanel(
        h4("Simulation Parameters"),
        numericInput("BW", "Body weight (kg)", value = 70, min = 1),
        numericInput("duration", "Duration", value = 24),
        selectInput("time_scale", "Time scale of simulation",
                    choices = c("minutes", "hours", "days", "weeks", "months", "years")),
        numericInput("time_step", "Time step", value = 0.1, min = 0.01, step = 0.01),
        selectInput("chemical", "PFAS",
                    choices = c("PFHpA", "PFOA", "PFNA", "PFDA", "PFBS",
                                "PFHxS", "PFOS", "DONA", "HFPO_DA", "PFBA", "PFHxA")),
        hr(),
        h4("Oral Exposure"),
        numericInput("n_ingestion", "Number of ingestion steps", value = 1, min = 1),
        uiOutput("ingestion_help"),
        uiOutput("ingestion_inputs"),
        selectInput("exp_type", "Exposure type",
                    choices = c("continuous", "pharmacokinetics")),
        helpText("Continuous: steady intake at the given rate (ng/kg/day).
                  Pharmacokinetics: each step is a single oral dose (ng/kg)
                  given at its start time."),
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
          condition = "input.exp_type == 'continuous'",
          plotOutput("fwd_ingestion_plot", height = "300px")
        )
      )
    )
  ),

  tabPanel("Reverse Dosimetry",
    sidebarLayout(
      sidebarPanel(
        h4("Simulation Parameters"),
        numericInput("rev_BW", "Body weight (kg)", value = 70, min = 1),
        numericInput("rev_duration", "Duration", value = 24),
        selectInput("rev_time_scale", "Time scale of simulation",
                    choices = c("minutes", "hours", "days", "weeks", "months", "years"),
                    selected = "years"),
        numericInput("rev_time_step", "Time step", value = 0.1, min = 0.01, step = 0.1),
        selectInput("rev_chemical", "PFAS",
                    choices = c("PFHpA", "PFOA", "PFNA", "PFDA", "PFBS",
                                "PFHxS", "PFOS", "DONA", "HFPO_DA", "PFBA", "PFHxA")),
        hr(),
        h4("Point of Departure"),
        numericInput("rev_POD", "POD (μg/L)", value = 1, min = 0),
        helpText("Target tissue concentration to back-calculate exposure from."),
        selectInput("rev_compartment", "Tissue",
                    choices = c("Plasma"         = "plasma",
                                "Liver"          = "liver",
                                "Adipose"        = "adipose",
                                "Brain"          = "brain",
                                "Gonads"         = "gonads",
                                "Gut (GI tract)" = "gut",
                                "Heart"          = "heart",
                                "Lung"           = "lung",
                                "Muscle"         = "muscle",
                                "Skin"           = "skin",
                                "Kidney (PTC)"   = "kidney",
                                "Pancreas"       = "pancreas")),
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
          condition = "input.rev_compartment != 'plasma'",
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
            tags$dt("Body weight (kg)"),
            tags$dd("Subject body weight used to scale physiological parameters
                     and to convert ingestion rates from ng/kg/day to absolute
                     amounts."),
            tags$dt("Duration & Time scale"),
            tags$dd("Total length of the simulation and its unit (e.g. 24 hours,
                     30 days, 1 year; a month is 30 days and a year 365 days).
                     The time step, in the same unit, controls output
                     resolution — smaller values give smoother curves but take
                     longer."),
            tags$dt("PFAS compound"),
            tags$dd("Selects the chemical-specific parameters (molecular weight,
                     free fraction, renal transporter kinetics, biliary
                     excretion, efflux and tissue partition coefficients) used
                     in the model.")
          ),

          h4("Oral Exposure"),
          p("Exposure is oral (dietary or drinking-water ingestion). Its unit
            depends on the exposure type: an intake rate in ng/kg body weight
            per day for continuous exposure, or a single dose in ng/kg body
            weight for pharmacokinetics. Times are given in the selected time
            scale."),

          h4("Exposure Types"),
          tags$dl(
            tags$dt("Continuous"),
            tags$dd("Ingestion (ng/kg/day) is applied as a constant, ongoing
                     intake rate. Each step sets a new rate from its start time
                     until the next step (or the end of the simulation).
                     Appropriate for chronic dietary exposure scenarios."),
            tags$dt("Pharmacokinetics"),
            tags$dd("Each dose (ng/kg body weight) is given once into the GI
                     lumen at its dose time and then absorbed. The dose does not
                     depend on the selected time scale. Use for single- or
                     repeated-dose pharmacokinetic studies.")
          ),

          h4("Multiple Exposure Steps"),
          p("Multiple steps model changing exposure over time. For example,
            define two continuous steps to simulate a high-exposure period
            followed by a clean period: the first with a non-zero intake
            starting at time 0, and the second with intake 0 starting at the
            switchover time. With pharmacokinetics, each step is a separate
            dose, e.g. repeated doses at fixed intervals."),

          h4("Outputs"),
          tags$ul(
            tags$li(tags$b("Plasma concentration:"), " total (free + protein-bound)
                    plasma concentration over time (μg/L)."),
            tags$li(tags$b("All compartments:"), " concentrations in all modelled
                    compartments (tissues, plasma, kidney proximal tubule
                    cells, filtrate and bone) on a log scale."),
            tags$li(tags$b("Cumulative elimination:"), " total PFAS excreted via
                    urine and faeces (μg)."),
            tags$li(tags$b("Ingestion rate:"), " (continuous only) the intake
                    rate applied to the model (μg per time-scale unit).")
          ),
          p("Results can be downloaded as a CSV file using the Download button."),

          hr(),
          h3("Reverse Dosimetry"),

          h4("Simulation Parameters"),
          tags$dl(
            tags$dt("Body weight (kg)"),
            tags$dd("Subject body weight, used to scale physiological parameters
                     and convert the estimated exposure back to ng/kg/day."),
            tags$dt("Duration & Time scale"),
            tags$dd("The simulation is run until the end of this period. The
                     estimated exposure is the constant intake that produces the
                     target tissue concentration at the end of the duration.
                     The time step controls the resolution of the output plots."),
            tags$dt("PFAS compound"),
            tags$dd("Selects the chemical-specific parameters used in the model.")
          ),

          h4("Estimated Exposure"),
          p("The optimizer searches for the constant daily oral intake
            (ng/kg/day) that results in the target tissue concentration at the
            end of the simulation duration."),

          h4("Point of Departure"),
          tags$dl(
            tags$dt("POD (μg/L)"),
            tags$dd("The target internal concentration (point of departure) in
                     the selected tissue. This is the value the model is
                     back-calculated from."),
            tags$dt("Tissue"),
            tags$dd("The compartment in which the POD is defined. Supported:
                     plasma (total plasma concentration), liver, adipose, brain,
                     gonads, gut (GI tract), heart, lung, muscle, skin,
                     kidney (proximal tubule cell concentration) and
                     pancreas.")
          ),

          h4("Outputs"),
          tags$ul(
            tags$li(tags$b("Results table:"), " estimated intake (ng/kg/day)
                    and relative error between the achieved and target
                    concentration."),
            tags$li(tags$b("Plasma concentration:"), " total plasma concentration
                    over time under the estimated exposure (μg/L)."),
            tags$li(tags$b("All compartments:"), " concentrations in all modelled
                    compartments on a log scale."),
            tags$li(tags$b("Cumulative elimination:"), " total PFAS excreted via
                    urine and faeces (μg)."),
            tags$li(tags$b("Target tissue:"), " (when tissue is not plasma)
                    concentration-time plot for the selected tissue, with the
                    POD shown as a dashed line.")
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

# ── Readable names for model outputs ──────────────────────────────────────────
comp_labels <- c(CA = "Plasma", CL = "Liver", CPTC = "Kidney (PTC)", Cfil = "Filtrate",
                 CAdi = "Adipose", CBra = "Brain", CGon = "Gonads", CHea = "Heart",
                 CLun = "Lung", CMus = "Muscle", CSki = "Skin", CSpl = "Spleen",
                 Cpan = "Pancreas", CGI = "GI tract", CR = "Bone")
route_labels <- c(Aurine = "Urine", Afeces = "Faeces")
route_colours <- c(Urine = "#3A86FF", Faeces = "#FF6B6B")
tissue_labels <- c(liver = "Liver", adipose = "Adipose", brain = "Brain",
                   gonads = "Gonads", gut = "GI tract", heart = "Heart",
                   lung = "Lung", muscle = "Muscle", skin = "Skin",
                   kidney = "Kidney (PTC)", pancreas = "Pancreas")

# ── Server ────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  # Labels and units follow the exposure type: a rate for continuous intake,
  # a single dose for pharmacokinetics
  observeEvent(input$exp_type, {
    updateNumericInput(session, "n_ingestion",
                       label = if (input$exp_type == "pharmacokinetics") "Number of doses"
                               else "Number of ingestion steps")
  })

  output$ingestion_help <- renderUI({
    if (input$exp_type == "pharmacokinetics") {
      helpText("Each dose (ng/kg body weight) is given once at its dose time,
                in the selected time scale.")
    } else {
      helpText("Each step defines a daily intake rate (ng/kg body weight/day)
                starting at a given time, in the selected time scale.")
    }
  })

  output$ingestion_inputs <- renderUI({
    n     <- input$n_ingestion
    is_pk <- input$exp_type == "pharmacokinetics"
    lapply(1:n, function(i) {
      # keep values the user already typed when the labels change
      amount <- isolate(input[[paste0("ingestion_", i)]])
      start  <- isolate(input[[paste0("ingestion_time_", i)]])
      fluidRow(
        column(6, numericInput(paste0("ingestion_", i),
                               if (is_pk) paste("Dose", i, "(ng/kg)")
                               else paste("Ingestion", i, "(ng/kg/day)"),
                               value = if (is.null(amount)) 10 else amount)),
        column(6, numericInput(paste0("ingestion_time_", i),
                               paste0(if (is_pk) "Dose time " else "Start time ",
                                      i, " (", input$time_scale, ")"),
                               value = if (is.null(start)) i - 1 else start))
      )
    })
  })

  fwd_data <- eventReactive(input$run, {
    n <- input$n_ingestion
    ingestion_vec      <- sapply(1:n, function(i) input[[paste0("ingestion_", i)]])
    ingestion_time_vec <- sapply(1:n, function(i) input[[paste0("ingestion_time_", i)]])

    forward_dosimetry(
      BW             = input$BW,
      duration       = input$duration,
      time_step      = input$time_step,
      chemical       = input$chemical,
      ingestion      = ingestion_vec,
      ingestion_time = ingestion_time_vec,
      admin_dose     = 0,
      admin_time     = 0,
      admin_type     = "oral",
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
                      "CMus", "CSki", "CSpl", "Cpan", "CPTC", "Cfil",
                      "CL", "CGI", "CA")
    df_long <- df %>%
      select(time, all_of(compartments)) %>%
      pivot_longer(-time, names_to = "compartment", values_to = "concentration") %>%
      mutate(compartment = comp_labels[compartment])

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
      pivot_longer(-time, names_to = "route", values_to = "amount") %>%
      mutate(route = route_labels[route])

    ggplot(df_long, aes(x = time, y = amount, colour = route, fill = route)) +
      geom_line(linewidth = 1) +
      geom_area(alpha = 0.08, position = "identity") +
      scale_colour_manual(values = route_colours) +
      scale_fill_manual(values  = route_colours) +
      labs(
        x      = paste0("Time (", input$time_scale, ")"),
        y      = "Cumulative amount (μg)",
        colour = NULL, fill = NULL,
        title  = "Cumulative elimination"
      ) +
      theme_pfas()
  })

  output$fwd_ingestion_plot <- renderPlot({
    df <- fwd_data()
    ggplot(df, aes(x = time, y = ingestion)) +
      geom_step(colour = "#E07B39", linewidth = 1) +
      geom_area(fill = "#E07B39", alpha = 0.08) +
      labs(
        x     = paste0("Time (", input$time_scale, ")"),
        y     = paste0("Ingestion rate (μg/", sub("s$", "", input$time_scale), ")"),
        title = "Ingestion rate over time"
      ) +
      theme_pfas()
  })

  output$dl_fwd <- downloadHandler(
    filename = function() paste0("forward_", Sys.Date(), ".csv"),
    content  = function(file) write.csv(fwd_data(), file, row.names = FALSE)
  )

  rev_results <- eventReactive(input$rev_run, {
    res <- .reverse_dosimetry(
      chemical    = input$rev_chemical,
      BW          = input$rev_BW,
      duration    = input$rev_duration,
      time_scale  = input$rev_time_scale,
      POD         = input$rev_POD,
      compartment = input$rev_compartment,
      admin_type  = "oral"
    )

    solution <- forward_dosimetry(
      BW             = input$rev_BW,
      duration       = input$rev_duration,
      time_step      = input$rev_time_step,
      chemical       = input$rev_chemical,
      ingestion      = res$exposure,
      ingestion_time = 0,
      admin_dose     = 0,
      admin_time     = 0,
      admin_type     = "oral",
      exp_type       = "continuous",
      time_scale     = input$rev_time_scale
    )

    list(result = res, solution = solution)
  })

  output$rev_results_table <- renderTable({
    res <- rev_results()$result
    data.frame(
      Chemical = input$rev_chemical,
      Tissue   = input$rev_compartment,
      `POD (μg/L)` = res$POD,
      Duration = paste(input$rev_duration, input$rev_time_scale),
      `Estimated intake (ng/kg/day)` = round(res$exposure, 4),
      `Relative error` = round(res$rel_error, 6),
      check.names = FALSE
    )
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
                      "CMus", "CSki", "CSpl", "Cpan", "CPTC", "Cfil",
                      "CL", "CGI", "CA")
    df_long <- df %>%
      select(time, all_of(compartments)) %>%
      pivot_longer(-time, names_to = "compartment", values_to = "concentration") %>%
      mutate(compartment = comp_labels[compartment])
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
      pivot_longer(-time, names_to = "route", values_to = "amount") %>%
      mutate(route = route_labels[route])
    ggplot(df, aes(x = time, y = amount, colour = route, fill = route)) +
      geom_line(linewidth = 1) +
      geom_area(alpha = 0.08, position = "identity") +
      scale_colour_manual(values = route_colours) +
      scale_fill_manual(values  = route_colours) +
      labs(x      = paste0("Time (", input$rev_time_scale, ")"),
           y      = "Cumulative amount (μg)",
           colour = NULL, fill = NULL,
           title  = "Cumulative elimination") +
      theme_pfas()
  })

  output$rev_tissue_plot <- renderPlot({
    req(input$rev_compartment != "plasma")
    col_map <- c(liver = "CL", adipose = "CAdi", brain = "CBra", gonads = "CGon",
                 gut = "CGI", heart = "CHea", lung = "CLun", muscle = "CMus",
                 skin = "CSki", kidney = "CPTC", pancreas = "Cpan")
    col  <- col_map[input$rev_compartment]
    df   <- rev_results()$solution
    ggplot(df, aes(x = time, y = .data[[col]])) +
      geom_line(colour = "#9B5DE5", linewidth = 1.1) +
      geom_area(fill = "#9B5DE5", alpha = 0.08) +
      geom_hline(yintercept = input$rev_POD, linetype = "dashed", colour = "grey40") +
      labs(x     = paste0("Time (", input$rev_time_scale, ")"),
           y     = "Concentration (μg/L)",
           title = paste0(input$rev_chemical, " — ",
                          tissue_labels[input$rev_compartment],
                          " concentration vs time (dashed line: POD)")) +
      theme_pfas()
  })

  output$dl_rev <- downloadHandler(
    filename = function() paste0("reverse_", Sys.Date(), ".csv"),
    content  = function(file) {
      res <- rev_results()$result
      write.csv(data.frame(
        chemical       = input$rev_chemical,
        tissue         = input$rev_compartment,
        POD_ug_L       = res$POD,
        duration       = input$rev_duration,
        time_scale     = input$rev_time_scale,
        exposure_ng_kg_day = res$exposure,
        relative_error = res$rel_error
      ), file, row.names = FALSE)
    }
  )
}

shinyApp(ui = ui, server = server)

