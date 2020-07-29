# This web app is used to convert mother and daughter layouts into echo instruction files
library(varhandle) # has check.numeric function
library(shinyalert) # pop-up error messages
library(platetools) # plate-view plot tools, only lightly used 
library(tidyverse)
library(shiny)

# new daughter layout function
df_to_layout <- function(df, layout_type) {
    df_m <-   set_names( df ,  c("type","row",as.numeric( df [1,-c(1,2)]))) %>%
        . [ -1 , -1] %>%
        reshape2::melt( . ,id.vars = "row") %>%
        mutate( . , well = as_vector(map2( . $row,  . $variable, paste0)) ) %>%
        set_names( . , c("row", "column", layout_type, "well"))
    df_m
}

make_layout <- function( filename ) { # from path to raw layout to a final fomatted layout file
    # read the layout file, and split each layout into an individual
    layout_list <- data.table::fread( filename, header = TRUE) %>%
        as_tibble() %>%
        split( . ,  . $Type)
    
    # put into a merge-able form
    layout <- df_to_layout(layout_list[[1]], names(layout_list)[[1]])[c(1,2,4)] # initialize the list
    for (i in c(1:length(layout_list))) {
        layout <- layout %>%
            mutate("var" =  as_vector(df_to_layout(layout_list[[i]], layout_type = names(layout_list)[[i]])[3] )) %>% # append the column of interest
            set_names(c(names(layout), names(layout_list)[[i]])) # rename based on the column of interest
    }
    layout <- layout %>%
        unite("condition", c(4:ncol(.)), remove = FALSE) %>% # create a unique column, used to define groups after averaging
        mutate_if(is.factor, as.character)
    
    layout
}

round_25 <- function( x ) {
    if ( x%%25 == 0 ) {
        return( x )
    } else {
        x + 25 - x%%25 
    }
}

convert_numerics <- function( vec ) {
    
    if(all(check.numeric(vec))){
        # convert the vector to numeric
        vec <- as.numeric(vec)
    }
    vec
}

add_DMSO_dil <- function(df_all, df_mother) {
    daughter_DMSO <- df_all %>% 
        filter(DMSO_vol > 0 )
    
    mother_DMSO <- df_mother %>% 
        filter(compound == "DMSO")
    
    n_DMSO_daughter <- daughter_DMSO %>% nrow() # number of DMSO-diluted wells in daugther
    n_DMSO_mother <- mother_DMSO %>% nrow()
    
    mother_reps <- ceiling(n_DMSO_daughter/n_DMSO_mother) 
    
    merge_vec <- rep(mother_DMSO$`Source Well`, times = mother_reps) %>%
        .[c(1:n_DMSO_daughter)]
    
    daughter_DMSO_out <- daughter_DMSO %>%
        mutate("Source Well" = merge_vec) %>%
        mutate(mother_vol = DMSO_vol)
}

make_platemap_plot <- function( data, plot_title = "", fill_var, alpha_var = NULL ) {
    
    # names(data)
    #  if(!alpha_var %in% names(data)) {alpha_var <- NULL}
    data <- data %>%
        mutate("-" = rep("", nrow(.)))
    
    fill_var <- enquo(fill_var)
    alpha_var <- enquo(alpha_var)
    
    
    p <- platetools::plate_map(data = data$well, well = data$well ) %>%
        mutate(well_f = well,
               well = as.character(well)) %>%
        select(-values) %>%
        
        left_join(data, by = "well") %>%
        mutate(well = well_f) %>%
        filter(!!fill_var != "Empty") %>%
        
        
        ggplot( . , aes_string(x = "Column", y = "Row")) +
        geom_point(data = expand.grid(seq(1, 24), seq(1, 16)),
                   aes_string(x = "Var1", y = "Var2"),
                   color = "grey90", fill = "white", shape = 22, size = 5-2, alpha = 0.1) +
        coord_fixed(ratio = (24.5 / 24) / (16.5 / 16), xlim = c(0.5, 24.5)) +
        scale_y_reverse(breaks = seq(1, 16), labels = LETTERS[1:16]) +
        scale_x_continuous(position = "top", breaks = seq(1, 24)) +
        xlab("") +
        ylab("") +
        theme_dark() +
        
       # geom_point(aes(fill = !!fill_var, alpha = !!alpha_var), colour = "gray20", shape = 22, size = 5) 
        geom_point(aes(fill = !!fill_var), colour = "gray20", shape = 22, size = 5) +
        labs(title = plot_title)
    
    if(all(check.numeric(data %>% select(!!fill_var) %>% as_vector()))){
        p <- p + 
            scale_fill_viridis_c(begin = 0.8, end = 0)
    } else {
        
        p <- p + 
            scale_fill_viridis_d(begin = 0.8, end = 0)
    }
    
    p
    
}

# server.R
server <- function(input, output) {
    values <- reactiveValues() # initialize the reactive container
    
    #### instruction template downloads
    output$download_echo_sample <- downloadHandler(
        filename = "sample_echo_instruction_file.csv",
        content = function(fname) {
            write.csv( read_csv("sample_echo_instruction_file.csv"), fname, row.names = FALSE)
        }
    )
    
    output$download_daughter_sample <- downloadHandler(
        filename = "sample_daughter_layout.csv",
        content = function(fname) {
            write.csv( read_csv("sample_daughter_layout.csv"), fname, row.names = FALSE)
        }
    )
    
    output$download_mother_sample <- downloadHandler(
        filename = "sample_mother_layout.csv",
        content = function(fname) {
            write.csv(read_csv("sample_mother_layout.csv"), fname, row.names = FALSE)
        }
    )
    
    ##### uploading the files 
    observeEvent(input$analyze_button, {
        req(input$mother_file$datapath)
        req(input$daughter_file$datapath)
        
        
        tryCatch({
            values$mother <- make_layout(input$mother_file$datapath) %>%
                            mutate_all(convert_numerics) %>%
                            select(well, compound, concentration) %>%
                            rename("Source Well" = "well",
                                   "mother_conc" = "concentration")
            
            values$daughter_raw <- make_layout(input$daughter_file$datapath) %>%
                                    mutate_all(convert_numerics) 

            
            values$all <- values$daughter_raw %>% 
                            select(well, compound, concentration, volume) %>%
                            mutate_at(c("concentration", "volume"), as.numeric) %>%
                            rename("Destination Well" = "well",
                                   "daughter_conc" = "concentration",
                                   "daughter_final_vol" = "volume") %>%
                            dplyr::filter(daughter_final_vol != 0) %>%
                            left_join( . , values$mother, by = "compound") %>%
                            mutate(mother_dil = (daughter_conc/mother_conc) * ( daughter_final_vol)) %>%
                            mutate(mother_vol = sapply(X = mother_dil, FUN = round_25)) %>%
                            mutate(DMSO_vol = daughter_final_vol - mother_vol) %>%
                            mutate("Source Plate Name" = rep("Source[1]", times = nrow(.)),
                                   "Destination Plate Name" = rep("Destination[1]", times = nrow(.)),
                                   "Destination Well X Offset"	= rep(NA, times = nrow(.)),
                                   "Destination Well Y Offset"	= rep(NA, times = nrow(.))
                            ) 
            
            values$final <- add_DMSO_dil(values$all, values$mother) %>%
                            bind_rows(values$all, . ) %>%
                            mutate(final_daughter_conc = (mother_vol * mother_conc)/daughter_final_vol)
            
            values$echo_instructions <- values$final %>%
                rename("Transfer Volume" = mother_vol) %>%
                select("Source Plate Name",	"Source Well",	"Destination Plate Name",	"Destination Well",	"Transfer Volume",	"Destination Well X Offset",	"Destination Well Y Offset")
            
            values$read_me <- tibble("Exp" = input$exp_num,
                "Date" = as.character(base::Sys.Date()),
                "Notes" = input$additional_notes,
                "Session_info" = capture.output(sessionInfo()))
            
           values$name_prefix <- paste0(input$exp_num, "--",as.character(base::Sys.Date()), "_")
            
        }, # end tryCatch function
        error = function(e){
            shinyalert("Oops, something went wrong!", "Please check the layout file formats and try again.")
            values$echo_instructions <- NULL
        }
        )
        
        output$data_table <-  DT::renderDataTable({ 
            
            req(values$echo_instructions)
            values$echo_instructions
        
            })
    })

    
    output$color_by <-renderUI({ # this is reactive by nature of being a render call? it can accept, therefore, rt(), which is a reactive expression. Can we
        req(values$daughter_raw)
        
        values$val_names <- values$daughter_raw %>% 
                            select(-c(row, column, well, condition)) %>% 
                            names()  %>% 
             c("-", .)
        
        varSelectInput("color_by", label = "Color plot by",  # as select input
                       data = values$daughter_raw %>% 
                           mutate("-" = rep("", nrow(.))) %>%
                              select( one_of(values$val_names )),
                       selected = "-"

        ) # this is able to take a reactive value
    }) 
    
    # output$alpha_by <-renderUI({ # this is reactive by nature of being a render call? it can accept, therefore, rt(), which is a reactive expression. Can we
    #     req(values$daughter_raw)
    #     
    #     values$val_names <- values$daughter_raw %>% 
    #         select(-c(row, column, well, condition)) %>% 
    #         names() 
    #     # %>% 
    #     #     c("-", .)
    #     
    #     varSelectInput("alpha_by", label = "Vary transparency by",  # as select input
    #                    data = values$daughter_raw %>% 
    #                        #mutate("-" = rep("", nrow(.))) %>%
    #                       # mutate("Don't vary transparency" = rep("", nrow(.))) %>%
    #                        select( one_of(values$val_names ) )
    #                    # ,
    #                    # selected = "-"
    #                    # selected = "Don't vary transparency"
    #     ) # this is able to take a reactive value
    # }) 
    
    # output$make_plot <-renderUI({ # this is reactive by nature of being a render call? it can accept, therefore, rt(), which is a reactive expression. Can we
    #     req(values$daughter_raw)
    #     actionButton("make_plot", "Make plot")
    # }) 
    
    # values$out_plot <- observeEvent(input$color_by,{
    #     req(values$daughter_raw)
    #     
    #     make_platemap_plot( data = values$daughter_raw, 
    #                         fill_var = !!input$color_by,
    #                         alpha_var = NULL)
    #     })
    # 
    values$plot_made <- TRUE
    # out_plot <- eventReactive( {input$color_by 
    #                             input$analyze_button}, { # only when the "update plot" button is clicked, update the plot
    out_plot <- eventReactive( input$color_by, { # only when the "update plot" button is clicked, update the plot
       print("input color by changes")
         req(values$daughter_raw)
        # plot <-  make_platemap_plot( data = values$daughter_raw, 
        #                              fill_var = !!input$color_by,
        #                              plot_title = paste0(input$exp_num, "--",as.character(base::Sys.Date()), ": Plate Map"),
        #                              alpha_var = NULL)
        tryCatch({
            print("making plot")
            values$plot_made <- TRUE
            plot <-  make_platemap_plot( data = values$daughter_raw, 
                                         fill_var = !!input$color_by,
                                         plot_title = paste0(input$exp_num, "--",as.character(base::Sys.Date()), ": Plate Map"),
                                         alpha_var = NULL)
        }, # end tryCatch function
        error = function(e){
            #shinyalert("Oops, something went wrong!", "Please check the layout file formats and try again.")
           values$plot_made <- FALSE
        })
         plot
    })
    
    output$plot <- renderPlot({
        if (values$plot_made == TRUE) {
            out_plot()   
        } else {
            NULL
        }
    })
                               
    
    #### downloads
   
    
    output$downloadData <- downloadHandler(

        
        filename = paste0(paste0(input$exp_num, "--",as.character(base::Sys.Date()), "_"), "echo_instructions.zip"),
        
        content = function(fname) {
            tmpdir <- tempdir()
            setwd(tempdir())
            print(tempdir())
            name_prefix <- paste0(input$exp_num, "--",as.character(base::Sys.Date()), "_")
            
            fs <- c(paste0(name_prefix, "plot.pdf"), 
                    paste0(name_prefix, "all_data.csv"),
                    paste0(name_prefix, "echo_instructions.csv"),
                    paste0(name_prefix, "read_me.csv")
                    )
            
            
            ggsave(paste0(name_prefix, "plot.pdf"), plot = out_plot())
            write.csv(values$final, file = paste0(name_prefix, "all_data.pdf"), sep =",")
            write.csv(values$echo_instructions, paste0(name_prefix, "echo_instructions.csv"), sep =",")
            write.csv(values$read_me, file = paste0(name_prefix, "read_me.csv"), sep =",")
            print (fs)

            zip(zipfile=fname, files=fs)
        },
        contentType = "application/zip"
    )
    
}


ui <- navbarPage(useShinyalert(),
                 tabPanel("Instructions",
                          fluidRow(
                              # left panel: interactive sliders
                              column(3),
                              column(6,
                                     tags$h1("Echo Instruction File Creation"), 
                                     tags$p("This app creates Echo instruction files. These files can be uploaded directly to the Echo software to produce a desired daughter plate from a given mother plate."),
                                     downloadButton("download_echo_sample", "Download example Echo instruction file", width = '50%',style="font-size: 14px; color: #00000; background-color: #fffff; border-color: #ffff"),
                                     tags$br(),
                                     tags$br(),

                                    tags$li( tags$strong("Step 1: Upload plate layouts."), 
                                                  "Upload layout files for both the mother and daughter plates. These layout files can have an arbitrary number of experimental variables, but require that both mother and daughter layouts have 'compound' and 'concentration' variables."), 
                                    tags$br(),     
                                    downloadButton("download_mother_sample", "Download example mother layout", width = '50%',style="font-size: 14px; color: #00000; background-color: #fffff; border-color: #ffff"),
                                    downloadButton("download_daughter_sample", "Download example daughter layout", width = '50%',style="font-size: 14px; color: #00000; background-color: #fffff; border-color: #ffff"),
                                         tags$br(),
                                        tags$br(),
                                    
                                     tags$li(tags$strong("Step 2: Create Echo instructions."),
                                             "Because this app was developed for the creation of dye daughter plates, its default process has a couple of idiosyncrasies when used more generally. I also haven't written error catchers yet. So keep an eye on the following:",
                                         tags$ul(
                                                tags$li("Both mother and daughter must have 'condition' and 'concentration' variables."),
                                                 tags$li("If the final concentration of a compound in the daughter is lower than in the mother, the app automatically calculates the diliutions and looks for 'DMSO' wells in the mother plate to bring the daughter well to its final volume. Each DMSO well in the mother will be used one-by-one for each well needing dilution in the daughter plate. All dilution volumes are rounded up to the nearest 25 nL.")
                                                )),
                                    tags$br(),
                                     tags$li(tags$strong("Step 3: Make plate plots."),
                                             "You can make a plate-view plot using the variables uploaded in the daughter layout."),
                                     tags$br(),
                                    
                                     tags$li(tags$strong("Step 4: Download."),
                                             "The 'Download all' button will download the following:",
                                             tags$ul(
                                                 tags$li("A csv of Echo instructions to create the uploaded daughter from the uploaed mother"),
                                                 tags$li("A csv with all relevant values (e.g. final concentrations, contents of each well)"),
                                                 tags$li("The currently-displayed plot")
                                            ),
                                             tags$br()
                                    
                              ),
                              tags$p("This is still pretty beta; if there are more features which are desperately needed, let me know! Thanks and good luck!")),
                              column(3)
                          )),
                 #"UCSF Dye Screen Processing",
                 tabPanel("Create Echo Instructions",
                          sidebarLayout(
                              sidebarPanel(
                                 textInput("exp_num", "Experiment number (for naming downloads)", "Exp000"),
                                 textAreaInput("additional_notes", "Additional notes","Add additional notes to self here.", height = 70),
                            
                                 tags$hr(),    
                            p("Screening files", style = "font-size: 16px;",align = "center"),

                                 fileInput("daughter_file", "Upload daughter layout (csv)",
                                           accept = c(
                                               "text/csv",
                                               "text/comma-separated-values,text/plain",
                                               ".csv")
                                 ),
                                 fileInput("mother_file", "Upload mother layout (csv)",
                                           accept = c(
                                               "text/csv",
                                               "text/comma-separated-values,text/plain",
                                               ".csv")
                                 ),

                                 
                            actionButton("analyze_button", "Create Echo Instructions"),
                            downloadButton('downloadData', 'Download results')
                          ),
                          
                          mainPanel(
                              tabsetPanel(
                                  tabPanel("Preview Echo file", 
                                           DT::dataTableOutput("data_table"),  
                                           style = "overflow-x: scroll;overflow-y: scroll;"),
                                  tabPanel("Make layout plot", 
                                           fluidRow(column(3,
                                                           uiOutput("color_by"),
                                           ),
                                           column(3, offset = 1,
                                                  uiOutput("alpha_by")
                                           ),
                                           column(3#,
                                                 # uiOutput("make_plot") ),
                                           ),
                                           plotOutput("plot")
                                           
                                           # sidebarPanel(
                                           #     uiOutput("color_by"),
                                           #     uiOutput("alpha_by")
                                           # )
                                         
                                           
                                        
                                          
                                           
                                           # tableOutput("table_set2"),  
                                           # #dataTableOutput("table_set2"), 
                                           # style = "overflow-x: scroll;")
                              )))

                          )
                        
                 
)
))

shinyApp(ui = ui, server = server)

