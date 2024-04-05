#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinydashboard)
library(shinydashboardPlus)
library(plotly)
library(DT)
library(htmlwidgets)
source("www/functions_2.R")
library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)
library(ggpubr)


# Define UI for application that draws a histogram
ui <- dashboardPage(
 # skin = "midnight",
  dashboardHeader(title = "EcB Omics"),
  dashboardSidebar(
    sidebarMenu(
      # Making the tabs
      menuItem("Proteomics Clustering", tabName = "Clustering", icon = icon("magnifying-glass-chart")),
      menuItem("vs Healthy", tabName = "Healthy", icon = icon("notes-medical")),
      menuItem("EcB Type", tabName = "EcB", icon = icon("bacteria")),
      menuItem("Mortality", tabName = "Mortality", icon = icon("skull-crossbones")),
      # Metabolomics Options # 
      menuItem("Metabolomics Clustering", tabName = "Clustering2", icon = icon("magnifying-glass-chart")),
      menuItem("vs Healthy", tabName = "Healthy2", icon = icon("notes-medical")),
      menuItem("EcB Type", tabName = "EcB2", icon = icon("bacteria")),
      menuItem("Mortality", tabName = "Mortality2", icon = icon("skull-crossbones"))
    ),
    sliderInput("range",
                "EFS Rank to Display", 1, 100, c(1,10))
    
    
  ),
  ### Proteomics ###
  dashboardBody(
    ### Comparison to Healthy###
    tabItems(
      tabItem(tabName = "Clustering",
              fluidPage(
                box(title = "Hierchical Clustering",
                    width= 12,
                    plotlyOutput("hc",width = "100%", height = "100%"))
              )),
      tabItem(tabName = "Healthy",
              h1("Healthy"),
              fluidRow(
                box(title = "Volcano plot",
                    plotlyOutput("healthy_vp")
                ),
                box(title = "ROC Curves",
                    plotOutput("healthy_roc", brush = "plot_brush")),
                box(title = "violin plots",
                    plotOutput("healthy_vip", brush = "plot_brush")),
                box(title = "statistics",
                    DTOutput('healthy_tbl'))
              )
      ),
      ### Comparison between E. faecalis and E. faecium ###
      tabItem(tabName = "EcB",
              h1("ECB"),
              fluidRow(
                box(title = "Volcano plot",
                    plotlyOutput("EcB_vp")),
                box(title = "ROC Curves",
                    plotOutput("EcB_roc", brush = "plot_brush")),
                box(title = "violin plots",
                    plotOutput("EcB_vip", brush = "plot_brush")),
                box(title = "statistics",
                    DTOutput('EcB_tbl'))
              )
                ),
      ### Comparison between Mortality and Survival ###
      tabItem(tabName = "Mortality",
              fluidRow(
                box(title = "Volcano plot",
                    plotlyOutput("mortality_vp")),
                box(title = "ROC Curves",
                    plotOutput("mortality_roc", brush = "plot_brush")),
                box(title = "violin plots",
                    plotOutput("mortality_vip", brush = "plot_brush")),
                box(title = "statistics",
                    DTOutput('Mortality_tbl'))
      )
    ),
      tabItem(tabName = "Clustering2",
              fluidPage(
                box(title = "Hierchical Clustering",
                    width= 12,
                    plotlyOutput("hc2",width = "100%", height = "100%"))
              )),
       tabItem(tabName = "Healthy2",
            h1("Healthy"),
            fluidRow(
              box(title = "Volcano plot",
                  plotlyOutput("healthy_vp2")
              ),
              box(title = "ROC Curves",
                  plotOutput("healthy_roc2", brush = "plot_brush")),
              box(title = "violin plots",
                  plotOutput("healthy_vip2", brush = "plot_brush")),
              box(title = "statistics",
                  DTOutput('healthy_tbl2'))
            )
    ),
    ### Comparison between E. faecalis and E. faecium ###
        tabItem(tabName = "EcB2",
            h1("ECB"),
            fluidRow(
              box(title = "Volcano plot",
                  plotlyOutput("EcB_vp2")),
              box(title = "ROC Curves",
                  plotOutput("EcB_roc2", brush = "plot_brush")),
              box(title = "violin plots",
                  plotOutput("EcB_vip2", brush = "plot_brush")),
              box(title = "statistics",
                  DTOutput('EcB_tbl2'))
            )
    ),
    ### Comparison between Mortality and Survival ###
        tabItem(tabName = "Mortality2",
            fluidRow(
              box(title = "Volcano plot",
                  plotlyOutput("mortality_vp2")),
              box(title = "ROC Curves",
                  plotOutput("mortality_roc2", brush = "plot_brush")),
              box(title = "violin plots",
                  plotOutput("mortality_vip2", brush = "plot_brush")),
              box(title = "statistics",
                  DTOutput('Mortality_tbl2')),
              
    
    )
    )
  )
  )
)



# Define server logic required to draw a histogram
server <- function(input, output) {
  
  ### Loading all the Data we need ### 
    ## Stats ##
    # Changing based on the type of data we are looking at. 
    healthy = read_csv("www/Infected_vs_Healthy.csv")
    EcB = read_csv("www/faecalis_v_faecium.csv")
    mortality = read_csv("www/death_vs_life.csv")
    
    
    ## Continuous data ##
    continuous = read_csv("www/protein_data.csv")
    

    ## EFS data ##
    efs_p_h = read_csv("www/healthy_prot_EFS_rank.csv")
    efs_p_e = read_csv("www/EcB_prot_EFS_rank.csv")
    efs_p_m = read_csv("www/mortality_prot_EFS_rank.csv")
  
  
    
  #### Proteomics #### 

  ### Overall Data Clustering ###
    
    # formating to get in heatmap mode
    wide = continuous %>% 
      pivot_wider(names_from = Feature, values_from = Abundance) 
    
    md = wide %>% 
      select(1:3)
    
    na_rm = wide %>% 
      select(4:length(.)) %>% 
      dplyr::select_if(~ !any(is.na(.))) 
    
    wide = bind_cols(md,na_rm) %>% 
      tibble::column_to_rownames("sample_name")
    
  
    
    output$hc = renderPlotly({
      
      heatmaply::heatmaply(wide,scale = "column",
                           fontsize_col = 5,
                           fontsize_row = 5) %>% 
        layout(height=800,width=1000)
      
      
      
      
    })

  
### Healthy Tab ### 
    ## Making the volcano plot
    output$healthy_vp = renderPlotly({
      
     p =  plot_ly(data = healthy, x = ~log2FC, y = ~(-log10(adj.pvalue)),hovertext = ~paste0(Feature,"<br>",Description)) %>% 
        add_markers(customdata = ~url) %>% 
        add_segments(x = min(healthy$log2FC[!is.infinite(healthy$log2FC)]), 
                     xend = max(healthy$log2FC[!is.infinite(healthy$log2FC)]),
                     y = -log10(0.05),
                     yend = -log10(0.05), 
                     line = list(dash = "dash"),
                     color = "red") %>% 
       layout(showlegend = FALSE,
              xaxis = list(title = "Log2FC (Infected / Healthy)"))
     
     
     onRender(
       p, "
  function(el) {
    el.on('plotly_click', function(d) {
      var url = d.points[0].customdata;
      window.open(url);
    });
  }
")
    })
    
    ## Showing the stats table
    
    healthy_stat = healthy %>% 
      select(Feature,uniprotid,Description,pvalue,adj.pvalue,log2FC) %>% 
      mutate(pvalue = round(-log10(pvalue),2),
             adj.pvalue = round(-log10(adj.pvalue),2),
             log2FC = round(log2FC,2)) %>% 
      arrange(-adj.pvalue)
    
    
    output$healthy_tbl = renderDT(healthy_stat)
    
    
    # ## Making the ROC plots 
    # 
    # efs_p_h_r = efs_p_h %>% 
    #   filter(between(rank,1,10)) %>% 
    #   pull(Feature)
    # 
    
    
    output$healthy_roc = renderPlot({
      
      filt = efs_p_h %>% 
        filter(between(rank,input$range[1],input$range[2])) %>% 
        pull(Feature)
      
      continuous_h = continuous %>% 
        filter(Feature %in% filt) %>% 
        dplyr::mutate(condition = case_when(condition %in% c("faecalis","faecium") ~ "infected",
                                            TRUE ~ "healthy"))
      
      roc = make_roc(continuous_h,filt)
      
      
      roc %>% ggplot(aes((1-specificity),sensitivity))+
         geom_point()+
          geom_line()+
         facet_wrap(~Feature)+
        ggprism::theme_prism()
    
    })
    
    
  
    
    
    ## Making the Healthy Violin Plots

    output$healthy_vip = renderPlot({

      filt = efs_p_h %>%
        filter(between(rank,input$range[1],input$range[2])) %>%
        pull(Feature)

      continuous_h = continuous %>%
        filter(Feature %in% filt) %>%
        dplyr::mutate(condition = case_when(condition %in% c("faecalis","faecium") ~ "infected",
                                            TRUE ~ "healthy"))

      vp = continuous_h %>%
        filter(Feature %in% filt)
      
      
      stat = vp %>% 
        dplyr::group_by(Feature) %>% 
        rstatix::t_test(Abundance ~ condition) %>% 
        rstatix::add_xy_position(x="condition",scales = "free") %>% 
        rstatix::add_significance()
  
      vp %>% ggboxplot(x = "condition",
                       y = "Abundance",
                       facet.by = "Feature",
                       scales = "free")+
        #facet_wrap(~condition, scales = "free") +
        stat_pvalue_manual(stat,label = "p.signif")+
        ggprism::theme_prism()
        

    })


### EcB Tab ###    
    ## Making the volcano plot
    output$EcB_vp = renderPlotly({
      
      p = plot_ly(data = EcB, x = ~log2FC, y = ~(-log10(adj.pvalue)),hovertext = ~paste0(Feature,"<br>",Description)) %>% 
        add_markers(customdata = ~url) %>% 
        add_segments(x = min(EcB$log2FC[!is.infinite(EcB$log2FC)]), 
                     xend = max(EcB$log2FC[!is.infinite(EcB$log2FC)]),
                     y = -log10(0.05),
                     yend = -log10(0.05), 
                     line = list(dash = "dash"),
                     color = "red") %>% 
        layout(showlegend = FALSE,
               xaxis = list(title = "Log2FC (Faecalis / Faecium)"))
      onRender(
        p, "
  function(el) {
    el.on('plotly_click', function(d) {
      var url = d.points[0].customdata;
      window.open(url);
    });
  }
")
    })
    
    ## Showing the stats table
    
    EcB_stat = EcB %>% 
      select(Feature,uniprotid,Description,pvalue,adj.pvalue,log2FC) %>% 
      mutate(pvalue = round(-log10(pvalue),2),
             adj.pvalue = round(-log10(adj.pvalue),2),
             log2FC = round(log2FC,2)) %>% 
      arrange(-adj.pvalue)
    
    
    output$EcB_tbl = renderDT(EcB_stat)
    
    
    ## Making the ROC curves 
    
    output$EcB_roc = renderPlot({
      
      filt = efs_p_e %>% 
        filter(between(rank,input$range[1],input$range[2])) %>% 
        pull(Feature)
      
      continuous_e = continuous %>% 
        filter(Feature %in% filt) %>% 
        dplyr::filter(condition != "healthy")
      
      roc = make_roc(continuous_e,filt)
      
      
      roc %>% ggplot(aes((1-specificity),sensitivity))+
        geom_point()+
        geom_line()+
        facet_wrap(~Feature)+
        ggprism::theme_prism()
    })
    
    output$EcB_vip = renderPlot({
      
      filt = efs_p_e %>%
        filter(between(rank,input$range[1],input$range[2])) %>%
        pull(Feature)
      
      continuous_e = continuous %>%
        filter(Feature %in% filt) %>%
        dplyr::filter(condition != "healthy")
      
      vp = continuous_e %>%
        filter(Feature %in% filt)
      
      
      stat = vp %>% 
        dplyr::group_by(Feature) %>% 
        rstatix::t_test(Abundance ~ condition) %>% 
        rstatix::add_xy_position(x="condition",scales = "free") %>% 
        rstatix::add_significance()
      
      vp %>% ggboxplot(x = "condition",
                       y = "Abundance",
                       facet.by = "Feature",
                       scales = "free")+
        #facet_wrap(~condition, scales = "free") +
        stat_pvalue_manual(stat,label = "p.signif")+
        ggprism::theme_prism()
      
      
    })
    
      

    
### Mortality Tab ###         
    ## Making the volcano plot
    output$mortality_vp = renderPlotly({
      
     p = plot_ly(data = mortality, x = ~log2FC, y = ~(-log10(adj.pvalue)),hovertext = ~paste0(Feature,"<br>",Description)) %>% 
        add_markers(customdata = ~url) %>% 
       add_segments(x = min(mortality$log2FC[!is.infinite(mortality$log2FC)]), 
                    xend = max(mortality$log2FC[!is.infinite(mortality$log2FC)]),
                    y = -log10(0.05),
                    yend = -log10(0.05), 
                    line = list(dash = "dash"),
                    color = "red") %>% 
       layout(showlegend = FALSE,
              xaxis = list(title = "Log2FC (Survival / Mortality)"))
     onRender(
       p, "
  function(el) {
    el.on('plotly_click', function(d) {
      var url = d.points[0].customdata;
      window.open(url);
    });
  }
")
    })
    
    ## Adding the stats table 
    mort_stat = mortality %>% 
      select(Feature,uniprotid,Description,pvalue,adj.pvalue,log2FC) %>% 
      mutate(pvalue = round(-log10(pvalue),2),
             adj.pvalue = round(-log10(adj.pvalue),2),
             log2FC = round(log2FC,2)) %>% 
      arrange(-adj.pvalue)
    
    
    output$Mortality_tbl = renderDT(mort_stat)
    
    
    output$mortality_roc = renderPlot({
      
      filt = efs_p_m %>% 
        filter(between(rank,input$range[1],input$range[2])) %>% 
        pull(Feature)
      
      continuous_m = continuous %>% 
        filter(Feature %in% filt) %>% 
        dplyr::filter(!is.na(death_during_admission)) %>% 
        mutate(condition = death_during_admission)
      
      roc = make_roc(continuous_m,filt)
      
      
      roc %>% ggplot(aes((1-specificity),sensitivity))+
        geom_point()+
        geom_line()+
        facet_wrap(~Feature)+
        ggprism::theme_prism()
    })
    
    output$mortality_vip = renderPlot({
      
      filt = efs_p_m %>%
        filter(between(rank,input$range[1],input$range[2])) %>%
        pull(Feature)
      
      continuous_m = continuous %>%
        filter(Feature %in% filt) %>%
        dplyr::filter(!is.na(death_during_admission)) %>% 
        mutate(condition = death_during_admission)
      
      vp = continuous_m %>%
        filter(Feature %in% filt)
      
      
      stat = vp %>% 
        dplyr::group_by(Feature) %>% 
        rstatix::t_test(Abundance ~ condition) %>% 
        rstatix::add_xy_position(x="condition",scales = "free") %>% 
        rstatix::add_significance()
      
      vp %>% ggboxplot(x = "condition",
                       y = "Abundance",
                       facet.by = "Feature",
                       scales = "free")+
        #facet_wrap(~condition, scales = "free") +
        stat_pvalue_manual(stat,label = "p.signif")+
        ggprism::theme_prism()
      
      
    })
    
    
    #### Metabolomics ####
    
    ## Loading the data 
    # stats 
    healthy_m = read_csv("www/healthy_vs_infected_metabolomics.csv")
    EcB_m = read_csv("www/faecalis_v_faecium_metabolomics.csv")
    mortality_m = read_csv("www/death_vs_life_metabolomics.csv")
    
    #Continuous
    continuous_m = read_csv("www/normalized_metabolomics.csv")
    
    
    # EFS data 
    m_efs_p_h = read_csv("www/infected_metabolomics_EFS.csv")
    m_efs_p_e = read_csv("www/EcB_metabolomics_EFS.csv")
    m_efs_p_m = read_csv("www/mortality_metabolomics_EFS.csv")
    
    ## Making the overall clustering plot
    #formating to get in heatmap mode
    
    wide2 = continuous_m %>% 
      mutate(Feature = paste0(Feature,":",compound_name_cleaned)) %>% 
      select(-compound_name_cleaned) %>% 
      pivot_wider(names_from = Feature, values_from = Abundance) 
    
    
    md2 = wide2 %>% 
      select(1:3) 
    
    na_rm2 = wide2 %>% 
      select(4:length(.)) %>% 
      dplyr::select_if(~ !any(is.na(.))) 
    
    wide2 = bind_cols(md2,na_rm2) %>% 
      tibble::column_to_rownames("sample_name")
    
    
    
    output$hc2 = renderPlotly({
      
      heatmaply::heatmaply(wide2,
                           scale = "column",
                           fontsize_col = 5,
                           fontsize_row = 5,
                           showticklabels = c(FALSE, TRUE)) %>% 
        layout(height=1200,width=1200)
      
      
    })
    
    
  
    ## Making the volcano plot
    output$healthy_vp2 = renderPlotly({
      
      p =  plot_ly(data = healthy_m, x = ~log2FC, y = ~(-log10(adj.pvalue)),hovertext = ~Feature,color = ~annotated) %>% 
        add_markers(customdata = ~url) %>% 
        add_segments(x = min(healthy_m$log2FC[!is.infinite(healthy_m$log2FC)]), 
                     xend = max(healthy_m$log2FC[!is.infinite(healthy_m$log2FC)]),
                     y = -log10(0.05),
                     yend = -log10(0.05), 
                     line = list(dash = "dash"),
                     color = "red") %>% 
        layout(xaxis = list(title = "Log2FC (Infected / Healthy)"))
      
      onRender(
        p, "
  function(el) {
    el.on('plotly_click', function(d) {
      var url = d.points[0].customdata;
      window.open(url);
    });
  }
")
    })
    
    ## Showing the stats table
    
    healthy_stat2 = healthy_m %>% 
      filter(annotated == TRUE) %>% 
      select(Feature,pvalue,adj.pvalue,log2FC) %>% 
      mutate(pvalue = round(-log10(pvalue),2),
             adj.pvalue = round(-log10(adj.pvalue),2),
             log2FC = round(log2FC,2)) %>% 
      arrange(-adj.pvalue)
    
    
    output$healthy_tbl2 = renderDT(healthy_stat2)
    
    ## Showing the ROC curves
    
    output$healthy_roc2 = renderPlot({
      
      filt = m_efs_p_h %>% 
        filter(between(rank,input$range[1],input$range[2])) %>% 
        pull(Feature)
      
      continuous2 = continuous_m %>% 
        filter(Feature %in% filt) %>% 
        dplyr::mutate(condition = case_when(condition %in% c("faecalis","faecium") ~ "infected",
                                            TRUE ~ "healthy"))
      
      roc = make_roc_m(continuous2,filt)
      
      
      roc %>% ggplot(aes((1-specificity),sensitivity))+
        geom_point()+
        geom_line()+
        facet_wrap(~compound_name,labeller = labeller(groupwrap = label_wrap_gen()))+
        ggprism::theme_prism()
      
    })
    
    
    ## Making the Healthy Violin Plots
    
    output$healthy_vip2 = renderPlot({
      
      filt = m_efs_p_h %>%
        filter(between(rank,input$range[1],input$range[2])) %>%
        pull(Feature)
      
      continuous_h = continuous_m %>%
        filter(Feature %in% filt) %>%
        dplyr::mutate(condition = case_when(condition %in% c("faecalis","faecium") ~ "infected",
                                            TRUE ~ "healthy"))
      
      vp = continuous_h %>%
        filter(Feature %in% filt)
      
      
      stat = vp %>% 
        dplyr::group_by(Feature,compound_name_cleaned) %>% 
        rstatix::t_test(Abundance ~ condition) %>% 
        rstatix::add_xy_position(x="condition",scales = "free") %>% 
        rstatix::add_significance()
      
      vp %>% ggboxplot(x = "condition",
                       y = "Abundance",
                       facet.by = "compound_name_cleaned",
                       scales = "free")+
        #facet_wrap(~condition, scales = "free") +
        stat_pvalue_manual(stat,label = "p.signif")+
        ggprism::theme_prism()
      
      
    })
    
    ### Faecalis vs Faecium
    
    ## Making the volcano plot
    output$EcB_vp2 = renderPlotly({
      
      p =  plot_ly(data = EcB_m, x = ~log2FC, y = ~(-log10(adj.pvalue)),hovertext = ~Feature,color = ~annotated) %>% 
        add_markers(customdata = ~url) %>% 
        add_segments(x = min(EcB_m$log2FC[!is.infinite(EcB_m$log2FC)]), 
                     xend = max(EcB_m$log2FC[!is.infinite(EcB_m$log2FC)]),
                     y = -log10(0.05),
                     yend = -log10(0.05), 
                     line = list(dash = "dash"),
                     color = "red") %>% 
        layout(xaxis = list(title = "Log2FC (Faecalis / Faecium)"))
      
      onRender(
        p, "
  function(el) {
    el.on('plotly_click', function(d) {
      var url = d.points[0].customdata;
      window.open(url);
    });
  }
")
    })
    
    ## Showing the stats table
    
    EcB_stat2 = EcB_m %>% 
      filter(annotated == TRUE) %>% 
      select(Feature,pvalue,adj.pvalue,log2FC) %>% 
      mutate(pvalue = round(-log10(pvalue),2),
             adj.pvalue = round(-log10(adj.pvalue),2),
             log2FC = round(log2FC,2)) %>% 
      arrange(-adj.pvalue)
    
    
    output$EcB_tbl2 = renderDT(EcB_stat2)
    
    ## Showing the ROC curves
    
    output$EcB_roc2 = renderPlot({
      
      filt = m_efs_p_e %>% 
        filter(between(rank,input$range[1],input$range[2])) %>% 
        pull(Feature)
      
      continuous2 = continuous_m %>% 
        filter(Feature %in% filt) %>% 
        dplyr::filter(condition != "healthy")
      
      roc = make_roc_m(continuous2,filt)
      
      
      roc %>% ggplot(aes((1-specificity),sensitivity))+
        geom_point()+
        geom_line()+
        facet_wrap(~compound_name,labeller = labeller(groupwrap = label_wrap_gen()))+
        ggprism::theme_prism()
      
    })
    
    
    ## Making the Healthy Violin Plots
    
    output$EcB_vip2 = renderPlot({
      
      filt = m_efs_p_e %>%
        filter(between(rank,input$range[1],input$range[2])) %>%
        pull(Feature)
      
      continuous_h = continuous_m %>%
        filter(Feature %in% filt) %>%
        dplyr::filter(condition != "healthy")
      
      vp = continuous_h %>%
        filter(Feature %in% filt)
      
      
      stat = vp %>% 
        dplyr::group_by(Feature,compound_name_cleaned) %>% 
        rstatix::t_test(Abundance ~ condition) %>% 
        rstatix::add_xy_position(x="condition",scales = "free") %>% 
        rstatix::add_significance()
      
      vp %>% ggboxplot(x = "condition",
                       y = "Abundance",
                       facet.by = "compound_name_cleaned",
                       scales = "free")+
        #facet_wrap(~condition, scales = "free") +
        stat_pvalue_manual(stat,label = "p.signif")+
        ggprism::theme_prism()
      
      
    })
    
  
    ### Mortality ### 
    
    ## Making the volcano plot
    output$mortality_vp2 = renderPlotly({
      
      p =  plot_ly(data = mortality_m, x = ~log2FC, y = ~(-log10(adj.pvalue)),hovertext = ~Feature,color = ~annotated) %>% 
        add_markers(customdata = ~url) %>% 
        add_segments(x = min(mortality_m$log2FC[!is.infinite(mortality_m$log2FC)]), 
                     xend = max(mortality_m$log2FC[!is.infinite(mortality_m$log2FC)]),
                     y = -log10(0.05),
                     yend = -log10(0.05), 
                     line = list(dash = "dash"),
                     color = "red") %>% 
        layout(xaxis = list(title = "Log2FC (Survival / Mortality)"))
      
      onRender(
        p, "
  function(el) {
    el.on('plotly_click', function(d) {
      var url = d.points[0].customdata;
      window.open(url);
    });
  }
")
    })
    
    ## Showing the stats table
    
    mortality_stat2 = mortality_m %>% 
      filter(annotated == TRUE) %>% 
      select(Feature,pvalue,adj.pvalue,log2FC) %>% 
      mutate(pvalue = round(-log10(pvalue),2),
             adj.pvalue = round(-log10(adj.pvalue),2),
             log2FC = round(log2FC,2)) %>% 
      arrange(-adj.pvalue)
    
    
    output$Mortality_tbl2 = renderDT(mortality_stat2)
    
    ## Showing the ROC curves
    
    output$mortality_roc2 = renderPlot({
      
      filt = m_efs_p_m %>% 
        filter(between(rank,input$range[1],input$range[2])) %>% 
        pull(Feature)
      
      continuous2 = continuous_m %>% 
        filter(Feature %in% filt) %>% 
        dplyr::filter(!is.na(death_during_admission)) %>% 
        dplyr::mutate(condition = death_during_admission)
      
      roc = make_roc_m(continuous2,filt)
      
      
      roc %>% ggplot(aes((1-specificity),sensitivity))+
        geom_point()+
        geom_line()+
        facet_wrap(~compound_name,labeller = labeller(groupwrap = label_wrap_gen()))+
        ggprism::theme_prism()
      
    })
    
    
    ## Making the Healthy Violin Plots
    
    output$mortality_vip2 = renderPlot({
      
      filt = m_efs_p_m %>%
        filter(between(rank,input$range[1],input$range[2])) %>%
        pull(Feature)
      
      continuous_h = continuous_m %>%
        filter(Feature %in% filt) %>%
        dplyr::filter(!is.na(death_during_admission)) %>% 
        dplyr::mutate(condition = death_during_admission)
      
      vp = continuous_h %>%
        filter(Feature %in% filt)
      
      
      stat = vp %>% 
        dplyr::group_by(Feature,compound_name_cleaned) %>% 
        rstatix::t_test(Abundance ~ condition) %>% 
        rstatix::add_xy_position(x="condition",scales = "free") %>% 
        rstatix::add_significance()
      
      vp %>% ggboxplot(x = "condition",
                       y = "Abundance",
                       facet.by = "compound_name_cleaned",
                       scales = "free")+
        #facet_wrap(~condition, scales = "free") +
        stat_pvalue_manual(stat,label = "p.signif")+
        ggprism::theme_prism()
      
      
    })
    
    
    
    
    
}

# Run the application 
shinyApp(ui = ui, server = server)
