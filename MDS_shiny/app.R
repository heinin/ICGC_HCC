#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(RColorBrewer)
#library(edgeR)
library(DESeq)
library(limma)
library(BiocManager)


#options(repos = c("CRAN" = "https://cran.rstudio.com/", "BioCsoft" = "https://bioconductor.org/packages/3.8/bioc", "BioCann" = "https://bioconductor.org/packages/3.8/data/annotation"))
#getOption("repos")

# Working directory
#setwd("/Users/hnatri/Dropbox (ASU)/ICGC_HCC/shinyMDS/MDS/")

# Defining colors
viralPalette <- brewer.pal(8, "Set1")
hbvColor <- viralPalette[1]
hcvColor <- viralPalette[2]
bothColor <- viralPalette[3]
neitherColor <- viralPalette[4]

sexTissuePalette <- brewer.pal(12, "Paired")
maleColor <- sexTissuePalette[4]
femaleColor <- sexTissuePalette[6]
tumorColor <- sexTissuePalette[2]
adjacentColor <- sexTissuePalette[8]
maleTumorColor <- sexTissuePalette[4]
maleAdjacentColor <- sexTissuePalette[3]
femaleTumorColor <- sexTissuePalette[6]
femaleAdjacentColor <- sexTissuePalette[5]

voomCounts <- readRDS("voomDGElist.Rdata")

genderColors <- ifelse(voomCounts$targets$sex=="M", maleColor,
                       ifelse(voomCounts$targets$sex=="F", femaleColor, "azure3"))

tissueColors <- ifelse(voomCounts$targets$tumor=="1", tumorColor,
                       ifelse(voomCounts$targets$tumor=="0", adjacentColor, "azure3"))

tissuegenderColors <- ifelse(voomCounts$targets$gender_tissue=="M_1", maleTumorColor,
                             ifelse(voomCounts$targets$gender_tissue=="M_0", maleAdjacentColor,
                                    ifelse(voomCounts$targets$gender_tissue=="F_1", femaleTumorColor, femaleAdjacentColor)))

libtypeColors <- ifelse(voomCounts$targets$libtype=="stranded", "pink1", "lightblue1")

infectionColors <- ifelse(voomCounts$targets$viral %in% c("HBV"), hbvColor,
                          ifelse(voomCounts$targets$viral %in% c("HCV"), hcvColor,
                                 ifelse(voomCounts$targets$viral %in% c("both"), bothColor, neitherColor)))

# Define UI for app that draws an MDS plot ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("MDS plots of the ICGC Japanese HCC data"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Input: Select dimension to plot on the X and Y axes, colors, and shapes----
      #selectInput('selection', 'Gene selection', c("common", "pairwise"), "common"),
      numericInput('ngenes', 'Number of top genes', 50, min = 1, max = 12534), 
      selectInput('xval', 'Dimension to plot on the X-axis', c(1, 2, 3, 4, 5), 1),
      selectInput('yval', 'Dimension to plot on the Y-axis', c(1, 2, 3, 4, 5), 2),
      selectInput('colorby', 'Color by', c("None", "Tissue type", "Sex", "Tissue and sex", "Infection", "Library type"), "Tissue and sex"),
      selectInput('legendpos', 'Legend position', c("topleft", "topright", "bottomleft", "bottomright"), "topleft"),
      textOutput(outputId = "test")
      #selectInput('shapeby', 'Shape by', c("None", "Tissue type", "Sex", "Tissue and sex", "Infection", "Library type"))
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Histogram ----
      plotOutput(outputId = "MDS", width = "800px", height = "800px")
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    
  xVal <- reactive({input$xval})
  yVal <- reactive({input$yval})
  colorby <- reactive({input$colorby})
  #shapeby <- input$shapeby
  selection <- reactive({input$selection})
  
  output$MDS <- renderPlot({
    
    colors <- if (input$colorby=="Tissue type") {
      tissueColors
    } else if (input$colorby=="Sex") {
      genderColors
    } else if (input$colorby=="Tissue and sex") {
      tissuegenderColors
    } else if (input$colorby=="Infection") {
      infectionColors
    } else if (input$colorby=="Library type") {
      libtypeColors
    } else { "black" }
    
    legends <- if (input$colorby=="Tissue type") {
      c("Adjacent", "Tumor")
    } else if (input$colorby=="Sex") {
      c("Male", "Female")
    } else if (input$colorby=="Tissue and sex") {
      c("Male adjacent", "Male tumor", "Female adjacent", "Female tumor")
    } else if (input$colorby=="Infection") {
      c("HBV", "HCV", "Both", "Neither")
    } else if (input$colorby=="Library type") {
      c("Stranded", "Nonstranded")
    } else { "black" }
    
    legend_colors <- if (input$colorby=="Tissue type") {
      c(adjacentColor, tumorColor)
    } else if (input$colorby=="Sex") {
      c(maleColor, femaleColor)
    } else if (input$colorby=="Tissue and sex") {
      c(maleAdjacentColor, maleTumorColor, femaleAdjacentColor, femaleTumorColor)
    } else if (input$colorby=="Infection") {
      c(hbvColor, hcvColor, bothColor, neitherColor)
    } else if (input$colorby=="Library type") {
      c("pink1", "lightblue1")
    } else { "black" }
    

    plotMDS(voomCounts, top = input$ngenes, ndim = 10, dim.plot = c(input$xval, input$yval), plot=TRUE, cex=2,
              pch = 15,
              col = colors,
              gene.selection = "common")
    
    legend(input$legendpos, pch=c(15), 
           col=legend_colors, legend=legends)
    
    })
  
  output$test <- renderText({"This Shiny App uses the plotMDS() function of 
    limma to produce MDS plots of the International Cancer Genome Consortium's 
    gene expression data on hepatitis-associated HCC. The gene selection method 
    used here is 'common'. See the limma manual for more information.
    
    Contact: heini.natri@gmail.com"})
}

# Run the application 
shinyApp(ui = ui, server = server)


