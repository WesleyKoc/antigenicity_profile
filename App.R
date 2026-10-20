library(ggplot2)
library(data.table)
library(shiny)
library(plotly)
library(cli)
library(shinycssloaders)
library(DT)
library(r3dmol)
library(bio3d)
library(Biostrings)
library(pwalign)
library(bslib)
library(msaR)
data(BLOSUM100)
source("R/config.R")
source("R/utils.R")
options(shiny.maxRequestSize = 1024 * 1024 * 2048) # 2 GB
#unlink("result", recursive = TRUE)

brrep <- function(pixels = 20) {
  div(style = sprintf("height: %dpx;", pixels))
}

# Define UI ----
ui <- fluidPage(
  theme = bs_theme(version = 5, bootswatch = "morph"),
  
  titlePanel(
    div(
      p("FluWatch ☢️ \n", tags$em('@TropMed')),
      tags$style(HTML(
        "body{margin-top: 50px;}"
      )),
      style = "min-height: 120px; display: flex; justify-content: center;
      align-items: center; width: 100%; font-family: 'Times New Roman'; font-size: 3.5vw"
        )
    ),
  
  # Main panel content first
  tabsetPanel(type = "tabs",
              tabPanel("Demo",
                       tags$style(HTML(
                         "body {
                           margin-right: 150px;
                           margin-left: 150px;
                           margin-bottom: 150px;
                         }"
                       )),
                       br(),
                       fluidRow(tags$p("Welcome to FluWatch! A tool for testing if your influenza virus strain works with current vaccine strains in circulation!"),
                         tags$a(href="https://www.youtube.com/watch?v=9XDaVNsC_4g", "Click here for a youtube demo of the app!")),
                       fluidRow(
                         column(
                           width = 12,
                           br(),
                           div(id="demo_report_section",
                             div(style = "font-size:17px; text-align:center; border:2px solid #A9A9A9; 
                               background-color:#F5F5F5; padding:15px; border-radius:10px; width:100%; margin:0 auto;", 
                                 textOutput("demo_report")),
                             )
                           )
                         ),
                       br(),
                       fluidRow(
                         column(
                           width = 12,
                           div(id="antigenicity_demo_section",
                               div(
                                 style = "width: 100%; margin: 0 auto;",
                                 plotlyOutput("plot_antigencity", height = "500px")
                                 )
                               )
                           )
                       ),
                       hr(),
                       fluidRow(
                         column(
                           width = 10,
                           offset=1,
                           div(id="demo_protein_model_section",
                               #div(style = "font-size: 17px;",style = "text-align: center;",textOutput("num_epitope")),
                               selectInput(inputId = "protein_file",
                                           label = "Choose Protein Model to Visualise Hemagglutinin Mutations On:",
                                           list("2009 H1N1 influenza virus" = "3LZG",
                                                "H3N2 A/Hong Kong/1/1968 influenza virus hemagglutinin" = "6CEX",
                                                "H5N1 influenza virus" = "2FK0",
                                                "H5N1 A/goose/Guangdong/1/96" = "4MHI",
                                                "H5N1 hemagglutinin mutant (N158D, N224K and Q226L) A/Viet Nam/1203/2004"="4N5Y",
                                                "H10N8 A/Jiangxi-Donghu/346/2013"="4WSX",
                                                "H17N10 A/little yellow-shouldered bat/Guatemala/060/2010"="4H32",
                                                "H1N1 A/South Carolina/1/1918"="4JUG",
                                                "H5N1 A/duck/Laos/3295/2006(H5N1)"="4JUL",
                                                "H7N9 A/Hong Kong/61/2016"="6D7C"
                                                ),
                                           width = 999),
                               r3dmolOutput("mol", height = "500px"),
                               br(),
                               plotlyOutput("plot_epitope", height = "300px", width = "100%")
                               )
                           )
                         ),
                      hr(),
                      br(),
                      fluidRow(
                        column(
                          width = 12,
                          dataTableOutput("segment_summary")
                          )
                        )),
              tabPanel("User's analysis",
                       tags$style(HTML(
                         "body {
                           margin-right: 150px;
                           margin-left: 150px;
                         }
                         .radio label {
                         white-space: nowrap;
                         }"
                       )),
                       br(),
                       fluidRow(
                         column(
                           width = 5,
                           br(),
                           radioButtons( 
                             inputId = "radio", 
                             label = "Analysis Type", 
                             choices = list( 
                               "Genome Assembly with Antigenicity Prediction" = 1, 
                               "Antigenicity Prediction Only" = 2
                             ),
                           )
                         ),
                         column(
                           width = 3,
                           conditionalPanel(
                             condition = "input.radio == 1",
                             fileInput(
                               "fastq_files", 
                               "Upload FASTQ File", 
                               multiple = FALSE,
                               accept = c(".fastq",".fastq.gz", ".fq",".fq.gz","application/gzip",
                                          "application/x-gzip")
                             )
                           ),
                           conditionalPanel(
                             condition = "input.radio == 2",
                             fileInput(
                               "fasta_files", 
                               "Upload FASTA File", 
                               multiple = FALSE,
                               accept = c(".fasta",".fasta.gz", ".fa",".fa.gz")
                             )
                           )
                         ),
                         column(
                           width = 5,
                           radioButtons( 
                             inputId = "subtype_choice", 
                             label = "Vaccine Type", 
                             choices = list( 
                               "H1N1" = "H1N1", 
                               "VIC" = "VIC",
                               "H3N2" = "H3N2"
                             ),
                           )
                         ),
                         column(
                           width = 5,
                           uiOutput("vax_strain_ui")
                         ),
                         column(
                           width = 2,
                           dateInput("date1", "Collected Date:", value = "2020-02-10")
                         ),
                         column(
                           width = 1,
                           br(), # adds spacing before button
                           actionButton("run_pipeline", "Start", class = "btn-primary")
                         )
                       ),
                       br(),
                       fluidRow(
                         column(
                           width = 12,
                           div(style = "font-size:17px; text-align:center; border:2px solid #A9A9A9; 
                               background-color:#F5F5F5; padding:15px; border-radius:10px; width:100%; 
                               margin:0 auto;", uiOutput("demo_report_user")),
                           br(),
                           div(
                             style = "width: 100%; margin: 0 auto;",
                             plotlyOutput("plot_antigencity_user", height = "500px"),
                           ),
                         )
                       ),
                       br(),
                       fluidRow(
                         column(
                           width=12,
                           h4("Amino Acid Differences/Mismatches"),
                           textOutput("mismatches")
                         )
                       ),
                       br(),
                       fluidRow(
                         column(
                           width=12,
                           msaROutput("alignment", height="500px", width = "100%")
                         )
                       ),
                       fluidRow(
                         column(
                           width = 10,
                           offset = 1,
                           uiOutput("protein_analysis_user"),
                           r3dmolOutput("mol_user", height = "500px"),
                           hr(),
                           plotlyOutput("plot_epitope_user", height = "500px", width = "100%"),
                         )
                       ),
                       br(),
                       fluidRow(
                         column(
                           br(),
                           width = 12,
                           dataTableOutput("segment_summary_usr")
                         )
                       )
                       ),
              tabPanel("Help (Eng)",
                       tags$style(HTML(
                         "body {
                           margin-right: 150px;
                           margin-left: 150px;
                         }
                         .radio label {
                         white-space: nowrap;
                         }"
                       )),
                       br(),
                       fluidRow(
                         tags$h1("Welcome to FluWatch, the influenza vaccine antigenicity detector!"),
                         hr(style = "border-top: 2px solid #337ab7;"),
                         tags$h2("User-input data guide"),
                         brrep(4),
                         tags$p("Start off by selecting if you wanna do both genome assembly and antigenicity prediction, or just antigenicity prediction."),
                         tags$p("If you're doing genome assembly and antigenicity prediction, you need to provide a FastQ file of your desired influenza strain's haemagglutinin protein DNA sequencing data."),
                         brrep(1),
                         tags$img(src = "app_images/genomeassem+antigen.png", style = "width: 65%; height: 65%"), 
                         brrep(10),
                         tags$p("If you're only doing antigenicity prediction, you need to provide a Fasta file of your desired influenza strain's haemagglutinin protein's DNA sequence"),
                         brrep(1),
                         tags$img(src = "app_images/antigenonly.png", style = "width: 65%; height: 65%"),
                         brrep(2),
                         tags$p("Either way, the upload should only take a few seconds"),
                         brrep(4),
                         hr(style = "border-top: 2px solid #337ab7;"),
                         brrep(4),
                         tags$p("Once the upload is complete, you need to select a vaccine strain to compare against"),
                         brrep(2),
                         tags$img(src="app_images/analysisinput.png", width="100%", height="100%"),
                         brrep(2),
                         tags$p("The selection starts with whether the inactivated vaccine targets H1N1, VIC, or H3N2, select as appropriate"),
                         tags$p("Then, the dropdown menu will show you the vaccine strains in the given target you selected available for comparing against"),
                         tags$p("The number of options per vaccine target will be updated over time as more data gets added to the database"),
                         tags$p("The default selection is H1N1 strain A/California/7/2009"),
                         brrep(4),
                         hr(style = "border-top: 2px solid #337ab7;"),
                         brrep(4),
                         tags$p("If you want, there is an option to change the collection date of your data"),
                         tags$img(src="app_images/dateselection.png", width="100%", height="100%"),
                         tags$p("It simply eases the comparison of the antigenicity against different vaccine strains, but it has little effect on the final result"),
                         brrep(4),
                         hr(style = "border-top: 2px solid #337ab7;"),
                         brrep(4),
                         tags$h2("Results interface guide"),
                         brrep(4),
                         tags$p("The results page is divided into 4 sections:"),
                         tags$p("1.  Antigenicity Plot"),
                         tags$p("2.  Mutant Epitope Molecular Visualisation Panel"),
                         tags$p("3.  Mutant Epitope Plot"),
                         tags$p("4.  User Input Statistic Summary"),
                         brrep(4),
                         hr(style = "border-top: 2px solid #337ab7;"),
                         brrep(4),
                         tags$h4("Antigenicity Plot"),
                         tags$p("This plot shows the predicted efficacy of your selected vaccine against your selected virus strain given by the HA gene sequence."),
                         tags$img(src="app_images/antigenicityplotuser.png", width="100%", height="100%"),
                         tags$p("Above the plot, the summary of how well your selected vaccine strain performs against your selected virus strain is indicated in percentage change."),
                         tags$p("In the plot, the red dotted line indicates the theoretically ideal level of performance your vaccine should perform at while each virus strain stored in the database is distributed across the X-axis by the collection date and season, and each virus' performance is indicated by the Y-axis."),
                         tags$p("The level of your red triangle indicates how well your selected vaccine was predicted to perform it against other viral strains across the years, so don't feel too disappointed if it's underperforming a bit in the model, most vaccines do anyway."),
                         tags$p("Hover over any of the points to see which specific influenza strain is being represented, as well as its collection date and the vaccine's predicted efficacy against said strain."),
                         brrep(4),
                         hr(style = "border-top: 2px solid #337ab7;"),
                         brrep(4),
                         tags$h4("Mutant Epitope Molecular Visualisation Panel"),
                         tags$img(src="app_images/usermoleculemodel.png", width="100%", height="100%"),
                         tags$p("This panel generates an interactive 3D model of existing hemagglutinin structures, where the locations of all mutations on your given virus strain (relative to the consensus sequence in the app's database) will also be mapped out. Only the mutations on any of the 5 canonical epitopes are displayed."),
                         tags$p("There are currently 10 models to select from, which you can select from the dropdown menu above the panel. Each option is indicated by the specific influenza strain followed by its official name."),
                         brrep(4),
                         hr(style = "border-top: 2px solid #337ab7;"),
                         brrep(4),
                         tags$h4("Mutant Epitope Plot"),
                         tags$img(src="app_images/userepitopeplot.png", width="100%", height="100%"),
                         tags$p("This panel shows the distribution of the total mutations found in your virus strain across the 5 canonical epitopes A through E. It's just a basic bar chart showing how many mutations are at each epitope. The more mutations you have, the greater the reduction in vaccine efficacy due to divergence from conditions supporting antigenicity"),
                         )
                       ),
              tabPanel("Help (Thai)",
                       tags$style(HTML(
                         "body {
                           margin-right: 150px;
                           margin-left: 150px;
                         }
                         .radio label {
                         white-space: nowrap;
                         }"
                       )),
                       br(),
                       fluidRow(
                         tags$h1("ยินดีต้อนรับสู่ FluWatch เครื่องมือตรวจสอบคุณสมบัติทางแอนติเจนของวัคซีนไข้หวัดใหญ่!"),
                         hr(style = "border-top: 2px solid #337ab7;"),
                         tags$h2("คู่มือสำหรับข้อมูลที่ผู้ใช้ป้อนเข้า"),
                         brrep(4),
                         tags$p("เริ่มต้นด้วยการเลือกก่อนว่าคุณต้องการทำการประกอบจีโนมและการทำนายแอนติเจนิกส์พร้อมกัน หรือต้องการทำนายแอนติเจนิกส์เพียงอย่างเดียว"),
                         tags$p("หากคุณกำลังดำเนินการประกอบจีโนมและการทำนายแอนติเจน คุณต้องจัดเตรียมไฟล์ FastQ ที่มีข้อมูลลำดับดีเอ็นเอของโปรตีนฮีแมกกลูตินินของสายพันธุ์ไวรัสไข้หวัดใหญ่ที่คุณต้องการ"),
                         brrep(1),
                         tags$img(src = "app_images/genomeassem+antigen.png", style = "width: 65%; height: 65%"), 
                         brrep(10),
                         tags$p("หากคุณต้องการทำนายแอนติเจนิกส์เพียงอย่างเดียว คุณต้องจัดเตรียมไฟล์ Fasta ที่มีลำดับดีเอ็นเอของโปรตีนฮีแมกกลูตินินของสายพันธุ์ไวรัสไข้หวัดใหญ่ที่คุณต้องการ"),
                         brrep(1),
                         tags$img(src = "app_images/antigenonly.png", style = "width: 65%; height: 65%"),
                         brrep(2),
                         tags$p("ไม่ว่าจะด้วยวิธีใด การอัปโหลดก็น่าจะใช้เวลาเพียงไม่กี่วินาทีเท่านั้น"),
                         brrep(4),
                         hr(style = "border-top: 2px solid #337ab7;"),
                         brrep(4),
                         tags$p("เมื่อการอัปโหลดเสร็จสมบูรณ์ คุณจำเป็นต้องเลือกสายพันธุ์วัคซีนเพื่อนำมาเปรียบเทียบ"),
                         brrep(2),
                         tags$img(src="app_images/analysisinput.png", width="100%", height="100%"),
                         brrep(2),
                         tags$p("ขั้นตอนการเลือกเริ่มต้นจากการระบุว่าวัคซีนชนิดเชื้อตายนั้นมุ่งเป้าไปที่สายพันธุ์ H1N1, VIC หรือ H3N2 โดยให้เลือกตามความเหมาะสม จากนั้นเมนูแบบดร็อปดาวน์จะแสดงสายพันธุ์วัคซีนที่มีอยู่ภายใต้กลุ่มเป้าหมายที่คุณเลือก เพื่อให้คุณนำไปเปรียบเทียบกันได้ ทั้งนี้ จำนวนตัวเลือกในแต่ละกลุ่มเป้าหมายจะมีการอัปเดตอยู่เรื่อยๆ ตามข้อมูลที่ถูกเพิ่มเข้ามาในฐานข้อมูล"),
                         tags$p("ตัวเลือกเริ่มต้นคือเชื้อ H1N1 สายพันธุ์ A/California/7/2009"),
                         brrep(4),
                         hr(style = "border-top: 2px solid #337ab7;"),
                         brrep(4),
                         tags$p("หากต้องการ คุณสามารถเลือกเปลี่ยนวันที่จัดเก็บข้อมูลของคุณได้"),
                         tags$img(src="app_images/dateselection.png", width="100%", height="100%"),
                         tags$p("สิ่งนี้เพียงแค่ช่วยให้การเปรียบเทียบคุณสมบัติทางแอนติเจนของวัคซีนที่พัฒนาจากสายพันธุ์ต่าง ๆ ทำได้ง่ายขึ้นเท่านั้น แต่แทบไม่มีผลกระทบใด ๆ ต่อผลลัพธ์สุดท้ายเลย"),
                         brrep(4),
                         hr(style = "border-top: 2px solid #337ab7;"),
                         brrep(4),
                         tags$h2("คู่มืออินเทอร์เฟซผลลัพธ์"),
                         brrep(4),
                         tags$p("หน้าแสดงผลลัพธ์แบ่งออกเป็น 4 ส่วน ได้แก่"),
                         tags$p("1. แผนภูมิแอนติเจนิก"),
                         tags$p("2. แผงภาพจำลองโมเลกุลของอีพิโทปกลายพันธุ์"),
                         tags$p("3. แผนภาพอีพิโทปกลายพันธุ์"),
                         tags$p("4. สรุปสถิติข้อมูลที่ผู้ใช้ป้อนเข้ามา"),
                         brrep(4),
                         hr(style = "border-top: 2px solid #337ab7;"),
                         brrep(4),
                         tags$h4("Antigenicity Plot"),
                         tags$p("แผนภาพนี้แสดงประสิทธิภาพที่คาดการณ์ไว้ของวัคซีนที่คุณเลือก ในการรับมือกับสายพันธุ์ไวรัสเฉพาะที่คุณระบุ โดยอ้างอิงจากลำดับยีนเฮแมกกลูตินิน"),
                         tags$img(src="app_images/antigenicityplotuser.png", width="100%", height="100%"),
                         tags$p("เหนือแผนภาพดังกล่าว จะมีการแสดงสรุปประสิทธิภาพของสายพันธุ์วัคซีนที่คุณเลือก เมื่อเทียบกับสายพันธุ์ไวรัสที่คุณเลือก โดยระบุเป็นค่าการเปลี่ยนแปลงในรูปแบบเปอร์เซ็นต์"),
                         tags$p("ในแผนภาพนี้ เส้นประสีแดงแสดงถึงระดับประสิทธิภาพในอุดมคติตามทฤษฎีที่วัคซีนควรทำได้ ส่วนไวรัสแต่ละสายพันธุ์ที่จัดเก็บอยู่ในฐานข้อมูลจะถูกกระจายอยู่บนแกน X ตามวันที่และฤดูกาลที่เก็บตัวอย่าง โดยประสิทธิภาพของไวรัสแต่ละสายพันธุ์จะแสดงไว้บนแกน Y"),
                         tags$p("ระดับของรูปสามเหลี่ยมสีแดงบ่งบอกว่าวัคซีนที่คุณเลือกไว้มีประสิทธิภาพในการรับมือกับไวรัสสายพันธุ์อื่น ๆ ตลอดช่วงหลายปีที่ผ่านมาได้ดีเพียงใด ดังนั้นหากผลลัพธ์ในแบบจำลองออกมาไม่ดีนัก ก็อย่าเพิ่งผิดหวังไปเลย เพราะโดยปกติแล้ววัคซีนส่วนใหญ่ก็มักจะมีประสิทธิภาพในระดับนี้อยู่แล้ว"),
                         tags$p("เลื่อนเมาส์ไปวางเหนือจุดใดจุดหนึ่งเพื่อดูว่าจุดนั้นแสดงถึงเชื้อไข้หวัดใหญ่สายพันธุ์ใด รวมถึงวันที่เก็บตัวอย่างเชื้อ และประสิทธิภาพของวัคซีนที่คาดการณ์ไว้ในการป้องกันเชื้อสายพันธุ์ดังกล่าว"),
                         brrep(4),
                         hr(style = "border-top: 2px solid #337ab7;"),
                         brrep(4),
                         tags$h4("แผงภาพจำลองโมเลกุลของอีพิโทปกลายพันธุ์"),
                         tags$img(src="app_images/usermoleculemodel.png", width="100%", height="100%"),
                         tags$p("ส่วนนี้จะสร้างแบบจำลอง 3 มิติเชิงโต้ตอบของโครงสร้างโปรตีนเฮแมกกลูตินิน (hemagglutinin) ที่มีอยู่ โดยจะแสดงตำแหน่งของการกลายพันธุ์ทั้งหมดที่พบในสายพันธุ์ไวรัสที่คุณระบุ (เมื่อเทียบกับลำดับมาตรฐานในฐานข้อมูลของแอปพลิเคชัน) ทั้งนี้ จะแสดงเฉพาะการกลายพันธุ์ที่เกิดขึ้นภายในบริเวณ epitope หลักทั้ง 5 ตำแหน่งเท่านั้น"),
                         tags$p("ปัจจุบันมีโมเดลให้เลือกทั้งหมด 10 แบบ โดยคุณสามารถเลือกได้จากเมนูแบบเลื่อนลงที่อยู่เหนือแผงควบคุม แต่ละตัวเลือกจะระบุสายพันธุ์ของเชื้อไข้หวัดใหญ่ไว้ ตามด้วยชื่อเรียกอย่างเป็นทางการของสายพันธุ์นั้นๆ"),
                         brrep(4),
                         hr(style = "border-top: 2px solid #337ab7;"),
                         brrep(4),
                         tags$h4("แผนภาพอีพิโทปกลายพันธุ์"),
                         tags$img(src="app_images/userepitopeplot.png", width="100%", height="100%"),
                         tags$p("แผนภาพนี้แสดงการกระจายตัวของการกลายพันธุ์ทั้งหมดที่พบในสายพันธุ์ไวรัสของคุณ โดยแบ่งตามตำแหน่งเอพิโทป (epitope) หลัก 5 ตำแหน่ง คือ A ถึง E ซึ่งนำเสนอในรูปแบบแผนภูมิแท่งพื้นฐานที่แสดงจำนวนการกลายพันธุ์ในแต่ละเอพิโทป ยิ่งมีการกลายพันธุ์มากเท่าใด ประสิทธิภาพของวัคซีนก็จะยิ่งลดลงมากขึ้นเท่านั้น เนื่องจากไวรัสมีการเปลี่ยนแปลงไปจากสภาวะที่เอื้อต่อการจดจำของภูมิคุ้มกัน (antigenicity)"),
                       )
              )
              )
)

# Define server logic ----
server <- function(input, output) {
  
  ## execute code
  result_ready <- reactiveVal(FALSE)
  observeEvent(input$run_pipeline, {
    
    #Assembly+antigenicity
    if (input$radio == "1") {
    req(!is.null(input$fastq_files))
    input_path <- input$fastq_files$datapath
    req(length(input_path) == 1, !is.na(input_path), nzchar(input_path))
    req(file.exists(input_path))
    message("Shiny is launching assembly + antigenicity pipeline with: ", input_path)
    status <- system2("conda",args = c("run", "--no-capture-output", "-n", "antigenicity_profile","bash", "Flu_assembler.sh", shQuote(input_path)))
    if (status != 0) {
      showNotification("Assembly pipeline failed; inspect the R console.", type = "error")
      return()
    }
    }
    
    #Antigenicity only
    else if (input$radio == "2") {
      req(!is.null(input$fasta_files))
      input_path <- input$fasta_files$datapath
      req(length(input_path) == 1, !is.na(input_path), nzchar(input_path))
      req(file.exists(input_path))
      message("Shiny is launching antigenicity pipeline with: ", input_path)
      dir.create("result/consensus", recursive = TRUE, showWarnings = FALSE)
      ok <- file.copy (
        from = input_path,
        to = "result/consensus/draft_segment_4.fasta",
        overwrite =TRUE
      )
      
      if (!ok) {
        showNotification("Sorry, FASTA file copy failed, please check ur server logs", type = "error")
        return()
      }
      
      stats_dir <- file.path(getwd(), "result", "stat")
      dir.create(stats_dir, recursive = TRUE, showWarnings = FALSE)
      stats_file <- file.path(stats_dir, "summary_statistics.tsv")
      
      if (!file.exists(stats_file)) {
        dummy_file <- data.frame(
          `genomic segment` = paste0("Segment_", 1:8),
          `Consensus Length` = NA_integer_,
          `Read count` = NA_integer_,
          `Averaged genomic depth` = NA_integer_,
          `N number` = NA_integer_,
          check.names = FALSE
        )
        fwrite(dummy_file, stats_file, sep = "\t")
      }
      
      message("Copied to: ", normalizePath("result/consensus/draft_segment_4.fasta"))
      message("File exists? ", file.exists("result/consensus/draft_segment_4.fasta"))
      #status <- system2("conda",args = c("run", "--no-capture-output", "-n", "antigenicity_profile","bash", "Flu_assembler.sh", shQuote(input_path)))
      #if (status != 0) {
      # showNotification("Assembly pipeline failed; inspect the R console.", type = "error")
      #  return()
      #}
    }
    
    result_ready(TRUE)
  }, ignoreInit = TRUE)

  

  ## check subtypes ##
  subtype <- reactive({
    
    file <- "result/consensus/draft_segment_4.fasta"
    req(file.exists(file))
    
    HA_seg <- readDNAStringSet("result/consensus/draft_segment_4.fasta") |> translate_AA()
    nr <- readAAStringSet("data/01_numbering_reference_strain/amino_acid.fasta")
    seqs <- c(HA_seg,nr)
    
    score_vec <- c()
    for (u in c(2:4)) {
      alm <- pwalign::pairwiseAlignment(seqs[[1]], seqs[[u]], substitutionMatrix=BLOSUM100) |> score()  
      score_vec <- c(score_vec,alm)
    }
    
    max_score <- which.max(score_vec)
    if (max_score == 1) {
        type <- "H1N1"      
    } else if (max_score == 2){
        type <- "VIC"
    } else{
        type <- "H3N2"
    }
    
    type
    
  })
  
  ## select virus strain ##
  output$vax_strain_ui <- renderUI({
    req(input$subtype_choice)
    path <- switch(
      input$subtype_choice,
      "H1N1" = "data/02_vaccine_strain/H1N1.fasta",
      "VIC" = "data/02_vaccine_strain/VIC.fasta",
      "H3N2" = "data/02_vaccine_strain/H3N2.fasta"
    )
    
    headers <- get_fasta_headers(path)
    
    selectInput(
      inputId = "vaccine_strain",
      label = "Vaccine Strain to Compare",
      choices = headers,
      width = "100%"
    )
  })
  
  ## prediction ##
  
  pred <- reactive({
    
    file <- "result/consensus/draft_segment_4.fasta"
    req(file.exists(file))
    dir.create("result/tmp", recursive = TRUE, showWarnings = FALSE)
    
    type <- input$subtype_choice
    
    sample_AA <- translate_AA(readDNAStringSet("result/consensus/draft_segment_4.fasta"))
  
    # Target
    vax_id <- input$vaccine_strain
    vacx <- readAAStringSet(paste0("data/02_vaccine_strain/",type,".fasta"))
    vax_id_clean <- sub("^>", "", vax_id)
    vaccine_strain <- vacx[vax_id_clean]
    names(vaccine_strain) <- "reference"
    
    # Circulating strain
    cir_strain <- readAAStringSet(paste0("data/04_circulating_strain/",type,".fasta"))
    cir_strain <- cir_strain[sample(length(cir_strain), 100)]
    
    temp_in <- "result/tmp/temp_seqs.fasta"
    temp_out <- "result/tmp/aligned.fasta"
    
    writeXStringSet(c(sample_AA,vaccine_strain,cir_strain), temp_in)
    
    muscle_status <- system2(
      "conda",
      args = c(
        "run", "--no-capture-output", "-n", "antigenicity_profile",
        "muscle", "-align", temp_in, "-output", temp_out
      ),
      stdout = TRUE,
      stderr = TRUE
    )
    message("muscle status: ", attr(muscle_status, "status"))
    message("muscle output:\n", paste(muscle_status, collapse = "\n"))
    
    if (!file.exists(temp_out)) {
      stop("MUSCLE did not produce aligned.fasta, check logs above.")
    }
    
    aln <- readAAStringSet("result/tmp/aligned.fasta") |> 
      as.matrix() |> 
      t() |> 
      as.data.frame()
    
    aln_filtered <- aln[aln$reference != "-",]
    aln_df <- aln_filtered[ , names(aln_filtered) != "reference"]
    
    ###
    if (type == "H1N1") {
      epitope <- H1N1_idx
    } else if(type == "H3N2"){
      epitope <- H3N2_idx
    } else{
      epitope <- VIC_idx
    }
    
    ###
    count_epitope_l <- list()
    epitope_idx_l <- list()
    cnt_ <- 1
    Ve <- c()
    
    for (k in colnames(aln_df)) {
        
        aln_tmp <- aln_filtered[unlist(epitope),c(k,"reference")]
        score <- ifelse(aln_tmp[,1] == aln_tmp[,2],0,1) 
        epi_idx <- aln_tmp[score == 1,] |> rownames() |> as.numeric()

        ###
        idx_sub <- epitope
        for (g in names(epitope)) {
          idx_sub[[g]] <- intersect(idx_sub[[g]],epi_idx)
        }
        
        datEpi <- sapply(idx_sub, length) |> as.data.frame()
        colnames(datEpi) <- "mutations"
        datEpi$epitopes <- rownames(datEpi)
        
        count_epitope_l[[cnt_]] <- datEpi
        epitope_idx_l[[cnt_]] <- idx_sub
        
        
        ##############################################################################
        # The simple linear model for antigenicity prediction
        # You might want to refine model later
        ##############################################################################
        Pepitope <- sum(score)/length(epitope)    
        if (type == "H1N1") {
          E <- -1.19*Pepitope+0.53

        }  else if(type == "H3N2"){
          E <- -2.47*Pepitope+0.47
          
        } else if(type == "VIC"){
          E <- -0.86*Pepitope+0.68
          
        }
        
        Ve <- c(Ve,E)
        cnt_ <- cnt_ + 1

    }
    
    
  
    list_name <- strsplit(colnames(aln_df),"\\|")
    idx <- 1
    
    vec_idx <- c()
    vec_date <- c()
    vec_id <- c()
    vec_clade <- c()
    
    for (v in list_name) {
      
      if (v[1] == "sample") {
        vec_idx <- c(vec_idx,idx)
        vec_date <- c(vec_date, as.character(input$date1))
        vec_id <- c(vec_id,v[1])
        vec_clade <- c(vec_clade,"unassigned")
        
      } else{
        if (nchar(v[3]) == 0) {
          print("pass")
        }else {
          vec_idx <- c(vec_idx,idx)
          vec_date <- c(vec_date,v[3])
          vec_id <- c(vec_id,v[1])
          vec_clade <- c(vec_clade,v[2])
          
        }
      }

      idx <- idx + 1
    }

    
    df_final <- data.frame(date=normalize_date(vec_date),E=Ve[vec_idx],type=paste0(vec_id,"|",vec_clade))

    count_epitope_filtered <- count_epitope_l[vec_idx]
    names(count_epitope_filtered) <- paste0(vec_id,"|",vec_clade)
    
    epitope_idx_filtered <- epitope_idx_l[vec_idx]
    names(epitope_idx_filtered) <- paste0(vec_id,"|",vec_clade)
    
    result <- list(ant=df_final,ecount=count_epitope_filtered,epi=epitope_idx_filtered)
    
    result
    
  })
  
  
  clicked_strain_user <- reactive({
    
    d <- event_data("plotly_click", source = "A_user")  # optional: add `source` for isolation
    if (is.null(d)) return("sample|unassigned") # default
    strain <- d$customdata
    strain
    
  })
  
  ###
  clicked_epitope_user <- reactive({
    d <- event_data("plotly_click", source = "B_user")  # optional: add `source` for isolation
    if (is.null(d)) return("A") # default
    epitope <- d$customdata[1]
    epitope
    
  })
  
  ###
  output$demo_report_user <- renderUI({
    
    req(result_ready()) 
    
    demo_df <- pred()$ant
    demo_df$group <- ifelse(demo_df$type == "sample|unassigned","samples","Circulating strains")
    
    E <- round(unique(0.53-demo_df[demo_df$type == clicked_strain_user(),2]),3)*100
    if (input$subtype_choice == "H1N1") {
      max_score = "0.53 (53%)"
    } else if (input$subtype_choice == "H3N2") {
      max_score = "0.47 (47%)"
    } else {
      max_score = "0.68 (68%)"
      }
    HTML(paste0("Vaccine efficacy of ",clicked_strain_user(),
           " compared with ", paste0(input$vaccine_strain)," strain is reduced by ",E," percent
           <br>Note: don't be alarmed by the low percentage, the maximum efficacy ie. a perfect match is ", max_score))
  })
  
  
  ##plot antigencity
  output$plot_antigencity_user <- renderPlotly({
    
    req(result_ready()) 
    demo_df <- pred()$ant
    write.csv(demo_df,"test.csv")
    demo_df$group <- ifelse(demo_df$type == "sample|unassigned","sample","Circulating strains from GISAID database")
    p <- demo_df |> ggplot(aes(date,E,
                               color=group,shape=group,
                               text = paste(
                                 "Date:", date,
                                 "<br>Name:", type,
                                 "<br>Vaccine efficacy:", round(E, 2)
                               ),customdata = type))+
      geom_point(position = position_dodge(width = .9),alpha=0.7,size=3)+
      geom_hline(yintercept=0.53, linetype="dashed", 
                 color = "red")+
      scale_colour_manual(values = c("grey","red"))+
      scale_x_date(
        date_breaks = "3 month",
        date_labels = "%Y-%m"
      )+
      labs(color="",shape="")+
      ylab("Vaccine efficacy ((u - v)/u)")+
      xlab("Date (Year-month)")+
      theme_minimal()+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            plot.title = element_text(hjust = 0.5))
    
    ggplotly(p, tooltip = "text",source = "A_user") |>
      layout(
        legend = list(
          orientation = "v",   # vertical legend
          x = 0.05,            # distance from left (0 = far left, 1 = far right)
          y = 0.05,            # distance from bottom (0 = bottom, 1 = top)
          xanchor = "left",
          yanchor = "bottom",
          bgcolor = "rgba(255,255,255,0.6)", # semi-transparent background
          bordercolor = "black",
          borderwidth = 1
        )
      )
    
  })
  
  
  output$plot_epitope_user <- renderPlotly({
    
    req(result_ready()) 
    
    strain <- clicked_strain_user()
    dx <- pred()$ecount[strain][[1]]
    print(dx)
    # Plot
    fig <- plot_ly(
      x = dx$epitopes,
      y = dx$mutations,
      customdata = dx$epitopes,
      type = 'bar',
      marker = list(color = '#4C78A8'),
      source = "B_user"
    ) %>%
      layout(
        xaxis = list(title = "Epitope"),
        yaxis = list(title = "Number of Mutations"),
        plot_bgcolor = '#ffffff',   # plot area white
        paper_bgcolor = '#ffffff',  # outside plot white
        font = list(size = 12),
        bargap = 0.3
      )
    
    fig
  })
  
  
  #Set up pairwise alignment logic
  
  alignment_prep <- reactive({
    
    req(input$subtype_choice, input$vaccine_strain)
    
    #Sample sequence
    file <- "result/consensus/draft_segment_4.fasta"
    req(file.exists(file))
    
    sample_AA <- translate_AA(readDNAStringSet("result/consensus/draft_segment_4.fasta"))
    names(sample_AA) <- "sample"
    
    #Target
    req(file.exists(file.path("data/02_vaccine_strain",paste0(input$subtype_choice, ".fasta"))))
    vacx <- readAAStringSet(paste0("data/02_vaccine_strain/",input$subtype_choice,".fasta"))
    vax_ids <-sub(">","",input$vaccine_strain)
    idx <- match(vax_ids, names(vacx))
    validate(need(!is.na(idx), paste("Selected vaccine strain (", vax_ids, ") not found")))
    reference_AA <- vacx[idx]
    names(reference_AA) <- "reference"
    
    #pairwise alignment
    alignment <- pwalign::pairwiseAlignment(
      pattern = sample_AA[[1]],
      subject = reference_AA[[1]],
      type = "global",
      substitutionMatrix = "BLOSUM62",
      gapOpening = 10,
      gapExtension = 0.5
    )
    
    aligndpair <- c(
      pwalign::alignedPattern(alignment),
      pwalign::alignedSubject(alignment)
    )
    
    names(aligndpair) <- c("User Input", "Selected Vaccine")
    
    aligndpair
    
  })
  
  #cir_strain = user-selected circulating strain sequence
  #vaccine_strain = user-selected vaccine strain
  #
  output$alignment <- renderMsaR({
    align <- alignment_prep()
    msaR(align,
         menu = TRUE,
         overviewbox = TRUE,
         colorscheme = "taylor",
         height = "500px",
         seqlogo = FALSE)
  })
  
  mismatches <- reactive({
    aligned <- alignment_prep()
    user_input <- strsplit(as.character(aligned[[1]]),"")[[1]]
    vaccine <- strsplit(as.character(aligned[[2]]),"")[[1]]
    vax_pos <- cumsum(vaccine != "-")
    changes <- vaccine[vaccine != user_input]
    if (length(changes) == 0){
      return("All amino acids identical")
    }
    paste0(vaccine[vaccine != user_input], vax_pos[vaccine != user_input], user_input[vaccine != user_input], collapse = ", ")
  })
  
  output$mismatches <- renderText({
    mismatches()
  })
  
  output$mol_user <- renderR3dmol({
    
    req(result_ready()) 

    #type <- subtype()
    epitope <- clicked_epitope_user()
    strain <- clicked_strain_user()
    
    idx <- pred()$epi[strain][[1]]
    
    
    # Read your legacy-format PDB file
    pdb_file <- file.path("data/03_tmp/",paste0(input$protein_file,".pdb"))
    pdb <- bio3d::read.pdb(pdb_file)
    pdb_trim <- trim.pdb(pdb, chain = c("A","B","C","D","E","F"))
    pdb_text <- paste(capture.output(write.pdb(pdb_trim)), collapse = "\n")  
    
    pdb_lines <- pdb_trim$atom
    pdb_text <- paste0(
      apply(pdb_lines, 1, function(row) {
        sprintf("ATOM  %5d %-4s %3s %1s%4d    %8.3f%8.3f%8.3f  1.00  0.00           %2s",
                as.integer(row["eleno"]), row["elety"], row["resid"],
                row["chain"], as.integer(row["resno"]),
                as.numeric(row["x"]), as.numeric(row["y"]), as.numeric(row["z"]),
                row["elety"])
      }),
      collapse = "\n"
    )
    

    r3dmol() %>%
      m_add_model(data = pdb_text, format = "pdb") %>% 
      m_set_style(sel = m_sel(chain = c("A","B","C","D","E","F")),
                  style = m_style_cartoon(color = "lightgray")) %>%
      m_zoom_to() %>%
      m_set_style(sel = m_sel(chain="A", resi = idx[[epitope]]),
                  style = m_style_sphere(color = "red",radius = 1)) %>%
      m_add_res_labels(
        m_sel(
          resi = idx[[epitope]],
          chain = "A"
        ),
        style = m_style_label(inFront = T,font = list(size = 6),
                              backgroundOpacity = 0.7,
                              fontColor = "white",showBackground = T,
                              alignment = "bottomRight")
      ) 
    
  })
  
  ## sunnary statistics
  output$segment_summary_usr <- renderDataTable({
    req(result_ready()) 
    mydf <- fread("result/stat/summary_statistics.tsv",header=T)
    datatable(mydf,
              extensions = 'Buttons',
              options = list(
                dom = 'Bfrtip',   # show buttons
                buttons = c('copy', 'csv', 'excel'),
                pageLength = 4
              ))
    
  })
  
  
  
  
  ##################################################################
  
  
  output$segment_summary <- renderDataTable({
    
    mydf <- data.frame(
      `genomic segment` = paste0("Segment_", 1:8),
      `Consensus Length` = sample(1000:2000, 8, replace = TRUE),
      `Read count` = sample(100:500, 8, replace = TRUE),
      `Averaged genomic depth` = sample(1000:5000, 8, replace = TRUE),
      `N number` = sample(10:50, 8, replace = TRUE)
    )
    
    
    datatable(mydf,
              extensions = 'Buttons',
              options = list(
                dom = 'Bfrtip',   # show buttons
                buttons = c('copy', 'csv', 'excel'),
                pageLength = 4
              ))
    
  })

  
  
  clicked_strain <- reactive({
    
    d <- event_data("plotly_click", source = "A")  # optional: add `source` for isolation
    if (is.null(d)) return("Demo_sample") # default
    strain <- d$customdata
    strain
    
  })
  
  output$num_epitope <- renderText({
    paste0(clicked_epitope()," epitope changes \n in ",clicked_strain())
  })
  
  
  output$demo_report <- renderText({
    
    demo_df <- fread("data/03_tmp/demo_antigenicty.tsv")
    demo_df$group <- ifelse(demo_df$type == "Demo_sample","Demo samples","Circulating strains")
    
    E <- round(unique(0.53-demo_df[demo_df$type == clicked_strain(),2]),3)*100
    
    paste0("Vaccine efficacy of ",clicked_strain(),"\n",
           "compared with A/Victoria/2570/2019 strain is reduced by ",E," percent")
  })
  
  
  clicked_epitope <- reactive({
    d <- event_data("plotly_click", source = "B")  # optional: add `source` for isolation
    if (is.null(d)) return("A") # default
    epitope <- d$customdata[1]
    epitope
    
  })
  
  
  output$plot_antigencity <- renderPlotly({
    
      demo_df <- fread("data/03_tmp/demo_antigenicty.tsv")
      demo_df$group <- ifelse(demo_df$type == "Demo_sample","Demo samples","Circulating strains from GISAID database")
      p <- demo_df |> ggplot(aes(date,E,
                                 color=group,shape=group,
                                 text = paste(
                                   "Date:", date,
                                   "<br>Name:", type,
                                   "<br>Vaccine efficacy:", round(E, 2)
                                 ),customdata = type))+
        geom_point(position = position_dodge(width = .9),alpha=0.7,size=3)+
        geom_hline(yintercept=0.53, linetype="dashed", 
                   color = "red")+
        scale_colour_manual(values = c("grey","red"))+
        scale_x_date(
          date_breaks = "3 month",
          date_labels = "%Y-%m"
        )+
        labs(color="",shape="")+
        ylab("Vaccine efficacy ((u - v)/u)")+
        xlab("Date (Year-month)")+
        theme_minimal()+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
              plot.title = element_text(hjust = 0.5))
      
      ggplotly(p, tooltip = "text",source = "A") |>
        layout(
          legend = list(
            orientation = "v",   # vertical legend
            x = 0.05,            # distance from left (0 = far left, 1 = far right)
            y = 0.05,            # distance from bottom (0 = bottom, 1 = top)
            xanchor = "left",
            yanchor = "bottom",
            bgcolor = "rgba(255,255,255,0.6)", # semi-transparent background
            bordercolor = "black",
            borderwidth = 1
          )
        )

  })
  
  
  output$plot_epitope <- renderPlotly({
    
    strain <- clicked_strain()
    epitopes <- c("A", "B", "C", "D", "E")
    mutations <- round(runif(5,min = 0,max = 5))
    
    # Plot
    fig <- plot_ly(
      x = epitopes,
      y = mutations,
      customdata = epitopes,
      type = 'bar',
      marker = list(color = '#4C78A8'),
      source = "B"
    ) %>%
      layout(
        xaxis = list(title = "Epitope"),
        yaxis = list(title = "Number of Mutations"),
        plot_bgcolor = '#ffffff',   # plot area white
        paper_bgcolor = '#ffffff',  # outside plot white
        font = list(size = 12),
        bargap = 0.3
      )
    
    fig
  })
  
  output$mol <- renderR3dmol({
    
    epitope <- clicked_epitope()
  
    # Read your legacy-format PDB file
    pdb_file <- file.path("data/03_tmp/",paste0(input$protein_file,".pdb"))
    pdb <- bio3d::read.pdb(pdb_file)
    pdb_trim <- trim.pdb(pdb, chain = c("A","B","C","D","E","F"))
    pdb_text <- paste(capture.output(write.pdb(pdb_trim)), collapse = "\n")  
  
    pdb_lines <- pdb_trim$atom
    pdb_text <- paste0(
      apply(pdb_lines, 1, function(row) {
        sprintf("ATOM  %5d %-4s %3s %1s%4d    %8.3f%8.3f%8.3f  1.00  0.00           %2s",
                as.integer(row["eleno"]), row["elety"], row["resid"],
                row["chain"], as.integer(row["resno"]),
                as.numeric(row["x"]), as.numeric(row["y"]), as.numeric(row["z"]),
                row["elety"])
      }),
      collapse = "\n"
    )
    
    
    # Visualize it
    idx <- list("A"=c(132,133,134,135),
                "B"=c(54,155,156,157,160),
                "C"=c(38,40,41,43,44,45),
                "D"=c(89,94,95,96,113,117,163),
                "E"=c(258,259,260,261,263,267))
    
    r3dmol() %>%
      m_add_model(data = pdb_text, format = "pdb") %>% 
      m_set_style(sel = m_sel(chain = c("A","B","C","D","E","F")),
                  style = m_style_cartoon(color = "lightgray")) %>%
      m_zoom_to() %>%
      m_set_style(sel = m_sel(chain="A", resi = idx[[epitope]]),
                  style = m_style_sphere(color = "red",radius = 1)) %>%
      m_add_res_labels(
        m_sel(
          resi = idx[[epitope]],
          chain = "A"
        ),
        style = m_style_label(inFront = T,font = list(size = 6),
                              backgroundOpacity = 0.7,
                              fontColor = "white",showBackground = T,
                              alignment = "bottomRight")
      )
      
  })
  
  output$protein_analysis_user <- renderUI({
    req(result_ready())
    
    selectInput(inputId = "protein_file",
                label = "Choose Protein Model to Visualise Mutations On:",
                list("2009 H1N1 influenza virus" = "3LZG",
                     "H3N2 A/Hong Kong/1/1968 influenza virus hemagglutinin" = "6CEX",
                     "H5N1 influenza virus" = "2FK0",
                     "H5N1 A/goose/Guangdong/1/96" = "4MHI",
                     "H5N1 hemagglutinin mutant (N158D, N224K and Q226L) A/Viet Nam/1203/2004"="4N5Y",
                     "H10N8 A/Jiangxi-Donghu/346/2013"="4WSX",
                     "H17N10 A/little yellow-shouldered bat/Guatemala/060/2010"="4H32",
                     "H1N1 A/South Carolina/1/1918"="4JUG",
                     "H5N1 A/duck/Laos/3295/2006(H5N1)"="4JUL",
                     "H7N9 A/Hong Kong/61/2016"="6D7C"),
                width = 999)
  })
  
  
}

# Run the app ----
shinyApp(ui = ui, server = server)
