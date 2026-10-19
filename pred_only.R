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
options(shiny.maxRequestSize = 1024 * 1024 * 2048)

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



#The code I used to test the translation function
HA_seg <- readDNAStringSet("data/05_example_genes/A:California:07:2009(H1N1).fasta") |> translate_AA()
as.character(HA_seg)