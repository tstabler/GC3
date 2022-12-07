#GC3 Python output Cleaning and Processing
#Author: Thomas C. Stabler


# Introduction ------------------------------------------------------------
#Determine directory
getwd()

#Update R (if neccesary) - GC3 requires latest version of R (currently run on v4.1.1)
version

#NOTE: Current version specific to detecting gene deletions. 
#User's can adjust deletion criteria according to their specifications. Search "deletions assignments"

#Symbols < and > show where user defined input is required. Do not include < or > in command

#GC3 R packages below to function. 
#If packages not installed, delete "#" before install.packages commands below 
#Import libraries 
# install.packages("tidyr")
# install.packages("ggplot2")
# install.packages("writexl")
# install.packages("dplyr")
# install.packages("readxl")
# install.packages("reshape2")
library(tidyr)
library(ggplot2)
library(writexl)
library(dplyr)
library(tibble)
library(readxl)
library(reshape2)

# User input --------------------------------------------------------------

#Once all necessary parameters are included, select all (ctrl a) and click run above. 

# Change working directory - User input required

directory="C:/Users/stabth/Documents/UMaryland/GC3/duplications/DD2_dup"
analysis="C:/Users/stabth/Documents/UMaryland/GC3/duplications/DD2_dup"


setwd(directory)


#Sliding window variables (if sliding window used)
larger_input="samples_coverage_chrm05.txt"

#Zoom files (if every coverage position extracted)
smaller_input="samples_coverage_chrm05_zoom.txt"

#ean coverage file (if standardization required)
#If user needs to standardized, then user needs to have run mean_covearage function in Python script
#and obtained mean_coverage.txt file
standardization="chrm05_mean_cov.txt"


#Target gene coordinates - Value needs to be in integer form
gene_start=957890 
gene_end=963095

# Intron coordinates (i.e. select for CDS positions only) - must be within target gene coordinates
# NOTE: GC3 currently currently marks coordinates of 1 intron on figures
# If no intron, then keep intron_start and intron_end to 0
# If more than 1 intron in gene, then additional conditions required ("Generate meta-database" section). 
# For example -> subset(gene_pos, (Position<intron_end | Position>intron_start | Position<intron2_end | Position>intron2_start))
intron_start=0
intron_end=0

#Flanking gene coordinates (optional)
#NOTE: make sure start/end positions are within range of coordinates provided in python script!
#If flanking gene information not desired, leave all variables = 0
downstream_start=0
downstream_end=0

upstream_start=0
upstream_end=0

# Zoom (e.g. exon) - if user wants to zoom on a particular area of interest (e.g. exon)
#If information not desired, leave variables = 0
zoom_start=0
zoom_end=0

#Subgroup sample set (optional) - Currently set to subset 2 groups
#Include all names of 1 subgroup in sample_change1 list ->
#(GC3 will automatically assign other samples to other subgroup)
#If more subgroups needed, then editing is required under the "GROUPED" sections below  
#Uncomment lines below to include subgroups (delete "#" at start of lines)
# sample_change1 <- c("<country_samplesID_1>", "<country_sampleID_2>",
#                        )
# subgroup1 <- "Group1"
# subgroup2 <- "Group2"

#When setting groups - plots visuals can vary - use limity variable to find best plot view
#limity value is for non-standardized plot and limity_stan is for standardized plot
#Highlight and Rerun only "GROUPED" sections for faster results.
# limity=50
# limity_stan=4

# Data cleaning/preparation - sliding window (i.e. larger_input) -----------------------------------------------
#Input files - sliding window (e.g. 500bp with 100bp jump)
if(exists("larger_input"))
{
  input1 <- readLines(larger_input)
  
  #Get rid of { and } and quotations
  input1 <- gsub("\\{", "", input1)
  input1 <- gsub("\\}", "", input1)
  input1 <- gsub("\\'", "", input1)
  
  #Separate using comma delimiter
  input1 <- unlist(strsplit(input1, ","))
  
  #Create data frame from txt file
  input2 <- as.data.frame(input1)
  
  #Create dataframe split into 3 columns by start, end and coverage values
  input3 <- input2 %>%
    separate(input1, c("Start_position", "End_position", "Coverage"), ":")
  
  #Using Start_position and end_position as template, move coverages to separate column with country/sample as header name
  #Extract names (keep columns with NA in End_position column - should be NA for only sample name rows)
  sample_names <- input3 %>% filter(is.na(End_position))
  
  sample_names <- sample_names %>% select(-one_of('End_position', 'Coverage'))
  colnames(sample_names)[colnames(sample_names)=="Start_position"] <- "Sample_name"
  
  sample_names <- gsub(".*samples_(.+)_alignments.*", "\\1", sample_names$Sample_name)
  
  
  sample_names
  
  #Delete rows that aren't integers in input3
  input4 <- na.omit(input3)
  
  #Get number of distinct observations in Position column to determine where to split dataframe
  
  unique_rows <- length(unique(input4[["Start_position"]]))
  
  #Create new column every 1528 rows (NOTE: this is based on predescribed methods - this may require adjustment depending on start/stop positions)
  input4 <- as.data.frame(matrix(input4$Coverage, nrow=unique_rows))
  
  #rename header names
  colnames(input4) <- sample_names
  
  #Create input5 that gets rid of Coverage column
  input5 <- input3 %>% select(-one_of('Coverage'))
  input5 <- na.omit(input5)
  input5[, c(1,2)] <- apply(input5[, c(1,2)], 2,
                            function(x) as.numeric(as.character(x)))
  #Get rid of duplicates
  input5 <- distinct(input5)
  
  #Combine input4 and input5 columns
  input6 <- cbind(input5, input4)
  
  #ake all values in input files numeric
  input6[] <- lapply(input6, as.numeric)
  input4[] <- lapply(input4, as.numeric)
  
  rm(input1, input2, input3, input5)
  
  #Write resulting file to excel
  write_xlsx(input6, "region_pos_only.xlsx", col_names = T, format_headers = T)
}

# target gene file - magnified file (coverage of every position)  --------------------------------------------------------------
#Import file
if(exists("smaller_input"))
{
  input1_zoom <- readLines(smaller_input)
  
  #Get rid of { and } and quotations
  input1_zoom <- gsub("\\{", "", input1_zoom)
  input1_zoom <- gsub("\\}", "", input1_zoom)
  input1_zoom <- gsub("\\'", "", input1_zoom)
  
  #Separate using comma delimiter
  input1_zoom <- unlist(strsplit(input1_zoom, ","))
  
  #Create data frame from txt file
  input2_zoom <- as.data.frame(input1_zoom)
  
  #Create dataframe split into 3 columns by start, end and coverage values
  input3_zoom <- input2_zoom %>%
    separate(input1_zoom, c("Position", "Coverage"), ":")
  
  # #Using "Position" as template, move coverages to separate column with country/sample as header name
  # #Extract names
  sample_names <- input3_zoom %>% filter(is.na(Coverage))
  
  sample_names <- sample_names %>% select(-one_of('Coverage'))
  colnames(sample_names)[colnames(sample_names)=="Position"] <- "Sample_name"
  
  #Extract word(s) between keywords samples and alignments
  sample_names <- gsub(".*samples_(.+)_alignments.*", "\\1", sample_names$Sample_name)
  
  sample_names <- gsub("samples_", "", sample_names)
  sample_names
  
  #Delete rows that aren't integers in input3
  input4_zoom <- na.omit(input3_zoom)
  
  #Get number of unique rows and split dataframe
  unique_rows_zoom <- length(unique(input4_zoom[["Position"]]))
  
  #Create new column based on unique rows 
  input4_zoom <- as.data.frame(matrix(input4_zoom$Coverage, nrow=unique_rows_zoom))
  
  #rename header names
  colnames(input4_zoom) <- sample_names
  
  #Create input5 that gets rid of Coverage column
  input5_zoom <- na.omit(input3_zoom)
  input5_zoom <- input5_zoom %>% select(-one_of('Coverage'))
  input5_zoom[] <- lapply(input5_zoom, as.numeric)
  
  #Get rid of duplicates
  input5_zoom <- distinct(input5_zoom)
  
  #Combine input4 and input5 columns
  input6_zoom <- cbind(input5_zoom, input4_zoom)
  
  #ake all values in input files numeric
  input6_zoom[] <- lapply(input6_zoom, as.numeric)
  input4_zoom[] <- lapply(input4_zoom, as.numeric)
  
  rm(input1_zoom, input2_zoom, input3_zoom, input5_zoom)
  
  #Subset data to only target gene positions
  gene_pos <- input6_zoom[ which(input6_zoom$Position>=gene_start & input6_zoom$Position<=gene_end) , ]
  
  #Write resulting gene_pos file to excel
  write_xlsx(gene_pos, "gene_pos_only.xlsx", col_names = T, format_headers = T)
}

# Standardization - process mean coverage file ----------------------------
if(exists("standardization"))
{
  input_stan <- readLines(standardization)
  
  #Get rid of { and } and quotations
  input_stan <- gsub("\\{", "", input_stan)
  input_stan <- gsub("\\}", "", input_stan)
  input_stan <- gsub("\\'", "", input_stan)
  
  
  #Separate using comma delimiter
  input_stan <- unlist(strsplit(input_stan, ","))
  
  #Create data frame from txt file
  input2_stan <- as.data.frame(input_stan)
  
  #Create dataframe split into 3 columns by start, end and coverage values
  input3_stan <- input2_stan %>%
    separate(input_stan, c("Start Position", "End Position", "Mean coverage"), ":")
  
  #Get rid of rows with sample names
  input4_stan <- input3_stan[-grep("_", input3_stan$`Start Position`),]
  
  #Cbind sample_names (from magnified file) to input4
  data_reshape <- cbind(sample_names, input4_stan)
  
  data_reshape[,2:ncol(data_reshape)] <- lapply(data_reshape[,2:ncol(data_reshape)], as.numeric)
  
  write_xlsx(data_reshape,"samples_mean_coverage.xlsx", col_names = T, format_headers = T)
  
  rm(input_stan, input2_stan, input3_stan, input4_stan)
}

# Additional coverage files (if optional inputs included) ---------------------------------------------------

#Create subsets for flanking genes (downstream and upstream)
if(downstream_start>0 | upstream_start>0)
{
  flankgene_down <- input6_zoom[ which(input6_zoom$Position>=downstream_start & input6_zoom$Position<=downstream_end) , ]
  flankgene_up <- input6_zoom[ which(input6_zoom$Position>=upstream_start & input6_zoom$Position<=upstream_end) , ]
}

#Standardized target gene positions
if(exists("data_reshape"))
{
  gene_pos_stan <- data.frame(rbind(gene_pos[,2:ncol(gene_pos), drop=F], data_reshape$`Mean coverage`))
  gene_pos_stan <- data.frame(lapply(gene_pos_stan, function(x) x/x[nrow(gene_pos_stan)]))
  gene_pos_stan <- head(gene_pos_stan, -1)

  gene_pos_stan <- cbind(gene_pos[,1,drop=F],gene_pos_stan)
  
  #Save intermediate file
  write_xlsx(gene_pos_stan, "gene_pos_standardized.xlsx", col_names = T, format_headers = T)
}

#Subset data according to smaller areas of interest (e.g. exon)
if(zoom_start>0)
{
  gene_zoom <- gene_pos[ which(gene_pos$Position>=zoom_start & gene_pos$Position<=zoom_end) , ]
}
#Create target gene CDS, excluding intron region
if(intron_start>0)
{
  gene_cds <- subset(gene_pos, (Position>intron_end | Position<intron_start))
} else {
  gene_cds <- gene_pos
}

# Generate meta-database ----------------------------------------------------------

meta_data <- data.frame(sample_names)
#Create country column (NOTE: this is based on format of "<country of origin>_<sample name>")
meta_data$country <- sub("_.*", "", meta_data$sample_names)
#Create sample column
meta_data$sample <- sub(".*?_", "", meta_data$sample_names)

#Create column of mean coverage
if(exists("data_reshape"))
{meta_data$regional_mean <- data_reshape[,2,drop=F]
}
meta_data$smaller_mean <- sapply(gene_pos[,2:ncol(gene_pos),drop=F],2,FUN=mean)

#Calculate proportional coverage of target gene CDS
meta_data$gene_cds_cov <- apply(gene_cds[,2:ncol(gene_cds),drop=F],2,function(x) mean(x>0))

#replace any NA of CDS column with 0
meta_data$gene_cds_cov[is.na(meta_data$gene_cds_cov)] <- 0

#Create column with deletions assignments
meta_data <- meta_data %>% mutate(deletion = case_when(gene_cds_cov>=0.75 ~ "No deletion",
                                                       gene_cds_cov<0.75 & gene_cds_cov>0 ~ "Parital deletion",
                                                       gene_cds_cov==0 ~ "Complete deletion"))

#Create counts of target gene positions with zero coverage
meta_data$del_0X_count <- apply(gene_pos[,2:ncol(gene_pos),drop=F],2,function(x) sum(x==0))

#Create column of proportion coverage of zoom positions
if(zoom_start>0)
{
  meta_data$gene_zoom_cov <- apply(gene_zoom[,2:ncol(gene_zoom),drop=F],2,function(x) mean(x>0))
}

#Create column(s) of mean coverage for upstream/downstream flanking genes
if(downstream_start>0)
{
  meta_data$down_flank_mean <- sapply(flankgene_down[,2:ncol(flankgene_down),drop=F],FUN=mean)
  meta_data$down_flank_cov <- apply(flankgene_down[,2:ncol(flankgene_down),drop=F],2,function(x) mean(x>0))
}
if(upstream_start>0)
{
  meta_data$up_flank_mean <- sapply(flankgene_up[,2:ncol(flankgene_up),drop=F],FUN=mean)
  meta_data$up_flank_cov <- apply(flankgene_up[,2:ncol(flankgene_up),drop=F],2,function(x) mean(x>0))
}

#Create column if deletions extending into upstream or downstream regions around target gene
if(downstream_start>0 | upstream_start>0)
{
  meta_data <- meta_data %>% mutate(gene_flank = case_when(down_flank_cov>=0.75 & up_flank_cov>=0.75 ~ "No flanking deletions",
                                                           down_flank_cov<0.75 & up_flank_cov>=0.75 ~ "Downstream deletion, no upstream deletion",
                                                           down_flank_cov>=0.75 & up_flank_cov<0.75 ~ "Upstream deletion, no downstream deletion",
                                                           down_flank_cov<0.75 & up_flank_cov<0.75 ~ "Upstream and downstream deletion"))
}
#Generate dataframe of each position with counts of positions w/ and without coverage (x>0)
cov_count1 <- apply(gene_pos[,2:ncol(gene_pos),drop=F],1,function(x) sum(x > 0))
cov_count0 <- apply(gene_pos[,2:ncol(gene_pos),drop=F],1,function(x) sum(x==0))
mean_cov <- apply(gene_pos[,2:ncol(gene_pos),drop=F],1,function(x) mean(x))
median_cov <- apply(gene_pos[,2:ncol(gene_pos),drop=F],1,function(x) median(x))

cov_count <- data.frame(gene_pos[,1,drop=F],cov_count1,cov_count0,mean_cov,median_cov)

#Write meta_data to excel sheet
write_xlsx(meta_data, "meta_data.xlsx", col_names = T, format_headers = T)
write_xlsx(cov_count,"gene_pos_descriptives.xlsx", col_names = T, format_headers = T)


# Regional plot - sliding window ------------------------------------------------------
if(exists("larger_input"))
{
  #elt sliding window files
  cov_melt <- melt(input6, id.vars = c('Start_position','End_position'), variable.name = 'samples')
  
  #Log scale all values
  cov_melt[,4] <- log10(cov_melt[,4])
  
  #Substitute all -Inf or negative number for 0
  cov_melt[cov_melt=="-Inf"] <-0
  cov_melt[cov_melt<0] <- 0
  cov_melt[cov_melt=="NaN"] <- 0
  
  #Plot melted database
  cov_plot <- ggplot(cov_melt, aes(x=End_position, y=value,colour=samples)) +
    geom_line(size=0.5) +
    geom_vline(xintercept=gene_start, color="black", size=1) +
    geom_vline(xintercept=gene_end, color="black", size=1) +
    theme(axis.title.x = element_text(size=15),
          axis.title.y = element_text(size=15),
          # legend.text = element_text(size=15),
          # legend.title = element_text(size=15),
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(size=1.5),
          axis.ticks.length = unit(0.1, "inch"),
          aspect.ratio = 1) +
    labs(x = "Chromosomal coordinates", y = bquote(log[10](Coverage)), colour = "Samples") +
    guides(color = guide_legend(override.aes = list(size=4)))
  
  cov_plot
  
  #SAVE
  ggsave("Regional_sliding_window.jpeg", 
         path = analysis,
         width = 7,
         height = 7,
         units = "in")
}
# target gene loci plot - UNGROUPED/UNSTANDARDIZED -------------------------------------------------------
if(exists("smaller_input"))
{
  #Melt target positions
  melt_gene <- melt(gene_pos, id.vars = 'Position', variable.name = 'samples')
  
  #Replace Na with 0 values
  melt_gene$value[is.na(melt_gene$value)] <- 0
  
  #gene plot
  plot_gene <- ggplot(melt_gene, aes(x=Position, y=value,colour=samples)) +
    geom_line() +
    labs(x = "target gene locus", y = "Coverage")
  if(intron_start>0&intron_end>0)
  {
    plot_gene <- plot_gene +
      geom_rect(aes(xmin = intron_start, xmax = intron_end, ymin=0, ymax=Inf), fill="tan1", color=NA, alpha=0.01) +
      theme(axis.title.x = element_text(size=15),
            axis.title.y = element_text(size=15),
            # legend.text = element_text(size=12),
            # legend.title = element_text(size=15),
            legend.position = "none",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.ticks = element_line(size=1.5),
            axis.ticks.length = unit(0.1, "inch"),
            aspect.ratio = 1) +
      guides(fill="none", colour=guide_legend(override.aes = list(size=6)))
  
  
    plot_gene
  
    #SAVE
    ggsave("gene_pos_plot_ungrouped_unstandardized.jpeg", 
           path = analysis,
           width = 7,
           height = 7,
           unit = "in")
  
    rm(melt_gene, plot_gene)
  }
  if(intron_start==0&intron_end==0)
  {
    plot_gene <- plot_gene +
      theme(axis.title.x = element_text(size=15),
            axis.title.y = element_text(size=15),
            # legend.text = element_text(size=12),
            # legend.title = element_text(size=15),
            legend.position = "none",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.ticks = element_line(size=1.5),
            axis.ticks.length = unit(0.1, "inch"),
            aspect.ratio = 1) +
      guides(fill="none", colour=guide_legend(override.aes = list(size=6)))
    
    
    plot_gene
    
    #SAVE
    ggsave("gene_pos_plot_ungrouped_unstandardized.jpeg", 
           path = analysis,
           width = 7,
           height = 7,
           unit = "in")
    
    rm(melt_gene, plot_gene)
  }
}
# target gene loci plot - GROUPED/UNSTANDARDIZED -------------------------------------------------------

if(exists("smaller_input") & exists("sample_change1"))
{
  
  #elt target gene positions
  melt_gene <- melt(gene_pos, id.vars = 'Position', variable.name = 'samples')
  
  #Replace Na with 0 values
  melt_gene$value[is.na(melt_gene$value)] <- 0
  
  #Group sample according to user input of sample names. If not on sample name or similar patter and relabel with group1 or group2 name
  melt_gene$samples <- ifelse(melt_gene$samples%in%sample_change1, subgroup1, subgroup2) 
  
  #gene plot
  plot_gene <- ggplot(melt_gene, aes(x=Position, y=value,colour=samples)) +
    geom_line(stat = "summary", fun = "median") +
    labs(x = "target gene locus", y = "edian Coverage",
         color = "Group") +
    scale_color_manual(name = "Group", 
                       labels=c("Group 1", "Group 2"),
                       values = c("#009E73", "#CC79A7"))
  
  if(intron_start>0&intron_end>0)
  {
    plot_gene <- plot_gene +
      geom_rect(aes(xmin = intron_start, xmax = intron_end, ymin=0, ymax=Inf), fill="tan1", color=NA, alpha=0.01) +
      theme(axis.title.x = element_text(size=15),
            axis.title.y = element_text(size=15),
            legend.text = element_text(size=12),
            legend.title = element_text(size=15),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.ticks = element_line(size=1.5),
            axis.ticks.length = unit(0.1, "inch"),
            aspect.ratio = 1) +
      guides(fill="none", colour=guide_legend(override.aes = list(size=6))) +
      coord_cartesian(ylim=c(0,limity))
    
    
    plot_gene
    
    #SAVE
    ggsave("gene_pos_plot_grouped_unstandardized.jpeg", 
           path = analysis,
           width = 7,
           height = 7,
           unit = "in")
    
    rm(plot_gene, melt_gene)
  }
  
  if(intron_start==0&intron_end==0)
  {
    plot_gene <- plot_gene +
      theme(axis.title.x = element_text(size=15),
            axis.title.y = element_text(size=15),
            legend.text = element_text(size=12),
            legend.title = element_text(size=15),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.ticks = element_line(size=1.5),
            axis.ticks.length = unit(0.1, "inch"),
            aspect.ratio = 1) +
      guides(fill="none", colour=guide_legend(override.aes = list(size=6))) +
      coord_cartesian(ylim=c(0,limity))
    
    
    plot_gene
    
    #SAVE
    ggsave("gene_pos_plot_grouped_unstandardized.jpeg", 
           path = analysis,
           width = 7,
           height = 7,
           unit = "in")
    
    rm(plot_gene, melt_gene)
  }
}

# target gene loci plot - UNGROUPED/STANDARDIZED -------------------------------------------------------
if(exists("smaller_input") & exists("gene_pos_stan"))
{
  #elt standardized target gene positions
  melt_gene <- melt(gene_pos_stan, id.vars = 'Position', variable.name = 'samples')
  
  #Replace Na with 0 values
  melt_gene$value[is.na(melt_gene$value)] <- 0
  
  #gene plot
  plot_gene <- ggplot(melt_gene, aes(x=Position, y=value,colour=samples)) +
    geom_line() +
    labs(x = "target gene locus", y = "Standardized Coverage")
  if(intron_start>0&intron_end>0)
  {
    plot_gene <- plot_gene +
      geom_rect(aes(xmin = intron_start, xmax = intron_end, ymin=0, ymax=Inf), fill="tan1", color=NA, alpha=0.01) +
      theme(axis.title.x = element_text(size=15),
            axis.title.y = element_text(size=15),
            # legend.text = element_text(size=12),
            # legend.title = element_text(size=15),
            legend.position = "none",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.ticks = element_line(size=1.5),
            axis.ticks.length = unit(0.1, "inch"),
            aspect.ratio = 1) +
      guides(fill="none", colour=guide_legend(override.aes = list(size=6)))
    
    
    plot_gene
    
    #SAVE
    ggsave("gene_pos_plot_ungrouped_standardized.jpeg", 
           path = analysis,
           width = 7,
           height = 7,
           unit = "in")
    
    rm(plot_gene, melt_gene)
  }
  
  if(intron_start==0&intron_end==0)
  {
    plot_gene <- plot_gene +
      theme(axis.title.x = element_text(size=15),
            axis.title.y = element_text(size=15),
            # legend.text = element_text(size=12),
            # legend.title = element_text(size=15),
            legend.position = "none",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.ticks = element_line(size=1.5),
            axis.ticks.length = unit(0.1, "inch"),
            aspect.ratio = 1) +
      guides(fill="none", colour=guide_legend(override.aes = list(size=6)))
    
    
    plot_gene
    
    #SAVE
    ggsave("gene_pos_plot_ungrouped_standardized.jpeg", 
           path = analysis,
           width = 7,
           height = 7,
           unit = "in")
    
    rm(plot_gene, melt_gene)
  }
}
# target gene loci plot - GROUPED/STANDARDIZED -------------------------------------------------------
if(exists("smaller_input") & exists("sample_change1") & exists("standardization"))
{
  
  #elt target gene positions
  melt_gene_group <- melt(gene_pos_stan, id.vars = 'Position', variable.name = 'samples')
  
  #Replace Na with 0 values
  melt_gene_group$value[is.na(melt_gene_group$value)] <- 0
  
  #Group sample according to user input of sample names. If not on sample name or similar patter and relabel with group1 or group2 name
  melt_gene_group$samples <- ifelse(melt_gene_group$samples%in%sample_change1, subgroup1, subgroup2) 
  
  #gene plot
  plot_gene <- ggplot(melt_gene_group, aes(x=Position, y=value,colour=samples)) +
    geom_line(stat = "summary", fun = "median") +
    labs(x = "target gene locus", y = "edian Standardized Coverage",
         color = "Group") +
    scale_color_manual(name = "Group", 
                       labels=c("Group 1", "Group 2"),
                       values = c("#009E73", "#CC79A7"))
  
  if(intron_start>0&intron_end>0)
  {
    plot_gene <- plot_gene +
      geom_rect(aes(xmin = intron_start, xmax = intron_end, ymin=0, ymax=Inf), fill="tan1", color=NA, alpha=0.01) +
      theme(axis.title.x = element_text(size=15),
            axis.title.y = element_text(size=15),
            legend.text = element_text(size=12),
            legend.title = element_text(size=15),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.ticks = element_line(size=1.5),
            axis.ticks.length = unit(0.1, "inch"),
            aspect.ratio = 1) +
      guides(fill="none", colour=guide_legend(override.aes = list(size=6))) +
      coord_cartesian(ylim=c(0,limity_stan))
    
    
    plot_gene
    
    #SAVE
    ggsave("gene_pos_plot_grouped_standardized.jpeg", 
           path = analysis,
           width = 7,
           height = 7,
           unit = "in")
    
    rm(plot_gene, melt_gene_group)
    
  }
  
  if(intron_start==0&intron_end==0)
  {
    plot_gene <- plot_gene +
      theme(axis.title.x = element_text(size=15),
            axis.title.y = element_text(size=15),
            legend.text = element_text(size=12),
            legend.title = element_text(size=15),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.ticks = element_line(size=1.5),
            axis.ticks.length = unit(0.1, "inch"),
            aspect.ratio = 1) +
      guides(fill="none", colour=guide_legend(override.aes = list(size=6))) +
      coord_cartesian(ylim=c(0,limity_stan))
    
    
    plot_gene
    
    #SAVE
    ggsave("gene_pos_plot_grouped_standardized.jpeg", 
           path = analysis,
           width = 7,
           height = 7,
           unit = "in")
    
    rm(plot_gene, melt_gene_group)
  }
}

# Zero coverage positions count -------------------------------------------
if(exists("cov_count"))
{
  #Filter for only cambodia samples
  stacked <- cov_count %>% select(c(contains("Position") | contains("count")))
  
  #elt
  stacked_melt_gene <- melt(stacked, id.vars = 'Position', variable.name = 'count')
  
  #Plot
  stacked_plot <- ggplot(stacked_melt_gene, aes(x=Position,y=(value/sum(value)), fill=count)) +
    geom_bar(position = "fill", stat="identity", width = 1) + 
    labs(x = "target gene locus", y = "Sample Proportion") + 
    scale_fill_manual(name = "Position Coverage", labels=c("\u22651X Coverage", "Zero Coverage"), values=c("lightsteelblue1", "firebrick4"))
  
  if(intron_start>0)
  {
    stacked_plot <- stacked_plot +
      geom_rect(aes(xmin = intron_start, xmax = intron_end, ymin=0, ymax=1), fill="tan1", alpha=0.01)
  }
  
  
  stacked_plot <- stacked_plot + scale_y_continuous(labels=scales::percent) +
    theme(axis.title.x = element_text(size=15),
          axis.title.y = element_text(size=15),
          legend.text = element_text(size=15),
          legend.title = element_text(size=15),
          aspect.ratio = 1)
  
  
  stacked_plot
  
  #SAVE
  ggsave("gene_coverage_count.jpeg", path = analysis)
}