
calculateCumulativeSizes <- function(filename) {
  # Read the file into a data frame
  ChrSize <- read.table(filename, header = FALSE)
  
  # Initialize the cumulative sizes
  ChrSize[1, 3] <- 0
  ChrSize[2, 3] <- ChrSize[1, 2]
  
  # Calculate cumulative sizes
  for (x in 2:(nrow(ChrSize) - 1)) {
    ChrSize[x + 1, 3] <- ChrSize[x, 3] + ChrSize[x, 2]
  }
  
  # Remove rows with missing values
  ChrSize <- na.omit(ChrSize)
  
  names(ChrSize)=c("Chromosome","Position_old","Position_new")
  # Return the resulting data frame
  return(ChrSize)
}


filterRawSnpCall <- function(rawsnpcall) {
  # Add a new column 'total' to rawsnpcall
  rawsnpcall <- rawsnpcall %>% mutate(total = Count_reference_allele + Count_alternative_allele)
  
  # Calculate mean and standard deviation of the 'total' column
  covpersnp <- mean(rawsnpcall$total)
  stdvcovpersnps <- sd(rawsnpcall$total)
  
  # Calculate lower and upper bounds for filtering
  lower <- covpersnp - 2 * stdvcovpersnps
  upper <- covpersnp + 2 * stdvcovpersnps
  
  # Filter rawsnpcall based on the lower and upper bounds
  rawsnpcallfiltered <- rawsnpcall %>% filter(lower < total, total < upper)
  
  return(rawsnpcallfiltered)
}



runGeneticStats <- function(data, chrSizeData) {
  data$gstat="NA"
  for (i in 1:length(data[, 1])) {
    data$gstat[i] <- G.test(x=c(data[i,3],data[i,4]),p = c(0.5,0.5))$p.value
  }
  
  snpfinalfinal <- merge(data, chrSizeData, by = c("Chromosome")) %>%
    mutate(newpos = Position + Position_new,
           gstatrunner = runner(as.numeric(gstat), f = mean, k = 10, lag = -5),
           propref = Count_reference_allele / total,
           propalt = Count_alternative_allele / total,
           proprefrunner = runner(propref, f = mean, k = 10, lag = -5),
           propaltrunner = runner(propalt, f = mean, k = 10, lag = -5))
  
  return(snpfinalfinal)
}



makePlotBSA <- function(data, ChrSize,removeUn=FALSE) {
  # Calculate LabelPos
  Chrom = ChrSize$Chromosome
  ChrSize<- ChrSize%>% mutate(LabelPos = Position_new + (Position_old / 2)) %>%
    select(LabelPos)
  vbreaksforplotLabel <- if (removeUn) as.vector(ChrSize$LabelPos[1:11]) else as.vector(ChrSize$LabelPos)
  vChr <- if (removeUn) as.vector(Chrom[1:11]) else as.vector(Chrom)
  data1 =  subset(data, Chromosome %in% vChr)
  
  # Extract breaks for x-axis
  # Create ggplot
  plot1 <- ggplot(data1) +
    geom_point(aes(x = newpos, y = -log10(gstatrunner), color = Chromosome), size=0.001) +
    theme_bw() +
    ylab("-log10(pvalue)") +
    xlab("Chromosome") +
    ylim(0, 50) +
    theme_bw() +
    labs(colour = "Chromosomes") +
    scale_x_continuous(breaks = vbreaksforplotLabel, labels = seq(1, length(vbreaksforplotLabel), 1)) + 
    scale_color_manual(values=c("black","darkgrey","black","darkgrey","black","darkgrey","black","darkgrey","black","darkgrey","black","darkgrey")) + 
    theme(legend.position = "none", axis.text = element_text(size = 20), axis.title = element_text(size = 24))
  plot1
}


makeLinePlotBSA <- function(data, ChrSize, removeUn = FALSE, setcolorline) {
  # Calculate LabelPos
  ChrSize<- ChrSize%>% mutate(LabelPos = Position_new + (Position_old / 2)) %>%
    select(LabelPos)
  
  # Extract breaks for x-axis
  vbreaksforplotLabel <- if (removeUn) as.vector(ChrSize$LabelPos[1:11]) else as.vector(ChrSize$LabelPos)
  
  # Create ggplot
  plot1 <- ggplot(data) +
    geom_line(aes(x = newpos, y = -log10(gstatrunner)), color = setcolorline) +
    theme_bw() +
    ylab("-log10(pvalue)") +
    xlab("Chromosome") +
    theme_bw() +
    ylim(0, 50) +
    labs(colour = "Chromosomes") +
    scale_x_continuous(breaks = vbreaksforplotLabel, labels = seq(1, 11, 1)) +
    theme(legend.position = "none", axis.text = element_text(size = 20), axis.title = element_text(size = 24))
  plot1
  
}


makeLinePlotBSAZoom <- function(data, ChrSize, removeUn = FALSE, setcolorline, chromosome) {
  # Calculate LabelPos
  ChrSize<- ChrSize%>% mutate(LabelPos = Position_new + (Position_old / 2)) %>%
    select(LabelPos)
  
  # Filter for the chromosome indicated by chromosome
  data = data %>% filter(Chromosome==chromosome)
  
  # Create ggplot
  plot1 <- ggplot(data) +
    geom_line(aes(x = newpos, y = -log10(gstatrunner)), color = setcolorline) +
    theme_bw() +
    ylab("-log10(pvalue)") +
    xlab(chromosome) +
    theme_bw() +
    ylim(0, 50) +
    labs(colour = "Chromosomes") +
    theme(legend.position = "none", axis.text = element_text(size = 20), axis.title = element_text(size = 24))
  plot1
  
}


makeLinePlotBSAAllelProp <- function(data, ChrSize, removeUn = FALSE, setcolorline) {
  # Calculate LabelPos
  ChrSize<- ChrSize%>% mutate(LabelPos = Position_new + (Position_old / 2)) %>%
    select(LabelPos)
  
  # Extract breaks for x-axis
  vbreaksforplotLabel <- if (removeUn) as.vector(ChrSize$LabelPos[1:11]) else as.vector(ChrSize$LabelPos)
  
  
  # Create ggplot
  plot1 <- ggplot(data) +
    geom_line(aes(x = newpos, y = proprefrunner), color = setcolorline) +
    theme_bw() +
    ylab("-log10(pvalue)") +
    xlab("Chromosome") +
    theme_bw() +
    ylim(0,1) +
    labs(colour = "Chromosomes") +
    scale_x_continuous(breaks = vbreaksforplotLabel, labels = seq(1, 11, 1)) +
    theme(legend.position = "none", axis.text = element_text(size = 20), axis.title = element_text(size = 24))
  plot1
  
}






runBSA <- function(input_file, chr_size_file) {
  # Suppress messages and warnings during function execution
  options(warn = -1)
  
  # Function 1: Read input file
  data <- read.table(input_file, h = TRUE)
  
  # Function 2: Filter raw SNP calls
  filtered_data <- filterRawSnpCall(data)
  
  # Function 3: Run genetic statistics
  final_data <- runGeneticStats(filtered_data, chr_size_file)
  return(final_data)
  
  # Reset warning settings
  options(warn = 0)
}

