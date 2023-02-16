# Author: Arvis Sulovari
# Date last edited: 02/19/2020
# Modifications are made to the original script for adaptability.
# Goal: convert color-coded FASTA sequences into MSA_style sequence composition plots


# Log stuff
log <- file(snakemake@log[[1]], open="wt")

# Define arguments (only 1 -- the region ID)
args <- commandArgs(TRUE)
REGION_NUM <- as.numeric(as.character(toString(args[1])))

# All haplotype names:
all_haps <- read.table("put/path/HERE/MEI", header = T)
all_haps <- as.character(all_haps$samples_haplo)
#get rid of ones in the meantime that are misbehaving
haps_keep <- read.table("./5ac_samples_unique", sep = "\t", header = T)
haps_keep <- as.character(haps_keep$sample)
all_haps <- all_haps[all_haps %in% haps_keep]


# Read in bed file with chr:start-end format region IDs. This will be placed as simple text at the bottom of plots.
all_strvntr <- snakemake@params[["motif"]]

top_3_colors <- c("#9E0142","#A90D44","#B41947")

# Function for determining the point where to cut the color strip

sum_neighours <- function(myarr=c(1:10),nr_neigh=2){
  sum_arr <- c()
  for (i in 1:(length(myarr)-(nr_neigh))){
    j=i+nr_neigh-1
    sum_arr <- c(sum_arr,sum(myarr[i:j]))
    i=i+nr_neigh-1
  }
  return(sum_arr)
}

###check longest pure track for each haplotype file
lpt_files <- vector()
for (i in 1:length(all_haps)){
  tryCatch({curr_lpt <- read.delim(paste0("../fasta/", all_haps[i],"_MUC5AC_VNTRflanks.fasta.lpt"))}, error=function(e){})
  lpt_files <- c(lpt_files, curr_lpt[,1])
}
names(lpt_files) <- all_haps
lpt_files <- lpt_files[rev(order(lpt_files))]

#make list of data for samples and check to make sure they aren't empty!
kmer_data <- list()
for (i in 1:length(all_haps)){
  tryCatch({curr_fa <- read.delim(paste0("../fasta2/", all_haps[i], "_MUC5AC_VNTRflanks.fasta"), header = T)}, error=function(e){})
  curr_fa <- as.character(curr_fa[,1])
  kmer_data[[i]] <- curr_fa
}
names(kmer_data) <- all_haps


info_data <- list()
for (i in 1:length(all_haps)){
  curr_info <- file.info(paste0("../fasta2/", all_haps[i],"_MUC5AC_VNTRflanks.fasta"))
  info_data[[i]] <- curr_info
}
names(info_data) <- all_haps



###order samples by size in info_data
sizes <- vector()
for (i in 1:length(all_haps)){
  curr_hap <- info_data[all_haps[i]]
  curr_hap <- curr_hap[[1]]$size
  sizes <- c(sizes, curr_hap)
}
names(sizes) <- all_haps
sizes <- sizes[rev(order(sizes))]

largest_file <- names(sizes[1])
largest_file_size <- unname(sizes[1])


###sort data by size order in lpt vector
myind <- unname(which((sort(sizes,decreasing = T))>0))
re_sorted_names <- names(sizes)
re_sorted_indeces <- c(order(sizes,decreasing = T))[myind]
size_arr <- unname(sizes)


##bring in VNTR coloring data
hprc_trf <- read.table("/net/eichler/vol27/projects/mucin/nobackups/airway_mucins_FIXED/MUC1/all_trf_compiled.txt", header = T)
hprc_trf <- hprc_trf[which(hprc_trf$PeriodSize == 60),]
hprc_trf <- hprc_trf[which(hprc_trf$CopyNumber >= 3.3),]

remove <- c("353", "357", "806", "830", "892", "895", "803", "866", "871", "879", "885", "797", "900")
for (i in 1:length(remove)){
    hprc_trf <- hprc_trf[-which(rownames(hprc_trf) == remove[i]),]
}

hprc_trf$primate <- c(rep("human", 101), rep("nhp", 4))

hprc_trf <- hprc_trf[which(hprc_trf$haplotype %in% names(kmer_data)),]

# specify colors for haplogroup
mycolors_tmp <- colorRampPalette(c("cornflowerblue", "chartreuse3", "darkorchid4"))(length(unique(hprc_trf$CopyNumber)))
colors_df <- as.numeric(names(table(hprc_trf$CopyNumber)))
colors_df <- data.frame(length = colors_df, color = mycolors_tmp)


color_vec <- vector()
for(i in 1:length(re_sorted_names)){
   curr_hap <- re_sorted_names[i]
   curr_length <- hprc_trf$CopyNumber[which(hprc_trf$haplotype == curr_hap)]
   curr_color <- as.character(colors_df$color[which(colors_df$length == curr_length)])
   names(curr_color) <- curr_hap
   color_vec <- c(color_vec, curr_color)
}


# START PLOTTING HERE
pdf(file = snakemake@output[['waterfall_pdf']], width = 18,height = 14)
par(mfrow=c(length(re_sorted_names),1))
par(mar=c(0,16,0,2))

tryCatch({  
  for(ind_file in re_sorted_names){
    print(ind_file)
    haplotype_name <- all_haps

    # Convert file names to actual array of colors, according to the content of column V1
    col_arr <- unname(unlist(kmer_data[ind_file]))
    
    # Read in files wiht the longest pure track values
    lpt <- unname(lpt_files[ind_file])


    current_index <- which(re_sorted_names==ind_file)
    current_index2 <- re_sorted_indeces[current_index]

    # Find position in the array where to insert the gap:
    ## Find all positions where the color is one of the top three colors (*NOTE: consider adjusting this based on the KMER size)
    col_indeces <- c(which(col_arr %in% top_3_colors))

    ## Find the last occurrence of *** 3 consecutive *** hot colors
    col_arr_tmp <- c(col_indeces[-1]-col_indeces[-length(col_indeces)])
    consec_difs <- sum_neighours(myarr = col_arr_tmp,nr_neigh = 2)
    index_to_cut <- which(consec_difs<=49)[length(which(consec_difs<=49))]+1
    index_to_cut <- col_indeces[index_to_cut]

    # Add both left and right alignments

    
    # Calculate gap size.
    #gap_size <- largest_file_size-length(col_arr)
    gap_size <- 370-length(col_arr)
    # Redefine left_coords and right_coords to match tht end of the VNTR region
    left_coords <- c(1:index_to_cut)  
    right_coords <- seq(index_to_cut+1,length(col_arr))

    # Plot each half of the barplots separartely, separated by the gap_size. Fix gap-size so that we break at the end of the VNTR
    # Left half of the plot
    barplot(c(rep(1,length(left_coords)),rep(0,length(left_coords))),col=as.character(col_arr),
            border = NA,space = 0,width = 1,xlim = c(0,largest_file_size),axes="FALSE")

    # Right half of the plot
    barplot(c(rep(0,length(left_coords)+gap_size),rep(1,length(right_coords))),col=c(rep("#ffffff",gap_size),as.character(col_arr)),
            border = NA,space = 0,width = 1,xlim = c(0,largest_file_size),axes="FALSE",xpd = T,add=T)
    
    # Add labels of each haplotype
    axis(side=2,at=0.5,labels = ind_file, cex.lab=1, cex.axis=1, cex.main=1, 
         cex.sub=1,tick = F,las=1,hadj = 0.7)

    #axis(side=2,at=0.5,labels = ind_file,col.axis=unname(color_vec[ind_file]),cex.lab=1, cex.axis=1, cex.main=1, 
         #cex.sub=1,tick = F,las=1,hadj = 0.7)
    
    axis(side=4,at=0.5,labels = lpt,cex.lab=1, cex.axis=1, cex.main=1, 
         cex.sub=3,tick = F,las=1,hadj = 1.5)
  }

}, error=function(e){})


tryCatch({

  region_name=toString(all_strvntr[REGION_NUM+1,]$REGION)
  motif_seq=toString(all_strvntr[REGION_NUM+1,]$MOTIF_SEQ)
  text(x = largest_file_size*0.5,y=0.7,labels = region_name,cex=1)
  text(x = largest_file_size*0.5,y=0.4,labels = paste0("MOTIF: ",motif_seq),cex=1)
}, error=function(e){})

  dev.off() 
