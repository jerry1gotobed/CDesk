suppressMessages(library(ggplot2))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(viridis))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ggsignif))
suppressMessages(library(VennDiagram))
suppressMessages(library(ggpointdensity))
suppressMessages(library(stringr))
suppressMessages(library(PRIMME))

sparse2full <- function(sparse.mat,resolution){
  sparse.mat <- as.matrix(sparse.mat)
  bins <- unique(c(sparse.mat[, 1], sparse.mat[, 2]))
  bins <- as.numeric(bins)
  bins <- bins[order(bins)]
  #bin.size <- min(diff(bins))
  bin.size = resolution
  min.bin <- min(bins)
  max.bin <- max(bins)
  # set sequence of rows/cols for full matrix
  cols <- seq(min.bin, max.bin, by = bin.size)
  if (tail(cols, 1) != max.bin) {
    cols <- c(cols, max.bin)
  }
  # set up matrix to be correct size and row/col names
  mat <- matrix(nrow = length(cols), ncol = length(cols))
  rownames(mat) <- cols
  colnames(mat) <- cols
  
  sparse.mat[,1] <- match(sparse.mat[,1], cols)
  sparse.mat[,2] <- match(sparse.mat[,2], cols)
  # populate matrix
  rcix <- sparse.mat[, 1:2]
  mat[rcix] <- sparse.mat[,3]
  rcix <- rcix[,c(2,1)]
  mat[rcix] <- sparse.mat[,3]
  
  mat[is.na(mat)] <- 0
  message(paste0("Matrix dimensions: ", dim(mat)[1], "x", dim(mat)[2]))
  return(mat)
}

TADCompare = function(cont_mat1,
                      cont_mat2,
                      resolution = "auto",
                      z_thresh = 2,
                      window_size = 15,
                      gap_thresh = .2,
                      pre_tads = NULL) {
  
  #Pulling out dimensions to test for matrix type
  row_test = dim(cont_mat1)[1]
  col_test = dim(cont_mat1)[2]
  
  if (row_test == col_test) {
    if (all(is.finite(cont_mat1)) == FALSE) {
      stop("Contact matrix 1 contains non-numeric entries")
    }
    
    if (all(is.finite(cont_mat2)) == FALSE) {
      stop("Contact matrix 2 contains non-numeric entries")
    }
  }
  
  if (col_test == 3) {
    
    
    #Convert sparse matrix to n x n matrix
    
    message("Converting to n x n matrix")
    
    if (nrow(cont_mat1) == 1) {
      stop("Matrix 1 is too small to convert to full")
    }
    
    if (nrow(cont_mat2) == 1) {
      stop("Matrix 2 is too small to convert to full")
    }
    
    cont_mat1 = sparse2full(cont_mat1,resolution)
    cont_mat2 = sparse2full(cont_mat2,resolution)
    
    if (all(is.finite(cont_mat1)) == FALSE) {
      stop("Contact matrix 1 contains non-numeric entries")
    }
    
    if (all(is.finite(cont_mat2)) == FALSE) {
      stop("Contact matrix 2 contains non-numeric entries")
    }
    if (resolution == "auto") {
      message("Estimating resolution")
      resolution = as.numeric(names(table(as.numeric(colnames(cont_mat1))-
                                            lag(
                                              as.numeric(
                                                colnames(cont_mat1)
                                              ))))[1]
      )
    }
    
  } else if (col_test-row_test == 3) {
    
    message("Converting to n x n matrix")
    
    #Find the start coordinates based on the second column of the
    #bed file portion of matrix
    
    start_coords = cont_mat1[,2]
    
    #Calculate resolution based on given bin size in bed file
    
    resolution = as.numeric(cont_mat1[1,3])-as.numeric(cont_mat1[1,2])
    
    #Remove bed file portion
    
    cont_mat1 = as.matrix(cont_mat1[,-c(seq_len(3))])
    cont_mat2 = as.matrix(cont_mat2[,-c(seq_len(3))])
    
    if (all(is.finite(cont_mat1)) == FALSE) {
      stop("Contact matrix 1 contains non-numeric entries")
    }
    
    if (all(is.finite(cont_mat2)) == FALSE) {
      stop("Contact matrix 2 contains non-numeric entries")
    }
    
    #Make column names correspond to bin start
    
    colnames(cont_mat1) = start_coords
    colnames(cont_mat2) = start_coords
    
    
  } else if (col_test!=3 & (row_test != col_test) & (col_test-row_test != 3)) {
    
    #Throw error if matrix does not correspond to known matrix type
    
    stop("Contact matrix must be sparse or n x n or n x (n+3)!")
    
  } else if ( (resolution == "auto") & (col_test-row_test == 0) ) {
    message("Estimating resolution")
    
    #Estimating resolution based on most common distance between loci
    
    resolution = as.numeric(names(table(as.numeric(colnames(cont_mat1))-
                                          lag(
                                            as.numeric(colnames(cont_mat1))
                                          )))[1])
  }
  
  #Make sure contact matrices only include shared columns
  coord_sum = list(colnames(cont_mat1), colnames(cont_mat2))
  
  #Only include shared columns in analysis
  shared_cols = Reduce(intersect, coord_sum)
  
  cont_mat1 = cont_mat1[colnames(cont_mat1) %in% shared_cols,
                        colnames(cont_mat1) %in% shared_cols]
  
  cont_mat2 = cont_mat2[colnames(cont_mat2) %in% shared_cols,
                        colnames(cont_mat2) %in% shared_cols]
  
  #Set maximize size of sliding window
  
  window_size = window_size
  
  #Remove full gaps from matrices
  
  non_gaps = which(colSums(cont_mat1) !=0 & (colSums(cont_mat2) !=0))
  
  #Remove gaps
  cont_mat1 = cont_mat1[non_gaps,non_gaps]
  cont_mat2 = cont_mat2[non_gaps,non_gaps]
  
  #Defining window size
  max_end = window_size
  max_size = window_size/ceiling(200000/resolution)
  min_size = ceiling(200000/resolution)
  Group_over = bind_rows()
  
  start = 1
  end = max_end
  
  end_loop = 0
  
  #If window is larger than matrix make it equal to matrix size
  if (end+window_size>nrow(cont_mat1)) {
    end = nrow(cont_mat1)
  }
  
  #Pre-allocate vectors
  point_dists1 = c()
  point_dists2 = c()
  Regions = c()
  
  while (end_loop == 0) {
    #Subsetting
    sub_filt1 = cont_mat1[seq(start,end,1), seq(start,end,1)]
    sub_filt2 = cont_mat2[seq(start,end,1), seq(start,end,1)]
    
    #Removing gap regions from sub_matrices
    
    Per_Zero1 = colSums(sub_filt1 !=0)/nrow(sub_filt1)
    Per_Zero2 = colSums(sub_filt2 !=0)/nrow(sub_filt2)
    
    #Remove columns with more zeros than threshold
    sub_gaps1 = Per_Zero1>gap_thresh
    sub_gaps2 = Per_Zero2>gap_thresh
    
    comp_rows = sub_gaps1 & sub_gaps2
    sub_filt1 = sub_filt1[ comp_rows, comp_rows]
    sub_filt2 = sub_filt2[ comp_rows, comp_rows]
    
    #Slide window to end if window size is less than 2
    if ( (length(sub_filt1) == 0) | (length(sub_filt1) == 1) ) {
      start = start+max_end
      end = end+max_end
      
    } else {
      
      #Getting degree matrices
      
      dr1 = rowSums(abs(sub_filt1))
      dr2 = rowSums(abs(sub_filt2))
      
      #Creating the normalized laplacian
      
      Dinvsqrt1 = diag((1/sqrt(dr1+2e-16)))
      Dinvsqrt2 = diag((1/sqrt(dr2+2e-16)))
      
      P_Part1 = crossprod(as.matrix(sub_filt1), Dinvsqrt1)
      sub_mat1 = crossprod(Dinvsqrt1, P_Part1)
      
      P_Part2 = crossprod(as.matrix(sub_filt2), Dinvsqrt2)
      sub_mat2 = crossprod(Dinvsqrt2, P_Part2)
      
      #Reading names
      
      colnames(sub_mat1) = colnames(sub_mat2) = colnames(sub_filt1)
      
      #Find gaps at 2mb and remove
      
      
      #Get first two eigenvectors
      
      Eigen1 = PRIMME::eigs_sym(sub_mat1, NEig = 2)
      
      eig_vals1 = Eigen1$values
      eig_vecs1 = Eigen1$vectors
      
      #Get order of eigenvalues from largest to smallest
      
      large_small1 = order(-eig_vals1)
      
      eig_vals1 = eig_vals1[large_small1]
      eig_vecs1 = eig_vecs1[,large_small1]
      
      #Repeat for matrix 2
      
      Eigen2 = eigs_sym(sub_mat2, NEig = 2)
      
      eig_vals2 = Eigen2$values
      eig_vecs2 = Eigen2$vectors
      
      #Get order of eigenvalues from largest to smallest
      
      large_small2 = order(-eig_vals2)
      
      eig_vals2 = eig_vals2[large_small2]
      eig_vecs2 = eig_vecs2[,large_small2]
      
      
      #Normalize the eigenvectors
      
      norm_ones = sqrt(dim(sub_mat1)[2])
      
      for (i in seq_len(dim(eig_vecs1)[2])) {
        eig_vecs1[,i] = (eig_vecs1[,i]/sqrt(sum(eig_vecs1[,i]^2)))  * norm_ones
        if (eig_vecs1[1,i] !=0) {
          eig_vecs1[,i] = -1*eig_vecs1[,i] * sign(eig_vecs1[1,i])
        }
      }
      
      
      for (i in seq_len(dim(eig_vecs2)[2])) {
        eig_vecs2[,i] = (eig_vecs2[,i]/sqrt(sum(eig_vecs2[,i]^2)))  * norm_ones
        if (eig_vecs2[1,i] !=0) {
          eig_vecs2[,i] = -1*eig_vecs2[,i] * sign(eig_vecs2[1,i])
        }
      }
      
      eps = 2.2204e-16
      
      n = dim(eig_vecs1)[1]
      k = dim(eig_vecs1)[2]
      
      #Project eigenvectors onto a unit circle
      
      vm1 = matrix(
        kronecker(rep(1,k), as.matrix(sqrt(rowSums(eig_vecs1^2)))),n,k
      )
      eig_vecs1 = eig_vecs1/vm1
      
      vm2 = matrix(
        kronecker(rep(1,k), as.matrix(sqrt(rowSums(eig_vecs2^2)))),n,k
      )
      eig_vecs2 = eig_vecs2/vm2
      
      #Get distance between points on circle
      
      point_dist1 = sqrt(
        rowSums( (eig_vecs1-rbind(NA,eig_vecs1[-nrow(eig_vecs1),]))^2)
      )
      point_dist2 = sqrt(
        rowSums( (eig_vecs2-rbind(NA,eig_vecs2[-nrow(eig_vecs2),]))^2)
      )
      
      #Remove NA entry at start of windows
      
      point_dists1 = c(point_dists1, point_dist1[-1])
      point_dists2 = c(point_dists2, point_dist2[-1])
      
      #Assign to regions based on column names
      Regions = c(Regions, colnames(sub_filt1)[-1])
      
    }
    
    #Test if we've reached end of matrix
    if (end == nrow(cont_mat1)) {
      end_loop = 1
    }
    
    #Set new start and end for window
    start = end
    end = end+max_end
    
    if ( (end + max_end) >nrow(cont_mat1)) {
      end = nrow(cont_mat1)
    }
    
    if (start == end | start>nrow(cont_mat1)) {
      end_loop = 1
    }
  }
  #Calculating the difference between log gaps
  
  dist_diff = (point_dists1)-(point_dists2)
  #Getting the z-scores
  sd_diff = (dist_diff-mean(dist_diff, na.rm = TRUE))/(sd(dist_diff,
                                                          na.rm = TRUE))
  
  TAD_Score1 = (point_dists1-mean(point_dists1, na.rm = TRUE))/
    (sd(point_dists1, na.rm = TRUE))
  
  TAD_Score2 = (point_dists2-mean(point_dists2, na.rm = TRUE))/
    (sd(point_dists2, na.rm = TRUE))
  
  #Get areas with high z-scores
  gaps = which(abs(sd_diff)>z_thresh)
  
  #Put differential regions into a data frame
  diff_loci = data.frame(Region = as.numeric(Regions)[gaps],
                         Gap_Score = sd_diff[gaps])
  
  #Return differential TAD boundaries
  Gap_Scores = data.frame(Boundary = as.numeric(Regions),
                          TAD_Score1 = TAD_Score1,
                          TAD_Score2 =TAD_Score2,
                          Gap_Score = sd_diff)
  TAD_Frame = data.frame(Boundary = as.numeric(Regions),
                         Gap_Score = sd_diff,
                         TAD_Score1,
                         TAD_Score2)
  
  #Assign labels to boundary type and identify which matrix has the boundary
  
  if(!is.null(pre_tads)) {
    
    pre_tads = lapply(pre_tads, as.data.frame)
    
    #pre_tads = bind_rows(pre_tads)
    TAD_Frame = TAD_Frame %>% 
      filter(Boundary %in% bind_rows(pre_tads)$end)  %>%
      mutate(Differential = ifelse(abs(Gap_Score)>z_thresh, "Differential",
                                   "Non-Differential"),
             Enriched_In = ifelse(Gap_Score>0, "Matrix 1", "Matrix 2")) %>%
      arrange(Boundary) %>%
      mutate(Bound_Dist =  abs(Boundary-lag(Boundary))/resolution) %>%
      mutate(Differential = ifelse( (Differential == "Differential") &
                                      (Bound_Dist<=5) & !is.na(Bound_Dist) & 
                                      (   Enriched_In!=lag(Enriched_In)) &
                                      (lag(Differential)=="Differential"),
                                    "Shifted", Differential)) %>% 
      mutate(Differential= ifelse(lead(Differential) == "Shifted", "Shifted",
                                  Differential)) %>%
      dplyr::select(-Bound_Dist)
    
    #Pull out non-shared boundaries
    
    
  } else {
    
    TAD_Frame = TAD_Frame %>%
      filter( (TAD_Score1>1.5) | TAD_Score2>1.5) %>%
      mutate(Differential = ifelse(abs(Gap_Score)>z_thresh, "Differential",
                                   "Non-Differential"),
             Enriched_In = ifelse(Gap_Score>0, "Matrix 1", "Matrix 2")) %>%
      arrange(Boundary) %>%
      mutate(Bound_Dist =  abs(Boundary-lag(Boundary))/resolution) %>%
      mutate(Differential = ifelse( (Differential == "Differential") &
                                      (Bound_Dist<=5) & !is.na(Bound_Dist) & 
                                      (   Enriched_In!=lag(Enriched_In)) &
                                      (lag(Differential)=="Differential"),
                                    "Shifted", Differential)) %>% 
      mutate(Differential= ifelse(lead(Differential) == "Shifted", "Shifted",
                                  Differential)) %>%
      dplyr::select(-Bound_Dist)
    
  }
  #Classifying merged-split
  TAD_Frame = TAD_Frame %>%
    mutate(Type = ifelse( (Differential == "Differential") &
                            (lag(Differential) == "Non-Differential") &
                            (lead(Differential) == "Non-Differential"),
                          ifelse(Enriched_In == "Matrix 1", "Split", "Merge"),
                          Differential))
  
  #Add up-down enrichment of TAD boundaries
  TAD_Frame = TAD_Frame %>%
    mutate(Type = ifelse( (TAD_Score1>1.5) &
                            (TAD_Score2>1.5) &
                            (Differential == "Differential"),
                          "Strength Change", Type))
  
  #Classify leftovers as complex
  
  TAD_Frame = TAD_Frame %>% mutate(Type = gsub("^Differential$",
                                               "Complex", Type))
  
  #Another step for pre-specified
  
  if (!is.null(pre_tads)) {
    #Pulling out shared ends by overlap
    shared_ends = ((TAD_Frame$Boundary  %in%
                      pre_tads[[1]]$end + TAD_Frame$Boundary  %in%
                      pre_tads[[2]]$end)==1)
    #Converting non-differential to non-overlap
    TAD_Frame = TAD_Frame %>% mutate(Type = ifelse(
      (shared_ends == TRUE)&(Type=="Non-Differential"),
      "Non-Overlap", Type))
  }
  #Redo for gap score frame as well
  
  #Assign labels to boundary type and identify which matrix has the boundary
  
  Gap_Scores = Gap_Scores %>%
    mutate(Differential = ifelse(abs(Gap_Score)>z_thresh, "Differential",
                                 "Non-Differential"),
           Enriched_In = ifelse(Gap_Score>0, "Matrix 1", "Matrix 2")) %>%
    arrange(Boundary) %>%
    mutate(Bound_Dist = pmin(abs(Boundary-lag(Boundary))/resolution,
                             abs((Boundary-lead(Boundary)))/resolution)) %>%
    
    mutate(Differential = ifelse( (Differential == "Differential") &
                                    (Bound_Dist<=5) & !is.na(Bound_Dist),
                                  "Shifted", Differential)) %>%
    dplyr::select(-Bound_Dist)
  
  #Classifying merged-split
  Gap_Scores = Gap_Scores %>%
    mutate(Type = ifelse( (Differential == "Differential") &
                            (lag(Differential) == "Non-Differential") &
                            (lead(Differential) == "Non-Differential"),
                          ifelse(Enriched_In == "Matrix 1", "Split", "Merge"),
                          Differential))
  
  #Add up-down enrichment of TAD boundaries
  Gap_Scores = Gap_Scores %>%
    mutate(Type = ifelse( (TAD_Score1>1.5) &
                            (TAD_Score2>1.5) &
                            (Differential == "Differential"),
                          "Strength Change", Type))
  
  #Classify leftovers as complex
  
  Gap_Scores = Gap_Scores %>% mutate(Type = gsub("^Differential$",
                                                 "Complex", Type))
  
  
  TAD_Sum = TAD_Frame %>% group_by(Type) %>% summarise(Count = n())
  
  #Fix double counting of shifted boundaries
  
  TAD_Sum = TAD_Sum %>% mutate(Count = ifelse(Type == "Shifted",
                                              Count/2,
                                              Count))
  
  Count_Plot = ggplot2::ggplot(TAD_Sum,
                               aes(x = 1,
                                   y = Count, fill = Type)) +
    geom_bar(stat="identity") + theme_bw(base_size = 24) +
    theme(axis.title.x = element_blank(), panel.grid = element_blank(),
          axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    labs(y = "Number of Boundaries")
  
  return(list(TAD_Frame =TAD_Frame,
              Boundary_Scores = Gap_Scores,
              Count_Plot = Count_Plot ))
}

sys_argv <- commandArgs(T)
input_dir <- sys_argv[1]
tags = sys_argv[2]
resolution = as.numeric(sys_argv[3])

tags = unlist(strsplit(tags, ","))

###dis_bdr##########
tmp_dir = file.path(input_dir,'insulation_boundary_di')
img_dir = file.path(input_dir,'img')
insulation_distance = list()
di_distance = list()

for (tag in tags){
  insulation_distance[[tag]] = read.table(file.path(tmp_dir,paste0(tag,'.insulation_distance2boundary.txt')),header = FALSE)
  insulation_distance[[tag]] = data.frame(insulation_distance[[tag]],sample=tag)
  colnames(insulation_distance[[tag]])[1:2] <- c("dis","insulation")
  di_distance[[tag]] = read.table(file.path(tmp_dir,paste0(tag,'.direction_distance2boundary.txt')),header = FALSE)
  di_distance[[tag]] = data.frame(di_distance[[tag]],sample=tag)
  colnames(di_distance[[tag]])[1:2] <- c("dis","di")
}

insulation_all <- do.call(rbind, insulation_distance)
di_all <- do.call(rbind, di_distance)

pdf(file.path(img_dir,'insulation_distance2boundary.pdf'),width = 10,height = 10)
ggplot(insulation_all)+geom_smooth(aes(x=dis,y=insulation,col=sample),se = FALSE)+
  theme_linedraw(base_line_size = I(0))+
  theme(text = element_text(size = 25),axis.text.x = element_text(angle = 45,hjust = 1),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  xlab("distance to TAD boundaries / bp")
while (!is.null(dev.list()))  dev.off()

pdf(file.path(img_dir,'direction_distance2boundary.pdf'),width = 10,height = 10)
ggplot(di_all)+geom_smooth(aes(x=dis,y=di,col=sample),se = FALSE)+
  theme_linedraw(base_line_size = I(0))+
  theme(text = element_text(size = 25),axis.text.x = element_text(angle = 45,hjust = 1),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  xlab("distance to TAD boundaries / bp")
while (!is.null(dev.list()))  dev.off()

# TAD compare
contact_dir = file.path(input_dir,'contact')
tad_dir = file.path(input_dir,'tad')

# Get common chroms
chr_lists <- list()
for (tag in tags) {
  pattern <- paste0("^", tag, "_chr[0-9XYM]{1,2}_contact\\.txt$") 
  files <- list.files(contact_dir, pattern = pattern, full.names = FALSE)
  chr_names <- str_extract(files, "chr[0-9XYM]{1,2}")
  chr_lists[[tag]] <- unique(na.omit(chr_names))
}
chroms = Reduce(intersect, chr_lists)

chrom_cmp = list()
if (length(tags) > 1){
  for (i in 1:(length(tags)-1)){
    for (j in (i+1):length(tags)){
      vs_files = paste0(tags[i],'_vs_',tags[j])
      for (chrom in chroms) {
        cat(paste0('Process ',chrom,'\n'))
        tryCatch({
          file1 = read.table(file.path(contact_dir, paste0(tags[i], '_', chrom, "_contact.txt")))
          file2 = read.table(file.path(contact_dir, paste0(tags[j], '_', chrom, "_contact.txt")))
          results <- TADCompare(file1, file2, resolution = resolution)
          temp = results$TAD_Frame
          temp$Type[(temp$Type == 'Strength Change') & (temp$TAD_Score2 > temp$TAD_Score1)] = 'strength-up'
          temp$Type[(temp$Type == 'Strength Change') & (temp$TAD_Score2 < temp$TAD_Score1)] = 'strength-down'
          chrom_cmp[[chrom]] = temp
        }, error = function(e) {
          message(paste("Error processing chromosome", chrom, ":", e$message))
          chrom_cmp[[chrom]] <- NULL  
        })
      }
      cmp = do.call(rbind,chrom_cmp)
      output_file = file.path(tad_dir,paste0(vs_files,'.tadcompare.csv'))
      write.csv(cmp, output_file, quote = FALSE, sep = ",", row.names = FALSE, col.names = TRUE)
    }
  }
}

