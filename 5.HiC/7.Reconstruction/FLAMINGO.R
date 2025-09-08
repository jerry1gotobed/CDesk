library(FLAMINGOrLite)

args = commandArgs(trailingOnly=TRUE)

hic_file <- args[1]
chr_num <- args[2]
domain_res = as.numeric(args[3])
frag_res = as.numeric(args[4])
normalization = args[5]
nThread = as.numeric(args[6])
output_dir = args[7]

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
setwd(output_dir)

available_norm = strawr::readHicNormTypes(hic_file)
if (!normalization %in% available_norm){
  cat('Available normalizations:')
  cat(available_norm)
  stop('Please specify the available normalization')
}
available_chr = strawr::readHicChroms(hic_file)$name
if (!chr_num %in% available_chr){
  cat('Available chromosomes:')
  cat(available_chr)
  stop('Please specify the available chromosome')
}
available_frag = strawr::readHicBpResolutions(hic_file)
if (!frag_res %in% available_frag){
  cat('Available fragment resolutions:')
  cat(available_frag)
  stop('Please specify the available fragment resolution')
}
available_domain = strawr::readHicBpResolutions(hic_file)
if (!domain_res %in% available_domain){
  cat('Available domain resolutions:')
  cat(available_domain)
  stop('Please specify the available domain resolution')
}

new_chr_num <- sub("chr","chchrr",chr_num)
## Note: chr name function has a sub "chr" to "", so using "chchrr" to avoid error.
res_wt <- flamingo_main(hic_data=hic_file,file_format='hic',domain_res=domain_res,frag_res=frag_res,chr_name=new_chr_num,normalization=normalization,nThread=nThread)

## output data to .vtk file
write.vtk(res_wt,lookup_table = res_wt$start,name = "out.vtk",opt_path = file.path(output_dir,"out.vtk"))

cat('Done, you can see the result now\n')
