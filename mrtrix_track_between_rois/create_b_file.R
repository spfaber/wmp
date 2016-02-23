## Sam Faber
## 10-23-2015
## convert bval/bvec to combined mrtrix .b file
##

## read in bval file and convert it to an unlabeled vector
bval <- read.table("dwi_data_b2000_aligned_trilin.bvals", sep = " ")
colnames(bval) <- NULL
bval <- unlist(bval)

## read in bvecs and convert them to a data frame
bvec <- read.table("dwi_data_b2000_aligned_trilin.bvecs", sep = " ")
bvec <- t(bvec)
bvec <- as.data.frame(bvec)

## merge bval/bvec into .b output format by concatenation
bvec$bval <- bval

## save text file out in .b format
write.table(bvec, "dwi_data_b2000_aligned_trilin.b", quote = FALSE, col.names = FALSE, row.names = FALSE)
