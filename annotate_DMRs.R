

annotate_DMRs <- function(dm_file){
#chr     start       end length nCG meanMethy1 meanMethy2 diff.Methy  areaStat
extraCols = c(length = 'integer', nCG = 'integer', meanMethy1 = 'numeric',meanMethy2='numeric',diff.Methy="numeric",areaStat="numeric")

}

 # chr      start        end     length        nCG meanMethy1 meanMethy2 diff.Methy 
  #"factor"  "integer"  "integer"  "integer"  "integer"  "numeric"  "numeric"  "numeric" 
  #areaStat 
 #"numeric" 

dm_regions = read_regions(con = dm_file, genome = 'hg19', extraCols = extraCols, format = 'bed',
