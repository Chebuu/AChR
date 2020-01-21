smi2ap <- function(smi) {
  return(sdf2ap(smiles2sdf(smi)))
}
getDescriptors <- function(x.mol) {
  return(cbind(extractDrugALOGP(x.mol),
        extractDrugApol(x.mol),
        extractDrugECI(x.mol),
        extractDrugTPSA(x.mol),
        extractDrugWeight(x.mol),
        extractDrugWienerNumbers(x.mol),
        extractDrugZagrebIndex(x.mol)))
}
getDescriptorsByBatch <- function(smis) { # for a list of groups, each group (mAg, nAn, nAg, mAn), return a data frame of descriptors for each molecule within the group
  return(sapply(smis, function(sets){
    sapply(sets, function(set){
      sapply(sets, function(x.mol){
        getDescriptors(x.mol)
      })
    })
  }))
}
removeNADescriptors <- function(descriptors) {
  apply(descriptors, 2, function(set){
    sapply(set, function(row){
      sapply(row, function(cell) {
        ifelse(is.na(cell), 0, cell)
      })
    })
  })
}


