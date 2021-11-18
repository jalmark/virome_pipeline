# 16.2.2021 
# Jalmari Kettunen
# Kehitysvaiheessa oleva R-skripti
# Mukaan taksonomia Artun NCBI-taulukoista
# Kaytetty analysoimaan naytteet gtex_vanhimmat_B 21.1.2021
##### Aseta nimet oikein kohdissa 1.1, 1.2 ja 1.3! ######
# Read count tiedostojen taytyy olla samassa hakemistossa keskenaan, ei muita tiedostoja!


library(stringr)
library(openxlsx)


# Lue Excel-taulukot taksonomiatiedoiksi
readRefseqs = function() {
  # 1.1 Kirjoita oikea polku tahan!
  setwd("~/postGradu/Artun_skripteja")
  
  print("Reading 1/4")
  virInfo1 = read.xlsx("virusInfoRefseq1.xlsx", colNames = FALSE)
  print("Reading 2/4")
  virInfo2 = read.xlsx("virusInfoRefseq2.xlsx", colNames = FALSE)
  print("Reading 3/4")
  virInfo3 = read.xlsx("virusInfoRefseq3.xlsx", colNames = FALSE)
  print("Reading 4/4")
  virInfo4 = read.xlsx("virusInfoRefseq4.xlsx", colNames = FALSE)
  virInfo = rbind(virInfo1,virInfo2,virInfo3,virInfo4)
  colnames(virInfo) = c("Accession","Release_Date","Species","Genus","Family","Length",
                        "Nuc._Completeness","Genotype","Segment","Country",
                        "Host","Isolation_Source","Collection_Date","BioSample","GenBank_Title")
  # Poista muistin saastamiseksi
  rm(virInfo1)
  rm(virInfo2)
  rm(virInfo3)
  rm(virInfo4)
  
  return(virInfo)
}






# Yhdist‰ read countit eri tiedostoista samaan DataFrameen ja yhdista taksonomiatietoihin
readSampleReads = function() {
  # 1.2 Kirjoita oikea polku tahan!
  setwd("~/postGradu/GTEx/GTEx_vanhimmat_testi")
  
  tulostaulukko = data.frame()
  
  for (k in seq(1,length(dir()))) {
    uusi_nayte = read.delim(dir()[k], header = FALSE)
    
    # Ota vain olennaiset sarakkeet: ID ja read count.
    uusi_nayte = uusi_nayte[,c(1,3)]
    
    # Muokkaa Virosaurus ID vastaamaan NCBI GenBankin accessionia ID:ta
    uusi_nayte[,1] = str_replace(uusi_nayte[,1], '\\:\\S+\\;',"")
    
    # Summaa sitten saman ID:n read countit
    uusi_nayte = aggregate(V3 ~ V1, uusi_nayte, sum)
    
    # Jos eka tiedosto, joka lisataan DataFrameen...
    if (k == 1) {
      colnames(uusi_nayte) = c('accession', dir()[k])
      tulostaulukko = uusi_nayte
    }
    # ...muulloin...
    else {
      tulostaulukko[[dir()[k]]] = rep(0,nrow(tulostaulukko))                        # Luodaan uusi sarake taynna nollia.
      for (i in 1:nrow(uusi_nayte)) {
        if (uusi_nayte[i,1] %in% tulostaulukko[,1]) {                               # Jos accession on jo tulostaulukossa, lis‰‰ readCount oikealle riville.
          correctRow = which(as.vector(tulostaulukko[,1]) == as.vector(uusi_nayte[i,1]))
          tulostaulukko[correctRow,ncol(tulostaulukko)] = uusi_nayte[i,2]
        }
        else {                                                                      # Jos ei, lis‰‰ uusi rivi.
          newRow = cbind(as.vector(uusi_nayte[i,1]), t(rep(0,ncol(tulostaulukko)-2)), uusi_nayte[i,2])
          colnames(newRow) = colnames(tulostaulukko)
          tulostaulukko = rbind(tulostaulukko,newRow)
        }
      }
    }
  }
  # Muuta read count sarakkeet numeerisiksi
  tulostaulukko[ , c(2:ncol(tulostaulukko))] <- apply(tulostaulukko[ , c(2:ncol(tulostaulukko))], 2, function(x) as.numeric(x))
  
  
  
  
  
  
  
  # Hylataan alle viiden readin kentat "kohinana"
  maxvektori = apply(tulostaulukko[,2:ncol(tulostaulukko)], 1, FUN = max)
  cutoff5 = tulostaulukko[which(maxvektori > 4),]
  for (i in 2:ncol(cutoff5)) {
    for (j in 1:nrow(cutoff5)) {
      if (cutoff5[j,i] < 5) {
        cutoff5[j,i] = 0
      }
    }
  }
  
  
  
  
  # Yhdista GenBankin taksonomiatiedot accession ID:n perusteella
  # Valmistele
  species = rep("NA", nrow(cutoff5))
  genus = rep("NA", nrow(cutoff5))
  family = rep("NA", nrow(cutoff5))
  genotype = rep("NA", nrow(cutoff5))
  segment = rep("NA", nrow(cutoff5))
  country = rep("NA", nrow(cutoff5))
  host = rep("NA", nrow(cutoff5))
  isol_source = rep("NA", nrow(cutoff5))
  biosample = rep("NA", nrow(cutoff5))
  genbank = rep("NA", nrow(cutoff5))
  
  
  
  
  # Hae tiedot jokaiselle DataFramen ID:lle
  for (i in 1:nrow(cutoff5)) {
    correctRow = which(virInfo["Accession"] == cutoff5[i,"accession"])
    if (length(correctRow) > 0) {
      species[i] = virInfo[correctRow, "Species"]
      genus[i] = virInfo[correctRow, "Genus"]
      family[i] = virInfo[correctRow, "Family"]
      genotype[i] = virInfo[correctRow, "Genotype"]
      segment[i] = virInfo[correctRow, "Segment"]
      country[i] = virInfo[correctRow, "Country"]
      host[i] = virInfo[correctRow, "Host"]
      isol_source[i] = virInfo[correctRow, "Isolation_Source"]
      biosample[i] = virInfo[correctRow, "BioSample"]
      genbank[i] = virInfo[correctRow, "GenBank_Title"]
    }
  }
  
  # Yhdista tiedot DataFrameen
  DFWithTaxonomy = cbind(cutoff5, species, genus, family, genotype, segment, country, host, isol_source, biosample, genbank)
  
  # Siirr‰ accession ID:t rivinimiksi
  row.names(DFWithTaxonomy) = DFWithTaxonomy[,1]
  DFWithTaxonomy = DFWithTaxonomy[,2:ncol(DFWithTaxonomy)]
  
  return(DFWithTaxonomy)
}



virInfo = readRefseqs()
result = readSampleReads()

# 1.3 Vaihda sarakenimet GTEx:n SUBJID:ksi.
# Varmista, etta ovat oikeassa jarkassa!
colnames[1:15] = c()

# Vie DataFrame tekstitiedostoon ja sit‰ kautta Exceliin
write.table(result, file="noname.txt", sep='\t', quote=FALSE)
# Poista virInfo, jos haluat saastaa tallennustilaa.
rm(virInfo)
# Tallenna workspace haluamallasi nimella.
save.image("noname.RData")







