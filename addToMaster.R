# Jalmari Kettunen
# Last modified 18.11.2021
# This R script is designed for processing of virome pipeline results.
# It takes output of Jalmarin_viromiDev.R as input.
# Includes 3 functions.










# GTEx:n metadatassa vasta-ainetestin arvo voi olla 96 (indeterminate) tai 97 (not performed)
# Halutaan, etta nama arvot lisataan masterfileen NA-arvoina.
# Tama funktio hoitaa homman.
checkTestResult = function(arvo) {
  if (arvo == 97 || arvo == 96) {
    return(NA)
  }
  else {
    return(arvo)
  }
}



# Normalisointifunktio
# Normalisoi per miljoonaa laadukasta readiparia.
# Funktio addToMaster() on riippuvainen tasta.
Normalisoi = function(taulu, laatureadiVektori) {
  NormFactors = (10^6)/laatureadiVektori
  # Tee normalisaatio vain read count sarakkeille.
  abundance = apply(taulu[,9:ncol(taulu)], 2, function(x){t = x*NormFactors; t})
  # Yhdistä sitten takaisin muuhun masterfileen
  taulu = cbind(taulu[,1:8],abundance)
  return(taulu)
}









# Funktio addToMaster on riippuvainen kirjastosta stringr, funktioista Normalisoi() ja checkTestResult() seka dataframesta phenotypesOf317.
# Funktiolla on 4 syotetta:
# 1) Jo olemassa oleva masterfile (dataframe) (alusta nollaksi ennen funktion kayttoa, jos ei viela olemassa).
# 2) batch on satsi (dataframe), joka halutaan lisata masterfileen.
# HUOM! Näytesarakkeiden (alkavat sanalla "GTEX") täytyy olla ensimmäisinä sarakkeina peräkkäin batch-objektissa!
# 3) qualityReadPairsPerSample on normalisointia varten: numeerinen vektori, jossa laadukkaiden readiparien maarat.
# HUOM! SUBJID:t taytyy olla samassa jarkassa (aakkosjarkassa) vektorissa ja satsin sarakenimissa. 
# 4) date on paivamaara, jolloin satsi on analysoitu (merkkijono, esim. '16.12.2020')

addToMaster = function(mastertaulu, batch, qualityReadPairsPerSample, date) {
  library(stringr)
  # sampleSize = kuinka monta naytetta satsissa on.
  sampleSize = length(which(startsWith(colnames(batch),"GTEX")))
  # Sisaanluvussa naytteiden nimet vaaristyvat. Korjaa. Jarjesta ne  myos aakkosjarjestykseen.
  colnames(batch)[1:sampleSize] = str_replace(colnames(batch)[1:sampleSize], "\\.", "\\-")
  batch = cbind( batch[,order(colnames(batch)[1:sampleSize])] , batch[,(sampleSize+1):ncol(batch)] )
  
  pvm = rep(date, sampleSize)
  # Kasitellaan jatkossa vain batchin transpoosia, ei batchia itseaan.
  temp = t(batch[,1:sampleSize])
  
  # Etsi fenotyyppitiedot
  group = c()
  age = c()
  sex = c()
  CMV_ab = c()
  EBV_IgG = c()
  EBV_IgM = c()
  
  for (i in 1:sampleSize) {
    correctRowNumber = which(phenotypesOf317$SUBJID == rownames(temp)[i])
    age = c(age, phenotypesOf317[correctRowNumber,]$AGE)
    sex = c(sex, phenotypesOf317[correctRowNumber,]$SEX)
    CMV_ab = c(CMV_ab, checkTestResult(phenotypesOf317[correctRowNumber,]$LBCMVTAB))
    EBV_IgG = c(EBV_IgG, checkTestResult(phenotypesOf317[correctRowNumber,]$LBEBVGAB))
    EBV_IgM = c(EBV_IgM, checkTestResult(phenotypesOf317[correctRowNumber,]$LBEBVMAB))
    if (phenotypesOf317[correctRowNumber,]$AGE > 59) {
      group = c(group, 1)
    } else {
      group = c(group, 2)
    }
  }
  
  
  
  
  
  # Rakenna read countit
  # Jos mastertaulu ei ole tyhja...
  if (length(dim(mastertaulu)) > 0) {
    # tempDf on aluksi nollataulu, joka vastaa mastertaulun read count sarakkeita, ei muita sarakkeita.
    tempDf = as.data.frame(matrix(0, nrow=sampleSize, ncol = (ncol(mastertaulu)-8) ),
                           row.names = rownames(temp))
    colnames(tempDf) = colnames(mastertaulu)[9:ncol(mastertaulu)]
    
    
    for (j in 1:ncol(temp)) {
      # Jos satsin accession on jo mastertaulussa, lisaa satsin sarake dataFrameen tempDf oikeaan sarakkeeseen.
      if (colnames(temp)[j] %in% colnames(mastertaulu)){
        tempDf[[ colnames(temp)[j] ]] = temp[,j]
      } else {
        # Jos ei, mastertauluun lisataan nollasarakkeita, jotta yhdistaminen onnistuu myohemmin.
        mastertaulu[[colnames(temp)[j]]] = rep(0, nrow(mastertaulu))
        tempDf[[ colnames(temp)[j] ]] = temp[,j]
      }
    }
    # Laajenna tempDf fenotyyppitiedoilla.
    tempDf = cbind(SUBJID = rownames(temp), Date=date, Group=group, age=age, sex=sex, CMV_ab=CMV_ab, EBV_IgG=EBV_IgG, EBV_IgM=EBV_IgM, 
                   tempDf)
    # Normalisoi
    tempDf = Normalisoi(tempDf,qualityReadPairsPerSample)
    
    # Yhdista tempDf ja mastertaulu.
    mastertaulu = rbind(mastertaulu, tempDf)
    
    
  } else {
    # Jos mastertaulu on tyhja, mastertaulu = eka lisattava satsi
    mastertaulu = data.frame(SUBJID = rownames(temp), Date=pvm, Group=group, age=age, sex=sex, CMV_ab=CMV_ab, EBV_IgG=EBV_IgG, EBV_IgM=EBV_IgM,
                             row.names = rownames(temp))
    mastertaulu = cbind(mastertaulu, temp)
    # Normalisoi
    mastertaulu = Normalisoi(mastertaulu,qualityReadPairsPerSample)
  }
  
  
  
  # Poista ilmeiset vaarat positiiviset virusreferenssit
  # Ks. tiedosto GTEx_masterfile_selite.txt (28.5.2021)
  poistettavat = c('KU746283', 'KU746280', 'AB972431', 'KF478765', 'JQ679013', 
                   'MN646692', 'MN646693', 'MN646694', 'MN646695', 'MF444119',
                   'S57428')
  mastertaulu = mastertaulu[,!colnames(mastertaulu) %in% poistettavat]
 
}