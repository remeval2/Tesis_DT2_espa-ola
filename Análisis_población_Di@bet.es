

### Preparación de los datos clínicos

```{r,eval=doNot}
clinica_i = read_excel("../data/Di@betes_uni_20211202_en-uso.xlsx")
clinica_apoyo = read_sav("../data/2018 06 25 DI@BET.ES cols_ORI.sav")
clinica_apoyo_2 = read.spss("../data/2018 06 25 DI@BET.ES cols_ORI.sav", to.data.frame = TRUE)

#Variable respuesta: diabetes tipo 2

 clinica_i$DM[clinica_i$DM == 0] <- NA
 clinica_i$DM[clinica_i$DM == 1] <- 0
 clinica_i$DM[clinica_i$DM == 2] <- 1
 
 #Definimos los grupos de variables, las variables "_s" indican los datos del seguimiento.
 
 #Lipidos

 clinica_i$`HDLc_mg/dL` = as.numeric(clinica_i_3$`HDLc_mg/dL`)
 clinica_i$`HDLc_mg/dL_sgto` = as.numeric(clinica_i_3$`HDLc_mg/dL_sgto`)
 
 clinica_i$`LDLc_mg/dL` = as.numeric(clinica_i_3$`LDLc_mg/dL`)
 clinica_i$`LDLc_mg/dL_sgto`= as.numeric(clinica_i_3$`LDLc_mg/dL_sgto`)
 
  clinica_i$`CT_mg/dL` = as.numeric(clinica_i_3$`CT_mg/dL`)
  clinica_i$`CT_mg/dL_sgto` = as.numeric(clinica_i_3$`CT_mg/dL_sgto`)
 
 clinica_i$`TG_mg/dL` = as.numeric(clinica_i_3$`TG_mg/dL`)
 clinica_i$`TG_mg/dL_sgto` = as.numeric(clinica_i_3$`TG_mg/dL_sgto`)
 
 ITI_lipidos = c("HDLc_mg/dL","LDLc_mg/dL","CT_mg/dL","TG_mg/dL")
 ITI_lipidos_s = c("HDLc_mg/dL_sgto","LDLc_mg/dL_sgto","TG_mg/dL_sgto","CT_mg/dL_sgto")
 
 # Variables de confusión
 
 BMI  = as.numeric(clinica_i_3$`IMC_Kg/m2`)
 BMI_s = as.numeric(clinica_i_3$`IMC_Kg/m2_sgto`)

 age = as.numeric(clinica_i$Edad)
 age_s = as.numeric(clinica_i$edad_sgto)

 Sex = as.factor(clinica_i$Genero)
 
 #Eliminamos los datos de pacientes que están tratados farmacológicamente
 
 for (i in 1:nrow(clinica_i)){
 if(clinica_i$Tratamiento_farmacologico_hipolipemiantes[i] == 2){
  clinica_i$`HDLc_mg/dL`[i] <- NA
  clinica_i$`LDLc_mg/dL`[i] <- NA
  clinica_i$`TG_mg/dL`[i] <- NA
  clinica_i$`CT_mg/dL`[i] <- NA
}
  if(clinica_i$Tratamiento_farmacologico_ADO[i] == 2 || clinica_i$Tratamiento_farmacologico_DM_ADO_paciente[i] == 2 || clinica_i$Tratamiento_farmacologico_INS[i] == 2) {
    clinica_i$`Glucosa_basal_mg/dL`[i] <- NA
    clinica_i$`Insulina_basal_microU/mL`[i] <- NA
  }
  
  if(clinica_i$Tratamiento_farmacologico_HTA[i] == 2){
    clinica_i$PAD_mmHg_media[i] <- NA
  }
}

for (i in 1:nrow(clinica_i)){
   if (clinica_i$Tratamiento_farmacologico_hipolipemiante_paciente_sgto[i] == 2 ){
  clinica_i$`HDLc_mg/dL_sgto`[i] <- NA
  clinica_i$`LDLc_mg/dL_sgto`[i] <- NA
  clinica_i$`TG_mg/dL_sgto`[i] <- NA
  clinica_i$`CT_mg/dL_sgto`[i] <- NA
 }
}

for (i in 1:nrow(clinica_i)){
  if(clinica_i$Tratamiento_farmacologico_DM_ADO_paciente_sgto[3] == "2" || clinica_i$Tratamiento_farmacologico_DM_ins_paciente_sgto[i] == "2" ||        clinica_i$Tratamiento_farmacologico_DM_ins_paciente_sgto[i] == "2") {
    clinica_i$`Glucosa_basal_mg/dL_sgto`[i] <- NA
    clinica_i$`Insulina_basal_mU/mL_sgto`[i] <- NA
  }
}
 for (i in 1:nrow(clinica_i)){ 
   if(clinica_i$Tratamiento_farmacologico_HTA_paciente_sgto[i] == "2"){
    clinica_i$PAS_mmHg_media_sgto[i] <- NA
  }
}

```

### Preparación de los datos genómicos

#### Lectura de archivos MAP y PED que contienen la información genética

```{r carga_datos}

map <- read.table("~/FISDM3_ITI_DI@BET.ES_V38_merged_output_only_variants_20230929.map")
ped <- read.table("~/FISDM3_ITI_DI@BET.ES_V38_merged_output_only_variants_20230929.ped")
ID <- ped$V1
ped_m <- ped[, -1:-6]  # Eliminamos columnas de información no genética
```

#### Convertimos los pares de alelos a formato tipo "./."

```{r formateo_genotipos}
dataframe_m <- data.frame(rep(0, nrow(ped_m)))

for (i in 2:ncol(ped_m)) {
  if (i %% 2 == 0) {
    n <- i - 1
    x <- paste(ped_m[, n], ped_m[, i], sep = "/")
    dataframe_m <- as.data.frame(cbind(dataframe_m, x))
  }
}
dataframe_m <- dataframe_m[, -1]
colnames(dataframe_m) <- data.frame(t(map$V2))
```

#### Anotación de los SNPs

```{r anotacion_biomart}
# Extraemos los identificadores numéricos para usar en biomart
num <- unlist(str_extract_all(colnames(dataframe_m), "\\d+"))
chr <- num[2 * (1:length(num)) - 1]  
pos <- num[2 * (1:length(num))] 
rep_pos <- pos
bio <- data.frame(chr, pos, rep_pos)

# Conectamos con Ensembl para anotación de SNPs
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                dataset = "hsapiens_gene_ensembl",
                host = 'https://www.ensembl.org')

# Inicializamos dataframe vacío para guardar anotaciones
df_anotacion <- data.frame()

for (i in 1:nrow(bio)) {
  total_anotacion <- getBM(attributes = c("refsnp_id","ensembl_gene_name","chr_name",
                                          "chrom_start","allele", "minor_allele_freq",
                                          "consequence_type_tv", "clinical_significance"), 
                            filters = c('chromosome_name','start','end'),
                            values = list(bio[i,1],bio[i,2],bio[i,3]),
                            mart = mart)

  if (nrow(total_anotacion) > 0) {
    total_anotacion <- cbind(bio[i, ], total_anotacion)
    df_anotacion <- rbind(df_anotacion, total_anotacion)
  }
}
```
## Función de corrección por BH
```{r}
rs_correction = function(x,dim,method="BH"){
  res = apply(x,dim,p.adjust,method=method)
  if(dim == 1) res = t(res)
res
}
```

### Análisis de la asociación de las SNPs con la diabetes tipo 2 y selección de modelo

```{r seleccion_snp_2_cat}
# Selección de variantes con solo dos categorías genotípicas
Variants <- dataframe_m[1:4754, ]

# SNPs con menos de 3 niveles excluyendo "0/0"
d_2 <- which(sapply(Variants, function(y) {
  y <- factor(y)
  nlevels(y) - is.element("0/0", levels(y)) * 1
}) < 3)
Variants_2 <- Variants[, d_2]

# SNPs con al menos 2 niveles válidos
d_1 <- which(sapply(Variants_2, function(y) {
  y <- factor(y)
  nlevels(y) - is.element("0/0", levels(y)) * 1
}) >= 2)
Variants_2 <- Variants_2[, d_1]
```

```{r analisis_3_cat}
# Selección de variantes con tres categorías genotípicas
d_3 <- which(sapply(Variants, function(y){
  y <- factor(y)
  nlevels(y) - is.element("0/0",levels(y)) * 1}) >= 3)
Variants_3 <- Variants[, d_3]

# Inicialización de matrices para almacenar resultados
Genetic_3 <- data.frame(SNP_name = colnames(Variants_3), 
                        codominant = 0, dominant = 0, recessive = 0, 
                        overdominant = 0, log_additive = 0, default = 0)
Genetic_3_g <- setNames(data.frame(matrix(ncol = length(colnames(Variants_3)), 
                                          nrow = nrow(Variants_3))), colnames(Variants_3))

models <- c("codominant", "dominant", "recesive", "overdominant", "log_additive")

for (i in 1:ncol(Variants_3)) {
  u0 <- as.character(Variants_2[, i])
  c <- strsplit(u0, "/")
  c <- as.data.frame(c)
  c1 <- unname(unlist(c[1,]))
  c2 <- unname(unlist(c[2,]))

  h <- "0"
  sorted_table <- sort(table(c(c1, c2)), decreasing = TRUE)
  dominant <- names(sorted_table[which(names(sorted_table) != h)][1])
  v0 <- unique(c(c1, c2))
  alternative <- v0[v0 != dominant & v0 != 0] 

  homozygous_dominant <- paste(dominant, dominant, sep = "/")
  heterozygous <- paste(alternative, dominant, sep = "/")
  heterozygous_1 <- paste(dominant, alternative, sep = "/")
  homozygous_alternative <- paste(alternative, alternative, sep = "/")

  tipo <- gsub(homozygous_dominant, 1, u0)
  tipo <- gsub(heterozygous, 0, tipo)
  tipo <- gsub(heterozygous_1, 0, tipo)
  tipo <- gsub(homozygous_alternative, 2, tipo)
  tipo <- gsub("0/0", NA, tipo)

  for (j in models) {
    if (j == "codominant") {
      tipo_1 <- factor(tipo, levels = 0:2,
                       labels = c(homozygous_alternative,
                                  homozygous_dominant, heterozygous))
      y <- table(tipo_1, Response)
      y1 <- chisq.test(y)
      Genetic_3[i, 2] <- y1$p.value
    }

    if (j == "dominant") {
      tipo_2 <- gsub(2, 1, tipo)
      tipo_2 <- factor(tipo_2, levels = 0:1,
                       labels = c(homozygous_dominant,
                                  paste(homozygous_dominant, heterozygous, sep = " ")))
      y <- table(tipo_2, Response)
      y1 <- chisq.test(y)
      Genetic_3[i, 3] <- y1$p.value
    }

    if (j == "recesive") {
      tipo_3 <- gsub(2, 0, tipo)
      tipo_3 <- factor(tipo_3, levels = 0:1,
                       labels = c(paste(homozygous_alternative, heterozygous, sep = " "),
                                  homozygous_dominant))
      y <- table(tipo_3, Response)
      y1 <- chisq.test(y)
      Genetic_3[i, 4] <- y1$p.value
    }

    if (j == "overdominant") {
      tipo_4 <- gsub(0, 1, tipo)
      tipo_4 <- factor(tipo_4, levels = 1:2,
                       labels = c(paste(homozygous_alternative, homozygous_dominant, sep = " "),
                                  heterozygous))
      y <- table(tipo_4, Response)
      y1 <- chisq.test(y)
      Genetic_3[i, 5] <- y1$p.value
    }

    if (j == "log_additive") {
      tipo_5 <- factor(tipo, levels = 0:2, ordered = TRUE,
                       labels = c(homozygous_alternative, homozygous_dominant, heterozygous))
      y <- table(tipo_5, Response)
      y1 <- CochranArmitageTest(y)
      Genetic_3[i, 6] <- y1$p.value
    }
  }

  Genetic_3[i, 7] <- NA
  model <- which.min(Genetic_3[i, ]) 
  Genetic_3_g[, i] <- switch(as.character(model), tipo_1, tipo_2, tipo_3, tipo_4, tipo_5)
}

correction_3 <- data.frame(rs_correction(Genetic_3, 1, method = "BH"))
correction_3[, 1] <- colnames(Variants_3)
```

Guardamos los resultados

```{r guardar_resultados}
# Combine SNPs with 2 and 3 categories
Genetic_2_3 <- rbind(Genetic_2, Genetic_3)
correction_variant <- data.frame(rs_correction(Genetic_2_3, 1, method = "BH"))
correction_variant[, 1] <- c(colnames(Variants_2), colnames(Variants_3))
Genetic <- data.frame(cbind(Genetic_2_g, Genetic_3_g))

# Save the results
save(Genetic, file = "~/Nextcloud/RebecaMeleroTesis/Hortega/data/Genetic_ITI_3.rda")
```

#### Merge clinical and genetic data

```{r}
orden = match(dataframe_i$ID,clinica_i$Paciente_ID)
ITI_m = data.frame(clinica_i[orden,],dataframe_i)
table(dataframe_i$ID == clinica_i$Paciente_ID[orden])

ITI_m = data.frame(ITI_m)
save(ITI_m, file = "~/Nextcloud/RebecaMeleroTesis/snpspanel/data/ITI_m_2.rda")
```

### Definition of the study population

```{r}
Name_variables = c("Gender(F/M)","Age","BMI", "HDL (mg/dL)", "CT (mg/dL)", "TG (mg/dL)", "LDL (mg/dL)" ,"BMI f","HDL (mg/dL) f" ,"CT (mg/dL) f", "TG (mg/dL) f", "LDL (mg/dL) f")

z = list(Gender = ITI_m$Sex ,
Age = ITI_m$age[ITI_m$DM == 1],
BMI = ITI_m$IMC_Kg.m2[ITI_m$DM == 1],
HDL = ITI_m$HDLc_mg.dL[ITI_m$DM == 1],
CT = ITI_m$CT_mg.dL[ITI_m$DM == 1],
TG = ITI_m$TG_mg.dL[ITI_m$DM == 1],
LDL = ITI_m$LDLc_mg.dL[ITI_m$DM == 1],
BMI_s = ITI_m$IMC_Kg.m2_sgto[ITI_m$DM_sgto == 1],
HDL_s = ITI_m$HDLc_mg.dL_sgto[ITI_m$DM_sgto== 1],
CT_s = ITI_m$CT_mg.dL_sgto[ITI_m$DM_sgto == 1],
TG_s = ITI_m$TG_mg.dL_sgto[ITI_m$DM_sgto == 1],
LDL_s = ITI_m$LDLc_mg.dL_sgto[ITI_m$DM_sgto == 1]

)

z0 = list(Gender = ITI_m$Sex,
Age_0 =  ITI_m$age[ITI_m$DM == 0],
BMI_0 = ITI_m$IMC_Kg.m2[ITI_m$DM == 0],
HDL_0 = ITI_m$HDLc_mg.dL[ITI_m$DM == 0],
CT_0 = ITI_m$CT_mg.dL[ITI_m$DM == 0],
TG_0 = ITI_m$TG_mg.dL[ITI_m$DM == 0],
LDL_0 = ITI_m$LDLc_mg.dL[ITI_m$DM == 0],
BMI_s0 = ITI_m$IMC_Kg.m2_sgto[ITI_m$DM_sgto == 0],
HDL_s0 = ITI_m$HDLc_mg.dL_sgto[ITI_m$DM_sgto== 0],
CT_s0 = ITI_m$CT_mg.dL_sgto[ITI_m$DM_sgto == 0],
TG_s0 = ITI_m$TG_mg.dL_sgto[ITI_m$DM_sgto == 0],
LDL_s0 = ITI_m$LDLc_mg.dL_sgto[ITI_m$DM_sgto == 0]
)

Gender = ITI_m$Sex
Response = ITI_m$DM

definitions = data.frame(Variables = Name_variables,
                         Diabetics = rep(0, length(Name_variables)),
                         Controls = rep(0, length(Name_variables)))
Aso_definitions = setNames(data.frame(matrix(data = NA, 
                                             nrow = 1, 
                                             ncol = length(Name_variables)))
                           ,Name_variables)


d = table(Gender,Response)
d = d[-1,]
x = d[1,1]
y = d[2,1]
xy = paste0(x,"/",y)
definitions$Controls[1] <- xy

x1= d[1,2]
y1 = d[2,2]
xy1 =  paste0(x1,"/",y1)
definitions$Diabetics[1] <- xy1

t = chisq.test(d)
Aso_definitions[1] <- t$p.value

for ( i in 2:length(z)){
x = round(mean(z[[i]], na.rm = TRUE),1)
y = round(sd(z[[i]], na.rm = TRUE),1)
xy = paste0(x,"(",y,")")
definitions$Diabetics[i] <- xy}

for ( i in 2:length(z0)){
x = round(mean(z0[[i]], na.rm = TRUE),1)
y = round(sd(z0[[i]], na.rm = TRUE),1)
xy = paste0(x,"(",y,")")
definitions$Controls[i] <- xy

t = t.test(
  x           = z[[i]],
  y           = z0[[i]],
  alternative = "two.sided",   
  mu          = 0,
  var.equal   = TRUE,
  conf.level  = 0.95
)
Aso_definitions[i]<- t$p.value
xtable(definitions)
}
definitions_p = data.frame(definitions, t(Aso_definitions[1,]), p.adjust(t(Aso_definitions[1,]), method = "BH"))
```

#### Analysis of the association between SNPs and lipid levels

```{r}
Variants = ITI_genetic[,125:420]
#Replace "0/0"  por NA
for ( i in colnames(Variants)){
 x = as.character(Variants[,i]) 
 x[x == "0/0"] <- NA
 Variants[,i] = factor(x)
}

d_2 =  which(sapply(Variants,function(y){
    y = factor(y)
    nlevels(y) - is.element("0/0",levels(y))*1}) < 3)

  Variants_2 = Variants[,d_2]
  
  d_1 = which(sapply(Variants_2,function(y){
    y = factor(y)
    nlevels(y) - is.element("0/0",levels(y))*1}) >= 2)
  
  Variants_2 =   Variants_2[,d_1]

 # Select variables with three categories:

  d_3 = which(sapply(Variants,function(y){
    y = factor(y)
    nlevels(y) - is.element("0/0",levels(y))*1}) >= 3)
  
  
  
 
 Aso_covariables_2 = data.frame(SNPs = colnames(Variants_2), 
                                p.valor_HDL = rep(0, ncol(Variants_2)), 
                                p.valor_CT = rep(0, ncol(Variants_2)), 
                                p.valor_TG = rep(0, ncol(Variants_2)), 
                                p.valor_LDL = rep(0, ncol(Variants_2)))

#f = function(Variants_2,data,j){
  
for (j in 1:ncol(Variants_2)) {
D_HDL = list()
D0_HDL = list()

D_CT = list()
D0_CT= list()

D_TG = list()
D0_TG = list()

D_LDL = list()
D0_LDL = list()

for ( i in 1:nrow(Variants_2)){
 z = levels(Variants_2[,j])
    x1 = z[1]
    x2 = z[2]
   
  if (is.na(Variants_2[i,j])){
  D_HDL[[i]] = NA
  D0_HDL[[i]] = NA 
  
  D_CT[[i]] = NA 
  D0_CT[[i]] = NA

  D_TG[[i]] = NA  
  D0_TG[[i]] = NA  
  
  D_LDL[[i]] = NA 
  D0_LDL[[i]] = NA 
      
  } else if (Variants_2[i,j] == x1){
  D1 = data$HDLc_mg.dL[i] 
  D_HDL[[i]] = D1
  D0_HDL[[i]] = NA 
  
  D2 = data$CT_mg.dL[i]
   D_CT[[i]] = D2
  D0_CT[[i]] = NA
 
  D3 = data$TG_mg.dL[i]
  D_TG[[i]] = D3
  D0_TG[[i]] = NA 
  
  D4 = data$LDLc_mg.dL[i]
  D_LDL[[i]] = D4
  D0_LDL[[i]] = NA 
  
   }  else if (Variants_2[i,j] == x2){
    D1 = data$HDLc_mg.dL[i] 
    D_HDL[[i]] = NA
    D0_HDL[[i]] =  D1
    
    D2 = data$CT_mg.dL[i]
    D_CT[[i]] = NA
    D0_CT[[i]] = D2
 
    D3 = data$TG_mg.dL[i]
    D_TG[[i]] = NA
    D0_TG[[i]] = D3
    
     D4 = data$LDLc_mg.dL[i]
     D_LDL[[i]] = NA
     D0_LDL[[i]] = D4 
   }
}
t1 = t.test(
  x           = unlist(D_CT),
  y           = unlist(D0_CT),
  alternative = "two.sided",
  mu          = 0,
  var.equal   = TRUE,
  conf.level  = 0.95
)
 Aso_covariables_2$p.valor_CT[j]<- t1$p.value
}
  t3 = t.test(
  x           = unlist(D_LDL),
  y           = unlist(D0_LDL) ,
  alternative = "two.sided",
  mu          = 0,
  var.equal   = TRUE,
  conf.level  = 0.95
)
Aso_covariables_2$p.valor_LDL[j]<- t3$p.value

 
 
 t = t.test(
  x           = unlist(D_HDL),
  y           = unlist(D0_HDL),
  alternative = "two.sided",
  mu          = 0,
  var.equal   = TRUE,
  conf.level  = 0.95
)
 Aso_covariables_2$p.valor_HDL[j]<- t$p.value
 
  t2 = t.test(
  x           = unlist(D_TG),
  y           = unlist(D0_TG) ,
  alternative = "two.sided",
  mu          = 0,
  var.equal   = TRUE,
  conf.level  = 0.95
)
 Aso_covariables_2$p.valor_TG[j]<- t2$p.value
 
   t3 = t.test(
  x           = unlist(D_LDL),
  y           = unlist(D0_LDL) ,
  alternative = "two.sided",
  mu          = 0,
  var.equal   = TRUE,
  conf.level  = 0.95
)
Aso_covariables_2$p.valor_LDL[j]<- t3$p.value
}
Aso_covariables_p.adjust = data.frame( Aso_covariables_2$SNPs, rs_correction(Aso_covariables_2,2, method = "BH"))


 View(Aso_covariables_2)
save(Aso_covariables_2, file = "~/Nextcloud/RebecaMeleroTesis/snpspanel/data/Aso_covariables_p.adjust_ITI.rda")

 definitions_p = data.frame(definitions, t(Aso_definitions[1,]), p.adjust(t(Aso_definitions[1,]), method = "BH"))
```

REGRESION LOGISTICA, VARIABLE RESPUESTA DIABETES A LOS 8 AÑOS:

```{r}
glmSNPs_2 = function(response,SNP,vars,data){ 

   ff = as.formula(paste0(response," ~ ",paste0("(",paste0("rs1240521563",collapse="+"),")","*","(",
                                               paste0(vars,collapse = "+"),")")))
data0 = data[,c(response,"rs1240521563",vars)]
data0 = data0[complete.cases(data0),]
data0[,vars] = as.numeric(as.character(data0[,vars]))
data0[,response] = as.numeric(as.character(data0[,response]))
fit0 = glm(ff,data = data0)
fit0.s = step(fit0,direction="both",trace=0)
fit0.s
}

b = coef(summary(fit0.s))
  x = data.frame(SNP = SNP[i],  coeficientes = b[,1],p.valores=b[,4])
```

```{r}
response = "DM"
vars_1 = c("IMC_Kg.m2","Sex","age")

#"HDLc_mg.dL", "CT_mg.dL", "TG_mg.dL",  "LDLc_mg.dL""LDLc_mg.dL"
vars =  "LDLc_mg.dL"
SNP = colnames(data[128:667])

for ( i in 1:ncol(data[128:667])){
 x = as.character(data[,i]) 
 x[x == "./."] <- NA
 data[,i] = factor(x)
}
```

```{r}
l2_1 = NULL

for ( i in 3:length(SNP)){
  res= glmSNPs_2(response = response ,SNP = SNP[i],
                 vars= vars ,data= data)
  b = coef(summary(res))
  x = data.frame(SNP = SNP[i],  coeficientes = b[,1],p.valores=b[,4])
 l2_1= rbind(l2_1,x)
}

```

```{r}
l3_1 = data.frame( SNP = l2_1[,1] ,coeficientes = l2_1[,2],p.valores=l2_1[,3],
                 p.valores.ajustados =  rs_correction(l2_1,2,method = "BH")[,3])  
     x2_1 = grep(":" , rownames(l3_1))
     x3_1 = l3_1[x2_1,]
     s =  which(x3_1$p.valores.ajustados < 0.05)
     x4_1 = x3_1[s,]
     SNPs_sig = unique(x4_1$SNP)
     x5_1 = x3_1[-s,]
     SNP_no_sig = unique(x5_1$SNP)
```

```{r}
Adjust_confuse_variable = function(response,SNP,vars,vars_1,data){
  data0 = data[,c(response,SNP,vars,vars_1)]
  data0 = data0[complete.cases(data0),]

 ff = as.formula(paste0(response," ~ ",paste0("(",paste0(vars,collapse = "+"),")", "+","(",
                                               paste0(vars_1,collapse = "+"),")")))
                                               
 ff_1 = as.formula(paste0(response," ~ ",paste0("(",paste0(SNP,collapse="+"),")","+","(",
                                               paste0(vars,collapse = "+"),")", "+","(",
                                             paste0(vars_1,collapse = "+"),")")))
                                          
 fit = glm(ff,data = data0)
 fit_1 = glm(ff_1 ,data = data0)
 x = anova(fit,fit_1, test = "Chisq")
 z = x$`Pr(>Chi)`
 z[!(is.na(z))]
}

Adjust_confuse_variable_v = Vectorize(Adjust_confuse_variable, vectorize.args = "SNP")

```

```{r}
z2 = Adjust_confuse_variable_v (response = response, SNP = SNPs_sig , vars = vars, vars_1 = vars_1, data = data)
z3 = p.adjust(z2, method = "BH")
d1 = which(z3 <= 0.05)
SNPs_ajustadas_significativas_1 = names(d1)
h = data.frame(z3)
h1 = data.frame(SNP = SNPs_sig, p.value = z2, adjusted.p.value = h$z3)
h1$adjusted.p.value
print(xtable(h1, digits = c(5,5,5,5)), include.rownames= FALSE)
```
```{r}
# Create the dataframe to be filled:

df = data.frame(pair=rep(1:4755,rep(2,4755)), respuesta = rep(0,4755), 
                   HDLc= rep(0,4755),CT= rep(0,4755),
                   TG = rep(0,4755),LDLc = rep(0,4755),
                   Age = rep(0,4755) , Sex = rep(0,4755), BMI = rep(0,4755),
                   SNP_1 = rep(0,4755))

# Define the response variable (y1) y las variables explicativas(x's), 
# teniendo en cuenta que las variables "_s" son de seguimiento.

y1 = data$DM
y1_s = data$DM_sgto
x1 = data$HDLc_mg.dL
x1_s = data$HDLc_mg.dL_sgto
x2 = data$CT_mg.dL
x2_s = data$CT_mg.dL_sgto
x3 = data$TG_mg.dL
x3_s = data$TG_mg.dL_sgto
x4 = data$LDLc_mg.dL
x4_s = data$LDLc_mg.dL_sgto

s = as.factor(data$Sex)
which(data$Sex == 0)
s_1 = as.factor(data$Sex)
a = data$age
a_s = data$age_s
b = data$IMC_Kg.m2
b_s = data$IMC_Kg.m2_sgto


df[,1] = rep(1:4755,rep(2,4755))
df[2*(1:4755)-1,2] = y1
df[2*(1:4755),2] = y1_s
df[2*(1:4755)-1,3] = x1
df[2*(1:4755),3] = x1_s
df[2*(1:4755)-1,4] = x2
df[2*(1:4755),4] = x2_s
df[2*(1:4755)-1,5] = x3
df[2*(1:4755),5] = x3_s
df[2*(1:4755)-1,6] = x4
df[2*(1:4755),6] = x4_s
df[2*(1:4755)-1,7] = a
df[2*(1:4755),7] = a_s
df[2*(1:4755)-1,8] = s
df[2*(1:4755),8] = s
df[2*(1:4755)-1,9] = b
df[2*(1:4755),9] = b_s


```

```{r}
tabla_gee_2 = function(tabla0_4_id, name_SNP){
  control <- glmerControl(optimizer = "bobyqa",
                        optCtrl = list(maxfun = 10000000),
                        tolPwrss = 1e-10, # Tolerancia para la convergencia en pwrss
                        restart_edge = FALSE) # Evita reinicios si está en el borde

  fit1 = glmer(response ~ Age + Sex + BMI + (LDLc * time * SNP_1) + (1|ID), 
               data = tabla0_4_id, family = binomial, control = control)
  #t = stepcAIC(fit1,direction="both",trace=0)
x = summary(fit1)
x1 = data.frame(x$coefficients, SNP = rep(name_SNP, nrow(x$coefficients)))
fit1_p = x1 
fit1_p
}

l2 = NULL
#table_pairs = function(tabla0_2 ,z ,duplicados_totales, SNP){
library("Matrix")
library("lme4")

SNPs = data[,128:667]

# Replace ./. por "0/0"
for ( i in 1:ncol(SNPs)){
 x = as.character(SNPs[,i]) 
 x[x == "./."] <- "0/0"
 SNPs[,i] = factor(x)
}

levels(as.factor(SNPs[,"rs139953107"]))
```

```{r}

for (u in 534:length(SNPs)){
name_SNP =colnames(SNPs)[u]
SNP = SNPs[,u]
#df[2*(1:4755) -1,10] = SNP
df[,10] = SNP
# Remove individuals containing NAs.

tabla0 = df[complete.cases(df),]

# Select individuals without NAs at both time points:
  
z = which(duplicated(tabla0$pair) == TRUE)
duplicados = as.numeric(rownames(tabla0[z,]))
mas_uno = rep(1,length(z))
duplicados_mas_uno = duplicados - mas_uno 
duplicados_totales = as.numeric(c(duplicados,duplicados_mas_uno))
duplicados_totales = duplicados_totales[order(duplicados_totales)]

tabla0_2 = df[duplicados_totales,]
n = nrow(tabla0_2)
n = n/2
b = rep(1:n,rep(2,n))
tabla0_2$pair = b

## Keep only complete cases for SNP.
tabla0_4 = tabla0_2[complete.cases(tabla0_2),]

time = rep(c(0,1), nrow(tabla0_4) + 1 /2 )

tabla0_4_id = data.frame(ID = as.numeric(tabla0_4$pair), response = as.factor(tabla0_4$respuesta), SNP_1 = as.factor(tabla0_4$SNP_1),  LDLc = as.numeric(tabla0_4$LDLc), time = as.factor(time), Age = as.numeric(tabla0_4$Age) , Sex = as.factor(tabla0_4$Sex), BMI = as.numeric(tabla0_4$BMI))

s = data.frame(table(tabla0_4_id$response, tabla0_4_id$SNP_1))
 if(all(s$Freq) != 0 & nrow(s) >= 4){
   print(u)
 tryCatch({fit1_p1 = tabla_gee_2(tabla0_4_id = tabla0_4_id, name_SNP = name_SNP)})
    l2 = rbind(l2,fit1_p1)
 }
  else { print(NULL)}
}

```


```{r}
l3_1 = data.frame( SNP = l2[,5] ,coeficientes = l2[,1],p.valores=l2[,4], p.valores.ajustados =  rs_correction(l2, dim = 2, method = "BH")[,4]) 
     x2_1 = grep(":" , rownames(l3_1))
     x3_1 = l3_1[x2_1,]
     s =  which(x3_1$p.valores.ajustados < 0.05)
     x4_1 = x3_1[s,]
     SNPs_sig = unique(x4_1$SNP)
     x5_1 = x3_1[-s,]
     SNP_no_sig = unique(x5_1$SNP)
     
     x4_2 = x4_1[x4_1[, 3] != 0, ]
     f = grep("time1:" , rownames(x4_2))
     x4_3 = x4_2[f,]
   
     f = grep(":time1:" , rownames(x4_2))
     y4_3 = x5_1[f,]
      
```
for (z in 2:length(unique(x4_3$SNP))){

w = unique(x4_3$SNP)[z]
SNP = SNPs[,w]

#df[2*(1:4755) -1,10] = SNP
df[,10] = SNP

for ( i in 1:length(df[,10])){
 x = as.character(df[i,10]) 
 x[x == "0/0"] <- NA
 df[i,10] = as.factor(x)
}
# Remove individuals containing NAs.

tabla0 = df[complete.cases(df),]

# Select individuals without NAs at both time points:
  
z = which(duplicated(tabla0$pair) == TRUE)
duplicados = as.numeric(rownames(tabla0[z,]))
mas_uno = rep(1,length(z))
duplicados_mas_uno = duplicados - mas_uno 
duplicados_totales = as.numeric(c(duplicados,duplicados_mas_uno))
duplicados_totales = duplicados_totales[order(duplicados_totales)]

tabla0_2 = df[duplicados_totales,]
n = nrow(tabla0_2)
n = n/2
b = rep(1:n,rep(2,n))
tabla0_2$pair = b
levels(tabla0$SNP_1)
 
## Keep only complete cases for SNP.
tabla0_4 = tabla0_2[complete.cases(tabla0_2),]

time = rep(c(0,1), nrow(tabla0_4) + 1/2)

tabla0_4_id = data.frame(ID = as.numeric(tabla0_4$pair), response = as.factor(tabla0_4$respuesta), SNP_1 = as.factor(tabla0_4$SNP_1),  LDLc = as.numeric(tabla0_4$LDLc), time = as.factor(time), Age = as.numeric(tabla0_4$Age) , Sex = as.factor(tabla0_4$Sex), BMI = as.numeric(tabla0_4$BMI))

 control <- glmerControl(optimizer = "bobyqa",
                        optCtrl = list(maxfun = 10000000),
                        tolPwrss = 1e-10, # Tolerancia para la convergencia en pwrss
                        restart_edge = FALSE) # Evita reinicios si está en el borde

  fit1 = glmer(response ~ Age + Sex + BMI +  (LDLc * time * SNP_1) + (1|ID), 
               data = tabla0_4_id, family = binomial, control = control)


d1 = "~/Downloads/"
vars = "LDL_rep2_" 
probabilidad_enfermedad = "diabetes"
  filename = paste0(d1,"/",probabilidad_enfermedad,vars,w ,".png")
  png(filename = filename)
ggpredict(fit1, terms = c("LDLc[all]","time","SNP_1" )) %>% plot()
print(plot)
}
