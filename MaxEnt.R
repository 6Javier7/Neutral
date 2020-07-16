library(spocc)
library(raster)
library(dismo)
library(rJava)
library(rgeos)
library(grDevices)
library(mapr)
library(ENMeval)
library(ROCR)
library(vcd)
library(sp)
library(maptools)

data(wrld_simpl)
clima <- getData("worldclim", var = 'bio', res = 2.5)

clima1 <- stack(clima)
bio1 <- clima1[[1]]

nrow(bio1)
extent(bio1)
head(values(bio1))
length(values(bio1))

#Registros
ocurrencias <- gbif(genus = "Coenogonium", species = "linkii")
write.csv(ocurrencias, "registros")
ocurrencias <- subset(ocurrencias,!is.na(lat) & !is.na(lon))
ocurrencias <- ocurrencias[!duplicated(ocurrencias[c("lat", "lon"),])]

table(ocurrencias$basisORecord)
preservados <- subset(ocurrencias, basisOfRecord == "PRESERVED_SPECIMEN")
coordinates(preservados) <- ~lon + lat #espaciales para graficar
coordinates(ocurrencias) <- ~lon + lat #espaciales para graficar
src <- CRS("+init=epsg:4326")
crs(preservados) <- src
crs(ocurrencias) <- src

plot(bio1[[1]])
plot(preservados, add = T) #Analizar datos extremos

condiciones <- extract(clima, preservados)
head(condiciones)

sindatos <- is.na(condiciones[,1])
table(sindatos)

mapa <- spTransform(wrld_simpl, crs(preservados))

sobrelapan <- over(preservados, mapa)

#para evitar el ruido es preferible solo un punto por celda
celdas <- cellFromXY(clima[[1]], preservados)
celdas1 <- cellFromXY(clima[[1]], ocurrencias)

datos <- preservados[!duplicated(celdas),]
datos1 <- ocurrencias[!duplicated(celdas),]

#Buffer
bufer <- buffer(datos, width = 400000) #bufer de 400km

#jugar con el mapa

areas <- crop(clima, extent(bufer)) #un rectangulo del area de los puntos
areas <- mask(areas, bufer) #el area exacta que coge todos los registros

set.seed(0728)
bg1 <- randomPoints(areas, 8000, p = datos)
bg <- sampleRandom(x = areas, size = 10000, na.rm = T, sp = T) 
plot(areas[[1]])
plot(bg, add = T)
plot(datos, add = T, col = "red")

set.seed(0728)
azar <- kfold(datos, k = 5)
ram <- sample(1:5, 1)
test <- datos[azar == ram,]
train <- datos[azar != ram,]

#se extraen las condiciones para la prueba, el training y el entorno
p <- extract(clima, train)
ptest <- extract(clima, test)
a <- extract(clima, bg)
a1<- extract(clima, bg1)

label <- c(rep(1, nrow(p)), rep(0, nrow(a))) #presente en el test y ausente en el entorno
combinado <- as.data.frame(rbind(p, a))

#model
mod <- maxent(x = combinado, p = label, args = c("responsecurves"))

#evaluaciones
evaltrain <- evaluate(p = ptest, a = a, model = mod)
plot(evaltrain, "ROC")

print(evaltrain)
response(mod)

#predicciones

ped1 <- predict(mod, areas)
ped2 <- predict(mod, clima)

#desviacion estandar

library(boot)

AUC1 <- function(data = list(ptest, a), i) {
  pres <- ptest[i]
  aus <- a
  combinado <- c(pres, aus)
  label <- c(rep(1, length(pres)), rep(0, length(aus)))
  predic <- prediction(combinado, label)
  perfor <- performance(predic, "auc")
  return(perfor@y.values[[1]])
}

AUC2 <- function(data = list(ptest, a, mod), i) {
  pres <- p[i,]
  aus <- a
  eval <- evaluate(p = pres, a = aus, model = mod)
  return(eval@auc)
}

bo1 <- boot(data = list(ptest, a), AUC1, 1000)
bo2 <- boot(data = list(ptest, a), AUC2, 1000)
bo3 <- boot(data = list(ptest, a1), AUC1, 1000)
bo4 <- boot(data = list(ptest, a1), AUC2, 1000)

bo1
bo3
bo2
bo4
plot(bo3)
plot(bo4)
boot.ci(bo3)
boot.ci(bo4)
plot(ped2)
