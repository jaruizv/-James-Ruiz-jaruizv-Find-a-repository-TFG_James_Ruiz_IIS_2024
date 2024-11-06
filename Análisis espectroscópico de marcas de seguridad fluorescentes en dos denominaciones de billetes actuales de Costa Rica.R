#Se importan las bibliotecas necesarias
library(tidyverse) # Conjunto de paquetes para manipulación de datos
library(photobiology)# Manipulación de espectros
library(caTools) # Cálculo numérico
library(lmtest) # Prueba de los supuestos del modelo
library(stats) # Paquete de funciones estadísticas
library(car) # Paquete para regresión aplicada
library(ARTool) # Paquete de estadística no paramétrica
library(lme4) # Modelo lineal de efectos mixtos
library(rcompanion) # Funciones para evaluación de datos extendida


# FUNCIONES #

mean_spectrum_INT <- function(df){
  # Descripción: La función determina el promedio ponderado integral de la intensidad
  # Entradas
  #  df = data frame del espectro
  # Salidas
  # promedio ponderado integral de la intensidad
  
  # Se emplea la regla del trapecio para determinar la integral
  mean = trapz(df$nm, df$INT) / (max(df$nm)-min(df$nm))
  return(mean)
}

centroid_wavelength_nm <- function(df){
  # Descripción: La función determina el promedio ponderado de la longitud de onda
  # Entradas
  #  df = data frame del espectro
  # Salidas
  # promedio ponderado integral de la longitud de onda
  centroid = sum(df$INT*df$nm)/sum(df$INT)
  return(centroid)
}

localpeaks <- function(df) {
  # Descripción: La función determina los picos máximos y su intensidad de cada espectro
  # Entradas
  #  df = data frame del espectro
  # Salidas
  # Longitud e intensidad del pico máximo de emisión
  localPeaks <- peaks(df,
                      strict = TRUE,
                      span = 5, 
                      x.var.name = "nm", 
                      y.var.name = "INT")
  # Se ordenan los máximos de mayor a menor según su intensidad de emisión
  locPeaks <- localPeaks[order(-localPeaks$INT),]
  return(locPeaks)
}


### DISEÑO EXPERIMENTAL EN BILLETES DE DOS MIL COLONES ###

#BILLETES DE 2 MIL: POSICIÓN DE LA MARCA FLUORESCENTE

par(mfrow = c(2, 2), cex.axis = 1.2, cex.lab = 1.2, cex.main = 1.2)

pos_X_esq_inf_izq = data.frame(Año =c("2020","2020","2020","2020","2020",
                                      "2021","2021","2021","2021","2021",
                                      "2022","2022","2022","2022","2022"))
pos_X_esq_inf_izq$Año = factor(pos_X_esq_inf_izq$Año, levels = c("2020", "2021", "2022"))
pos_X_esq_inf_izq$Posición = c(65.37,
                               65.22,
                               66.15,
                               67.30,
                               67.32,
                               67.57,
                               66.85,
                               66.85,
                               66.22,
                               65.71,
                               65.54,
                               68.67,
                               67.46,
                               65.94,
                               66.00)
par(mar = c(2.5, 5, 2, 0))
boxplot(data = pos_X_esq_inf_izq, Posición ~ Año, ylab="Posición en X / (mm)", xlab = NULL, main = "Esquina inferior izquierda", col = rgb(0.5, 0.5, 1, alpha = 0.5))


pos_X_esq_inf_der = data.frame(Año =c("2020","2020","2020","2020","2020",
                                      "2021","2021","2021","2021","2021",
                                      "2022","2022","2022","2022","2022"))
pos_X_esq_inf_der$Año = factor(pos_X_esq_inf_der$Año, levels = c("2020", "2021", "2022"))
pos_X_esq_inf_der$Posición = c(79.07,
                               79.18,
                               80.42,
                               81.44,
                               82.44,
                               82.25,
                               81.35,
                               81.52,
                               80.00,
                               79.21,
                               79.04,
                               84.33,
                               82.58,
                               80.62,
                               80.10)
par(mar = c(2.5, 3, 2, 2))
boxplot(data = pos_X_esq_inf_der, Posición ~ Año, ylab = NULL, xlab = NULL, main = "Esquina inferior derecha", col = rgb(0.5, 0.5, 1, alpha = 0.5))


pos_Y_esq_inf_izq = data.frame(Año =c("2020","2020","2020","2020","2020",
                                      "2021","2021","2021","2021","2021",
                                      "2022","2022","2022","2022","2022"))
pos_Y_esq_inf_izq$Año = factor(pos_Y_esq_inf_izq$Año, levels = c("2020", "2021", "2022"))
pos_Y_esq_inf_izq$Posición = c(29.59,
                               28.84,
                               28.04,
                               28.48,
                               27.54,
                               27.96,
                               28.54,
                               28.12,
                               29.07,
                               28.57,
                               28.75,
                               27.33,
                               27.14,
                               28.37,
                               28.13)
par(mar = c(4, 5, 0.5, 0))
boxplot(data = pos_Y_esq_inf_izq, Posición ~ Año, ylab="Posición en Y / (mm)", col = rgb(0.5, 0.5, 1, alpha = 0.5))

pos_Y_esq_inf_der = data.frame(Año =c("2020","2020","2020","2020","2020",
                                      "2021","2021","2021","2021","2021",
                                      "2022","2022","2022","2022","2022"))
pos_Y_esq_inf_der$Año = factor(pos_Y_esq_inf_der$Año, levels = c("2020", "2021", "2022"))
pos_Y_esq_inf_der$Posición = c(29.74,
                               29.07,
                               27.96,
                               27.96,
                               27.38,
                               27.92,
                               28.54,
                               28.04,
                               28.93,
                               28.43,
                               28.54,
                               27.42,
                               26.62,
                               28.37,
                               27.90)
par(mar = c(4, 3, 0.5, 2))
boxplot(data = pos_Y_esq_inf_der, Posición ~ Año, ylab = NULL, col = rgb(0.5, 0.5, 1, alpha = 0.5))

# Se define un nivel de significancia de 0.05
# Prueba de normalidad
shapiro.test(pos_X_esq_inf_izq$Posición) # p-value = 0.444 (Los datos son normales)
shapiro.test(pos_Y_esq_inf_izq$Posición) # p-value = 0.984 (Los datos son normales)
shapiro.test(pos_X_esq_inf_der$Posición) # p-value = 0.2932 (Los datos son normales)
shapiro.test(pos_Y_esq_inf_der$Posición) # p-value = 0.9411 (Los datos son normales)

# Prueba de igualdad de varianzas
leveneTest(pos_X_esq_inf_izq$Posición, pos_X_esq_inf_izq$Año, center = "mean") # p-value = 0.1939 (Varianzas iguales)
leveneTest(pos_Y_esq_inf_izq$Posición, pos_Y_esq_inf_izq$Año, center = "mean") # p-value = 0.4313 (Varianzas iguales)
leveneTest(pos_X_esq_inf_der$Posición, pos_X_esq_inf_der$Año, center = "mean") # p-value = 0.3175 (Varianzas iguales)
leveneTest(pos_Y_esq_inf_der$Posición, pos_Y_esq_inf_der$Año, center = "mean") # p-value = 0.1232 (Varianzas iguales)

# ANOVA de un factor (Factor temporal)
pos_X_esq_inf_izq.lm <- lm(Posición ~ Año, data = pos_X_esq_inf_izq)
anova(pos_X_esq_inf_izq.lm) # valor-p = 0.7708 (Medias iguales)

pos_Y_esq_inf_izq.lm <- lm(Posición ~ Año, data = pos_Y_esq_inf_izq)
anova(pos_Y_esq_inf_izq.lm) # valor-p = 0.3594 (Medias iguales)

pos_X_esq_inf_der.lm <- lm(Posición ~ Año, data = pos_X_esq_inf_der)
anova(pos_X_esq_inf_der.lm) # valor-p = 0.7343 (Medias iguales)

pos_Y_esq_inf_der.lm <- lm(Posición ~ Año, data = pos_Y_esq_inf_der)
anova(pos_Y_esq_inf_der.lm) # valor-p = 0.3429 (Medias iguales)

# Análisis de los residuos

# Normalidad de los residuos
shapiro.test(pos_X_esq_inf_izq.lm$residuals) # valor-p = 0.1991 (Residuos normales)
shapiro.test(pos_Y_esq_inf_izq.lm$residuals) # valor-p = 0.939 (Residuos normales)
shapiro.test(pos_X_esq_inf_der.lm$residuals) # valor-p = 0.7272 (Residuos normales)
shapiro.test(pos_Y_esq_inf_der.lm$residuals) # valor-p = 0.7575 (Residuos normales)

# Varianza constante
leveneTest(pos_X_esq_inf_izq.lm$residuals, pos_X_esq_inf_izq$Año, center = "mean") # valor-p = 0.1939 (Varianza constante)
leveneTest(pos_Y_esq_inf_izq.lm$residuals, pos_Y_esq_inf_izq$Año, center = "mean") # valor-p = 0.4313 (Varianza constante)
leveneTest(pos_X_esq_inf_der.lm$residuals, pos_X_esq_inf_der$Año, center = "mean") # valor-p = 0.3175 (Varianza constante)
leveneTest(pos_Y_esq_inf_der.lm$residuals, pos_Y_esq_inf_der$Año, center = "mean") # valor-p = 0.1232 (Varianza constante)

# Independencia de los residuos
bgtest(pos_X_esq_inf_izq.lm) # valor-p = 0.2303 (Residuos independientes)
bgtest(pos_Y_esq_inf_izq.lm) # valor-p = 0.7864 (Residuos independientes)
bgtest(pos_X_esq_inf_der.lm) # valor-p = 0.2462 (Residuos independientes)
bgtest(pos_Y_esq_inf_der.lm) # valor-p = 0.5506 (Residuos independientes)

dev.off()
# Verificación gráfica
par(mfrow = c(2, 2))
plot(pos_X_esq_inf_izq.lm)
plot(pos_Y_esq_inf_izq.lm)
plot(pos_X_esq_inf_der.lm)
plot(pos_Y_esq_inf_der.lm)

# IMPORTACIÓN DE DATOS ESPECTROSCÓPICOS#

# Las siguientes líneas leen los archivos .CSV y los convierte en data frames (df)

# Se define el directorio de trabajo
setwd("C:\\Users\\james\\Google Drive\\TEC\\2024\\II SEMESTRE 2024\\TFG\\Espectros\\Espectros Art2\\DOE 2 MIL")

# Se lee un archivo .csv con los códigos empleados para guardar cada espectro de emisión medido
# Dichos códigos fueron generados a partir del diseño factorial aleatorizado, el cual fue obtenido de Minitab
codes <- read.csv("Codigos_DOE.csv", header = TRUE, colClasses = c("Codigo.DOE" = "character"))
code_names = c()
variable_names = c()
rep <- c(1,2,3)
for (i in seq_along(codes)) {
  for (j in seq(1,3)){
    code_names <- c(code_names, paste(codes[[i]], "-", rep[j], sep = ""))
    variable_names <- c(variable_names, paste("e", codes[[i]], "_", rep[j], sep = ""))
  }
}

# Se leen y se guardan los espectros de acuerdo a los códigos generados anteriormente
for (i in seq_along(code_names)) {
  file_name <- paste(code_names[[i]], ".csv", sep = "")
  assign(variable_names[[i]], read.csv(file_name, skip = 1, header = TRUE))
}

# ANALISIS DE LOS ESPECTROS #

# Los resultados de los espectros se guardan en un data frame para su análisis
results_df <- data.frame(code_name = character(),
                         OrdenEst = numeric(),
                         OrdenCorrida = numeric(),
                         TipoPt = numeric(),
                         Bloques = numeric(),
                         Rep = numeric(),
                         lambda_max = numeric(),
                         lambda_ponderada = numeric(),
                         INT_max = numeric(),
                         INT_ponderada = numeric(),
                         stringsAsFactors = FALSE)

for (i in seq_along(variable_names)) {
  # Se obtiene el dataframe a analizar
  data_to_analyze <- get(variable_names[[i]])
  
  # Se determina los picos locales de máxima emisión
  peaks <- localpeaks(data_to_analyze)
  
  # Se promedian integralmente la intensidad medida
  mean_intensity <- mean_spectrum_INT(data_to_analyze)
  
  # Se promedian ponderadamente las longitudes de onda a partir de la intensidad
  centroid_wavelength <- centroid_wavelength_nm(data_to_analyze)
  
  # Se extraen los números de los códigos del DOE
  parts <- strsplit(code_names[[i]], "-")[[1]]
  numbers <- strsplit(parts[1], "")[[1]]
  
  # Se crea una nueva fila con los resultados de los picos en cada espectro
  new_row <- data.frame(
    code_name = code_names[[i]],
    OrdenEst = paste(numbers[1:2], collapse = ""),
    OrdenCorrida = paste(numbers[3:4], collapse = ""),
    TipoPt = paste(numbers[5:5], collapse = ""),
    Bloques = paste(numbers[6:6], collapse = ""),
    Rep = paste(parts[2:2], collapse = ""),
    lambda_max = peaks[1, "nm"],
    lambda_ponderada = centroid_wavelength,
    INT_max = peaks[1, "INT"],
    INT_ponderada = mean_intensity)
  
  # Se añade la nueva fila al dataframe de resultados
  results_df <- rbind(results_df, new_row)
}

# Se ordena por code_name
results_df <- results_df[order(results_df$code_name, decreasing = FALSE),]

# Se añaden las variables de los factores al dataframe de resultado
A <- c("2020", "2021", "2022") #Año de circulación
B <- c("A", "B", "C") # Posición
C <- c("1", "2", "3") # Orden de medición

# Se crea una combinación de los factores
combinations <- expand.grid(C= C, B = B, A = A)

# Se replican las combinaciones 5 veces (debido al factor bloque)
factors_df <- combinations[rep(seq_len(nrow(combinations)), times = 5), ]

# Se añaden las columnas al results_df
results_df$A = factors_df$A
results_df$B = factors_df$B

# Se añade un factor de muestras (1:15)
results_df$Muestras = rep(1:15, each = 9)

# Se guarda el results_df como un archivo ".csv"
write.csv(results_df, "Results_df.csv", row.names = FALSE)


# ANALISIS DISEÑO EXPERIMENTAL # alpha = 0.05 (95% de confianza)

# Histograma con gráfico de caja para cada variable de respuesta (Guardar figuras a 550x550)

dev.off()
# Lambda max
layout(matrix(1, nrow = 1))
par(mar = c(4, 4, 1, 1), cex = 1.2) 
# Gráfico principal: Histograma
hist(results_df$lambda_max, main = NULL, xlab = "Longitud de onda de máxima emisión / (nm)", 
     xlim = c(400, 700), ylab = "Densidad", ylim = c(-0.001, 0.025), freq = FALSE)
lines(density(results_df$lambda_max, bw = "nrd"), col = "blue", lwd = 2)
box()
# Gráfico secundario: diagrama de caja
par(new = TRUE, mar = c(4, 4, 18.5, 1))
boxplot(results_df$lambda_max, horizontal = TRUE, axes = FALSE, 
        ylim = c(400, 700), col = rgb(0.5, 0.5, 1, alpha = 0.5), at = -0.002)

dev.off()
# Lambda ponderada
layout(matrix(1, nrow = 1))
par(mar = c(4, 4, 1, 1), cex = 1.2) 
# Gráfico principal: Histograma
hist(results_df$lambda_ponderada, main = NULL, xlab = "Longitud de onda ponderada / (nm)", 
     xlim = c(515, 555), ylab = "Densidad", ylim = c(-0.01, 0.15), freq = FALSE, breaks = 12)
lines(density(results_df$lambda_ponderada, bw = "nrd"), col = "blue", lwd = 2)
box()
# Gráfico secundario: diagrama de caja
par(new = TRUE, mar = c(4, 4, 18.5, 1))
boxplot(results_df$lambda_ponderada, horizontal = TRUE, axes = FALSE, 
        ylim = c(515, 555), col = rgb(0.5, 0.5, 1, alpha = 0.5), at = -0.002)

dev.off()
# INT max
layout(matrix(1, nrow = 1))
par(mar = c(4, 4, 1, 1), cex = 1.2)
# Gráfico principal: Histograma
hist(results_df$INT_max, main = NULL, xlab = "Intensidad de emisión máxima / (CPS)",
     xlim = c(0,70000), ylab = "Densidad", ylim = c(-4e-06, 5e-05), freq= FALSE)
lines(density(results_df$INT_max, bw = "nrd"), col = "blue", lwd = 2)
box()
# Gráfico secundario: diagrama de caja
par(new = TRUE, mar = c(4, 4, 18.5, 1))
boxplot(results_df$INT_max, horizontal = TRUE, axes = FALSE, 
        ylim = c(0,70000), col = rgb(0.5, 0.5, 1, alpha = 0.5), at = -0.002)

dev.off()
# INT ponderada
layout(matrix(1, nrow = 1))
par(mar = c(4, 4, 1, 1), cex = 1.2) 
# Gráfico principal: Histograma
hist(results_df$INT_ponderada, main = NULL, xlab = "Intensidad de emisión ponderada / (CPS)",
     xlim = c(0,20000), ylab = "Densidad", ylim = c(-0.000008, 0.00012), freq= FALSE)
lines(density(results_df$INT_ponderada, bw = "nrd"), col = "blue", lwd = 2)
box()
# Gráfico secundario: diagrama de caja
par(new = TRUE, mar = c(4, 4, 18.5, 1))
boxplot(results_df$INT_ponderada, horizontal = TRUE, axes = FALSE, 
        ylim = c(0,20000), col = rgb(0.5, 0.5, 1, alpha = 0.5), at = -0.002)


dev.off()
# Gráficos de caja (Guardar 1200x675)

# Factores
results_df$A = factor(results_df$A, levels = c("2020", "2021", "2022"))
results_df$B = factor(results_df$B, levels = c("A", "B", "C"))

# Lambda max
par(cex = 1.6, mar = c(4,4,1,1))
boxplot(results_df$lambda_max ~ results_df$A * results_df$B, ylab="Longitud de onda de máxima emisión / (nm)", xlab = "Año de circulación y posición", main = NULL, col = rgb(0.5, 0.5, 1, alpha = 0.5))

# Lambda ponderada
par(cex = 1.6, mar = c(4,4,1,1))
boxplot(results_df$lambda_ponderada ~ results_df$A * results_df$B, ylab="Longitud de onda ponderada / (nm)", xlab = "Año de circulación y posición", main = NULL, col = rgb(0.5, 0.5, 1, alpha = 0.5))

# INT max
par(cex = 1.6, mar = c(4,4,1,1))
boxplot(results_df$INT_max ~ results_df$A * results_df$B, ylab="Intensidad de emisión máxima / (CPS)", xlab = "Año de circulación y posición", main = NULL, col = rgb(0.5, 0.5, 1, alpha = 0.5))

# INT ponderada
par(cex = 1.6, mar = c(4,4,1,1))
boxplot(results_df$INT_ponderada ~ results_df$A * results_df$B, ylab="Intensidad ponderada / (CPS)", xlab = "Año de circulación y posición", main = NULL, col = rgb(0.5, 0.5, 1, alpha = 0.5))


# Prueba de los supuestos del ANOVA
# Prueba de normalidad de los datos
shapiro.test(results_df$lambda_max) # p-value = 6.183e-13 (Los datos NO son normales)
shapiro.test(results_df$lambda_ponderada) # p-value = 4.756e-06 (Los datos NO son normales)
shapiro.test(results_df$INT_max) # p-value = 0.002533 (Los datos NO son normales)
shapiro.test(results_df$INT_ponderada) # p-value = 0.004931 (Los datos NO son normales)

# Prueba de Levene para igualdad de varianzas
leveneTest(lambda_max ~ A * B, data = results_df) # p-value = 0.08888 (Varianzas iguales)
leveneTest(lambda_ponderada ~ A * B, data = results_df) # p_value = 0.002728 (Varianzas diferentes)
leveneTest(INT_max ~ A * B, data = results_df) # p_value = 6.145e-05 (Varianzas diferentes)
leveneTest(INT_ponderada ~ A * B, data = results_df) # p_value = 0.2009 (Varianzas iguales)

# Dado que no se cumple con los supuestos de normalidad e igualdad de varianza, se emplea estadística no paramétrica
# Aligned Ranks Transformation ANOVA (ART ANOVA)

# Se definen los factores del experimento, de forma que se vean como variables
results_df$A <- factor(results_df$A) #Año de circulación
A <- factor(results_df$A)
results_df$B <- factor(results_df$B) #Posición
B <- factor(results_df$B)
results_df$Bloques <- factor(results_df$Bloques) # Bloque
Bloques <- factor(results_df$Bloques)
results_df$Muestras <- factor(results_df$Muestras) # Muestras
Muestras <- factor(results_df$Muestras)

# Se realiza un modelo lineal donde las variables "lambda_prom" y "INT_pond_prom" dependen de los factores A, B y su interacción
lambda_max.mod = art(lambda_max ~ (A * B + (1 | Bloques) + (1 | Bloques:Muestras)), data = results_df)
lambda_ponderada.mod = art(lambda_ponderada ~ (A * B + (1 | Bloques) + (1 | Bloques:Muestras)), data = results_df)
INT_max.mod = art(INT_max ~ (A * B + (1 | Bloques) + (1 | Bloques:Muestras)), data = results_df)
INT_ponderada.mod = art(INT_ponderada ~ (A * B + (1 | Bloques) + (1 | Bloques:Muestras)), data = results_df)

# Se verifica que el procedimiento ART es adecuado
summary(lambda_max.mod)
summary(lambda_ponderada.mod)
summary(INT_max.mod)
summary(INT_ponderada.mod)

# Se realiza la tabla ANOVA para cada variable de respuesta y se determina el eta parcial cuadrado
ART_lambda_max = anova(lambda_max.mod)
ART_lambda_max$part.eta.sq = with(anova(lambda_max.mod), `F` * `Df` / (`F` * `Df` + `Df.res`))
ART_lambda_max

ART_lambda_pond = anova(lambda_ponderada.mod)
ART_lambda_pond$part.eta.sq = with(anova(lambda_ponderada.mod), `F` * `Df` / (`F` * `Df` + `Df.res`))
ART_lambda_pond

ART_INT_max = anova(INT_max.mod)
ART_INT_max$part.eta.sq = with(anova(INT_max.mod), `F` * `Df` / (`F` * `Df` + `Df.res`))
ART_INT_max

ART_INT_pond = anova(INT_ponderada.mod)
ART_INT_pond$part.eta.sq = with(anova(INT_ponderada.mod), `F` * `Df` / (`F` * `Df` + `Df.res`))
ART_INT_pond

# Comparación entre pares para los factores significativos
# Lambda max
art.con(lambda_max.mod, ~ B)
cldList(p.value ~ contrast, data = as.data.frame(art.con(lambda_max.mod, ~ B)))
# Lambda ponderada
art.con(lambda_ponderada.mod, ~ B)
cldList(p.value ~ contrast, data = as.data.frame(art.con(lambda_ponderada.mod, ~ B)))
art.con(lambda_max.mod, ~ A*B)
cldList(p.value ~ contrast, data = as.data.frame(art.con(lambda_ponderada.mod, ~ A*B)))
# INT max
art.con(INT_max.mod, ~ B)
cldList(p.value ~ contrast, data = as.data.frame(art.con(INT_max.mod, ~ B)))
# INT ponderada
art.con(INT_ponderada.mod, ~ B)
cldList(p.value ~ contrast, data = as.data.frame(art.con(INT_ponderada.mod, ~ B)))
art.con(INT_ponderada.mod, ~ A*B)
cldList(p.value ~ contrast, data = as.data.frame(art.con(INT_ponderada.mod, ~ A*B)))

dev.off()
# Se realizan las gráficas de las interacciones
par(mfrow = c(2, 2), oma = c(4, 1, 0, 0), mar = c(4, 4.5, 0, 1), cex = 1.15, cex = 1.15, cex.lab=1.3)
interaction.plot(results_df$B,results_df$A,results_df$lambda_max,
                 legend = FALSE,
                 xlab = "",        
                 ylab = expression(paste(lambda[max], " / (nm)")),         
                 lty = 1:3,
                 col = c("blue", "red", "black"),
                 lwd = 2)
interaction.plot(results_df$B,results_df$A,results_df$lambda_ponderada,
                 legend = FALSE,
                 xlab = "",        
                 ylab = expression(paste(lambda[ponderada], " / (nm)")),         
                 lty = 1:3,
                 col = c("blue", "red", "black"),
                 lwd = 2)
interaction.plot(results_df$B,results_df$A,results_df$INT_max,
                 legend = FALSE,
                 xlab = "Posición",        
                 ylab = expression(paste(I[max], " / (CPS)")),         
                 lty = 1:3,
                 col = c("blue", "red", "black"),
                 lwd = 2)
interaction.plot(results_df$B,results_df$A,results_df$INT_ponderada,
                 legend = FALSE,
                 xlab = "Posición",        
                 ylab = expression(paste(I[ponderada], " / (CPS)")),         
                 lty = 1:3,
                 col = c("blue", "red", "black"),
                 lwd = 2)
par(xpd = NA)
legend("bottomleft", inset = c(-0.9,-0.7), title = "Año de circulación", 
       legend = c("2020", "2021", "2022"),
       col = c("blue", "red", "black"),
       lty = 1:3, horiz = TRUE, cex = 1)



### DISEÑO EXPERIMENTAL EN BILLETES DE DIEZ MIL COLONES ###

#BILLETES DE 10 MIL: POSICIÓN DE LA MARCA FLUORESCENTE
dev.off()
par(mfrow = c(2, 2), cex.axis = 1.2, cex.lab = 1.2, cex.main = 1.2)

pos_X_esq_inf_izq = data.frame(Numeración =c("Baja","Baja","Baja","Baja","Baja",
                                             "Media","Media","Media","Media","Media",
                                             "Alta","Alta","Alta","Alta","Alta"))
pos_X_esq_inf_izq$Numeración = factor(pos_X_esq_inf_izq$Numeración, levels = c("Baja", "Media", "Alta"))
pos_X_esq_inf_izq$Posición = c(74.67,
                               73.00,
                               73.69,
                               73.79,
                               72.84,
                               73.29,
                               72.81,
                               74.53,
                               73.54,
                               73.58,
                               74.52,
                               75.45,
                               75.03,
                               75.22,
                               72.30)
par(mar = c(2.5, 5, 2, 0))
boxplot(data = pos_X_esq_inf_izq, Posición ~ Numeración, ylab="Posición en X / (mm)", main = "Esquina inferior izquierda", col = rgb(0.5, 1, 0.5, alpha = 0.5))

pos_X_esq_inf_der = data.frame(Numeración =c("Baja","Baja","Baja","Baja","Baja",
                                             "Media","Media","Media","Media","Media",
                                             "Alta","Alta","Alta","Alta","Alta"))
pos_X_esq_inf_der$Numeración = factor(pos_X_esq_inf_der$Numeración, levels = c("Baja", "Media", "Alta"))
pos_X_esq_inf_der$Posición = c(90.11,
                               87.95,
                               88.96,
                               88.70,
                               87.42,
                               88.04,
                               87.36,
                               90.23,
                               88.62,
                               88.91,
                               90.13,
                               91.68,
                               91.01,
                               91.04,
                               86.86)
par(mar = c(2.5, 3, 2, 2))
boxplot(data = pos_X_esq_inf_der, Posición ~ Numeración, ylab="Posición en X / (mm)", main = "Esquina inferior derecha", col = rgb(0.5, 1, 0.5, alpha = 0.5))


pos_Y_esq_inf_izq = data.frame(Numeración =c("Baja","Baja","Baja","Baja","Baja",
                                             "Media","Media","Media","Media","Media",
                                             "Alta","Alta","Alta","Alta","Alta"))
pos_Y_esq_inf_izq$Numeración = factor(pos_Y_esq_inf_izq$Numeración, levels = c("Baja", "Media", "Alta"))
pos_Y_esq_inf_izq$Posición = c(26.96,
                               27.55,
                               27.46,
                               27.00,
                               28.38,
                               28.79,
                               28.34,
                               27.72,
                               27.90,
                               27.48,
                               28.13,
                               26.95,
                               27.08,
                               27.47,
                               28.80)
par(mar = c(4, 5, 0.5, 0))
boxplot(data = pos_Y_esq_inf_izq, Posición ~ Numeración, ylab="Posición en Y / (mm)", col = rgb(0.5, 1, 0.5, alpha = 0.5))

pos_Y_esq_inf_der = data.frame(Numeración =c("Baja","Baja","Baja","Baja","Baja",
                                             "Media","Media","Media","Media","Media",
                                             "Alta","Alta","Alta","Alta","Alta"))
pos_Y_esq_inf_der$Numeración = factor(pos_Y_esq_inf_der$Numeración, levels = c("Baja", "Media", "Alta"))
pos_Y_esq_inf_der$Posición = c(26.83,
                               27.29,
                               27.72,
                               26.96,
                               28.30,
                               29.00,
                               28.36,
                               27.64,
                               27.73,
                               27.47,
                               28.52,
                               26.86,
                               26.98,
                               27.51,
                               28.78)
par(mar = c(4, 3, 0.5, 2))
boxplot(data = pos_Y_esq_inf_der, Posición ~ Numeración, ylab="Posición en Y / (mm)", col = rgb(0.5, 1, 0.5, alpha = 0.5))


# Verificación de los supuestos del ANOVA de un factor

# Prueba de normalidad
shapiro.test(pos_X_esq_inf_izq$Posición) # p-value = 0.6155 (Los datos son normales)
shapiro.test(pos_Y_esq_inf_izq$Posición) # p-value = 0.1795 (Los datos son normales)
shapiro.test(pos_X_esq_inf_der$Posición) # p-value = 0.5658 (Los datos son normales)
shapiro.test(pos_Y_esq_inf_der$Posición) # p-value = 0.2751 (Los datos son normales)

# Prueba de igualdad de varianzas
leveneTest(pos_X_esq_inf_izq$Posición, pos_X_esq_inf_izq$Numeración, center = "mean") # p-value = 0.4362 (Varianzas iguales)
leveneTest(pos_Y_esq_inf_izq$Posición, pos_Y_esq_inf_izq$Numeración, center = "mean") # p-value = 0.4771 (Varianzas iguales)
leveneTest(pos_X_esq_inf_der$Posición, pos_X_esq_inf_der$Numeración, center = "mean") # p-value = 0.5153 (Varianzas iguales)
leveneTest(pos_Y_esq_inf_der$Posición, pos_Y_esq_inf_der$Numeración, center = "mean") # p-value = 0.3332 (Varianzas iguales)

# ANOVA de un factor (Factor temporal)
pos_X_esq_inf_izq.lm <- lm(Posición ~ Numeración, data = pos_X_esq_inf_izq)
anova(pos_X_esq_inf_izq.lm) # valor-p = 0.2253 (Medias iguales)

pos_Y_esq_inf_izq.lm <- lm(Posición ~ Numeración, data = pos_Y_esq_inf_izq)
anova(pos_Y_esq_inf_izq.lm) # valor-p = 0.3771 (Medias iguales)

pos_X_esq_inf_der.lm <- lm(Posición ~ Numeración, data = pos_X_esq_inf_der)
anova(pos_X_esq_inf_der.lm) # valor-p = 0.185 (Medias iguales)

pos_Y_esq_inf_der.lm <- lm(Posición ~ Numeración, data = pos_Y_esq_inf_der)
anova(pos_Y_esq_inf_der.lm) # valor-p = 0.4179 (Medias iguales)

# Análisis de los residuos

# Normalidad de los residuos
shapiro.test(pos_X_esq_inf_izq.lm$residuals) # valor-p = 0.1442 (Residuos normales)
shapiro.test(pos_Y_esq_inf_izq.lm$residuals) # valor-p = 0.2401 (Residuos normales)
shapiro.test(pos_X_esq_inf_der.lm$residuals) # valor-p = 0.1847 (Residuos normales)
shapiro.test(pos_Y_esq_inf_der.lm$residuals) # valor-p = 0.09644 (Residuos normales)

# Varianza constante
leveneTest(pos_X_esq_inf_izq.lm$residuals, pos_X_esq_inf_izq$Numeración, center = "mean") # valor-p = 0.4362 (Varianza constante)
leveneTest(pos_Y_esq_inf_izq.lm$residuals, pos_Y_esq_inf_izq$Numeración, center = "mean") # valor-p = 0.4771 (Varianza constante)
leveneTest(pos_X_esq_inf_der.lm$residuals, pos_X_esq_inf_der$Numeración, center = "mean") # valor-p = 0.5153 (Varianza constante)
leveneTest(pos_Y_esq_inf_der.lm$residuals, pos_Y_esq_inf_der$Numeración, center = "mean") # valor-p = 0.3332 (Varianza constante)

# Independencia de los residuos
bgtest(pos_X_esq_inf_izq.lm) # valor-p = 0.2698 (Residuos independientes)
bgtest(pos_Y_esq_inf_izq.lm) # valor-p = 0.798 (Residuos independientes)
bgtest(pos_X_esq_inf_der.lm) # valor-p = 0.4769 (Residuos independientes)
bgtest(pos_Y_esq_inf_der.lm) # valor-p = 0.8325 (Residuos independientes)

dev.off()
# Verificación gráfica
par(mfrow = c(2, 2))
plot(pos_X_esq_inf_izq.lm)
plot(pos_Y_esq_inf_izq.lm)
plot(pos_X_esq_inf_der.lm)
plot(pos_Y_esq_inf_der.lm)

# IMPORTACIÓN DE DATOS #

# Las siguientes líneas leen los archivos .CSV y los convierte en data frames (df)

# Se define el directorio de trabajo
setwd("C:\\Users\\james\\Google Drive\\TEC\\2024\\II SEMESTRE 2024\\TFG\\Espectros\\Espectros Art2\\DOE 10 MIL")


# Se lee un archivo .csv con los códigos empleados para guardar cada espectro de emisión medido
# Dichos códigos fueron generados a partir del diseño factorial aleatorizado, el cual fue obtenido de Minitab
codes <- read.csv("Codigos_DOE.csv", header = TRUE, colClasses = c("Codigo.DOE" = "character"))
code_names = c()
variable_names = c()
rep <- c(1,2,3)
for (i in seq_along(codes)) {
  for (j in seq(1,3)){
    code_names <- c(code_names, paste(codes[[i]], "-", rep[j], sep = ""))
    variable_names <- c(variable_names, paste("e", codes[[i]], "_", rep[j], sep = ""))
  }
}

# Se leen y se guardan los espectros de acuerdo a los códigos generados anteriormente
for (i in seq_along(code_names)) {
  file_name <- paste(code_names[[i]], ".csv", sep = "")
  assign(variable_names[[i]], read.csv(file_name, skip = 1, header = TRUE))
}

# ANALISIS DE LOS ESPECTROS #

# Los resultados de los espectros se guardan en un data frame para su análisis
results_df <- data.frame(code_name = character(),
                         OrdenEst = numeric(),
                         OrdenCorrida = numeric(),
                         TipoPt = numeric(),
                         Bloques = numeric(),
                         Rep = numeric(),
                         lambda_max = numeric(),
                         lambda_ponderada = numeric(),
                         INT_max = numeric(),
                         INT_ponderada = numeric(),
                         stringsAsFactors = FALSE)

for (i in seq_along(variable_names)) {
  # Se obtiene el dataframe usando get()
  data_to_analyze <- get(variable_names[[i]])
  
  # Se determina los picos locales de máxima emisión
  peaks <- localpeaks(data_to_analyze)
  
  # Se promedian integralmente la intensidad medida
  mean_intensity <- mean_spectrum_INT(data_to_analyze)
  
  # Se promedian ponderadamente las longitudes de onda a partir de la intensidad
  centroid_wavelength <- centroid_wavelength_nm(data_to_analyze)
  
  # Se extraen los números de los códigos del DOE
  parts <- strsplit(code_names[[i]], "-")[[1]]
  numbers <- strsplit(parts[1], "")[[1]]
  
  # Se crea una nueva fila con los resultados de los picos en cada espectro
  new_row <- data.frame(
    code_name = code_names[[i]],
    OrdenEst = paste(numbers[1:2], collapse = ""),
    OrdenCorrida = paste(numbers[3:4], collapse = ""),
    Bloques = paste(numbers[5:5], collapse = ""),
    Rep = paste(parts[2:2], collapse = ""),
    lambda_max = peaks[1, "nm"],
    lambda_ponderada = centroid_wavelength,
    INT_max = peaks[1, "INT"],
    INT_ponderada = mean_intensity
  )
  
  # Se añade la nueva fila al dataframe de resultados
  results_df <- rbind(results_df, new_row)
}

# Se ordena por code_name
results_df <- results_df[order(results_df$code_name, decreasing = FALSE),]

# Se añaden las variables de los factores al dataframe de resultado
A <- c("Baja", "Media", "Alta") # Numeración de la serie
B <- c("A", "B", "C") # Posición
C <- c("1", "2", "3") # Orden de medición

# Se crea una combinación de los factores
combinations <- expand.grid(C= C, B = B, A = A)

# Se replican las combinaciones 5 veces (debido al factor bloque)
factors_df <- combinations[rep(seq_len(nrow(combinations)), times = 5), ]

# Se añaden las columnas al results_df
results_df$A = factors_df$A
results_df$B = factors_df$B

# Se añade un factor de muestras (1:15)
results_df$Muestras = rep(1:15, each = 9)

# Se guarda el results_df como un archivo ".csv"
write.csv(results_df, "Results_df.csv", row.names = FALSE)


# ANALISIS DISEÑO EXPERIMENTAL # alpha = 0.05 (95% de confianza)


# Histograma con gráfico de caja para cada variable de respuesta (Guardar figuras a 550x550)

dev.off()
# Lambda max
layout(matrix(1, nrow = 1))
par(mar = c(4, 4, 1, 1), cex = 1.2) 
# Gráfico principal: Histograma
hist(results_df$lambda_max, main = NULL, xlab = "Longitud de onda de máxima emisión / (nm)", 
     xlim = c(450, 675), ylab = "Densidad", ylim = c(-0.0023, 0.05), freq = FALSE)
lines(density(results_df$lambda_max, bw = "nrd"), col = "darkgreen", lwd = 2)
box()
# Gráfico secundario: diagrama de caja
par(new = TRUE, mar = c(4, 4, 18.5, 1))
boxplot(results_df$lambda_max, horizontal = TRUE, axes = FALSE, 
        ylim = c(475, 675), col = rgb(0.5, 1, 0.5, alpha = 0.5), at = -0.002)

dev.off()
# Lambda ponderada
layout(matrix(1, nrow = 1))
par(mar = c(4, 4, 1, 1), cex = 1.2) 
# Gráfico principal: Histograma
hist(results_df$lambda_ponderada, main = NULL, xlab = "Longitud de onda ponderada / (nm)", 
     xlim = c(530, 560), ylab = "Densidad", ylim = c(-0.013, 0.2), freq = FALSE, breaks = 12)
lines(density(results_df$lambda_ponderada, bw = "nrd"), col = "darkgreen", lwd = 2)
box()
# Gráfico secundario: diagrama de caja
par(new = TRUE, mar = c(4, 4, 18.5, 1))
boxplot(results_df$lambda_ponderada, horizontal = TRUE, axes = FALSE, 
        ylim = c(530, 560), col = rgb(0.5, 1, 0.5, alpha = 0.5), at = -0.002)

dev.off()
# INT max
layout(matrix(1, nrow = 1))
par(mar = c(4, 4, 1, 1), cex = 1.2) 
# Gráfico principal: Histograma
hist(results_df$INT_max, main = NULL, xlab = "Intensidad de emisión máxima / (CPS)",
     xlim = c(0,35000), ylab = "Densidad", ylim = c(-7e-06, 0.0001), freq= FALSE)
lines(density(results_df$INT_max, bw = "nrd"), col = "darkgreen", lwd = 2)
box()
# Gráfico secundario: diagrama de caja
par(new = TRUE, mar = c(4, 4, 18.5, 1))
boxplot(results_df$INT_max, horizontal = TRUE, axes = FALSE, 
        ylim = c(0,35000), col = rgb(0.5, 1, 0.5, alpha = 0.5), at = -0.002)

dev.off()
# INT ponderada
layout(matrix(1, nrow = 1))
par(mar = c(4, 4, 1, 1), cex = 1.2) 
# Gráfico principal: Histograma
hist(results_df$INT_ponderada, main = NULL, xlab = "Intensidad de emisión ponderada / (CPS)",
     xlim = c(0,12000), ylab = "Densidad", ylim = c(-0.00002, 0.0003), freq= FALSE)
lines(density(results_df$INT_ponderada, bw = "nrd"), col = "darkgreen", lwd = 2)
box()
# Gráfico secundario: diagrama de caja
par(new = TRUE, mar = c(4, 4, 18.5, 1))
boxplot(results_df$INT_ponderada, horizontal = TRUE, axes = FALSE, 
        ylim = c(0,12000), col = rgb(0.5, 1, 0.5, alpha = 0.5), at = -0.002)

dev.off()
# Gráficos de caja (Guardar 1200x675)

# Factores
results_df$A = factor(results_df$A, levels = c("Baja", "Media", "Alta"))
results_df$B = factor(results_df$B, levels = c("A", "B", "C"))

# Lambda max
par(cex = 1.55, mar = c(4,4,1,1))
boxplot(results_df$lambda_max ~ results_df$A * results_df$B, ylab="Longitud de onda de máxima emisión / (nm)", xlab = "Numeración y posición", main = NULL, col = rgb(0.5, 1, 0.5, alpha = 0.5))

# Lambda ponderada
par(cex = 1.55, mar = c(4,4,1,1))
boxplot(results_df$lambda_ponderada ~ results_df$A * results_df$B, ylab="Longitud de onda ponderada / (nm)", xlab = "Numeración y posición", main = NULL, col = rgb(0.5, 1, 0.5, alpha = 0.5))

# INT max
par(cex = 1.55, mar = c(4,4,1,1))
boxplot(results_df$INT_max ~ results_df$A * results_df$B, ylab="Intensidad de emisión máxima / (CPS)", xlab = "Numeración y posición", main = NULL, col = rgb(0.5, 1, 0.5, alpha = 0.5))

# INT ponderada
par(cex = 1.55, mar = c(4,4,1,1))
boxplot(results_df$INT_ponderada ~ results_df$A * results_df$B, ylab="Intensidad ponderada / (CPS)", xlab = "Numeración y posición", main = NULL, col = rgb(0.5, 1, 0.5, alpha = 0.5))


# Prueba de los supuestos del ANOVA
# Prueba de normalidad de los datos
shapiro.test(results_df$lambda_max) # p-value = 1.671e-15 (Los datos NO son normales)
shapiro.test(results_df$lambda_ponderada) # p-value = 4.518e-05 (Los datos NO son normales)
shapiro.test(results_df$INT_max) # p-value = 0.01616 (Los datos NO son normales)
shapiro.test(results_df$INT_ponderada) # p-value = 0.001185 (Los datos NO son normales)

# Prueba de Levene para igualdad de varianzas
leveneTest(lambda_max ~ A * B, data = results_df) # p-value = 0.2465 (Varianzas iguales)
leveneTest(lambda_ponderada ~ A * B, data = results_df) # p_value = 0.001765 (Varianzas diferentes)
leveneTest(INT_max ~ A * B, data = results_df) # p_value = 0.02882 (Varianzas diferentes)
leveneTest(INT_ponderada ~ A * B, data = results_df) # p_value = 0.1723 (Varianzas iguales)

# Dado que no se cumple con los supuestos de normalidad e igualdad de varianza (en la mitad de los casos), se emplea estadística no paramétrica
# Aligned Ranks Transformation ANOVA (ART ANOVA)

# Se definen los factores del experimento, de forma que se vean como variables
results_df$A <- factor(results_df$A) #Año de circulación
A <- factor(results_df$A)
results_df$B <- factor(results_df$B) #Posición
B <- factor(results_df$B)
results_df$Bloques <- factor(results_df$Bloques) # Bloque
Bloques <- factor(results_df$Bloques)
results_df$Muestras <- factor(results_df$Muestras) # Muestras
Muestras <- factor(results_df$Muestras)

# Se realiza un modelo lineal donde las variables "lambda_prom" y "INT_pond_prom" dependen de los factores A, B y su interacción
lambda_max.mod = art(lambda_max ~ (A * B + (1 | Bloques) + (1 |Bloques:Muestras)), data = results_df)
lambda_ponderada.mod = art(lambda_ponderada ~ (A * B + (1 | Bloques) + (1 | Bloques:Muestras)), data = results_df)
INT_max.mod = art(INT_max ~ (A * B + (1 | Bloques) + (1 | Bloques:Muestras)), data = results_df)
INT_ponderada.mod = art(INT_ponderada ~ (A * B + (1 | Bloques) + (1 | Bloques:Muestras)), data = results_df)

# Se verifica que el procedimiento ART es adecuado
summary(lambda_max.mod)
summary(lambda_ponderada.mod)
summary(INT_max.mod)
summary(INT_ponderada.mod)

# Se realiza la tabla ANOVA para cada variable de respuesta y se determina el eta parcial cuadrado
ART_lambda_max = anova(lambda_max.mod)
ART_lambda_max$part.eta.sq = with(anova(lambda_max.mod), `F` * `Df` / (`F` * `Df` + `Df.res`))
ART_lambda_max

ART_lambda_pond = anova(lambda_ponderada.mod)
ART_lambda_pond$part.eta.sq = with(anova(lambda_ponderada.mod), `F` * `Df` / (`F` * `Df` + `Df.res`))
ART_lambda_pond

ART_INT_max = anova(INT_max.mod)
ART_INT_max$part.eta.sq = with(anova(INT_max.mod), `F` * `Df` / (`F` * `Df` + `Df.res`))
ART_INT_max

ART_INT_pond = anova(INT_ponderada.mod)
ART_INT_pond$part.eta.sq = with(anova(INT_ponderada.mod), `F` * `Df` / (`F` * `Df` + `Df.res`))
ART_INT_pond

# Comparación entre pares para los factores significativos
# Lambda max
art.con(lambda_max.mod, ~ B)
cldList(p.value ~ contrast, data = as.data.frame(art.con(lambda_max.mod, ~ B)))
# Lambda ponderada
art.con(lambda_ponderada.mod, ~ B)
cldList(p.value ~ contrast, data = as.data.frame(art.con(lambda_ponderada.mod, ~ B)))
# INT max
art.con(INT_max.mod, ~ B)
cldList(p.value ~ contrast, data = as.data.frame(art.con(INT_max.mod, ~ B)))
# INT ponderada
art.con(INT_ponderada.mod, ~ B)
cldList(p.value ~ contrast, data = as.data.frame(art.con(INT_ponderada.mod, ~ B)))
art.con(INT_ponderada.mod, ~ A*B)
cldList(p.value ~ contrast, data = as.data.frame(art.con(INT_ponderada.mod, ~ A*B)))

dev.off()
# Se realizan las gráficas de las interacciones
par(mfrow = c(2, 2), oma = c(4, 1, 0, 0), mar = c(4, 4.5, 0, 1), cex = 1.15, cex.lab=1.3)
interaction.plot(results_df$B,results_df$A,results_df$lambda_max,
                 legend = FALSE,
                 xlab = "",        
                 ylab = expression(paste(lambda[max], " / (nm)")),
                 lty = 1:3,
                 col = c("blue", "red", "black"),
                 lwd = 2)
interaction.plot(results_df$B,results_df$A,results_df$lambda_ponderada,
                 legend = FALSE,
                 xlab = "",        
                 ylab = expression(paste(lambda[ponderada], " / (nm)")),         
                 lty = 1:3,
                 col = c("blue", "red", "black"),
                 lwd = 2)
interaction.plot(results_df$B,results_df$A,results_df$INT_max,
                 legend = FALSE,
                 xlab = "Posición",        
                 ylab = expression(paste(I[max], " / (CPS)")),         
                 lty = 1:3,
                 col = c("blue", "red", "black"),
                 lwd = 2)
interaction.plot(results_df$B,results_df$A,results_df$INT_ponderada,
                 legend = FALSE,
                 xlab = "Posición",        
                 ylab = expression(paste(I[ponderada], " / (CPS)")),         
                 lty = 1:3,
                 col = c("blue", "red", "black"),
                 lwd = 2)
par(xpd = NA)
legend("bottomleft", inset = c(-0.9,-0.7), title = "Numeración", 
       legend = c("Baja", "Media", "Alta"),
       col = c("blue", "red", "black"),
       lty = 1:3, horiz = TRUE, cex = 1)
