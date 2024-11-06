# Autor: James Andrés Ruiz Vásquez
# Código empleado para el desarrollo del artículo:
# "Estudios espectroscópicos en billetes de Costa Rica actualmente en circulación"
# Versión: 1.0
# Fecha: 05/11/2024

#Se importan las bibliotecas necesarias
library(tidyverse) # Conjunto de paquetes para manipulación de datos
library(ggpmisc)
library(photobiology) # Manipulación de espectros
library(plotly) # Figuras 3D
library(cooltools) # Wavelength to RGB


### ESPECTROS DE EMISIÓN/EXCITACIÓN UV-VIS ###


# Para espectros en 3D, las columnas están categorizadas de la siguiente forma:
# V1 -> EM Wavelength / nm [Longitud de onda de emisión / (nm)]
# V2 -> EX Wavelength / nm [Longitud de onda de excitación / (nm)]
# V3 -> S1c / R1c / CPS/μA [Intensidad / (CPS/μA) ]
# NOTA: Para hacer las figuras 3D, se debe usar los archivos "Contour" no los "Graph"

# Para espectros en 2D, las columnas están categorizadas de la siguiente forma:
# V1 -> Wavelength / nm [Longitud de onda (EM ó EX) / (nm)]
# V2 -> S1c / R1c / CPS/μA [Intensidad / (CPS/μA)]
# V3 -> Wavelength / nm [Longitud de onda (EM ó EX) / (nm)]
# V4 -> S1 / CPS [Intensidad / (CPS)]


# IMPORTACIÓN DE DATOS #

# Las siguientes líneas leen los archivos .dat y los convierte en data frames (df)


# Figuras 3D en dos posiciones diferentes: A Y B

# Billete de 1 Mil
Bill1mil_3D_PosA_Cont <- read.table("C:\\Users\\james\\Google Drive\\TEC\\2024\\II SEMESTRE 2024\\TFG\\Espectros\\Espectros Art1\\ESPECTROS EM EX\\Billete 1 Mil\\3D\\PosA\\Bill1Mil 3D (01)_Contour_S1c_R1c.dat", header = FALSE, skip = 2)
Bill1mil_3D_PosB_Cont <- read.table("C:\\Users\\james\\Google Drive\\TEC\\2024\\II SEMESTRE 2024\\TFG\\Espectros\\Espectros Art1\\ESPECTROS EM EX\\Billete 1 Mil\\3D\\PosB\\Bill1Mil 1 (02)_Contour_S1c_R1c.dat", header = FALSE, skip = 2)

# Billete de 2 Mil
Bill2mil_3D_PosA_Cont <- read.table("C:\\Users\\james\\Google Drive\\TEC\\2024\\II SEMESTRE 2024\\TFG\\Espectros\\Espectros Art1\\ESPECTROS EM EX\\Billete 2 Mil\\3D\\PosA\\Bill2Mil 3D A (01)_Contour_S1c_R1c.dat", header  <-  FALSE, skip = 2)
Bill2mil_3D_PosB_Cont <- read.table("C:\\Users\\james\\Google Drive\\TEC\\2024\\II SEMESTRE 2024\\TFG\\Espectros\\Espectros Art1\\ESPECTROS EM EX\\Billete 2 Mil\\3D\\PosB\\Bill2Mil 3D B (01)_Contour_S1c_R1c.dat", header = FALSE, skip = 2)

# Billete de 5 Mil
Bill5mil_3D_PosA_Cont <- read.table("C:\\Users\\james\\Google Drive\\TEC\\2024\\II SEMESTRE 2024\\TFG\\Espectros\\Espectros Art1\\ESPECTROS EM EX\\Billete 5 Mil\\3D\\PosA\\Bill5Mil 3D A (01)_Contour_S1c_R1c.dat", header = FALSE, skip = 2)
Bill5mil_3D_PosB_Cont <- read.table("C:\\Users\\james\\Google Drive\\TEC\\2024\\II SEMESTRE 2024\\TFG\\Espectros\\Espectros Art1\\ESPECTROS EM EX\\Billete 5 Mil\\3D\\PosB\\Bill5Mil 3D B (01)_Contour_S1c_R1c.dat", header = FALSE, skip = 2)

# Billete de 10 Mil
Bill10mil_3D_PosA_Cont <- read.table("C:\\Users\\james\\Google Drive\\TEC\\2024\\II SEMESTRE 2024\\TFG\\Espectros\\Espectros Art1\\ESPECTROS EM EX\\Billete 10 Mil\\3D\\PosA\\Billete10mil (01)_Contour_S1c_R1c.dat", header = FALSE, skip = 2)
Bill10mil_3D_PosB_Cont <- read.table("C:\\Users\\james\\Google Drive\\TEC\\2024\\II SEMESTRE 2024\\TFG\\Espectros\\Espectros Art1\\ESPECTROS EM EX\\Billete 10 Mil\\3D\\PosB\\Billete10mil (01)_Contour_S1c_R1c.dat", header = FALSE, skip = 2)
# Filtrado de datos para emisiones a partir de 400 nm
Bill10mil_3D_PosA_Cont <- Bill10mil_3D_PosA_Cont %>% filter(V1 >= 400)
Bill10mil_3D_PosB_Cont <- Bill10mil_3D_PosB_Cont %>% filter(V1 >= 400)

# Billete de 20 Mil
Bill20mil_3D_PosA_Cont <- read.table("C:\\Users\\james\\Google Drive\\TEC\\2024\\II SEMESTRE 2024\\TFG\\Espectros\\Espectros Art1\\ESPECTROS EM EX\\Billete 20 Mil\\3D\\PosA\\Billete20mil (01)_Contour_S1c_R1c.dat", header = FALSE, skip = 2)
Bill20mil_3D_PosB_Cont <- read.table("C:\\Users\\james\\Google Drive\\TEC\\2024\\II SEMESTRE 2024\\TFG\\Espectros\\Espectros Art1\\ESPECTROS EM EX\\Billete 20 Mil\\3D\\PosB\\Billete20mil (01)_Contour_S1c_R1c.dat", header = FALSE, skip = 2)
# Filtrado de datos para emisiones a partir de 400 nm
Bill20mil_3D_PosA_Cont <- Bill20mil_3D_PosA_Cont %>% filter(V1 >= 400)
Bill20mil_3D_PosB_Cont <- Bill20mil_3D_PosB_Cont %>% filter(V1 >= 400)


# Espectros de Emisión con excitación a 365nm para cada denominación en dos posiciones diferentes: A Y B

# Billete de 1 Mil
Bill1milEmExc365nm_PosA <- read.table("C:\\Users\\james\\Google Drive\\TEC\\2024\\II SEMESTRE 2024\\TFG\\Espectros\\Espectros Art1\\ESPECTROS EM EX\\Billete 1 Mil\\EmExc365nm\\EmExc365nmA (01)_Graph.dat", header = FALSE, skip = 2)
Bill1milEmExc365nm_PosB <- read.table("C:\\Users\\james\\Google Drive\\TEC\\2024\\II SEMESTRE 2024\\TFG\\Espectros\\Espectros Art1\\ESPECTROS EM EX\\Billete 1 Mil\\EmExc365nm\\EmExc365nmB (01)_Graph.dat", header = FALSE, skip = 2)

# Billete de 2 Mil
Bill2milEmExc365nm_PosA <- read.table("C:\\Users\\james\\Google Drive\\TEC\\2024\\II SEMESTRE 2024\\TFG\\Espectros\\Espectros Art1\\ESPECTROS EM EX\\Billete 2 Mil\\EmExc365nm\\EmExc365nmA (01)_Graph.dat", header = FALSE, skip = 2)
Bill2milEmExc365nm_PosB <- read.table("C:\\Users\\james\\Google Drive\\TEC\\2024\\II SEMESTRE 2024\\TFG\\Espectros\\Espectros Art1\\ESPECTROS EM EX\\Billete 2 Mil\\EmExc365nm\\EmExc365nmB (01)_Graph.dat", header = FALSE, skip = 2)

# Billete de 5 Mil
Bill5milEmExc365nm_PosA <- read.table("C:\\Users\\james\\Google Drive\\TEC\\2024\\II SEMESTRE 2024\\TFG\\Espectros\\Espectros Art1\\ESPECTROS EM EX\\Billete 5 Mil\\EmExc365nm\\EmExc365nmA (01)_Graph.dat", header = FALSE, skip = 2)
Bill5milEmExc365nm_PosB <- read.table("C:\\Users\\james\\Google Drive\\TEC\\2024\\II SEMESTRE 2024\\TFG\\Espectros\\Espectros Art1\\ESPECTROS EM EX\\Billete 5 Mil\\EmExc365nm\\EmExc365nmB (01)_Graph.dat", header = FALSE, skip = 2)

# Billete de 10 Mil
Bill10milEmExc365nm_PosA <- Bill10mil_3D_PosA_Cont[Bill10mil_3D_PosA_Cont$V2 == 365, ]
Bill10milEmExc365nm_PosB <- Bill10mil_3D_PosB_Cont[Bill10mil_3D_PosB_Cont$V2 == 365, ]

# Billete de 20 Mil
Bill20milEmExc365nm_PosA <- Bill20mil_3D_PosA_Cont[Bill20mil_3D_PosA_Cont$V2 == 365, ]
Bill20milEmExc365nm_PosB <- Bill20mil_3D_PosB_Cont[Bill20mil_3D_PosB_Cont$V2 == 365, ]


# Espectros de Excitación con longitud de onda de emisión características, únicamente en posición A

# Billete de 1 Mil
Bill1milExcEm536nm <- read.table("C:\\Users\\james\\Google Drive\\TEC\\2024\\II SEMESTRE 2024\\TFG\\Espectros\\Espectros Art1\\ESPECTROS EM EX\\Billete 1 Mil\\Exc\\PosA\\ExcEm536nmA (01)_Graph.dat", header = FALSE, skip = 2)
Bill1milExcEm610nm <- read.table("C:\\Users\\james\\Google Drive\\TEC\\2024\\II SEMESTRE 2024\\TFG\\Espectros\\Espectros Art1\\ESPECTROS EM EX\\Billete 1 Mil\\Exc\\PosA\\ExcEm610nmA (01)_Graph.dat", header = FALSE, skip = 2)
# Filtrado de datos para excitaciones entre 320 nm y 400 nm
Bill1milExcEm536nm <- Bill1milExcEm536nm[Bill1milExcEm536nm$V1 >= 320,]
Bill1milExcEm610nm <- Bill1milExcEm610nm[Bill1milExcEm610nm$V1 >= 320,]

# Billete de 2 Mil
Bill2milExcEm510nm <- read.table("C:\\Users\\james\\Google Drive\\TEC\\2024\\II SEMESTRE 2024\\TFG\\Espectros\\Espectros Art1\\ESPECTROS EM EX\\Billete 2 Mil\\Exc\\PosA\\ExcEm510nmA (01)_Graph.dat", header = FALSE, skip = 2)
Bill2milExcEm610nm <- read.table("C:\\Users\\james\\Google Drive\\TEC\\2024\\II SEMESTRE 2024\\TFG\\Espectros\\Espectros Art1\\ESPECTROS EM EX\\Billete 2 Mil\\Exc\\PosA\\ExcEm610nmA (01)_Graph.dat", header = FALSE, skip = 2)
Bill2milExcEm510nm <- Bill2milExcEm510nm[Bill2milExcEm510nm$V1 >= 320,]
Bill2milExcEm610nm <- Bill2milExcEm610nm[Bill2milExcEm610nm$V1 >= 320,]

# Billete de 5 Mil
Bill5milExcEm510nm <- read.table("C:\\Users\\james\\Google Drive\\TEC\\2024\\II SEMESTRE 2024\\TFG\\Espectros\\Espectros Art1\\ESPECTROS EM EX\\Billete 5 Mil\\Exc\\PosA\\ExcEm510nmA (01)_Graph.dat", header = FALSE, skip = 2)
Bill5milExcEm610nm <- read.table("C:\\Users\\james\\Google Drive\\TEC\\2024\\II SEMESTRE 2024\\TFG\\Espectros\\Espectros Art1\\ESPECTROS EM EX\\Billete 5 Mil\\Exc\\PosA\\ExcEm610nmA (01)_Graph.dat", header = FALSE, skip = 2)
Bill5milExcEm510nm <- Bill5milExcEm510nm[Bill5milExcEm510nm$V1 >= 320,]
Bill5milExcEm610nm <- Bill5milExcEm610nm[Bill5milExcEm610nm$V1 >= 320,]

# Billete de 10 Mil
Bill10milExcEm510nm <- Bill10mil_3D_PosA_Cont[Bill10mil_3D_PosA_Cont$V1 == 510, ]
Bill10milExcEm610nm <- Bill10mil_3D_PosA_Cont[Bill10mil_3D_PosB_Cont$V1 == 610, ]

# Billete de 20 Mil
Bill20milExcEm510nm <- Bill20mil_3D_PosA_Cont[Bill20mil_3D_PosA_Cont$V1 == 510, ]
Bill20milExcEm610nm <- Bill20mil_3D_PosA_Cont[Bill20mil_3D_PosA_Cont$V1 == 610, ]


# FUNCIONES #

# Se realizan las gráficas con ggplot2, para ello se crea la siguiente función:

plot_3D <- function(dfCont) {
  # Descripción: Función para generar los espectros en 3D
  # Entradas:
  # dfCont = dataframe del contorno
  # Salidas:
  # Gráfico 3D de la intensidad en función de la longitud de onda de EM y EX.
  
  # Se escala la intensidad en unidades arbitrarias
  int_norm <- dfCont$V3 / max(dfCont$V3)
  dfCont$V3 <- int_norm
  
  # Se determinan los valores de frontera en "x" y en "y"
  x_min <- min(dfCont$V1)-2
  x_max <- max(dfCont$V1)+2
  y_min <- min(dfCont$V2)-3
  y_max <- max(dfCont$V2)+3
  
  # Se crean filas adicionales con NA en "z" en las fronteras
  na_rows <- data.frame(
    V1 = c(rep(x_min, 2), rep(x_max, 2)),
    V2 = c(y_min, y_max, y_min, y_max),
    V3 = NA
  )
  
  # Se unen las filas de NA al dataframe original
  df_with_na <- rbind(dfCont, na_rows)
  
  # Se genera la matriz de color VIS (400 nm a 700 nm)
  lambda = seq(400, 700, length.out = 100) # Rango visible
  colors <- sapply(lambda, wavelength2col)
  
  # Se crea un colorscale personalizado para plotly
  custom_colorscale <- lapply(seq_along(colors), function(i) {
    c((i - 1) / (length(colors) - 1), colors[i])
  })
  
  plot_ly(df_with_na, 
          x = ~df_with_na$V2, 
          y = ~df_with_na$V1, 
          z = ~df_with_na$V3, 
          type = "mesh3d", 
          intensity = ~df_with_na$V1,  # Color basado en los valores "x"
          colorscale = custom_colorscale,
          showscale = FALSE,
          opacity = 1,
          flatshading = TRUE
  ) %>%
    layout(scene = list(
      xaxis = list(title = "Longitud de onda excitación / (nm)"),
      yaxis = list(title = "Longitud de onda emisión / (nm)"),
      zaxis = list(title = "Intensidad / (u.a)")
    ))
}


plot_spectrum <- function(df, type) {
  # Descripción: Función para generar espectros de emisión o excitación en 2D
  # Entradas:
  # df = dataframe
  # type = tipo de gráfico
  #        0 = Espectro de emisión regular
  #        1 = Espectro de emisión a partir de datos medidos en espectro 3D
  #        2 = Espectro de excitación a partir de datos medidos en espectro 3D
  #        3 = Espectro de excitación regular
  # Salidas:
  # Gráfico del espectro de emisión o excitación en 2D
  
  # TIPO 0: Espectro de emisión regular
  if (type==0) {
    # Se escala la intensidad en unidades arbitrarias
    int_norm <- df$V2 / max(df$V2)
    df$V2 <- int_norm
    
    # La siguiente línea grafica el espectro
    graph <- ggplot(mapping = aes(x=df$V1, y=df$V2)) +
      geom_line(size = 0.75) +
      stat_peaks(col = "red", span=NULL) +
      labs(
        x = "Longitud de onda / (nm)",
        y = "Intensidad / (u.a)"
      ) +
      theme_bw() +
      theme(
        legend.text=element_text(size=12),
        legend.title = element_text(size = 12),
        legend.position = c(.9, .85),
        legend.background = element_rect(linetype = "solid", colour = "black"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12)
      ) +
      scale_x_continuous(breaks = seq(min(df$V1), max(df$V1), by = 50))
    # La siguiente línea determina los picos máximos locales
    localPeaks <- peaks(x = df, span = 7, x.var.name = "V1", y.var.name = "V2")
    locPeaks <- localPeaks[order(-localPeaks$V2),]
    print("Picos máximos locales (Resaltados en rojo)")
    print(locPeaks)
    absPeak <- df[which.max(df$V2),]
    print("Pico máximo global")
    print(absPeak)}
  
  # TIPO 1: Espectro de emisión modificado
  else if (type==1) {
    # Se escala la intensidad en unidades arbitrarias
    int_norm <- df$V3 / max(df$V3)
    df$V3 <- int_norm
    
    # La siguiente línea grafica el espectro
    graph <- ggplot(mapping = aes(x=df$V1, y=df$V3)) +
      geom_line(size = 0.75) +
      stat_peaks(col = "red", span=NULL) +
      labs(
        x = "Longitud de onda / (nm)",
        y = "Intensidad / (u.a)"
      ) +
      theme_bw() +
      theme(
        legend.text=element_text(size=12),
        legend.title = element_text(size = 12),
        legend.position = c(.9, .85),
        legend.background = element_rect(linetype = "solid", colour = "black"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12)
      ) +
      scale_x_continuous(breaks = seq(min(df$V1), max(df$V1), by = 50))
    # La siguiente línea determina los picos máximos locales
    localPeaks <- peaks(x = df, span = 7, x.var.name = "V1", y.var.name = "V3")
    locPeaks <- localPeaks[order(-localPeaks$V3),]
    print("Picos máximos locales (Resaltados en rojo)")
    print(locPeaks)
    absPeak <- df[which.max(df$V3),]
    print("Pico máximo global")
    print(absPeak)}
  
  # TIPO 2: Espectro de excitación modificado
  else if (type==2) {
    # Se escala la intensidad en unidades arbitrarias
    int_norm <- df$V3 / max(df$V3)
    df$V3 <- int_norm
    
    # La siguiente línea grafica el espectro
    graph <- ggplot(mapping = aes(x=df$V2, y=df$V3)) +
      geom_line(size = 0.75) +
      stat_peaks(col = "red", span=NULL) +
      labs(
        x = "Longitud de onda / (nm)",
        y = "Intensidad / (u.a)"
      ) +
      theme_bw() +
      theme(
        legend.text=element_text(size=12),
        legend.title = element_text(size = 12),
        legend.position = c(.9, .85),
        legend.background = element_rect(linetype = "solid", colour = "black"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12)
      ) +
      scale_x_continuous(breaks = seq(min(df$V2), max(df$V2), by = 5))
    # La siguiente línea determina los picos máximos locales
    localPeaks <- peaks(x = df, span = 7, x.var.name = "V2", y.var.name = "V3")
    locPeaks <- localPeaks[order(-localPeaks$V3),]
    print("Picos máximos locales (Resaltados en rojo)")
    print(locPeaks)
    absPeak <- df[which.max(df$V3),]
    print("Pico máximo global")
    print(absPeak)}
  
  # TIPO 3: Espectro de excitación regular
  else if (type==3) {
    # Se escala la intensidad en unidades arbitrarias
    int_norm <- df$V2 / max(df$V2)
    df$V2 <- int_norm
    
    # La siguiente línea grafica el espectro
    graph <- ggplot(mapping = aes(x=df$V1, y=df$V2)) +
      geom_line(size = 0.75) +
      stat_peaks(col = "red", span=NULL) +
      labs(
        x = "Longitud de onda / (nm)",
        y = "Intensidad / (u.a)"
      ) +
      theme_bw() +
      theme(
        legend.text=element_text(size=12),
        legend.title = element_text(size = 12),
        legend.position = c(.9, .85),
        legend.background = element_rect(linetype = "solid", colour = "black"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12)
      ) +
      scale_x_continuous(breaks = seq(min(df$V1), max(df$V1), by = 10))
    # La siguiente línea determina los picos máximos locales
    localPeaks <- peaks(x = df, span = 7, x.var.name = "V1", y.var.name = "V2")
    locPeaks <- localPeaks[order(-localPeaks$V2),]
    print("Picos máximos locales (Resaltados en rojo)")
    print(locPeaks)
    absPeak <- df[which.max(df$V2),]
    print("Pico máximo global")
    print(absPeak)
    
  }
  print(graph)}


plot_espectra <- function(l, type){
  # Descripción: Función para graficar varios espectros en conjunto en un mismo plot
  # Usados para comparar espectros de emisión en dos posiciones diferentes
  # Entradas:
  # l = lista de dataframes de los espectros a graficar
  # type = tipo de gráfico
  #        0 = Espectro de emisión regular
  #        1 = Espectro de emisión a partir de datos medidos en espectro 3D
  #        2 = Espectro de excitación 
  # Salidas:
  # Gráfico con los espectros superpuestos
  
  # Vector con los colores de las gráficas
  color_val <- c("red","blue", "green")
  
  # Se itera, de forma que se añada un identificador para cada dataframe
  for (i in seq_along(l)) {
    l[[i]]$id <- paste(LETTERS[i])  # Se agrega un identificador
  }
  
  # Se combina todos los dataframes en uno solo
  data_combinada <- do.call(rbind, l)
  
  # TIPO 0: Espectros de emisión regulares 
  if (type==0) {
    # Se itera, de forma que se añada un identificador para cada dataframe
    for (i in seq_along(l)) {
      l[[i]]$id <- paste(LETTERS[i])  # Se agrega un identificador
    }
    
    # Se combina todos los dataframes en uno solo
    data_combinada <- do.call(rbind, l)
    
    # Se escala la intensidad en unidades arbitrarias
    int_norm <- data_combinada$V2 / max(data_combinada$V2)
    data_combinada$V2 <- int_norm
    
    # Se crea la gráfica
    graph <- ggplot(data_combinada, aes(x = V1, y = V2, color = id)) +
      geom_line(size = 0.75) +
      scale_color_manual(values = color_val, name = "Posición" ) +
      labs(
        x = "Longitud de onda / (nm)",
        y = "Intensidad / (u.a)"
      ) +
      theme_bw() +
      theme(
        legend.text=element_text(size=12),
        legend.title = element_text(size = 12),
        legend.position = c(.9, .85),
        legend.background = element_rect(linetype = "solid", colour = "black"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12)
      ) +
      scale_x_continuous(breaks = seq(400, 700, by = 50))}
  
  # TIPO 1: Espectros de emisión modificados
  else if (type==1) {
    # Se itera, de forma que se añada un identificador para cada dataframe
    for (i in seq_along(l)) {
      l[[i]]$id <- paste(LETTERS[i])  # Se agrega un identificador
    }
    
    # Se combina todos los dataframes en uno solo
    data_combinada <- do.call(rbind, l)
    
    # Se escala la intensidad en unidades arbitrarias
    int_norm <- data_combinada$V3 / max(data_combinada$V3)
    data_combinada$V3 <- int_norm
    
    # Se crea la gráfica
    graph <- ggplot(data_combinada, aes(x = V1, y = V3, color = id)) +
      geom_line(size = 0.75) +
      scale_color_manual(values = color_val, name = "Posición" ) +
      labs(
        x = "Longitud de onda / (nm)",
        y = "Intensidad / (u.a)"
      ) +
      theme_bw() +
      theme(
        legend.text=element_text(size=12),
        legend.title = element_text(size = 12),
        legend.position = c(.9, .85),
        legend.background = element_rect(linetype = "solid", colour = "black"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12)
      ) +
      scale_x_continuous(breaks = seq(400, 700, by = 50))}
  
  # TIPO 2: Espectros de excitación regulares
  else if (type==2) {
    # Se itera, de forma que se añada un identificador para cada dataframe
    lambda <- c("510 nm", "610 nm")
    for (i in seq_along(l)) {
      l[[i]]$id <- paste(lambda[i])  # Se agrega un identificador
    }
    
    # Se combina todos los dataframes en uno solo
    data_combinada <- do.call(rbind, l)
    
    # Se escala la intensidad en unidades arbitrarias
    int_norm <- data_combinada$V2 / max(data_combinada$V2)
    data_combinada$V2 <- int_norm
    
    # La siguiente línea grafica el espectro
    graph <- ggplot(data_combinada, aes(x = V1, y = V2, color = id)) +
      geom_line(size = 0.75) +
      scale_color_manual(values = color_val, name = "Longitud de onda de emisión" ) +
      labs(
        x = "Longitud de onda / (nm)",
        y = "Intensidad / (u.a)"
      ) +
      theme_bw() +
      theme(
        legend.text=element_text(size=12),
        legend.title = element_text(size = 12),
        legend.position = c(.725, .85),
        legend.background = element_rect(linetype = "solid", colour = "black"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12)
      ) +
      scale_x_continuous(breaks = seq(min(data_combinada$V1), max(data_combinada$V1), by = 10))}
  
  # TIPO 3: Espectros de excitación modificados
  else if (type==3) {
    # Se itera, de forma que se añada un identificador para cada dataframe
    lambda <- c("510 nm", "610 nm")
    for (i in seq_along(l)) {
      l[[i]]$id <- paste(lambda[i])  # Se agrega un identificador
    }
    
    # Se combina todos los dataframes en uno solo
    data_combinada <- do.call(rbind, l)
    
    # Se escala la intensidad en unidades arbitrarias
    int_norm <- data_combinada$V3 / max(data_combinada$V3)
    data_combinada$V3 <- int_norm
    
    # La siguiente línea grafica el espectro
    graph <- ggplot(data_combinada, aes(x = V2, y = V3, color = id)) +
      geom_line(size = 0.75) +
      scale_color_manual(values = color_val, name = "Longitud de onda de emisión" ) +
      labs(
        x = "Longitud de onda / (nm)",
        y = "Intensidad / (u.a)"
      ) +
      theme_bw() +
      theme(
        legend.text=element_text(size=12),
        legend.title = element_text(size = 12),
        legend.position = c(.725, .85),
        legend.background = element_rect(linetype = "solid", colour = "black"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12)
      ) +
      scale_x_continuous(breaks = seq(min(data_combinada$V2), max(data_combinada$V2), by = 10))}
  
  return(graph)
}

# GENERACIÓN DE FIGURAS #

# Billete 1 Mil
# 3D
plot_3D(Bill1mil_3D_PosA_Cont)
plot_3D(Bill1mil_3D_PosB_Cont)
# EmExc365nm
plot_spectrum(Bill1milEmExc365nm_PosA,0)
plot_spectrum(Bill1milEmExc365nm_PosB,0)
plot_espectra(list(Bill1milEmExc365nm_PosA, Bill1milEmExc365nm_PosB),0)
# ExcEm 536nm, 610nm
plot_spectrum(Bill1milExcEm536nm,3)
plot_spectrum(Bill1milExcEm610nm,3)
plot_espectra(list(Bill1milExcEm536nm, Bill1milExcEm610nm), 2)

# Billete 2 Mil
# 3D
plot_3D(Bill2mil_3D_PosA_Cont)
plot_3D(Bill2mil_3D_PosB_Cont)
# EmExc365nm
plot_spectrum(Bill2milEmExc365nm_PosA,0)
plot_spectrum(Bill2milEmExc365nm_PosB,0)
plot_espectra(list(Bill2milEmExc365nm_PosA, Bill2milEmExc365nm_PosB),0)
# ExcEm 510nm, 610nm
plot_spectrum(Bill2milExcEm510nm,3)
plot_spectrum(Bill2milExcEm610nm,3)
plot_espectra(list(Bill2milExcEm510nm, Bill2milExcEm610nm), 2)

# Billete 5 Mil
# 3D
plot_3D(Bill5mil_3D_PosA_Cont)
plot_3D(Bill5mil_3D_PosB_Cont)
# EmExc365nm
plot_spectrum(Bill5milEmExc365nm_PosA,0)
plot_spectrum(Bill5milEmExc365nm_PosB,0)
plot_espectra(list(Bill5milEmExc365nm_PosA, Bill5milEmExc365nm_PosB),0)
# ExcEm 510nm, 610nm
plot_spectrum(Bill5milExcEm510nm,3)
plot_spectrum(Bill5milExcEm610nm,3)
plot_espectra(list(Bill5milExcEm510nm, Bill5milExcEm610nm), 2)

# Billete 10 Mil
# 3D
plot_3D(Bill10mil_3D_PosA_Cont)
plot_3D(Bill10mil_3D_PosB_Cont)
# EmExc365nm
plot_spectrum(Bill10milEmExc365nm_PosA,1)
plot_spectrum(Bill10milEmExc365nm_PosB,1)
plot_espectra(list(Bill10milEmExc365nm_PosA, Bill10milEmExc365nm_PosB),1)
# ExcEm 510nm, 610nm
plot_spectrum(Bill10milExcEm510nm,2)
plot_spectrum(Bill10milExcEm610nm,2)
plot_espectra(list(Bill10milExcEm510nm, Bill10milExcEm610nm), 3)

# Billete 20 Mil
# 3D
plot_3D(Bill20mil_3D_PosA_Cont)
plot_3D(Bill20mil_3D_PosB_Cont)
# EmExc365nm
plot_spectrum(Bill20milEmExc365nm_PosA,1)
plot_spectrum(Bill20milEmExc365nm_PosB,1)
plot_espectra(list(Bill20milEmExc365nm_PosA, Bill20milEmExc365nm_PosB),1)
# ExcEm 510nm, 610nm
plot_spectrum(Bill20milExcEm510nm,2)
plot_spectrum(Bill20milExcEm610nm,2)
plot_espectra(list(Bill20milExcEm510nm, Bill20milExcEm610nm), 3)


### ESPECTROS IR ###


# IMPORTACIÓN DE DATOS #

# Las siguientes líneas leen los archivos .csv y los convierte en data frames (df)

# Billete de 1 Mil
Bill1Mil_VentInfDer = read.csv("C:\\Users\\james\\Google Drive\\TEC\\2024\\II SEMESTRE 2024\\TFG\\Espectros\\Espectros Art1\\ESPECTROS IR\\Billete 1 Mil\\Ventana transparente\\Bill1Mil VentInfDer.CSV", header = FALSE)
Bill1Mil_VentSupDer = read.csv("C:\\Users\\james\\Google Drive\\TEC\\2024\\II SEMESTRE 2024\\TFG\\Espectros\\Espectros Art1\\ESPECTROS IR\\Billete 1 Mil\\Ventana transparente\\Bill1Mil VentSupDer.CSV", header = FALSE)
Bill1Mil_Mapa = read.csv("C:\\Users\\james\\Google Drive\\TEC\\2024\\II SEMESTRE 2024\\TFG\\Espectros\\Espectros Art1\\ESPECTROS IR\\Billete 1 Mil\\Bill1Mil Mapa Iridiscente.CSV", header = FALSE)
Bill1Mil_MarcaF = read.csv("C:\\Users\\james\\Google Drive\\TEC\\2024\\II SEMESTRE 2024\\TFG\\Espectros\\Espectros Art1\\ESPECTROS IR\\Billete 1 Mil\\Bill1Mil MarcaFluor.CSV", header = FALSE)

# Billete de 2 Mil
Bill2Mil_VentInfDer = read.csv("C:\\Users\\james\\Google Drive\\TEC\\2024\\II SEMESTRE 2024\\TFG\\Espectros\\Espectros Art1\\ESPECTROS IR\\Billete 2 Mil\\Ventana transparente\\Bill2Mil VentInfDer.CSV", header = FALSE)
Bill2Mil_VentSupDer = read.csv("C:\\Users\\james\\Google Drive\\TEC\\2024\\II SEMESTRE 2024\\TFG\\Espectros\\Espectros Art1\\ESPECTROS IR\\Billete 2 Mil\\Ventana transparente\\Bill2Mil VentSupDer.CSV", header = FALSE)
Bill2Mil_Mapa = read.csv("C:\\Users\\james\\Google Drive\\TEC\\2024\\II SEMESTRE 2024\\TFG\\Espectros\\Espectros Art1\\ESPECTROS IR\\Billete 2 Mil\\Bill2Mil Mapa Iridiscente.CSV", header = FALSE)
Bill2Mil_MarcaF = read.csv("C:\\Users\\james\\Google Drive\\TEC\\2024\\II SEMESTRE 2024\\TFG\\Espectros\\Espectros Art1\\ESPECTROS IR\\Billete 2 Mil\\Bill2Mil MarcaFluor.CSV", header = FALSE)

# Billete de 5 Mil
Bill5Mil_VentInfDer = read.csv("C:\\Users\\james\\Google Drive\\TEC\\2024\\II SEMESTRE 2024\\TFG\\Espectros\\Espectros Art1\\ESPECTROS IR\\Billete 5 Mil\\Ventana transparente\\Bill5Mil VentInfDer.CSV", header = FALSE)
Bill5Mil_VentSupDer = read.csv("C:\\Users\\james\\Google Drive\\TEC\\2024\\II SEMESTRE 2024\\TFG\\Espectros\\Espectros Art1\\ESPECTROS IR\\Billete 5 Mil\\Ventana transparente\\Bill5Mil VentSupDer.CSV", header = FALSE)
Bill5Mil_Mapa = read.csv("C:\\Users\\james\\Google Drive\\TEC\\2024\\II SEMESTRE 2024\\TFG\\Espectros\\Espectros Art1\\ESPECTROS IR\\Billete 5 Mil\\Bill5Mil Mapa Iridiscente.CSV", header = FALSE)
Bill5Mil_MarcaF = read.csv("C:\\Users\\james\\Google Drive\\TEC\\2024\\II SEMESTRE 2024\\TFG\\Espectros\\Espectros Art1\\ESPECTROS IR\\Billete 5 Mil\\Bill5Mil MarcaFluor.CSV", header = FALSE)

# Billete de 10 Mil 
Bill10Mil_VentInfDer = read.csv("C:\\Users\\james\\Google Drive\\TEC\\2024\\II SEMESTRE 2024\\TFG\\Espectros\\Espectros Art1\\ESPECTROS IR\\Billete 10 Mil\\Ventana transparente\\Bill10Mil VentInfDer.CSV", header = FALSE)
Bill10Mil_VentSupDer = read.csv("C:\\Users\\james\\Google Drive\\TEC\\2024\\II SEMESTRE 2024\\TFG\\Espectros\\Espectros Art1\\ESPECTROS IR\\Billete 10 Mil\\Ventana transparente\\Bill10Mil VentSupDer.CSV", header = FALSE)
Bill10Mil_Mapa = read.csv("C:\\Users\\james\\Google Drive\\TEC\\2024\\II SEMESTRE 2024\\TFG\\Espectros\\Espectros Art1\\ESPECTROS IR\\Billete 10 Mil\\Bill10Mil Mapa Iridiscente.CSV", header = FALSE)
Bill10Mil_MarcaF = read.csv("C:\\Users\\james\\Google Drive\\TEC\\2024\\II SEMESTRE 2024\\TFG\\Espectros\\Espectros Art1\\ESPECTROS IR\\Billete 10 Mil\\Bill10Mil MarcaFluor.CSV", header = FALSE)

# Billete de 20 Mil 
Bill20Mil_VentInfDer = read.csv("C:\\Users\\james\\Google Drive\\TEC\\2024\\II SEMESTRE 2024\\TFG\\Espectros\\Espectros Art1\\ESPECTROS IR\\Billete 20 Mil\\Ventana transparente\\Bill20Mil VentInfDer.CSV", header = FALSE)
Bill20Mil_VentSupDer = read.csv("C:\\Users\\james\\Google Drive\\TEC\\2024\\II SEMESTRE 2024\\TFG\\Espectros\\Espectros Art1\\ESPECTROS IR\\Billete 20 Mil\\Ventana transparente\\Bill20Mil VentSupDer.CSV", header = FALSE)
Bill20Mil_Mapa = read.csv("C:\\Users\\james\\Google Drive\\TEC\\2024\\II SEMESTRE 2024\\TFG\\Espectros\\Espectros Art1\\ESPECTROS IR\\Billete 20 Mil\\Bill20Mil Mapa Iridiscente.CSV", header = FALSE)
Bill20Mil_MarcaF = read.csv("C:\\Users\\james\\Google Drive\\TEC\\2024\\II SEMESTRE 2024\\TFG\\Espectros\\Espectros Art1\\ESPECTROS IR\\Billete 20 Mil\\Bill20Mil MarcaFluor.CSV", header = FALSE)

# FUNCIONES #

plot_ir_spectrum <- function(df) {
  # Descripción: Función para graficar espectros IR en 2D , así como determinar
  # los picos locales de absorbancia
  # Entradas:
  # df = dataframes del espectro a graficar
  # type = tipo de gráfico
  #        0 = Espectro de emisión regular
  #        1 = Espectro de emisión a partir de datos medidos en espectro 3D
  # Salidas:
  # Gráfico del espectro IR
  # Picos locales de absorbancia máxima
  
    graph = ggplot(mapping = aes(x=df$V1, y=df$V2)) +
      geom_line(size = 0.75) +
      labs(
        x = "Número de onda / (cm^-1)",
        y = "Absorbancia"
      ) +
      theme_bw() +
      theme(
        legend.text=element_text(size=12),
        legend.title = element_text(size = 12),
        legend.position = c(.85, .85),
        legend.background = element_rect(linetype = "solid", colour = "black"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12)
      ) +
      scale_x_reverse(breaks = seq(500, 4000, by = 500))
    
    # La siguiente línea determina los picos máximos locales
    #localPeaks = peaks(x = df, strict = TRUE, span = 11, x.var.name = "V1", y.var.name = "V2")
    #localPeaks = spikes(x = df, z.threshold = 0.2, x.var.name = "V1", y.var.name = "V2")
    localPeaks = get_peaks(x = df$V1, y = df$V2, strict = TRUE, span = 11)
    locPeaks = localPeaks[order(-localPeaks$y),]
    print("Picos máximos locales")
    print(locPeaks$x)
    absPeak = df[which.max(df$V2),]
    print("Pico máximo global")
    print(absPeak)
    print(graph)}

# Función para varios espectros IR en conjunto en un mismo plot
plot_ir_espectra <- function(l,mod){
  # Donde l: lista de dataframes de los espectros IR a graficar
  #     mod: modo del gráfico:
  #          1 es para comparar espectros en una misma ventana)
  #          2 es para comparar espectros entre billetes (Mapa/MarcaFluor)
  
  # Vector con los colores de las gráficas
  color_val <- c("red","blue", "yellow", "green", "orange")
  
  if (mod == 1) { # MOD 1: COMPARAR EN UNA MISMA VENTANA
  # Se itera, de forma que se añada un identificador para cada dataframe
    pos <- c("Inferior derecha", "Superior derecha")
  for (i in seq_along(l)) {
    l[[i]]$id <- paste(pos[i])  # Se agrega un identificador
  }
  
  # Se combina todos los dataframes en uno solo
  data_combinada <- do.call(rbind, l)

    # Se crea la gráfica
    graph <- ggplot(data_combinada, aes(x = V1, y = V2, color = id)) +
      geom_line(size = 0.75) +
      scale_color_manual(values = color_val, name = "Ventana" ) +
      labs(
        x = "Número de onda / (cm^-1)",
        y = "Absorbancia"
      ) +
      theme_bw() +
      theme(
        legend.text=element_text(size=10),
        legend.title = element_text(size = 10),
        legend.position = c(0.19, 0.84),
        legend.background = element_rect(linetype = "solid", colour = "black"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12)
      ) +
      scale_x_reverse(breaks = seq(500, 4000, by = 500))}
  
  else if (mod == 2) { # MOD 2: COMPARAR ENTRE BILLETES
    # Se itera, de forma que se añada un identificador para cada dataframe
    denominacion = c(1,2,5,10,20)
    for (i in seq_along(l)) {
      l[[i]]$id <- paste(denominacion[i], "Mil")  # Se agrega un identificador
    }
    
    # Se combina todos los dataframes en uno solo
    data_combinada <- do.call(rbind, l)
    
    data_combinada$id <- factor(data_combinada$id, levels = c("1 Mil", "2 Mil", "5 Mil", "10 Mil", "20 Mil"))
    
    # Se crea la gráfica
    graph <- ggplot(data_combinada, aes(x = V1, y = V2, color = id)) +
      geom_line(size = 0.75) +
      scale_color_manual(values = c("1 Mil" = "red",
                                    "2 Mil" = "blue1",
                                    "5 Mil" = "yellow2",
                                    "10 Mil" = "green3",
                                    "20 Mil" = "orange4"), 
                         name = "Billete", 
                         breaks = c("1 Mil", 
                                    "2 Mil", 
                                    "5 Mil", 
                                    "10 Mil", 
                                    "20 Mil")) +
      labs(
        x = "Número de onda / (cm^-1)",
        y = "Absorbancia"
      ) +
      theme_bw() +
      theme(
        legend.text=element_text(size=10),
        legend.title = element_text(size = 10),
        legend.position = c(0.093, 0.761),
        legend.background = element_rect(linetype = "solid", colour = "black"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12)
      ) +
      scale_x_reverse(breaks = seq(500, 4000, by = 500))
  }
  return(graph)
}

    
# GENERACIÓN DE FIGURAS #

# Delimitación de 500 cm^-1 a 4000 cm^-1
# Marca Fluorescente
Bill1Mil_MarcaF <- Bill1Mil_MarcaF[Bill1Mil_MarcaF$V1 >= 500, ]
Bill2Mil_MarcaF <- Bill2Mil_MarcaF[Bill2Mil_MarcaF$V1 >= 500, ]
Bill5Mil_MarcaF <- Bill5Mil_MarcaF[Bill5Mil_MarcaF$V1 >= 500, ]
Bill10Mil_MarcaF <- Bill10Mil_MarcaF[Bill10Mil_MarcaF$V1 >= 500, ]
Bill20Mil_MarcaF <- Bill20Mil_MarcaF[Bill20Mil_MarcaF$V1 >= 500, ]

# Billete de 1 Mil
plot_ir_spectrum(Bill1Mil_VentInfDer)
plot_ir_spectrum(Bill1Mil_VentSupDer)
plot_ir_espectra(list(Bill1Mil_VentInfDer, Bill1Mil_VentSupDer), 1)
plot_ir_spectrum(Bill1Mil_Mapa)
plot_ir_spectrum(Bill1Mil_MarcaF)

# Billete de 2 Mil
plot_ir_spectrum(Bill2Mil_VentInfDer)
plot_ir_spectrum(Bill2Mil_VentSupDer)
plot_ir_espectra(list(Bill2Mil_VentInfDer, Bill2Mil_VentSupDer), 1)
plot_ir_spectrum(Bill2Mil_Mapa)
plot_ir_spectrum(Bill2Mil_MarcaF)

# Billete de 5 Mil
plot_ir_spectrum(Bill5Mil_VentInfDer)
plot_ir_spectrum(Bill5Mil_VentSupDer)
plot_ir_espectra(list(Bill5Mil_VentInfDer, Bill5Mil_VentSupDer), 1)
plot_ir_spectrum(Bill5Mil_Mapa)
plot_ir_spectrum(Bill5Mil_MarcaF)

# Billete de 10 Mil 
plot_ir_spectrum(Bill10Mil_VentInfDer)
plot_ir_spectrum(Bill10Mil_VentSupDer)
plot_ir_espectra(list(Bill10Mil_VentInfDer, Bill10Mil_VentSupDer), 1)
plot_ir_spectrum(Bill10Mil_Mapa)
plot_ir_spectrum(Bill10Mil_MarcaF)

# Billete de 20 Mil 
plot_ir_spectrum(Bill20Mil_VentInfDer)
plot_ir_spectrum(Bill20Mil_VentSupDer)
plot_ir_espectra(list(Bill20Mil_VentInfDer, Bill20Mil_VentSupDer), 1)
plot_ir_spectrum(Bill20Mil_Mapa)
plot_ir_spectrum(Bill20Mil_MarcaF)

# Comparaciones entre billetes

plot_ir_espectra(list(Bill1Mil_MarcaF,
                      Bill2Mil_MarcaF, 
                      Bill5Mil_MarcaF, 
                      Bill10Mil_MarcaF, 
                      Bill20Mil_MarcaF), 2)

# Mapa Iridiscente
plot_ir_espectra(list(Bill1Mil_Mapa, 
                      Bill2Mil_Mapa, 
                      Bill5Mil_Mapa, 
                      Bill10Mil_Mapa, 
                      Bill20Mil_Mapa), 2)
