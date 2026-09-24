# =============================================================================
# IIO422 Hidrología - Ayudantía 04: Reproducibilidad de un análisis hidrológico
# Autora: Fernanda Ravanal Medina
#
# Instrucciones: Ay04.pdf [lámina 7/7], paso 2
# Pruebas:       Ay03a-homogeneidad_y_estacionariedad.pdf (J. Morán, 17-04-2026)
#
# Uso, desde la raíz del repositorio en Git Bash:
#   Rscript Rscript/Ravanal_Medina.R <Código_BNA>
#   ej.: Rscript Rscript/Ravanal_Medina.R 9123001
# En Positron o RStudio: abrir la carpeta del repo, fijar `codigo.manual`
# (sección 1) y usar Source.
#
# Marcas de los comentarios:
#   "2.1)", "2.2)"   -> sub-ítems del paso 2 de la lámina 7 de la Ay04
#   "[Ay03a lám. N]" -> lámina de la Ay03a de donde sale la prueba o el código
#   "+)"             -> agregado mío, no lo pide la lámina
#
# El archivo está en UTF-8. Si en RStudio se ven mal las tildes:
#   File > Reopen with Encoding > UTF-8
# =============================================================================


# 0) Directorio de trabajo [Ay04 lám. 7, paso 0] ------------------------------
# En la Ay03a el ayudante fija main.drty con la ruta de su PC [Ay03a lám. 13].
# +) Aquí se toma la carpeta desde donde se corre, para que el script funcione
#    igual en cualquier computador que clone el repositorio.
main.drty <- getwd()
if (!dir.exists(file.path(main.drty, "data"))) {
  stop("Corre el script desde la raíz del repositorio (donde están data/ y Rscript/).")
}
setwd(main.drty)

library(trend)      # Pettitt, Buishand, SNHT   [Ay03a lám. 11]
library(DescTools)  # Von Neumann               [Ay03a lám. 11]
library(tseries)    # ADF y KPSS                [Ay03a lám. 21]


# 1) Estación asignada (Código BNA, el mismo nombre de la rama) ---------------
codigo.manual <- "8317001"   # mi estación; en Git Bash se puede pasar otra como argumento

args <- commandArgs(trailingOnly = TRUE)
codigo <- if (length(args) >= 1) args[1] else codigo.manual

disponibles <- sub("^precip_cr2met_day_(.*)\\.csv$", "\\1",
                   list.files("data", pattern = "^precip_cr2met_day_.*\\.csv$"))
if (!nzchar(codigo)) {
  stop("Falta el Código BNA. Disponibles: ", paste(disponibles, collapse = ", "))
}
if (!codigo %in% disponibles) {
  stop("No existe data/precip_cr2met_day_", codigo, ".csv. Disponibles: ",
       paste(disponibles, collapse = ", "))
}

dir.create("results", showWarnings = FALSE)
salida <- function(nombre) file.path("results", paste0(codigo, "_", nombre))


# 2.1) Análisis exploratorio de la base de datos asignada ---------------------
data.pp <- read.csv(file.path("data", paste0("precip_cr2met_day_", codigo, ".csv")),
                    check.names = FALSE)
# La 5a columna se llama igual que el código BNA; se renombra para trabajarla.
stopifnot(names(data.pp)[5] == codigo)
names(data.pp)[5] <- "pp"
data.pp$date <- as.Date(data.pp$date)

# +) Controles antes de analizar: fechas continuas, sin duplicados, sin
#    valores negativos. Si algo falla, el script se detiene aquí.
stopifnot(!anyNA(data.pp$date),
          !any(duplicated(data.pp$date)),
          all(diff(data.pp$date) == 1),
          all(data.pp$pp >= 0, na.rm = TRUE))
n.na <- sum(is.na(data.pp$pp))

# +) Umbral de día lluvioso. 1 mm es un criterio usual; se deja como parámetro.
umbral.lluvia <- 1   # mm

# Totales anuales: solo años calendario completos y sin faltantes.
# (La serie termina el 30 de abril de 2020, así que 2020 queda fuera.)
dias.anio <- tapply(!is.na(data.pp$pp), data.pp$year, sum)
dias.esperados <- sapply(as.integer(names(dias.anio)), function(y)
  as.numeric(as.Date(paste0(y + 1, "-01-01")) - as.Date(paste0(y, "-01-01"))))
anios.completos <- as.integer(names(dias.anio))[dias.anio == dias.esperados]

sub.anual <- data.pp[data.pp$year %in% anios.completos, ]
anual <- data.frame(year = anios.completos,
                    pp   = as.numeric(tapply(sub.anual$pp, sub.anual$year, sum)))

# Totales mensuales (todos los meses de la serie están completos).
mensual <- aggregate(pp ~ year + month, data = data.pp, FUN = sum)
mensual <- mensual[order(mensual$year, mensual$month), ]
stopifnot(nrow(mensual) == length(unique(format(data.pp$date, "%Y-%m"))))

# +) Anomalías mensuales: total del mes menos el promedio de ese mes en toda la
#    serie. Quita el ciclo estacional, que de otro modo "engaña" a las pruebas.
clima.mes <- tapply(mensual$pp, mensual$month, mean)
mensual$anom <- mensual$pp - clima.mes[as.character(mensual$month)]

# Resumen exploratorio
dia.max  <- data.pp[which.max(data.pp$pp), ]
anio.max <- anual[which.max(anual$pp), ]
anio.min <- anual[which.min(anual$pp), ]
mes.max  <- mensual[which.max(mensual$pp), ]
mes.min  <- mensual[which.min(mensual$pp), ]
meses <- c("ene", "feb", "mar", "abr", "may", "jun",
           "jul", "ago", "sep", "oct", "nov", "dic")

# +) Datos faltantes por año: si hubiera NA, aquí se ve dónde se concentran.
na.anio <- data.frame(year = as.integer(names(dias.anio)),
                      dias = as.integer(table(data.pp$year)),
                      na   = as.integer(tapply(is.na(data.pp$pp), data.pp$year, sum)))
write.csv(na.anio, salida("na_por_anio.csv"), row.names = FALSE)

# +) Media diaria en los años hidrológicos 1990-2009 (1-abr-1990 a 31-mar-2010),
#    el mismo período con que CAMELS-CL calcula su atributo p_mean_cr2met.
#    Sirve para contrastar la serie contra un valor publicado.
periodo.camels <- data.pp$date >= as.Date("1990-04-01") &
                  data.pp$date <= as.Date("2010-03-31")
media.camels <- mean(data.pp$pp[periodo.camels], na.rm = TRUE)

exploratorio <- data.frame(
  variable = c("Código BNA", "Fecha inicial", "Fecha final", "Días de registro",
               "Días sin dato", "Días con pp >= 1 mm (%)",
               "Días con pp = 0 (%)", "Precipitación diaria media (mm)",
               "Mínimo diario (mm)",
               "Máximo diario (mm)", "Fecha del máximo diario",
               "Mínimo mensual (mm)", "Mes del mínimo mensual",
               "Máximo mensual (mm)", "Mes del máximo mensual",
               "Media diaria abr-1990 a mar-2010 (mm)",
               "Años completos", "Precipitación anual media (mm)",
               "Desviación estándar anual (mm)", "Coeficiente de variación anual (-)",
               "Año más lluvioso", "Total del año más lluvioso (mm)",
               "Año más seco", "Total del año más seco (mm)",
               "Mes más lluvioso (media)", "Mes más seco (media)"),
  valor = c(codigo,
            format(min(data.pp$date)), format(max(data.pp$date)),
            nrow(data.pp), n.na,
            round(100 * mean(data.pp$pp >= umbral.lluvia, na.rm = TRUE), 1),
            round(100 * mean(data.pp$pp == 0, na.rm = TRUE), 1),
            round(mean(data.pp$pp, na.rm = TRUE), 2),
            round(min(data.pp$pp, na.rm = TRUE), 2),
            round(dia.max$pp, 1), format(dia.max$date),
            round(mes.min$pp, 1), sprintf("%d-%02d", mes.min$year, mes.min$month),
            round(mes.max$pp, 1), sprintf("%d-%02d", mes.max$year, mes.max$month),
            round(media.camels, 4),
            paste0(nrow(anual), " (", min(anual$year), "-", max(anual$year), ")"),
            round(mean(anual$pp), 1), round(sd(anual$pp), 1),
            round(sd(anual$pp) / mean(anual$pp), 3),
            anio.max$year, round(anio.max$pp, 1),
            anio.min$year, round(anio.min$pp, 1),
            meses[which.max(clima.mes)], meses[which.min(clima.mes)])
)
print(exploratorio, row.names = FALSE)

write.csv(exploratorio, salida("exploratorio.csv"), row.names = FALSE)
write.csv(anual, salida("anual.csv"), row.names = FALSE)
write.csv(mensual, salida("mensual.csv"), row.names = FALSE)

# Figuras exploratorias
png(salida("fig1_serie_diaria.png"), width = 1800, height = 700, res = 200)
plot(data.pp$date, data.pp$pp, type = "h", col = "steelblue",
     xlab = "Fecha", ylab = "Precipitación diaria (mm)",
     main = paste("Estación", codigo, "- precipitación diaria CR2MET"))
invisible(dev.off())

png(salida("fig2_totales_anuales.png"), width = 1800, height = 800, res = 200)
barplot(anual$pp, names.arg = anual$year, las = 2, cex.names = 0.7,
        col = "steelblue", border = NA,
        ylab = "Precipitación anual (mm)",
        main = paste("Estación", codigo, "- totales anuales"))
abline(h = mean(anual$pp), lty = 2, col = "firebrick")
legend("topright", legend = "Media", lty = 2, col = "firebrick", bty = "n")
invisible(dev.off())

png(salida("fig3_ciclo_mensual.png"), width = 1600, height = 900, res = 200)
boxplot(pp ~ month, data = mensual, names = meses, col = "lightsteelblue",
        xlab = "Mes", ylab = "Precipitación mensual (mm)",
        main = paste("Estación", codigo, "- distribución de totales mensuales"))
invisible(dev.off())

# +) Serie mensual: totales de cada mes en el tiempo, con la media móvil de
#    12 meses para ver si el nivel cambia.
fecha.mes <- as.Date(sprintf("%d-%02d-15", mensual$year, mensual$month))
movil12   <- stats::filter(mensual$pp, rep(1 / 12, 12), sides = 2)
png(salida("fig5_serie_mensual.png"), width = 1800, height = 800, res = 200)
plot(fecha.mes, mensual$pp, type = "l", col = "steelblue",
     xlab = "Fecha", ylab = "Precipitación mensual (mm)",
     main = paste("Estación", codigo, "- serie de totales mensuales"))
lines(fecha.mes, movil12, col = "firebrick", lwd = 2)
legend("topright", legend = c("Total mensual", "Media móvil 12 meses"),
       col = c("steelblue", "firebrick"), lwd = c(1, 2), bty = "n")
invisible(dev.off())

# +) Histogramas. El diario usa solo días con lluvia (pp > 0): con todos los
#    días, la barra del cero tapa el resto de la distribución.
png(salida("fig6_histogramas.png"), width = 1800, height = 800, res = 200)
old.par <- par(mfrow = c(1, 2))
hist(data.pp$pp[data.pp$pp > 0], breaks = 50, col = "lightsteelblue",
     border = "white", xlab = "Precipitación diaria (mm)", ylab = "Frecuencia",
     main = "Días con pp > 0")
hist(mensual$pp, breaks = 30, col = "lightsteelblue", border = "white",
     xlab = "Precipitación mensual (mm)", ylab = "Frecuencia",
     main = "Totales mensuales")
par(old.par)
invisible(dev.off())


# 2.2) Análisis de homogeneidad y estacionariedad -----------------------------
# Se aplican a los totales anuales (serie principal) y, como complemento, a las
# anomalías mensuales. +) No se usan los datos diarios: la mayoría de los días
# son cero y la serie diaria está muy autocorrelacionada, lo que distorsiona
# las pruebas.
#
# Buishand y SNHT calculan el p-valor por simulación, por eso se fija la semilla.
# El ayudante usa m = 50 para que corra rápido [Ay03a lám. 15]. +) Aquí se usa
# m = 20000 (el valor por defecto de trend) para que el p-valor sea estable.
set.seed(20260924)
m.sim <- 20000
alfa  <- 0.05

aplicar.tests <- function(x, etiquetas, nombre.serie) {
  # Pettitt: cambio en la mediana                          [Ay03a lám. 12-13]
  pt <- pettitt.test(x)
  # Rango de Buishand: cambio en la media                  [Ay03a lám. 14-15]
  br <- br.test(x, m = m.sim)
  # SNHT: cambio en la media                               [Ay03a lám. 16-17]
  sn <- snh.test(x, m = m.sim)
  # Von Neumann: aleatoriedad global, sin punto de quiebre [Ay03a lám. 18-19]
  vn <- VonNeumannTest(x)
  # ADF (H0: raíz unitaria) y KPSS (H0: estacionaria)      [Ay03a lám. 21-22]
  # +) tseries interpola el p-valor en una tabla y lo recorta a sus bordes
  #    (ADF: 0.01-0.99; KPSS: 0.01-0.10). Se guarda el aviso para declararlo.
  capturar <- function(expr) {
    aviso <- ""
    res <- withCallingHandlers(expr, warning = function(w) {
      aviso <<- conditionMessage(w); invokeRestart("muffleWarning")
    })
    list(res = res, aviso = aviso)
  }
  adf <- capturar(adf.test(x))
  kp  <- capturar(kpss.test(x, null = "Level"))

  cambio <- function(K) etiquetas[as.integer(K)]
  decide.homog <- function(p) ifelse(p < alfa, "Rechaza H0: no homogénea",
                                     "No rechaza H0: homogénea")
  data.frame(
    serie       = nombre.serie,
    prueba      = c("Pettitt", "Buishand", "SNHT", "Von Neumann", "ADF", "KPSS"),
    H0          = c("Sin cambio en la mediana", "Sin cambio en la media",
                    "Sin cambio en la media", "Serie aleatoria (homogénea)",
                    "Raíz unitaria (no estacionaria)", "Estacionaria en nivel"),
    estadistico = c(unname(pt$statistic), unname(br$statistic), unname(sn$statistic),
                    unname(vn$statistic[1]), unname(adf$res$statistic),
                    unname(kp$res$statistic)),
    p_valor     = c(pt$p.value, br$p.value, sn$p.value, vn$p.value,
                    adf$res$p.value, kp$res$p.value),
    p_recortado = c("", "", "", "", adf$aviso, kp$aviso),
    ultimo_antes_del_cambio = c(cambio(pt$estimate), cambio(br$estimate),
                                cambio(sn$estimate), NA, NA, NA),
    decision    = c(decide.homog(c(pt$p.value, br$p.value, sn$p.value, vn$p.value)),
                    ifelse(adf$res$p.value < alfa,
                           "Rechaza H0: estacionaria", "No rechaza H0: no estacionaria"),
                    ifelse(kp$res$p.value < alfa,
                           "Rechaza H0: no estacionaria", "No rechaza H0: estacionaria")),
    stringsAsFactors = FALSE
  )
}

tests.anual <- aplicar.tests(anual$pp, as.character(anual$year), "Total anual")
tests.mensual <- aplicar.tests(mensual$anom,
                               sprintf("%d-%02d", mensual$year, mensual$month),
                               "Anomalía mensual")
tests <- rbind(tests.anual, tests.mensual)
tests$p_recortado <- ifelse(nzchar(tests$p_recortado), "sí", "no")
print(tests[, c("serie", "prueba", "estadistico", "p_valor",
                "ultimo_antes_del_cambio", "decision")], row.names = FALSE)
write.csv(tests, salida("tests.csv"), row.names = FALSE)

# +) Figura del quiebre de Pettitt en los totales anuales, con la media de cada
#    tramo. La prueba dice si hay quiebre; esta figura muestra de cuánto es.
K <- as.integer(pettitt.test(anual$pp)$estimate)
antes   <- anual$pp[seq_len(K)]
despues <- anual$pp[(K + 1):nrow(anual)]
png(salida("fig7_quiebre_pettitt.png"), width = 1800, height = 800, res = 200)
plot(anual$year, anual$pp, type = "b", pch = 16, col = "steelblue",
     xlab = "Año", ylab = "Precipitación anual (mm)",
     main = paste("Estación", codigo, "- quiebre más probable según Pettitt"))
abline(v = anual$year[K] + 0.5, lty = 2, col = "firebrick")
segments(anual$year[1], mean(antes), anual$year[K], mean(antes),
         lwd = 2, col = "darkorange")
segments(anual$year[K + 1], mean(despues), anual$year[nrow(anual)], mean(despues),
         lwd = 2, col = "darkorange")
legend("topright", legend = c("Quiebre más probable", "Media de cada tramo"),
       lty = c(2, 1), lwd = c(1, 2), col = c("firebrick", "darkorange"), bty = "n")
invisible(dev.off())

write.csv(data.frame(ultimo_anio_antes = anual$year[K],
                     media_antes = mean(antes), media_despues = mean(despues),
                     diferencia = mean(despues) - mean(antes),
                     diferencia_pct = 100 * (mean(despues) - mean(antes)) / mean(antes)),
          salida("quiebre_pettitt.csv"), row.names = FALSE)

# +) Registro de versiones, para que otro pueda reproducir el resultado.
writeLines(capture.output(sessionInfo()), salida("sessionInfo.txt"))
cat("\nListo. Resultados en results/", codigo, "_*\n", sep = "")
