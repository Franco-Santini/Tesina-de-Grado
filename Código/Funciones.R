# Librerias
library(dplyr)
library(lubridate)
library(fpp3)
library(ggplot2)
library(ggrepel)
library(MASS) # Transformación Box-Cox

# Lectura de la base de datos
pasajeros <- readxl::read_excel("./Datos/solo pasajeros 2006-2024.xlsx")

# Creo funciones para los gráficos
g1_evol_serie <- function(df){
  df |> 
    ggplot() +
    aes(x = mes_anio, y = pasajeros) + 
    geom_line() +
    geom_point(color = "dodgerblue2") +
    labs(x = "Fecha", y = "Pasajeros") +
    theme_bw()
}

g2_boxplot_disp <- function(df){
  df |> 
    ggplot() +
    aes(x = Anio, y = pasajeros, group = Anio) +
    geom_boxplot(color = "black", fill = "dodgerblue2") +
    theme_bw()
}

g3_estacionalidad_anual <- function(df){
  df |> 
    ggplot() +
    aes(x = Mes, y = pasajeros, group = Anio, color = factor(Anio)) +
    geom_line() +
    geom_text_repel(data = subset(df |> 
                                    mutate(anio = factor(year(mes_anio)),
                                           mes = factor(month(mes_anio))), mes == 1),
                    aes(color = anio, label = anio),
                    fontface = "bold",
                    size = 2.5,
                    # direction = "y",
                    segment.size = 0,
                    segment.alpha = 0,
                    nudge_x = -0.25,
                    box.padding = 0,
                    force = 0.4
    ) +
    scale_x_discrete(limits = 1:12, labels = c("Ene", "Feb", "Mar", "Abr", "May", "Jun", 
                                               "Jul", "Ago", "Sep", "Oct", "Nov", "Dic")) +
    labs(color = "Año") +
    theme_bw() +
    theme(legend.position = "none")
}

g4_acf <- function(df, conf_limit){
  autocorrelacion_2 <- acf(df$pasajeros, lag.max = 80, plot = F)
  datos_autocorrelacion <- data.frame(
    acf = autocorrelacion_2$acf,
    lag = autocorrelacion_2$lag
  )
  
  alpha <- conf_limit
  conf.lims <- c(-1,1)*qnorm((1 + alpha)/2)/sqrt(autocorrelacion_2$n.used)
  
  datos_autocorrelacion |> 
    ggplot() +
    aes(x = lag, y = acf) +
    geom_segment(aes(x = lag, xend = lag, y = 0, yend = acf), linewidth = 1) +
    geom_point(size = 1.5, color = "dodgerblue2") +
    geom_hline(yintercept=conf.lims, lty=2, col='blue') +
    geom_hline(yintercept = 0) +
    scale_x_continuous(breaks = seq(0, 80, 4)) +
    labs(x = "Rezago", y = expression(rho[k])) +
    theme_bw()
}

g5_pacf <- function(df, conf_limit){
  pautocorrelacion_2 <- pacf(df$pasajeros, lag.max = 80, plot = F)
  datos_pautocorrelacion <- data.frame(
    pacf = pautocorrelacion_2$acf,
    lag = pautocorrelacion_2$lag
  )
  
  alpha <- conf_limit
  conf.lims2 <- c(-1,1)*qnorm((1 + alpha)/2)/sqrt(pautocorrelacion_2$n.used)
  
  datos_pautocorrelacion |> 
    ggplot() +
    aes(x = lag, y = pacf) +
    geom_segment(aes(x = lag, xend = lag, y = 0, yend = pacf), linewidth = 1) +
    geom_point(size = 1.5, color = "dodgerblue2") +
    geom_hline(yintercept=conf.lims2, lty=2, col='blue') +
    geom_hline(yintercept = 0) +
    scale_x_continuous(breaks = seq(0, 80, 4)) +
    labs(y = expression(Phi[kk]), x = "Rezago") +
    theme_bw()
}

g6_graf_pronosticos <- function(estrategia, anio, niv_conf){
  df_filtrado <- datos_series |> filter(mes_anio <= yearmonth(paste("Dec", anio)))
  
  estrategia |> 
    filter(.id == paste("Pronóstico de", anio)) |> 
    autoplot(level = niv_conf) +
    geom_line(data = df_filtrado, aes(y = pasajeros)) +
    theme_bw()
}

# Creo funciones para el ajuste de los modelos
arima_ets_sol1 <- function(df){
  df |> 
    stretch_tsibble(.step = 12, .init = 168) |> 
    model(
      ets = ETS(log(pasajeros), ic = c("aicc", "aic", "bic")),
      arima = ARIMA(log(pasajeros))
    ) |> 
    forecast(h = 12) |> 
    mutate(.id = paste("Pronóstico de", .id + 2019))
}

m_int_sol2 <- function(df){
  # Construir el conjunto de datos con la variable de intervención
  datos_series_stretch <- df |> 
    mutate(pandemia = if_else((mes_anio >= yearmonth("2020 mar") & mes_anio <= yearmonth("2021 jul")), 1, 0),
           post_pandemia = if_else(mes_anio > yearmonth("2021 jul"), 1, 0)) |> 
    stretch_tsibble(.step = 12, .init = 168)
  
  # Ajustamos los modelos usando las covariables
  ajustes <- datos_series_stretch |> 
    model(
      arima1 = ARIMA(log(pasajeros)),
      arima2 = ARIMA(log(pasajeros) ~ pandemia),
      arima3 = ARIMA(log(pasajeros) ~ pandemia + post_pandemia)
    ) |> suppressWarnings()
  
  ajustes_largos <- ajustes |> 
    pivot_longer(
      cols = starts_with("arima"),
      names_to = ".model",
      values_to = "arima"
    )
  
  ajustes_final <- ajustes_largos |> 
    filter(!is_null_model(arima)) |> 
    slice_max(order_by = .model, n = 1, by = .id) |> 
    dplyr::select(.id, arima)
  
  # Creamos el objeto fable
  bind_rows(
    ajustes_final |> filter(.id == 1) |> m_int_sol2_fc(h = 12),
    ajustes_final |> filter(.id == 2) |> m_int_sol2_fc(h = 12),
    ajustes_final |> filter(.id == 3) |> m_int_sol2_fc(h = 12),
    ajustes_final |> filter(.id == 4) |> m_int_sol2_fc(h = 12),
    ajustes_final |> filter(.id == 5) |> m_int_sol2_fc(h = 12),
    ajustes_final |> filter(.id == 6) |> m_int_sol2_fc(h = 12)
  ) |> 
    as_fable(
      index = mes_anio, key = c(.id, .model), response = "pasajeros", distribution = pasajeros
    ) |> 
    mutate(.id = paste("Pronóstico de", .id + 2019))
}

m_int_sol2_fc <- function(fit, h = h){
  datos_entrenamiento <- fit$arima[[1]]$data
  id <- fit$.id[1]
  datos_nuevos <- new_data(datos_entrenamiento, n = h) |> 
    mutate(pandemia = if_else((mes_anio >= yearmonth("2020 mar") & mes_anio <= yearmonth("2021 jul")), 1, 0),
           post_pandemia = if_else(mes_anio > yearmonth("2021 jul"), 1, 0))
  
  fit |> 
    dplyr::select(-.id) |> 
    dplyr::select(arima) |> 
    forecast(new_data = datos_nuevos) |> 
    mutate(.id = id)
}

stm_sol3 <- function(df){
  df |> 
    stretch_tsibble(.step = 12, .init = 168) |>
    mutate(
      pasajeros = if_else((mes_anio >= yearmonth("2020 mar") & mes_anio <= yearmonth("2021 jul")), NA_real_, pasajeros)
    ) |> 
    model(
      arima = ARIMA(log(pasajeros))
    ) |> 
    forecast(h = 12) |> 
    mutate(.id = paste("Pronóstico de", .id + 2019))
}

wmhb_sol4 <- function(df){
  # Se reconstruye el período de mar 2020 - nove 2022 con el promedio de los últimos 3 años
  promedio_3a <- df |>
    filter(
      mes_anio >= yearmonth("2020 March") - 3 * 12,
      mes_anio <= yearmonth("2020 Feb")
    ) |>
    as_tibble() |>
    summarise(ave = mean(pasajeros), .by = "Mes")
  
  # Se completa el dataframe con la imputacion en el periodo de la pandemia
  datos_series_wmhb <- df |> 
    left_join(promedio_3a, by = "Mes") |> 
    mutate(pasajeros = if_else((mes_anio >= yearmonth("2020 Mar") & mes_anio <= yearmonth("2021 Jul")),
                              ave, pasajeros)) |> 
    dplyr::select(-ave)
  
  # Se ajusta el modelo en los distintos períodos
  datos_series_wmhb |> 
    stretch_tsibble(.step = 12, .init = 168) |> 
    model(
      arima = ARIMA(log(pasajeros))
    ) |> 
    forecast(h = 12) |> 
    mutate(.id = paste("Pronóstico de", .id + 2019))
}

ensemble_sol5 <- function(object) {
  lapply(object, function(x) {
    x |> as_tibble() |> mutate(.model = NULL)
  }) |> 
    bind_rows() |> 
    mutate(pasajeros = generate(pasajeros, 2000)) |> 
    group_by(.id, mes_anio) |> 
    summarise(pasajeros = distributional::dist_sample(list(c(unlist(pasajeros))))) |> 
    ungroup() |> 
    as_fable(index = mes_anio, key = .id,
             response = "Pasajeros", distribution = pasajeros) |> 
    suppressWarnings()
}
