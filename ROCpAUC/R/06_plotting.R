# ==============================================================================
# 06_plotting.R
# ==============================================================================
# Visualisation helpers for ROC curves and bootstrap diagnostics.
#
# Functions in this file:
#   plot_ROC_curves         - ROC curve overlay and/or individual plots with
#                             shaded partial-area polygons.
#   plot_bootstrap_results  - CI panel and bootstrap distribution boxplots.
#
# Both functions are internal; they are called by the exported *AUC and *boot
# functions when plot = TRUE. They use ggplot2 and gridExtra throughout.
# ==============================================================================





# ------------------------------------------------------------------------------
# plot_ROC_curves
# ------------------------------------------------------------------------------
# Descripcion:
#   Representa las curvas ROC para multiples variables predictoras y resalta
#   la region de area parcial de interes con un poligono sombreado.
#   Admite modos de visualizacion combinado (superposicion) e individual (paginado).
#   Para graficos individuales incluye:
#     - Dominio TVP (axis="y"): lineas de referencia NLR0, TPR0 y FPR0,
#       etiqueta semantica del tipo de curva (BpNLR / Proper / Improper).
#     - Dominio TFP (axis="x"): lineas horizontales en TPR1 y TPR2 interpolados,
#       etiqueta semantica del tipo de curva (Bounded TPR / Proper / Improper).
#   Para el grafico combinado la leyenda incluye el tipo de cada curva.
#
# Parametros:
#   results        - Lista nombrada con $name.variable (caracter) y
#                    $results (lista de listas, cada una con una matriz $ROC_points)
#   low.value      - Limite inferior de la region de interes (TFP o TVP)
#   up.value       - Limite superior de la region de interes (TFP o TVP)
#   plot_type      - "combined": un grafico de superposicion; "individual": uno por variable;
#                    "both": ambos tipos de salida
#   axis           - "x" para dominio TFP (corte vertical);
#                    "y" para dominio TVP (corte horizontal)
#   plots_per_page - Numero maximo de paneles individuales por pagina
# ------------------------------------------------------------------------------

plot_ROC_curves = function(results, low.value, up.value,
                           plot_type      = c("combined", "individual", "both"),
                           axis           = c("x", "y"),
                           plots_per_page = 2) {
  
  plot_type = match.arg(plot_type)
  axis      = match.arg(axis)
  
  # Construimos un unico data.frame largo con todos los puntos ROC y una columna Variable
  df_all = do.call(rbind, lapply(seq_along(results$name.variable), function(i) {
    data.frame(
      FPR      = results$results[[i]]$ROC_points[, 1],
      TPR      = results$results[[i]]$ROC_points[, 2],
      Variable = results$name.variable[i]
    )
  }))
  
  # Un color por variable (coherente entre los graficos combinado e individual)
  palette = scales::hue_pal()(length(unique(df_all$Variable)))
  names(palette) = unique(df_all$Variable)
  
  # -- Auxiliar: construir un poligono de sombreado cerrado para una variable --
  # El poligono sigue la curva ROC en la region de interes (borde superior)
  # y se cierra a lo largo de la linea base (borde inferior: TVP = 0 para dominio TFP,
  # TVP = low.value para dominio TVP).
  make_shade = function(df_v) {
    fpr_roc = df_v$FPR
    tpr_roc = df_v$TPR
    
    if (axis == "x") {
      slice = tryCatch(
        portion_ROC_FPR(up.value, low.value, fpr_roc, tpr_roc),
        error = function(e) NULL
      )
      if (is.null(slice) || nrow(slice) < 2) return(NULL)
      fpr_s = slice[, 1]
      tpr_s = slice[, 2]
      # Cerramos hacia abajo hasta TVP = 0 (linea base para el pAUC del dominio TFP)
      return(data.frame(
        FPR = c(fpr_s, rev(fpr_s)),
        TPR = c(tpr_s, rep(0, length(tpr_s)))
      ))
      
    } else {
      slice = tryCatch(
        portion_ROC_TPR(up.value, low.value, fpr_roc, tpr_roc),
        error = function(e) NULL
      )
      if (is.null(slice) || nrow(slice) < 2) return(NULL)
      fpr_s = slice[, 1]
      tpr_s = slice[, 2]
      # Cerramos hacia abajo hasta TVP = low.value (parte inferior de la banda horizontal)
      return(data.frame(
        FPR = c(fpr_s,             rev(fpr_s)),
        TPR = c(tpr_s,             rep(low.value, length(tpr_s)))
      ))
    }
  }
  
  # -- Clasificacion del tipo de curva en la region de interes ------------------
  # Para axis="x": usa classification_Tp -> "Proper", "Bounded TPR" o "Improper"
  # Para axis="y": usa portion_ROC_TPR + shapepROC, identico a FpA, para garantizar
  #   que la frontera inferior se interpola exactamente en TPR0. El filtrado crudo
  #   por indice no interpola la frontera, produciendo un NLR0 incorrecto y
  #   clasificaciones erroneas (p. ej. SPNCRNA.1080 aparecia como "Proper").
  classify_curve = function(fpr_v, tpr_v) {
    if (axis == "x") {
      idx = fpr_v >= low.value & fpr_v <= up.value
      if (sum(idx) < 2) return("unknown")
      type = classification_Tp(fpr_v[idx], tpr_v[idx])
      if      (type[1]) return("Bounded TPR")
      else if (type[2]) return("Proper")
      else              return("Improper")
    } else {
      # Replicamos exactamente el pipeline de NpA (linea 2516 de NpAUC exportada):
      #   portion = portion_ROC_TPR(up.value, low.value, roc[,1], roc[,2])
      #   NpA(portion[,2], portion[,1])   <- tpr.proc=col2, fpr.proc=col1
      # Y dentro de NpA se llama classification_Fp(tpr.proc, fpr.proc).
      # Usar shapepROC era incorrecto porque shapepROC elimina el ultimo punto
      # mientras classification_Fp filtra is.finite(NLR), lo que produce
      # clasificaciones distintas para curvas impropias cerca de (1,1).
      portion = portion_ROC_TPR(up.value, low.value, fpr_v, tpr_v)
      fpr_s   = portion[, 1]
      tpr_s   = portion[, 2]
      if (length(fpr_s) < 2) return("unknown")
      # classification_Fp(tpr.proc, fpr.proc): mismo orden que NpA
      tipo = classification_Fp(tpr_s, fpr_s)
      if      (tipo[1]) return("BpNLR")    # bounded: RVN acotada por RVN0
      else if (tipo[3]) return("Improper") # improper: RVN > 1 en algun punto
      else              return("pProp")    # partial: RVN <= 1 pero > RVN0 (NpAUC = NA)
    }
  }
  
  curve_types = setNames(
    vapply(unique(df_all$Variable), function(vn) {
      dv = df_all[df_all$Variable == vn, ]
      classify_curve(dv$FPR, dv$TPR)
    }, character(1)),
    unique(df_all$Variable)
  )
  
  # Color semantico de la etiqueta de tipo:
  #   verde  (#2E7D32) -> Proper / Bounded TPR / BpNLR
  #   naranja (#E65100) -> Improper
  #   gris   (#607D8B) -> indeterminado
  label_colour = function(tipo) {
    switch(tipo,
           "Proper"      = "#2E7D32",
           "Bounded TPR" = "#2E7D32",
           "BpNLR"       = "#2E7D32",
           "pProp"       = "#1565C0",   # azul -> parcialmente propia (NpAUC indefinido)
           "Improper"    = "#E65100",
           "#607D8B"
    )
  }
  
  # -- Auxiliar: interpolar TPR en un FPR dado ----------------------------------
  interp_tpr = function(fpr_v, tpr_v, fpr_target) {
    if (fpr_target <= min(fpr_v)) return(tpr_v[which.min(fpr_v)])
    if (fpr_target >= max(fpr_v)) return(tpr_v[which.max(fpr_v)])
    i_up  = min(which(fpr_v >= fpr_target))
    i_low = max(which(fpr_v <= fpr_target))
    if (i_up == i_low) return(tpr_v[i_up])
    t = (fpr_target - fpr_v[i_low]) / (fpr_v[i_up] - fpr_v[i_low])
    tpr_v[i_low] + t * (tpr_v[i_up] - tpr_v[i_low])
  }
  
  # -- Auxiliar: interpolar FPR en un TPR dado ----------------------------------
  # Cuando hay empate exacto en tpr_target (varios puntos con TPR == tpr_target),
  # se elige el punto con menor FPR (el más cercano al origen en el eje X).
  # Esto garantiza que fpr0 sea el punto de corte correcto con la frontera
  # horizontal y evita división 0/0 (NaN) en el cálculo de la línea NLR0.
  interp_fpr = function(fpr_v, tpr_v, tpr_target) {
    if (tpr_target <= min(tpr_v)) return(fpr_v[which.min(tpr_v)])
    if (tpr_target >= max(tpr_v)) return(fpr_v[which.max(tpr_v)])
    # Empate exacto: devolvemos el punto de menor FPR (más cercano al origen)
    exact = which(tpr_v == tpr_target)
    if (length(exact) > 0) return(min(fpr_v[exact]))
    i_up  = min(which(tpr_v >  tpr_target))   # estrictamente mayor
    i_low = max(which(tpr_v <  tpr_target))   # estrictamente menor
    t = (tpr_target - tpr_v[i_low]) / (tpr_v[i_up] - tpr_v[i_low])
    fpr_v[i_low] + t * (fpr_v[i_up] - fpr_v[i_low])
  }
  
  # -- Graficos individuales (uno por variable, paginados) ----------------------
  plot_list = lapply(unique(df_all$Variable), function(vn) {
    dv    = df_all[df_all$Variable == vn, ]
    shade = make_shade(dv)
    col   = palette[vn]
    tipo  = curve_types[[vn]]
    lcol  = label_colour(tipo)
    
    p = ggplot2::ggplot(dv, ggplot2::aes(x = FPR, y = TPR)) +
      ggplot2::geom_abline(slope = 1, intercept = 0,
                           linetype = "dashed", colour = "grey60") +
      ggplot2::geom_line(colour = col, linewidth = 1.2) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::labs(
        title = paste("ROC Curve:\n", vn),
        x     = "False Positive Ratio (FPR)",
        y     = "True Positive Ratio (TPR)"
      ) +
      ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
    
    # Poligono sombreado del area parcial (color de la curva)
    if (!is.null(shade))
      p = p + ggplot2::geom_polygon(
        data    = shade,
        mapping = ggplot2::aes(x = FPR, y = TPR),
        fill = col, alpha = 0.22, inherit.aes = FALSE
      )
    
    if (axis == "y") {
      # Dominio TVP: lineas NLR0, TPR0 y FPR0 (sin etiquetas)
      fpr0  = interp_fpr(dv$FPR, dv$TPR, low.value)
      nlr0  = (1 - low.value) / (1 - fpr0)
      tpr_nlr_end = min(max(low.value + nlr0 * (1 - fpr0), 0), 1)
      nlr_line = data.frame(FPR = c(fpr0, 1), TPR = c(low.value, tpr_nlr_end))
      
      p = p +
        # Frontera inferior de la banda: linea horizontal en TPR0
        ggplot2::geom_hline(yintercept = low.value,
                            linetype = "dashed", colour = "grey40", linewidth = 0.7) +
        # Frontera vertical en FPR0
        ggplot2::geom_vline(xintercept = fpr0,
                            linetype = "dashed", colour = "grey40", linewidth = 0.7) +
        # Segmento vertical en x=1 representando 1-FPR0
        ggplot2::annotate("segment",
                          x = 1, xend = 1, y = low.value, yend = 1,
                          linetype = "solid", colour = "grey50", linewidth = 0.6) +
        # Linea NLR0 desde (FPR0, TPR0) hasta el borde derecho
        ggplot2::geom_line(
          data    = nlr_line,
          mapping = ggplot2::aes(x = FPR, y = TPR),
          linetype = "dotted", colour = "grey30", linewidth = 0.8,
          inherit.aes = FALSE
        )
      
    } else {
      # Dominio TFP: lineas horizontales en TPR1 y TPR2 (sin etiquetas)
      tpr1 = interp_tpr(dv$FPR, dv$TPR, low.value)
      tpr2 = interp_tpr(dv$FPR, dv$TPR, up.value)
      
      p = p +
        # Fronteras verticales de la banda: FPR1 y FPR2
        ggplot2::geom_vline(xintercept = c(low.value, up.value),
                            linetype = "dashed", colour = "grey40", linewidth = 0.7) +
        # Lineas horizontales en TPR1 y TPR2
        ggplot2::geom_hline(yintercept = tpr1,
                            linetype = "dotted", colour = "grey40", linewidth = 0.7) +
        ggplot2::geom_hline(yintercept = tpr2,
                            linetype = "dotted", colour = "grey40", linewidth = 0.7)
    }
    
    # Etiqueta del tipo de curva en esquina inferior derecha con color semántico
    p = p + ggplot2::annotate(
      "label",
      x          = 0.97, y = 0.03,
      label      = paste("Type:", tipo),
      hjust      = 1,    vjust = 0,
      size       = 3,
      fill       = lcol,
      colour     = "white",
      fontface   = "bold",
      alpha      = 0.85
    )
    
    return(p)
  })
  
  # -- Grafico de superposicion combinado ---------------------------------------
  if (plot_type %in% c("combined", "both")) {
    
    # Etiquetas de leyenda enriquecidas: "variable (tipo en region)"
    legend_labels = setNames(
      paste0(names(curve_types), " (", curve_types, ")"),
      names(curve_types)
    )
    
    p_comb = ggplot2::ggplot(
      df_all, ggplot2::aes(x = FPR, y = TPR, colour = Variable)
    ) +
      ggplot2::geom_abline(slope = 1, intercept = 0,
                           linetype = "dashed", colour = "grey60") +
      ggplot2::geom_line(linewidth = 1.1) +
      ggplot2::scale_colour_manual(values = palette, labels = legend_labels) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::labs(
        title  = "ROC Curves \u2014 All Variables",
        x      = "False Positive Ratio (FPR)",
        y      = "True Positive Ratio (TPR)",
        colour = "Variable (type in region)"
      ) +
      ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
    
    if (axis == "x")
      p_comb = p_comb +
      ggplot2::geom_vline(xintercept = c(low.value, up.value), linetype = "dotted")
    else
      p_comb = p_comb +
      ggplot2::geom_hline(yintercept = c(low.value, up.value), linetype = "dotted")
    
    print(p_comb)
  }
  
  # -- Graficos individuales paginados ------------------------------------------
  if (plot_type %in% c("individual", "both")) {
    for (j in seq(1, length(plot_list), by = plots_per_page)) {
      grp = plot_list[j:min(j + plots_per_page - 1L, length(plot_list))]
      gridExtra::grid.arrange(grobs = grp, ncol = min(plots_per_page, 3L))
    }
  }
  
  invisible(NULL)
}


# ------------------------------------------------------------------------------
# plot_bootstrap_results
# ------------------------------------------------------------------------------
# Descripcion:
#   Genera graficos de diagnostico para los resultados de los intervalos de confianza
#   bootstrap. Hay dos tipos de visualizacion disponibles (seleccionados mediante plot_type):
#
#   "ci"           : Grafico de intervalos de confianza que muestra la estimacion original
#                    (circulo naranja), la media bootstrap (triangulo azul) y
#                    los bigotes del intervalo de confianza. La diferencia entre la
#                    estimacion original y la media bootstrap es el sesgo bootstrap,
#                    que debe ser pequeno en relacion con el ancho del IC.
#
#   "distribution" : Diagrama de caja de la distribucion bootstrap para cada variable,
#                    superpuesto con valores individuales de replica en jitter y la
#                    estimacion original como un diamante naranja.
#
#   "both"         : Representa ambos graficos (por defecto).
#
# Parametros:
#   result_boot   - Objeto devuelto por boot::boot()
#   name.variable - Vector de caracteres con los nombres de los predictores
#   index_label   - Nombre del indice (p.ej., "MCpAUC"), usado en las etiquetas de los ejes
#   lwr_vals      - Vector numerico de limites inferiores del IC (uno por variable)
#   upr_vals      - Vector numerico de limites superiores del IC (uno por variable)
#   plot_type     - Uno de "both", "ci", "distribution"
# ------------------------------------------------------------------------------

plot_bootstrap_results = function(result_boot, name.variable, index_label,
                                  lwr_vals, upr_vals,
                                  plot_type = c("both", "ci", "distribution")) {
  
  plot_type = match.arg(plot_type)
  
  # Solo representamos las variables que tienen una estimacion puntual finita y al menos
  # algunas replicas bootstrap no NA
  valid = which(
    !is.na(result_boot$t0) &
      !apply(result_boot$t, 2, function(col) all(is.na(col) | is.nan(col)))
  )
  if (length(valid) == 0) {
    message("No hay variables validas para representar.")
    return(invisible(NULL))
  }
  
  # -- Grafico de IC + sesgo ----------------------------------------------------
  if (plot_type %in% c("both", "ci")) {
    t_means = colMeans(result_boot$t[, valid, drop = FALSE], na.rm = TRUE)
    df_ic = data.frame(
      Variable = factor(name.variable[valid], levels = name.variable[valid]),
      Estimate = result_boot$t0[valid],
      BootMean = t_means,
      bias     = t_means - result_boot$t0[valid],
      lwr      = lwr_vals[valid],
      upr      = upr_vals[valid]
    )
    
    p_ci = ggplot2::ggplot(df_ic, ggplot2::aes(x = Variable)) +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = lwr, ymax = upr),
        width = 0.25, colour = "grey40"
      ) +
      ggplot2::geom_point(
        ggplot2::aes(y = Estimate, shape = "Original estimate"),
        size = 3.5, colour = "#D55E00"
      ) +
      ggplot2::geom_point(
        ggplot2::aes(y = BootMean, shape = "Bootstrap mean"),
        size = 3.5, colour = "#0072B2"
      ) +
      ggplot2::scale_shape_manual(
        name   = NULL,
        values = c("Original estimate" = 16, "Bootstrap mean" = 17)
      ) +
      ggplot2::scale_y_continuous(limits = c(0, 1)) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(legend.position = "bottom") +
      ggplot2::labs(
        title    = paste("Bootstrap Confidence Intervals \u2014", index_label),
        subtitle = "Orange circle = original estimate  |  Blue triangle = bootstrap mean",
        y        = index_label,
        x        = NULL
      )
    print(p_ci)
  }
  
  # -- Diagrama de caja de la distribucion bootstrap ----------------------------
  if (plot_type %in% c("both", "distribution")) {
    df_dist = as.data.frame(result_boot$t[, valid, drop = FALSE])
    colnames(df_dist) = name.variable[valid]
    
    df_long = tidyr::pivot_longer(
      df_dist,
      cols      = tidyselect::everything(),
      names_to  = "Variable",
      values_to = index_label
    )
    df_long$Variable = factor(df_long$Variable, levels = name.variable[valid])
    
    # Estimaciones puntuales originales como diamantes naranjas sobre los diagramas de caja
    df_t0 = data.frame(
      Variable = factor(name.variable[valid], levels = name.variable[valid]),
      t0       = result_boot$t0[valid]
    )
    
    print(
      ggplot2::ggplot(df_long, ggplot2::aes(
        x = .data[["Variable"]], y = .data[[index_label]]
      )) +
        ggplot2::geom_boxplot(outlier.shape = NA, fill = "#56B4E9", alpha = 0.7) +
        ggplot2::geom_jitter(width = 0.15, alpha = 0.35, size = 0.9) +
        ggplot2::geom_point(
          data        = df_t0,
          mapping     = ggplot2::aes(x = Variable, y = t0),
          shape       = 18, size = 4.5, colour = "#D55E00",
          inherit.aes = FALSE
        ) +
        ggplot2::scale_y_continuous(limits = c(0, 1)) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::labs(
          title    = paste("Bootstrap Distribution \u2014", index_label),
          subtitle = "Orange diamond = original estimate (t(0))",
          y        = paste(index_label, "(bootstrap)"),
          x        = NULL
        )
    )
  }
  
  invisible(NULL)
}
