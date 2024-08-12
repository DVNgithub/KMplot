###--------------------------------
#--- importing packages
###--------------------------------
usethis::use_package("survival")
usethis::use_package("ggpubr")
usethis::use_package("ggplot2")
usethis::use_package("survminer")
usethis::use_package("gridExtra")
usethis::use_package("grid")

#' KM plot with only two subgroups
#'
#' @param dat data plot
#' @param trt treatments/arms
#' @param tte time to event variable
#' @param ttecen censor variable
#' @param group subgroup variable
#' @param strat whether using stratification
#' @param ttime.inc incremate interval time
#' @param switch.leg just for special case, in which legends switching
#' @param censor.status whether censor status displaying
#' @param legend.type for different legends format to display
#' @param level.grp subgroup levels
#' @param legend.lb legend labs
#' @param legend.title.h legend titles
#' @param cols KM curves colors
#' @param legend.pos legend position
#' @param legend.adj legend adjusted position
#' @param expand.val for expand x-axis
#' @param legend.title size of legend title
#' @param legend.text size of legend text
#' @param axis.title size of axis title
#' @param axis.text size of axis text
#' @param yname y-axis lab
#' @param margin.plot margins of plot
#' @param margin.tb margings of risk of table
#' @param xlim.lo x-axis limit
#' @param xminlab length of risk table legend
#' @param xmaxlab spaces between risk table legend and table
#' @param heights.plot the heights of plot
#' @param n.round statistics round number
#' @param p.round statistics round number
#' @param n.space the lists for adjusting spaces between values of legends
#'
#' @return KM plot with two curves
#'
#' @export
#'
#' @examples
#' # start example
#' attach(datplot)
#'
#' library(survival)
#' library(survminer)
#' library(dplyr)
#' library(ggplot2)
#' library(gridExtra)
#' library(grid)
#'
#' df.plot <- mutate(datplot, grp2=trt)
#'
#' p1 <- KMplot.onlytwogrp(dat=df.plot,
#'                          tte="timevar",
#'                          ttecen="cenvar",
#'                          group = "grp2",
#'                          yname="Invasive Disease-Free Survival",
#'                          ttime.inc=3, legend.type=3,
#'                          legend.title.h="           Arm                 Median      HR (95% CI)      p-value",
#'                          level.grp=c("Arm1", "Arm2"),
#'                          legend.lb=c('Arm1', 'Arm2'),
#'                          cols=c('red', "blue"),
#'                          legend.pos = c(0, 0),
#'                          legend.adj = c(-0.07, -0.3),
#'                          legend.title=11, legend.text=10,
#'                          axis.title=11, axis.text=10,
#'                          margin.plot=c(2, 2.5, 0, 1.5),
#'                          margin.tb=c(0, 2.1, 3, 0.7),
#'                          xminlab=-2.5, xmaxlab =-1, xlim.lo = -0.5,
#'                          heights.plot = c(2, 0.8),
#'                          n.space = list(c(17,16,17, 17,13)))
#' # different format
#' p2 <- KMplot.onlytwogrp(dat=df.plot,
#'                        tte="timevar",
#'                        ttecen="cenvar",
#'                        group = "grp2",
#'                        legend.type=1,
#'                        yname="Invasive Disease-Free Survival",
#'                        ttime.inc=3,
#'                        legend.title.h="           Arm                 Median      HR (95% CI)",
#'                        level.grp=c("Arm1", "Arm2"),
#'                        legend.lb=c('Arm1', 'Arm2'),
#'                        cols=c('red', "blue"),
#'                        legend.pos = c(0, 0),
#'                        legend.adj = c(-0.07, -0.3),
#'                        legend.title=11, legend.text=10,
#'                        axis.title=11, axis.text=10,
#'                        margin.plot=c(2, 2.5, 0, 1.5),
#'                        margin.tb=c(0, 2.1, 3, 0.7),
#'                        xminlab=-2.5, xmaxlab =-1, xlim.lo = -0.5,
#'                        heights.plot = c(2, 0.8),
#'                        n.space = list(c(17, 16, 17, 13)))
#'
#' # output
#' pdf(file=file.path(getwd(), "f_km_plot_2grp.pdf"), width=11, height=8.5)
#' grid.newpage()
#' print(p1, vp = viewport(x = unit(0.50, "npc"), y = unit(0.52, "npc"), width=0.90, height=0.85, just = "centre"))
#' grid.newpage()
#' print(p2, vp = viewport(x = unit(0.50, "npc"), y = unit(0.52, "npc"), width=0.90, height=0.85, just = "centre"))
#' dev.off()
#' #' # end example
KMplot.onlytwogrp <- function(dat,
                              trt = "trt",                              # treatments/arms variable
                              tte, ttecen, group,                       # TTE, censor, and the groups for comparison
                              strat = NA,
                              ttime.inc=4,                              # TTE break interval
                              switch.leg=TRUE,                          # for switching legends of KMs
                              censor.status=FALSE,
                              legend.type=1,
                              level.grp=c("Abema+EDT", "EDT"),
                              legend.lb=c('Abema+EDT', 'EDT'),
                              legend.title.h="Legend head",
                              cols=c('red', 'blue', 'orange', 'green'),
                              legend.pos = c(0, 0),
                              legend.adj = c(-0.1, -0.3),
                              expand.val = c(0, 0),
                              legend.title=11, legend.text=10,
                              axis.title=11, axis.text=10,
                              yname="",
                              margin.plot=c(1, 1, 1, 1),
                              margin.tb=c(1, 1, 1, 1),
                              xlim.lo = 0,
                              xminlab=-2,                               # length of risk table legend
                              xmaxlab = -1,                             # space between risk table legend and table
                              heights.plot = c(1.7, 0.75),
                              n.round=2, p.round=3,
                              n.space=list(c(14,5,7,15,5,7))) {

  names(dat)[match(trt, names(dat))] <- "trt"
  names(dat)[match(tte, names(dat))] <- "time"
  names(dat)[match(ttecen, names(dat))] <- "ttecen"
  names(dat)[match(group, names(dat))] <- "group"

  dat.use <- mutate(dat,
                    GRP = factor(group, levels = level.grp),
                    time = time,
                    trt = trt,
                    evnt = 1 - ttecen)

  if(switch.leg) {
    dat.lv <- mutate(dat.use,
                     GRP=factor(group, levels = level.grp))
  } else {
    dat.lv <- dat.use
  }
  fit.plot <- survfit(Surv(time, evnt) ~ GRP, data = dat.lv, conf.type = 'log-log')
  x <- data.frame(summary(fit.plot)$table)
  meds <- quantile(fit.plot)

  if (is.na(strat)){
    hr1 <- coxph(Surv(time, evnt) ~ trt, data = dat.use, ties = 'exact')
  } else {
    form1 <- paste0('Surv(time, evnt) ~ trt + ', strat)
    hr1 <- coxph(as.formula(form1), data = dat.use, ties = 'exact')
  }
  hr <- format_num(summary(hr1)$conf.int[1], n.round)
  ci <- paste0('(', format_num(summary(hr1)$conf.int[3], n.round), ', ', format_num(summary(hr1)$conf.int[4], n.round), ')')

  # get Log−rank p-value
  p.val <- ifelse(summary(hr1)[["coefficients"]][,"Pr(>|z|)"] < 0.0001, "<.0001",
                  round(summary(hr1)[["coefficients"]][,"Pr(>|z|)"], p.round))

  legend.lb.txt <- tte.legend.form(x=x, subgrpComp=list(level.grp),
                                   hr=hr, ci=ci, pval=p.val,
                                   legend.type=legend.type,
                                   n.space = n.space)

  p <- ggsurvplot(fit.plot, data = dat.use,
                  size = 1,
                  censor=censor.status, censor.shape = "|", censor.size = 4,
                  palette = cols, conf.int = F,
                  #linetype = ltypes,
                  risk.table = TRUE, fontsize = 4, tables.theme = clean_theme(),
                  xlab = "Time (months)", ylab = yname,
                  surv.scale = "percent",
                  break.time.by = ttime.inc, break.y.by = 0.1,
                  xlim=c(xlim.lo, (max(dat$time)+ttime.inc)),
                  risk.table.y.text =FALSE,
                  legend.title = legend.title.h,
                  legend.labs = legend.lb.txt)

  p$plot <- p$plot +
    scale_x_continuous(breaks = seq(0, (max(dat.use$time)+ttime.inc), ttime.inc), expand = expand.val) +
    scale_y_continuous(breaks=seq(0, 10, 0.1), limits = c(0, 1),
                       labels = scales::percent_format(accuracy=5L)) +
    geom_hline(yintercept = .5, color = 'grey28', linetype = 'dashed') +
    theme(legend.position = legend.pos,
          legend.justification = legend.adj,
          legend.key.width=unit(1.2, "cm"),
          legend.title = element_text(size = legend.title, face = "bold"),
          legend.text = element_text(size = legend.text, face = "bold"),
          axis.title.y=element_text(size = axis.title, face = "bold"),
          axis.title.x=element_text(size = axis.title, face = "bold"),
          axis.text.x =element_text(size = axis.text, face = "bold"),
          axis.text.y =element_text(size = axis.text, vjust=0.5, hjust=0.5, face = "bold"),
          axis.ticks=element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          text=element_text(size = axis.title, face = "bold"),
          plot.margin = unit(margin.plot, "cm"))


  # update risk table
  p$table <- p$table +
    ggtitle("Patients at risk") +
    theme(plot.title = element_text(hjust = 0.05, size = 11),
          plot.margin = unit(margin.tb, "cm"),
          axis.text.y=element_blank(),
          axis.title.y=element_blank())

  for(i in 1:length(legend.lb)) {
    p$table <- p$table +
      annotate("segment",
               x = xminlab, xend = xmaxlab, y = i, yend = i,
               colour = rev(cols)[i], linewidth = 1) +
      coord_cartesian(clip = "off")
  }

  p.out <- ggarrange(p$plot, p$table, heights = heights.plot, ncol = 1, nrow = 2)

  return(p.out)
}#end



#' KM plot with more than two subgroups
#'
#' @param dat data plot
#' @param trt treatments/arms
#' @param tte time to event variable
#' @param ttecen censor variable
#' @param group subgroup variable
#' @param strat whether using stratification
#' @param ttime.inc incremate interval time
#' @param switch.leg just for special case, in which legends switching
#' @param level.grp subgroup levels
#' @param subgrpName subgroup names
#' @param subgrpComp subgroups for comparison
#' @param censor.status whether censor status displaying
#' @param censor.size size of censor
#' @param titles title of plot
#' @param legend.type for different legends format to display
#' @param legend.head legend titles
#' @param legend.pos legend position
#' @param legend.adj adjusted legend position
#' @param expand.val for x-axis expanding
#' @param font.legend.txt size of legend text
#' @param font.legend.title size of legend title
#' @param color.val KM curves colors
#' @param linetype KM curve line types
#' @param ylab.use y-axis lab
#' @param y.seq y-axis sequence
#' @param y.limit y-axis limit
#' @param xlim.lo x-axis limit
#' @param xminlab length of risk table legend
#' @param xmaxlab spaces between risk table legend and table
#' @param special.case specical case for plotting risk lines legend
#' @param hjust.val spaces adjust for risk lines of legends if special.case is TRUE
#' @param axis.text size of axis text
#' @param axis.title size of axis title
#' @param font.tb size of text of table
#' @param n.round statistics round number
#' @param p.round statistics round number
#' @param heights.plot heights of plot
#' @param margin.plot margins of plot
#' @param margin.tb margins of table
#' @param n.space the lists for adjusting spaces between values of legends
#'
#' @return KM curves
#'
#' @export
#'
#' @examples
#' # start example
#' attach(datplot)
#'
#' library(survival)
#' library(survminer)
#' library(dplyr)
#' library(ggplot2)
#' library(gridExtra)
#' library(grid)
#'
#' df.plot <- mutate(datplot2,
#'                   timevar=IDFS1_TIME, cenvar=IDFS1_CNSR,
#'                   grp=paste(grp.pos, trt, sep="_"),
#'                   trt=factor(trt, levels = c("EDT", "Abema+EDT"))) %>%
#'            select(-group)
#'
#' p3 <- KMplot.w.risk(dat=df.plot,
#'                    trt = "trt",
#'                    tte="timevar", ttecen="cenvar", group="grp",
#'                    ttime.inc = 3,
#'                    subgrpName=c("ER+/PR+_Abema+EDT", "ER+/PR+_EDT",
#'                                 "ER/PR_other_Abema+EDT", "ER/PR_other_EDT"),
#'                    subgrpComp=list(c("ER+/PR+_Abema+EDT", "ER+/PR+_EDT"),
#'                                    c("ER/PR_other_Abema+EDT", "ER/PR_other_EDT")),
#'                    censor.status = FALSE,
#'                    legend.type = 2,
#'                    legend.head="               Patient Group                          HR (95% CI)     P-value",
#'                    legend.pos = c(0, 0),
#'                    legend.adj = c(-0.05, -0.1),
#'                    font.legend.txt = 10,
#'                    font.legend.title = 11,
#'                    color.val=c('red', 'blue', 'red', 'blue'),
#'                    linetype = c('solid','solid','dashed','dashed'),
#'                    #linetype=c(1,1,2,2),
#'                    ylab.use="Progression-Free Survival",
#'                    y.seq = seq(0, 10, 0.1),
#'                    y.limit = c(0, 1),
#'                    xminlab = -5,
#'                    xmaxlab = -2.2,
#'                    xlim.lo = -1,
#'                    axis.text = 12, axis.title = 13,
#'                    font.tb = 4,
#'                    n.round=3, p.round = 4,
#'                    special.case = TRUE,
#'                    heights.plot = c(2, 0.85),
#'                    margin.plot=c(1, 1, 0, 1),
#'                    margin.tb=c(0, 0.4, 2, 0.7),
#'                    n.space=list(c(20, 15, 29),
#'                                 c(13, 15, 24)))
#' # output
#' pdf(file=file.path(getwd(), "f_km_plot_pair.comp.pdf"), width=11, height=8.5)
#' grid.newpage()
#' print(p3, vp = viewport(x = unit(0.52, "npc"), y = unit(0.54, "npc"), width=0.90, height=0.9, just = "centre"))
#' dev.off()
#' #' # end example
KMplot.w.risk <- function(dat,
                          trt = "trt",                              # treatments/arms variable
                          tte, ttecen, group,                       # TTE, censor, and the groups for comparison
                          strat = NA,
                          ttime.inc = 4,                            # TTE break interval
                          switch.leg = FALSE,                       # for switching legends of KMs
                          level.grp = NULL,                         # if need to switch legend
                          subgrpName = c("grp1", "grp2", "grp3", "grp4"),           # these are from "group"
                          subgrpComp = list(c("grp1", "grp2"), c("grp3", "grp4")),  # paired comparison
                          censor.status = FALSE,                    # if displaying censor
                          censor.size = 1,
                          titles = "",
                          legend.type = 1,                          # type of legend to show on the plot
                          legend.head = "Patient Group  HR (95% CI)   P-value",     # legend title
                          legend.pos = c(0, 0),
                          legend.adj = c(-0.05, -0.1),
                          expand.val = c(0, 0),
                          font.legend.txt = 10,
                          font.legend.title = 11,
                          color.val = c(1, 1, 1, 1),
                          linetype = c('solid','solid','dashed','dashed'),
                          ylab.use = "",
                          y.seq = seq(0, 10, 0.1),
                          y.limit = c(0, 1),
                          xlim.lo = 0,
                          xminlab = -5,                            # length of risk table legend
                          xmaxlab = -1,                            # space between risk table legend and table
                          special.case = FALSE,                    # specical case for plotting risk lines legend
                          hjust.val = 0.05,                        # space adjust for risk lines of legends if special.case is TRUE
                          axis.text = 12, axis.title = 13,
                          font.tb = 4,
                          n.round=2, p.round = 3,
                          heights.plot = c(2, 0.85),
                          margin.plot = c(1, 1, 0.1, 1.2),
                          margin.tb = c(0.1, 0.2, 2, 2.1),
                          n.space = list(c(5, 5, 5, 5, 5),        # these lists for adjust spaces between values of legends
                                         c(5, 5, 5, 5, 5),
                                         c(5, 5, 5, 5, 5))) {

  names(dat)[match(trt, names(dat))] <- "trt"
  names(dat)[match(tte, names(dat))] <- "time"
  names(dat)[match(ttecen, names(dat))] <- "ttecen"
  names(dat)[match(group, names(dat))] <- "group"

  dat.use <- mutate(dat,
                    GRP = factor(group, levels = subgrpName),
                    time = time,
                    trt = trt,
                    evnt = 1 - ttecen)

  if(switch.leg) {
    dat.lv <- mutate(dat.use,
                     GRP=factor(group, levels = unlist(level.grp)))
  } else {
    dat.lv <- dat.use
  }

  fit.plot <- survfit(Surv(time, evnt) ~ GRP, data = dat.lv, conf.type = "log-log")
  x <- data.frame(summary(fit.plot)$table)
  meds <- quantile(fit.plot)

  #---------------------------------------------------
  ###--- compute stat to generate legends
  #---------------------------------------------------
  # for more than 2 subgroups (4, 6)
  if (TRUE) {
    if (is.na(strat)) {
      form1 <- "Surv(time, evnt) ~ trt"
    } else {
      form1 <- paste0("Surv(time, evnt) ~ trt + ", strat)
    }

    hr <- NULL; ci <- NULL; pval.logr <- NULL
    for(i in 1:length(subgrpComp)) {
      df.sub <- filter(dat.use, GRP %in% subgrpComp[[i]]) %>%
        mutate(trt=factor(GRP, levels=subgrpComp[[i]]))
      val <- coxph(as.formula(form1), data = df.sub, ties = "exact")
      hr0 <- format_num(summary(val)$conf.int[1], n.round)
      ci0 <- paste0("(", format_num(summary(val)$conf.int[3], n.round), ", ", format_num(summary(val)$conf.int[4], n.round), ")")
      # for log-rank
      fit.diff <- survdiff(Surv(time, evnt) ~ trt, data = filter(dat.use, GRP %in% subgrpComp[[i]]))
      pval0 <- pchisq(summary(val)[["sctest"]]["test"], df=1, lower.tail=FALSE)
      pval <- ifelse(pval0 < 0.0001, "<.0001",
                     format_num(pval0, p.round))
      if(is.null(hr)) hr <- hr0 else hr <- c(hr, hr0)
      if(is.null(ci)) ci <- ci0 else ci <- c(ci, ci0)
      if(is.null(pval.logr)) pval.logr <- pval else pval.logr <- c(pval.logr, pval)
    }

    # create legends
    if(switch.leg) {
      legend.txt <- tte.legend.form(x=x, subgrpComp=level.grp, hr=hr, ci=ci, pval=pval.logr,
                                    legend.type=legend.type, n.space = n.space)
    } else {
      legend.txt <- tte.legend.form(x=x, subgrpComp=subgrpComp, hr=hr, ci=ci, pval=pval.logr,
                                    legend.type=legend.type, n.space = n.space)
    }
  }#end legends

  #---------------------------------------------------
  ###--- do plot
  #---------------------------------------------------
  p <- ggsurvplot(fit.plot,
                  data = dat.use, size = 1, censor.shape = "|", censor.size = censor.size,
                  conf.int = F,
                  linetype=linetype,
                  risk.table = TRUE,
                  fontsize = font.tb, tables.theme = clean_theme(),
                  xlab = "Time (months)", ylab = ylab.use, title=titles,
                  surv.scale = "percent",
                  break.time.by = ttime.inc, break.y.by = 0.1,
                  xlim = c(xlim.lo, (max(dat.use$time)+ttime.inc)),
                  risk.table.y.text =FALSE,
                  legend = legend.adj,
                  legend.title = legend.head,
                  legend.labs = legend.txt)

  p$plot <- p$plot +
    scale_color_manual(name=legend.head, values = color.val) +
    scale_linetype_manual(name=legend.head, values = linetype) +
    scale_x_continuous(breaks = seq(0, (max(dat.use$time)+ttime.inc), ttime.inc), expand = expand.val) +
    scale_y_continuous(breaks=y.seq, limits = y.limit,
                       labels = scales::percent_format(accuracy=5L)) +
    geom_hline(yintercept = .5, color = 'grey28', linetype = 'dashed') +
    theme(legend.position = legend.pos,
          legend.justification = legend.adj,
          legend.key.width=unit(1.2, "cm"),
          legend.title = element_text(size = font.legend.title, face = "bold"),
          legend.text = element_text(size = font.legend.txt, face = "bold"),
          axis.title.y=element_text(margin=margin(0,15,0,0), size = axis.title, face = "bold"),
          axis.title.x=element_text(size = axis.title, face = "bold"),
          axis.text.x =element_text(size = axis.text, face = "bold"),
          axis.text.y =element_text(size = axis.text, vjust=0.5, hjust=0.5, face = "bold"),
          #axis.ticks=element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          text=element_text(size = axis.title, face = "bold"),
          plot.margin = unit(margin.plot, "cm"))


  # update risk table
  p$table <- p$table +
    ggtitle("Patients at risk") +
    theme(plot.title = element_text(hjust = 0.05, size = 11),
          plot.margin = unit(margin.tb, "cm"),
          axis.text.y=element_blank(),
          axis.title.y=element_blank())

  if(special.case) {
    index.val <- seq(1, length(linetype), 1)
    for(i in index.val) {
      yval <- rev(index.val)[i]
      if(linetype[i] %in% "solid") {
        p$table <- p$table +
          annotate("segment",
                   x = xminlab, xend = xmaxlab, y = yval, yend = yval,
                   colour = color.val[i], linewidth = 1.1)
      } else {
        p$table <- p$table +
          annotate("text", x = xminlab, y = yval, label = "--",
                   vjust=0.5, hjust=hjust.val, size = 9, color=color.val[i], fontface = "bold") +
          coord_cartesian(clip = "off")
      }
    }
  } else {
    for(i in 1:length(legend.txt)) {
      p$table <- p$table +
        annotate("segment",
                 x = xminlab, xend = xmaxlab, y = i, yend = i,
                 colour = rev(color.val)[i], linewidth = 1.1)
    }
  }

  p.out <- ggarrange(p$plot, p$table, heights = heights.plot, ncol = 1, nrow = 2)

  return(p.out)

}#end


#' KM plot with more than two subgroups
#'
#' @param dat data plot
#' @param subcohana subcohort analysis
#' @param trt.var treatments/arms
#' @param tte.var time to event variable
#' @param cen.var censor variable
#' @param cov.add whether using covariate
#' @param cov.var covarite varibale
#' @param STRAT whether using stratification
#' @param pval.type p-value type. Default is "Wald"
#' @param COMP there are only two groups comparison
#' @param COMPGRP paired comparison
#' @param INTER whether using interaction model
#' @param pairedgrp giving paired subgrp for comparison
#' @param pairedgrplv giving paired subgrp levels
#' @param conf.level confident level
#' @param conf.type confident type
#' @param method method using for model. Default is 'efron'
#' @param ttime.inc incremate interval time
#' @param cols KM curves colors
#' @param ltypes KM curves line types
#' @param main title of plot
#' @param xname x-axis lab
#' @param yname y-axis lab
#' @param level.grp subgroup levels
#' @param legend.head legend title
#' @param legend.type for different legends format to display
#' @param legend.pos legend position
#' @param expand.val for x-axis expanding
#' @param legend.adj legend adjusted position
#' @param legend.text size of legend text
#' @param legend.title size of legend title
#' @param axis.text size of axis text
#' @param axis.title size of axis title
#' @param fontsize size of text of table
#' @param round.hr statistics round number
#' @param round.p statistics round number
#' @param n.round statistics round number
#' @param margin.plot margins of plot
#' @param margin.tb margins of table
#' @param heights.plot height of plot
#' @param xminlab length of risk table legend
#' @param xmaxlab space between risk table legend and table
#' @param xlim.lo x-axis limit
#' @param n.space the lists for adjusting spaces between values of legends
#'
#' @return KM curves
#'
#' @export
#'
#' @examples
#' # start example
#' attach(datplot)
#'
#' library(survival)
#' library(survminer)
#' library(dplyr)
#' library(ggplot2)
#' library(gridExtra)
#' library(grid)
#'
#' df.plot <- mutate(datplot,
#'                  grp=factor(grp, levels=c("Arm1_Detected", "Arm2_Detected",
#'                                           "Arm1_No-Detected", "Arm2_No-Detected"))) %>%
#'            select(-trt) %>%
#'            mutate(trt=trtn)
#'
#' p1 <- KM.w.table(dat=df.plot,
#'                 trt.var = "grp",
#'                 tte.var = "timevar", cen.var = "cenvar",
#'                 ttime.inc = 3,
#'                 cov.add=FALSE, cov.var=NULL, STRAT=FALSE,
#'                 pval.type="Log-rank",
#'                 COMP=FALSE,
#'                 COMPGRP=TRUE,
#'                 INTER=FALSE,
#'                 pairedgrp=list(c("Arm1_Detected", "Arm2_Detected"),
#'                                c("Arm1_No-Detected", "Arm2_No-Detected")),
#'                 pairedgrplv=list(c("Arm1_Detected", "Arm2_Detected"),
#'                                  c("Arm1_No-Detected", "Arm2_No-Detected")),
#'                 method="exact",
#'                 cols=c("red", "blue", "green", "orange"),
#'                 ltypes=c(1,1,1,1),
#'                 main="KM plots", xname="",
#'                 yname='Progression Free Survival',
#'                 level.grp=NULL,
#'                 legend.head="Patient Group                       Median    HR (95% CI)     P-value",
#'                 legend.type = 4,
#'                 legend.pos = c(0, 0),
#'                 legend.adj = c(-0.05, -0.1),
#'                 legend.text=10, legend.title=12,
#'                 axis.text=12, axis.title=13,
#'                 fontsize=4.5,
#'                 round.hr=2, round.p=3, n.round=2,
#'                 heights.plot = c(2, 0.85),
#'                 margin.plot=c(1, 1, 0, 1),
#'                 margin.tb=c(0, 0.45, 2, 0.35),
#'                 xminlab = -2.5, xmaxlab = -1, xlim.lo = -1,
#'                 n.space=list(c(16, 12, 21, 16, 6),
#'                              c(10, 12, 20, 10, 6)))
#'
#' # for interaction model
#' df.plot.int <- dplyr::mutate(df.plot, trt=1-trtn)
#'
#' p2 <- KM.w.table(dat=df.plot.int,
#'                 trt.var = "grp",
#'                 tte.var = "timevar", cen.var = "cenvar",
#'                 ttime.inc = 3,
#'                 cov.add=FALSE, cov.var=NULL, STRAT=FALSE,
#'                 pval.type="Log-rank",
#'                 COMP=FALSE,
#'                 COMPGRP=TRUE,
#'                 INTER=TRUE,
#'                 pairedgrp=list(c("Arm1_Detected", "Arm2_Detected"),
#'                                c("Arm1_No-Detected", "Arm2_No-Detected")),
#'                 pairedgrplv=list(c("Arm1_Detected", "Arm2_Detected"),
#'                                  c("Arm1_No-Detected", "Arm2_No-Detected")),
#'                 method="exact",
#'                 cols=c("red", "blue", "green", "orange"),
#'                 ltypes=c(1,1,1,1),
#'                 main="KM plots with interaction model", xname="",
#'                 yname='Progression Free Survival',
#'                 level.grp=NULL,
#'                 legend.head="Patient Group                       Median    HR (95% CI)     P-value",
#'                 legend.type = 4,
#'                 legend.pos = c(0, 0),
#'                 legend.adj = c(-0.05, -0.1),
#'                 legend.text=10, legend.title=12,
#'                 axis.text=12, axis.title=13,
#'                 fontsize=4.5,
#'                 round.hr=2, round.p=3, n.round=2,
#'                 heights.plot = c(2, 0.85),
#'                 margin.plot=c(1, 1, 0, 1),
#'                 margin.tb=c(0, 0.45, 2, 0.35),
#'                 xminlab = -2.5, xmaxlab = -1, xlim.lo = -1,
#'                 n.space=list(c(16, 12, 21, 16, 6),
#'                              c(10, 12, 20, 10, 6)))
#'
#'
#' pdf(file=file.path(getwd(), "f_km_plot_pair.comp.pdf"), width=11, height=8.5)
#' grid.newpage()
#' print(p1, vp = viewport(x = unit(0.5, "npc"), y = unit(0.5, "npc"), width=0.95, height=0.9, just = "centre"))
#' grid.newpage()
#' print(p2, vp = viewport(x = unit(0.5, "npc"), y = unit(0.5, "npc"), width=0.95, height=0.9, just = "centre"))
#' dev.off()
#'
# end example
KM.w.table <- function(dat, subcohana=FALSE,
                       trt.var, tte.var, cen.var,
                       cov.add=FALSE, cov.var=NULL, STRAT=FALSE,
                       pval.type="Wald",
                       COMP=FALSE,         # there are only two groups
                       COMPGRP=TRUE,       # paired comparison
                       INTER=FALSE,        # for interaction
                       pairedgrp=list(),   # giving paired subgrp for comparison
                       pairedgrplv=list(), # giving paired subgrp levels
                       conf.level=0.95, conf.type='log-log', method='efron',
                       ttime.inc=2,
                       cols=c("red", "blue", "orange", "green"),
                       ltypes=c(1,1,1,1),
                       main="KM plots", xname="",
                       yname='Progression Free Survival',
                       level.grp=NULL,
                       legend.head="",     # for legend title if apply
                       legend.type = 4,
                       legend.pos = c(1,1),
                       expand.val = c(0, 0),
                       legend.adj = c(1.02,1.015), # for top/righ c(0.3, 0.15) is for bottom/left
                       legend.text=10, legend.title=12,
                       axis.text=12, axis.title=13,
                       fontsize=4.5,
                       round.hr=2, round.p=3, n.round=3,
                       margin.plot=c(20, 20, 30, 85),
                       margin.tb=c(0, 20, 30, 2),
                       heights.plot = c(2, 0.4),
                       xminlab=-2.5, xmaxlab = -1, xlim.lo = 0,
                       n.space=list(c(5, 5, 5, 5, 5),
                                    c(5, 5, 5, 5, 5),
                                    c(5, 5, 5, 5, 5))) {

  # rename variables
  names(dat)[match(trt.var, names(dat))] <- "grp"
  names(dat)[match(tte.var, names(dat))] <- "timevar"
  names(dat)[match(cen.var, names(dat))] <- "cenvar"
  dat$evnt <- 1- as.numeric(dat$cenvar)

  # called misc function for computing stat
  ret <- Km.w.table.misc(dat=dat, subcohana=subcohana,
                         conf.level=conf.level, ttime.inc=ttime.inc,
                         conf.type=conf.type, method=method,
                         level.grp=level.grp, STRAT=STRAT,
                         cov.add = cov.add, cov.var = cov.var,
                         COMP=COMP, COMPGRP=COMPGRP, INTER=INTER,
                         pval.type=pval.type, legend.type=legend.type,
                         pairedgrp=pairedgrp, pairedgrplv=pairedgrplv,
                         round.hr=round.hr, round.p=round.p, n.round=n.round,
                         n.space = n.space)


  #--- plot KM
  fit <- ret$km

  # plot
  legend.lb <- ret$legend.lb

  p <- ggsurvplot(fit, data = dat,
                  size = 1, censor.shape = "|", censor.size = 4,
                  palette = cols, conf.int = F,
                  linetype = ltypes,
                  risk.table = TRUE, fontsize = fontsize, tables.theme = clean_theme(),
                  xlab = xname, ylab = yname, title=main,
                  surv.scale = "percent",
                  break.time.by = ttime.inc, break.y.by = 0.1,
                  xlim = c(xlim.lo, (max(dat$timevar)+ttime.inc)),
                  surv.plot.height = 0.5,
                  risk.table.height = 0.6,
                  risk.table.y.text =FALSE,
                  legend.title = legend.head,
                  legend.labs = legend.lb)

  p$plot <- p$plot +
    scale_x_continuous(breaks = seq(0, (max(dat$timevar)+ttime.inc), ttime.inc), expand = expand.val) +
    scale_y_continuous(breaks=seq(0, 10, 0.1), limits = c(0, 1),
                       labels = scales::percent_format(accuracy=5L)) +
    geom_hline(yintercept = .5, color = 'grey28', linetype = 'dashed') +
    theme(legend.position = legend.pos,
          legend.justification = legend.adj,
          legend.key.width=unit(1.5, "cm"),
          legend.title = element_text(size = legend.title, face = "bold"),
          legend.text = element_text(size = legend.text, face = "bold"),
          axis.title.y=element_text(size = axis.title, face = "bold"),
          axis.title.x=element_text(size = axis.title, face = "bold"),
          axis.text.x =element_text(size = axis.text, face = "bold"),
          axis.text.y =element_text(size = axis.text, vjust=0.5, hjust=0.5, face = "bold"),
          axis.ticks=element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          text=element_text(size = axis.title, face = "bold"),
          plot.margin = unit(margin.plot, "cm"))


  # update risk table
  p$table <- p$table +
    ggtitle("Patients at risk") +
    theme(plot.title = element_text(hjust = 0.05, size = 11),
          plot.margin = unit(margin.tb, "cm"),
          axis.text.y=element_blank(),
          axis.title.y=element_blank())

  for(i in 1:length(legend.lb)) {
    p$table <- p$table +
      annotate("segment",
               x = xminlab, xend = xmaxlab, y = i, yend = i,
               colour = rev(cols)[i], linewidth = 1) +
      coord_cartesian(clip = "off")
  }

  p.out <- ggarrange(p$plot, p$table, heights = heights.plot, ncol = 1, nrow = 2)

  return(p.out)

}#end

###----------------------------------------------------------------
#
###--- misc functions
#
###----------------------------------------------------------------

#' @export
format_num <- function(x, num.sigf=2, width=1, num.space=1, trun=FALSE) {
  num <- formatC(format(round(x, num.sigf), nsmall=num.sigf),
                 width=width, flag=paste(rep("", num.space), collapse=" "))
  if(trun==TRUE) {
    num <- gsub(" ", "", num)
  }
  return(num)
}#end format_num

#--- function to replace value
#' @export
replace_val <- function(x, replaceVal=NA, byVal=NULL) {
  id <- which(x %in% replaceVal)
  x[id] <- byVal
  x
}#end

# function to map Lilly color with chart R color
#' @export
color.code <- function(str.color) {
  colz <- NULL
  for(i in 1:length(str.color)) {
    colz[i] <- str.color[i]
    if(str.color[i] %in% "red")     {colz[i] <- '#D52B1E'}
    if(str.color[i] %in% "blue")    {colz[i] <- '#00A1DE'}
    if(str.color[i] %in% "green")   {colz[i] <- '#00AF3F'}
    if(str.color[i] %in% "lightgreen")  {colz[i] <- '#C2EDCE'}
    if(str.color[i] %in% "darkgreen")   {colz[i] <- '#275E37'}
    if(str.color[i] %in% "yellow")  {colz[i] <- '#FED100'}
    if(str.color[i] %in% "orange")  {colz[i] <- '#FF6D22'}
    if(str.color[i] %in% c('grey', 'gray'))           {colz[i] <- '#A59D95'}
    if(str.color[i] %in% c('lightgrey', 'lightgray')) {colz[i] <- '#D5D2CA'}
    if(str.color[i] %in% "brown")    {colz[i] <- '#4E2E2D'}

  }#end
  return(colz)
}#end

#' @export
tte.get.result <- function(fit, confint.level,
                           pval.type="Wald",
                           n.round=3, round.p=3) {
  est.val <- exp(summary(fit)$coef[1])
  CI <- summary(fit, conf.int=confint.level)$conf.int
  ci.low <- CI[,3]
  ci.up <- CI[,4]
  if(pval.type == "Wald") {
    pval <- summary(fit)$waldtest[3]
  } else if(pval.type == "Ratio") {
    pval <- summary(fit)$logtest[3]
  } else if(pval.type == "Log-rank") {
    pval <- summary(fit)$sctest[3]
  } else {
    pval <- NULL
  }

  result <- data.frame(HR=format_num(est.val, n.round),
                       ci.lo=format_num(ci.low, n.round),
                       ci.up=format_num(ci.up, n.round),
                       pval=ifelse(pval < 0.001, "<.001", format_num(pval, round.p)))
  return(result)
}#end

#' @export
surv.med.comp <- function(dat, round.n=2) {
  # compute n and min, max by grp
  tb <- dplyr::group_by(dat, grp) %>%
    summarise(n=n(),
              min=min(timevar),
              max=max(timevar)) %>%
    ungroup() %>% data.frame()

  # make sure the levels would be consistant
  fit <- survfit(Surv(timevar, cenvar) ~ grp, data = dat, conf.type="log-log")
  x <- data.frame(summary(fit)$table)
  rownames(x) <- gsub("grp=\\b", "", rownames(x))

  val.med <- quantile(fit, probs = c(0.25, 0.5, 0.75))$quantile
  rownames(val.med) <- gsub("grp=\\b", "", rownames(val.med))
  val.lo <- quantile(fit, probs = c(0.25, 0.5, 0.75))$lower
  rownames(val.lo) <- gsub("grp=\\b", "", rownames(val.lo))
  val.up <- quantile(fit, probs = c(0.25, 0.5, 0.75))$upper
  rownames(val.up) <- gsub("grp=\\b", "", rownames(val.up))

  x2 <- t(x) %>% data.frame() %>%
    mutate(subgrp=colnames(x))
  val.med2 <- t(val.med) %>% data.frame() %>%
    mutate(subgrp=colnames(val.med))
  val.lo2 <- t(val.lo) %>% data.frame() %>%
    mutate(subgrp=colnames(val.lo))
  val.up2 <- t(val.up) %>% data.frame() %>%
    mutate(subgrp=colnames(val.up))

  ret <- data.frame(subgrp=c("Number (%) of Subjects with Events",
                             "Number (%) of Subjects Censored",
                             "Observed range (min, max)",
                             "25th percentile (95% CI)",
                             "Median (95% CI)",
                             "75th percentile (95% CI)"),
                    grp1=c(paste0(x2[x2$subgrp %in% "events", 1], " ( ",
                                  round(x2[x2$subgrp %in% "events", 1]/tb$n[1]*100, round.n), ")"),
                           paste0((tb$n[1] - x2[x2$subgrp %in% "events", 1]), " ( ",
                                  round((tb$n[1] - x2[x2$subgrp %in% "events", 1])/tb$n[1]*100, round.n), ")"),
                           paste0("( ", round(tb$min[1], round.n), " - ", round(tb$max[1], round.n), ")"),
                           paste0(round(val.med2[1,1], round.n), " ( ",
                                  round(val.lo2[1,1], round.n), " - ", round(val.up2[1,1], round.n), ")"),
                           paste0(round(val.med2[2,1], round.n), " ( ",
                                  round(val.lo2[2,1], round.n), " - ", round(val.up2[2,1], round.n), ")"),
                           paste0(round(val.med2[3,1], round.n), " ( ",
                                  round(val.lo2[3,1], round.n), " - ", round(val.up2[3,1], round.n), ")")

                    ),
                    grp2=c(paste0(x2[x2$subgrp %in% "events", 2], " ( ",
                                  round(x2[x2$subgrp %in% "events", 2]/tb$n[2]*100, round.n), ")"),
                           paste0((tb$n[2] - x2[x2$subgrp %in% "events", 2]), " ( ",
                                  round((tb$n[2] - x2[x2$subgrp %in% "events", 2])/tb$n[2]*100, round.n), ")"),
                           paste0("( ", round(tb$min[2], round.n), " - ", round(tb$max[2], round.n), ")"),
                           paste0(round(val.med2[1,2], round.n), " ( ",
                                  round(val.lo2[1,2], round.n), " - ", round(val.up2[1,2], round.n), ")"),
                           paste0(round(val.med2[2,2], round.n), " ( ",
                                  round(val.lo2[2,2], round.n), " - ", round(val.up2[2,2], round.n), ")"),
                           paste0(round(val.med2[3,2], round.n), " ( ",
                                  round(val.lo2[3,2], round.n), " - ", round(val.up2[3,2], round.n), ")")
                    ),

                    grp3=c(paste0(x2[x2$subgrp %in% "events", 3], " ( ",
                                  round(x2[x2$subgrp %in% "events", 3]/tb$n[3]*100, round.n), ")"),
                           paste0((tb$n[3] - x2[x2$subgrp %in% "events", 3]), " ( ",
                                  round((tb$n[3] - x2[x2$subgrp %in% "events", 3])/tb$n[3]*100, round.n), ")"),
                           paste0("( ", round(tb$min[3], round.n), " - ", round(tb$max[3], round.n), ")"),
                           paste0(round(val.med2[1,3], round.n), " ( ",
                                  round(val.lo2[1,3], round.n), " - ", round(val.up2[1,3], round.n), ")"),
                           paste0(round(val.med2[2,3], round.n), " ( ",
                                  round(val.lo2[2,3], round.n), " - ", round(val.up2[2,3], round.n), ")"),
                           paste0(round(val.med2[3,3], round.n), " ( ",
                                  round(val.lo2[3,3], round.n), " - ", round(val.up2[3,3], round.n), ")")
                    ),

                    grp4=c(paste0(x2[x2$subgrp %in% "events", 4], " ( ",
                                  round(x2[x2$subgrp %in% "events", 4]/tb$n[4]*100, round.n), ")"),
                           paste0((tb$n[4] - x2[x2$subgrp %in% "events", 4]), " ( ",
                                  round((tb$n[4] - x2[x2$subgrp %in% "events", 4])/tb$n[4]*100, round.n), ")"),
                           paste0("( ", round(tb$min[4], round.n), " - ", round(tb$max[4], round.n), ")"),
                           paste0(round(val.med2[1,4], round.n), " ( ",
                                  round(val.lo2[1,4], round.n), " - ", round(val.up2[1,4], round.n), ")"),
                           paste0(round(val.med2[2,4], round.n), " ( ",
                                  round(val.lo2[2,4], round.n), " - ", round(val.up2[2,4], round.n), ")"),
                           paste0(round(val.med2[3,4], round.n), " ( ",
                                  round(val.lo2[3,4], round.n), " - ", round(val.up2[3,4], round.n), ")"))
  )

  return(ret)
}#end

#' @export
surv.summ <- function(datIn,
                      grpval=c("Exon 19 Shedding Ram", "Exon 19 Shedding Placebo"),
                      grplv=c("Exon 19 Shedding Placebo", "Exon 19 Shedding Ram"),
                      conf.level=0.95, cov.add=FALSE,
                      cov.var="",
                      conf.type='log-log', method='efron',
                      round.p=3) {
  # subset subgrp data
  dat <- filter(datIn, grp %in% grpval) %>%
    mutate(grp=factor(as.character(grp), levels = grplv))


  # get log-rank p-val --------------------------------------------------------
  # for both cases: multivariates and stratificaiton of factor covariates
  if(cov.add) {
    form.cov <- paste0("Surv(timevar, cenvar) ~ grp + ", paste(cov.var, collapse = " + "))
    form.strat <- paste0("Surv(timevar, cenvar) ~ grp + ", paste(paste0("strata(", cov.var, ")"), collapse = " + "))
  } else {
    form.cov <- "Surv(timevar, cenvar) ~ grp"
    form.strat <- "Surv(timevar, cenvar) ~ grp"
  }
  fit.diff.cov <- survdiff(as.formula(eval(expression(form.cov))), data=dat)
  fit.diff.strat <- survdiff(as.formula(eval(expression(form.strat))), data=dat)
  p.logr.cov <- 1 - pchisq(fit.diff.cov$chisq, (length(fit.diff.cov$n) - 1))
  p.logr.strat <- 1 - pchisq(fit.diff.strat$chisq, (length(fit.diff.strat$n) - 1))

  if(p.logr.cov < 0.0001) {
    p.val.cov <- "< 0.0001"
  } else {
    p.val.cov <- format(p.logr.cov, nsmall=round.p, digits=1, scientific=FALSE)
  }#end
  if(p.logr.strat < 0.0001) {
    p.val.strat <- "< 0.0001"
  } else {
    p.val.strat <- format(p.logr.strat, nsmall=round.p, digits=1, scientific=FALSE)
  }#end


  # compute HR (CIs), pvalue if comparison of two subgroups ---------------------------------------
  if(cov.add) {
    form.full <- paste0("Surv(timevar, cenvar) ~ grp + ", paste(paste0("strata(", cov.var, ")"), collapse = " + "))
    form.rd <- paste0("Surv(timevar, cenvar) ~ 1 + ", paste(paste0("strata(", cov.var, ")"), collapse = " + "))

    cox.full <- coxph(as.formula(eval(expression(form.full))), data=dat, method=method)
    cox.rd <- coxph(as.formula(eval(expression(form.rd))), data = dat, method=method)
    mul <- summary(cox.full)[["logtest"]]["test"]
    pval <- ifelse(pchisq(mul, df=1, lower.tail=FALSE) < 0.0001, "< 0.0001",
                   format_num(pchisq(mul, df=1, lower.tail=FALSE), round.p))

  } else { # none covariates
    form.full <- "Surv(timevar, cenvar) ~ grp"
    form.rd <- "Surv(timevar, cenvar) ~ 1"
    cox.full <- coxph(as.formula(eval(expression(form.full))), data=dat, method=method)
    mul <- summary(cox.full)[["logtest"]]["test"]
    pval <- ifelse(pchisq(mul, df=1, lower.tail=FALSE) < 0.0001, "< 0.0001",
                   format_num(pchisq(mul, df=1, lower.tail=FALSE), round.p))

  }#end
  val <- tte.get.result(fit=cox.full, confint.level=conf.level,
                        pval.type="Wald", n.round=2, round.p=3)


  ret <- data.frame(subgrp=c("Treatment Effect", "Log-rank p-value", "HR (95% CI)"),
                    hr=c("", "", paste0(val$HR, " ( ", val$ci.lo, " - ", val$ci.up, ")")),
                    pval=c("", p.val.strat, val$pval) )


  # for return results ------------------------------------------------------
  return(ret)
}#end

#' @export
tte.legend.form <- function(x, subgrpComp,
                            hr, ci, pval,
                            legend.type=1,
                            n.space=list(c(5, 5, 5, 5, 5),
                                         c(5, 5, 5, 5, 5),
                                         c(5, 5, 5, 5, 5))) {
  rownames(x) <- gsub("GRP=\\b", "", rownames(x))

  if (legend.type==1) { # for the case legend including median, HR
    legend.txt <- NULL
    for(i in 1:length(subgrpComp)) {
      grpname <- subgrpComp[[i]]
      n.sp <- n.space[[i]]
      val <- c(paste0(grpname[1], paste(rep(" ", n.sp[1]), collapse = ""),
                      format_num(x$median[rownames(x)%in%grpname[1]], 2),
                      paste(rep(" ", n.sp[2]), collapse = ""), hr[i]),
               paste0(grpname[2], paste(rep(" ", n.sp[3]), collapse = ""),
                      format_num(x$median[rownames(x)%in%grpname[2]], 2),
                      paste(rep(" ", n.sp[4]), collapse = ""), ci[i]))
      if(is.null(legend.txt)) legend.txt <- val else legend.txt <- c(legend.txt, val)
    }
  } else if(legend.type==2) { # for the case legend including HR, log-rank pvalue
    legend.txt <- NULL
    for(i in 1:length(subgrpComp)) {
      grpname <- subgrpComp[[i]]
      n.sp <- n.space[[i]]
      val <- c(paste0(grpname[1], paste(rep(" ", n.sp[1]), collapse = ""), hr[i],
                      paste(rep(" ", n.sp[2]), collapse = ""), pval[i]),
               paste0(grpname[2], paste(rep(" ", n.sp[3]), collapse = ""), ci[i]))
      if(is.null(legend.txt)) legend.txt <- val else legend.txt <- c(legend.txt, val)
    }
  } else if(legend.type==10) { # for special case, where diplaying N, Event, Median (CIs)
    legend.txt <- NULL
    for(i in 1:length(subgrpComp)) {
      grpname <- subgrpComp[[i]]
      n.sp <- n.space[[i]]
      val <- c(paste0(grpname[1], paste(rep(" ", n.sp[1]), collapse = ""),
                      x$records[rownames(x)%in%grpname[1]], paste(rep(" ", n.sp[2]), collapse = ""),
                      x$events[rownames(x)%in%grpname[1]], paste(rep(" ", n.sp[3]), collapse = ""),
                      format_num(x$median[rownames(x)%in%grpname[1]], 2), " (",
                      format_num(x[rownames(x)%in%grpname[1], 8], 2), ", ",
                      format_num(x[rownames(x)%in%grpname[1], 9], 2), ")"),
               paste0(grpname[2], paste(rep(" ", n.sp[4]), collapse = ""),
                      x$records[rownames(x)%in%grpname[2]], paste(rep(" ", n.sp[5]), collapse = ""),
                      x$events[rownames(x)%in%grpname[2]], paste(rep(" ", n.sp[6]), collapse = ""),
                      format_num(x$median[rownames(x)%in%grpname[2]], 2), " (",
                      format_num(x[rownames(x)%in%grpname[2], 8], 2), ", ",
                      format_num(x[rownames(x)%in%grpname[2], 9], 2), ")"))
      if(is.null(legend.txt)) legend.txt <- val else legend.txt <- c(legend.txt, val)
    }

  } else { # for the case legend including median, HR, log-rank pvalue
    legend.txt <- NULL
    for(i in 1:length(subgrpComp)) {
      grpname <- subgrpComp[[i]]
      n.sp <- n.space[[i]]
      val <- c(paste0(grpname[1], paste(rep(" ", n.sp[1]), collapse = ""),
                      format_num(x$median[rownames(x)%in%grpname[1]], 2),
                      paste(rep(" ", n.sp[2]), collapse = ""), hr[i],
                      paste(rep(" ", n.sp[3]), collapse = ""), pval[i]),
               paste0(grpname[2], paste(rep(" ", n.sp[4]), collapse = ""),
                      format_num(x$median[rownames(x)%in%grpname[2]], 2),
                      paste(rep(" ", n.sp[5]), collapse = ""), ci[i]))
      if(is.null(legend.txt)) legend.txt <- val else legend.txt <- c(legend.txt, val)
    }
  }

  return(legend.txt)

}#end

# function to compute log-rank pvalue, HR (CI)
# n, median, and form legends
# also computing number at risk table
#' @export
Km.w.table.misc <- function(dat, subcohana=FALSE,
                            conf.level=0.95, ttime.inc=2,
                            conf.type='log-log', method='efron',
                            level.grp=NULL,
                            COMP=FALSE,         # there are only two groups
                            COMPGRP=TRUE,       # paired comparison
                            INTER=FALSE,        # for interaction model
                            pairedgrp=list(),   # giving paired subgrp for comparison
                            pairedgrplv=list(), # giving paired subgrp levels
                            cov.var="", cov.add=TRUE, STRAT=FALSE,
                            pval.type="Wald",
                            round.hr=2, round.p=3, n.round=3,
                            legend.type=4,
                            n.space=list(c(5, 5, 5, 5, 5),
                                         c(5, 5, 5, 5, 5),
                                         c(5, 5, 5, 5, 5))) {

  # compute KM
  km <- survfit(Surv(timevar,  evnt) ~ grp, data=dat, conf.type=conf.type, conf.int=conf.level)
  x <- data.frame(summary(km)$table)
  rownames(x) <- gsub("grp=\\b", "", rownames(x))

  ciname <- paste0(100*conf.level,"% CI")

  # get log-rank p-val --------------------------------------------------------
  # for both cases: multivariates and stratificaiton of factor covariates
  if(cov.add) {
    form.cov <- paste0("Surv(timevar, evnt) ~ grp + ", paste(cov.var, collapse = " + "))
    form.strat <- paste0("Surv(timevar, evnt) ~ grp + ", paste(paste0("strata(", cov.var, ")"), collapse = " + "))
  } else {
    form.cov <- "Surv(timevar, evnt) ~ grp"
    form.strat <- "Surv(timevar, evnt) ~ grp"
  }
  fit.diff.cov <- survdiff(as.formula(eval(expression(form.cov))), data=dat)
  fit.diff.strat <- survdiff(as.formula(eval(expression(form.strat))), data=dat)
  p.logr.cov <- 1 - pchisq(fit.diff.cov$chisq, (length(fit.diff.cov$n) - 1))
  p.logr.strat <- 1 - pchisq(fit.diff.strat$chisq, (length(fit.diff.strat$n) - 1))

  if(p.logr.cov < 0.0001) {
    p.val.cov <- "< 0.0001"
  } else {
    p.val.cov <- format(p.logr.cov, nsmall=round.p, digits=1, scientific=FALSE)
  }#end
  if(p.logr.strat < 0.0001) {
    p.val.strat <- "< 0.0001"
  } else {
    p.val.strat <- format(p.logr.strat, nsmall=round.p, digits=1, scientific=FALSE)
  }#end


  # compute HR (CIs), pvalue if comparison of two subgroups ---------------------------------------
  if(cov.add) {
    if(STRAT) {
      form.full <- paste0("Surv(timevar, evnt) ~ grp + ", paste(paste0("strata(", cov.var, ")"), collapse = " + "))
      form.rd <- paste0("Surv(timevar, evnt) ~ 1 + ", paste(paste0("strata(", cov.var, ")"), collapse = " + "))

      cox.full <- coxph(as.formula(eval(expression(form.full))), data=dat, method=method)
      cox.rd <- coxph(as.formula(eval(expression(form.rd))), data = dat, method=method)
      mul <- summary(cox.full)[["logtest"]]["test"]
      pval <- ifelse(pchisq(mul, df=1, lower.tail=FALSE) < 0.0001, "< 0.0001",
                     format_num(pchisq(mul, df=1, lower.tail=FALSE), round.p))

    } else {
      form.full <- paste0("Surv(timevar, evnt) ~ grp + ", paste(cov.var, collapse = " + "))
      form.rd <- paste0("Surv(timevar, evnt) ~ 1 + ", paste(cov.var, collapse = " + "))
      cox.full <- coxph(as.formula(eval(expression(form.full))), data=dat, method=method)
      cox.rd <- coxph(as.formula(eval(expression(form.rd))), data = dat, method=method)
      mul <- summary(cox.full)[["logtest"]]["test"] - summary(cox.rd)[["logtest"]]["test"]
      pval <- ifelse(pchisq(mul, df=1, lower.tail=FALSE) < 0.0001, "< 0.0001",
                     format_num(pchisq(mul, df=1, lower.tail=FALSE), round.p))
    }#end

  } else { # none covariates
    form.full <- "Surv(timevar, evnt) ~ grp"
    form.rd <- "Surv(timevar, evnt) ~ 1"
    cox.full <- coxph(as.formula(eval(expression(form.full))), data=dat, method=method)
    mul <- summary(cox.full)[["logtest"]]["test"]
    pval <- ifelse(pchisq(mul, df=1, lower.tail=FALSE) < 0.0001, "< 0.0001",
                   format_num(pchisq(mul, df=1, lower.tail=FALSE), round.p))
  }#end

  # for HR (CIs) overall
  if (COMP) {  # there are only two subgroups
    hr <- summary(cox.full)[["coefficients"]][,"exp(coef)"]
    ci.lo <- summary(cox.full)[["conf.int"]][,3]
    ci.hi <- summary(cox.full)[["conf.int"]][,4]
    p.val <- summary(cox.full)[["coefficients"]][,"Pr(>|z|)"]
    HR.CI <- paste0(format_num(hr, round.hr), " (",
                    format_num(ci.lo, round.hr), "-", format_num(ci.hi, round.hr), ")",
                    ", [", format_num(p.val, round.p), "]")
  } else {  # more than two subgroups, not good for comparison
    hr <- NA; ci.lo <- NA; ci.hi <- NA
    HR.CI <- ""
  }

  if (COMPGRP) {  # paired two subgroups
    # pairedgrp=list(c("Arm1_Detected", "Arm1_No-Detected"),
    #                c("Arm2_Detected", "Arm2_No-Detected"))
    # pairedgrplv=list(c("Arm1_Detected", "Arm1_No-Detected"),
    #                  c("Arm2_Detected", "Arm2_No-Detected"))
    if(subcohana) {
      # using cch function in the survival package
      hr.result <- NULL
      for(i in 1:length(pairedgrp)) {
        dat.comp <- dplyr::filter(dat, grp %in% pairedgrp[[i]]) %>%
          dplyr::mutate(grp=factor(grp, levels = pairedgrplv[[i]]))
        if(any(grepl("TR_", pairedgrplv[[i]]))) {
          dat.comp$subcohort <- as.logical(dat.comp$subcohort)
          fit <- cch(Surv(timevar, evnt) ~ grp, data=dat.comp,
                     id=~SUBJID, subcoh=~subcohort, cohort.size=unique(dat.comp$N))

          vari <- fit$var
          z.norm <- qnorm(1-(1-conf.level)/2, lower.tail=T)
          b <- sqrt(vari)
          HR <- exp(summary(fit)$coef[1])
          ci.lo <- HR*exp(-z.norm*b)
          ci.up <- HR*exp(z.norm*b)

          hr.out <- data.frame(HR=format_num(HR, round.hr),
                               ci.lo=format_num(ci.lo, round.hr), ci.up=format_num(ci.up, round.hr),
                               pval=format_num(summary(fit)$coeff[[4]], round.p))

        } else {
          fit <- coxph(as.formula(eval(expression(form.full))), data=dat.comp, method=method)
          hr.out <- tte.get.result(fit=fit, confint.level=conf.level,
                                   n.round=n.round, round.p=round.p, pval.type=pval.type)
        }

        if(is.null(hr.result)) hr.result <- hr.out else hr.result <- rbind(hr.result, hr.out)
      }#end

    } else {
      hr.result <- NULL
      for(i in 1:length(pairedgrp)) {
        dat.comp <- dplyr::filter(dat, grp %in% pairedgrp[[i]]) %>%
          dplyr::mutate(grp=factor(grp, levels = pairedgrplv[[i]]))
        fit <- coxph(as.formula(eval(expression(form.full))), data=dat.comp, method=method)
        hr.out <- tte.get.result(fit=fit, confint.level=conf.level,
                                 n.round=n.round, round.p=round.p, pval.type=pval.type)
        if(is.null(hr.result)) hr.result <- hr.out else hr.result <- rbind(hr.result, hr.out)
      }#end
    }#end

    rownames(hr.result) <- 1:nrow(hr.result)

    # for interaction model: need to convert trt and mrk as 0/1 numeric
    # before running the function
    if (INTER) {
      df.int <- dplyr::mutate(dat,
                              #trt=ifelse(trt %in% "RAM", 1, 0),
                              int=mrk*trt)
      fit.full <- coxph(Surv(timevar, evnt) ~ mrk + trt + int, data=df.int, method=method)
      fit.rd <- coxph(Surv(timevar, evnt) ~ mrk + trt, data=df.int, method=method)

      # get HR (CIs)
      vari <- fit.full$var  # var matrix
      z.norm <- qnorm(1-(1-conf.level)/2,lower.tail=T)
      CI <- summary(fit.full, conf.int=conf.level)$conf.int

      # For High level
      HR.hi <- exp(sum(fit.full$coef[c("trt", "int")]))  # for trt and int
      b <- ifelse(HR.hi %in% c(NA, "Inf"), NA, sqrt(vari[2, 2]+vari[3, 3]+2*vari[2, 3]))
      HR.hi.l <- HR.hi*exp(-z.norm*b)
      HR.hi.u <- HR.hi*exp(z.norm*b)

      # For Low level
      HR.lo <- exp(sum(fit.full$coef["trt"]))
      HR.lo.l <- CI["trt",3]
      HR.lo.u <- CI["trt",4]

      hr.result <- data.frame(HR=c(format_num(HR.hi, round.hr) , format_num(HR.lo, round.hr)),
                              ci.lo=c(format_num(HR.hi.l, n.round), format_num(HR.lo.l, n.round)),
                              ci.up=c(format_num(HR.hi.u, n.round), format_num(HR.lo.u, n.round)),
                              pval=hr.result$pval)

    }#end

  } else {
    hr.result <- NULL
  }

  # form legends ------------------------------------------------------------
  ret.lb <- tte.legend.form(x=x, subgrpComp=pairedgrp, legend.type=legend.type,
                            hr=hr.result$HR,
                            ci=paste0("(", hr.result$ci.lo, ", ", hr.result$ci.up, ")"),
                            pval=hr.result$pval,
                            n.space=n.space)

  # for return results ------------------------------------------------------
  return(list(df.sur=dat, df.km=km, legend.lb=ret.lb,
              p.diff.cov=p.val.cov, p.diff.strat=p.val.strat,
              pval=pval, hr=hr, hr.lo=ci.lo, hr.hi=ci.hi, HR.CI=HR.CI,
              x.tb=x, km=km, HR.paired=hr.result))

}#end

#-------------------------------------------------------------------

