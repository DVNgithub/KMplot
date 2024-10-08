% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/KMplot_func.R
\name{KM.w.table}
\alias{KM.w.table}
\title{KM plot with more than two subgroups}
\usage{
KM.w.table(
  dat,
  subcohana = FALSE,
  trt.var,
  tte.var,
  cen.var,
  cov.add = FALSE,
  cov.var = NULL,
  STRAT = FALSE,
  pval.type = "Wald",
  COMP = FALSE,
  COMPGRP = TRUE,
  INTER = FALSE,
  pairedgrp = list(),
  pairedgrplv = list(),
  conf.level = 0.95,
  conf.type = "log-log",
  method = "efron",
  ttime.inc = 2,
  cols = c("red", "blue", "orange", "green"),
  ltypes = c(1, 1, 1, 1),
  main = "KM plots",
  xname = "",
  yname = "Progression Free Survival",
  level.grp = NULL,
  legend.head = "",
  legend.type = 4,
  legend.pos = c(1, 1),
  expand.val = c(0, 0),
  legend.adj = c(1.02, 1.015),
  legend.text = 10,
  legend.title = 12,
  axis.text = 12,
  axis.title = 13,
  fontsize = 4.5,
  round.hr = 2,
  round.p = 3,
  n.round = 3,
  margin.plot = c(20, 20, 30, 85),
  margin.tb = c(0, 20, 30, 2),
  heights.plot = c(2, 0.4),
  xminlab = -2.5,
  xmaxlab = -1,
  xlim.lo = 0,
  n.space = list(c(5, 5, 5, 5, 5), c(5, 5, 5, 5, 5), c(5, 5, 5, 5, 5))
)
}
\arguments{
\item{dat}{data plot}

\item{subcohana}{subcohort analysis}

\item{trt.var}{treatments/arms}

\item{tte.var}{time to event variable}

\item{cen.var}{censor variable}

\item{cov.add}{whether using covariate}

\item{cov.var}{covarite varibale}

\item{STRAT}{whether using stratification}

\item{pval.type}{p-value type. Default is "Wald"}

\item{COMP}{there are only two groups comparison}

\item{COMPGRP}{paired comparison}

\item{INTER}{whether using interaction model}

\item{pairedgrp}{giving paired subgrp for comparison}

\item{pairedgrplv}{giving paired subgrp levels}

\item{conf.level}{confident level}

\item{conf.type}{confident type}

\item{method}{method using for model. Default is 'efron'}

\item{ttime.inc}{incremate interval time}

\item{cols}{KM curves colors}

\item{ltypes}{KM curves line types}

\item{main}{title of plot}

\item{xname}{x-axis lab}

\item{yname}{y-axis lab}

\item{level.grp}{subgroup levels}

\item{legend.head}{legend title}

\item{legend.type}{for different legends format to display}

\item{legend.pos}{legend position}

\item{expand.val}{for x-axis expanding}

\item{legend.adj}{legend adjusted position}

\item{legend.text}{size of legend text}

\item{legend.title}{size of legend title}

\item{axis.text}{size of axis text}

\item{axis.title}{size of axis title}

\item{fontsize}{size of text of table}

\item{round.hr}{statistics round number}

\item{round.p}{statistics round number}

\item{n.round}{statistics round number}

\item{margin.plot}{margins of plot}

\item{margin.tb}{margins of table}

\item{heights.plot}{height of plot}

\item{xminlab}{length of risk table legend}

\item{xmaxlab}{space between risk table legend and table}

\item{xlim.lo}{x-axis limit}

\item{n.space}{the lists for adjusting spaces between values of legends}
}
\value{
KM curves
}
\description{
KM plot with more than two subgroups
}
\examples{
# start example
attach(datplot)

library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(grid)

df.plot <- mutate(datplot,
                 grp=factor(grp, levels=c("Arm1_Detected", "Arm2_Detected",
                                          "Arm1_No-Detected", "Arm2_No-Detected"))) \%>\%
           select(-trt) \%>\%
           mutate(trt=trtn)

p1 <- KM.w.table(dat=df.plot,
                trt.var = "grp",
                tte.var = "timevar", cen.var = "cenvar",
                ttime.inc = 3,
                cov.add=FALSE, cov.var=NULL, STRAT=FALSE,
                pval.type="Log-rank",
                COMP=FALSE,
                COMPGRP=TRUE,
                INTER=FALSE,
                pairedgrp=list(c("Arm1_Detected", "Arm2_Detected"),
                               c("Arm1_No-Detected", "Arm2_No-Detected")),
                pairedgrplv=list(c("Arm1_Detected", "Arm2_Detected"),
                                 c("Arm1_No-Detected", "Arm2_No-Detected")),
                method="exact",
                cols=c("red", "blue", "green", "orange"),
                ltypes=c(1,1,1,1),
                main="KM plots", xname="",
                yname='Progression Free Survival',
                level.grp=NULL,
                legend.head="Patient Group                       Median    HR (95\% CI)     P-value",
                legend.type = 4,
                legend.pos = c(0, 0),
                legend.adj = c(-0.05, -0.1),
                legend.text=10, legend.title=12,
                axis.text=12, axis.title=13,
                fontsize=4.5,
                round.hr=2, round.p=3, n.round=2,
                heights.plot = c(2, 0.85),
                margin.plot=c(1, 1, 0, 1),
                margin.tb=c(0, 0.45, 2, 0.35),
                xminlab = -2.5, xmaxlab = -1, xlim.lo = -1,
                n.space=list(c(16, 12, 21, 16, 6),
                             c(10, 12, 20, 10, 6)))

# for interaction model
df.plot.int <- dplyr::mutate(df.plot, trt=1-trtn)

p2 <- KM.w.table(dat=df.plot.int,
                trt.var = "grp",
                tte.var = "timevar", cen.var = "cenvar",
                ttime.inc = 3,
                cov.add=FALSE, cov.var=NULL, STRAT=FALSE,
                pval.type="Log-rank",
                COMP=FALSE,
                COMPGRP=TRUE,
                INTER=TRUE,
                pairedgrp=list(c("Arm1_Detected", "Arm2_Detected"),
                               c("Arm1_No-Detected", "Arm2_No-Detected")),
                pairedgrplv=list(c("Arm1_Detected", "Arm2_Detected"),
                                 c("Arm1_No-Detected", "Arm2_No-Detected")),
                method="exact",
                cols=c("red", "blue", "green", "orange"),
                ltypes=c(1,1,1,1),
                main="KM plots with interaction model", xname="",
                yname='Progression Free Survival',
                level.grp=NULL,
                legend.head="Patient Group                       Median    HR (95\% CI)     P-value",
                legend.type = 4,
                legend.pos = c(0, 0),
                legend.adj = c(-0.05, -0.1),
                legend.text=10, legend.title=12,
                axis.text=12, axis.title=13,
                fontsize=4.5,
                round.hr=2, round.p=3, n.round=2,
                heights.plot = c(2, 0.85),
                margin.plot=c(1, 1, 0, 1),
                margin.tb=c(0, 0.45, 2, 0.35),
                xminlab = -2.5, xmaxlab = -1, xlim.lo = -1,
                n.space=list(c(16, 12, 21, 16, 6),
                             c(10, 12, 20, 10, 6)))


pdf(file=file.path(getwd(), "f_km_plot_pair.comp.pdf"), width=11, height=8.5)
grid.newpage()
print(p1, vp = viewport(x = unit(0.5, "npc"), y = unit(0.5, "npc"), width=0.95, height=0.9, just = "centre"))
grid.newpage()
print(p2, vp = viewport(x = unit(0.5, "npc"), y = unit(0.5, "npc"), width=0.95, height=0.9, just = "centre"))
dev.off()

}
