library(ggplot2)
library(RColorBrewer)

source("code/experimentFunctions.R")

#####

experiment <- "yeast_random"

makeAUClinePlotYeastTCPA(experiment, c(0.4,0.8))
ggsave(file.path("output", experiment, "plots", paste0(experiment, "_AUC.pdf")), width = 15, height = 15, units = "cm")

makeROCcurvePlotYeastTCPA(experiment, c(0,1), c(0,1), pointMethodSet = "final")
ggsave(file.path("output", experiment, "plots", paste0(experiment, "_ROC.pdf")), width = 15, height = 22.5, units = "cm")

#####

experiment <- "yeast_row-wise"

makeAUClinePlotYeastTCPA(experiment, c(0.4,0.8))
ggsave(file.path("output", experiment, "plots", paste0(experiment, "_AUC.pdf")), width = 15, height = 15, units = "cm")

makeROCcurvePlotYeastTCPA(experiment, c(0,1), c(0,1), pointMethodSet = "final")
ggsave(file.path("output", experiment, "plots", paste0(experiment, "_ROC.pdf")), width = 15, height = 22.5, units = "cm")

#####

experiment <- "yeast_row-wise_GIES"

makeROCcurvePlotYeastTCPA(experiment, c(0,1), c(0,1), pointMethodSet = "final")
ggsave(file.path("output", experiment, "plots", paste0(experiment, "_ROC.pdf")), width = 15, height = 22.5, units = "cm")

#####

experiment <- "TCPA"

makeAUClinePlotYeastTCPA(experiment, c(0.4,0.8))
ggsave(file.path("output", experiment, "plots", paste0(experiment, "_AUC.pdf")), width = 15, height = 15, units = "cm")

makeROCcurvePlotYeastTCPA(experiment, c(0,1), c(0,1), pointMethodSet = "final")
ggsave(file.path("output", experiment, "plots", paste0(experiment, "_ROC.pdf")), width = 15, height = 22.5, units = "cm")

#####

experiment <- "cellLine"

makeAUCbarPlotCellLine(c(0, 0.95))
ggsave(file.path("output", experiment, "plots", paste0(experiment, "_AUC.pdf")), width=15, height=15, units = "cm")
