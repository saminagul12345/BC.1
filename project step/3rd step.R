
 #3rd step

rm(list = ls())

options(stringsAsFactors = F)

suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(stringr))
suppressMessages(library(ggsci))
suppressMessages(library(dplyr))
suppressMessages(library(circlize))


## Load data
load('rdata/step_1_TCGA_EXP.Rdata')
load('rdata/step_1_TCGA_clinc.Rdata')

##  mRNAsi and clinical features
names(TCGA_clinc)
TCGA_mRNAsi <- TCGA_clinc[order(TCGA_clinc$mRNAsi),]
ha1 = HeatmapAnnotation(df = data.frame(Age = TCGA_mRNAsi$Age,
                                        Stage = TCGA_mRNAsi$Stage,
                                        OS = TCGA_mRNAsi$Status),
                        col = list(Age = colorRamp2(c(30, 90), c("white", "firebrick")),
                                   Stage = c('I'='#FFCC99','II'='#CC9999','III'='#996699','IV'='#993333'),
                                   OS = c('Alive' = "#99CC66", 'Dead' = "#336699",'NA' = 'grey')),
                                   annotation_height = unit(c(0.5, 0.5), "cm"))
ha2 = HeatmapAnnotation(barplot = anno_barplot(TCGA_mRNAsi$mRNAsi,bar_width = 1,
                                               gp = gpar(fill = 'blue',col = 'blue')),
                        annotation_height = unit(4, "cm"))

zero_row_mat = matrix(nrow = 0, ncol = nrow(TCGA_mRNAsi))
colnames(zero_row_mat) = TCGA_mRNAsi$submitter

ht = Heatmap(zero_row_mat, top_annotation = ha2, bottom_annotation = ha1,
             show_column_names = F,show_row_names = F)
draw(ht, padding = unit(c(2, 20, 2, 2), "mm"))

pdf(file="fig/step_2_mRNAsi_ann.pdf",width = 7,height = 4)
draw(ht, padding = unit(c(2, 20, 2, 2), "mm"))
decorate_annotation("Age", {grid.text("Age", unit(-2, "mm"), just = "right")})
decorate_annotation("OS", {grid.text("OS", unit(-2, "mm"), just = "right")})
decorate_annotation("Stage", {grid.text("Stage", unit(-2, "mm"), just = "right")})
decorate_annotation("barplot", {
  grid.text("mRNAsi", unit(-10, "mm"), just = "bottom", rot = 90)
})
dev.off()


## Add clinical data to the immune score
TCGA_imm <- TCGA_clinc[order(TCGA_clinc$ImmuneScore),]
names(TCGA_imm)
ha1 = HeatmapAnnotation(df = data.frame(Age = TCGA_clinc$Age,
                                        Stage = TCGA_clinc$Stage,
                                        OS = TCGA_clinc$State),
                        col = list(Age = colorRamp2(c(30, 90), c("white", "firebrick")),
                                   Stage = c('I'='#FFCC99','II'='#CC9999','III'='#996699','IV'='#993333'),
                                   OS = c('Alive' = "#99CC66", 'Dead' = "#336699",'NA' = 'grey')),
                        annotation_height = unit(c(0.5, 0.5), "cm"))
ha2 = HeatmapAnnotation(barplot = anno_barplot(TCGA_imm$ImmuneScore,bar_width = 1,
                                               gp = gpar(fill = 'blue',col = 'blue')),
                        annotation_height = unit(4, "cm"))

zero_row_mat = matrix(nrow = 0, ncol = nrow(TCGA_imm))
colnames(zero_row_mat) = TCGA_imm$submitter

ht = Heatmap(zero_row_mat, top_annotation = ha2, bottom_annotation = ha1,
             show_column_names = F,show_row_names = F)
draw(ht, padding = unit(c(2, 20, 2, 2), "mm"))

pdf(file="fig/step_2_immunescore_ann.pdf",width = 7,height = 4)
draw(ht, padding = unit(c(2, 20, 2, 2), "mm"))
decorate_annotation("Age", {grid.text("Age", unit(-2, "mm"), just = "right")})
decorate_annotation("OS", {grid.text("OS", unit(-2, "mm"), just = "right")})
decorate_annotation("Stage", {grid.text("Stage", unit(-2, "mm"), just = "right")})
decorate_annotation("barplot", {
  grid.text("ImmuneScore", unit(-10, "mm"), just = "bottom", rot = 90)
})
dev.off()


## OCLR with clinical features
names(TCGA_clinc)
table(TCGA_clinc$Status)


for (k in c('mRNAsi','ImmuneScore')) {
  for (i in c('Status','Age2','Stage2')) {
    color1 = c('#ff9900','#146eb4')
    order1 <- names(table(TCGA_clinc[,i]))
    comparison = list(order1)
    
    boxp2 = ggboxplot(TCGA_clinc, size = 1, bxp.errorbar = F,
                      x=i, y=k,  fill = i,
                      order = order1) + 
      stat_compare_means(comparisons = comparison, method = "wilcox.test") +
      theme(text = element_text(),legend.position="none")
    print(boxp2)
    pdf(file=paste0("fig/step_2_",k,"_",i,".pdf"),width = 3,height = 4)
    print(boxp2)
    dev.off()
  }
}
