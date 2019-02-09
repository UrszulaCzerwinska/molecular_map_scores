GSE72056_melanoma_single_cell <- read.delim("../Data/GSE72056_melanoma_single_cell.txt", row.names=NULL)


GSE72056_melanoma_single_cell.w=data.frame(GSE72056_melanoma_single_cell)
headers = GSE72056_melanoma_single_cell.w[,1]
text = GSE72056_melanoma_single_cell.w[1:3,-1]
GSE72056_melanoma_single_cell.num = GSE72056_melanoma_single_cell.w[-c(1:3),-1]
GSE72056_melanoma_single_cell.num.t=t(GSE72056_melanoma_single_cell.num)
GSE72056_melanoma_single_cell.num.t.text = data.frame(t(text),GSE72056_melanoma_single_cell.num.t)
colnames(GSE72056_melanoma_single_cell.num.t.text) = headers


GSE72056_melanoma_single_cell.w.t=GSE72056_melanoma_single_cell.num.t.text
colnames(GSE72056_melanoma_single_cell.w.t)[2:3]=c("malignant","cell_type")


Macrophages.sinCell=subset(GSE72056_melanoma_single_cell.w.t,cell_type == "3" & malignant == "1" )
Macrophages.sinCell[1:10,1:10]

Macrophages.sinCell[c("CY88CD45POS_2_G03_S459_comb","cy53.1.CD45.pos.2.D10.S1006.comb","CY88CD45POS_2_D06_S426_comb","CY89A_CD45_POS_10_E06_S246_comb","cy53.1.CD45.pos.2.B12.S984.comb","CY88CD45POS_7_G04_S268_comb"), 1:10]
drp.rows=c("CY88CD45POS_2_G03_S459_comb","cy53.1.CD45.pos.2.D10.S1006.comb","CY88CD45POS_2_D06_S426_comb","CY89A_CD45_POS_10_E06_S246_comb","cy53.1.CD45.pos.2.B12.S984.comb","CY88CD45POS_7_G04_S268_comb")
Macrophages.sinCell.filtr1<- Macrophages.sinCell[!(row.names(Macrophages.sinCell) %in% drp.rows),]


Macrophages.sinCell.filtr1$tumor=as.factor(Macrophages.sinCell.filtr1$tumor)
Macrophages.sinCell.filtr1.n=Macrophages.sinCell.filtr1[,4:23689]


drop.rows2=c("cy60_1_cd_45_pos_4_C08_S32_comb","cy84_Primary_CD45_pos_C08_S416_comb")


Macrophages.sinCell.filtr2=Macrophages.sinCell.filtr1.n[!(row.names(Macrophages.sinCell.filtr1.n) %in% drop.rows2),]


global.mean=mean(as.matrix(Macrophages.sinCell.filtr2))
rows.mean=apply(as.matrix(Macrophages.sinCell.filtr2), 1, mean)
cols.mean=apply(as.matrix(Macrophages.sinCell.filtr2), 2, mean)


m1=apply(as.matrix(Macrophages.sinCell.filtr2), 2, function(x) x-rows.mean)
m2=apply(m1, 1, function(x) x-cols.mean)
m3=t(m2)+global.mean
m3[1:10,1:10]

Macrophages.dc=m3


drop.rows3=c("cy84_Primary_CD45_pos_E01_S433_comb","cy60_1_cd_45_pos_4_E01_S49_comb")
Macrophages.sinCell.filtr3=Macrophages.sinCell.filtr2[!(row.names(Macrophages.sinCell.filtr2) %in% drop.rows3),]
dim(Macrophages.sinCell.filtr2)
dim(Macrophages.sinCell.filtr3)

global.mean=mean(as.matrix(Macrophages.sinCell.filtr3))
rows.mean=apply(as.matrix(Macrophages.sinCell.filtr3), 1, mean)
cols.mean=apply(as.matrix(Macrophages.sinCell.filtr3), 2, mean)

m1=apply(as.matrix(Macrophages.sinCell.filtr3), 2, function(x) x-rows.mean)
m2=apply(m1, 1, function(x) x-cols.mean)
m3=t(m2)+global.mean

Macrophages.dc_2=m3
Macrophages.dc_2[1:10,1:10]

Macrophages.dc.t = t(Macrophages.dc_2)
rownames(Macrophages.dc.t) = headers[4:length(headers)]
Macrophages.dc.t = Macrophages.dc.t[!duplicated(rownames(Macrophages.dc.t)),]
# write.table(data.frame(GENES = rownames(Macrophages.dc.t), Macrophages.dc.t), "../Data/preprocessed2/macro_dc_outrem.txt", row.names=FALSE, col.names = TRUE, quote=FALSE, sep="\t")
