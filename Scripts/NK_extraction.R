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




NK.sinCell=subset(GSE72056_melanoma_single_cell.w.t,cell_type == "6" & malignant == "1" )
dim(NK.sinCell)



#NK.sinCell[c("cy80.Cd45.pos.Pd1.neg.S293.E05.S293.comb","cy72.CD45.pos.H08.S956.comb","cy82.CD45.pos.3.A07.S7.comb","Cy67.CD45pos.S2.C4_S28"),1:10]

drp.rows=c("cy80.Cd45.pos.Pd1.neg.S293.E05.S293.comb","cy53.1.CD45.pos.2.D10.S1006.comb","CY88CD45POS_2_D06_S426_comb","CY89A_CD45_POS_10_E06_S246_comb","cy53.1.CD45.pos.2.B12.S984.comb","CY88CD45POS_7_G04_S268_comb")
NK.sinCell.filtr1<- NK.sinCell[!(row.names(NK.sinCell) %in% drp.rows),]


NK.sinCell.filtr1$tumor=as.factor(NK.sinCell.filtr1$tumor)
NK.sinCell.filtr1.n=NK.sinCell.filtr1[,4:23689]

drp.rows2=c("cy60_1_cd_45_pos_4_E08_S56_comb","cy84_Primary_CD45_pos_E08_S440_comb","cy84_Primary_CD45_pos_C07_S415_comb","cy60_1_cd_45_pos_4_C07_S31_comb","cy60_1_cd_45_pos_4_G03_S75_comb","cy84_Primary_CD45_pos_G03_S459_comb","cy60_1_cd_45_pos_4_C11_S35_comb","cy84_Primary_CD45_pos_C11_S419_comb")
NK.sinCell.filtr2<- NK.sinCell.filtr1.n[!(row.names(NK.sinCell.filtr1.n) %in% drp.rows2),]

#double center
global.mean=mean(as.matrix(NK.sinCell.filtr2))
rows.mean=apply(as.matrix(NK.sinCell.filtr2), 1, mean)
cols.mean=apply(as.matrix(NK.sinCell.filtr2), 2, mean)

length(rows.mean)
length(cols.mean)

m1=apply(as.matrix(NK.sinCell.filtr2), 2, function(x) x-rows.mean)
m2=apply(m1, 1, function(x) x-cols.mean)
m3=t(m2)+global.mean

NK.dc=m3


# write.table(t(NK.dc), "nk_dc_outrem.txt", row.names=TRUE, col.names = TRUE, quote=FALSE, sep="\t")
# write.table(t(NK.dc), "nk_dc_outrem_num.txt", row.names=FALSE, col.names = FALSE, quote=FALSE, sep="\t")
# write.table(row.names(t(NK.dc)), "nk_dc_outrem_ids.txt", row.names=FALSE, col.names = FALSE, quote=FALSE, sep="\t")
# write.table(row.names(NK.dc), "nk_dc_outrem_samples.txt", row.names=FALSE, col.names = FALSE, quote=FALSE, sep="\t")


NK.dc.t = t(NK.dc)
rownames(NK.dc.t) = headers[4:length(headers)]
NK.dc.t = NK.dc.t[!duplicated(rownames(NK.dc.t)),]

#write.table(data.frame(GENES = rownames(NK.dc.t),NK.dc.t), "../Data/preprocessed2/nk_dc_outrem.txt", row.names=FALSE, col.names = TRUE, quote=FALSE, sep="\t")
