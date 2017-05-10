last <- function(x) { 
	return( x[length(x)] ) 
}

first_entry <- function(x){
	return(x[1])
}

dendrogram_ordering <- function(dataset, weights, filename) {

				hc <- hclust(dist(dataset),method="single")
				dd <- as.dendrogram(hc)
				dd.reorder <- reorder(dd,weights,agglo.FUN=first_entry)
				the_labels <- labels(dd.reorder)
				#print(colnames(drug_data_key)[metric])
				write(the_labels,file=filename,ncolumns=1)
				#labels(dd.reorder) <- leaf_names[the_labels]
}


order_this_dendrogram <- function() {
	require(graphics)
		library(dendextend)

		drug_data <- read.table('/home/scratch/testoutput/Tox_Res_Paper/collated_data.tsv')
		drug_data_key <- read.csv('/home/scratch/testoutput/Tox_Res_Paper/collated_data_key.csv')
		scaled <- read.table('/home/scratch/testoutput/Tox_Res_Paper/scaled_metrics.tsv')
		weights <- as.vector(drug_data[,2])
		drug_names <- drug_data[,1]
		leaf_names <- paste(drug_names,weights)
		dendrogram_ordering(scaled[,1],weights, 'order_APD90')
		dendrogram_ordering(scaled[,6],weights, 'order_INa_EADs')
		dendrogram_ordering(scaled[,7],weights, 'order_ICaL_EADs')
		dendrogram_ordering(scaled[,8],weights, 'order_IKr_EADs')
		dendrogram_ordering(scaled[,c(1,6,7,8)],weights, 'order_APD90_EADs')
		dendrogram_ordering(scaled[,c(6,7,8)],weights, 'order_EADs')
		dendrogram_ordering(scaled[,c(1,6)],weights, 'order_APD90_ICaL')
		dendrogram_ordering(scaled[,c(2,3)],weights, 'order_Grandi_LS')
		dendrogram_ordering(scaled[,c(4,5)],weights, 'order_OHara_LS')
		dendrogram_ordering(drug_data[,11],weights, 'order_hERG_cmax')
}
