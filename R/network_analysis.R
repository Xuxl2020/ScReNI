###source('/data2/jwang/Retina/RetReg/CellType/Programs/network_analysis.R')

network_analysis <- function (regulatory_relationships, Kmeans_result, TFFDR1 = 10,
    TFFDR2 = 10, ModuleFDR = 0.05){
    
    source('/data2/jwang/Retina/RetReg/CellType/Programs/reconstruct_network_part.R')


    if (!"Correlation" %in% colnames(regulatory_relationships)) {
        stop("regulatory_relationships should contain Correlation column")
    }
    if (!"TF" %in% colnames(regulatory_relationships)) {
        stop("regulatory_relationships should contain TF column")
    }
    if (!"Target" %in% colnames(regulatory_relationships)) {
        stop("regulatory_relationships should contain Target column")
    }
    if (!"KmeansGroup" %in% colnames(Kmeans_result)) {
        stop("Kmeans_result should contain KmeansGroup column")
    }
    if (length(regulatory_relationships$TF[!regulatory_relationships$TF %in%
        rownames(Kmeans_result)]) > 0) {
        stop("rownames of Kmeans_result should contain all TF in regulatory_relationships")
    }
    if (length(regulatory_relationships$Target[!regulatory_relationships$Target %in%
        rownames(Kmeans_result)]) > 0) {
        stop("rownames of Kmeans_result should contain all Target in regulatory_relationships")
    }
    TFs_list <- get_Enriched_TFs(regulatory_relationships, Kmeans_result,
        TFFdrThr1 = TFFDR1)
    TFs_list <- get_regulation_of_TFs_to_modules(TFs_list, TFFDR2)
    TFs_list <- get_partial_regulations(TFs_list)
    TFs_list <- merge_Module_Regulations(TFs_list, Kmeans_result,
        ModuleThr1 = ModuleFDR)
    return(TFs_list)
}

