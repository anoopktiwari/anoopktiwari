#################################
#################################



# taking a quick look of columns of df
glimpse(df)

# to remove perticular values from data frame
output_name <- df[apply(df, 1, function(row) all(row !=1)), ]
                        
# to write df as required format
write.table(genes_counts1, file = "15k_genecount.gct", row.names=T, sep="\t", quote=F)
                        
                        
# to convert row.names into columns
cleared_tpm2 <- tibble::rownames_to_column(cleared_tpm, "Name")

# add blank column                        
add_column(cleared_tpm2, add_column="Dummy") -> cleared_tpm2                       
                        
# relocate columns by their names
cleared_tpm2 %>% relocate(add_column, .after = Name) -> cleared_tpm2
                        
                        
                        
