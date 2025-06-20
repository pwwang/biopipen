library(SeuratWrappers)
library(Seurat)
library(purrr)
library(stringr)
library(biopipen.utils)

infile = {{in.infile | r}}
outfile = {{out.outfile | r}}
envs = {{envs.alra_args | r}}

log <- get_logger()

log$info("Loading Seurat object")
sobj <- read_obj(infile)
assay <- DefaultAssay(sobj)

# https://github.com/mojaveazure/seurat-disk/issues/102
# https://github.com/simoncmo/shared_seurat_scripts/blob/main/function_seurat_janitor.R
# Try to fix the issue with SCTModel
log$info("Trying to fix SCTModel issue (see mojaveazure/seurat-disk#102)")
# --------------------------------------------------------------------------
# Handle missing median_umi
fix_median_umi = function(SCTModel_obj){
	err_message = ''
	tryCatch({ test <- methods::validObject(SCTModel_obj) },
		error = function(error_message) {
			err_message <<- as.character(error_message)
	})
	missing_medium_umi = stringr::str_detect(err_message, 'median_umi')

	if(missing_medium_umi){
		message('Missing medium_umi, calculate again from cell.attributes$umi')
		slot(SCTModel_obj, 'median_umi') = median(SCTModel_obj@cell.attributes$umi)
	}
	return(SCTModel_obj)
}

# Cleaning empty objects
# General purpose
clean_seurat_obj_list = function(obj_list, attirbute_to_check){
	if(missing(attirbute_to_check)) {stop("Need attributes to check for cleaning")}
	# Object type
	obj_type = class(obj_list[[1]])[[1]]

	# Count
	obj_size = unlist(purrr::map(obj_list, function(object){
		nrow(slot(object, attirbute_to_check))
	}))

	# Remove empty
	if(length(obj_size ==0)  != 0 ){
		message(str_glue('Removing {length(obj_size ==0)} empty object from the {obj_type} object list'))
		obj_list = obj_list[obj_size!=0]
		message(str_glue('{length(obj_list)} {obj_type} object(s) left'))
	}
	obj_list
}

# for SCTModel.list slot
clean_seurat_SCTModel_list = function(sct_model_list){
	clean_seurat_obj_list(obj_list = sct_model_list, attirbute_to_check = 'cell.attributes')
}

fix_seurat_SCT = function(obj){
	# Check first
	if(!'SCT' %in% Assays(obj)){
		message('SCT assay not found. Nothing to fix')
		return(obj)
	}

	# Model list
	sct_model_list = obj$SCT@SCTModel.list
	# 1. clean SCTModel list
	sct_model_list = clean_seurat_SCTModel_list(sct_model_list)

	# 2. fix missing median_umi
	sct_model_list = map(sct_model_list, function(sct_model){
		fix_median_umi(sct_model)
	})

	# Add back and retrun
	obj$SCT@SCTModel.list = sct_model_list

	return(obj)
}
# --------------------------------------------------------------------------
sobj = fix_seurat_SCT(sobj)

log$info("Imputing expression values, using ALRA")
envs$object <- sobj
sobj = do_call(RunALRA, envs)
envs$object <- NULL
gc()

log$info("Renaming assays")
sobj = RenameAssays(sobj, assay.name = assay, new.assay.name = "RAW")
sobj = RenameAssays(sobj, assay.name = "alra", new.assay.name = assay)
DefaultAssay(sobj) <- assay

sobj@misc$impute_method = "alra"

log$info("Saving Seurat object")
save_obj(sobj, outfile)

# choosek_plot_file = file.path(dirname(outfile), "choosek.png")
# png(choosek_plot_file, width = 1200, height = 1000, res = 100)
# p = ALRAChooseKPlot(sobj)
# print(p)
# dev.off()
