#' Diagnoses Codes for Thrombotic Events (ICD10)
#' 
#' This table contains ICD10 codes for thrombo-embolic events (e.g. pulmonary embolism, myocardial 
#' infarction, deep vein thromboembolism, arterial thromboembolism).
#' This table is an abbreviated version derived from the main table thromb_ref
#' 
#' The headers in the table are as follows:
#' 
#' - icd_code - shortened ICD code (see below)
#' - full_code - full ICD code (see below)
#' - type - diagnosis category as above
#' - description - full test description
#' 
#' Due to the variations in the way ICD codes may be stored, we have adhered to the following conventions:
#' - icd_code will store the broad category of the diagnosis
#' - Full ICD diagnosis codes will have the major category and subcategory separated by a period.
#' For example, the ICD-10 diagnosis "Septic pulmonary embolism with acute cor pulmonale" with full code I260.1
#' will have icd_code = "I260" and full_code = "I260.1".
#' Please ensure that your data tables store diagnoses codes as their broad categories. (e.g. I260.1 will be 
#' stored as "I260").
#' In future versions, this convention may be changed depending on the granularity required for analysis.
#' 
#' @name thromb_icd10_ref
#' 
#' @docType data
#' 
#' @usage data(thromb_icd10_ref)
#' 
#' @keywords datasets
#' 
"thromb_icd10_ref"