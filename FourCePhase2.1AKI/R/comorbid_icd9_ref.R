#' Diagnoses Codes for Common Comorbidities and Thrombotic Events (ICD10)
#' 
#' This table contains ICD10 codes for diagnoses pertaining to the following broad categories:
#' 1) Cardiovascular Risk Factors (Diabetes mellitus, hypertension, hyperlipidemia, ischemic heart
#'    disease)
#' 2) Chronic kidney disease / End-stage Renal Failure
#' 3) Malignancies (solid tumours, myeloproliferative and lymphoproliferative disorders)
#' 3) Thrombo-embolic events (e.g. pulmonary embolism, myocardial infarction, deep vein thromboembolism,
#'    arterial thromboembolism)
#' 4) Pulmonary diseases (e.g. asthma, chronic obstructive pulmonary disease, pulmonary tuberculosis, 
#'    interstitial lung diseases)
#' 5) Rheumatological conditions (e.g. rheumatoid arthritis, systemic lupus erythematosus, systemic
#'    sclerosis, dermatomyositis, polymyositis, Sjogren's syndrome)
#'    
#' This table is a truncated version derived from the main table comorbid_ref.
#' 
#' The headers in the table are as follows:
#' 
#' - icd_code - shortened ICD code (see below)
#' - full_code - full ICD code (see below)
#' - comorbid_type - diagnosis category as above
#' 
#' Due to the variations in the way ICD codes may be stored, we have adhered to the following conventions:
#' - icd_code will store the broad category of the diagnosis
#' - Full ICD diagnosis codes will have the major category and subcategory separated by a period.
#' For example, the ICD-9 diagnosis "Chronic kidney disease, Stage V" with full code 585.5 will have 
#' icd_code = "585" and full_code = "585.5".
#' Please ensure that your data tables store diagnoses codes as their broad categories. (e.g. 585.5 will be 
#' stored as "585").
#' In future versions, this convention may be changed depending on the granularity required for analysis.
#' 
#' @name comorbid_icd9_ref
#' 
#' @docType data
#' 
#' @usage data(comorbid_icd9_ref)
#' 
#' @keywords datasets
#' 
"comorbid_icd9_ref"