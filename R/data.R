#' Blood pressure data
#' 
#' Systolic and diastolic blood pressure was measured in a sample of 384 subjects using a handle mercury sphygmomanometer device and an automatic device. The blood pressure was simultaneously measured twice by each instrument, thus every subject had four measurements, two by each method. 
#' 
#'
#' @format A data frame with the following columns:
#' \describe{
#'   \item{ID}{Subject identifier}
#'   \item{SIS}{Systolic blood pressure in mmHg}
#'   \item{DIA}{Diastolic blood pressure in mmHg}
#'   \item{METODE}{Device identifier}
#'   \item{NM}{Identifier of replicates}
#'   \item{ALTURA}{Height in cm}
#'   \item{EDAD}{Age in years}
#'   \item{FRECUENC}{Heart rate}
#'   \item{INFOR_AR}{Have the subject been informed about he is hypertense?}
#'   \item{PESO}{Weight in Kg}
#'   \item{SEXO}{Gender. 1 for Male. 2 for Female}
#'   \item{TA}{Was the subject's blood pressure measured the last year? 1=Yes, 2=No, 9=Unknown}
#'   \item{TNSI_MED}{Does the subject receive treatment for hypertension? 1=Yes, 2=No, 3=Doubtful, 8=not applicable, 9=Insufficient data.}
#' }
"bpres"

#' Blood draw data
#' 
#' Plasma cortisol area under curve (AUC) was calculated from the trapezoidal rule over the 12-h period of the hourly blood draws. The subjects were required to repeat the process in five visits. The aim of the agreement study was to assess how well the plasma cortisol AUC from hourly measurements agreed with plasma cortisol AUC that was measured every two hours.
#' 
#'
#' @format A data frame with the following columns:
#' \describe{
#'   \item{SUBJ}{Subject identifier}
#'   \item{VNUM}{Visit number}
#'   \item{AUC}{Area under the curve}
#'   \item{MET}{Device identifier}
#' }
"bdaw"


#' Body fat data
#' 
#' Percentage body fat was estimated from skinfold calipers and DEXA on a cohort of 90 adolescent girls. Skinfold caliper and DEXA measurements were taken at ages 12.5, 13 and 13.5. The objective was to determine the amount of agreement between the skinfold caliper and DEXA measurements of percentage body fat.
#' 
#'
#' @format A data frame with the following columns:
#' \describe{
#'   \item{SUBJECT}{Subject identifier}
#'   \item{VISITNO}{Visit number}
#'   \item{BF}{Percentage body fat}
#'   \item{MET}{Device identifier}
#' }
"bfat"



