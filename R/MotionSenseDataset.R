#' Motion Sense Dataset
#'
#' Measurements of user acceleration and pitch attitude collected by
#' smartphone accelerometer and gyroscope sensors from 24 participants
#' performing four activities.
#'
#' @format A list with two numeric matrices, each with 200 rows (time points)
#'   and 96 columns (24 participants x 4 activities):
#' \describe{
#'   \item{user_acceleration}{user acceleration}
#'   \item{pitch_attitude}{pitch attitude}
#' }
#' @docType data
#' @keywords datasets
#' @name motion_sense_data
#' @usage data("motion_sense_data")
#' @examples
#' data("motion_sense_data")
#' dim(motion_sense_data$user_acceleration)
NULL

#' Electrical Power Dataset
#'
#' Daily curves of household electric power and voltage, starting on
#' 17 December 2006, with 288 readings per day (5-minute intervals).
#'
#' @format A list with two data frames, `power` and `voltage`, each with 1440
#'   rows (days) and 291 columns: `year`, `month` and `day`, followed by the
#'   288 readings of that day in columns `X1` to `X288`.
#' @docType data
#' @keywords datasets
#' @name electrical_power_data
#' @usage data("electrical_power_data")
#' @examples
#' data("electrical_power_data")
#' dim(electrical_power_data$power)
NULL
