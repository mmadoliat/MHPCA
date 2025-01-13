#' @title Define an Image Class for Multidimensional Non Functional Data
#'
#' @description
#' The `img` class represents image data as a specialized form of `nfd`. It can
#' convert images from the `imager` package into the `img` class, ensuring they
#' are stored as matrices or arrays with the appropriate class attributes.
#'
#' @param image An image object, preferably from the `imager` package.
#'
#' @return An object of class `img`, `nfd`, and the original class of the image.
#'
#' @export
#'
#' @examples
#' library(imager)
#' img_cimg <- load.example("parrots")  # Example image from imager
#' img_obj <- img(img_cimg)
#' class(img_obj)
img <- function(image) {
  # Check if the input is from the 'imager' package
  if (inherits(image, "cimg")) {
    # Convert 'cimg' object to an array
    image_data <- as.array(image)
  } else if (is.matrix(image) || is.array(image)) {
    image_data <- image
  } else {
    stop("Input must be a matrix, array, or an image object from the 'imager' package.")
  }
  
  # Assign the appropriate classes
  class(image_data) <- c("img", "nfd", class(image_data))
  
  return(image_data)
}
