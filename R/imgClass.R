#' @title Define a Grayscale Image Class for Multidimensional Non-Functional Data
#'
#' @description
#' The `img` class represents grayscale image data as a specialized form of `nfd`.
#' It converts a 2D matrix or a grayscale image from the `imager` package into a
#' 2-dimensional matrix where each row corresponds to a pixel and the single column
#' represents the intensity value.
#'
#' @param image A grayscale image object, preferably from the `imager` package (`cimg` class), or a 2D matrix/array.
#'
#' @return An object of class `img`, `nfd`, and the original class of the image.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(imager)
#' img_cimg <- load.example("lena")  # Example grayscale image from imager
#' img_obj <- img(img_cimg)
#' print(img_obj)
#' class(img_obj)
#' }
img <- function(image) {
  
  # Handle input based on its class
  if (inherits(image, "cimg")) {
    # Ensure the image is grayscale
    if (dim(image)[4] != 1) {
      stop("Input 'cimg' object must be a grayscale image with a single channel.")
    }
    # Convert 'cimg' to a matrix
    image_data <- as.matrix(image)
  } else if (is.matrix(image)) {
    image_data <- image
  } else if (is.array(image)) {
    # Ensure the array is 2D for grayscale
    if (length(dim(image)) != 2) {
      stop("Array input must be 2-dimensional for grayscale images.")
    }
    image_data <- image
  } else {
    stop("Input must be a 2D matrix, a 2D array, or an image object from the 'imager' package.")
  }
  
  # Get original dimensions
  orig_dim <- dim(image_data)
  
  # Vectorize the image: (height * width) x 1
  height <- orig_dim[1]
  width <- orig_dim[2]
  vectorized_data <- matrix(as.vector(image_data), nrow = 1, ncol = height * width)
  
  # Assign column name
  rownames(vectorized_data) <- "Intensity"
  
  # Assign the appropriate classes
  class(vectorized_data) <- c("img", "nfd", class(vectorized_data))
  
  # Attach original dimensions as an attribute
  attr(vectorized_data, "original_dimensions") <- orig_dim
  
  return(vectorized_data)
}