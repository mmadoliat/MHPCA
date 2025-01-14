#' @title Define a Grayscale Image Class for Multidimensional Non-Functional Data
#'
#' @description
#' The `img` class represents grayscale image data as a specialized form of `nfd`.
#' It converts one or multiple 2D matrices or grayscale images from the `imager` package into
#' a 2-dimensional matrix where each row corresponds to a vectorized image.
#'
#' @param ... A grayscale image object (from `imager` package as `cimg`, or a 2D matrix/array),
#'            or a single list containing such objects.
#'
#' @return An object of class `img`, `nfd`, `matrix`, and `array`. If multiple images are provided,
#'         the returned matrix will have multiple rows, each representing a vectorized image.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(imager)
#'
#' # Single image as separate argument
#' img_cimg1 <- load.example("lena") %>% grayscale()
#' img_obj1 <- img(img_cimg1)
#' print(img_obj1)
#' class(img_obj1)
#'
#' # Multiple images as separate arguments
#' img_cimg2 <- load.example("boats") %>% grayscale()
#' img_cimg3 <- load.example("peppers") %>% grayscale()
#' img_obj2 <- img(img_cimg1, img_cimg2, img_cimg3)
#' print(img_obj2)
#' class(img_obj2)
#'
#' # Multiple images as a list
#' img_list <- list(img_cimg1, img_cimg2, img_cimg3)
#' img_obj3 <- img(img_list)
#' print(img_obj3)
#' class(img_obj3)
#' }
img <- function(...) {
  # Capture all arguments
  args <- list(...)
  
  # If a single list is provided, treat it as list of images
  if (length(args) == 1 && is.list(args[[1]])) {
    img_list <- args[[1]]
  } else {
    img_list <- args
  }
  
  # Check that there is at least one image
  if (length(img_list) == 0) {
    stop("At least one image must be provided to the 'img' function.")
  }
  
  vectorized_images <- list()
  original_dim <- NULL
  
  # Iterate over each image
  for (i in seq_along(img_list)) {
    image <- img_list[[i]]
    
    # Handle input based on its class
    if (inherits(image, "cimg")) {
      dims <- dim(image)
      if (length(dims) != 4 || dims[3] != 1 || dims[4] != 1) {
        stop(paste0("Input 'cimg' object at position ", i, " must be a 2D grayscale image with a single channel."))
      }
      # Convert 'cimg' to a matrix
      image_data <- as.matrix(image)
    } else if (is.matrix(image)) {
      # Ensure the matrix is 2D
      if (length(dim(image)) != 2) {
        stop(paste0("Matrix input at position ", i, " must be 2-dimensional for grayscale images."))
      }
      image_data <- image
    } else if (is.array(image)) {
      # Ensure the array is 2D for grayscale
      if (length(dim(image)) != 2) {
        stop(paste0("Array input at position ", i, " must be 2-dimensional for grayscale images."))
      }
      image_data <- image
    } else {
      stop(paste0("Input at position ", i, " must be a 2D matrix, a 2D array, or an image object from the 'imager' package."))
    }
    
    # Get original dimensions
    current_dim <- dim(image_data)
    
    # On first image, store the dimensions
    if (is.null(original_dim)) {
      original_dim <- current_dim
    } else {
      # Check if current image has the same dimensions as previous images
      if (!all(current_dim == original_dim)) {
        stop(paste0("All images must have the same dimensions. Image at position ", i, 
                    " has dimensions ", paste(current_dim, collapse = "x"), 
                    " which differs from previous images with dimensions ", 
                    paste(original_dim, collapse = "x"), "."))
      }
    }
    
    # Vectorize the image: 1 x (height * width)
    height <- current_dim[1]
    width <- current_dim[2]
    num_pixels <- height * width
    vectorized_data <- matrix(as.vector(image_data), nrow = 1, ncol = num_pixels)
    
    # Assign column names if it's the first image
    if (i == 1) {
      colnames(vectorized_data) <- paste0("Pixel_", 1:num_pixels)
    }
    
    # Assign row name
    rownames(vectorized_data) <- paste0("Image_", i)
    
    # Append to list
    vectorized_images[[i]] <- vectorized_data
  }
  
  # Combine all vectorized images into a single matrix
  combined_matrix <- do.call(rbind, vectorized_images)
  
  # Assign the appropriate classes
  class(combined_matrix) <- c("img", "nfd", class(combined_matrix))
  
  # Attach original dimensions as an attribute
  attr(combined_matrix, "original_dimensions") <- original_dim
  
  return(combined_matrix)
}