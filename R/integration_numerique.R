# Numerical integration using the trapezoidal rule
#' Title Numerical Integration using the Trapezoidal Rule
#' @param f function to integrate
#' @param a lower bound of integration
#' @param b upper bound of integration
#' @param n number of subdivisions
#' @return value of the integral
#' @export
#'
#' @examples
#' f <- function(x) 3*(x^3)
#' integrationtrapeze(f,0,1,10)
integrationtrapeze<-function(f,a,b,n){
  # Calculate the width of the subintervals
  largdesousinter<-(b-a)/n
  # Sum the function values at the points in the interval
  somdesval<-(f(a)+f(b))/2
  for (i in 1:(n-1)) {
    somdesval<-somdesval+f(a+i*largdesousinter)
  }
  # Calculate the integral by multiplying the sum of values by the width of the subintervals
  return(largdesousinter*somdesval)
}

#Numerical integration using Simpsons rule
#' Title Numerical Integration using Simpsons Rule
#'
#' @param f function to integrate
#' @param a lower bound of integration
#' @param b upper bound of integration
#' @param n number of subdivisions (must be even)
#'
#' @return The approximate value of the integral
#' @export
#'
#' @examples
#' f <- function(x) sin(x)
#' integrationsimpson(f,0,pi,32)
integrationsimpson <- function(f, a, b, n) {
  # Check if n is even
  if (n %% 2 != 0) stop("n doit etre pair")
  # Calculate the width of the subdivisions
  largdessubd <- (b - a) / n
  # Division points
  ptsdediv <- seq(a, b, length.out = n + 1)
  y <- f(ptsdediv)
  # Calculate the integral using Simpsons rule
  integral <- largdessubd / 3 * (y[1] + y[n + 1] + 4 * sum(y[seq(2, n, by = 2)]) + 2 * sum(y[seq(3, n - 1, by = 2)]))
  return(integral)
}

#Integration numerique par la methode de Boole
#' Title Integration numerique par la methode de Boole
#'
#' @param f function to integrate
#' @param a lower bound of integration
#' @param b upper bound of integration
#' @param n  number of subdivisions
#'
#' @return The approximate value of the integral
#' @export
#'
#' @examples
#' f <- function(x) x^2 * sin(x)
#' resultat <- integrationboole(f, 0, pi, 16)
integrationboole <- function(f, a, b, n) {
  # Check if n is a multiple of 4
  if (n %% 4 != 0) stop("n doit etre multiple de 4")
  # Calculate the width of the subintervals
  h <- (b - a) / n
  # Sequence of division points
  x <- seq(a, b, by = h)
  # Evaluate the function at the division points
  y <- f(x)
  # Calculate the integral using Booles rule
  integral <- 2 * h / 45 * (7 * sum(y[c(1, n + 1)]) + 32 * sum(y[seq(2, n, by = 2)]) + 12 * sum(y[seq(3, n - 1, by = 4)]) + 14 * sum(y[seq(5, n - 3, by = 4)]))
  return(integral)
}

#Numerical integration using the Gauss-Legendre method
#Function to obtain the Gauss-Legendre nodes and weights for a given order
gaussLegendreNodesWeights <- function(n) {
  # x represents the vector of abscissas of the Gauss points (nodes) on the interval [-1,1]
  # w represents the vector of weights associated with each node x
  if (n == 2) {
    x <- c(-0.5773502692, 0.5773502692)
    w <- c(1.0, 1.0)
  } else if (n == 3) {
    x <- c(-0.7745966692, 0.0, 0.7745966692)
    w <- c(0.5555555556, 0.8888888889, 0.5555555556)
  } else if (n == 4) {
    x <- c(-0.8611363116, -0.3399810436, 0.3399810436, 0.8611363116)
    w <- c(0.3478548451, 0.6521451549, 0.6521451549, 0.3478548451)
  } else if (n == 5) {
    x <- c(-0.9061798459, -0.5384693101, 0.0, 0.5384693101, 0.9061798459)
    w <- c(0.2369268851, 0.4786286705, 0.5688888889, 0.4786286705, 0.2369268851)
  } else {
    stop("Order n must be 2, 3, 4, or 5")
  }
  return(list(x = x, w = w))
}

# Function for Gauss-Legendre integration for a given order
#' Title Numerical Integration using the Gauss-Legendre Method
#'
#' @param f function to integrate
#' @param a lower bound of integration
#' @param b upper bound of integration
#' @param n  number of subdivisions
#'
#' @return The approximate value of the integral
#' @export
#'
#' @examples
#' f <- function(x) exp(-x^2)
#' resultat <- integrationgauss(f, -1, 1, 5)
integrationgauss <- function(f, a, b, n) {
  # Obtain the nodes and weights for the given order n
  gauss <- gaussLegendreNodesWeights(n)
  x <- gauss$x
  w <- gauss$w
  # Initialize the integral
  integral <- 0
  # Calculate the integral using the Gauss-Legendre nodes and weights
  for (i in 1:n) {
    integral <- integral + w[i] * f(((b - a) * x[i] + b + a) / 2)
  }
  # Final adjustment for the interval [a, b]
  integral <- (b - a) / 2 * integral
  return(integral)
}
