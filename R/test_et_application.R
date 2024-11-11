# Define the function to integrate
f <- function(x) exp(-x^2)

# Integration bounds
a <- 0
b <- 1

# Number of subdivisions
n <- 10

# Calculate integrals using different methods
integral_trapeze <- integrationtrapeze(f, a, b, n)
integral_simpson <- integrationsimpson(f, a, b, n)
# Boole's method requires n to be a multiple of 4
n_boole <- 12  # Adjust n to be multiple of 4
integral_boole <- integrationboole(f, a, b, n_boole)
integral_gauss <- integrationgauss(f, a, b, 5)

# Print the results
cat("Integral using Trapezoidal Rule:", integral_trapeze, "\n")
cat("Integral using Simpson's Rule:", integral_simpson, "\n")
cat("Integral using Boole's Rule:", integral_boole, "\n")
cat("Integral using Gauss-Legendre Method:", integral_gauss, "\n")

