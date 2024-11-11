# Probleme : Calcul du Transfert de Chaleur

# Contexte
# Nous avons une paroi plane homogene de materiau isolant avec une conductivite thermique constante k.
# La temperature a l'interieur de la paroi varie le long de son epaisseur en reponse a une difference
# de temperature entre les deux surfaces. Nous voulons calculer le flux de chaleur qui traverse cette paroi.

# Modele Mathematique
# Le flux de chaleur q a travers une paroi plane est donne par la loi de Fourier pour la conduction thermique :
# q = -k * dT/dx
# ou :
# k est la conductivite thermique du materiau (W/m·K)
# T est la temperature (K)
# x est la position a travers l'epaisseur de la paroi (m)

# Pour une paroi de longueur L avec une temperature T1 d'un cote et T2 de l'autre,
# nous devons integrer la temperature pour determiner le gradient de temperature a chaque point de la paroi
# et en deduire le flux thermique global.

# Etapes
# 1. Definir les parametres du probleme :
#    - Longueur de la paroi L
#    - Temperatures T1 et T2
#    - Conductivite thermique k

# 2. Formuler le gradient de temperature :
#    - Si la temperature varie lineairement entre T1 et T2, nous pouvons approximativement modeliser
#      la temperature le long de x comme :
#      T(x) = T1 + ((T2 - T1) / L) * x

# 3. Calculer l'integrale :
#    - Integrer la temperature T(x) pour trouver le flux thermique total.

# Define the parameters
k <- 0.5  # Conductivity in W/m·K
L <- 0.1  # Length of the wall in meters
T1 <- 350  # Temperature at x=0 in Kelvin
T2 <- 300  # Temperature at x=L in Kelvin

# Define the temperature function T(x)
T <- function(x) T1 + ((T2 - T1) / L) * x

# Define the derivative function dT/dx
dTdx <- function(x) (T2 - T1) / L

# Number of subdivisions
n <- 10



# Calculate the heat flux using different integration methods
heat_flux_trapeze <- (-k * integrationtrapeze(dTdx, 0, L, n))
heat_flux_gauss <- -k * integrationgauss(dTdx, 0, L, 5)

# Print the results
cat("Heat flux using Trapezoidal Rule:", heat_flux_trapeze, "W\n")
cat("Heat flux using Gauss-Legendre Method:", heat_flux_gauss, "W\n")
