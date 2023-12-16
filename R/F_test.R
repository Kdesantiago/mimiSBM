#' moyenne d’un vecteur
#' Une fonction pour faire une moyenne en enlevant les valeurs manquantes
#'
#' @param x un vecteur numerique
#'
#' @return la fonction renvoie la moyenne d'un vecteur
#' @import magrittr
#' @importFrom stats na.omit
#' @examples
#' moyenne(c(4,5))
#' @export
moyenne <- function(x){
  x <- x %>% na.omit()
  sum(x)/length(x)
}
#usethis::use_package("stats")
#usethis::use_package("magrittr")


#' Fonction qui ne sert à rien
#'
#' @return NULL
#' @export
#'
test <-function(){
  print("Ceci est un test")
}
