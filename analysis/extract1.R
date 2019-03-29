# Convenience function to simplify extracting single parameters from stanfit objects
extract1 <- function(object, par)
{
  extract(object, par)[[1]]
}
