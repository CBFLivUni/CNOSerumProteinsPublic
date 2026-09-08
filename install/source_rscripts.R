
SourceExternalScripts <- function(dir, pattern, ...) {
  
  require(here)
  
  # Get list of files to source
  files2source <- here::here(dir,
                         list.files(pattern=pattern, path=dir, recursive=TRUE, ...))
  print(files2source)
  
  # Source extra scripts
  for (script in files2source) { source(script) }
  
}
