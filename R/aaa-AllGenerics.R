#' @rdname coseq
#' @export
setGeneric("coseq", function(object, ...) standardGeneric("coseq"))

#' @rdname coseqHelpers
#' @export
setGeneric("coseqFullResults", function(object, ...) standardGeneric("coseqFullResults"))

#' @rdname compareARI
#' @export
setGeneric("compareARI", function(object, ...) standardGeneric("compareARI"))

#' @rdname coseqHelpers
#' @export
setGeneric("clusters", function(object, ...) standardGeneric("clusters"))

#' @rdname coseqHelpers
#' @export
setGeneric("likelihood", function(object,...) standardGeneric("likelihood"))

#' @rdname coseqHelpers
#' @export
setGeneric("nbCluster", function(object,...) standardGeneric("nbCluster"))

#' @rdname coseqHelpers
#' @export
setGeneric("proba", function(object,...) standardGeneric("proba"))

#' @rdname coseqHelpers
#' @export
setGeneric("ICL", function(object,...) standardGeneric("ICL"))

#' @rdname coseqHelpers
#' @export
setGeneric("profiles", function(object,...) standardGeneric("profiles"))

#' @rdname coseqHelpers
#' @export
setGeneric("tcounts", function(object,...) standardGeneric("tcounts"))

#' @rdname coseqHelpers
#' @export
setGeneric("transformationType", function(object,...) standardGeneric("transformationType"))

#' @rdname coseqHelpers
#' @export
setGeneric("model", function(object,...) standardGeneric("model"))

#' @rdname coseqHelpers
setGeneric("DDSEextract", function(object,...) standardGeneric("DDSEextract"))

#' @rdname coseqHelpers
setGeneric("Djumpextract", function(object,...) standardGeneric("Djumpextract"))