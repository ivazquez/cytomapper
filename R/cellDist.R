#' Computes the mean pixel distance between each cell and a specified cell class.
#'
#' For each object (e.g. cell) identified by segmentation, the
#' \code{cellDist} function computes the distance between this object and
#' a specified region (e.g. cell-type, compartment, etc.)
#'
#' @param object
#' @param mask
#' @param img_id
#' @param cell_id
#' @param pattern logical vector indicating which cells are part of the pattern
#' @param expand_by
#'
#' @return A \linkS4class{SingleCellExperiment} object (see details)
#'
#' @section
#'
#' @examples
#'
#'
#' @seealso
#' \code{\link{distmap}}, for underlying computations
#'
#' @author Nils Eling (\email{nils.eling@@dqbm.uzh.ch}),
#'
#' @export
#' @importFrom EBImage computeFeatures.basic distmap
#' @importFrom BiocParallel SerialParam bpmapply
measureObjects <- function(object,
                           mask,
                           img_id,
                           cell_id,
                           pattern_slot,
                           pattern_name,
                           metric = c("min", "mean"),
                           BPPARAM = SerialParam()) {

    # Validity checks
    .valid.mask(mask, img_id = img_id)
    .valid.matchObjects.plotCells(object, mask, img_id)

    # Check for valid pattern

    # Define metric
    metric <- match.arg(metric)

    cur_out <- bplapply(mcols(mask)[,img_id], function(x){

      cur_object <- object[,colData(object)[,img_id] == x]
      cur_mask <- mask[mcols(mask)[,img_id] == x][[1]]
      cur_mask_bin <- cur_mask

      # Build binary mask of pattern - inside pattern
      cur_ind <- colData(cur_object)[,pattern_slot] == pattern_name
      cur_mask_bin[!(cur_mask_bin %in% colData(cur_object)[cur_ind,cell_id])] <- 0
      cur_mask_bin[cur_mask_bin != 0] <- 1

      if (expand_by > 0) {
          kern <- makeBrush(expand_by * 2 + 1, shape = "diamond")
          cur_mask_bin <- dilate(cur_mask_bin, kern)
      }

      cur_distmap <- distmap(cur_mask_bin)

      # Compute features
      if (metric == "mean") {
          cur_dist_1 <- computeFeatures.basic(cur_mask, cur_distmap)
          cur_dist_1 <- cur_dist_1[,"b.mean"]
      } else if (metric == "min") {
         cur_dist_1 <- computeFeatures.basic(cur_mask, cur_distmap, basic.quantiles = 0)
         cur_dist_1 <- cur_dist_1[,"b.q0"]
      }

      # Build binary mask of pattern - outside pattern
      cur_mask_bin <- cur_mask
      cur_ind <- colData(cur_object)[,pattern_slot] == pattern_name
      cur_mask_bin[cur_mask_bin %in% colData(cur_object)[cur_ind,cell_id]] <- 0
      cur_mask_bin[cur_mask_bin > 0] <- 1

      if (expand_by > 0) {
        kern <- makeBrush(expand_by * 2 + 1, shape = "diamond")
        cur_mask_bin <- dilate(cur_mask_bin, kern)
      }

      cur_distmap <- distmap(cur_mask_bin)

      # Compute features
      if (metric == "mean") {
        cur_dist_2 <- computeFeatures.basic(cur_mask, cur_distmap)
        cur_dist_2 <- cur_dist_2[,"b.mean"]
      } else if (metric == "min") {
        cur_dist_2 <- computeFeatures.basic(cur_mask, cur_distmap, basic.quantiles = 0)
        cur_dist_2 <- cur_dist_2[,"b.q0"]
      }

      cur_dist <- cur_dist_2 - cur_dist_1

      names(cur_dist) <- paste0(x, "_", names(cur_dist))

      return(cur_dist)

    }, BPPARAM = BPPARAM)

    cur_out <- do.call(c, cur_out)

    # Correct ordering
    final_out <- rep(NA, ncol(object))
    names(final_out) <- paste(colData(object)[,img_id], colData(object)[,cell_id], sep = "_")

    final_out <-

    cur_out <- cur_out

    colData(object)[,paste0("distance_to_", pattern_slot, "_", pattern_name)] <- cur_out

    return(sce)
}
