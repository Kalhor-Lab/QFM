# library(TreeTools)
# RootedSPRSwapMod <- function (parent, child, nEdge = length(parent), nNode = nEdge/2L,
#           edgeToBreak = NULL, mergeEdge = NULL)
# {
#         if (nEdge < 5)
#                 return(SPRWarning(parent, child, "Too few tips to rearrange."))
#         rootNode <- parent[1]
#         rootEdges <- parent == rootNode
#         breakable <- !logical(nEdge) & !rootEdges
#         if (!is.null(edgeToBreak) && edgeToBreak == -1) {
#                 notDuplicateRoot <- NonDuplicateRoot(parent, child, nEdge)
#                 return(unique(unlist(lapply(which(breakable), AllSPR,
#                                             parent = parent, child = child, nEdge = nEdge, notDuplicateRoot = notDuplicateRoot),
#                                      recursive = FALSE)))
#         }
#         rightSide <- DescendantEdges(1, parent, child, nEdge)
#         leftSide <- !rightSide
#         nEdgeRight <- which(rootEdges)[2] - 1
#         nEdgeLeft <- nEdge - nEdgeRight
#         if (nEdgeRight < 4) {
#                 if (nEdgeLeft < 4)
#                         return(SPRWarning(parent, child, "No rearrangement possible with this root position."))
#                 breakable <- breakable & !rightSide
#                 rightHalfOfLeftSide <- DescendantEdges(nEdgeRight + 2L,
#                                                        parent, child, nEdge)
#                 leftHalfOfLeftSide <- leftSide & !rightHalfOfLeftSide &
#                         !rootEdges
#                 if (sum(rightHalfOfLeftSide) == 1)
#                         breakable[nEdgeRight + 3] <- FALSE
#                 if (sum(leftHalfOfLeftSide) == 1)
#                         breakable[nEdgeRight + 2] <- FALSE
#         }
#         else {
#                 if (nEdgeLeft < 4) {
#                         breakable <- breakable & rightSide
#                 }
#                 else {
#                         rightHalfOfLeftSide <- DescendantEdges(nEdgeRight +
#                                                                        2L, parent, child, nEdge)
#                         leftHalfOfLeftSide <- leftSide & !rightHalfOfLeftSide &
#                                 !rootEdges
#                         if (sum(rightHalfOfLeftSide) == 1)
#                                 breakable[nEdgeRight + 3] <- FALSE
#                         if (sum(leftHalfOfLeftSide) == 1)
#                                 breakable[nEdgeRight + 2] <- FALSE
#                 }
#                 rightHalfOfRightSide <- DescendantEdges(2L, parent, child,
#                                                         nEdge)
#                 leftHalfOfRightSide <- rightSide & !rightHalfOfRightSide &
#                         !rootEdges
#                 if (sum(rightHalfOfRightSide) == 1)
#                         breakable[3] <- FALSE
#                 if (sum(leftHalfOfRightSide) == 1)
#                         breakable[2] <- FALSE
#         }
#         if (is.null(edgeToBreak)) {
#                 edgeToBreak <- SampleOne(which(breakable))
#         }
#         else {
#                 if (!breakable[edgeToBreak])
#                         return(SPRWarning(parent, child, paste("Nowhere to regraft if pruning on edge",
#                                                                edgeToBreak)))
#                 if (edgeToBreak > nEdge)
#                         return(SPRWarning(parent, child, "edgeToBreak > nEdge"))
#                 if (edgeToBreak < 1)
#                         return(SPRWarning(parent, child, "edgeToBreak < 1"))
#         }
#         brokenEdge <- seq_along(parent) == edgeToBreak
#         brokenEdge.parentNode <- parent[edgeToBreak]
#         brokenEdge.childNode <- child[edgeToBreak]
#         edgesCutAdrift <- DescendantEdges(edgeToBreak, parent, child,
#                                           nEdge)
#         edgesOnAdriftSegment <- edgesCutAdrift | brokenEdge
#         brokenEdgeParent <- child == brokenEdge.parentNode
#         brokenEdgeSister <- parent == brokenEdge.parentNode & !brokenEdge
#         brokenEdgeDaughters <- parent == brokenEdge.childNode
#         nearBrokenEdge <- brokenEdge | brokenEdgeSister | brokenEdgeParent |
#                 brokenEdgeDaughters
#         if (!is.null(mergeEdge)) {
#                 if (mergeEdge > nEdge)
#                         return(SPRWarning(parent, child, "mergeEdge value > number of edges"))
#                 if (length(mergeEdge) != 1)
#                         return(SPRWarning(parent, child, paste0("mergeEdge value ",
#                                                                 paste(mergeEdge, collapse = "|"), " invalid; must be NULL or a vector of length 1\n")))
#                 if (nearBrokenEdge[mergeEdge])
#                         return(SPRWarning(parent, child, "Selected mergeEdge will not change tree topology."))
#                 if (DescendantEdges(edgeToBreak, parent, child, nEdge)[mergeEdge])
#                         stop("mergeEdge is within pruned subtree")
#         }
#         else {
#                 edgesOnThisSide <- if (rightSide[edgeToBreak])
#                         rightSide
#                 else leftSide
#                 mergeEdge <- which(edgesOnThisSide & !nearBrokenEdge &
#                                            !edgesOnAdriftSegment)
#                 nCandidates <- length(mergeEdge)
#                 if (nCandidates > 1)
#                         mergeEdge <- SampleOne(mergeEdge, len = nCandidates)
#         }
#         parent[brokenEdgeSister] <- parent[brokenEdgeParent]
#         parent[brokenEdgeParent] <- parent[mergeEdge]
#         parent[mergeEdge] <- brokenEdge.parentNode
#         return(list(parent, child))
# }
