# Tests for ArchR 
# change in random number generation in R3.6, this ensures tests will pass under older and newer Rs
# (this was implemented in Seurat's testthat)
library(ArchR)
library(testthat)
library(Matrix)
suppressWarnings(RNGversion(vstr = "3.5.3"))
set.seed(1)

#Context
context("test_workflow")

#Genome
addArchRGenome("hg19test2")

#Arrows
arrows <- file.path(system.file("testdata", package="ArchR"), "PBSmall.arrow")

#Project
proj <- ArchRProject(arrows, outputDirectory = "Test")

#LSI
proj <- addIterativeLSI(proj, dimsToUse = 1:5, varFeatures=1000, iterations = 2, force=TRUE)

features <- getInitialLSIFeatures(proj)
# Test getInitialLSIFeature
test_that("getInitialLSIFeatures returns expected output", {
  
  expect_s4_class(features, "DataFrame")
  expect_true(all(c("seqnames", "start", "idx", "rowSums") %in% names(features)))
  expect_gt(nrow(features), 100)

  # Only check 'end' if the matrix is TileMatrix
  rds_path <- file.path(proj@projectMetadata$outputDirectory, "IterativeLSI", "Save-LSI-Iteration-1.rds")
  iteration <- readRDS(rds_path)
  if (iteration$LSI$useMatrix == "TileMatrix") {
    expect_true("end" %in% names(features))
    expect_equal(features$end, features$start + iteration$LSI$tileSize)
  }
})

test_that("getInitialLSIFeatures returns error if file is nonexistent"){
    expect_error(getInitialLSIFeatures(proj, iterationName = "nonexistent.rds"),
    sprintf("No saved iteration 1 at %s/IterativeLSI/%s. Run `runIterativeLSI()` with saveIterations=TRUE!.", proj@projectMetadata$outputDirectory, rds_path)
  )
    
})

test_that("getInitialLSIFeatures output works in addIterativeLSI", {
  proj2 <- ArchR::addIterativeLSI(
    ArchRProj = proj,
    useMatrix = "TileMatrix",
    name = "IterativeLSI_test_from_features",
    varFeatures = features
  )
  expect_s4_class(proj2, "ArchRProject")
  expect_true("IterativeLSI_test_from_features" %in% names(ArchR::getReducedDims(proj2)))
  assign("proj2", proj2, envir = .GlobalEnv)
})


test_that("IterativeLSI with getInitialLSIFeatures reproduces original LSI embedding", {
  skip_if_not(exists("proj2"), "proj2 not available from previous test")
  # Get original IterativeLSI result
  original_lsi <- ArchR::getReducedDims(proj, reduction = "IterativeLSI")

  # Get new LSI result
  reproduced_lsi <- ArchR::getReducedDims(proj2, reduction = "IterativeLSI_from_features")

  # Ensure the same cells are in both matrices
  common_cells <- intersect(rownames(original_lsi), rownames(reproduced_lsi))
  original_lsi <- original_lsi[common_cells, , drop = FALSE]
  reproduced_lsi <- reproduced_lsi[common_cells, , drop = FALSE]

  # Compare matrices: must be numerically equal within tolerance
  expect_equal(
    reproduced_lsi,
    original_lsi,
    tolerance = 1e-6,
    scale = 1,
    info = "Recomputed LSI should match original"
  )
})

#Clusters
proj <- addClusters(proj, force=TRUE, dimsToUse = 1:5)

#Cell type
proj$CellType <- substr(getCellNames(proj), 9, 9)

#Check Clusters
cM <- confusionMatrix(proj$CellType, proj$Clusters)

test_that("Test Cluster Purity...", {
	expect_equal(all(apply(cM / Matrix::rowSums(cM), 1, max) > 0.9), TRUE)
})

#UMAP
proj <- addUMAP(proj, nNeighbors = 40, dimsToUse = 1:5, minDist=0.1, force=TRUE)

#Plot UMAP
p1 <- plotEmbedding(proj, name = "CellType", size = 3)
p2 <- plotEmbedding(proj, name = "Clusters", size = 3)
plotPDF(p1, p2, name = "UMAP", addDOC = FALSE, ArchRProj = proj)
 
#Group Coverages
proj <- addGroupCoverages(proj)

#Custom Peak Calling
proj <- addReproduciblePeakSet(
    ArchRProj = proj, 
    groupBy = "Clusters", 
    peakMethod = "tiles",
    minCells = 20,
    cutOff = 0.1
)

#Macs2 Peak Calling
pathToMacs2 <- findMacs2()
proj <- addReproduciblePeakSet(
    ArchRProj = proj, 
    groupBy = "Clusters", 
    pathToMacs2 = pathToMacs2,
    minCells = 20,
    cutOff = 0.1
)

#Add Peak Matrix
proj <- addPeakMatrix(proj)

#Motif Annotations
proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbpTest", name = "Motif", force=TRUE)

#Motif Deviations
proj <- addBgdPeaks(proj)
proj <- addDeviationsMatrix(
  ArchRProj = proj, 
  peakAnnotation = "Motif",
  force = TRUE
)

#Check Matrices
se <- getMatrixFromProject(proj, "MotifMatrix")

pval_CEBP <- t.test(
	assays(se)$z["CEBPB_1", se$CellType=="T"],
	assays(se)$z["CEBPB_1", se$CellType=="M"]
)$p.value

pval_EOMES <- t.test(
	assays(se)$z["EOMES_6", se$CellType=="T"],
	assays(se)$z["EOMES_6", se$CellType=="M"]
)$p.value

pval_PAX <- t.test(
	assays(se)$z["PAX5_5", se$CellType=="M"],
	assays(se)$z["PAX5_5", se$CellType=="B"]
)$p.value

test_that("Test Motif Deviations...", {
	expect_equal(pval_CEBP < 0.01, TRUE)
	expect_equal(pval_EOMES < 0.1, TRUE)
	expect_equal(pval_PAX < 0.01, TRUE)
})

################################################
# Clear
################################################

files <- list.files()
files <- files[!grepl("\\.R", files)]
for(i in seq_along(files)){
	if(dir.exists(files[i])){
		unlink(files[i], recursive=TRUE)
	}else if(file.exists(files[i])){
		file.remove(files[i])
	}
}

