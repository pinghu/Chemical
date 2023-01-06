#for python https://www.rdkit.org/docs/GettingStartedInPython.html#fingerprinting-and-molecular-similarity
# for similarity search https://github.com/milvus-io/milvus
#https://my.oschina.net/u/4209276/blog/5191465
#https://www.reddit.com/r/chemistry/comments/q6jn4z/i_built_a_molecular_structure_similarity_search/


#https://www.bioconductor.org/packages/devel/bioc/vignettes/ChemmineR/inst/doc/ChemmineR.html

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("ChemmineR")

library("ChemmineR") # Loads the package
library(help="ChemmineR") # Lists all functions and classes 
vignette("ChemmineR") # Opens this PDF manual from R 

data(sdfsample) 
sdfset <- sdfsample
sdfset # Returns summary of SDFset 

sdfset[1:4] 
sdfset[[1]]
view(sdfset[1:4]) # Returns summarized content of many SDFs, not printed here 
as(sdfset[1:4], "list") # Returns complete content of many SDFs, not printed here 
sdfset <- read.SDFset("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/Samples/sdfsample.sdf") 
header(sdfset[1:4]) # Not printed here
header(sdfset[[1]])
atomblock(sdfset[1:4])
atomblock(sdfset[[1]])[1:4,] 
bondblock(sdfset[1:4])
bondblock(sdfset[[1]])[1:4,] 
datablock(sdfset[1:4])
datablock(sdfset[[1]])[1:4]
cid(sdfset)[1:4]
sdfid(sdfset)[1:4]
unique_ids <- makeUnique(sdfid(sdfset))
cid(sdfset) <- unique_ids 
blockmatrix <- datablock2ma(datablocklist=datablock(sdfset)) # Converts data block to matrix 
numchar <- splitNumChar(blockmatrix=blockmatrix) # Splits to numeric and character matrix 
numchar[[1]][1:2,1:2] # Slice of numeric matrix 
numchar[[2]][1:2,10:11] # Slice of character matrix 
propma <- data.frame(MF=MF(sdfset), MW=MW(sdfset), atomcountMA(sdfset))
propma[1:4, ]
datablock(sdfset) <- propma 
datablock(sdfset[1])
grepSDFset("650001", sdfset, field="datablock", mode="subset") # Returns summary view of matches. Not printed here.
grepSDFset("650001", sdfset, field="datablock", mode="index") 
write.SDF(sdfset[1:4], file="sub.sdf", sig=TRUE) 
#adjust plot margins
par(mar = c(1, 1, 1, 1))
plot(sdfset[1:4], print=FALSE) # Plots structures to R graphics device 
sdf.visualize(sdfset[1:4]) # Compound viewing in web browser 
#http://chemmine.ucr.edu/myCompounds/addCompounds/
#[1] "http://chemmine.ucr.edu/ChemmineR/showJob/812c6daf-d68a-46c2-9233-869b2b25f414"

#Structure similarity searching and clustering:
apset <- sdf2ap(sdfset) # Generate atom pair descriptor database for searching 
data(apset) # Load sample apset data provided by library. 
cmp.search(apset, apset[1], type=3, cutoff = 0.3, quiet=TRUE) # Search apset database with single compound. 
cmp.cluster(db=apset, cutoff = c(0.65, 0.5), quiet=TRUE)[1:4,] # Binning clustering using variable similarity cutoffs. 

####Error 
#BiocManager::install("ChemmineOB", force=TRUE)
#library(ChemmineOB)
#convertFormatFile("SML","SDF","mycompound.sml","mycompound.sdf")
#sdfset=read.SDFset("mycompound.sdf")
propOB(sdfset[1])
fingerprintOB(sdfset,"FP2")
smartsSearchOB(sdfset[1:5],"[!$(*#*)&!D1]-!@[!$(*#*)&!D1]",uniqueMatches=FALSE)

library(fmcsR)
data(fmcstest) # Loads test sdfset object 
test <- fmcs(fmcstest[1], fmcstest[2], au=2, bu=1) # Searches for MCS with mismatches 
plotMCS(test) # Plots both query compounds with MCS in color 
cmp.similarity(apset[1],
               apset[2])
cid(sdfset) <- sdfid(sdfset)
fpset <- fp2bit(sdfset, type=1) 
fpset <- fp2bit(sdfset, type=2) 
fpset <- fp2bit(sdfset, type=3) 
fpset
fpSim(fpset[1], fpset[2]) 
fpSim(fpset["650065"], fpset, method="Tanimoto", cutoff=0.6, top=6) 
cid(sdfset) <-
  cid(apset) # Assure compound name consistency among objects. 

plot(sdfset[names(cmp.search(apset, apset["650065"], type=2, cutoff=4, quiet=TRUE))], print=FALSE) 
similarities <- cmp.search(apset, apset[1], type=3, cutoff = 10)
sdf.visualize(sdfset[similarities[,1]]) 
#[1] "http://chemmine.ucr.edu/ChemmineR/showJob/a5d51d1f-4545-4019-83df-84e73f5df903"

#Searching PubChem
#Get Compound SDF from PubChem by Id

compounds <- pubchemCidToSDF(c(111,123))
compounds 

#Get Compound SDF from PubChem by InChIkey -- did not work
inchikeys <- c(
  "ZFUYDSOHVJVQNB-FZERPYLPSA-N",
  "KONGRWVLXLWGDV-BYGOPZEFSA-N",
  "AANKDJLVHZQCFG-WLIQWNBFSA-N",
  "SNFRINMTRPQQLE-JQWAAABSSA-N"
)
# You should only have 2 SDF returned, 2 other not found
inchikey_query <- pubchemInchikey2sdf(inchikeys)
inchikey_query$sdf_set

# successful queries
inchikey_query_index <- inchikey_query$sdf_index[inchikey_query$sdf_index != 0]

# get CID of these queries
inchikey_query_cid <- cid(inchikey_query$sdf_set[inchikey_query_index])
names(inchikey_query_cid) <- names(inchikey_query_index)
inchikey_query_cid

inchis <- c(
  "InChI=1S/C15H26O/c1-9(2)11-6-5-10(3)15-8-7-14(4,16)13(15)12(11)15/h9-13,16H,5-8H2,1-4H3/t10-,11+,12-,13+,14+,15-/m1/s1", 
  "InChI=1S/C3H8/c1-3-2/h3H2,1-2H3", 
  "InChI=1S/C15H20Br2O2/c1-2-12(17)13-7-3-4-8-14-15(19-13)10-11(18-14)6-5-9-16/h3-4,6,9,11-15H,2,7-8,10H2,1H3/t5-,11-,12+,13+,14-,15-/m1/s1",
  "InChI=abc"
)
pubchemInchi2cid(inchis)

#Search a SMILES Query in PubChem ### did not work because the URL link is not good anymore
#should use https://pubchem.ncbi.nlm.nih.gov/#query=CC(%3DO)OC1%3DCC%3DCC%3DC1C(%3DO)O

compounds <- searchString("CC(=O)OC1=CC=CC=C1C(=O)O") 
compounds
#
data(sdfsample); 
sdfset <- sdfsample[1] 
compounds <- searchSim(sdfset) 
compounds 
listCMTools()
toolDetails("ChEMBL Fingerprint Search")

job1 <- launchCMTool("pubchemID2SDF", 2244)
status(job1)
result1 <- result(job1)

job2 <- launchCMTool('ChEMBL Fingerprint Search', result1, 'Similarity Cutoff'=0.95, 'Max Compounds Returned'=200)
result2 <- result(job2)
job3 <- launchCMTool("pubchemID2SDF", result2)
result3 <- result(job3)

job4 <- launchCMTool("OpenBabel Descriptors", result3)
result4 <- result(job4)
result4[1:10,] # show first 10 lines of result

sessionInfo()


#https://cran.r-project.org/web/packages/rcdk/vignettes/using-rcdk.html
#this is simple not very detail usage, good for create new molecules. 