library(velocyto.R)
library(future)
library(promises)
plan(multiprocess)

show1 <- readRDS("clytia_show.Rds")
velo <- readRDS("clytia_velocity.Rds")

cell_attrs <- readRDS("clytia_cell_attrs.Rds")
emb <- as.matrix(cell_attrs[,c("tSNE1", "tSNE2")])
rownames(emb) <- cell_attrs$cell_names
Sys.time()
s1 <- future({
  show1 <- readRDS("clytia_show.Rds")
  velo <- readRDS("clytia_velocity.Rds")
  show.velocity.on.embedding.cor(emb = emb, 
                               vel = velo, show.grid.flow = TRUE, 
                               arrow.scale = 3, grid.n = 45, cc = show1$cc,
                               cex = 1, xlab = "tSNE1", ylab = "tSNE2")})
show2 <- value(s1)
Sys.time()
Sys.time()
show.velocity.on.embedding.cor(emb = emb, 
                               vel = velo, show.grid.flow = TRUE, 
                               arrow.scale = 3, grid.n = 45, cc = show1$cc,
                               cex = 1, xlab = "tSNE1", ylab = "tSNE2")
Sys.time()