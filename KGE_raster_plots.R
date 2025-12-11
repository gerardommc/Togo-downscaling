
library(terra)

r <- rast("KGE_all-downscales.tif")

rdf <- as.data.frame(r, xy = T)

tasmax <- data.frame(rdf[, c("x", "y")], rdf[,grep("tasmax", names(rdf))])
tasmin <- data.frame(rdf[, c("x", "y")], rdf[,grep("tasmin", names(rdf))])
tas <- data.frame(rdf[, c("x", "y")], rdf[,grep("tas_", names(rdf))])
hurs <- data.frame(rdf[, c("x", "y")], rdf[,grep("hurs", names(rdf))])
prec <- data.frame(rdf[, c("x", "y")], rdf[,grep("pr_", names(rdf))])

names(tasmax)[3:7] <- substr(names(tasmax)[3:7], 8, nchar(names(tasmax)[3:7]))
names(tasmin)[3:7] <- substr(names(tasmin)[3:7], 8, nchar(names(tasmin)[3:7]))
names(tas)[3:7] <- substr(names(tas)[3:7], 5, nchar(names(tas)[3:7]))
names(prec)[3:7] <- substr(names(prec)[3:7], 4, nchar(names(prec)[3:7]))
names(hurs)[3:7] <- substr(names(hurs)[3:7], 6, nchar(names(hurs)[3:7]))

tasmax.l <- reshape2::melt(tasmax, id.vars = c("x", "y"))
tasmin.l <- reshape2::melt(tasmin, id.vars = c("x", "y"))
tas.l <- reshape2::melt(tas, id.vars = c("x", "y"))
prec.l <- reshape2::melt(prec, id.vars = c("x", "y"))
hurs.l <- reshape2::melt(hurs, id.vars = c("x", "y"))


library(ggplot2)
library(viridis)

dir.create("KGE-plots")

pdf("KGE-plots/KGE-validation-maps.pdf", width = 6, height = 9)
ggplot(tasmax.l) + geom_tile(aes(x = x, y = y, fill = value)) +
  facet_wrap(~ variable) + coord_equal() +
  scale_fill_viridis(option = "A") +
  labs(fill = "KGE", title = "Maximum temperature") +
  theme_light()

ggplot(tasmin.l) + geom_tile(aes(x = x, y = y, fill = value)) +
  facet_wrap(~ variable) + coord_equal() +
  scale_fill_viridis(option = "A") +
  labs(fill = "KGE", title = "Minimum temperature") +
  theme_light()

ggplot(tas.l) + geom_tile(aes(x = x, y = y, fill = value)) +
  facet_wrap(~ variable) + coord_equal() +
  scale_fill_viridis(option = "A") +
  labs(fill = "KGE", title = "Mean temperature") +
  theme_light()

ggplot(prec.l) + geom_tile(aes(x = x, y = y, fill = value)) +
  facet_wrap(~ variable) + coord_equal() +
  scale_fill_viridis(option = "A") +
  labs(fill = "KGE", title = "Precipitation") +
  theme_light()

ggplot(hurs.l) + geom_tile(aes(x = x, y = y, fill = value)) +
  facet_wrap(~ variable) + coord_equal() +
  scale_fill_viridis(option = "A") +
  labs(fill = "KGE", title = "Relative humidity") +
  theme_light()
dev.off()
