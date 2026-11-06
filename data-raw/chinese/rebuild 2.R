# Run from package root. Source: Stellarium v0.22.2 skycultures/chinese.
# These input data and their adaptation are CC BY-SA 4.0; see inst/COPYRIGHTS.
source('R/importConstellations.R')
x <- importConstellations('data-raw/chinese/constellationship.fab',set='chinese')
n <- readLines('data-raw/chinese/constellation_names.eng.fab',encoding='UTF-8')
n <- n[grepl('^[0-9]+[[:space:]]',n)]
ids <- sub('^([0-9]+).*','\\1',n)
pinyin <- sub('^[^\"]*\"([^\"]*)\".*','\\1',n)
english <- sub('.*_\\(\"(.*)\"\\).*','\\1',n)
x$name3 <- pinyin[match(x$name1,ids)]
x$name2 <- x$name4 <- paste0('CHN',x$name1)
x$name1 <- english[match(x$name1,ids)]
order <- c(9,11,4,6,26,24,8,5,16,17,28,23,18,2,12,14,25,15,1,34,3,10,7,13,27,30,29,31)
m <- x[match(sprintf('CHN%03d',order),x$name2),];m$set<-'chinese_mansions'
load('data/constellations.rda')
constellations <- rbind(subset(constellations,set=='modern'),x,m)
rownames(constellations) <- NULL
stopifnot(nrow(x)==318L,nrow(m)==28L)
save(constellations,file='data/constellations.rda',compress='xz')
