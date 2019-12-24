

# Nicholas Ducharme-Barth
# 18/12/2018
# create coast polygon for pacifc basin and WCPFC convention area

# load packages
	library(raster)
	library(sp)
	library(rgeos)
	library(rgdal)

# load the coastline data
	coast = readOGR("T:\\NicholasDB\\gshhg-shp-2\\GSHHS_shp\\l","GSHHS_l_L1")
	coast = spTransform(coast, CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

	# rotate to center on the Pacific
	coast = maptools::nowrapSpatialPolygons(maptools::nowrapRecenter(coast))

# load the WCPFC convention area	
	load("T:\\NicholasDB\\WCPFC.skj.shp.RData")

# clip coast to pacific basin area
# define bounding box
	x1 = c(260,260,90,90)
	y1 = c(-70,70,70,-70)
	c1 = cbind(x1, y1)
	r1 = rbind(c1, c1[1, ])  # join
	pb.bb = SpatialPolygons(list(Polygons(list(Polygon(r1)), ID = "bb")),proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
	pb.coast = gIntersection(coast,pb.bb)

# clip coast to wcpo area
# define bounding box
	x1 = c(210,210,110,110)
	y1 = c(-50,50,50,-50)
	c1 = cbind(x1, y1)
	r1 = rbind(c1, c1[1, ])  # join
	wcpo.bb = SpatialPolygons(list(Polygons(list(Polygon(r1)), ID = "bb")),proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
	wcpo.coast = gIntersection(coast,wcpo.bb)

# clip coast to wcpfc area
	wcpfc.skj.coast = gIntersection(coast,WCPFC.skj.shp)

# save
	save(pb.coast,wcpo.coast,wcpfc.skj.coast,file="T:\\NicholasDB\\coast.shp.RData")