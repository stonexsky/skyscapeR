library(skyscapeR)

# Use an empty sky to isolate the atmosphere and sample the final display pixels.
sample.sky <- function(solar.altitude) {
  sky <- list(light=list(sun=c(180,solar.altitude), k=.2, calibrated.twilight=TRUE))
  field <- skyscapeR:::.shared_radiance(sky,c(0,360,0,90),4,3,colour=TRUE)
  skyscapeR:::.shared_pixels(field,sky)[1,1]
}

# These twilight and sunrise pixels are from 2.0.9009.
stopifnot(sample.sky(-6) == "#344557", sample.sky(0) == "#526D7E")
for (altitude in c(15,45,75)) {
  rgb <- grDevices::col2rgb(sample.sky(altitude))
  stopifnot(rgb[3] > rgb[2], rgb[2] > rgb[1], rgb[3]-rgb[1] > 65)
}
