add_library(module OBJECT 
				abs.cxx add.cxx billow.cxx blend.cxx cache.cxx checkerboard.cxx clamp.cxx const.cxx
				curve.cxx cylinders.cxx displace.cxx exponent.cxx invert.cxx max.cxx min.cxx
				modulebase.cxx multiply.cxx perlin.cxx power.cxx ridgedmulti.cxx rotatepoint.cxx
				scalebias.cxx scalepoint.cxx select.cxx spheres.cxx terrace.cxx translatepoint.cxx
				turbulence.cxx voronoi.cxx)

install(FILES 
				abs.h add.h billow.h blend.h cache.h checkerboard.h clamp.h const.h curve.h
				cylinders.h displace.h exponent.h invert.h max.h min.h module.h modulebase.h
				multiply.h perlin.h power.h ridgedmulti.h rotatepoint.h scalebias.h scalepoint.h 
				select.h spheres.h terrace.h translatepoint.h turbulence.h voronoi.h
				DESTINATION include/noise/module)
