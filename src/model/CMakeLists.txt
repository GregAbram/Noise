add_library(model OBJECT cylinder.cxx line.cxx plane.cxx sphere.cxx)

install(FILES cylinder.h line.h model.h plane.h sphere.h DESTINATION include/noise/model)
