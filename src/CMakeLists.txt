project(Noise)

add_subdirectory(model)
add_subdirectory(module)

add_library(noise 
						latlon.cxx noisegen.cxx
						$<TARGET_OBJECTS:model>
						$<TARGET_OBJECTS:module>)

install(TARGETS noise DESTINATION lib)
install(FILES basictypes.h exception.h interp.h latlon.h mathconsts.h misc.h
				noise.h noisegen.h vectortable.h DESTINATION include/noise)
