########################################################
# cmake file for building Marlin example Package
# @author Gerald Grenier, IPNL (modifications)
# @author Jan Engels, Desy IT (original)
########################################################
### LIBRARY #################################################################
# include directories
INCLUDE_DIRECTORIES( ./include )

# add library
AUX_SOURCE_DIRECTORY( ./src library_sources )
ADD_SHARED_LIBRARY( ${PROJECT_NAME} ${library_sources} )
INSTALL_SHARED_LIBRARY( ${PROJECT_NAME} DESTINATION lib )


