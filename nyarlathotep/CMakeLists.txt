simple_plugin(Nyarlathotep "module"
      larcore_Geometry_Geometry_service
      larcorealg_Geometry
      lardataobj_RecoBase
      larsim_MCCheater_PhotonBackTracker_service
      larsim_MCCheater_BackTracker_service
      nusimdata_SimulationBase
      ${SIMULATIONBASE}
      ${ART_FRAMEWORK_CORE}
      ${ART_FRAMEWORK_PRINCIPAL}
      ${ART_FRAMEWORK_SERVICES_REGISTRY}
      ${ART_FRAMEWORK_SERVICES_OPTIONAL}
      ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
      ${ART_PERSISTENCY_COMMON}
      ${ART_PERSISTENCY_PROVENANCE}
      ${ART_UTILITIES}
      ${MF_MESSAGELOGGER}
      ${MF_UTILITIES}
      ${CETLIB}
      ${ROOT_BASIC_LIB_LIST}
      )

    #simple_plugin(Shoggoth "module"
    #  larcore_Geometry_Geometry_service
    #  larcorealg_Geometry
    #  larsim_MCCheater_PhotonBackTracker_service
    #  larsim_MCCheater_BackTracker_service
    #  ${SIMULATIONBASE}
    #  ${ART_FRAMEWORK_CORE}
    #  ${ART_FRAMEWORK_PRINCIPAL}
    #  ${ART_FRAMEWORK_SERVICES_REGISTRY}
    #  ${ART_FRAMEWORK_SERVICES_OPTIONAL}
    #  ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
    #  ${ART_PERSISTENCY_COMMON}
    #  ${ART_PERSISTENCY_PROVENANCE}
    #  ${ART_UTILITIES}
    #  ${MF_MESSAGELOGGER}
    #  ${MF_UTILITIES}
    #  ${CETLIB}
    #  ${ROOT_BASIC_LIB_LIST}
    #  )

    #simple_plugin(pbtDebug "module"
    #      larcore_Geometry_Geometry_service
    #      larcorealg_Geometry
    #      larsim_MCCheater_PhotonBackTracker_service
    #      larsim_MCCheater_BackTracker_service
    #      nusimdata_SimulationBase
    #      ${SIMULATIONBASE}
    #      ${ART_FRAMEWORK_CORE}
    #      ${ART_FRAMEWORK_PRINCIPAL}
    #      ${ART_FRAMEWORK_SERVICES_REGISTRY}
    #      ${ART_FRAMEWORK_SERVICES_OPTIONAL}
    #      ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
    #      ${ART_PERSISTENCY_COMMON}
    #      ${ART_PERSISTENCY_PROVENANCE}
    #      ${ART_UTILITIES}
    #      ${MF_MESSAGELOGGER}
    #      ${MF_UTILITIES}
    #      ${CETLIB}
    #      ${ROOT_BASIC_LIB_LIST}
    #      )

    #simple_plugin(Shoggoth "module"
    #  larcore_Geometry_Geometry_service
    #  larcorealg_Geometry
    #  larsim_MCCheater_PhotonBackTracker_service
    #  larsim_MCCheater_BackTracker_service
    #  ${SIMULATIONBASE}
    #  ${ART_FRAMEWORK_CORE}
    #  ${ART_FRAMEWORK_PRINCIPAL}
    #  ${ART_FRAMEWORK_SERVICES_REGISTRY}
    #  ${ART_FRAMEWORK_SERVICES_OPTIONAL}
    #  ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
    #  ${ART_PERSISTENCY_COMMON}
    #  ${ART_PERSISTENCY_PROVENANCE}
    #  ${ART_UTILITIES}
    #  ${MF_MESSAGELOGGER}
    #  ${MF_UTILITIES}
    #  ${CETLIB}
    #  ${ROOT_BASIC_LIB_LIST}
    #  )

install_headers()
install_fhicl()
#install_source( EXTRAS ADDITIONAL_NOTES.txt README.txt )
install_source()
install_scripts()
