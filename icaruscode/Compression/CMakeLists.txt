include_directories(/icarus/app/users/cerati/zstd-1.0.0/lib)

link_directories(/icarus/app/users/cerati/zstd-1.0.0/lib)

simple_plugin( ICARUSTPCCompress module
  ${LARDATAOBJ}
  lardataobj_RawData
  ${ART_UTILITIES}
  ${FHICLCPP}
  ${MF_MESSAGELOGGER}
  ${MF_UTILITIES}
  ${ART_FRAMEWORK_CORE}
  ${ROOT_BASIC_LIB_LIST}
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  ${ART_ROOT_IO_TFILE_SUPPORT}
  ${ROOT_CORE}
  ${ART_ROOT_IO_TFILESERVICE_SERVICE}
  zstd
)
install_headers()
install_fhicl()
install_source()
