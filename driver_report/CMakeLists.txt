set(driver_report_c_SRC Report.cpp)

PREPEND(driver_report_c_SRC ${CMAKE_CURRENT_SOURCE_DIR} ${driver_report_c_SRC})

target_sources(driver_report_c PRIVATE ${driver_report_c_SRC})

set(driver_report_fortran_wrapper_SRC ecrad_driver_report.f90)

PREPEND(driver_report_fortran_wrapper_SRC ${CMAKE_CURRENT_SOURCE_DIR} ${driver_report_fortran_wrapper_SRC})

target_sources(driver_report_fortran_wrapper PRIVATE ${driver_report_fortran_wrapper_SRC})

add_custom_target(dir-update-driver-report-format)
add_dependencies(update-driver-report-format dir-update-driver-report-format)
add_dependencies(dir-update-driver-report-format create-python-virtual-environment)
set(DRV_REPORT_XSD ${CMAKE_CURRENT_SOURCE_DIR}/driver_report.xsd)
set(DRV_REPORT_PY_SERIALIZATION ${CMAKE_CURRENT_SOURCE_DIR}/driver_report_serialization.py)
add_custom_command(
   TARGET dir-update-driver-report-format
   COMMAND ${PY_VENV_DIR}/bin/python3 ${PY_VENV_DIR}/bin/generateDS.py --silence -f -o ${DRV_REPORT_PY_SERIALIZATION} ${DRV_REPORT_XSD}
   COMMENT "Generating driver report python serialization code"
   )

