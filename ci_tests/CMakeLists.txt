set(CI_TESTS_RUNTIME_DIR ${CMAKE_BINARY_DIR}/ci_tests)
set(CI_TESTS_RUNTIME_DIR ${CI_TESTS_RUNTIME_DIR} PARENT_SCOPE)
set(ECRAD_CFG_DIR ${CMAKE_SOURCE_DIR}/data)
set(ECRAD_CI_TESTS_DATA_DIR ${CMAKE_SOURCE_DIR}/ci_tests/data)

file(MAKE_DIRECTORY ${CI_TESTS_RUNTIME_DIR})
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config.py.template
               ${CI_TESTS_RUNTIME_DIR}/config.py)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/tests_runner.py ${CI_TESTS_RUNTIME_DIR}/tests_runner.py COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/combine_block_len_test_results.py ${CI_TESTS_RUNTIME_DIR}/combine_block_len_test_results.py COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/plot_block_len_test_results.py ${CI_TESTS_RUNTIME_DIR}/plot_block_len_test_results.py COPYONLY)

add_custom_target(update-cpu-reg-test-format)
add_dependencies(update-test-report-formats update-cpu-reg-test-format)
add_dependencies(update-cpu-reg-test-format create-python-virtual-environment)
set(CPU_REG_REPORT_TEST_XSD ${CMAKE_CURRENT_SOURCE_DIR}/cpu_regression_test_report.xsd)
set(CPU_REG_REPORT_TEST_PY_SERIALIZATION ${CMAKE_CURRENT_SOURCE_DIR}/cpu_reg_test_report_serialization.py)
add_custom_command(
   TARGET update-cpu-reg-test-format
   COMMAND ${PY_VENV_DIR}/bin/python3 ${PY_VENV_DIR}/bin/generateDS.py --silence -f -o ${CPU_REG_REPORT_TEST_PY_SERIALIZATION} ${CPU_REG_REPORT_TEST_XSD}
   COMMENT "Generating cpu regression report python serialization code"
   )

