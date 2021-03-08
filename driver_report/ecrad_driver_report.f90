! Author:  Mikhail Zhigun
! Email:   mikhail.zhigun@meteoswiss.ch

module ecrad_driver_report
    use ISO_C_BINDING
interface

 function c_create_report(file_path, file_path_len) result(res)&
        bind(C,name="createReport")
        import C_INTPTR_T, C_CHAR, C_INT
      character(len=1,kind=C_CHAR), dimension(*), intent(in) :: file_path
      integer(kind=C_INT), value, intent(in) :: file_path_len
      integer(C_INTPTR_T) :: res
 end function c_create_report
        
 subroutine delete_report(report) &
        bind(C,name="deleteReport")
        import C_INTPTR_T
      integer(C_INTPTR_T),value, intent(in) :: report
 end subroutine delete_report
		
 subroutine save_report(report) &
        bind(C,name="saveReport")
        import :: C_INTPTR_T, C_CHAR, C_INT
      integer(C_INTPTR_T),value, intent(in) :: report
 end subroutine save_report
		
 subroutine c_start_report_step(report, step_name, name_len) &
        bind(C,name="startReportStep")
        import :: C_INTPTR_T, C_CHAR, C_INT
      integer(C_INTPTR_T),value, intent(in) :: report
      character(len=1,kind=C_char), dimension(*), intent(in) :: step_name
      integer(kind=C_int), value, intent(in) :: name_len
 end subroutine c_start_report_step
		
 subroutine finish_current_report_step(report, step_result) &
        bind(C,name="finishCurrentReportStep")
        import :: C_INTPTR_T, C_INT
      integer(C_INTPTR_T),value, intent(in) :: report
      integer(kind=C_INT), value, intent(in) :: step_result
 end subroutine finish_current_report_step

 subroutine c_fail_current_report_step(report, error_msg, error_msg_len) &
        bind(C,name="failCurrentReportStep")
        import :: C_INTPTR_T, C_INT, C_CHAR
      integer(C_INTPTR_T),value, intent(in) :: report
      character(len=1,kind=C_char), dimension(*), intent(in) :: error_msg
      integer(kind=C_INT), value, intent(in) :: error_msg_len
 end subroutine c_fail_current_report_step
end interface
public
    integer(kind=C_INT), parameter :: C_TRUE  = 1
    integer(kind=C_INT), parameter :: C_FALSE  = 0

	contains

 function create_report(file_path) result(report)
      integer(C_INTPTR_T) :: report
      character(*), intent(in) :: file_path
    report = c_create_report(file_path, len(file_path))
 end function create_report
	
 subroutine start_report_step(report, step_name)
      integer(C_INTPTR_T),value, intent(in) :: report
      character(*), intent(in) :: step_name
    call c_start_report_step(report, step_name, len(step_name))
 end subroutine start_report_step

 subroutine fail_current_report_step(report, error_msg)
      integer(C_INTPTR_T),value, intent(in) :: report
      character(*), intent(in) :: error_msg
    call c_fail_current_report_step(report, error_msg, len(error_msg))
 end subroutine fail_current_report_step

end module ecrad_driver_report
