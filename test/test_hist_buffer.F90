program test_hist_buffer
   use ISO_FORTRAN_ENV,  only: REAL64, REAL32, INT32, INT64
   use hist_msg_handler, only: hist_log_messages, hist_add_error, error_length
   use hist_msg_handler, only: ERROR, VERBOSE
   use hist_buffer,      only: hist_buffer_t
   use hist_field,       only: hist_field_info_t, hist_get_field
   use hist_api,         only: hist_new_field, hist_new_buffer
   use hist_api,         only: hist_field_accumulate, hist_buffer_norm_value

   implicit none

   integer                           :: index
   integer, parameter                :: num_cols = 5
   character(len=error_length)       :: errmsg
   class(hist_field_info_t), pointer :: my_fields => NULL()
   class(hist_field_info_t), pointer :: fld_ptr => NULL()
   class(hist_buffer_t),     pointer :: buffer => NULL()
   class(hist_buffer_t),     pointer :: buff_ptr => NULL()
   type(hist_log_messages)           :: errors
   real(REAL32)                      :: field32(5), test_val32
   real(REAL64)                      :: field64(5), test_val64

   errors%message_level = VERBOSE
   ! Create some fields and some buffers
   my_fields => hist_new_field('U', 'eastward_wind', 'Meridional Wind',       &
        'm s-1', 'real', errors=errors)
   my_fields%next => hist_new_field('T', 'temperature', 'Temperature',        &
        'K', 'real', errors=errors)
   call hist_new_buffer(my_fields, (/ num_cols /), REAL32, 1, 'lst', 1,       &
        buffer, errors=errors)
   call hist_new_buffer(my_fields, (/ num_cols /), REAL32, 1, 'min', 1,       &
        buffer, errors=errors)
   call hist_new_buffer(my_fields, (/ num_cols /), REAL32, 1, 'max', 1,       &
        buffer, errors=errors)
   call hist_new_buffer(my_fields, (/ num_cols /), REAL32, 1, 'avg', 1,       &
        buffer, errors=errors)
   fld_ptr => my_fields%next
   call hist_new_buffer(fld_ptr, (/ num_cols /), REAL64, 1, 'lst', 1,         &
        buffer, errors=errors)
   ! Put some data into the buffers
   do index = 1, num_cols
      field32(index) = 2.0_real32 * real(index, real32)
      field64(index) = (2.0_real64 * real(index, real64)) - 1.0_real64
   end do
   ! Do some accumulation
   do index = 1, 2
      if (errors%num_errors() == 0) then
         call hist_field_accumulate(my_fields, field32, 1, logger=errors)
         if (errors%num_errors() > 0) then
            call  errors%add_stack_frame(ERROR, __FILE__, __LINE__ - 2,       &
                 subname='my_fields32 accumulate')
         else
            call hist_field_accumulate(fld_ptr, field64, 1, logger=errors)
            if (errors%num_errors() > 0) then
               call  errors%add_stack_frame(ERROR, __FILE__, __LINE__ - 2,    &
                    subname='my_fields64 accumulate')
            end if
         end if
      end if
      field32(:) = field32(:) * 2.0_real32
      field64(:) = field64(:) * 3.0_real64
   end do
   ! Read out and check answers
   buff_ptr => my_fields%buffers
   do
      if (associated(buff_ptr) .and. (errors%num_errors() == 0)) then
         field32(:) = 0.0_real32
         call hist_buffer_norm_value(buff_ptr, field32, logger=errors)
         if (errors%num_errors() > 0) then
            call  errors%add_stack_frame(ERROR, __FILE__, __LINE__ - 3,       &
                 subname='my_fields accumulate')
            exit
         else
            ! Check answers
            select case(buff_ptr%buffer_type())
            case ('buff_1dreal32_lst', 'buff_1dreal32_max')
               do index = 1, num_cols
                  test_val32 = 4.0_real32 * real(index, real32)
                  if (field32(index) /= test_val32) then
                     write(errmsg, '(a,i0,2(a,f8.3))') 'field(', index,       &
                          ') = ', field32(index), ' /= ', test_val32
                     call hist_add_error(buff_ptr%buffer_type()//' check',    &
                          trim(errmsg), errors=errors)
                  end if
               end do
            case ('buff_1dreal32_min')
               do index = 1, num_cols
                  test_val32 = 2.0_real32 * real(index, real32)
                  if (field32(index) /= test_val32) then
                     write(errmsg, '(a,i0,2(a,f8.3))') 'field(', index,       &
                          ') = ', field32(index), ' /= ', test_val32
                     call hist_add_error(buff_ptr%buffer_type()//' check',    &
                          trim(errmsg), errors=errors)
                  end if
               end do
            case ('buff_1dreal32_avg')
               do index = 1, num_cols
                  test_val32 = 3.0_real32 * real(index, real32)
                  if (field32(index) /= test_val32) then
                     write(errmsg, '(a,i0,2(a,f8.2))') 'field(', index,       &
                          ') = ', field32(index), ' /= ', test_val32
                     call hist_add_error(buff_ptr%buffer_type()//' check',    &
                          trim(errmsg), errors=errors)
                  end if
               end do
            case default
               call hist_add_error('answer check, unknown buffer type',       &
                    buff_ptr%buffer_type(), errors=errors)
            end select
            buff_ptr => buff_ptr%next
         end if
      else
         exit
      end if
   end do
   field64(:) = 0.0_real64
   buff_ptr => fld_ptr%buffers
   call hist_buffer_norm_value(buff_ptr, field64, logger=errors)
   if (errors%num_errors() > 0) then
      call  errors%add_stack_frame(ERROR, __FILE__, __LINE__ - 2,             &
           subname='my_fields accumulate')
   else
      ! Check answers
      if (buff_ptr%buffer_type() == 'buff_1dreal64_lst') then
         do index = 1, num_cols
            test_val64 = 6.0_real64 * real(index, real64) - 3.0_real64
            if (field64(index) /= test_val64) then
               write(errmsg, '(a,i0,2(a,f8.3))') 'field(', index,             &
                    ') = ', field64(index), ' /= ', test_val64
               call hist_add_error('hist_buff_1dreal64_lst_t check',          &
                    trim(errmsg), errors=errors)
            end if
         end do
      else
         call hist_add_error('answer check, unknown buffer type',             &
              buff_ptr%buffer_type(), errors=errors)
      end if
   end if

   if (errors%num_errors() > 0) then
      write(6, '(a,i0,a)') 'FAIL, ', errors%num_errors(), ' errors found'
      call errors%output(6)
      flush(6)
      STOP 1
   else
      STOP 0
   end if

end program test_hist_buffer
