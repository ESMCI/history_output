module hist_api

   implicit none
   private

   ! Public API interfaces
   public :: hist_new_field         ! Allocate a hist_field_info_t object
   public :: hist_new_buffer        ! Create a new field buffer
   public :: hist_field_accumulate  ! Accumulate a new field state in all buffs
   public :: hist_field_norm_value  ! Grab the normalized value from field buffer

   ! Interfaces for public interfaces
   interface hist_field_accumulate
      module procedure hist_field_accumulate_1d
      module procedure hist_field_accumulate_2d
   end interface hist_field_accumulate

   interface hist_field_norm_value
      module procedure hist_field_norm_value_1d
      module procedure hist_field_norm_value_2d
   end interface hist_field_norm_value

CONTAINS

   !#######################################################################

   function hist_new_field(diag_name_in, std_name_in, long_name_in, units_in, &
        type_in, decomp_in, dimensions, acc_flag, num_levels, field_shape, fill_value, &
        sampling_seq, flag_xyfill, mixing_ratio, dim_bounds, mdim_sizes, beg_dims, end_dims,  &
        cell_methods, errors) result(new_field)
      use hist_msg_handler, only: hist_have_error, hist_log_messages, ERROR
      use hist_field,       only: hist_field_initialize, hist_field_info_t
      use ISO_FORTRAN_ENV,  only: REAL64

      type(hist_field_info_t), pointer                 :: new_field
      character(len=*),                  intent(in)    :: diag_name_in
      character(len=*),                  intent(in)    :: std_name_in
      character(len=*),                  intent(in)    :: long_name_in
      character(len=*),                  intent(in)    :: units_in
      character(len=*),                  intent(in)    :: type_in
      integer,                           intent(in)    :: decomp_in
      integer,                           intent(in)    :: dimensions(:)
      character(len=*),                  intent(in)    :: acc_flag
      integer,                           intent(in)    :: field_shape(:)
      integer,                           intent(in)    :: num_levels
      real(kind=REAL64),                 intent(in)    :: fill_value
      character(len=*),        optional, intent(in)    :: sampling_seq
      logical,                 optional, intent(in)    :: flag_xyfill
      character(len=*),        optional, intent(in)    :: mixing_ratio
      integer,                 optional, intent(in)    :: dim_bounds(:,:)
      integer,                 optional, intent(in)    :: mdim_sizes(:)
      integer,                 optional, intent(in)    :: beg_dims(:)
      integer,                 optional, intent(in)    :: end_dims(:)
      character(len=*),        optional, intent(in)    :: cell_methods
      type(hist_log_messages), optional, intent(inout) :: errors

      integer                     :: astat
      character(len=128)          :: errmsg
      character(len=*), parameter :: subname = 'hist_new_field'

      if (.not. hist_have_error(errors=errors)) then
         allocate(new_field, stat=astat)
         if ((astat /= 0) .and. present(errors)) then
            call errors%new_error(subname//' Unable to allocate <new_field>')
         end if
      end if
      if (.not. hist_have_error(errors=errors)) then
         call hist_field_initialize(new_field, diag_name_in, std_name_in,     &
              long_name_in, units_in, type_in, decomp_in, dimensions, acc_flag,  &
              num_levels, field_shape, sampling_seq=sampling_seq, flag_xyfill=flag_xyfill,   &
              fill_value=fill_value, mixing_ratio=mixing_ratio, dim_bounds=dim_bounds, &
              mdim_sizes=mdim_sizes, beg_dims=beg_dims, end_dims=end_dims, &
              cell_methods=cell_methods, errmsg=errmsg)
         if (hist_have_error(errors=errors)) then
            call errors%add_stack_frame(ERROR, __FILE__, __LINE__ - 3,        &
                 subname=subname)
         end if
      end if
   end function hist_new_field

   !#######################################################################

   subroutine hist_new_buffer(field, buff_shape, horiz_axis_ind,   &
        accum_type, output_vol, errors, block_ind, block_sizes)
      use hist_msg_handler, only: hist_log_messages, hist_add_error
      use hist_msg_handler, only: hist_add_alloc_error, ERROR
      use hist_hashable,    only: hist_hashable_t
      use hist_buffer,      only: hist_buffer_t, buffer_factory
      use hist_buffer,      only: hist_accum_lst, hist_accum_min, hist_accum_max
      use hist_buffer,      only: hist_accum_avg, hist_accum_var
      use hist_field,       only: hist_field_info_t

      ! Dummy arguments
      class(hist_field_info_t), pointer                 :: field
      integer,                            intent(in)    :: buff_shape(:)
      integer,                            intent(in)    :: horiz_axis_ind
      character(len=*),                   intent(in)    :: accum_type
      integer,                            intent(in)    :: output_vol
      type(hist_log_messages),  optional, intent(inout) :: errors
      integer,                  optional, intent(in)    :: block_ind
      integer,                  optional, intent(in)    :: block_sizes(:)

      ! Local variables
      integer                              :: rank
      integer                              :: line_loc
      integer                              :: accum_val
      integer                              :: shape_idx
      integer,                 allocatable :: beg_dims(:), end_dims(:)
      integer,                 allocatable :: buffer_shape(:)
      character(len=8)                     :: kind_string
      character(len=3)                     :: accum_string
      character(len=16)                    :: bufftype_string
      class(hist_hashable_t),  pointer     :: field_base
      class(hist_buffer_t),    pointer     :: buff_ptr
      class(hist_buffer_t),    pointer     :: buffer
      character(len=:),        allocatable :: type_str
      integer,                 parameter   :: max_rank = 2
      character(len=*),        parameter   :: subname = 'hist_new_buffer'

      ! Initialize output and local variables
      nullify(buffer)
      nullify(field_base)
      nullify(buff_ptr)
      rank = SIZE(buff_shape, 1)
      !! Some sanity checks
      ! We can select on the field's type string but not its kind string
      ! because we do not know the kind value for the kind string
      if (associated(field)) then
         type_str = field%type()
      else
         type_str = 'unknown'
      end if
      ! Check horiz_axis_ind
      if ((horiz_axis_ind < 1) .or. (horiz_axis_ind > rank)) then
         call hist_add_error(subname, 'horiz_axis_ind outside of ',           &
              errstr2='valid range, [1, ', errint2=rank, errstr3=']',         &
              errors=errors)
      end if
      ! Check for (proper) block structured buffer
      if (present(block_ind) .and. present(block_sizes)) then
         if ((block_ind < 1) .or. (block_ind > rank)) then
            call hist_add_error(subname, 'block_ind outside of ',             &
                 errstr2='valid range, [1, ', errint2=rank, errstr3=']',      &
                 errors=errors)
         else if (block_ind == horiz_axis_ind) then
            call hist_add_error(subname, 'block_ind cannot be the same ',     &
                 errstr2='as horiz_axis_ind', errors=errors)
         end if
      else if (present(block_ind)) then
         call hist_add_error(subname,                                         &
              'block_sizes required if block_ind is present', errors=errors)
      else if (present(block_sizes)) then
         call hist_add_error(subname,                                         &
              'block_ind required if block_sizes is present', errors=errors)
      end if ! No else, we just do not have a blocked buffer
      ! Check accumulation type
      select case(trim(accum_type))
      case ('I', 'i', 'lst')
         accum_string = 'lst'
         accum_val = hist_accum_lst
      case ('A', 'a', 'avg')
         accum_string = 'avg'
         accum_val = hist_accum_avg
      case ('M', 'm', 'min')
         accum_string = 'min'
         accum_val = hist_accum_min
      case ('X', 'x', 'max')
         accum_string = 'max'
         accum_val = hist_accum_max
      case ('S', 's', 'var')
         accum_string = 'var'
         accum_val = hist_accum_var
      case default
         call hist_add_error(subname,                                         &
              "Unknown accumulation operator type, '",                        &
              errstr2=trim(accum_type)//"'", errors=errors)
      end select
      ! Only real values handled right now
      if (trim(type_str) /= 'real') then
         call hist_add_error(subname, 'buffer type ', errstr2=type_str,       &
                 errstr3=' unsupported', errors=errors)
      end if
      ! We now know what sort of buffer we need
      ! Based on rank
      select case (rank)
      case (1)
         bufftype_string = 'real_1'
      case (2)
         bufftype_string = 'real_2'
      case default
         ! Over max rank currently handled
         call hist_add_error(subname, 'buffers have a max rank of ',          &
              errint1=max_rank, errors=errors)
      end select
      buffer => buffer_factory(trim(bufftype_string), logger=errors)
      line_loc = __LINE__ - 1
      if (associated(buffer)) then
         field_base => field
         beg_dims = field%beg_dims()
         end_dims = field%end_dims()
         allocate(buffer_shape(size(buff_shape,1)))
         buffer_shape(1) = end_dims(1) - beg_dims(1) + 1
         if (size(buffer_shape) > 1) then
            buffer_shape(2) = buff_shape(2)
         end if
         call buffer%initialize(field_base, output_vol, horiz_axis_ind,       &
              accum_val, buffer_shape, block_sizes, block_ind, &
              logger=errors)
         ! Add this buffer to its field (field should be there if buffer is)
         if (associated(field%buffers)) then
            buff_ptr => field%buffers
            field%buffers => buffer
            buffer%next => buff_ptr
         else
            field%buffers => buffer
         end if
      else
         call hist_add_error(subname, 'buffer ('//trim(bufftype_string),      &
              errstr2=') not created', errors=errors)
         if (present(errors)) then
            call errors%add_stack_frame(ERROR, __FILE__, line_loc,            &
                 subname=subname)
         end if
      end if

   end subroutine hist_new_buffer

   !#######################################################################

   subroutine hist_field_accumulate_1d(field, data, cols_or_block,      &
        cole, logger)
      use hist_msg_handler, only: hist_log_messages, hist_have_error, ERROR
      use hist_msg_handler, only: hist_add_message, VERBOSE
      use hist_field,       only: hist_field_info_t
      use hist_buffer,      only: hist_buff_1d_t, hist_buffer_t
      use ISO_FORTRAN_ENV,  only: REAL64

      ! Dummy arguments
      class(hist_field_info_t), pointer,  intent(inout) :: field
      real(REAL64),                       intent(in)    :: data(:)
      integer,                            intent(in)    :: cols_or_block
      integer,                  optional, intent(in)    :: cole
      type(hist_log_messages),  optional, intent(inout) :: logger
      ! Local variables
      class(hist_buffer_t), pointer     :: buff_ptr
      class(hist_buff_1d_t), pointer    :: buff
      character(len=*), parameter :: subname = 'hist_field_accumulate_1d'

      if (associated(field)) then
         buff_ptr => field%buffers
         do
            if (associated(buff_ptr) .and.                                    &
                 (.not. hist_have_error(errors=logger))) then
               select type(buff_ptr)
               class is (hist_buff_1d_t)
                  buff => buff_ptr
                  call buff%accumulate(data, cols_or_block, field%flag_xyfill(), &
                          field%fill_value(), cole, logger)
                  if (hist_have_error(errors=logger)) then
                     call  logger%add_stack_frame(ERROR, __FILE__, __LINE__ - 3, &
                          subname=subname)
                     exit
                  else
                     call hist_add_message(subname, VERBOSE,                     &
                          "Accumulated data for",                                &
                          msgstr2=trim(field%diag_name())//", Buffer type, ",    &
                          msgstr3=trim(buff%buffer_type()),                  &
                          logger=logger)
                  end if
               end select
               buff_ptr => buff_ptr%next
            else
               exit
            end if
         end do
      end if ! No else, it is legit to pass in a null pointer

   end subroutine hist_field_accumulate_1d

   !#######################################################################

   subroutine hist_field_accumulate_2d(field, data, cols_or_block,      &
        cole, logger)
      use hist_msg_handler, only: hist_log_messages, hist_have_error, ERROR
      use hist_msg_handler, only: hist_add_message, VERBOSE
      use hist_field,       only: hist_field_info_t
      use hist_buffer,      only: hist_buff_2d_t, hist_buffer_t
      use ISO_FORTRAN_ENV,  only: REAL64

      ! Dummy arguments
      class(hist_field_info_t), pointer,  intent(inout) :: field
      real(REAL64),                       intent(in)    :: data(:,:)
      integer,                            intent(in)    :: cols_or_block
      integer,                  optional, intent(in)    :: cole
      type(hist_log_messages),  optional, intent(inout) :: logger
      ! Local variables
      class(hist_buffer_t), pointer     :: buff_ptr
      class(hist_buff_2d_t), pointer    :: buff
      character(len=*),     parameter   :: subname = 'hist_field_accumulate_2d'

      if (associated(field)) then
         buff_ptr => field%buffers
         do
            if (associated(buff_ptr) .and.                                    &
                 (.not. hist_have_error(errors=logger))) then
               select type (buff_ptr)
               class is (hist_buff_2d_t)
                  buff => buff_ptr
                  call buff%accumulate(data, cols_or_block, field%flag_xyfill(), &
                          field%fill_value(), cole, logger)
                  if (hist_have_error(errors=logger)) then
                     call  logger%add_stack_frame(ERROR, __FILE__, __LINE__ - 3, &
                          subname=subname)
                     exit
                  else
                     call hist_add_message(subname, VERBOSE,                     &
                          "Accumulated data for",                                &
                          msgstr2=trim(field%diag_name())//", Buffer type, ",    &
                          msgstr3=trim(buff%buffer_type()),                  &
                          logger=logger)
                  end if
               end select
               buff_ptr => buff_ptr%next
            else
               exit
            end if
         end do
      end if ! No else, it is legit to pass in a null pointer

   end subroutine hist_field_accumulate_2d

   !#######################################################################

   subroutine hist_field_norm_value_1d(field, norm_val, logger)
      use hist_buffer,      only: hist_buffer_t, hist_buff_1d_t
      use hist_msg_handler, only: hist_log_messages, hist_have_error, ERROR
      use hist_msg_handler, only: hist_add_message, VERBOSE
      use hist_field,       only: hist_field_info_t
      use ISO_FORTRAN_ENV,  only: REAL64

      ! Dummy arguments
      class(hist_field_info_t),          intent(inout) :: field
      real(REAL64),                      intent(inout) :: norm_val(:)
      type(hist_log_messages), optional, intent(inout) :: logger
      ! Local variables
      class(hist_buffer_t), pointer  :: buff_ptr
      class(hist_buff_1d_t), pointer :: buff
      character(len=*), parameter    :: subname = 'hist_field_norm_value_1d'

      buff_ptr => field%buffers
      if (associated(buff_ptr) .and.                                  &
         (.not. hist_have_error(errors=logger))) then
         select type(buff_ptr)
         class is (hist_buff_1d_t)
            buff => buff_ptr
            call buff%norm_value(norm_val, field%flag_xyfill(), &
                    field%fill_value(), logger=logger)
            if (hist_have_error(errors=logger)) then
               call logger%add_stack_frame(ERROR, __FILE__, __LINE__-3,  &
                  subname=subname)
            else
               call hist_add_message(subname, VERBOSE,                   &
                  "Accumulated data for",                                &
                  msgstr2=trim(field%diag_name())//", Buffer type, ",    &
                  msgstr3=trim(buff%buffer_type()),                  &
                  logger=logger)
            end if
         end select
      end if

   end subroutine hist_field_norm_value_1d

   !#######################################################################

   subroutine hist_field_norm_value_2d(field, norm_val, logger)
      use hist_buffer,      only: hist_buffer_t, hist_buff_2d_t
      use hist_msg_handler, only: hist_log_messages, hist_have_error, ERROR
      use hist_msg_handler, only: hist_add_message, VERBOSE
      use hist_field,       only: hist_field_info_t
      use ISO_FORTRAN_ENV,  only: REAL64

      ! Dummy arguments
      class(hist_field_info_t),          intent(inout) :: field
      real(REAL64),                      intent(inout) :: norm_val(:,:)
      type(hist_log_messages), optional, intent(inout) :: logger
      ! Local variables
      class(hist_buffer_t), pointer  :: buff_ptr
      class(hist_buff_2d_t), pointer :: buff
      character(len=*), parameter    :: subname = 'hist_field_norm_value_2d'

      buff_ptr => field%buffers
      if (associated(buff_ptr) .and.                                  &
         (.not. hist_have_error(errors=logger))) then
         select type(buff_ptr)
         class is (hist_buff_2d_t)
            buff => buff_ptr
            call buff%norm_value(norm_val, field%flag_xyfill(), &
                    field%fill_value(), logger=logger)
            if (hist_have_error(errors=logger)) then
              call logger%add_stack_frame(ERROR, __FILE__, __LINE__-3,  &
                 subname=subname)
            else
               call hist_add_message(subname, VERBOSE,                   &
                  "Accumulated data for",                                &
                  msgstr2=trim(field%diag_name())//", Buffer type, ",    &
                  msgstr3=trim(buff_ptr%buffer_type()),                  &
                  logger=logger)
            end if
         end select
      end if

   end subroutine hist_field_norm_value_2d

end module hist_api
