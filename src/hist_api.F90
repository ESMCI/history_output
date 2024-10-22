module hist_api

   use ISO_FORTRAN_ENV, only: REAL64, REAL32, INT32, INT64
   use hist_buffer,     only: hist_buff_1dreal32_t
   use hist_buffer,     only: hist_buff_1dreal64_t
   use hist_buffer,     only: hist_buff_2dreal32_t
   use hist_buffer,     only: hist_buff_2dreal64_t

   implicit none
   private

   ! Public API interfaces
   public :: hist_new_field         ! Allocate a hist_field_info_t object
   public :: hist_new_buffer        ! Create a new field buffer
   public :: hist_field_accumulate  ! Accumulate a new field state in all buffs
   public :: hist_buffer_accumulate ! Accumulate a new field state
   public :: hist_buffer_norm_value ! Return current normalized field state

   ! Interfaces for public interfaces
   interface hist_field_accumulate
      module procedure hist_field_accumulate_1d_REAL32
      module procedure hist_field_accumulate_1d_REAL64
      module procedure hist_field_accumulate_2d_REAL32
      module procedure hist_field_accumulate_2d_REAL64
   end interface hist_field_accumulate

   interface hist_buffer_accumulate
      module procedure hist_buffer_accumulate_1d_REAL32
      module procedure hist_buffer_accumulate_1d_REAL64
      module procedure hist_buffer_accumulate_2d_REAL32
      module procedure hist_buffer_accumulate_2d_REAL64
   end interface hist_buffer_accumulate

   interface hist_buffer_norm_value
      module procedure hist_buffer_norm_value_1dreal32
      module procedure hist_buffer_norm_value_1dreal64
      module procedure hist_buffer_norm_value_2dreal32
      module procedure hist_buffer_norm_value_2dreal64
   end interface hist_buffer_norm_value

CONTAINS

   !#######################################################################

   function hist_new_field(diag_name_in, std_name_in, long_name_in, units_in, &
        type_in, decomp_in, dimensions, acc_flag, num_levels, field_shape, sampling_seq,  &
        flag_xyfill, mixing_ratio, dim_bounds, mdim_sizes, beg_dims, end_dims,  &
        cell_methods, errors) result(new_field)
      use hist_msg_handler, only: hist_have_error, hist_log_messages, ERROR
      use hist_field,       only: hist_field_initialize, hist_field_info_t

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
              mixing_ratio=mixing_ratio, dim_bounds=dim_bounds, mdim_sizes=mdim_sizes, &
              beg_dims=beg_dims, end_dims=end_dims, cell_methods=cell_methods, errmsg=errmsg)
         if (hist_have_error(errors=errors)) then
            call errors%add_stack_frame(ERROR, __FILE__, __LINE__ - 3,        &
                 subname=subname)
         end if
      end if
   end function hist_new_field

   !#######################################################################

   subroutine hist_new_buffer(field, buff_shape, buff_kind, horiz_axis_ind,   &
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
      integer,                            intent(in)    :: buff_kind
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
      select case (type_str)
      case ('integer')
         select case (buff_kind)
         case (INT32)
            kind_string = 'int32'
         case (INT64)
            kind_string = 'int64'
         case default
            kind_string = ''
         end select
      case ('real')
         select case(buff_kind)
         case (REAL32)
            kind_string = 'real32'
         case (REAL64)
            kind_string = 'real64'
         case default
            kind_string = ''
         end select
      case default
         kind_string = ''
         call hist_add_error(subname, "type, '"//type_str,                    &
              errstr2="' is not supported", errors=errors)
      end select
      if ((len_trim(kind_string) == 0) .and. present(errors)) then
         call errors%new_error("kind = ", errint1=buff_kind,                  &
              errstr2=" is not supported for type "//type_str, subname=subname)
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
      ! We now know what sort of buffer we need
      ! First, sort by rank
      select case (rank)
      case (1)
         ! sort by kind (already checked above)
         if (buff_kind == REAL32) then
            bufftype_string = 'real32_1'
         else if (buff_kind == REAL64) then
            bufftype_string = 'real64_1'
         end if
      case (2)
         ! sort by kind
         if (buff_kind == REAL32) then
            bufftype_string = 'real32_2'
         else if (buff_kind == REAL64) then
            bufftype_string = 'real64_2'
         end if
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
              accum_val, buffer_shape, block_sizes, block_ind, logger=errors)
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

   subroutine hist_field_accumulate_1d_REAL32(field, data, cols_or_block,      &
        cole, logger)
      use hist_msg_handler, only: hist_log_messages, hist_have_error, ERROR
      use hist_msg_handler, only: hist_add_message, VERBOSE
      use hist_field,       only: hist_field_info_t
      use hist_buffer,      only: hist_buffer_t

      ! Dummy arguments
      class(hist_field_info_t), pointer,  intent(inout) :: field
      real(REAL32),                       intent(in)    :: data(:)
      integer,                            intent(in)    :: cols_or_block
      integer,                  optional, intent(in)    :: cole
      type(hist_log_messages),  optional, intent(inout) :: logger
      ! Local variables
      class(hist_buffer_t), pointer     :: buff_ptr
      character(len=:),     allocatable :: buff_typestr
      character(len=*), parameter :: subname = 'hist_field_accumulate_1d_REAL32'

      if (associated(field)) then
         buff_ptr => field%buffers
         do
            if (associated(buff_ptr) .and.                                    &
                 (.not. hist_have_error(errors=logger))) then
               call hist_buffer_accumulate(buff_ptr, data, cols_or_block,     &
                    cole=cole, logger=logger)
               if (hist_have_error(errors=logger)) then
                  call  logger%add_stack_frame(ERROR, __FILE__, __LINE__ - 3, &
                       subname=subname)
                  exit
               else
                  call hist_add_message(subname, VERBOSE,                     &
                       "Accumulated data for",                                &
                       msgstr2=trim(field%diag_name())//", Buffer type, ",    &
                       msgstr3=trim(buff_ptr%buffer_type()),                  &
                       logger=logger)
                  buff_ptr => buff_ptr%next
               end if
            else
               exit
            end if
         end do
      end if ! No else, it is legit to pass in a null pointer

   end subroutine hist_field_accumulate_1d_REAL32

   !#######################################################################

   subroutine hist_field_accumulate_1d_REAL64(field, data, cols_or_block,      &
        cole, logger)
      use hist_msg_handler, only: hist_log_messages, hist_have_error, ERROR
      use hist_msg_handler, only: hist_add_message, VERBOSE
      use hist_field,       only: hist_field_info_t
      use hist_buffer,      only: hist_buffer_t

      ! Dummy arguments
      class(hist_field_info_t), pointer,  intent(inout) :: field
      real(REAL64),                       intent(in)    :: data(:)
      integer,                            intent(in)    :: cols_or_block
      integer,                  optional, intent(in)    :: cole
      type(hist_log_messages),  optional, intent(inout) :: logger
      ! Local variables
      class(hist_buffer_t), pointer     :: buff_ptr
      character(len=:),     allocatable :: buff_typestr
      character(len=*), parameter :: subname = 'hist_field_accumulate_1d_REAL64'

      if (associated(field)) then
         buff_ptr => field%buffers
         do
            if (associated(buff_ptr) .and.                                    &
                 (.not. hist_have_error(errors=logger))) then
               call hist_buffer_accumulate(buff_ptr, data, cols_or_block,     &
                    cole=cole, logger=logger)
               if (hist_have_error(errors=logger)) then
                  call  logger%add_stack_frame(ERROR, __FILE__, __LINE__ - 3, &
                       subname=subname)
                  exit
               else
                  call hist_add_message(subname, VERBOSE,                     &
                       "Accumulated data for",                                &
                       msgstr2=trim(field%diag_name())//", Buffer type, ",    &
                       msgstr3=trim(buff_ptr%buffer_type()),                  &
                       logger=logger)
                  buff_ptr => buff_ptr%next
               end if
            else
               exit
            end if
         end do
      end if ! No else, it is legit to pass in a null pointer

   end subroutine hist_field_accumulate_1d_REAL64

   !#######################################################################

   subroutine hist_field_accumulate_2d_REAL32(field, data, cols_or_block,      &
        cole, logger)
      use hist_msg_handler, only: hist_log_messages, hist_have_error, ERROR
      use hist_msg_handler, only: hist_add_message, VERBOSE
      use hist_field,       only: hist_field_info_t
      use hist_buffer,      only: hist_buffer_t

      ! Dummy arguments
      class(hist_field_info_t), pointer,  intent(inout) :: field
      real(REAL32),                       intent(in)    :: data(:,:)
      integer,                            intent(in)    :: cols_or_block
      integer,                  optional, intent(in)    :: cole
      type(hist_log_messages),  optional, intent(inout) :: logger
      ! Local variables
      class(hist_buffer_t), pointer     :: buff_ptr
      character(len=:),     allocatable :: buff_typestr
      character(len=*),     parameter   :: subname = 'hist_field_accumulate_2d_REAL32'

      if (associated(field)) then
         buff_ptr => field%buffers
         do
            if (associated(buff_ptr) .and.                                    &
                 (.not. hist_have_error(errors=logger))) then
               call hist_buffer_accumulate(buff_ptr, data, cols_or_block,     &
                    cole=cole, logger=logger)
               if (hist_have_error(errors=logger)) then
                  call  logger%add_stack_frame(ERROR, __FILE__, __LINE__ - 3, &
                       subname=subname)
                  exit
               else
                  call hist_add_message(subname, VERBOSE,                     &
                       "Accumulated data for",                                &
                       msgstr2=trim(field%diag_name())//", Buffer type, ",    &
                       msgstr3=trim(buff_ptr%buffer_type()),                  &
                       logger=logger)
                  buff_ptr => buff_ptr%next
               end if
            else
               exit
            end if
         end do
      end if ! No else, it is legit to pass in a null pointer

   end subroutine hist_field_accumulate_2d_REAL32

   !#######################################################################

   subroutine hist_field_accumulate_2d_REAL64(field, data, cols_or_block,      &
        cole, logger)
      use hist_msg_handler, only: hist_log_messages, hist_have_error, ERROR
      use hist_msg_handler, only: hist_add_message, VERBOSE
      use hist_field,       only: hist_field_info_t
      use hist_buffer,      only: hist_buffer_t

      ! Dummy arguments
      class(hist_field_info_t), pointer,  intent(inout) :: field
      real(REAL64),                       intent(in)    :: data(:,:)
      integer,                            intent(in)    :: cols_or_block
      integer,                  optional, intent(in)    :: cole
      type(hist_log_messages),  optional, intent(inout) :: logger
      ! Local variables
      class(hist_buffer_t), pointer     :: buff_ptr
      character(len=:),     allocatable :: buff_typestr
      character(len=*),     parameter   :: subname = 'hist_field_accumulate_2d_REAL64'

      if (associated(field)) then
         buff_ptr => field%buffers
         do
            if (associated(buff_ptr) .and.                                    &
                 (.not. hist_have_error(errors=logger))) then
               call hist_buffer_accumulate(buff_ptr, data, cols_or_block,     &
                    cole=cole, logger=logger)
               if (hist_have_error(errors=logger)) then
                  call  logger%add_stack_frame(ERROR, __FILE__, __LINE__ - 3, &
                       subname=subname)
                  exit
               else
                  call hist_add_message(subname, VERBOSE,                     &
                       "Accumulated data for",                                &
                       msgstr2=trim(field%diag_name())//", Buffer type, ",    &
                       msgstr3=trim(buff_ptr%buffer_type()),                  &
                       logger=logger)
                  buff_ptr => buff_ptr%next
               end if
            else
               exit
            end if
         end do
      end if ! No else, it is legit to pass in a null pointer

   end subroutine hist_field_accumulate_2d_REAL64

   !#######################################################################

   subroutine hist_buffer_accumulate_1d_REAL32(buffer, field, cols_or_block,   &
        cole, logger)
      use hist_msg_handler, only: hist_log_messages, hist_add_error
      use hist_buffer,      only: hist_buffer_t

      ! Dummy arguments
      class(hist_buffer_t),    target,   intent(inout) :: buffer
      real(REAL32),                      intent(in)    :: field(:)
      integer,                           intent(in)    :: cols_or_block
      integer,                 optional, intent(in)    :: cole
      type(hist_log_messages), optional, intent(inout) :: logger
      ! Local variables
      class(hist_buff_1dreal32_t), pointer          :: buff32
      class(hist_buff_1dreal64_t), pointer          :: buff64
      character(len=:),                 allocatable :: buff_typestr
      character(len=*), parameter :: subname = 'hist_buffer_accumulate_1d_REAL32'

      select type(buffer)
      class is (hist_buff_1dreal32_t)
         buff32 => buffer
         call buff32%accumulate(field, cols_or_block, cole, logger)
      class is (hist_buff_1dreal64_t)
         ! Do we want to accumulate 32bit data into 64bit buffers?
         buff64 => buffer
         call buff64%accumulate(real(field, REAL64), cols_or_block, cole,     &
              logger)
      class default
         buff_typestr = buffer%buffer_type()
         call hist_add_error(subname, "unsupported buffer type, '",           &
              errstr2=buff_typestr, errstr3="'", errors=logger)
      end select

   end subroutine hist_buffer_accumulate_1d_REAL32

   !#######################################################################

   subroutine hist_buffer_accumulate_1d_REAL64(buffer, field, cols_or_block,   &
        cole, logger)
      use hist_msg_handler, only: hist_log_messages, hist_add_error
      use hist_buffer,      only: hist_buffer_t

      ! Dummy arguments
      class(hist_buffer_t),    target,   intent(inout) :: buffer
      real(REAL64),                      intent(in)    :: field(:)
      integer,                           intent(in)    :: cols_or_block
      integer,                 optional, intent(in)    :: cole
      type(hist_log_messages), optional, intent(inout) :: logger
      ! Local variables
      class(hist_buff_1dreal32_t), pointer          :: buff32
      class(hist_buff_1dreal64_t), pointer          :: buff64
      character(len=:),                 allocatable :: buff_typestr
      character(len=*), parameter :: subname = 'hist_buffer_accumulate_1d_REAL64'

      select type(buffer)
      class is (hist_buff_1dreal32_t)
         buff32 => buffer
         call buff32%accumulate(real(field, REAL32), cols_or_block, cole, logger)
      class is (hist_buff_1dreal64_t)
         buff64 => buffer
         call buff64%accumulate(field, cols_or_block, cole, logger)
      class default
         buff_typestr = buffer%buffer_type()
         call hist_add_error(subname, "unsupported buffer type, '",           &
              errstr2=buff_typestr, errstr3="'", errors=logger)
      end select

   end subroutine hist_buffer_accumulate_1d_REAL64

   !#######################################################################

   subroutine hist_buffer_accumulate_2d_REAL32(buffer, field, cols_or_block,   &
        cole, logger)
      use hist_msg_handler, only: hist_log_messages, hist_add_error
      use hist_buffer,      only: hist_buffer_t

      ! Dummy arguments
      class(hist_buffer_t),    target,   intent(inout) :: buffer
      real(REAL32),                      intent(in)    :: field(:,:)
      integer,                           intent(in)    :: cols_or_block
      integer,                 optional, intent(in)    :: cole
      type(hist_log_messages), optional, intent(inout) :: logger
      ! Local variables
      type(hist_buff_2dreal32_t), pointer          :: buff32
      type(hist_buff_2dreal64_t), pointer          :: buff64
      character(len=:),                allocatable :: buff_typestr
      character(len=*), parameter :: subname = 'hist_buffer_accumulate_2d_REAL32'

      select type(buffer)
      class is (hist_buff_2dreal32_t)
         buff32 => buffer
         call buff32%accumulate(field, cols_or_block, cole, logger)
      class is (hist_buff_2dreal64_t)
         buff64 => buffer
         call buff64%accumulate(real(field, REAL64), cols_or_block, cole, logger)
      class default
         buff_typestr = buffer%buffer_type()
         call hist_add_error(subname, "unsupported buffer type, '",           &
              errstr2=buff_typestr, errstr3="'", errors=logger)
      end select

   end subroutine hist_buffer_accumulate_2d_REAL32

   !#######################################################################

   subroutine hist_buffer_accumulate_2d_REAL64(buffer, field, cols_or_block,   &
        cole, logger)
      use hist_msg_handler, only: hist_log_messages, hist_add_error
      use hist_buffer,      only: hist_buffer_t

      ! Dummy arguments
      class(hist_buffer_t),    target,   intent(inout) :: buffer
      real(REAL64),                      intent(in)    :: field(:,:)
      integer,                           intent(in)    :: cols_or_block
      integer,                 optional, intent(in)    :: cole
      type(hist_log_messages), optional, intent(inout) :: logger
      ! Local variables
      type(hist_buff_2dreal32_t), pointer          :: buff32
      type(hist_buff_2dreal64_t), pointer          :: buff64
      character(len=:),                allocatable :: buff_typestr
      character(len=*), parameter :: subname = 'hist_buffer_accumulate_2d_REAL64'

      select type(buffer)
      class is (hist_buff_2dreal32_t)
         buff32 => buffer
         call buff32%accumulate(real(field, REAL32), cols_or_block, cole, logger)
      class is (hist_buff_2dreal64_t)
         buff64 => buffer
         call buff64%accumulate(field, cols_or_block, cole, logger)
      class default
         buff_typestr = buffer%buffer_type()
         call hist_add_error(subname, "unsupported buffer type, '",           &
              errstr2=buff_typestr, errstr3="'", errors=logger)
      end select

   end subroutine hist_buffer_accumulate_2d_REAL64

   !#######################################################################

   subroutine hist_buffer_norm_value_1dreal32(buffer, norm_val, default_val,  &
        logger)
      use hist_msg_handler, only: hist_log_messages, hist_add_error
      use hist_buffer,      only: hist_buffer_t

      ! Dummy arguments
      class(hist_buffer_t),    target,   intent(inout) :: buffer
      real(REAL32),                      intent(inout) :: norm_val(:)
      real(REAL32),            optional, intent(in)    :: default_val
      type(hist_log_messages), optional, intent(inout) :: logger
      ! Local variables
      class(hist_buff_1dreal32_t), pointer          :: buff32
      class(hist_buff_1dreal64_t), pointer          :: buff64
      real(REAL64),                     allocatable :: norm_val64(:)
      character(len=:),                 allocatable :: buff_typestr
      character(len=*), parameter :: subname = 'hist_buffer_norm_value_1dreal32'

      select type(buffer)
      class is (hist_buff_1dreal32_t)
         buff32 => buffer
         call buff32%norm_value(norm_val, default_val=default_val,            &
              logger=logger)
      class is (hist_buff_1dreal64_t)
         ! Truncate 64bit buffer into 32bit output
         buff64 => buffer
         allocate(norm_val64(size(norm_val, 1)))
         if (present(default_val)) then
            call buff64%norm_value(norm_val64,                                &
                 default_val=REAL(default_val, REAL64), logger=logger)
         else
            call buff64%norm_value(norm_val64, logger=logger)
         end if
         norm_val(:) = REAL(norm_val64(:), REAL32)
      class default
         buff_typestr = buffer%buffer_type()
         call hist_add_error(subname, "unsupported buffer type, '",           &
              errstr2=buff_typestr, errstr3="'", errors=logger)
      end select

   end subroutine hist_buffer_norm_value_1dreal32

   !#######################################################################

   subroutine hist_buffer_norm_value_1dreal64(buffer, norm_val, default_val,  &
      logger)
      use hist_msg_handler, only: hist_log_messages, hist_add_error
      use hist_buffer,      only: hist_buffer_t

      ! Dummy arguments
      class(hist_buffer_t),    target,   intent(inout) :: buffer
      real(REAL64),                      intent(inout) :: norm_val(:)
      real(REAL64),            optional, intent(in)    :: default_val
      type(hist_log_messages), optional, intent(inout) :: logger
      ! Local variables
      type(hist_buff_1dreal32_t), pointer          :: buff32
      type(hist_buff_1dreal64_t), pointer          :: buff64
      real(REAL32),                    allocatable :: norm_val32(:)
      character(len=:),                allocatable :: buff_typestr
      character(len=*), parameter :: subname = 'hist_buffer_norm_value_1dreal64'

      select type(buffer)
      class is (hist_buff_1dreal32_t)
         ! Do we want to read out 32bit buffers into 64bit data?
         buff32 => buffer
         allocate(norm_val32(size(norm_val, 1)))
         if (present(default_val)) then
            call buffer%norm_value(norm_val32,                                &
                 default_val=REAL(default_val, REAL32), logger=logger)
         else
            call buffer%norm_value(norm_val32, logger=logger)
         end if
         norm_val(:) = REAL(norm_val32(:), REAL64)
      class is (hist_buff_1dreal64_t)
         buff64 => buffer
         call buffer%norm_value(norm_val, default_val=default_val,            &
              logger=logger)
      class default
         buff_typestr = buffer%buffer_type()
         call hist_add_error(subname, "unsupported buffer type, '",           &
              errstr2=buff_typestr, errstr3="'", errors=logger)
      end select

   end subroutine hist_buffer_norm_value_1dreal64

   !#######################################################################

   subroutine hist_buffer_norm_value_2dreal32(buffer, norm_val, default_val,  &
        logger)
      use hist_msg_handler, only: hist_log_messages, hist_add_error
      use hist_buffer,      only: hist_buffer_t

      ! Dummy arguments
      class(hist_buffer_t),    target,   intent(inout) :: buffer
      real(REAL32),                      intent(inout) :: norm_val(:,:)
      real(REAL32),            optional, intent(in)    :: default_val
      type(hist_log_messages), optional, intent(inout) :: logger
      ! Local variables
      class(hist_buff_2dreal32_t), pointer          :: buff32
      class(hist_buff_2dreal64_t), pointer          :: buff64
      real(REAL64),                     allocatable :: norm_val64(:,:)
      character(len=:),                 allocatable :: buff_typestr
      character(len=*), parameter :: subname = 'hist_buffer_norm_value_2dreal32'

      select type(buffer)
      class is (hist_buff_2dreal32_t)
         buff32 => buffer
         call buff32%norm_value(norm_val, default_val=default_val,            &
              logger=logger)
      class is (hist_buff_2dreal64_t)
         ! Truncate 64bit buffer into 32bit output
         buff64 => buffer
         allocate(norm_val64(size(norm_val, 1), size(norm_val, 2)))
         if (present(default_val)) then
            call buff64%norm_value(norm_val64,                                &
                 default_val=REAL(default_val, REAL64), logger=logger)
         else
            call buff64%norm_value(norm_val64, logger=logger)
         end if
         norm_val(:,:) = REAL(norm_val64(:,:), REAL32)
      class default
         buff_typestr = buffer%buffer_type()
         call hist_add_error(subname, "unsupported buffer type, '",           &
              errstr2=buff_typestr, errstr3="'", errors=logger)
      end select

   end subroutine hist_buffer_norm_value_2dreal32

   !#######################################################################

   subroutine hist_buffer_norm_value_2dreal64(buffer, norm_val, default_val,  &
      logger)
      use hist_msg_handler, only: hist_log_messages, hist_add_error
      use hist_buffer,      only: hist_buffer_t

      ! Dummy arguments
      class(hist_buffer_t),    target,   intent(inout) :: buffer
      real(REAL64),                      intent(inout) :: norm_val(:,:)
      real(REAL64),            optional, intent(in)    :: default_val
      type(hist_log_messages), optional, intent(inout) :: logger
      ! Local variables
      type(hist_buff_2dreal32_t), pointer          :: buff32
      type(hist_buff_2dreal64_t), pointer          :: buff64
      real(REAL32),                    allocatable :: norm_val32(:,:)
      character(len=:),                allocatable :: buff_typestr
      character(len=*), parameter :: subname = 'hist_buffer_norm_value_1dreal64'

      select type(buffer)
      class is (hist_buff_2dreal32_t)
         ! Do we want to read out 32bit buffers into 64bit data?
         buff32 => buffer
         allocate(norm_val32(size(norm_val, 1), size(norm_val, 2)))
         if (present(default_val)) then
            call buffer%norm_value(norm_val32,                                &
                 default_val=REAL(default_val, REAL32), logger=logger)
         else
            call buffer%norm_value(norm_val32, logger=logger)
         end if
         norm_val(:,:) = REAL(norm_val32(:,:), REAL64)
      class is (hist_buff_2dreal64_t)
         buff64 => buffer
         call buffer%norm_value(norm_val, default_val=default_val,            &
              logger=logger)
      class default
         buff_typestr = buffer%buffer_type()
         call hist_add_error(subname, "unsupported buffer type, '",           &
              errstr2=buff_typestr, errstr3="'", errors=logger)
      end select

   end subroutine hist_buffer_norm_value_2dreal64

end module hist_api
