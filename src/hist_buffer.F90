module hist_buffer
   use ISO_FORTRAN_ENV, only: REAL64, REAL32, INT32, INT64
   use hist_hashable,   only: hist_hashable_t

   implicit none
   private

   ! Public interfaces
   public :: buffer_factory

   ! Accumulation types -- the integers and array positions below must match
   integer, parameter, public :: hist_accum_lst = 1 ! last sample
   integer, parameter, public :: hist_accum_min = 2 ! minimum sample
   integer, parameter, public :: hist_accum_max = 3 ! maximum sample
   integer, parameter, public :: hist_accum_avg = 4 ! sample average
   integer, parameter, public :: hist_accum_var = 5 ! sample standard deviation

   integer, parameter         :: as_len = 36
   character(len=as_len), parameter, public :: accum_strings(5) = (/          &
        'last sampled value                  ',                               &
        'minimum of sampled values           ',                               &
        'maximum of sampled values           ',                               &
        'average of sampled values           ',                               &
        'standard deviation of sampled values' /)

   character(len=3), parameter, public :: accum_abbrev(5) =                   &
        (/ 'lst', 'min', 'max', 'avg', 'var' /)

   ! Time sampling flag indices
   !!XXgoldyXX: Todo: decide on sampling types

   type, abstract, public :: hist_buffer_t
      ! hist_buffer_t is an abstract base class for hist_outfld buffers
      class(hist_hashable_t), pointer              :: field_info => NULL()
      integer,                             private :: vol = -1 ! For host output
      integer,                             private :: horiz_axis_ind = 0
      integer,                             private :: rank = 0
      integer,                             private :: accum_type = 0
      integer,                allocatable, private :: field_shape(:)
      integer,                allocatable, private :: num_samples(:)
      integer,                allocatable, private :: block_begs(:)
      integer,                allocatable, private :: block_ends(:)
      character(len=:),       allocatable, private :: buff_type
      class(hist_buffer_t),   pointer              :: next => NULL()
   contains
      procedure                              :: field  => get_field_info
      procedure                              :: volume => get_volume
      procedure                              :: horiz_axis_index
      procedure                              :: init_buffer
      procedure                              :: accum_string
      procedure                              :: buffer_type
      procedure                              :: clear => hist_buff_clear
      procedure                              :: check_status
      procedure                              :: has_blocks
      procedure(hist_buff_init),    deferred :: initialize
   end type hist_buffer_t

   type, public, extends(hist_buffer_t) :: hist_buff_1dreal32_inst_t
      real(REAL32), pointer :: data(:) => NULL()
   CONTAINS
      procedure :: clear => buff_1dreal32_inst_clear
      procedure :: accumulate => buff_1dreal32_inst_accum
      procedure :: norm_value => buff_1dreal32_inst_value
      procedure :: initialize => init_buff_inst_1dreal32
   end type hist_buff_1dreal32_inst_t

   type, public, extends(hist_buff_1dreal32_inst_t) :: hist_buff_1dreal32_avg_t
   CONTAINS
      procedure :: accumulate => buff_1dreal32_avg_accum
      procedure :: norm_value => buff_1dreal32_avg_value
      procedure :: initialize => init_buff_avg_1dreal32
   end type hist_buff_1dreal32_avg_t

   type, public, extends(hist_buffer_t) :: hist_buff_1dreal64_inst_t
      real(REAL64), pointer :: data(:) => NULL()
   CONTAINS
      procedure :: clear => buff_1dreal64_inst_clear
      procedure :: accumulate => buff_1dreal64_inst_accum
      procedure :: norm_value => buff_1dreal64_inst_value
      procedure :: initialize => init_buff_inst_1dreal64
   end type hist_buff_1dreal64_inst_t

   ! Abstract interfaces for hist_buffer_t class
   abstract interface
      subroutine hist_buff_sub_log(this, logger)
         use hist_msg_handler, only: hist_log_messages
         import                   :: hist_buffer_t
         class(hist_buffer_t),              intent(inout) :: this
         type(hist_log_messages), optional, intent(inout) :: logger
      end subroutine hist_buff_sub_log
   end interface

   abstract interface
      subroutine hist_buff_init(this, field_in, volume_in, horiz_axis_in,     &
           accum_type_in, shape_in, block_sizes_in, block_ind_in, logger)
         use hist_msg_handler, only: hist_log_messages
         import                   :: hist_buffer_t
         import                   :: hist_hashable_t
         class(hist_buffer_t),              intent(inout) :: this
         class(hist_hashable_t),  pointer                 :: field_in
         integer,                           intent(in)    :: volume_in
         integer,                           intent(in)    :: horiz_axis_in
         integer,                           intent(in)    :: accum_type_in
         integer,                           intent(in)    :: shape_in(:)
         integer,                 optional, intent(in)    :: block_sizes_in(:)
         integer,                 optional, intent(in)    :: block_ind_in
         type(hist_log_messages), optional, intent(inout) :: logger
      end subroutine hist_buff_init
   end interface

CONTAINS

   !#######################################################################

   function get_field_info(this)
      class(hist_buffer_t), intent(in) :: this
      class(hist_hashable_t), pointer  :: get_field_info

      get_field_info => this%field_info
   end function get_field_info

   !#######################################################################

   integer function get_volume(this)
      class(hist_buffer_t), intent(in) :: this

      get_volume = this%vol
   end function get_volume

   !#######################################################################

   integer function horiz_axis_index(this)
      class(hist_buffer_t), intent(in) :: this

      horiz_axis_index = this%horiz_axis_ind
   end function horiz_axis_index
   !#######################################################################

   logical function check_status(this, logger, filename, line)
      ! Check to see if this buffer is properly initialized
      use hist_msg_handler, only: hist_log_messages, hist_add_error, ERROR

      ! Dummy arguments
      class(hist_buffer_t),              intent(inout) :: this
      type(hist_log_messages), optional, intent(inout) :: logger
      character(len=*),        optional, intent(in)    :: filename
      integer,                 optional, intent(in)    :: line
      ! Local variable
      character(len=*), parameter :: subname = 'check_status'

      check_status = .true.
      if ( (this%horiz_axis_index() < 1)      .or.                            &
           (.not. allocated(this%field_shape)) .or.                           &
           (.not. allocated(this%num_samples))) then
         check_status = .false.
         call hist_add_error(subname,                                         &
              "buffer not properly initialized '", errors=logger)
         if (present(filename) .and. present(line) .and. present(logger)) then
            call logger%add_stack_frame(ERROR, filename, line)
         end if
      end if

   end function check_status

   !#######################################################################

   logical function has_blocks(this)
      ! Dummy argument
      class(hist_buffer_t),   intent(inout) :: this

      has_blocks = allocated(this%block_begs) .and. allocated(this%block_ends)

   end function has_blocks

   !#######################################################################

   subroutine hist_buff_clear(this, logger)
      use hist_msg_handler, only: hist_log_messages, hist_add_alloc_error

      ! Dummy arguments
      class(hist_buffer_t),              intent(inout) :: this
      type(hist_log_messages), optional, intent(inout) :: logger
      ! Local variables
      integer                     :: aerr
      character(len=*), parameter :: subname = 'hist_buff_clear'

      if (this%check_status(logger, __FILE__, __LINE__ + 1)) then
         this%num_samples = 0
      end if

   end subroutine hist_buff_clear

   !#######################################################################

   subroutine init_buffer(this, field_in, volume_in, horiz_axis_in,           &
        accum_type_in, shape_in, block_sizes_in, block_ind_in, logger)
      use hist_msg_handler, only: hist_log_messages, hist_add_alloc_error
      use hist_msg_handler, only: hist_add_error

      ! Dummy arguments
      class(hist_buffer_t),              intent(inout) :: this
      class(hist_hashable_t),  pointer                 :: field_in
      integer,                           intent(in)    :: volume_in
      integer,                           intent(in)    :: horiz_axis_in
      integer,                           intent(in)    :: accum_type_in
      integer,                           intent(in)    :: shape_in(:)
      integer,                 optional, intent(in)    :: block_sizes_in(:)
      integer,                 optional, intent(in)    :: block_ind_in
      type(hist_log_messages), optional, intent(inout) :: logger
      ! Local variables
      integer                     :: astat
      integer                     :: hsize
      character(len=*), parameter :: subname = 'init_buffer'

      ! Sanity check
      this%rank = SIZE(shape_in, 1)
      if ((horiz_axis_in < 1) .or. (horiz_axis_in > this%rank)) then
         call hist_add_error(subname,                                         &
              "horiz_axis_in must be between 1 and '", errint1=this%rank,     &
              errors=logger)
      else
         this%field_info => field_in
         this%vol = volume_in
         this%horiz_axis_ind = horiz_axis_in
         this%accum_type = accum_type_in
         allocate(this%field_shape(size(shape_in, 1)), stat=astat)
         if (astat == 0) then
            this%field_shape(:) = shape_in(:)
         else
            call hist_add_alloc_error('field_shape', __FILE__, __LINE__ - 4,  &
                 subname=subname, errors=logger)
         end if
         ! Allocate num_samples
         hsize = this%field_shape(this%horiz_axis_ind)
         allocate(this%num_samples(hsize), stat=astat)
         if (astat /= 0) then
            call hist_add_alloc_error('num_samples', __FILE__, __LINE__ - 2,  &
                 subname=subname, errors=logger)
         end if
      end if
   end subroutine init_buffer

   !#######################################################################

   function accum_string(this) result(ac_str)
      class(hist_buffer_t), intent(in)   :: this
      character(len=as_len), allocatable :: ac_str

      ac_str = accum_strings(this%accum_type)

   end function accum_string

   !#######################################################################

   function buffer_type(this) result(bt_str)
      class(hist_buffer_t), intent(in) :: this
      character(len=:), allocatable    :: bt_str

      if (allocated(this%buff_type)) then
         bt_str = this%buff_type
      else
         bt_str = 'UNKNOWN'
      end if

   end function buffer_type

   !#######################################################################

   subroutine buff_1dreal32_inst_clear(this, logger)
      use hist_msg_handler, only: hist_log_messages, hist_add_alloc_error

      ! Dummy arguments
      class(hist_buff_1dreal32_inst_t),   intent(inout) :: this
      type(hist_log_messages), optional, intent(inout) :: logger
      ! Local variables
      integer                     :: aerr
      character(len=*), parameter :: subname = 'buff_1dreal32_inst_clear'

      call hist_buff_clear(this, logger)
      if (.not. associated(this%data)) then
         allocate(this%data(this%field_shape(1)), stat=aerr)
         if (aerr /= 0) then
            call hist_add_alloc_error('data', __FILE__, __LINE__ - 2,         &
                 subname=subname, errors=logger)
         end if
      end if
      this%data = 0.0_REAL32

   end subroutine buff_1dreal32_inst_clear

   !#######################################################################

   subroutine init_buff_inst_1dreal32(this, field_in, volume_in,              &
        horiz_axis_in, accum_type_in, shape_in, block_sizes_in, block_ind_in, &
        logger)
      use hist_msg_handler, only: hist_log_messages, hist_add_alloc_error

      ! Dummy arguments
      class(hist_buff_1dreal32_inst_t),  intent(inout) :: this
      class(hist_hashable_t),  pointer                 :: field_in
      integer,                           intent(in)    :: volume_in
      integer,                           intent(in)    :: horiz_axis_in
      integer,                           intent(in)    :: accum_type_in
      integer,                           intent(in)    :: shape_in(:)
      integer,                 optional, intent(in)    :: block_sizes_in(:)
      integer,                 optional, intent(in)    :: block_ind_in
      type(hist_log_messages), optional, intent(inout) :: logger

      call init_buffer(this, field_in, volume_in, horiz_axis_in,              &
           accum_type_in, shape_in, block_sizes_in, block_ind_in, logger=logger)
      call this%clear(logger=logger)
      this%buff_type = 'buff_1dreal32_'//accum_abbrev(accum_type_in)

   end subroutine init_buff_inst_1dreal32

   !#######################################################################

   subroutine buff_1dreal32_inst_accum(this, field, cols_or_block, cole, logger)
      use hist_msg_handler, only: hist_log_messages
      ! Dummy arguments
      class(hist_buff_1dreal32_inst_t),  intent(inout) :: this
      real(REAL32),                      intent(in)    :: field(:)
      integer,                           intent(in)    :: cols_or_block
      integer,                 optional, intent(in)    :: cole
      type(hist_log_messages), optional, intent(inout) :: logger
      ! Local variables
      integer      :: col_beg_use
      integer      :: col_end_use
      integer      :: ind1
      real(REAL32) :: fld_val

      if (this%has_blocks()) then
         ! For a blocked field, <cols_or_block> is a block index
         col_beg_use = this%block_begs(cols_or_block)
         col_end_use = this%block_ends(cols_or_block)
      else
         ! Non blocked, <cols_or_block> is the first column index
         col_beg_use = cols_or_block
         if (present(cole)) then
            col_end_use = cole
         else
            col_end_use = col_beg_use +                                       &
                 this%field_shape(this%horiz_axis_ind) - 1
         end if
      end if

      select case (this%accum_type)
      case (hist_accum_lst)
         this%data(col_beg_use:col_end_use) = field(:)
         this%num_samples(col_beg_use:col_end_use) = 1
      case (hist_accum_min)
         do ind1 = col_beg_use, col_end_use
            fld_val = field(ind1 - col_beg_use + 1)
            if (this%num_samples(ind1) == 0) then
               this%data(ind1) = fld_val
            else if (fld_val < this%data(ind1)) then
               this%data(ind1) = fld_val
            end if ! No else, we already have the minimum value for this col
            this%num_samples(ind1) = 1
         end do
      case (hist_accum_max)
         do ind1 = col_beg_use, col_end_use
            fld_val = field(ind1 - col_beg_use + 1)
            if (this%num_samples(ind1) == 0) then
               this%data(ind1) = fld_val
            else if (fld_val > this%data(ind1)) then
               this%data(ind1) = fld_val
            end if ! No else, we already have the maximum value for this col
            this%num_samples(ind1) = 1
         end do
      end select

   end subroutine buff_1dreal32_inst_accum

   !#######################################################################

   subroutine buff_1dreal32_inst_value(this, norm_val, default_val, logger)
      use hist_msg_handler, only: hist_log_messages

      ! Dummy arguments
      class(hist_buff_1dreal32_inst_t),  intent(inout) :: this
      real(REAL32),                      intent(inout) :: norm_val(:)
      real(REAL32),            optional, intent(in)    :: default_val
      type(hist_log_messages), optional, intent(inout) :: logger
      ! Local variable
      integer :: ind1
      integer :: nacc

      do ind1 = 1, this%field_shape(1)
         nacc = this%num_samples(ind1)
         if (nacc > 0) then
            norm_val(ind1) = this%data(ind1)
         else if (present(default_val)) then
            norm_val(ind1) = default_val
         end if
      end do

   end subroutine buff_1dreal32_inst_value

   !#######################################################################

   subroutine init_buff_avg_1dreal32(this, field_in, volume_in,               &
        horiz_axis_in, accum_type_in, shape_in, block_sizes_in, block_ind_in, &
        logger)
      use hist_msg_handler, only: hist_log_messages, hist_add_alloc_error

      ! Dummy arguments
      class(hist_buff_1dreal32_avg_t),   intent(inout) :: this
      class(hist_hashable_t),  pointer                 :: field_in
      integer,                           intent(in)    :: volume_in
      integer,                           intent(in)    :: horiz_axis_in
      integer,                           intent(in)    :: accum_type_in
      integer,                           intent(in)    :: shape_in(:)
      integer,                 optional, intent(in)    :: block_sizes_in(:)
      integer,                 optional, intent(in)    :: block_ind_in
      type(hist_log_messages), optional, intent(inout) :: logger

      call init_buffer(this, field_in, volume_in, horiz_axis_in,              &
           accum_type_in, shape_in, block_sizes_in, block_ind_in, logger=logger)
      call this%clear(logger=logger)
      this%buff_type = 'buff_1dreal32_'//accum_abbrev(accum_type_in)

   end subroutine init_buff_avg_1dreal32

   !#######################################################################

   subroutine buff_1dreal32_avg_accum(this, field, cols_or_block, cole, logger)
      use hist_msg_handler, only: hist_log_messages
      ! Dummy arguments
      class(hist_buff_1dreal32_avg_t),   intent(inout) :: this
      real(REAL32),                      intent(in)    :: field(:)
      integer,                           intent(in)    :: cols_or_block
      integer,                 optional, intent(in)    :: cole
      type(hist_log_messages), optional, intent(inout) :: logger
      ! Local variables
      integer :: col_beg_use
      integer :: col_end_use
      integer :: index

      if (this%has_blocks()) then
         ! For a blocked field, <cols_or_block> is a block index
         col_beg_use = this%block_begs(cols_or_block)
         col_end_use = this%block_ends(cols_or_block)
      else
         ! Non blocked, <cols_or_block> is the first column index
         col_beg_use = cols_or_block
         if (present(cole)) then
            col_end_use = cole
         else
            col_end_use = col_beg_use +                                       &
                 this%field_shape(this%horiz_axis_ind) - 1
         end if
      end if

      do index = col_beg_use, col_end_use
         this%data(index) = this%data(index) + field(index - col_beg_use + 1)
         this%num_samples(index) = this%num_samples(index) + 1
      end do

   end subroutine buff_1dreal32_avg_accum

   !#######################################################################

   subroutine buff_1dreal32_avg_value(this, norm_val, default_val, logger)
      use hist_msg_handler, only: hist_log_messages

      ! Dummy arguments
      class(hist_buff_1dreal32_avg_t),   intent(inout) :: this
      real(REAL32),                      intent(inout) :: norm_val(:)
      real(REAL32),            optional, intent(in)    :: default_val
      type(hist_log_messages), optional, intent(inout) :: logger
      ! Local variable
      integer :: ind1
      integer :: nacc

      do ind1 = 1, this%field_shape(1)
         nacc = this%num_samples(ind1)
         if (nacc > 0) then
            norm_val(ind1) = this%data(ind1) / real(nacc, REAL32)
         else if (present(default_val)) then
            norm_val(ind1) = default_val
         end if
      end do

   end subroutine buff_1dreal32_avg_value

   !#######################################################################

   subroutine buff_1dreal64_inst_clear(this, logger)
      use hist_msg_handler, only: hist_log_messages, hist_add_alloc_error

      ! Dummy arguments
      class(hist_buff_1dreal64_inst_t),  intent(inout) :: this
      type(hist_log_messages), optional, intent(inout) :: logger
      ! Local variables
      integer                     :: aerr
      character(len=*), parameter :: subname = 'buff_1dreal64_inst_clear'

      call hist_buff_clear(this, logger)
      if (.not. associated(this%data)) then
         allocate(this%data(this%field_shape(1)), stat=aerr)
         if (aerr /= 0) then
            call hist_add_alloc_error('data', __FILE__, __LINE__ - 1,         &
                 subname=subname, errors=logger)
         end if
      end if
      this%data = 0.0_REAL64

   end subroutine buff_1dreal64_inst_clear

   !#######################################################################

   subroutine init_buff_inst_1dreal64(this, field_in, volume_in,              &
        horiz_axis_in, accum_type_in, shape_in, block_sizes_in, block_ind_in, &
        logger)
      use hist_msg_handler, only: hist_log_messages

      class(hist_buff_1dreal64_inst_t),   intent(inout) :: this
      class(hist_hashable_t),  pointer                 :: field_in
      integer,                           intent(in)    :: volume_in
      integer,                           intent(in)    :: horiz_axis_in
      integer,                           intent(in)    :: accum_type_in
      integer,                           intent(in)    :: shape_in(:)
      integer,                 optional, intent(in)    :: block_sizes_in(:)
      integer,                 optional, intent(in)    :: block_ind_in
      type(hist_log_messages), optional, intent(inout) :: logger

      call init_buffer(this, field_in, volume_in, horiz_axis_in,              &
           accum_type_in, shape_in, block_sizes_in, block_ind_in, logger=logger)
      call this%clear(logger=logger)
      this%buff_type = 'buff_1dreal64_'//accum_abbrev(accum_type_in)

   end subroutine init_buff_inst_1dreal64

   !#######################################################################

   subroutine buff_1dreal64_inst_accum(this, field, cols_or_block, cole, logger)
      use hist_msg_handler, only: hist_log_messages
      ! Dummy arguments
      class(hist_buff_1dreal64_inst_t),  intent(inout) :: this
      real(REAL64),                      intent(in)    :: field(:)
      integer,                           intent(in)    :: cols_or_block
      integer,                 optional, intent(in)    :: cole
      type(hist_log_messages), optional, intent(inout) :: logger
      ! Local variables
      integer      :: col_beg_use
      integer      :: col_end_use
      integer      :: ind1
      real(REAL64) :: fld_val

      if (this%has_blocks()) then
         ! For a blocked field, <cols_or_block> is a block index
         col_beg_use = this%block_begs(cols_or_block)
         col_end_use = this%block_ends(cols_or_block)
      else
         ! Non blocked, <cols_or_block> is the first column index
         col_beg_use = cols_or_block
         if (present(cole)) then
            col_end_use = cole
         else
            col_end_use = col_beg_use +                                       &
                 this%field_shape(this%horiz_axis_ind) - 1
         end if
      end if


      select case (this%accum_type)
      case (hist_accum_lst)
         this%data(col_beg_use:col_end_use) = field(:)
         this%num_samples(col_beg_use:col_end_use) = 1
      case (hist_accum_min)
         do ind1 = col_beg_use, col_end_use
            fld_val = field(ind1 - col_beg_use + 1)
            if (this%num_samples(ind1) == 0) then
               this%data(ind1) = fld_val
            else if (fld_val < this%data(ind1)) then
               this%data(ind1) = fld_val
            end if ! No else, we already have the minimum value for this col
            this%num_samples(ind1) = 1
         end do
      case (hist_accum_max)
         do ind1 = col_beg_use, col_end_use
            fld_val = field(ind1 - col_beg_use + 1)
            if (this%num_samples(ind1) == 0) then
               this%data(ind1) = fld_val
            else if (fld_val > this%data(ind1)) then
               this%data(ind1) = fld_val
            end if ! No else, we already have the maximum value for this col
            this%num_samples(ind1) = 1
         end do
      end select

   end subroutine buff_1dreal64_inst_accum

   !#######################################################################

   subroutine buff_1dreal64_inst_value(this, norm_val, default_val, logger)
      use hist_msg_handler, only: hist_log_messages, ERROR, VERBOSE
      ! Dummy arguments
      class(hist_buff_1dreal64_inst_t),  intent(inout) :: this
      real(REAL64),                      intent(inout) :: norm_val(:)
      real(REAL64),            optional, intent(in)    :: default_val
      type(hist_log_messages), optional, intent(inout) :: logger
      ! Local variable
      integer :: ind1
      integer :: nacc

      do ind1 = 1, this%field_shape(1)
         nacc = this%num_samples(ind1)
         if (nacc > 0) then
            norm_val(ind1) = this%data(ind1)
         else if (present(default_val)) then
            norm_val(ind1) = default_val
         end if
      end do

      norm_val(:) = this%data(:)

   end subroutine buff_1dreal64_inst_value

   !#######################################################################

   function buffer_factory(buffer_type, logger) result(newbuf)
      ! Create a new buffer based on <buffer_type>.
      ! <buffer_type> has a format typekind_rank_accum
      ! Where:
      ! <typekind> is a lowercase string representation of a
      !    supported kind from the ISO_FORTRAN_ENV module.
      ! <rank> is the rank of the buffer (no leading zeros)
      ! <accum> is the accumulation type, one of:
      !    lst: Store the last value collected
      !    avg: Accumulate running average
      !    var: Accumulate standard deviation
      !    min: Accumulate smallest value
      !    max: Accumulate largest value

      use hist_msg_handler, only: hist_log_messages, ERROR, VERBOSE
      use hist_msg_handler, only: hist_add_error

      ! Arguments
      class(hist_buffer_t),    pointer                 :: newbuf
      character(len=*),                  intent(in)    :: buffer_type
      type(hist_log_messages), optional, intent(inout) :: logger

      ! Local variables
      character(len=*),                parameter :: subname = 'buffer_factory'
      ! For buffer                     allocation
      integer                                    :: aerr
      type(hist_buff_1dreal32_inst_t), pointer   :: real32_1_inst => NULL()
      type(hist_buff_1dreal32_avg_t),  pointer   :: real32_1_avg  => NULL()
      type(hist_buff_1dreal64_inst_t), pointer   :: real64_1_inst => NULL()


      nullify(newbuf)
      ! Create new buffer
      select case (trim(buffer_type))
      case ('real32_1_lst', 'real32_1_min', 'real32_1_max')
         allocate(real32_1_inst, stat=aerr)
         if (aerr == 0) then
            newbuf => real32_1_inst
         end if
      case ('real32_1_avg', 'real32_1_var')
         allocate(real32_1_avg, stat=aerr)
         if (aerr == 0) then
            newbuf => real32_1_avg
         end if
      case ('real64_1_lst', 'real64_1_min', 'real64_1_max')
         allocate(real64_1_inst, stat=aerr)
         if (aerr == 0) then
            newbuf => real64_1_inst
         end if
      case default
         call hist_add_error(subname,                                         &
              "Invalid or unsupported buffer type, '",                        &
              errstr2=trim(buffer_type), errstr3="'", errors=logger)
      end select
   end function buffer_factory

end module hist_buffer
